"""Structure fixing / cleanup (PDBFixer integration -- planned in the
original spec as `prot clean` but never implemented until now).

Unlike Scwrl4/MODELLER/Foldseek, PDBFixer (built on OpenMM) is free and
pip-installable (`pip install -e ".[fix]"`) -- no license, no separate
binary to hunt down. So this module treats it as a normal optional
extra rather than an external-tool-only wrapper, and raises a clear
PDBFixerNotAvailableError (with the install command, not a purchase/
registration pointer) if the extra isn't installed.

What this adds that `prot model`/`prot mutate` don't cover:
- Missing *atoms* within an existing residue (e.g. a sidechain
  incomplete due to weak density) -- prot model/prot mutate only ever
  handle whole missing residues or full identity changes, never
  "this residue is here but incomplete."
- Nonstandard-residue normalization (e.g. MSE -> MET).
- Heterogen removal (waters, ions, ligands).
- Hydrogen addition at a given pH.

What it does NOT replace:
- `prot model gaps` for missing-residue *detection*. PDBFixer's
  findMissingResidues() only sees gaps recorded in the file's SEQRES
  header (or an explicitly supplied sequence) -- a numbering-only gap
  with no SEQRES is invisible to it. Verified directly against this
  project's own gapped.pdb test fixture (residues 3-5 missing, inferred
  purely from the 2 -> 6 numbering jump): PDBFixer's
  findMissingResidues() reports no gap there at all, while
  `prot model gaps` (which only looks at numbering, no SEQRES needed)
  correctly reports it. Use `prot model gaps` for detection; this
  module's add_missing_residues option can optionally *fill* gaps
  PDBFixer does know about (from SEQRES), with more realistic
  template-based geometry than `prot model loop`'s straight-line
  placeholder trace -- but still not a substitute for real loop
  modeling/refinement for anything you'll trust downstream.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


class PDBFixerNotAvailableError(RuntimeError):
    pass


def _require_pdbfixer():
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        raise PDBFixerNotAvailableError(
            "PDBFixer (and OpenMM) are not installed. Install them with "
            "`pip install -e \".[fix]\"` (or `pip install pdbfixer openmm`) -- "
            "both are free/open-source and pip-installable, no license needed."
        ) from exc
    return PDBFixer, PDBFile


def _residue_label(residue) -> str:
    return f"{residue.chain.id}/{residue.name}{residue.id}"


@dataclass
class FixAnalysis:
    missing_residues: dict[str, list[str]]  # "chain_id/insertion_point" -> [resnames]
    missing_atoms: dict[str, list[str]]  # residue label -> [atom names]
    missing_terminals: dict[str, list[str]]  # residue label -> [atom names, usually OXT]
    nonstandard_residues: list[tuple[str, str]]  # (residue label, replacement resname)


def analyze(pdb_path: str | Path) -> FixAnalysis:
    """Report what PDBFixer would find, without changing anything --
    for deciding what a subsequent `fix()` call should actually do."""
    PDBFixer, _ = _require_pdbfixer()
    fixer = PDBFixer(filename=str(pdb_path))

    fixer.findMissingResidues()
    missing_residues = {
        f"{fixer.topology._chains[chain_idx].id}/insert_at_{insert_idx}": resnames
        for (chain_idx, insert_idx), resnames in fixer.missingResidues.items()
    }

    fixer.findMissingAtoms()
    missing_atoms = {
        _residue_label(residue): [a.name for a in atoms]
        for residue, atoms in fixer.missingAtoms.items()
    }
    missing_terminals = {
        _residue_label(residue): names
        for residue, names in fixer.missingTerminals.items()
    }

    fixer.findNonstandardResidues()
    nonstandard = [(_residue_label(residue), new_name) for residue, new_name in fixer.nonstandardResidues]

    return FixAnalysis(
        missing_residues=missing_residues, missing_atoms=missing_atoms,
        missing_terminals=missing_terminals, nonstandard_residues=nonstandard,
    )


@dataclass
class FixReport:
    residues_added: dict[str, list[str]] = field(default_factory=dict)
    atoms_added: dict[str, list[str]] = field(default_factory=dict)
    nonstandard_replaced: list[tuple[str, str]] = field(default_factory=list)
    heterogens_removed: bool = False
    hydrogens_added_at_ph: float | None = None


def fix(
    pdb_path: str | Path,
    output_path: str | Path,
    add_missing_residues: bool = False,
    add_missing_atoms: bool = True,
    replace_nonstandard: bool = True,
    remove_heterogens: str | None = None,  # None | "water" | "all"
    add_hydrogens_ph: float | None = None,
) -> FixReport:
    """Run the requested PDBFixer repair steps and write the result to
    `output_path`. Every step is opt-in/opt-out via the arguments above;
    nothing happens that wasn't asked for.

    add_missing_residues=False (the default) means PDBFixer's own
    missing-residue detection is *not* acted on even if it finds
    something (e.g. via SEQRES) -- only atom-level completions run. Set
    it True to also insert whole missing residues (PDBFixer's own
    template-based geometry, applied only where its detection actually
    found something -- see the module docstring for its SEQRES
    limitation).
    """
    PDBFixer, PDBFile = _require_pdbfixer()
    fixer = PDBFixer(filename=str(pdb_path))
    report = FixReport()

    fixer.findMissingResidues()
    if not add_missing_residues:
        fixer.missingResidues = {}
    else:
        report.residues_added = {
            f"{fixer.topology._chains[ci].id}/insert_at_{ii}": names
            for (ci, ii), names in fixer.missingResidues.items()
        }

    if replace_nonstandard:
        fixer.findNonstandardResidues()
        report.nonstandard_replaced = [
            (_residue_label(r), name) for r, name in fixer.nonstandardResidues
        ]
        fixer.replaceNonstandardResidues()

    if remove_heterogens is not None:
        keep_water = remove_heterogens == "water"
        fixer.removeHeterogens(keepWater=keep_water)
        report.heterogens_removed = True

    if add_missing_atoms:
        fixer.findMissingAtoms()
        report.atoms_added = {
            _residue_label(r): [a.name for a in atoms] for r, atoms in fixer.missingAtoms.items()
        }
        fixer.addMissingAtoms()

    if add_hydrogens_ph is not None:
        fixer.addMissingHydrogens(add_hydrogens_ph)
        report.hydrogens_added_at_ph = add_hydrogens_ph

    output_path = Path(output_path)
    with open(output_path, "w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle)

    return report
