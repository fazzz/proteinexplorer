"""Structure validation (spec-adjacent -- not in the original draft, added
after discussing MolProbity integration).

MolProbity's real value (Ramachandran favored/allowed/outlier zones,
rotamer outlier classification, CaBLAM) rests on statistically-calibrated
reference distributions from large curated structure sets. Reproducing
those thresholds from memory risks quietly getting them wrong and
presenting fabricated statistics as validated science -- the same concern
that kept a real rotamer library out of mutate.py's cb_only fallback. So
this module does NOT attempt Ramachandran/rotamer outlier classification
itself. What it does provide, dependency-free and without needing any
calibrated reference data:

- clashes(): steric (van der Waals) overlap between non-bonded atoms,
  reusing pocket.py's VDW_RADII. This is real, checkable geometry, not a
  statistical judgment call -- two atoms nearer than the sum of their
  vdW radii (minus a tolerance) are physically clashing, full stop.
  Heavy-atom-only input (the usual case for a crystal structure) means
  hydrogen-involving clashes are invisible here; consider running
  `prot fix apply --add-hydrogens <pH>` first for a more complete (if
  still not clash-optimized) picture.
- bond_geometry(): backbone bond lengths/angles compared against
  standard idealized covalent values (textbook bond chemistry, not
  empirical statistics) -- N-CA/CA-C/C-N lengths and the three backbone
  angles, each with a generous tolerance.

For the real thing -- Ramachandran/rotamer outliers, clashscore,
CaBLAM -- molprobity() wraps an external MolProbity installation
(Phenix's `phenix.molprobity`, or a standalone `mmtbx.molprobity`/
`molprobity.molprobity`), same external-tool-only pattern as
predict.py/search.py: no license-gate issue like Scwrl4/MODELLER, but
still a compiled/packaged tool this environment doesn't have, so a clear
error with an install pointer when it's missing. Output parsing is
deliberately shallow (the full MolProbity report format varies by
version/install method and isn't something to guess at) -- the raw
report text is always returned; a couple of commonly-present summary
numbers are extracted opportunistically, when they're findable, without
assuming a specific format.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Structure import Structure

from proteinexplorer.models import ResidueCategory, classify_residue
from proteinexplorer.pocket import VDW_RADII, _DEFAULT_VDW

CLASH_TOLERANCE = 0.4  # Angstrom, the same default Probe/MolProbity use


def _vdw(atom) -> float:
    return VDW_RADII.get((atom.element or "").strip().upper(), _DEFAULT_VDW)


def _residue_label(residue) -> str:
    return f"{residue.get_parent().id}/{residue.resname}{residue.id[1]}"


# --- Steric clashes ------------------------------------------------------

@dataclass
class Clash:
    atom_a: str  # "chain/resname+resnum:atomname"
    atom_b: str
    distance: float
    overlap: float  # how far into each other's vdW spheres, in Angstrom


def _protein_or_ligand_atoms(structure: Structure):
    model = next(iter(structure))
    return [
        a for a in model.get_atoms()
        if classify_residue(a.get_parent().resname, a.get_parent().id[0]) is not ResidueCategory.WATER
    ]


def _disulfide_bonded_residue_pairs(atoms: list, cutoff: float = 2.5) -> set[tuple[int, int]]:
    """Residue-object-id pairs connected by a real disulfide bond (SG-SG
    within a typical S-S bond distance) -- excluded from clash detection
    since they're legitimately, covalently close, not clashing. Same
    distance criterion as contact.py's find_disulfide_bonds."""
    sg_atoms = [
        a for a in atoms
        if a.get_name() == "SG" and a.get_parent().resname.strip().upper() == "CYS"
    ]
    pairs = set()
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            d = float(np.linalg.norm(sg_atoms[i].coord - sg_atoms[j].coord))
            if d <= cutoff:
                res_i, res_j = sg_atoms[i].get_parent(), sg_atoms[j].get_parent()
                pairs.add(tuple(sorted((id(res_i), id(res_j)))))
    return pairs


def clashes(
    structure: Structure,
    atoms: list | None = None,
    tolerance: float = CLASH_TOLERANCE,
    exclude_adjacent_residues: bool = True,
) -> list[Clash]:
    """Non-bonded steric overlaps. Excludes atom pairs within the same
    residue, between sequence-adjacent residues on the same chain (their
    backbone connection legitimately brings atoms close -- without an
    explicit bond topology, distinguishing "bonded" from "clashing"
    there isn't reliable, so adjacent-residue pairs are conservatively
    skipped rather than risking false positives), and between residues
    connected by a real disulfide bond.
    """
    all_atoms = atoms if atoms is not None else _protein_or_ligand_atoms(structure)
    if len(all_atoms) < 2:
        return []

    disulfide_pairs = _disulfide_bonded_residue_pairs(all_atoms)
    ns = NeighborSearch(all_atoms)
    max_vdw = max((_vdw(a) for a in all_atoms), default=_DEFAULT_VDW)
    max_reach = 2 * max_vdw - tolerance

    seen = set()
    found: list[Clash] = []
    for atom in all_atoms:
        nearby = ns.search(atom.coord, max_reach, level="A")
        for other in nearby:
            if other is atom:
                continue
            key = tuple(sorted((id(atom), id(other))))
            if key in seen:
                continue
            seen.add(key)

            res_a, res_b = atom.get_parent(), other.get_parent()
            if res_a is res_b:
                continue
            if tuple(sorted((id(res_a), id(res_b)))) in disulfide_pairs:
                continue
            if exclude_adjacent_residues:
                chain_a, chain_b = res_a.get_parent(), res_b.get_parent()
                if chain_a is chain_b and abs(res_a.id[1] - res_b.id[1]) <= 1:
                    continue

            distance = float(np.linalg.norm(atom.coord - other.coord))
            allowed = _vdw(atom) + _vdw(other) - tolerance
            if distance < allowed:
                found.append(
                    Clash(
                        atom_a=f"{_residue_label(res_a)}:{atom.get_name()}",
                        atom_b=f"{_residue_label(res_b)}:{other.get_name()}",
                        distance=distance,
                        overlap=allowed - distance,
                    )
                )
    found.sort(key=lambda c: -c.overlap)
    return found


# --- Backbone bond geometry ------------------------------------------

# Standard idealized covalent bond lengths (Angstrom) and backbone valence
# angles (degrees) -- textbook peptide geometry (Engh & Huber-style
# consensus values), not empirically-fit statistics.
IDEAL_BOND_LENGTHS = {"N-CA": 1.458, "CA-C": 1.525, "C-N": 1.329}
BOND_LENGTH_TOLERANCE = 0.12

IDEAL_BOND_ANGLES = {"N-CA-C": 111.2, "CA-C-N": 116.6, "C-N-CA": 121.7}
BOND_ANGLE_TOLERANCE = 8.0


@dataclass
class BondOutlier:
    kind: str  # e.g. "N-CA" or "N-CA-C"
    residue: str
    value: float
    ideal: float
    deviation: float


def _protein_residue_sequence(chain):
    return [
        r for r in chain
        if classify_residue(r.resname, r.id[0]) is ResidueCategory.PROTEIN
    ]


def bond_geometry(structure: Structure) -> list[BondOutlier]:
    """Backbone bond lengths/angles more than a generous tolerance away
    from standard idealized covalent values. A handful of outliers is
    normal in any real structure (strain, refinement noise); a cluster of
    them in one region is worth a closer look."""
    model = next(iter(structure))
    outliers: list[BondOutlier] = []

    for chain in model:
        residues = _protein_residue_sequence(chain)
        for residue in residues:
            if not all(name in residue for name in ("N", "CA", "C")):
                continue
            n, ca, c = residue["N"].coord, residue["CA"].coord, residue["C"].coord

            d_n_ca = float(np.linalg.norm(ca - n))
            d_ca_c = float(np.linalg.norm(c - ca))
            label = _residue_label(residue)
            if abs(d_n_ca - IDEAL_BOND_LENGTHS["N-CA"]) > BOND_LENGTH_TOLERANCE:
                outliers.append(BondOutlier("N-CA", label, d_n_ca, IDEAL_BOND_LENGTHS["N-CA"], d_n_ca - IDEAL_BOND_LENGTHS["N-CA"]))
            if abs(d_ca_c - IDEAL_BOND_LENGTHS["CA-C"]) > BOND_LENGTH_TOLERANCE:
                outliers.append(BondOutlier("CA-C", label, d_ca_c, IDEAL_BOND_LENGTHS["CA-C"], d_ca_c - IDEAL_BOND_LENGTHS["CA-C"]))

            v1 = n - ca
            v2 = c - ca
            cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            angle_n_ca_c = float(np.degrees(np.arccos(np.clip(cos_t, -1.0, 1.0))))
            if abs(angle_n_ca_c - IDEAL_BOND_ANGLES["N-CA-C"]) > BOND_ANGLE_TOLERANCE:
                outliers.append(
                    BondOutlier("N-CA-C", label, angle_n_ca_c, IDEAL_BOND_ANGLES["N-CA-C"],
                                angle_n_ca_c - IDEAL_BOND_ANGLES["N-CA-C"])
                )

        # inter-residue: C(i)-N(i+1) length, CA(i)-C(i)-N(i+1) and
        # C(i)-N(i+1)-CA(i+1) angles
        for prev_res, next_res in zip(residues, residues[1:]):
            if not all(name in prev_res for name in ("CA", "C")):
                continue
            if "N" not in next_res:
                continue
            c_i = prev_res["C"].coord
            n_next = next_res["N"].coord
            d_c_n = float(np.linalg.norm(n_next - c_i))
            label = _residue_label(prev_res)
            if abs(d_c_n - IDEAL_BOND_LENGTHS["C-N"]) > BOND_LENGTH_TOLERANCE:
                outliers.append(BondOutlier("C-N", label, d_c_n, IDEAL_BOND_LENGTHS["C-N"], d_c_n - IDEAL_BOND_LENGTHS["C-N"]))

            ca_i = prev_res["CA"].coord
            v1 = ca_i - c_i
            v2 = n_next - c_i
            cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            angle_ca_c_n = float(np.degrees(np.arccos(np.clip(cos_t, -1.0, 1.0))))
            if abs(angle_ca_c_n - IDEAL_BOND_ANGLES["CA-C-N"]) > BOND_ANGLE_TOLERANCE:
                outliers.append(
                    BondOutlier("CA-C-N", label, angle_ca_c_n, IDEAL_BOND_ANGLES["CA-C-N"],
                                angle_ca_c_n - IDEAL_BOND_ANGLES["CA-C-N"])
                )

    return sorted(outliers, key=lambda o: -abs(o.deviation))


# --- External: MolProbity ---------------------------------------------

class MolProbityNotAvailableError(RuntimeError):
    pass


_MOLPROBITY_BINARY_CANDIDATES = ["phenix.molprobity", "mmtbx.molprobity", "molprobity.molprobity"]


def molprobity_binary() -> str | None:
    for candidate in _MOLPROBITY_BINARY_CANDIDATES:
        found = shutil.which(candidate)
        if found is not None:
            return found
    return None


@dataclass
class MolProbityResult:
    binary: str
    raw_output: str
    summary: dict[str, str] = field(default_factory=dict)


def molprobity(pdb_path: str | Path, timeout: int = 900) -> MolProbityResult:
    """Run an external MolProbity installation for the real,
    statistically-calibrated validation report (Ramachandran/rotamer
    outliers, clashscore, CaBLAM, ...). Output format varies by
    install/version, so parsing is deliberately shallow: the full report
    text is always returned in `raw_output`; a few common summary lines
    are pulled into `summary` opportunistically when a recognizable
    "<label>: <value>" pattern is present, without assuming the rest of
    the report's structure.
    """
    binary = molprobity_binary()
    if binary is None:
        raise MolProbityNotAvailableError(
            "No MolProbity installation found on PATH (tried: "
            + ", ".join(_MOLPROBITY_BINARY_CANDIDATES) + "). MolProbity ships as "
            "part of Phenix (https://phenix-online.org/) or as a standalone build "
            "(https://github.com/rlabduke/MolProbity) -- there is no dependency-free "
            "substitute for its calibrated Ramachandran/rotamer/clashscore analysis; "
            "see `prot valid clashes`/`prot valid geometry` for what this package can "
            "check without it."
        )

    result = subprocess.run([binary, str(pdb_path)], capture_output=True, text=True, timeout=timeout)
    output = result.stdout + result.stderr

    summary = {}
    for line in output.splitlines():
        match = re.match(r"\s*([A-Za-z][A-Za-z0-9 _/-]{2,40}?)\s*[:=]\s*([-+]?\d+\.?\d*)\s*$", line)
        if match:
            summary[match.group(1).strip()] = match.group(2).strip()

    return MolProbityResult(binary=binary, raw_output=output, summary=summary)
