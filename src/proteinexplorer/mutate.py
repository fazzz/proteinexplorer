"""Point mutation / residue editing (spec section "Mutation").

Two independent side-chain construction methods, same two-tier pattern as
secondary.py (dssp/geometric) and pocket.py (no external tool available):

- scwrl4(): wraps an external Scwrl4 binary for full, energetically
  reasonable side-chain placement (dead-end elimination over a rotamer
  library). This is the real thing and should be preferred whenever
  Scwrl4 is installed.
- cb_only(): a dependency-free fallback that renames the residue, keeps
  the original backbone (N/CA/C/O) untouched, and places an idealized
  virtual C-beta using the standard backbone-only C-beta reconstruction
  formula. It does NOT attempt to build the rest of the side chain --
  going beyond C-beta needs residue-specific bond lengths/angles and
  rotamer statistics that are out of scope for a dependency-free
  fallback, so this is intentionally a partial result. This is honestly
  reported to the caller (`MutationResult.note`), not disguised as a
  complete side chain.

mutate_residue() picks scwrl4 when available and falls back to cb_only,
mirroring secondary_structure()'s auto/dssp/geometric selection.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB.Structure import Structure

from proteinexplorer import io as pio
from proteinexplorer.models import PROTEIN_RESIDUES, STANDARD_AMINO_ACIDS

ONE_TO_THREE: dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


class MutationError(ValueError):
    pass


class Scwrl4NotAvailableError(RuntimeError):
    pass


def normalize_resname(target: str) -> str:
    """Accept either a 1-letter or 3-letter amino acid code and return the
    canonical 3-letter PDB resname."""
    target = target.strip().upper()
    if len(target) == 1:
        resname = ONE_TO_THREE.get(target)
        if resname is None:
            raise MutationError(f"Unknown one-letter amino acid code: {target!r}")
        return resname
    if target in STANDARD_AMINO_ACIDS:
        return target
    raise MutationError(f"Unknown amino acid: {target!r} (expected a standard 1- or 3-letter code)")


@dataclass
class MutationResult:
    chain_id: str
    resseq: int
    original_resname: str
    new_resname: str
    method: str  # "scwrl4" or "cb_only"
    atoms_placed: list[str]
    note: str


def _find_residue(structure: Structure, chain_id: str, resseq: int):
    model = next(iter(structure))
    if chain_id not in model:
        raise MutationError(f"No chain {chain_id!r} in this structure")
    chain = model[chain_id]
    for residue in chain:
        if residue.id[1] == resseq and residue.id[0] == " ":
            return chain, residue
    raise MutationError(f"No standard residue at {chain_id}/{resseq}")


# --- Built-in fallback: backbone + idealized virtual C-beta ------------

def _virtual_cb(n_coord: np.ndarray, ca_coord: np.ndarray, c_coord: np.ndarray) -> np.ndarray:
    """Idealized C-beta position from backbone N/CA/C only, using the
    standard virtual-C-beta reconstruction formula widely used in protein
    structure tools (e.g. for placing C-beta on glycine or an incomplete
    side chain). Chirality follows the usual PDB N->CA->C atom ordering
    (L-amino acid convention); bond length/angle plausibility is covered
    by tests, but absolute chirality has not been cross-checked against
    an experimental reference structure in this environment (no network
    access to fetch one) -- flagged here for anyone relying on this for
    chirality-sensitive work.
    """
    b1 = ca_coord - n_coord
    b2 = c_coord - ca_coord
    a = np.cross(b1, b2)
    cb = -0.58273431 * a + 0.56802827 * b1 - 0.54067466 * b2 + ca_coord
    return cb


def cb_only(structure: Structure, chain_id: str, resseq: int, target_resname: str) -> MutationResult:
    """Rename a residue and rebuild it as backbone (N/CA/C/O, unchanged)
    plus an idealized C-beta (omitted for GLY). Does not build the rest
    of the side chain -- see module docstring."""
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Residue import Residue

    chain, residue = _find_residue(structure, chain_id, resseq)
    original_resname = residue.resname.strip().upper()

    required = ("N", "CA", "C")
    missing = [name for name in required if name not in residue]
    if missing:
        raise MutationError(
            f"Cannot mutate {chain_id}/{original_resname}{resseq}: missing backbone atom(s) {missing}"
        )

    n_coord, ca_coord, c_coord = (residue[name].coord for name in required)
    o_atom = residue["O"] if "O" in residue else None

    new_residue = Residue((" ", resseq, " "), target_resname, residue.segid)
    placed: list[str] = []

    def add(name: str, coord: np.ndarray, element: str):
        new_residue.add(Atom(name, coord, 0.0, 1.0, " ", name, 0, element=element))
        placed.append(name)

    add("N", n_coord, "N")
    add("CA", ca_coord, "C")
    add("C", c_coord, "C")
    if o_atom is not None:
        add("O", o_atom.coord, "O")

    note = "built-in fallback: backbone kept, C-beta idealized from N/CA/C."
    if target_resname != "GLY":
        cb_coord = _virtual_cb(n_coord, ca_coord, c_coord)
        add("CB", cb_coord, "C")
        note += " Full side chain beyond C-beta not built -- install Scwrl4 for a complete rotamer."
    else:
        note += " Glycine has no C-beta."

    chain.detach_child(residue.id)
    chain.add(new_residue)

    return MutationResult(
        chain_id=chain_id, resseq=resseq, original_resname=original_resname,
        new_resname=target_resname, method="cb_only", atoms_placed=placed, note=note,
    )


# --- Scwrl4 wrapper ------------------------------------------------------

def scwrl4_binary() -> str | None:
    return shutil.which("Scwrl4") or shutil.which("scwrl4")


def scwrl4(
    structure: Structure,
    pdb_path: str | Path,
    chain_id: str,
    resseq: int,
    target_resname: str,
) -> MutationResult:
    """Mutate one residue via an external Scwrl4 binary.

    Scwrl4 repacks side chains for the whole input structure it's given;
    to keep `prot mutate` surgical (only the target residue changes), this
    writes a copy of the structure with the target residue's side chain
    stripped back to its backbone and renamed, runs Scwrl4 on that copy,
    and then grafts *only* the rebuilt residue's atoms back onto an
    otherwise-untouched copy of the original structure -- every other
    residue keeps its original coordinates exactly.
    """
    binary = scwrl4_binary()
    if binary is None:
        raise Scwrl4NotAvailableError(
            "Scwrl4 executable not found on PATH. Install Scwrl4 (requires a "
            "license from the Dunbrack lab, UNC) to build a complete, "
            "energetically-optimized side chain; the built-in fallback "
            "(--method cb_only) places backbone + C-beta only."
        )

    chain, residue = _find_residue(structure, chain_id, resseq)
    original_resname = residue.resname.strip().upper()
    required = ("N", "CA", "C", "O")
    missing = [name for name in required if name not in residue]
    if missing:
        raise MutationError(
            f"Cannot mutate {chain_id}/{original_resname}{resseq}: missing backbone atom(s) {missing}"
        )

    stripped = pio.load_structure(pdb_path, structure_id="scwrl4_input")
    _, stripped_residue = _find_residue(stripped, chain_id, resseq)
    for atom in list(stripped_residue):
        if atom.get_name() not in required:
            stripped_residue.detach_child(atom.get_id())
    stripped_residue.resname = target_resname

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_in = Path(tmpdir) / "input.pdb"
        tmp_out = Path(tmpdir) / "output.pdb"
        pio.save_structure(stripped, tmp_in, fmt="pdb")

        result = subprocess.run(
            [binary, "-i", str(tmp_in), "-o", str(tmp_out)],
            capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0 or not tmp_out.exists():
            raise RuntimeError(
                f"Scwrl4 failed (exit {result.returncode}): {result.stderr or result.stdout}"
            )

        repacked = pio.load_structure(tmp_out, structure_id="scwrl4_output")

    _, repacked_residue = _find_residue(repacked, chain_id, resseq)
    placed = [atom.get_name() for atom in repacked_residue]

    chain.detach_child(residue.id)
    repacked_residue.detach_parent()
    chain.add(repacked_residue)

    return MutationResult(
        chain_id=chain_id, resseq=resseq, original_resname=original_resname,
        new_resname=target_resname, method="scwrl4", atoms_placed=placed,
        note="Full side chain built by Scwrl4 (dead-end elimination over its rotamer library).",
    )


# --- Unified entry point ------------------------------------------------

def mutate_residue(
    structure: Structure,
    pdb_path: str | Path,
    chain_id: str,
    resseq: int,
    target: str,
    method: str = "auto",
) -> MutationResult:
    """Apply a point mutation in place (the passed-in `structure` object is
    modified) and return a MutationResult describing what was done.

    method="auto" (default): use Scwrl4 if installed, otherwise cb_only.
    method="scwrl4": require Scwrl4, raise Scwrl4NotAvailableError if missing.
    method="cb_only": always use the dependency-free backbone+C-beta fallback.
    """
    target_resname = normalize_resname(target)
    if target_resname not in PROTEIN_RESIDUES:
        raise MutationError(f"{target_resname} is not a standard amino acid")

    if method == "cb_only":
        return cb_only(structure, chain_id, resseq, target_resname)

    if method == "scwrl4":
        return scwrl4(structure, pdb_path, chain_id, resseq, target_resname)

    if method == "auto":
        if scwrl4_binary() is not None:
            return scwrl4(structure, pdb_path, chain_id, resseq, target_resname)
        return cb_only(structure, chain_id, resseq, target_resname)

    raise ValueError(f"Unknown method: {method!r} (expected 'auto', 'scwrl4', or 'cb_only')")
