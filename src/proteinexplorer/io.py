"""Structure file I/O: PDB / mmCIF read & write, and summarization.

Thin wrapper around Bio.PDB's parsers/writers. Format is inferred from the
file extension unless explicitly given.
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import Literal

from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser
from Bio.PDB.Structure import Structure

from proteinexplorer.models import ChainSummary, StructureSummary, classify_residue

StructureFormat = Literal["pdb", "mmcif"]

_PDB_SUFFIXES = {".pdb", ".ent", ".pdb1"}
_MMCIF_SUFFIXES = {".cif", ".mmcif"}


class UnknownFormatError(ValueError):
    pass


def infer_format(path: str | Path) -> StructureFormat:
    p = Path(path)
    suffixes = p.suffixes
    # handle .pdb.gz / .cif.gz
    stem_suffix = suffixes[-2] if len(suffixes) >= 2 and suffixes[-1] == ".gz" else (
        suffixes[-1] if suffixes else ""
    )
    if stem_suffix in _PDB_SUFFIXES:
        return "pdb"
    if stem_suffix in _MMCIF_SUFFIXES:
        return "mmcif"
    raise UnknownFormatError(
        f"Cannot infer structure format from filename: {p.name} "
        f"(expected one of {sorted(_PDB_SUFFIXES | _MMCIF_SUFFIXES)}, optionally .gz)"
    )


def load_structure(
    path: str | Path,
    structure_id: str = "structure",
    fmt: StructureFormat | None = None,
) -> Structure:
    """Parse a PDB or mmCIF file (optionally gzip-compressed) into a
    Bio.PDB.Structure.Structure object."""
    p = Path(path)
    fmt = fmt or infer_format(p)

    if p.suffix == ".gz":
        opener = gzip.open
    else:
        opener = open

    parser = PDBParser(QUIET=True) if fmt == "pdb" else MMCIFParser(QUIET=True)
    with opener(p, "rt") as handle:  # type: ignore[arg-type]
        structure = parser.get_structure(structure_id, handle)
    return structure


def save_structure(
    structure: Structure,
    path: str | Path,
    fmt: StructureFormat | None = None,
) -> None:
    """Write a Bio.PDB structure to a PDB or mmCIF file."""
    p = Path(path)
    fmt = fmt or infer_format(p)

    io_writer = PDBIO() if fmt == "pdb" else MMCIFIO()
    io_writer.set_structure(structure)
    io_writer.save(str(p))


def copy_raw(src: str | Path, dst: str | Path) -> None:
    """Copy the original structure file byte-for-byte (used on import, so
    the project always retains an untouched copy of what was imported)."""
    shutil.copyfile(src, dst)


def header_info(structure) -> dict:
    """Extract selected fields from Bio.PDB's parsed header dict (PDB
    REMARK/HEADER records or mmCIF equivalents), when present."""
    header = getattr(structure, "header", {}) or {}
    method = header.get("structure_method") or None
    if method == "unknown":
        method = None
    return {
        "name": header.get("name") or None,
        "structure_method": method,
        "resolution": header.get("resolution"),
        "deposition_date": header.get("deposition_date") or None,
    }



def summarize(structure: Structure) -> StructureSummary:
    """Build a StructureSummary without needing to keep the full object
    around (used for `prot status` / `prot info`)."""
    n_models = len(structure)
    model = next(iter(structure), None)

    chains: list[ChainSummary] = []
    n_atoms = 0
    hetero_resnames: set[str] = set()
    has_altloc = False

    if model is not None:
        for chain in model:
            categories: dict[str, int] = {}
            n_residues = 0
            for residue in chain:
                hetero_flag = residue.id[0]
                category = classify_residue(residue.resname, hetero_flag)
                categories[category.value] = categories.get(category.value, 0) + 1
                n_residues += 1
                if hetero_flag not in (" ", "W"):
                    hetero_resnames.add(residue.resname.strip())
                for atom in residue:
                    n_atoms += 1
                    if atom.is_disordered():
                        has_altloc = True
            chains.append(
                ChainSummary(chain_id=chain.id, n_residues=n_residues, categories=categories)
            )

    return StructureSummary(
        n_models=n_models,
        chains=chains,
        n_atoms=n_atoms,
        hetero_resnames=sorted(hetero_resnames),
        has_altloc=has_altloc,
    )
