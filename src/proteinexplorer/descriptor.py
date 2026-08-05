"""Structure descriptors (spec section "Descriptor").

Everything here works with Biopython + numpy only, no external binaries.
Secondary structure composition delegates to secondary.py, which prefers
an external DSSP binary and falls back to a dependency-free phi/psi-based
classifier when DSSP isn't installed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure

from proteinexplorer.models import ResidueCategory, classify_residue

# Approximate atomic weights (g/mol) for elements commonly seen in PDB files.
# Hydrogens are usually absent from crystal structures, so this yields a
# heavy-atom molecular weight unless the input already has hydrogens added.
ATOMIC_WEIGHTS: dict[str, float] = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06,
    "P": 30.974, "SE": 78.971, "NA": 22.990, "K": 39.098, "CL": 35.45,
    "CA": 40.078, "MG": 24.305, "ZN": 65.38, "MN": 54.938, "FE": 55.845,
    "CU": 63.546, "NI": 58.693, "CO": 58.933, "CD": 112.414, "BR": 79.904,
}

# Kyte-Doolittle-style hydrophobic residue set (positive hydropathy index).
HYDROPHOBIC_RESIDUES: frozenset[str] = frozenset(
    {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "CYS", "PRO", "GLY"}
)

_DISULFIDE_SG_DISTANCE = 2.5  # Angstrom, generous cutoff around the ~2.05 A S-S bond
_CONTACT_CA_DISTANCE = 8.0  # Angstrom, standard CA-CA contact-map cutoff


@dataclass
class StructureDescriptors:
    molecular_weight: float
    n_atoms: int
    n_residues: int
    n_chains: int
    n_ligands: int
    n_waters: int
    sasa_total: float | None
    radius_of_gyration: float | None
    contact_density: float | None
    hydrophobic_ratio: float | None
    disulfide_count: int
    secondary_structure: dict[str, float] | None
    secondary_structure_method: str | None
    secondary_structure_error: str | None


def _protein_residues(structure: Structure):
    model = next(iter(structure))
    for chain in model:
        for residue in chain:
            if classify_residue(residue.resname, residue.id[0]) is ResidueCategory.PROTEIN:
                yield residue


def molecular_weight(structure: Structure) -> float:
    """Sum of atomic weights over all atoms present in the structure
    (heavy-atom MW if the input has no hydrogens, as is typical for
    crystallographic PDB/mmCIF files)."""
    total = 0.0
    model = next(iter(structure))
    for atom in model.get_atoms():
        element = (atom.element or "").strip().upper()
        total += ATOMIC_WEIGHTS.get(element, 0.0)
    return total


def sasa_total(structure: Structure) -> float:
    """Total solvent-accessible surface area (Angstrom^2) via the
    Shrake-Rupley algorithm built into Biopython (no external tool)."""
    sr = ShrakeRupley()
    sr.compute(structure, level="A")
    model = next(iter(structure))
    return float(sum(atom.sasa for atom in model.get_atoms() if atom.sasa is not None))


def radius_of_gyration(structure: Structure) -> float | None:
    """Mass-weighted radius of gyration over all atoms in the first model."""
    model = next(iter(structure))
    coords: list[np.ndarray] = []
    masses: list[float] = []
    for atom in model.get_atoms():
        element = (atom.element or "").strip().upper()
        mass = ATOMIC_WEIGHTS.get(element)
        if mass is None:
            continue
        coords.append(atom.coord)
        masses.append(mass)
    if not coords:
        return None

    coords_arr = np.array(coords)
    masses_arr = np.array(masses)
    total_mass = masses_arr.sum()
    com = (coords_arr * masses_arr[:, None]).sum(axis=0) / total_mass
    sq_dev = ((coords_arr - com) ** 2).sum(axis=1)
    rg_sq = (masses_arr * sq_dev).sum() / total_mass
    return float(np.sqrt(rg_sq))


def contact_density(structure: Structure, cutoff: float = _CONTACT_CA_DISTANCE) -> float | None:
    """Number of residue-residue CA-CA contacts within `cutoff` Angstrom,
    normalized by residue count. Only protein residues with a CA atom are
    considered."""
    ca_coords = []
    for residue in _protein_residues(structure):
        if "CA" in residue:
            ca_coords.append(residue["CA"].coord)
    n = len(ca_coords)
    if n < 2:
        return None
    coords = np.array(ca_coords)
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    n_contacts = int(np.sum((dist < cutoff) & (dist > 0)) / 2)
    return n_contacts / n


def hydrophobic_ratio(structure: Structure) -> float | None:
    """Fraction of protein residues that are hydrophobic (Kyte-Doolittle
    positive-hydropathy set)."""
    total = 0
    hydrophobic = 0
    for residue in _protein_residues(structure):
        total += 1
        if residue.resname.strip().upper() in HYDROPHOBIC_RESIDUES:
            hydrophobic += 1
    if total == 0:
        return None
    return hydrophobic / total


def disulfide_count(structure: Structure, cutoff: float = _DISULFIDE_SG_DISTANCE) -> int:
    """Count CYS-CYS pairs whose SG atoms are within `cutoff` Angstrom
    (each pair counted once)."""
    sg_atoms = []
    for residue in _protein_residues(structure):
        if residue.resname.strip().upper() == "CYS" and "SG" in residue:
            sg_atoms.append(residue["SG"].coord)
    n = len(sg_atoms)
    if n < 2:
        return 0
    coords = np.array(sg_atoms)
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    count = int(np.sum((dist < cutoff) & (dist > 0)) / 2)
    return count


class DSSPNotAvailableError(RuntimeError):
    pass


def compute_descriptors(
    structure: Structure,
    pdb_path: str | Path,
    category_totals: dict[str, int],
) -> StructureDescriptors:
    from proteinexplorer import secondary as sec

    model = next(iter(structure))
    n_atoms = sum(1 for _ in model.get_atoms())
    n_residues = sum(1 for _ in model.get_residues())
    n_chains = sum(1 for _ in model)

    try:
        sasa = sasa_total(structure)
    except Exception:
        sasa = None

    ss_composition: dict[str, float] | None = None
    ss_method: str | None = None
    ss_error: str | None = None
    try:
        residues, ss_method = sec.secondary_structure(structure, pdb_path=pdb_path, method="auto")
        ss_composition = sec.composition(residues)
    except Exception as exc:  # pragma: no cover - defensive
        ss_error = f"secondary structure assignment failed: {exc}"

    return StructureDescriptors(
        molecular_weight=molecular_weight(structure),
        n_atoms=n_atoms,
        n_residues=n_residues,
        n_chains=n_chains,
        n_ligands=category_totals.get(ResidueCategory.LIGAND.value, 0),
        n_waters=category_totals.get(ResidueCategory.WATER.value, 0),
        sasa_total=sasa,
        radius_of_gyration=radius_of_gyration(structure),
        contact_density=contact_density(structure),
        hydrophobic_ratio=hydrophobic_ratio(structure),
        disulfide_count=disulfide_count(structure),
        secondary_structure=ss_composition,
        secondary_structure_method=ss_method,
        secondary_structure_error=ss_error,
    )
