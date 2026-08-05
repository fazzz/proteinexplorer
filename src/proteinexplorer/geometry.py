"""Geometry analysis (spec section "Geometry").

All functions take already-resolved atom lists (from selection.select) or
raw coordinates, and stay independent of the CLI layer so they're directly
usable from `prot geometry` commands as well as from later commands
(contact, pocket, compare, ...) that need the same primitives.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from proteinexplorer.descriptor import ATOMIC_WEIGHTS


def _coords(atoms) -> np.ndarray:
    return np.array([atom.coord for atom in atoms], dtype=float)


def _masses(atoms) -> np.ndarray:
    return np.array(
        [ATOMIC_WEIGHTS.get((a.element or "").strip().upper(), 0.0) for a in atoms],
        dtype=float,
    )


class GeometryError(ValueError):
    pass


def centroid(atoms) -> np.ndarray:
    """Unweighted geometric center (a.k.a. center of geometry)."""
    if not atoms:
        raise GeometryError("Selection is empty")
    return _coords(atoms).mean(axis=0)


def center_of_mass(atoms) -> np.ndarray:
    coords = _coords(atoms)
    masses = _masses(atoms)
    if not atoms:
        raise GeometryError("Selection is empty")
    total_mass = masses.sum()
    if total_mass == 0:
        raise GeometryError("Selection has no recognized element masses")
    return (coords * masses[:, None]).sum(axis=0) / total_mass


# --- Distance / angle / dihedral -------------------------------------------

def distance(atoms_a, atoms_b) -> float:
    """Distance between two selections. If both selections are single
    atoms this is the atom-atom distance; otherwise it is the distance
    between the two selections' centroids (covers atom-residue,
    residue-residue, and chain-chain distances with one implementation)."""
    point_a = atoms_a[0].coord if len(atoms_a) == 1 else centroid(atoms_a)
    point_b = atoms_b[0].coord if len(atoms_b) == 1 else centroid(atoms_b)
    return float(np.linalg.norm(np.asarray(point_a) - np.asarray(point_b)))


def _point_of(atoms) -> np.ndarray:
    return np.asarray(atoms[0].coord if len(atoms) == 1 else centroid(atoms))


def angle(atoms_a, atoms_b, atoms_c) -> float:
    """Angle (degrees) at vertex B formed by B->A and B->C. Works for a
    3-atom bond angle as well as an arbitrary 3-point angle between
    selection centroids."""
    a, b, c = _point_of(atoms_a), _point_of(atoms_b), _point_of(atoms_c)
    v1 = a - b
    v2 = c - b
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_theta)))


def dihedral(atoms_a, atoms_b, atoms_c, atoms_d) -> float:
    """Torsion angle (degrees, -180..180) defined by four points/selections
    A-B-C-D. Covers arbitrary 4-point torsions as well as phi/psi/omega/chi
    once the relevant atoms have been picked out."""
    p0, p1, p2, p3 = (
        _point_of(atoms_a), _point_of(atoms_b), _point_of(atoms_c), _point_of(atoms_d)
    )
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1_norm = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1_norm) * b1_norm
    w = b2 - np.dot(b2, b1_norm) * b1_norm

    x = np.dot(v, w)
    y = np.dot(np.cross(b1_norm, v), w)
    return float(np.degrees(np.arctan2(y, x)))


# --- Backbone / sidechain torsions ------------------------------------------

# Side-chain chi-angle atom definitions (standard 4-atom sets per residue).
# Only residues with at least one chi angle are listed.
CHI_ATOMS: dict[str, list[tuple[str, str, str, str]]] = {
    "ARG": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"),
            ("CB", "CG", "CD", "NE"), ("CG", "CD", "NE", "CZ")],
    "ASN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "ASP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "CYS": [("N", "CA", "CB", "SG")],
    "GLN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"),
            ("CB", "CG", "CD", "OE1")],
    "GLU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"),
            ("CB", "CG", "CD", "OE1")],
    "HIS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")],
    "ILE": [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")],
    "LEU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "LYS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"),
            ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")],
    "MET": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD"),
            ("CB", "CG", "SD", "CE")],
    "PHE": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "PRO": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "SER": [("N", "CA", "CB", "OG")],
    "THR": [("N", "CA", "CB", "OG1")],
    "TRP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "TYR": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "VAL": [("N", "CA", "CB", "CG1")],
}


@dataclass
class ResidueTorsions:
    chain_id: str
    resseq: int
    resname: str
    phi: float | None
    psi: float | None
    omega: float | None
    chi: list[float]


def _protein_residue_list(chain):
    from proteinexplorer.models import ResidueCategory, classify_residue
    return [
        r for r in chain
        if classify_residue(r.resname, r.id[0]) is ResidueCategory.PROTEIN
    ]


def backbone_torsions(chain, resseq: int) -> ResidueTorsions:
    """phi/psi/omega for one residue, using its sequence neighbors in the
    same chain (None where a neighbor is missing, e.g. chain termini)."""
    residues = _protein_residue_list(chain)
    index = next((i for i, r in enumerate(residues) if r.id[1] == resseq), None)
    if index is None:
        raise GeometryError(f"No protein residue with resid {resseq} in chain {chain.id}")

    residue = residues[index]
    prev_res = residues[index - 1] if index > 0 else None
    next_res = residues[index + 1] if index < len(residues) - 1 else None

    def has(res, name):
        return res is not None and name in res

    phi = psi = omega = None
    if has(prev_res, "C") and all(a in residue for a in ("N", "CA", "C")):
        phi = dihedral([prev_res["C"]], [residue["N"]], [residue["CA"]], [residue["C"]])
    if has(next_res, "N") and all(a in residue for a in ("N", "CA", "C")):
        psi = dihedral([residue["N"]], [residue["CA"]], [residue["C"]], [next_res["N"]])
    if has(next_res, "CA") and all(a in residue for a in ("CA", "C")) and has(next_res, "N"):
        omega = dihedral([residue["CA"]], [residue["C"]], [next_res["N"]], [next_res["CA"]])

    chi_defs = CHI_ATOMS.get(residue.resname.strip().upper(), [])
    chi_values = []
    for a1, a2, a3, a4 in chi_defs:
        if all(name in residue for name in (a1, a2, a3, a4)):
            chi_values.append(
                dihedral([residue[a1]], [residue[a2]], [residue[a3]], [residue[a4]])
            )

    return ResidueTorsions(
        chain_id=chain.id, resseq=resseq, resname=residue.resname,
        phi=phi, psi=psi, omega=omega, chi=chi_values,
    )


# --- Coordinate analysis -----------------------------------------------------

@dataclass
class BoundingBox:
    min: np.ndarray
    max: np.ndarray

    @property
    def size(self) -> np.ndarray:
        return self.max - self.min


def bounding_box(atoms) -> BoundingBox:
    coords = _coords(atoms)
    if coords.size == 0:
        raise GeometryError("Selection is empty")
    return BoundingBox(min=coords.min(axis=0), max=coords.max(axis=0))


@dataclass
class PlaneFit:
    point: np.ndarray
    normal: np.ndarray
    rms_deviation: float


def fit_plane(atoms) -> PlaneFit:
    """Least-squares plane fit via SVD. The normal is the singular vector
    with the smallest singular value; rms_deviation is the RMS distance of
    the points from the fitted plane."""
    coords = _coords(atoms)
    if len(coords) < 3:
        raise GeometryError("Plane fitting needs at least 3 atoms")
    point = coords.mean(axis=0)
    centered = coords - point
    _, singular_values, vt = np.linalg.svd(centered)
    normal = vt[-1]
    normal = normal / np.linalg.norm(normal)
    deviations = centered @ normal
    rms = float(np.sqrt(np.mean(deviations ** 2)))
    return PlaneFit(point=point, normal=normal, rms_deviation=rms)


@dataclass
class PrincipalAxes:
    eigenvalues: np.ndarray  # descending order
    eigenvectors: np.ndarray  # columns are the axes, matching eigenvalues order


def principal_axes(atoms, mass_weighted: bool = True) -> PrincipalAxes:
    coords = _coords(atoms)
    if len(coords) < 2:
        raise GeometryError("Principal axes need at least 2 atoms")
    weights = _masses(atoms) if mass_weighted else np.ones(len(atoms))
    if mass_weighted and weights.sum() == 0:
        weights = np.ones(len(atoms))
    center = (coords * weights[:, None]).sum(axis=0) / weights.sum()
    centered = coords - center
    cov = (centered * weights[:, None]).T @ centered / weights.sum()
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = np.argsort(eigenvalues)[::-1]
    return PrincipalAxes(eigenvalues=eigenvalues[order], eigenvectors=eigenvectors[:, order])


@dataclass
class MomentOfInertia:
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray


def moment_of_inertia(atoms) -> MomentOfInertia:
    """Mass-weighted moment-of-inertia tensor eigen-decomposition about the
    center of mass."""
    coords = _coords(atoms)
    masses = _masses(atoms)
    if len(coords) < 2:
        raise GeometryError("Moment of inertia needs at least 2 atoms")
    if masses.sum() == 0:
        raise GeometryError("Selection has no recognized element masses")
    com = (coords * masses[:, None]).sum(axis=0) / masses.sum()
    centered = coords - com

    tensor = np.zeros((3, 3))
    for pos, mass in zip(centered, masses):
        x, y, z = pos
        tensor[0, 0] += mass * (y ** 2 + z ** 2)
        tensor[1, 1] += mass * (x ** 2 + z ** 2)
        tensor[2, 2] += mass * (x ** 2 + y ** 2)
        tensor[0, 1] -= mass * x * y
        tensor[0, 2] -= mass * x * z
        tensor[1, 2] -= mass * y * z
    tensor[1, 0] = tensor[0, 1]
    tensor[2, 0] = tensor[0, 2]
    tensor[2, 1] = tensor[1, 2]

    eigenvalues, eigenvectors = np.linalg.eigh(tensor)
    order = np.argsort(eigenvalues)
    return MomentOfInertia(eigenvalues=eigenvalues[order], eigenvectors=eigenvectors[:, order])


def radius_of_gyration(atoms) -> float:
    coords = _coords(atoms)
    masses = _masses(atoms)
    if masses.sum() == 0:
        raise GeometryError("Selection has no recognized element masses")
    com = (coords * masses[:, None]).sum(axis=0) / masses.sum()
    sq_dev = ((coords - com) ** 2).sum(axis=1)
    return float(np.sqrt((masses * sq_dev).sum() / masses.sum()))


def distance_matrix(groups: list) -> np.ndarray:
    """Pairwise centroid distance matrix for a list of atom-selections
    (e.g. one entry per residue) -- covers atom/residue/chain distance
    matrices depending on what each group contains."""
    if len(groups) < 2:
        raise GeometryError("Need at least 2 selections for a distance matrix")
    points = np.array([_point_of(g) for g in groups])
    diff = points[:, None, :] - points[None, :, :]
    return np.sqrt((diff ** 2).sum(axis=-1))


def rmsd(atoms_a, atoms_b, fit: bool = True) -> float:
    """RMSD between two atom lists of equal length and matching order.

    fit=True (default) superposes atoms_b onto atoms_a first (Kabsch, via
    Bio.PDB.Superimposer) before computing RMSD -- appropriate for
    comparing two independently-solved structures/conformations.
    fit=False computes the raw coordinate RMSD with no superposition --
    appropriate when both selections already share the same reference
    frame (e.g. two chains of one asymmetric unit already related by a
    crystallographic operation applied upstream).
    """
    if len(atoms_a) != len(atoms_b):
        raise GeometryError(
            f"Selections must have the same number of atoms for RMSD "
            f"({len(atoms_a)} vs {len(atoms_b)})"
        )
    if not atoms_a:
        raise GeometryError("Selections are empty")

    if not fit:
        coords_a = _coords(atoms_a)
        coords_b = _coords(atoms_b)
        return float(np.sqrt(np.mean(((coords_a - coords_b) ** 2).sum(axis=1))))

    from Bio.PDB.Superimposer import Superimposer

    sup = Superimposer()
    sup.set_atoms(list(atoms_a), list(atoms_b))
    return float(sup.rms)
