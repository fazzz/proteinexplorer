"""Pocket / cavity analysis (spec section "Pocket Analysis").

No external pocket-detection tool (fpocket, P2Rank, ...) is available in
this environment, so this module implements a simplified, dependency-free
grid-based cavity detector in the spirit of LIGSITE:

1. Lay a 3D grid over the region of interest.
2. Mark grid points that fall inside any atom's van der Waals radius (+ a
   solvent probe) as "occupied"; everything else is "empty".
3. For each empty point, cast rays along 7 fixed axis directions (each
   checked in both the + and - sense). An axis counts as "enclosed" if
   protein is hit on *both* sides within a search distance. Points with
   enough enclosed axes are cavity/pocket points (LIGSITE's core idea).
4. Cluster adjacent pocket points into pockets (grid connected components).
5. For each pocket, report volume, an approximate surface area, the
   lining residues, their hydrophobic fraction, and a simple heuristic
   druggability score.

This is a real, working approximation -- not a reimplementation of
fpocket's actual (more sophisticated, alpha-sphere-based) algorithm, and
the druggability score is a hand-built heuristic, not a trained model.
Both are called out as such in the docstrings/CLI output.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Structure import Structure

from proteinexplorer.descriptor import HYDROPHOBIC_RESIDUES
from proteinexplorer.models import ResidueCategory, classify_residue

VDW_RADII: dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "SE": 1.90, "ZN": 1.39, "MG": 1.73, "CA": 1.97,
    "FE": 1.56, "NA": 2.27, "K": 2.75, "CL": 1.75,
}
_DEFAULT_VDW = 1.70
PROBE_RADIUS = 1.4  # water probe

# 7 fixed axis directions (unit vectors); each is checked in both the +
# and - sense, so this covers the classic LIGSITE 7-direction scheme.
_RAW_AXES = [
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1),
]
AXIS_DIRECTIONS = [np.array(v, dtype=float) / np.linalg.norm(v) for v in _RAW_AXES]


def _vdw(atom) -> float:
    return VDW_RADII.get((atom.element or "").strip().upper(), _DEFAULT_VDW)


def _residue_label(residue) -> str:
    chain_id = residue.get_parent().id
    return f"{chain_id}/{residue.resname}{residue.id[1]}"


class PocketGrid:
    """Occupancy/enclosure testing against a fixed atom set, reused across
    every grid point so the NeighborSearch tree is only built once."""

    def __init__(self, atoms: list, max_ray_length: float = 6.0, ray_step: float = 1.5):
        self.atoms = atoms
        self.max_ray_length = max_ray_length
        self.ray_step = ray_step
        self._search = NeighborSearch(atoms) if atoms else None
        self._max_vdw = max((_vdw(a) for a in atoms), default=_DEFAULT_VDW)

    def is_occupied(self, point: np.ndarray) -> bool:
        if self._search is None:
            return False
        nearby = self._search.search(point, self._max_vdw + PROBE_RADIUS, level="A")
        for atom in nearby:
            if np.linalg.norm(atom.coord - point) <= _vdw(atom) + PROBE_RADIUS:
                return True
        return False

    def enclosed_axis_count(self, point: np.ndarray) -> int:
        count = 0
        steps = np.arange(self.ray_step, self.max_ray_length + 1e-6, self.ray_step)
        for axis in AXIS_DIRECTIONS:
            pos_hit = any(self.is_occupied(point + axis * d) for d in steps)
            neg_hit = any(self.is_occupied(point - axis * d) for d in steps)
            if pos_hit and neg_hit:
                count += 1
        return count


@dataclass
class Pocket:
    id: int
    n_grid_points: int
    volume: float  # Angstrom^3
    surface_area: float  # Angstrom^2, rough estimate
    centroid: np.ndarray
    residues: list[str]
    hydrophobicity: float
    druggability_score: float


def _grid_indices_in_box(low: np.ndarray, high: np.ndarray, spacing: float):
    n = np.maximum(np.round((high - low) / spacing).astype(int), 1)
    return n


def _protein_like_atoms(structure: Structure):
    """Everything except water: the receptor's own excluded volume for
    cavity detection (protein, nucleic acid, ions, and any bound ligand)."""
    model = next(iter(structure))
    return [
        a for a in model.get_atoms()
        if classify_residue(a.get_parent().resname, a.get_parent().id[0]) is not ResidueCategory.WATER
    ]


def find_pockets(
    structure: Structure,
    atoms: list | None = None,
    spacing: float = 1.5,
    padding: float = 3.0,
    max_ray_length: float = 6.0,
    ray_step: float = 1.5,
    min_enclosed_axes: int = 5,
    min_pocket_points: int = 3,
    lining_distance: float = 4.5,
    max_grid_points: int = 200_000,
) -> list[Pocket]:
    """Detect cavities via the LIGSITE-style grid method described in the
    module docstring.

    `atoms` restricts which atoms bound the search grid (e.g. a selection
    around a region of interest); it defaults to every non-water atom in
    the structure. All atoms (regardless of `atoms`) are still used to
    determine occupancy/enclosure, so a restricted region is still judged
    against its real surroundings.
    """
    all_atoms = _protein_like_atoms(structure)
    if not all_atoms:
        return []

    region_atoms = atoms if atoms is not None else all_atoms
    coords = np.array([a.coord for a in region_atoms])
    low = coords.min(axis=0) - padding
    high = coords.max(axis=0) + padding

    grid = PocketGrid(all_atoms, max_ray_length=max_ray_length, ray_step=ray_step)

    counts = np.maximum(np.round((high - low) / spacing).astype(int), 1)
    total_points = int(np.prod(counts + 1))
    if total_points > max_grid_points:
        raise ValueError(
            f"Grid would have {total_points} points (limit {max_grid_points}); "
            f"restrict the search with `atoms=`/--selection or increase `spacing`."
        )

    pocket_point_indices: dict[tuple[int, int, int], np.ndarray] = {}
    for i in range(counts[0] + 1):
        for j in range(counts[1] + 1):
            for k in range(counts[2] + 1):
                point = low + np.array([i, j, k]) * spacing
                if grid.is_occupied(point):
                    continue
                if grid.enclosed_axis_count(point) >= min_enclosed_axes:
                    pocket_point_indices[(i, j, k)] = point

    clusters = _cluster_grid_points(pocket_point_indices)

    all_ns = NeighborSearch(all_atoms)
    pockets: list[Pocket] = []
    for pocket_id, index_group in enumerate(clusters, start=1):
        if len(index_group) < min_pocket_points:
            continue
        points = np.array([pocket_point_indices[idx] for idx in index_group])
        volume = len(points) * (spacing ** 3)
        surface_area = _estimate_surface_area(index_group, spacing)
        centroid = points.mean(axis=0)

        lining_residues = set()
        for point in points:
            for atom in all_ns.search(point, lining_distance, level="A"):
                residue = atom.get_parent()
                category = classify_residue(residue.resname, residue.id[0])
                if category is ResidueCategory.PROTEIN:
                    lining_residues.add(residue)

        residue_labels = sorted(_residue_label(r) for r in lining_residues)
        if lining_residues:
            hydrophobic_count = sum(
                1 for r in lining_residues if r.resname.strip().upper() in HYDROPHOBIC_RESIDUES
            )
            hydrophobicity = hydrophobic_count / len(lining_residues)
        else:
            hydrophobicity = 0.0

        pockets.append(
            Pocket(
                id=pocket_id,
                n_grid_points=len(points),
                volume=volume,
                surface_area=surface_area,
                centroid=centroid,
                residues=residue_labels,
                hydrophobicity=hydrophobicity,
                druggability_score=_druggability_score(volume, hydrophobicity),
            )
        )

    pockets.sort(key=lambda p: p.volume, reverse=True)
    for new_id, p in enumerate(pockets, start=1):
        p.id = new_id
    return pockets


def _cluster_grid_points(pocket_point_indices: dict[tuple[int, int, int], np.ndarray]) -> list[list[tuple[int, int, int]]]:
    """6-connectivity connected components over the pocket grid indices."""
    remaining = set(pocket_point_indices.keys())
    clusters: list[list[tuple[int, int, int]]] = []
    neighbors6 = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

    while remaining:
        start = next(iter(remaining))
        stack = [start]
        remaining.discard(start)
        cluster = [start]
        while stack:
            cur = stack.pop()
            for dx, dy, dz in neighbors6:
                nb = (cur[0] + dx, cur[1] + dy, cur[2] + dz)
                if nb in remaining:
                    remaining.discard(nb)
                    stack.append(nb)
                    cluster.append(nb)
        clusters.append(cluster)
    return clusters


def _estimate_surface_area(index_group: list[tuple[int, int, int]], spacing: float) -> float:
    """Rough surface-area estimate: count grid-point faces that border a
    non-pocket cell, times the area of one grid face."""
    members = set(index_group)
    neighbors6 = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    exposed_faces = 0
    for idx in index_group:
        for dx, dy, dz in neighbors6:
            nb = (idx[0] + dx, idx[1] + dy, idx[2] + dz)
            if nb not in members:
                exposed_faces += 1
    return exposed_faces * (spacing ** 2)


def _druggability_score(volume: float, hydrophobicity: float) -> float:
    """Simple 0-1 heuristic (NOT a trained/validated model): rewards
    volumes in the typical small-molecule-druggable range (peaking around
    600 A^3) and higher hydrophobic lining, since real binding pockets are
    usually moderately sized, mostly-enclosed hydrophobic cavities."""
    ideal_volume = 600.0
    volume_score = float(np.exp(-((volume - ideal_volume) / 900.0) ** 2))
    return float(np.clip(0.6 * volume_score + 0.4 * hydrophobicity, 0.0, 1.0))
