from pathlib import Path

import numpy as np
import pytest
from Bio.PDB.Atom import Atom

from proteinexplorer import io as pio
from proteinexplorer import pocket as pk


def _atom(name, coord, element="C"):
    return Atom(name, np.array(coord, dtype=float), 20.0, 1.0, " ", name, 1, element=element)


# --- PocketGrid primitives -----------------------------------------------

def test_is_occupied_near_and_far():
    grid = pk.PocketGrid([_atom("C1", (0, 0, 0))])
    assert grid.is_occupied(np.array([1.0, 0.0, 0.0]))
    assert not grid.is_occupied(np.array([5.0, 0.0, 0.0]))


def test_is_occupied_empty_atom_list():
    grid = pk.PocketGrid([])
    assert not grid.is_occupied(np.array([0.0, 0.0, 0.0]))


def test_enclosed_axis_count_full_cage_is_seven():
    walls = [_atom("W", axis * sign * 4.0) for axis in pk.AXIS_DIRECTIONS for sign in (1, -1)]
    grid = pk.PocketGrid(walls)
    assert grid.enclosed_axis_count(np.array([0.0, 0.0, 0.0])) == 7
    assert not grid.is_occupied(np.array([0.0, 0.0, 0.0]))


def test_enclosed_axis_count_no_walls_is_zero():
    grid = pk.PocketGrid([_atom("W", (100.0, 100.0, 100.0))])
    assert grid.enclosed_axis_count(np.array([0.0, 0.0, 0.0])) == 0


def test_enclosed_axis_count_partial_cage_covers_at_least_primary_axes():
    primary_only = [
        _atom("W", np.array(v, dtype=float) * 4.0)
        for v in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    ]
    grid = pk.PocketGrid(primary_only)
    count = grid.enclosed_axis_count(np.array([0.0, 0.0, 0.0]))
    assert 3 <= count <= 7


# --- Clustering ------------------------------------------------------------

def test_cluster_grid_points_single_component():
    points = {(0, 0, 0): None, (1, 0, 0): None, (1, 1, 0): None}
    clusters = pk._cluster_grid_points(points)
    assert len(clusters) == 1
    assert set(clusters[0]) == set(points.keys())


def test_cluster_grid_points_two_components():
    points = {(0, 0, 0): None, (1, 0, 0): None, (50, 50, 50): None}
    clusters = pk._cluster_grid_points(points)
    assert len(clusters) == 2


def test_cluster_grid_points_empty():
    assert pk._cluster_grid_points({}) == []


# --- Surface area / druggability heuristics ---------------------------

def test_surface_area_single_point_has_six_exposed_faces():
    area = pk._estimate_surface_area([(0, 0, 0)], spacing=1.5)
    assert area == pytest.approx(6 * 1.5 ** 2)


def test_surface_area_two_adjacent_points_have_ten_exposed_faces():
    # a 1x1x2 block has 2*(1+1+2)=10 exposed unit faces
    area = pk._estimate_surface_area([(0, 0, 0), (0, 0, 1)], spacing=1.0)
    assert area == pytest.approx(10.0)


def test_druggability_score_peaks_near_ideal_volume():
    near_ideal = pk._druggability_score(600.0, hydrophobicity=1.0)
    tiny = pk._druggability_score(5.0, hydrophobicity=1.0)
    huge = pk._druggability_score(20000.0, hydrophobicity=1.0)
    assert near_ideal > tiny
    assert near_ideal > huge


def test_druggability_score_rewards_hydrophobicity():
    hydrophobic = pk._druggability_score(600.0, hydrophobicity=1.0)
    polar = pk._druggability_score(600.0, hydrophobicity=0.0)
    assert hydrophobic > polar


def test_druggability_score_bounded_0_to_1():
    assert 0.0 <= pk._druggability_score(0.0, 0.0) <= 1.0
    assert 0.0 <= pk._druggability_score(1e6, 1.0) <= 1.0


# --- find_pockets() end-to-end -------------------------------------------

@pytest.fixture()
def cage_structure():
    path = Path(__file__).parent / "data" / "cage.pdb"
    return pio.load_structure(path, structure_id="cage")


def test_find_pockets_detects_central_cavity(cage_structure):
    pockets = pk.find_pockets(
        cage_structure, spacing=1.0, padding=1.0,
        max_ray_length=6.0, ray_step=1.5,
        min_enclosed_axes=5, min_pocket_points=1,
    )
    central = [p for p in pockets if np.linalg.norm(p.centroid) < 0.5]
    assert len(central) == 1
    pocket = central[0]
    assert pocket.n_grid_points >= 1
    assert pocket.volume > 0
    assert len(pocket.residues) == 14  # all 14 cage "wall" residues line it
    assert pocket.hydrophobicity == pytest.approx(1.0)  # all walls are LEU
    assert 0.0 <= pocket.druggability_score <= 1.0


def test_find_pockets_sorted_by_volume_descending(cage_structure):
    pockets = pk.find_pockets(
        cage_structure, spacing=1.0, padding=1.0,
        max_ray_length=6.0, ray_step=1.5,
        min_enclosed_axes=5, min_pocket_points=1,
    )
    volumes = [p.volume for p in pockets]
    assert volumes == sorted(volumes, reverse=True)


def test_find_pockets_no_pockets_when_threshold_too_strict(cage_structure):
    pockets = pk.find_pockets(
        cage_structure, spacing=1.0, padding=1.0,
        max_ray_length=6.0, ray_step=1.5,
        min_enclosed_axes=8,  # impossible: only 7 axes exist
        min_pocket_points=1,
    )
    assert pockets == []


def test_find_pockets_empty_structure_returns_empty(tiny_pdb: Path):
    # tiny.pdb's residues are scattered and far apart, with no enclosed
    # cavity anywhere near the default detection parameters.
    structure = pio.load_structure(tiny_pdb, "t")
    pockets = pk.find_pockets(structure, spacing=1.5, padding=2.0)
    assert pockets == []


def test_find_pockets_grid_size_limit_raises(cage_structure):
    with pytest.raises(ValueError):
        pk.find_pockets(cage_structure, spacing=0.1, padding=5.0, max_grid_points=100)
