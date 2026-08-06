from pathlib import Path

import numpy as np
import pytest
from Bio.PDB.Atom import Atom

from proteinexplorer import cluster as clu


def _atom_group(base_coords, jitter=0.0, seed=0):
    rng = np.random.default_rng(seed)
    atoms = []
    for i, coord in enumerate(base_coords):
        c = np.array(coord, dtype=float) + rng.normal(scale=jitter, size=3)
        atoms.append(Atom(f"CA{i}", c, 20.0, 1.0, " ", f"CA{i}", i, element="C"))
    return atoms


BASE = [(0, 0, 0), (3.8, 0, 0), (7.6, 0, 0), (11.4, 0, 0), (15.2, 0, 0)]
# an "L-shaped" chain with the same bond spacing but a very different
# overall fold -- NOT reducible to BASE by rotation+translation, so
# fit-RMSD between the two groups stays large (unlike a pure translation,
# which Kabsch superposition would remove entirely).
BASE_BENT = [(0, 0, 0), (3.8, 0, 0), (7.6, 0, 0), (7.6, 3.8, 0), (7.6, 7.6, 0)]


@pytest.fixture()
def two_group_conformations():
    # 3 near-identical straight conformations ("group A") + 3 near-identical
    # bent conformations ("group B")
    group_a = [_atom_group(BASE, jitter=0.1, seed=i) for i in range(3)]
    group_b = [_atom_group(BASE_BENT, jitter=0.1, seed=100 + i) for i in range(3)]
    labels = ["a1", "a2", "a3", "b1", "b2", "b3"]
    return labels, group_a + group_b


# --- pairwise_rmsd_matrix -----------------------------------------------

def test_pairwise_rmsd_matrix_shape_and_symmetry(two_group_conformations):
    _, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    assert matrix.shape == (6, 6)
    assert np.allclose(matrix, matrix.T)
    assert np.allclose(np.diag(matrix), 0.0)


def test_pairwise_rmsd_matrix_within_group_smaller_than_across(two_group_conformations):
    _, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    within_a = matrix[0, 1]
    across = matrix[0, 3]
    assert within_a < across


def test_pairwise_rmsd_matrix_mismatched_atom_counts_raises():
    a = _atom_group(BASE)
    b = _atom_group(BASE[:3])
    with pytest.raises(clu.ClusterError):
        clu.pairwise_rmsd_matrix([a, b])


def test_pairwise_rmsd_matrix_needs_at_least_two():
    with pytest.raises(clu.ClusterError):
        clu.pairwise_rmsd_matrix([_atom_group(BASE)])


# --- greedy ------------------------------------------------------------

def test_greedy_separates_two_groups(two_group_conformations):
    labels, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    result = clu.greedy(labels, matrix, threshold=2.0)
    assert result.method == "greedy"
    assert len(result.clusters) == 2
    membership = {label: c.id for c in result.clusters for label in c.member_labels}
    assert membership["a1"] == membership["a2"] == membership["a3"]
    assert membership["b1"] == membership["b2"] == membership["b3"]
    assert membership["a1"] != membership["b1"]


def test_greedy_large_threshold_merges_everything(two_group_conformations):
    labels, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    result = clu.greedy(labels, matrix, threshold=1000.0)
    assert len(result.clusters) == 1
    assert set(result.clusters[0].member_labels) == set(labels)


def test_greedy_zero_threshold_splits_everything(two_group_conformations):
    labels, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    result = clu.greedy(labels, matrix, threshold=0.0001)
    assert len(result.clusters) == 6


def test_greedy_representative_is_a_member(two_group_conformations):
    labels, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    result = clu.greedy(labels, matrix, threshold=2.0)
    for c in result.clusters:
        assert c.representative_label in c.member_labels
        assert c.representative_index in c.member_indices


# --- hierarchical --------------------------------------------------------

def test_hierarchical_separates_two_groups_by_n_clusters(two_group_conformations):
    labels, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    result = clu.hierarchical(labels, matrix, n_clusters=2)
    assert result.method == "hierarchical"
    assert len(result.clusters) == 2
    membership = {label: c.id for c in result.clusters for label in c.member_labels}
    assert membership["a1"] == membership["a2"] == membership["a3"]
    assert membership["b1"] == membership["b2"] == membership["b3"]
    assert membership["a1"] != membership["b1"]


def test_hierarchical_separates_two_groups_by_distance_threshold(two_group_conformations):
    labels, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    result = clu.hierarchical(labels, matrix, distance_threshold=1.5)
    assert len(result.clusters) == 2


def test_hierarchical_requires_exactly_one_criterion(two_group_conformations):
    labels, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)
    with pytest.raises(clu.ClusterError):
        clu.hierarchical(labels, matrix, n_clusters=2, distance_threshold=5.0)
    with pytest.raises(clu.ClusterError):
        clu.hierarchical(labels, matrix)


def test_hierarchical_raises_clean_error_without_scipy(two_group_conformations, monkeypatch):
    import builtins
    labels, groups = two_group_conformations
    matrix = clu.pairwise_rmsd_matrix(groups)

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith("scipy"):
            raise ImportError("mocked: scipy not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(clu.ClusterExtraNotAvailableError):
        clu.hierarchical(labels, matrix, n_clusters=2)
