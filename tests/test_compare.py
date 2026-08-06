from pathlib import Path

import numpy as np
import pytest

from proteinexplorer import compare as cmp
from proteinexplorer import io as pio


@pytest.fixture()
def tiny_pair():
    a = pio.load_structure(Path(__file__).parent / "data" / "tiny.pdb", "a")
    b = pio.load_structure(Path(__file__).parent / "data" / "tiny.pdb", "b")
    return a, b


@pytest.fixture()
def cage_pair():
    a = pio.load_structure(Path(__file__).parent / "data" / "cage.pdb", "a")
    b = pio.load_structure(Path(__file__).parent / "data" / "cage.pdb", "b")
    return a, b


# --- RMSD ------------------------------------------------------------------

def test_rmsd_identical_structures_is_zero(tiny_pair):
    a, b = tiny_pair
    value, n = cmp.rmsd(a, b)
    assert value == pytest.approx(0.0, abs=1e-6)
    assert n == 4


def test_rmsd_perturbed_structure_is_nonzero(tiny_pair):
    a, b = tiny_pair
    b[0]["A"][2]["CA"].coord = b[0]["A"][2]["CA"].coord + np.array([5.0, 0.0, 0.0])
    value, n = cmp.rmsd(a, b)
    assert value > 0.5
    assert n == 4


def test_rmsd_no_fit_differs_from_fit(tiny_pair):
    a, b = tiny_pair
    b[0]["A"][2]["CA"].coord = b[0]["A"][2]["CA"].coord + np.array([5.0, 0.0, 0.0])
    fit_value, _ = cmp.rmsd(a, b, fit=True)
    nofit_value, _ = cmp.rmsd(a, b, fit=False)
    assert fit_value != pytest.approx(nofit_value)


def test_rmsd_no_common_residues_raises():
    a = pio.load_structure(Path(__file__).parent / "data" / "tiny.pdb", "a")
    b = pio.load_structure(Path(__file__).parent / "data" / "tiny.pdb", "b")
    for chain in b[0]:
        chain.id = f"Z{chain.id}"  # guarantee disjoint (chain_id, resseq) labels
    with pytest.raises(cmp.CompareError):
        cmp.rmsd(a, b)


# --- TM-score fallback ---------------------------------------------------

def test_tm_score_fallback_identical_is_one(tiny_pair):
    a, b = tiny_pair
    result = cmp.tm_score_fallback(a, b)
    assert result.score == pytest.approx(1.0, abs=1e-6)
    assert result.method == "fallback"
    assert result.n_residues == 4


def test_tm_score_fallback_perturbed_is_lower(tiny_pair):
    a, b = tiny_pair
    b[0]["A"][2]["CA"].coord = b[0]["A"][2]["CA"].coord + np.array([20.0, 0.0, 0.0])
    result = cmp.tm_score_fallback(a, b)
    assert 0.0 <= result.score < 0.9


def test_tm_score_bounded_0_to_1(tiny_pair):
    a, b = tiny_pair
    b[0]["A"][2]["CA"].coord = b[0]["A"][2]["CA"].coord + np.array([1000.0, 0.0, 0.0])
    result = cmp.tm_score_fallback(a, b)
    assert 0.0 <= result.score <= 1.0


def test_tm_score_auto_falls_back_without_tmalign(tiny_pair, monkeypatch):
    a, b = tiny_pair
    monkeypatch.setattr("shutil.which", lambda name: None)
    result = cmp.tm_score(a, b, structure_a_path="a.pdb", structure_b_path="b.pdb", method="auto")
    assert result.method == "fallback"


def test_tm_score_explicit_tmalign_requires_paths(tiny_pair):
    a, b = tiny_pair
    with pytest.raises(ValueError):
        cmp.tm_score(a, b, method="tmalign")


def test_tm_score_unknown_method_raises(tiny_pair):
    a, b = tiny_pair
    with pytest.raises(ValueError):
        cmp.tm_score(a, b, method="bogus")


# --- Secondary structure similarity -----------------------------------

def test_secondary_structure_similarity_identical_is_one(tiny_pair):
    a, b = tiny_pair
    score, n = cmp.secondary_structure_similarity(a, b, method="geometric")
    assert score == pytest.approx(1.0)
    assert n == 4


# --- Contact similarity --------------------------------------------------

def test_contact_similarity_identical_is_one(tiny_pair):
    a, b = tiny_pair
    jaccard, shared, union = cmp.contact_similarity(a, b)
    assert jaccard == pytest.approx(1.0)
    assert shared == union


def test_contact_similarity_differs_when_perturbed(tiny_pair):
    a, b = tiny_pair
    # move a residue far away so its contacts disappear in b
    b[0]["A"][3]["CA"].coord = b[0]["A"][3]["CA"].coord + np.array([100.0, 0.0, 0.0])
    jaccard, shared, union = cmp.contact_similarity(a, b)
    assert jaccard < 1.0


# --- Pocket overlap ------------------------------------------------------

def test_pocket_overlap_identical_structures(cage_pair):
    a, b = cage_pair
    jaccard, res_a, res_b = cmp.pocket_overlap(
        a, b, spacing=1.0, padding=1.0, min_pocket_points=1
    )
    assert jaccard == pytest.approx(1.0)
    assert set(res_a) == set(res_b)
    assert len(res_a) == 14


def test_pocket_overlap_missing_pocket_raises(cage_pair, tiny_pdb: Path):
    a, _ = cage_pair
    tiny = pio.load_structure(tiny_pdb, "tiny")
    with pytest.raises(cmp.CompareError):
        cmp.pocket_overlap(a, tiny, spacing=1.0, padding=1.0, min_pocket_points=1)


# --- Ligand comparison -----------------------------------------------------

def test_ligand_comparison_identical_structures(tiny_pair):
    a, b = tiny_pair
    result = cmp.ligand_comparison(a, b)
    assert result.common_resnames == ["LIG"]
    assert result.only_in_a == []
    assert result.only_in_b == []
    assert result.rmsd_by_resname["LIG"] == pytest.approx(0.0, abs=1e-6)


def test_ligand_comparison_no_common_ligand(tiny_pdb: Path, cage_pair):
    tiny = pio.load_structure(tiny_pdb, "tiny")
    cage, _ = cage_pair
    result = cmp.ligand_comparison(tiny, cage)
    assert result.common_resnames == []
    assert "LIG" in result.only_in_a


def test_ligand_comparison_perturbed_ligand_rmsd_nonzero(tiny_pair):
    a, b = tiny_pair
    # move the LIG ligand in b
    lig_b = next(r for r in b[0]["B"] if r.resname == "LIG")
    for atom in lig_b:
        atom.coord = atom.coord + np.array([3.0, 0.0, 0.0])
    result = cmp.ligand_comparison(a, b, fit=False)
    assert result.rmsd_by_resname["LIG"] == pytest.approx(3.0, abs=1e-6)
