from pathlib import Path

import numpy as np
import pytest

from proteinexplorer import io as pio
from proteinexplorer import model as mdl
from proteinexplorer.mutate import MutationError


@pytest.fixture()
def gapped_structure():
    path = Path(__file__).parent / "data" / "gapped.pdb"
    return pio.load_structure(path, structure_id="gapped")


# --- find_gaps -------------------------------------------------------

def test_find_gaps_detects_missing_range(gapped_structure):
    gaps = mdl.find_gaps(gapped_structure)
    assert len(gaps) == 1
    gap = gaps[0]
    assert gap.chain_id == "A"
    assert gap.prev_resseq == 2
    assert gap.next_resseq == 6
    assert gap.length == 3


def test_find_gaps_none_for_contiguous_chain(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    assert mdl.find_gaps(structure) == []


# --- fill_loop_linear --------------------------------------------------

def test_fill_loop_linear_adds_placeholder_residues(gapped_structure):
    result = mdl.fill_loop_linear(gapped_structure, "A", 3, 5)
    assert result.chain_id == "A"
    assert result.residues_added == ["ALA3", "ALA4", "ALA5"]
    assert "Crude placeholder" in result.note


def test_fill_loop_linear_closes_the_gap(gapped_structure):
    mdl.fill_loop_linear(gapped_structure, "A", 3, 5)
    assert mdl.find_gaps(gapped_structure) == []


def test_fill_loop_linear_ca_trace_is_between_anchors(gapped_structure):
    prev_ca = gapped_structure[0]["A"][2]["CA"].coord.copy()
    next_ca = gapped_structure[0]["A"][6]["CA"].coord.copy()
    mdl.fill_loop_linear(gapped_structure, "A", 3, 5)

    chain = gapped_structure[0]["A"]
    for resseq in (3, 4, 5):
        ca = chain[resseq]["CA"].coord
        # every interpolated CA must lie on the straight segment between anchors
        assert prev_ca[1] == pytest.approx(next_ca[1])  # sanity: same y/z in fixture
        t = (ca[0] - prev_ca[0]) / (next_ca[0] - prev_ca[0])
        assert 0.0 < t < 1.0
        expected = prev_ca + (next_ca - prev_ca) * t
        assert np.allclose(ca, expected, atol=1e-6)


def test_fill_loop_linear_residues_are_sequence_ordered(gapped_structure):
    mdl.fill_loop_linear(gapped_structure, "A", 3, 5)
    resseqs = [r.id[1] for r in gapped_structure[0]["A"]]
    assert resseqs == sorted(resseqs)


def test_fill_loop_linear_with_explicit_sequence(gapped_structure):
    result = mdl.fill_loop_linear(gapped_structure, "A", 3, 5, sequence="GVL")
    assert result.residues_added == ["GLY3", "VAL4", "LEU5"]


def test_fill_loop_linear_wrong_sequence_length_raises(gapped_structure):
    with pytest.raises(MutationError):
        mdl.fill_loop_linear(gapped_structure, "A", 3, 5, sequence="GV")


def test_fill_loop_linear_missing_anchor_raises(gapped_structure):
    with pytest.raises(MutationError):
        mdl.fill_loop_linear(gapped_structure, "A", 10, 12)


def test_fill_loop_linear_unknown_chain_raises(gapped_structure):
    with pytest.raises(MutationError):
        mdl.fill_loop_linear(gapped_structure, "Z", 3, 5)


def test_fill_loop_linear_end_before_start_raises(gapped_structure):
    with pytest.raises(MutationError):
        mdl.fill_loop_linear(gapped_structure, "A", 5, 3)


# --- homology_model (external tool only, not installed here) --------------

def test_homology_model_raises_without_modeller_package(tmp_path):
    with pytest.raises(mdl.ModellerNotAvailableError):
        mdl.homology_model(
            alignment_pir_path=tmp_path / "align.pir",
            template_codes=["template"],
            target_code="target",
            template_search_dir=tmp_path,
            output_dir=tmp_path / "out",
        )
