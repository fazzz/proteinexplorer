from pathlib import Path

import pytest

from proteinexplorer import geometry as geom
from proteinexplorer import io as pio
from proteinexplorer import secondary as sec


@pytest.fixture()
def eight_residue_structure():
    path = Path(__file__).parent / "data" / "eight_residues.pdb"
    return pio.load_structure(path, structure_id="eight")


# --- Ramachandran region classification -----------------------------------

def test_alpha_region_canonical_point():
    assert sec._in_alpha_region(-57.0, -47.0)


def test_alpha_region_excludes_beta_point():
    assert not sec._in_alpha_region(-120.0, 130.0)


def test_beta_region_canonical_point():
    assert sec._in_beta_region(-120.0, 130.0)


def test_beta_region_excludes_alpha_point():
    assert not sec._in_beta_region(-57.0, -47.0)


# --- Run smoothing -----------------------------------------------------

def _ss(chain_id, resseq, code, resname="GLY"):
    return sec.ResidueSS(chain_id=chain_id, resseq=resseq, resname=resname, code=code, phi=None, psi=None, method="geometric")


def test_smooth_runs_keeps_long_helix():
    residues = [_ss("A", i, "H") for i in range(1, 6)]  # 5-residue helix run
    smoothed = sec._smooth_runs(residues, min_helix_run=4, min_strand_run=2)
    assert all(r.code == "H" for r in smoothed)


def test_smooth_runs_demotes_short_helix():
    residues = [_ss("A", i, "H") for i in range(1, 3)]  # 2-residue "helix"
    smoothed = sec._smooth_runs(residues, min_helix_run=4, min_strand_run=2)
    assert all(r.code == "C" for r in smoothed)


def test_smooth_runs_keeps_long_strand():
    residues = [_ss("A", i, "E") for i in range(1, 4)]
    smoothed = sec._smooth_runs(residues, min_helix_run=4, min_strand_run=2)
    assert all(r.code == "E" for r in smoothed)


def test_smooth_runs_demotes_single_strand():
    residues = [_ss("A", 1, "E")]
    smoothed = sec._smooth_runs(residues, min_helix_run=4, min_strand_run=2)
    assert smoothed[0].code == "C"


def test_smooth_runs_mixed_sequence():
    residues = (
        [_ss("A", i, "H") for i in range(1, 6)]  # kept (5 >= 4)
        + [_ss("A", i, "C") for i in range(6, 8)]
        + [_ss("A", 8, "E")]  # isolated single strand -> demoted (1 < 2)
        + [_ss("A", 9, "C")]
        + [_ss("A", i, "E") for i in range(10, 13)]  # kept (3 >= 2)
    )
    smoothed = sec._smooth_runs(residues, min_helix_run=4, min_strand_run=2)
    codes = [r.code for r in smoothed]
    assert codes == ["H", "H", "H", "H", "H", "C", "C", "C", "C", "E", "E", "E"]


def test_smooth_runs_does_not_bridge_chain_boundary():
    residues = [_ss("A", 1, "H"), _ss("A", 2, "H"), _ss("B", 1, "H"), _ss("B", 2, "H")]
    smoothed = sec._smooth_runs(residues, min_helix_run=4, min_strand_run=2)
    # each chain only contributes a 2-residue run -> both demoted, not merged into one run of 4
    assert all(r.code == "C" for r in smoothed)


def test_smooth_runs_empty_list():
    assert sec._smooth_runs([], min_helix_run=4, min_strand_run=2) == []


# --- geometric() end-to-end with monkeypatched torsions ---------------------

def test_geometric_classifies_and_smooths_helix(eight_residue_structure, monkeypatch):
    # Residues 2-6 (5 residues) sit in the alpha region -> should survive
    # smoothing as a helix; residues 1, 7, 8 lack a full phi/psi pair or
    # sit outside any region -> coil.
    phi_psi = {
        1: (None, -47.0),
        2: (-57.0, -47.0),
        3: (-57.0, -47.0),
        4: (-57.0, -47.0),
        5: (-57.0, -47.0),
        6: (-57.0, -47.0),
        7: (-57.0, None),
        8: (None, None),
    }

    def fake_backbone_torsions(chain, resseq):
        phi, psi = phi_psi[resseq]
        return geom.ResidueTorsions(chain_id=chain.id, resseq=resseq, resname="GLY", phi=phi, psi=psi, omega=None, chi=[])

    monkeypatch.setattr(sec.geom, "backbone_torsions", fake_backbone_torsions)

    result = sec.geometric(eight_residue_structure)
    codes = {r.resseq: r.code for r in result}
    assert codes[1] == "C"
    for i in range(2, 7):
        assert codes[i] == "H"
    assert codes[7] == "C"
    assert codes[8] == "C"


def test_geometric_on_tiny_pdb_is_all_coil(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    result = sec.geometric(structure)
    assert all(r.code == "C" for r in result)


# --- DSSP wrapper ------------------------------------------------------

def test_dssp_raises_clean_error_without_binary(tiny_pdb: Path, monkeypatch):
    structure = pio.load_structure(tiny_pdb, "t")
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(sec.DSSPNotAvailableError):
        sec.dssp(structure, tiny_pdb)


# --- Unified entry point -------------------------------------------------

def test_secondary_structure_geometric_method_explicit(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    residues, method = sec.secondary_structure(structure, method="geometric")
    assert method == "geometric"
    assert len(residues) == 4


def test_secondary_structure_dssp_method_requires_pdb_path(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(ValueError):
        sec.secondary_structure(structure, pdb_path=None, method="dssp")


def test_secondary_structure_auto_falls_back_to_geometric_without_dssp(tiny_pdb: Path, monkeypatch):
    structure = pio.load_structure(tiny_pdb, "t")
    monkeypatch.setattr("shutil.which", lambda name: None)
    residues, method = sec.secondary_structure(structure, pdb_path=tiny_pdb, method="auto")
    assert method == "geometric"


def test_secondary_structure_unknown_method_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(ValueError):
        sec.secondary_structure(structure, method="bogus")


# --- composition ---------------------------------------------------------

def test_composition_fractions_sum_to_one():
    residues = [_ss("A", 1, "H"), _ss("A", 2, "H"), _ss("A", 3, "C"), _ss("A", 4, "E")]
    comp = sec.composition(residues)
    assert comp["H"] == pytest.approx(0.5)
    assert comp["C"] == pytest.approx(0.25)
    assert comp["E"] == pytest.approx(0.25)
    assert sum(comp.values()) == pytest.approx(1.0)


def test_composition_empty_list():
    assert sec.composition([]) == {}
