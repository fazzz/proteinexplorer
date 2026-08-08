from pathlib import Path

import pytest

from proteinexplorer import fix as fx
from proteinexplorer import io as pio
from proteinexplorer import model as mdl


@pytest.fixture()
def nonstandard_pdb() -> Path:
    return Path(__file__).parent / "data" / "nonstandard.pdb"


# --- analyze() ---------------------------------------------------------

def test_analyze_finds_missing_atoms(tiny_pdb: Path):
    result = fx.analyze(tiny_pdb)
    assert result.missing_atoms == {"B/VAL1": ["CB", "CG1", "CG2"]}


def test_analyze_finds_missing_terminals(tiny_pdb: Path):
    result = fx.analyze(tiny_pdb)
    assert "OXT" in result.missing_terminals["A/SER3"]
    assert "OXT" in result.missing_terminals["B/VAL1"]


def test_analyze_finds_nonstandard_residue(nonstandard_pdb: Path):
    result = fx.analyze(nonstandard_pdb)
    assert result.nonstandard_residues == [("A/MSE2", "MET")]


def test_analyze_no_missing_residues_without_seqres(tiny_pdb: Path):
    # tiny.pdb has no SEQRES, so PDBFixer's own detection (unlike
    # `prot model gaps`) finds nothing here even though there's no
    # numbering gap to find anyway in this particular fixture.
    result = fx.analyze(tiny_pdb)
    assert result.missing_residues == {}


def test_pdbfixer_gap_detection_blind_spot_vs_model_gaps():
    # The documented differentiator: PDBFixer needs SEQRES to see a
    # numbering-only gap; `prot model gaps` doesn't.
    gapped_path = Path(__file__).parent / "data" / "gapped.pdb"
    analysis = fx.analyze(gapped_path)
    assert analysis.missing_residues == {}

    structure = pio.load_structure(gapped_path, "g")
    gaps = mdl.find_gaps(structure)
    assert len(gaps) == 1
    assert gaps[0].length == 3


# --- fix() ---------------------------------------------------------------

def test_fix_adds_missing_atoms(tiny_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    report = fx.fix(tiny_pdb, out, add_missing_atoms=True)
    assert report.atoms_added == {"B/VAL1": ["CB", "CG1", "CG2"]}
    assert out.exists()

    text = out.read_text()
    assert " CB  VAL B" in text
    assert " CG1 VAL B" in text


def test_fix_skips_missing_atoms_when_disabled(tiny_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    report = fx.fix(tiny_pdb, out, add_missing_atoms=False)
    assert report.atoms_added == {}
    text = out.read_text()
    assert " CG1 VAL B" not in text


def test_fix_does_not_add_whole_residues_by_default(tmp_path: Path):
    gapped_path = Path(__file__).parent / "data" / "gapped.pdb"
    out = tmp_path / "fixed.pdb"
    report = fx.fix(gapped_path, out, add_missing_residues=False)
    assert report.residues_added == {}
    # still just 4 residues (1, 2, 6, 7) -- nothing inserted
    structure = pio.load_structure(out, "g", fmt="pdb")
    assert sum(1 for _ in structure[0].get_residues()) == 4


def test_fix_replaces_nonstandard_residue(nonstandard_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    report = fx.fix(nonstandard_pdb, out, replace_nonstandard=True, add_missing_atoms=False)
    assert report.nonstandard_replaced == [("A/MSE2", "MET")]
    text = out.read_text()
    assert "MSE" not in text
    assert "MET" in text


def test_fix_skips_nonstandard_replacement_when_disabled(nonstandard_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    report = fx.fix(nonstandard_pdb, out, replace_nonstandard=False, add_missing_atoms=False)
    assert report.nonstandard_replaced == []
    assert "MSE" in out.read_text()


def test_fix_remove_heterogens_water_keeps_water(nonstandard_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    report = fx.fix(nonstandard_pdb, out, remove_heterogens="water", add_missing_atoms=False, replace_nonstandard=False)
    assert report.heterogens_removed is True
    text = out.read_text()
    assert "HOH" in text


def test_fix_remove_heterogens_all_removes_water(nonstandard_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    fx.fix(nonstandard_pdb, out, remove_heterogens="all", add_missing_atoms=False, replace_nonstandard=False)
    text = out.read_text()
    assert "HOH" not in text


def test_fix_remove_heterogens_none_by_default(nonstandard_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    report = fx.fix(nonstandard_pdb, out, add_missing_atoms=False, replace_nonstandard=False)
    assert report.heterogens_removed is False
    assert "HOH" in out.read_text()


def test_fix_add_hydrogens(tiny_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    report = fx.fix(tiny_pdb, out, add_missing_atoms=True, add_hydrogens_ph=7.0)
    assert report.hydrogens_added_at_ph == pytest.approx(7.0)
    structure = pio.load_structure(out, "g", fmt="pdb")
    elements = {a.element for a in structure[0].get_atoms()}
    assert "H" in elements


def test_fix_no_hydrogens_by_default(tiny_pdb: Path, tmp_path: Path):
    out = tmp_path / "fixed.pdb"
    fx.fix(tiny_pdb, out, add_missing_atoms=True)
    structure = pio.load_structure(out, "g", fmt="pdb")
    elements = {a.element for a in structure[0].get_atoms()}
    assert "H" not in elements
