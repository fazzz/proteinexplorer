from pathlib import Path

import pytest

from proteinexplorer import assembly as asm
from proteinexplorer import io as pio


@pytest.fixture()
def a1a8o_path() -> Path:
    path = Path(__file__).parent.parent / "examples" / "1a8o" / "1A8O.pdb"
    if not path.exists():
        pytest.skip("examples/1a8o/1A8O.pdb not present")
    return path


# --- list_assemblies -----------------------------------------------------

def test_list_assemblies_finds_documented_dimer(a1a8o_path: Path):
    infos = asm.list_assemblies(a1a8o_path)
    assert len(infos) == 1
    info = infos[0]
    assert info.name == "1"
    assert info.oligomeric_details == "DIMERIC"
    assert info.chains_involved == ["A"]
    assert info.n_operators == 2


def test_list_assemblies_no_assembly_documented(tiny_pdb: Path):
    infos = asm.list_assemblies(tiny_pdb)
    assert infos == []


# --- generate_assembly ---------------------------------------------------

def test_generate_assembly_doubles_chains_and_atoms(a1a8o_path: Path, tmp_path: Path):
    out = tmp_path / "dimer.pdb"
    result = asm.generate_assembly(a1a8o_path, out)
    assert result.assembly_name == "1"
    assert result.chains_before == ["A"]
    assert result.chains_after == ["A1", "A2"]
    assert result.n_atoms_after == 1288
    assert out.exists()


def test_generate_assembly_output_loads_via_biopython_pipeline(a1a8o_path: Path, tmp_path: Path):
    # The whole point: gemmi is used only for this one gap, and the
    # result must be usable by every other (Bio.PDB-based) module
    # without any special-casing.
    out = tmp_path / "dimer.pdb"
    asm.generate_assembly(a1a8o_path, out)

    structure = pio.load_structure(out, "dimer")
    summary = pio.summarize(structure)
    assert summary.n_chains == 2
    assert summary.n_atoms == 1288
    assert summary.category_totals()["protein"] == 140


def test_generate_assembly_mmcif_output(a1a8o_path: Path, tmp_path: Path):
    out = tmp_path / "dimer.cif"
    result = asm.generate_assembly(a1a8o_path, out)
    assert out.exists()
    assert result.n_atoms_after == 1288
    structure = pio.load_structure(out, "dimer", fmt="mmcif")
    assert pio.summarize(structure).n_chains == 2


def test_generate_assembly_no_documented_assembly_raises(tiny_pdb: Path, tmp_path: Path):
    with pytest.raises(ValueError, match="no documented"):
        asm.generate_assembly(tiny_pdb, tmp_path / "out.pdb")


def test_generate_assembly_unknown_how_raises(a1a8o_path: Path, tmp_path: Path):
    with pytest.raises(ValueError, match="how"):
        asm.generate_assembly(a1a8o_path, tmp_path / "out.pdb", how="bogus")


def test_generate_assembly_duplicate_naming_keeps_original_chain_id(a1a8o_path: Path, tmp_path: Path):
    out = tmp_path / "dimer_dup.pdb"
    result = asm.generate_assembly(a1a8o_path, out, how="duplicate")
    assert result.chains_after == ["A", "A"]
