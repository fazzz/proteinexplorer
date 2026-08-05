from pathlib import Path

import pytest

from proteinexplorer import io as pio
from proteinexplorer.io import UnknownFormatError


def test_infer_format_pdb():
    assert pio.infer_format("foo.pdb") == "pdb"
    assert pio.infer_format("foo.ent") == "pdb"
    assert pio.infer_format("foo.pdb.gz") == "pdb"


def test_infer_format_mmcif():
    assert pio.infer_format("foo.cif") == "mmcif"
    assert pio.infer_format("foo.mmcif") == "mmcif"


def test_infer_format_unknown_raises():
    with pytest.raises(UnknownFormatError):
        pio.infer_format("foo.txt")


def test_load_structure_and_summarize(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, structure_id="tiny")
    summary = pio.summarize(structure)

    assert summary.n_models == 1
    assert summary.n_chains == 2
    assert summary.n_residues == 8  # 3 (chain A) + 5 (chain B incl. ZN/LIG/HOHx2)
    assert summary.n_atoms == 24
    assert summary.hetero_resnames == ["LIG", "ZN"]

    totals = summary.category_totals()
    assert totals["protein"] == 4
    assert totals["water"] == 2
    assert totals["ion"] == 1
    assert totals["ligand"] == 1


def test_save_and_reload_roundtrip_pdb_to_mmcif(tiny_pdb: Path, tmp_path: Path):
    structure = pio.load_structure(tiny_pdb, structure_id="tiny")
    out_cif = tmp_path / "out.cif"
    pio.save_structure(structure, out_cif, fmt="mmcif")

    reloaded = pio.load_structure(out_cif, structure_id="tiny", fmt="mmcif")
    summary = pio.summarize(reloaded)
    assert summary.n_atoms == 24
    assert summary.category_totals()["water"] == 2
