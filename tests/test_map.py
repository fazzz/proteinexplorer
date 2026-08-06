from pathlib import Path

import pytest

from proteinexplorer import io as pio
from proteinexplorer import map as map_mod


@pytest.fixture()
def tiny_structure(tiny_pdb: Path):
    return pio.load_structure(tiny_pdb, "t")


# --- _parse_label ---------------------------------------------------

def test_parse_label_chain_resnum():
    assert map_mod._parse_label("A/41") == ("A", 41)


def test_parse_label_chain_resname_resnum():
    assert map_mod._parse_label("A/GLY41") == ("A", 41)


def test_parse_label_no_digits_raises():
    with pytest.raises(map_mod.MapError):
        map_mod._parse_label("A/GLY")


# --- generate_group_script -------------------------------------------

def test_generate_group_script_pymol_contains_selection_and_color():
    groups = [map_mod.ResidueGroup(label="site1", residues=["A/10", "A/11"], color="red")]
    script = map_mod.generate_group_script(groups, tool="pymol")
    assert "color red" in script
    assert "chain A and resi 10" in script
    assert "chain A and resi 11" in script


def test_generate_group_script_chimerax():
    groups = [map_mod.ResidueGroup(label="site1", residues=["A/10"], color="blue")]
    script = map_mod.generate_group_script(groups, tool="chimerax")
    assert "/A:10" in script
    assert "blue" in script


def test_generate_group_script_vmd():
    groups = [map_mod.ResidueGroup(label="site1", residues=["A/10", "A/11"], color="green")]
    script = map_mod.generate_group_script(groups, tool="vmd")
    assert "atomselect" in script
    assert "resid 10 11" in script


def test_generate_group_script_unknown_tool_raises():
    groups = [map_mod.ResidueGroup(label="s", residues=["A/1"], color="red")]
    with pytest.raises(map_mod.MapError):
        map_mod.generate_group_script(groups, tool="bogus")


def test_generate_group_script_empty_raises():
    with pytest.raises(map_mod.MapError):
        map_mod.generate_group_script([], tool="pymol")


def test_assign_colors_cycles_palette():
    groups = map_mod.assign_colors({"a": ["A/1"], "b": ["A/2"]}, palette=["red", "blue"])
    assert groups[0].color == "red"
    assert groups[1].color == "blue"


# --- pocket / mutation / domain convenience wrappers ----------------

def test_pocket_map_script_one_color_per_pocket():
    class FakePocket:
        def __init__(self, id, residues):
            self.id = id
            self.residues = residues

    pockets = [FakePocket(1, ["A/1", "A/2"]), FakePocket(2, ["A/5"])]
    script = map_mod.pocket_map_script(pockets, tool="pymol")
    assert "pocket_1" in script
    assert "pocket_2" in script


def test_pocket_map_script_skips_empty_pockets():
    class FakePocket:
        def __init__(self, id, residues):
            self.id = id
            self.residues = residues

    pockets = [FakePocket(1, [])]
    with pytest.raises(map_mod.MapError):
        map_mod.pocket_map_script(pockets, tool="pymol")


def test_mutation_map_script_single_group():
    script = map_mod.mutation_map_script(["A/50", "B/12"], tool="pymol", color="magenta")
    assert "color magenta" in script
    assert "chain A and resi 50" in script
    assert "chain B and resi 12" in script


def test_domain_map_script_expands_range():
    domains = [map_mod.DomainRange(label="PF00062", chain_id="A", start=10, end=12)]
    script = map_mod.domain_map_script(domains, tool="pymol")
    assert "resi 10" in script
    assert "resi 11" in script
    assert "resi 12" in script


# --- write_bfactors / spectrum_script ---------------------------------

def test_write_bfactors_sets_values(tiny_structure):
    map_mod.write_bfactors(tiny_structure, {"A/1": 0.9, "A/2": 0.1})
    chain = tiny_structure[0]["A"]
    assert chain[1]["CA"].bfactor == pytest.approx(0.9)
    assert chain[2]["CA"].bfactor == pytest.approx(0.1)


def test_write_bfactors_uses_default_for_unlisted_residues(tiny_structure):
    map_mod.write_bfactors(tiny_structure, {"A/1": 0.9}, default=-1.0)
    chain = tiny_structure[0]["A"]
    assert chain[3]["CA"].bfactor == pytest.approx(-1.0)


def test_spectrum_script_pymol():
    assert "spectrum b" in map_mod.spectrum_script("pymol")


def test_spectrum_script_chimerax():
    assert "color byattribute" in map_mod.spectrum_script("chimerax")


def test_spectrum_script_unknown_tool_raises():
    with pytest.raises(map_mod.MapError):
        map_mod.spectrum_script("bogus")
