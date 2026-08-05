from pathlib import Path

import numpy as np
import pytest

from proteinexplorer import io as pio
from proteinexplorer import mutate as mut


# --- normalize_resname ---------------------------------------------------

def test_normalize_resname_one_letter():
    assert mut.normalize_resname("v") == "VAL"
    assert mut.normalize_resname("A") == "ALA"


def test_normalize_resname_three_letter():
    assert mut.normalize_resname("val") == "VAL"
    assert mut.normalize_resname("GLY") == "GLY"


def test_normalize_resname_unknown_one_letter_raises():
    with pytest.raises(mut.MutationError):
        mut.normalize_resname("X")


def test_normalize_resname_unknown_three_letter_raises():
    with pytest.raises(mut.MutationError):
        mut.normalize_resname("XYZ")


# --- cb_only ---------------------------------------------------------

def _bond_angle(a, vertex, b):
    v1 = a - vertex
    v2 = b - vertex
    cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cos_t, -1.0, 1.0)))


def test_cb_only_renames_residue_and_keeps_backbone(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    orig_n = structure[0]["A"][1]["N"].coord.copy()
    orig_ca = structure[0]["A"][1]["CA"].coord.copy()

    result = mut.cb_only(structure, "A", 1, "VAL")

    assert result.original_resname == "ALA"
    assert result.new_resname == "VAL"
    assert result.method == "cb_only"
    assert set(result.atoms_placed) == {"N", "CA", "C", "O", "CB"}

    new_residue = structure[0]["A"][1]
    assert new_residue.resname == "VAL"
    assert np.allclose(new_residue["N"].coord, orig_n)
    assert np.allclose(new_residue["CA"].coord, orig_ca)


def test_cb_only_glycine_has_no_cb(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    result = mut.cb_only(structure, "A", 1, "GLY")
    assert "CB" not in result.atoms_placed
    assert "CB" not in structure[0]["A"][1]
    assert "no C-beta" in result.note.lower() or "no c-beta" in result.note.lower()


def test_cb_only_note_mentions_fallback_limitation(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    result = mut.cb_only(structure, "A", 1, "TRP")
    assert "Scwrl4" in result.note


def test_cb_only_virtual_cb_bond_length_and_angles_are_plausible(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    n = structure[0]["A"][1]["N"].coord.copy()
    ca = structure[0]["A"][1]["CA"].coord.copy()
    c = structure[0]["A"][1]["C"].coord.copy()

    mut.cb_only(structure, "A", 1, "LEU")
    cb = structure[0]["A"][1]["CB"].coord

    bond_length = np.linalg.norm(cb - ca)
    assert 1.4 <= bond_length <= 1.7  # ~1.53 A for a Csp3-Csp3 bond

    angle_n = _bond_angle(n, ca, cb)
    angle_c = _bond_angle(c, ca, cb)
    assert 95 <= angle_n <= 125
    assert 95 <= angle_c <= 125


def test_cb_only_missing_backbone_atom_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    structure[0]["A"][1].detach_child("C")
    with pytest.raises(mut.MutationError):
        mut.cb_only(structure, "A", 1, "VAL")


def test_cb_only_unknown_chain_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(mut.MutationError):
        mut.cb_only(structure, "Z", 1, "VAL")


def test_cb_only_unknown_resid_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(mut.MutationError):
        mut.cb_only(structure, "A", 999, "VAL")


# --- Scwrl4 wrapper (binary not installed in this environment) -----------

def test_scwrl4_binary_not_found():
    assert mut.scwrl4_binary() is None


def test_scwrl4_raises_clean_error_without_binary(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(mut.Scwrl4NotAvailableError):
        mut.scwrl4(structure, tiny_pdb, "A", 1, "VAL")


# --- mutate_residue (unified entry point) ---------------------------------

def test_mutate_residue_auto_falls_back_to_cb_only(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    result = mut.mutate_residue(structure, tiny_pdb, "A", 1, "VAL", method="auto")
    assert result.method == "cb_only"
    assert structure[0]["A"][1].resname == "VAL"


def test_mutate_residue_explicit_scwrl4_raises_without_binary(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(mut.Scwrl4NotAvailableError):
        mut.mutate_residue(structure, tiny_pdb, "A", 1, "VAL", method="scwrl4")


def test_mutate_residue_accepts_one_letter_target(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    result = mut.mutate_residue(structure, tiny_pdb, "A", 1, "v", method="cb_only")
    assert result.new_resname == "VAL"


def test_mutate_residue_rejects_non_standard_target(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(mut.MutationError):
        mut.mutate_residue(structure, tiny_pdb, "A", 1, "MSE", method="cb_only")


def test_mutate_residue_unknown_method_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(ValueError):
        mut.mutate_residue(structure, tiny_pdb, "A", 1, "VAL", method="bogus")
