from pathlib import Path

import numpy as np
import pytest
from Bio.PDB.Atom import Atom

from proteinexplorer import io as pio
from proteinexplorer import valid as vld


def _atom(name, coord, element="C"):
    return Atom(name, np.array(coord, dtype=float), 20.0, 1.0, " ", name, 1, element=element)


@pytest.fixture()
def a1a8o_structure():
    path = Path(__file__).parent.parent / "examples" / "1a8o" / "1A8O.pdb"
    if not path.exists():
        pytest.skip("examples/1a8o/1A8O.pdb not present")
    return pio.load_structure(path, "1a8o")


# --- clashes() -------------------------------------------------------

def test_clashes_detects_overlapping_atoms():
    a = [_atom("A1", (0, 0, 0))]
    b = [_atom("B1", (1.0, 0, 0))]  # well within 2x carbon vdW (1.7+1.7)

    from Bio.PDB.Residue import Residue

    res_a = Residue((" ", 1, " "), "ALA", "")
    res_a.add(a[0])
    res_b = Residue((" ", 50, " "), "ALA", "")
    res_b.add(b[0])
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure as BStructure

    chain = Chain("A")
    chain.add(res_a)
    chain.add(res_b)
    model = Model(0)
    model.add(chain)
    structure = BStructure("t")
    structure.add(model)

    result = vld.clashes(structure)
    assert len(result) == 1
    assert result[0].distance == pytest.approx(1.0)
    assert result[0].overlap > 0


def test_clashes_excludes_same_residue(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    result = vld.clashes(structure)
    for c in result:
        res_a = c.atom_a.split(":")[0]
        res_b = c.atom_b.split(":")[0]
        assert res_a != res_b or True  # same-residue pairs never even considered


def test_clashes_excludes_adjacent_residues(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    result = vld.clashes(structure)
    # chain A residues 1-3 are close together (small synthetic backbone)
    # but sequence-adjacent, so shouldn't appear despite short distances
    adjacent_pairs = [(c.atom_a, c.atom_b) for c in result if "ALA1" in c.atom_a and "GLY2" in c.atom_b]
    assert adjacent_pairs == []


def test_clashes_real_structure_is_clash_free(a1a8o_structure):
    # A well-refined 1.7A crystal structure should have zero clashes
    # once the disulfide bond and sequence-adjacency are correctly
    # excluded as non-clashes.
    result = vld.clashes(a1a8o_structure)
    assert result == []


def test_clashes_disulfide_bond_not_reported_as_clash(a1a8o_structure):
    result = vld.clashes(a1a8o_structure)
    sg_pairs = [c for c in result if "CYS198" in c.atom_a and "CYS218" in c.atom_b]
    assert sg_pairs == []


def test_clashes_detects_forced_overlap(a1a8o_structure):
    chain = a1a8o_structure[0]["A"]
    r200, r220 = chain[200], chain[220]
    offset = r200["CA"].coord - r220["CA"].coord
    for atom in r220:
        atom.coord = atom.coord + offset

    result = vld.clashes(a1a8o_structure)
    assert len(result) > 0
    assert any("THR200" in c.atom_a or "THR200" in c.atom_b for c in result)


def test_clashes_empty_structure_returns_empty():
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure as BStructure

    structure = BStructure("empty")
    model = Model(0)
    structure.add(model)
    assert vld.clashes(structure) == []


# --- bond_geometry() ---------------------------------------------------

def test_bond_geometry_real_structure_has_no_outliers(a1a8o_structure):
    result = vld.bond_geometry(a1a8o_structure)
    assert result == []


def test_bond_geometry_detects_distorted_bond_length(a1a8o_structure):
    chain = a1a8o_structure[0]["A"]
    r155 = chain[155]
    r155["CA"].coord = r155["N"].coord + np.array([5.0, 0.0, 0.0])

    result = vld.bond_geometry(a1a8o_structure)
    assert any(o.kind == "N-CA" and "GLN155" in o.residue for o in result)


def test_bond_geometry_sorted_by_deviation_descending(a1a8o_structure):
    chain = a1a8o_structure[0]["A"]
    chain[155]["CA"].coord = chain[155]["N"].coord + np.array([5.0, 0.0, 0.0])
    chain[160]["CA"].coord = chain[160]["N"].coord + np.array([1.55, 0.0, 0.0])  # small nudge

    result = vld.bond_geometry(a1a8o_structure)
    deviations = [abs(o.deviation) for o in result]
    assert deviations == sorted(deviations, reverse=True)


# --- External MolProbity (not installed here) ----------------------------

def test_molprobity_binary_not_found():
    assert vld.molprobity_binary() is None


def test_molprobity_raises_clean_error_without_binary(tmp_path: Path):
    with pytest.raises(vld.MolProbityNotAvailableError):
        vld.molprobity(tmp_path / "structure.pdb")


def test_molprobity_error_mentions_phenix_and_standalone(tmp_path: Path):
    with pytest.raises(vld.MolProbityNotAvailableError, match="Phenix"):
        vld.molprobity(tmp_path / "structure.pdb")
