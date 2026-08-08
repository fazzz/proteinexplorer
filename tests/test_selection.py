from pathlib import Path

import pytest

from proteinexplorer import io as pio
from proteinexplorer import selection as sel
from proteinexplorer.selection import SelectionSyntaxError


def _atom_names(atoms):
    return sorted(f"{a.get_parent().get_parent().id}/{a.get_parent().resname}{a.get_parent().id[1]}/{a.get_name()}" for a in atoms)


def test_select_all(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    atoms = sel.select(structure, "all")
    assert len(atoms) == 24


def test_select_chain(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    atoms = sel.select(structure, "chain A")
    assert len(atoms) == 15  # ALA(5)+GLY(4)+SER(6) in chain A


def test_select_protein_category(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    atoms = sel.select(structure, "protein")
    residues = {a.get_parent().resname for a in atoms}
    assert residues == {"ALA", "GLY", "SER", "VAL"}


def test_select_water(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    atoms = sel.select(structure, "water")
    assert len(atoms) == 2


def test_select_ligand_and_ion(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    assert len(sel.select(structure, "ligand")) == 2  # LIG has 2 atoms (C1, O1)
    assert len(sel.select(structure, "ion")) == 1  # ZN


def test_select_resid_range(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    atoms = sel.select(structure, "chain A and resid 1:2")
    resnames = {a.get_parent().resname for a in atoms}
    assert resnames == {"ALA", "GLY"}


def test_select_resname(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    atoms = sel.select(structure, "resname ALA")
    assert all(a.get_parent().resname == "ALA" for a in atoms)
    assert len(atoms) == 5


def test_select_atom_name(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    atoms = sel.select(structure, "atom CA")
    assert all(a.get_name() == "CA" for a in atoms)
    assert len(atoms) == 4  # 4 protein residues, each with a CA


def test_select_backbone_and_sidechain_are_complementary_for_protein(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    protein = set(sel.select(structure, "protein"))
    backbone = set(sel.select(structure, "protein and backbone"))
    sidechain = set(sel.select(structure, "protein and sidechain"))
    assert backbone | sidechain == protein
    assert backbone & sidechain == set()


def test_select_not(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    all_atoms = set(sel.select(structure, "all"))
    water = set(sel.select(structure, "water"))
    not_water = set(sel.select(structure, "not water"))
    assert not_water == all_atoms - water


def test_select_and_or_precedence(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    # "protein or water and ion" should parse as "protein or (water and ion)"
    a = set(sel.select(structure, "protein or water and ion"))
    b = set(sel.select(structure, "protein or (water and ion)"))
    assert a == b


def test_select_parentheses(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    a = set(sel.select(structure, "(protein or water) and chain B"))
    expected_names = {"VAL", "HOH"}
    assert {atom.get_parent().resname for atom in a} <= expected_names


def test_select_within(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    # everything within 3 A of the ZN ion should at least include the ion itself
    atoms = sel.select(structure, "within 3 ion")
    ion_atoms = set(sel.select(structure, "ion"))
    assert ion_atoms <= set(atoms)


def test_select_unknown_keyword_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(SelectionSyntaxError):
        sel.select(structure, "bogus")


def test_select_unbalanced_parens_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(SelectionSyntaxError):
        sel.select(structure, "(protein")


def test_select_returns_canonical_order_regardless_of_file_order():
    # Reproduces a real bug: Bio.PDB's writer groups HETATM records after
    # ATOM records, so a hetero-flagged residue like MSE (selenomethionine,
    # classified as PROTEIN) can sit in a different relative position in
    # the parsed atom list after a save/reload round-trip even though its
    # residue number is unchanged. Two independently-parsed structures
    # covering the same residues must still select() in the same
    # (chain, resid) order for things like geometry.rmsd/cluster's
    # pairwise_rmsd_matrix (which pair atoms positionally) to stay correct.
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Structure import Structure

    def build(order):
        structure = Structure("t")
        model = Model(0)
        chain = Chain("A")
        residues = {
            1: ("ALA", " "),
            2: ("MSE", "H_MSE"),
            3: ("GLY", " "),
        }
        for resseq in order:
            resname, hetflag = residues[resseq]
            residue = Residue((hetflag, resseq, " "), resname, "")
            residue.add(Atom("CA", (float(resseq), 0.0, 0.0), 20.0, 1.0, " ", "CA", resseq, element="C"))
            chain.add(residue)
        model.add(chain)
        structure.add(model)
        return structure

    # "pre-round-trip" order (sequence order) vs "post-round-trip" order
    # (HETATM/MSE relocated to the end, as Bio.PDB's PDBIO would write it)
    structure_a = build([1, 2, 3])
    structure_b = build([1, 3, 2])

    atoms_a = sel.select(structure_a, "protein and atom CA")
    atoms_b = sel.select(structure_b, "protein and atom CA")
    labels_a = [a.get_parent().id[1] for a in atoms_a]
    labels_b = [a.get_parent().id[1] for a in atoms_b]
    assert labels_a == labels_b == [1, 2, 3]


def test_select_empty_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    with pytest.raises(SelectionSyntaxError):
        sel.select(structure, "")
