from proteinexplorer.models import ResidueCategory, classify_residue, is_backbone_atom


def test_classify_standard_amino_acid():
    assert classify_residue("ALA", " ") is ResidueCategory.PROTEIN


def test_classify_nonstandard_amino_acid():
    assert classify_residue("MSE", "H_MSE") is ResidueCategory.PROTEIN


def test_classify_water_by_hetero_flag():
    assert classify_residue("HOH", "W") is ResidueCategory.WATER


def test_classify_water_by_name_even_without_w_flag():
    assert classify_residue("HOH", "H_HOH") is ResidueCategory.WATER


def test_classify_nucleic():
    assert classify_residue("DA", "H_DA") is ResidueCategory.NUCLEIC
    assert classify_residue("G", "H_G") is ResidueCategory.NUCLEIC


def test_classify_ion():
    assert classify_residue("ZN", "H_ZN") is ResidueCategory.ION


def test_classify_unknown_hetero_is_ligand():
    assert classify_residue("LIG", "H_LIG") is ResidueCategory.LIGAND


def test_backbone_atom_protein():
    assert is_backbone_atom("CA", ResidueCategory.PROTEIN)
    assert not is_backbone_atom("CB", ResidueCategory.PROTEIN)


def test_backbone_atom_nucleic():
    assert is_backbone_atom("P", ResidueCategory.NUCLEIC)
    assert not is_backbone_atom("N1", ResidueCategory.NUCLEIC)


def test_backbone_atom_not_applicable_to_ligand():
    assert not is_backbone_atom("CA", ResidueCategory.LIGAND)
