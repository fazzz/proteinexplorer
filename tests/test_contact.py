from pathlib import Path

import pytest

from proteinexplorer import contact as ct
from proteinexplorer import io as pio


@pytest.fixture()
def contacts_structure():
    path = Path(__file__).parent / "data" / "contacts.pdb"
    return pio.load_structure(path, structure_id="contacts")


def test_disulfide_bonds(contacts_structure):
    bonds = ct.find_disulfide_bonds(contacts_structure)
    assert len(bonds) == 1
    bond = bonds[0]
    assert {bond.residue_a, bond.residue_b} == {"A/CYS1", "A/CYS2"}
    assert bond.distance == pytest.approx(2.05, abs=0.01)


def test_disulfide_bonds_respects_cutoff(contacts_structure):
    assert ct.find_disulfide_bonds(contacts_structure, cutoff=1.0) == []


def test_salt_bridges(contacts_structure):
    bridges = ct.find_salt_bridges(contacts_structure)
    assert len(bridges) == 1
    b = bridges[0]
    assert b.basic_residue == "A/ARG10"
    assert b.acidic_residue == "A/ASP11"
    assert b.distance == pytest.approx(3.67, abs=0.05)


def test_salt_bridges_respects_cutoff(contacts_structure):
    assert ct.find_salt_bridges(contacts_structure, cutoff=1.0) == []


def test_hydrophobic_contacts_includes_leu_val(contacts_structure):
    contacts = ct.find_hydrophobic_contacts(contacts_structure)
    pairs = {frozenset((c.residue_a, c.residue_b)) for c in contacts}
    assert frozenset({"A/LEU40", "A/VAL41"}) in pairs


def test_hydrophobic_contacts_excludes_charged_residues(contacts_structure):
    contacts = ct.find_hydrophobic_contacts(contacts_structure)
    for c in contacts:
        assert "ARG" not in c.residue_a and "ARG" not in c.residue_b
        assert "ASP" not in c.residue_a and "ASP" not in c.residue_b


def test_pipi_interaction_parallel_stack(contacts_structure):
    interactions = ct.find_pipi_interactions(contacts_structure)
    assert len(interactions) == 1
    p = interactions[0]
    assert {p.residue_a, p.residue_b} == {"A/PHE30", "A/PHE31"}
    assert p.centroid_distance == pytest.approx(4.5, abs=0.01)
    assert p.stack_type == "parallel"
    assert p.plane_angle == pytest.approx(0.0, abs=1.0)


def test_pipi_respects_cutoff(contacts_structure):
    assert ct.find_pipi_interactions(contacts_structure, cutoff=1.0) == []


def test_cationpi_interaction(contacts_structure):
    interactions = ct.find_cationpi_interactions(contacts_structure)
    assert len(interactions) == 1
    c = interactions[0]
    assert c.cation_residue == "A/LYS20"
    assert c.aromatic_residue == "A/PHE21"
    assert c.distance == pytest.approx(5.0, abs=0.05)


def test_cationpi_respects_cutoff(contacts_structure):
    assert ct.find_cationpi_interactions(contacts_structure, cutoff=1.0) == []


def test_hydrogen_bonds_found(contacts_structure):
    bonds = ct.find_hydrogen_bonds(contacts_structure)
    pairs = {frozenset((b.donor_residue, b.acceptor_residue)) for b in bonds}
    assert frozenset({"A/ALA50", "A/GLY51"}) in pairs
    assert frozenset({"A/ARG10", "A/ASP11"}) in pairs
    for b in bonds:
        assert b.distance <= 3.5


def test_hydrogen_bonds_exclude_same_residue(contacts_structure):
    bonds = ct.find_hydrogen_bonds(contacts_structure)
    for b in bonds:
        assert b.donor_residue != b.acceptor_residue


def test_contact_map_ca_mode(contacts_structure):
    cm = ct.contact_map(contacts_structure, mode="ca", cutoff=8.0)
    assert len(cm.labels) == 12  # 12 protein residues in the fixture
    assert cm.matrix.shape == (12, 12)
    assert not cm.matrix.diagonal().any()
    assert (cm.matrix == cm.matrix.T).all()


def test_contact_map_heavy_mode_more_sensitive_than_ca(contacts_structure):
    cm_ca = ct.contact_map(contacts_structure, mode="ca", cutoff=6.0)
    cm_heavy = ct.contact_map(contacts_structure, mode="heavy", cutoff=6.0)
    # heavy-atom min-distance contacts should be >= CA-CA contacts at the same cutoff
    assert cm_heavy.matrix.sum() >= cm_ca.matrix.sum()


def test_interaction_network_aggregates_all_types(contacts_structure):
    edges = ct.interaction_network(contacts_structure)
    kinds = {e.kind for e in edges}
    assert "disulfide" in kinds
    assert "salt_bridge" in kinds
    assert "hydrophobic" in kinds
    assert any(k.startswith("pi_pi") for k in kinds)
    assert "cation_pi" in kinds
    assert "hbond" in kinds


def test_interaction_network_edges_have_positive_distance(contacts_structure):
    edges = ct.interaction_network(contacts_structure)
    assert all(e.value > 0 for e in edges)


def test_no_contacts_on_isolated_tiny_structure(tiny_pdb: Path):
    # tiny.pdb's residues are all far apart / lack the relevant sidechain
    # atoms, so none of the specialized interaction finders should error.
    structure = pio.load_structure(tiny_pdb, "t")
    assert ct.find_disulfide_bonds(structure) == []
    assert ct.find_pipi_interactions(structure) == []
    assert ct.find_cationpi_interactions(structure) == []
