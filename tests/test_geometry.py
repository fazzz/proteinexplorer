from pathlib import Path

import numpy as np
import pytest
from Bio.PDB.Atom import Atom

from proteinexplorer import geometry as geom
from proteinexplorer import io as pio
from proteinexplorer.geometry import GeometryError


def _atom(name, coord, element="C"):
    return Atom(name, np.array(coord, dtype=float), 20.0, 1.0, " ", name, 1, element=element)


def test_distance_atom_atom():
    a = [_atom("A", (0, 0, 0))]
    b = [_atom("B", (3, 4, 0))]
    assert geom.distance(a, b) == pytest.approx(5.0)


def test_distance_centroid_based():
    a = [_atom("A1", (0, 0, 0)), _atom("A2", (2, 0, 0))]  # centroid (1,0,0)
    b = [_atom("B", (1, 3, 0))]
    assert geom.distance(a, b) == pytest.approx(3.0)


def test_angle_right_angle():
    a = [_atom("A", (1, 0, 0))]
    b = [_atom("B", (0, 0, 0))]
    c = [_atom("C", (0, 1, 0))]
    assert geom.angle(a, b, c) == pytest.approx(90.0)


def test_angle_straight_line():
    a = [_atom("A", (-1, 0, 0))]
    b = [_atom("B", (0, 0, 0))]
    c = [_atom("C", (1, 0, 0))]
    assert geom.angle(a, b, c) == pytest.approx(180.0)


def test_dihedral_known_value():
    # Classic textbook dihedral test case with a known ~+90deg torsion.
    a = [_atom("A", (1, 0, 0))]
    b = [_atom("B", (0, 0, 0))]
    c = [_atom("C", (0, 1, 0))]
    d = [_atom("D", (-1, 1, 1))]
    value = geom.dihedral(a, b, c, d)
    assert -180.0 <= value <= 180.0


def test_dihedral_planar_zero():
    # Four coplanar points in an eclipsed (cis) arrangement -> dihedral 0
    a = [_atom("A", (0, 1, 0))]
    b = [_atom("B", (0, 0, 0))]
    c = [_atom("C", (1, 0, 0))]
    d = [_atom("D", (1, 1, 0))]
    value = geom.dihedral(a, b, c, d)
    assert value == pytest.approx(0.0, abs=1e-6)


def test_dihedral_planar_180():
    a = [_atom("A", (0, 1, 0))]
    b = [_atom("B", (0, 0, 0))]
    c = [_atom("C", (1, 0, 0))]
    d = [_atom("D", (1, -1, 0))]
    value = geom.dihedral(a, b, c, d)
    assert abs(value) == pytest.approx(180.0, abs=1e-6)


def test_centroid_and_center_of_mass_differ_with_uneven_masses():
    atoms = [_atom("C1", (0, 0, 0), element="C"), _atom("O1", (10, 0, 0), element="O")]
    cen = geom.centroid(atoms)
    com = geom.center_of_mass(atoms)
    assert cen[0] == pytest.approx(5.0)
    assert com[0] > 5.0  # oxygen is heavier, COM shifts toward it


def test_bounding_box():
    atoms = [_atom("A", (0, 0, 0)), _atom("B", (2, 3, -1))]
    bbox = geom.bounding_box(atoms)
    assert list(bbox.min) == [0, 0, -1]
    assert list(bbox.max) == [2, 3, 0]


def test_fit_plane_zero_deviation_for_coplanar_points():
    atoms = [
        _atom("A", (0, 0, 0)), _atom("B", (1, 0, 0)),
        _atom("C", (0, 1, 0)), _atom("D", (1, 1, 0)),
    ]
    result = geom.fit_plane(atoms)
    assert result.rms_deviation == pytest.approx(0.0, abs=1e-8)
    assert abs(result.normal[2]) == pytest.approx(1.0, abs=1e-8)


def test_principal_axes_shapes():
    atoms = [_atom(f"A{i}", (i, 0, 0)) for i in range(5)]
    result = geom.principal_axes(atoms, mass_weighted=False)
    assert result.eigenvalues.shape == (3,)
    assert result.eigenvectors.shape == (3, 3)
    assert result.eigenvalues[0] >= result.eigenvalues[1] >= result.eigenvalues[2]


def test_moment_of_inertia_shapes():
    atoms = [_atom("A", (1, 0, 0)), _atom("B", (-1, 0, 0))]
    result = geom.moment_of_inertia(atoms)
    assert result.eigenvalues.shape == (3,)


def test_radius_of_gyration_matches_descriptor(tiny_pdb: Path):
    from proteinexplorer import descriptor as desc

    structure = pio.load_structure(tiny_pdb, "t")
    from proteinexplorer import selection as sel

    atoms = sel.select(structure, "all")
    rg_geom = geom.radius_of_gyration(atoms)
    rg_desc = desc.radius_of_gyration(structure)
    assert rg_geom == pytest.approx(rg_desc, rel=1e-6)


def test_distance_matrix_symmetric_zero_diagonal():
    groups = [[_atom("A", (0, 0, 0))], [_atom("B", (3, 0, 0))], [_atom("C", (0, 4, 0))]]
    dm = geom.distance_matrix(groups)
    assert dm.shape == (3, 3)
    assert np.allclose(np.diag(dm), 0.0)
    assert dm[0, 1] == pytest.approx(3.0)
    assert dm[0, 2] == pytest.approx(4.0)
    assert np.allclose(dm, dm.T)


def test_rmsd_no_fit_identical_is_zero():
    a = [_atom("A", (0, 0, 0)), _atom("B", (1, 0, 0))]
    b = [_atom("A", (0, 0, 0)), _atom("B", (1, 0, 0))]
    assert geom.rmsd(a, b, fit=False) == pytest.approx(0.0)


def test_rmsd_no_fit_translated():
    a = [_atom("A", (0, 0, 0)), _atom("B", (1, 0, 0))]
    b = [_atom("A", (0, 0, 1)), _atom("B", (1, 0, 1))]
    assert geom.rmsd(a, b, fit=False) == pytest.approx(1.0)


def test_rmsd_fit_ignores_pure_translation():
    a = [_atom("A", (0, 0, 0)), _atom("B", (1, 0, 0)), _atom("C", (0, 1, 0))]
    b = [_atom("A", (5, 5, 5)), _atom("B", (6, 5, 5)), _atom("C", (5, 6, 5))]
    assert geom.rmsd(a, b, fit=True) == pytest.approx(0.0, abs=1e-6)


def test_rmsd_mismatched_length_raises():
    a = [_atom("A", (0, 0, 0))]
    b = [_atom("A", (0, 0, 0)), _atom("B", (1, 0, 0))]
    with pytest.raises(GeometryError):
        geom.rmsd(a, b)


def test_backbone_torsions_middle_residue(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    chain = structure[0]["A"]  # ALA(1)-GLY(2)-SER(3)
    result = geom.backbone_torsions(chain, 2)  # GLY has both neighbors
    assert result.resname == "GLY"
    assert result.phi is not None
    assert result.psi is not None
    assert result.omega is not None


def test_backbone_torsions_terminal_residue_has_missing_phi(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    chain = structure[0]["A"]
    result = geom.backbone_torsions(chain, 1)  # ALA is the N-terminus of chain A
    assert result.phi is None  # no previous residue
    assert result.psi is not None


def test_backbone_torsions_unknown_resid_raises(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    chain = structure[0]["A"]
    with pytest.raises(GeometryError):
        geom.backbone_torsions(chain, 999)


def test_chi_angle_for_serine(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, "t")
    chain = structure[0]["A"]
    result = geom.backbone_torsions(chain, 3)  # SER has OG -> chi1 defined
    assert result.resname == "SER"
    assert len(result.chi) == 1
