from pathlib import Path

from proteinexplorer import descriptor as desc
from proteinexplorer import io as pio


def _load(tiny_pdb: Path):
    structure = pio.load_structure(tiny_pdb, structure_id="tiny")
    summary = pio.summarize(structure)
    return structure, summary.category_totals()


def test_molecular_weight_positive(tiny_pdb: Path):
    structure, _ = _load(tiny_pdb)
    mw = desc.molecular_weight(structure)
    assert mw > 0


def test_sasa_total_positive(tiny_pdb: Path):
    structure, _ = _load(tiny_pdb)
    sasa = desc.sasa_total(structure)
    assert sasa > 0


def test_radius_of_gyration_positive(tiny_pdb: Path):
    structure, _ = _load(tiny_pdb)
    rg = desc.radius_of_gyration(structure)
    assert rg is not None
    assert rg > 0


def test_hydrophobic_ratio_range(tiny_pdb: Path):
    structure, _ = _load(tiny_pdb)
    ratio = desc.hydrophobic_ratio(structure)
    assert ratio is not None
    assert 0.0 <= ratio <= 1.0
    # chain A: ALA, GLY, SER; chain B: VAL -> 3 of 4 protein residues hydrophobic
    assert ratio == 3 / 4


def test_disulfide_count_zero_when_no_cys(tiny_pdb: Path):
    structure, _ = _load(tiny_pdb)
    assert desc.disulfide_count(structure) == 0


def test_contact_density_reasonable(tiny_pdb: Path):
    structure, _ = _load(tiny_pdb)
    density = desc.contact_density(structure)
    assert density is not None
    assert density >= 0


def test_secondary_structure_falls_back_to_geometric_without_dssp(tiny_pdb: Path, monkeypatch):
    structure, _ = _load(tiny_pdb)
    monkeypatch.setattr("shutil.which", lambda name: None)
    from proteinexplorer import secondary as sec

    residues, method = sec.secondary_structure(structure, pdb_path=tiny_pdb, method="auto")
    assert method == "geometric"
    assert len(residues) == 4


def test_compute_descriptors_end_to_end(tiny_pdb: Path):
    structure, totals = _load(tiny_pdb)
    result = desc.compute_descriptors(structure, tiny_pdb, totals)
    assert result.n_atoms == 24
    assert result.n_ligands == 1
    assert result.n_waters == 2
    # DSSP isn't installed in the test environment -> auto falls back to
    # the geometric classifier, so a composition is still produced.
    assert result.secondary_structure is not None
    assert result.secondary_structure_method == "geometric"
    assert result.secondary_structure_error is None
