from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, tiny_pdb: Path, name: str = "tiny") -> str:
    result = runner.invoke(cli, ["import", str(tiny_pdb), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_geometry_distance_atom_atom(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(
        cli,
        ["geometry", "distance", struct_id,
         "chain A and resid 1 and atom CA", "chain A and resid 2 and atom CA"],
    )
    assert result.exit_code == 0
    assert "A" in result.output  # "... A" unit suffix


def test_geometry_distance_invalid_selection(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["geometry", "distance", struct_id, "bogus", "protein"])
    assert result.exit_code != 0
    assert "Invalid selection" in result.output


def test_geometry_distance_empty_selection(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["geometry", "distance", struct_id, "resname NOPE", "protein"])
    assert result.exit_code != 0
    assert "matched no atoms" in result.output


def test_geometry_angle(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(
        cli,
        ["geometry", "angle", struct_id,
         "chain A and resid 1 and atom N",
         "chain A and resid 1 and atom CA",
         "chain A and resid 1 and atom C"],
    )
    assert result.exit_code == 0
    assert "deg" in result.output


def test_geometry_dihedral(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(
        cli,
        ["geometry", "dihedral", struct_id,
         "chain A and resid 1 and atom N",
         "chain A and resid 1 and atom CA",
         "chain A and resid 1 and atom C",
         "chain A and resid 2 and atom N"],
    )
    assert result.exit_code == 0
    assert "deg" in result.output


def test_geometry_backbone_torsions(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["geometry", "backbone-torsions", struct_id, "--chain", "A", "--resid", "2"])
    assert result.exit_code == 0
    assert "phi=" in result.output


def test_geometry_backbone_torsions_unknown_chain(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["geometry", "backbone-torsions", struct_id, "--chain", "Z", "--resid", "2"])
    assert result.exit_code != 0


def test_geometry_rmsd_self_is_zero(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(
        cli, ["geometry", "rmsd", struct_id, struct_id, "--selection", "chain A and backbone"]
    )
    assert result.exit_code == 0
    assert "0.000" in result.output


def test_geometry_coords(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["geometry", "coords", struct_id, "chain A"])
    assert result.exit_code == 0
    assert "centroid" in result.output
    assert "principal axes" in result.output


def test_geometry_distmatrix(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(
        cli,
        ["geometry", "distmatrix", struct_id,
         "chain A and resid 1 and atom CA",
         "chain A and resid 2 and atom CA",
         "chain A and resid 3 and atom CA"],
    )
    assert result.exit_code == 0
    assert "[1]" in result.output and "[3]" in result.output


def test_geometry_distmatrix_needs_two_selections(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["geometry", "distmatrix", struct_id, "chain A"])
    assert result.exit_code != 0
