from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str) -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_valid_clashes_none_found(project_dir: Path):
    runner = CliRunner()
    path = Path(__file__).parent.parent / "examples" / "1a8o" / "1A8O.pdb"
    if not path.exists():
        import pytest
        pytest.skip("examples/1a8o/1A8O.pdb not present")
    struct_id = _import(runner, path, "1a8o")
    result = runner.invoke(cli, ["valid", "clashes", struct_id])
    assert result.exit_code == 0
    assert "No clashes found" in result.output


def test_valid_geometry_none_found(project_dir: Path):
    runner = CliRunner()
    path = Path(__file__).parent.parent / "examples" / "1a8o" / "1A8O.pdb"
    if not path.exists():
        import pytest
        pytest.skip("examples/1a8o/1A8O.pdb not present")
    struct_id = _import(runner, path, "1a8o")
    result = runner.invoke(cli, ["valid", "geometry", struct_id])
    assert result.exit_code == 0
    assert "No bond geometry outliers found" in result.output


def test_valid_molprobity_fails_cleanly_without_binary(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["valid", "molprobity", struct_id])
    assert result.exit_code != 0
    assert "Phenix" in result.output
