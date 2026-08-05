from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str = "s") -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_secondary_default_auto_falls_back_to_geometric(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["secondary", struct_id])
    assert result.exit_code == 0
    assert "method=geometric" in result.output
    assert "composition:" in result.output


def test_secondary_explicit_geometric(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["secondary", struct_id, "--method", "geometric"])
    assert result.exit_code == 0
    assert "method=geometric" in result.output


def test_secondary_explicit_dssp_fails_cleanly_without_binary(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["secondary", struct_id, "--method", "dssp"])
    assert result.exit_code != 0
    assert "DSSP" in result.output


def test_secondary_chain_filter(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["secondary", struct_id, "--chain", "A"])
    assert result.exit_code == 0
    assert "A:" in result.output
    assert "B:" not in result.output


def test_secondary_unknown_chain(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["secondary", struct_id, "--chain", "Z"])
    assert result.exit_code != 0
