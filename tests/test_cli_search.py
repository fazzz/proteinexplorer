from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str) -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_search_foldseek_fails_cleanly_without_binary(project_dir: Path, tiny_pdb: Path, tmp_path: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["search", "foldseek", struct_id, "--target-dir", str(tmp_path)])
    assert result.exit_code != 0
    assert "foldseek" in result.output


def test_search_foldseek_requires_exactly_one_target(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["search", "foldseek", struct_id])
    assert result.exit_code != 0
    assert "exactly one" in result.output


def test_search_foldseek_against_project_no_other_structures(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["search", "foldseek", struct_id, "--against-project"])
    assert result.exit_code != 0
    assert "No other structures" in result.output


def test_search_foldseek_against_project_reaches_binary_check(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    id1 = _import(runner, tiny_pdb, "s1")
    _import(runner, tiny_pdb, "s2")
    result = runner.invoke(cli, ["search", "foldseek", id1, "--against-project"])
    assert result.exit_code != 0
    assert "foldseek" in result.output


def test_search_createdb_fails_cleanly_without_binary(project_dir: Path, tmp_path: Path):
    runner = CliRunner()
    (tmp_path / "structures").mkdir()
    result = runner.invoke(cli, ["search", "createdb", str(tmp_path / "structures"), str(tmp_path / "db" / "targetDB")])
    assert result.exit_code != 0
    assert "foldseek" in result.output
