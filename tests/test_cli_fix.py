from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str) -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_fix_report(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["fix", "report", struct_id])
    assert result.exit_code == 0
    assert "B/VAL1" in result.output
    assert "CB" in result.output


def test_fix_apply_saves_new_structure(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["fix", "apply", struct_id])
    assert result.exit_code == 0
    assert "Saved as 'tiny_fixed'" in result.output

    status = runner.invoke(cli, ["status"])
    assert "tiny_fixed" in status.output


def test_fix_apply_original_unchanged(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    runner.invoke(cli, ["fix", "apply", struct_id])
    info = runner.invoke(cli, ["info", struct_id])
    assert "protein: 4" in info.output


def test_fix_apply_replace_nonstandard(project_dir: Path):
    runner = CliRunner()
    path = Path(__file__).parent / "data" / "nonstandard.pdb"
    struct_id = _import(runner, path, "ns")
    result = runner.invoke(cli, ["fix", "apply", struct_id, "--remove-heterogens", "all"])
    assert result.exit_code == 0
    assert "nonstandard residues replaced: 1" in result.output
    assert "heterogens removed" in result.output


def test_fix_apply_custom_name(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["fix", "apply", struct_id, "--name", "my_fixed"])
    assert result.exit_code == 0
    assert "Saved as 'my_fixed'" in result.output


def test_fix_report_unknown_structure_fails_cleanly(project_dir: Path):
    result = CliRunner().invoke(cli, ["fix", "report", "nope"])
    assert result.exit_code != 0
