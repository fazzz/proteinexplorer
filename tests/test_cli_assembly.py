from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str) -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def _a1a8o_path() -> Path:
    return Path(__file__).parent.parent / "examples" / "1a8o" / "1A8O.pdb"


def test_assembly_list(project_dir: Path):
    path = _a1a8o_path()
    if not path.exists():
        import pytest
        pytest.skip("examples/1a8o/1A8O.pdb not present")
    runner = CliRunner()
    struct_id = _import(runner, path, "1a8o")
    result = runner.invoke(cli, ["assembly", "list", struct_id])
    assert result.exit_code == 0
    assert "DIMERIC" in result.output


def test_assembly_list_none_documented(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["assembly", "list", struct_id])
    assert result.exit_code == 0
    assert "No biological assembly documented" in result.output


def test_assembly_generate_saves_new_structure(project_dir: Path):
    path = _a1a8o_path()
    if not path.exists():
        import pytest
        pytest.skip("examples/1a8o/1A8O.pdb not present")
    runner = CliRunner()
    struct_id = _import(runner, path, "1a8o")
    result = runner.invoke(cli, ["assembly", "generate", struct_id])
    assert result.exit_code == 0
    assert "['A1', 'A2']" in result.output
    assert "Saved as '1a8o_assembly'" in result.output

    status = runner.invoke(cli, ["status"])
    assert "1a8o_assembly" in status.output


def test_assembly_generate_no_assembly_fails_cleanly(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["assembly", "generate", struct_id])
    assert result.exit_code != 0
    assert "no documented" in result.output
