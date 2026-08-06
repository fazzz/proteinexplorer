from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def test_view_fails_cleanly_without_binary(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    struct_id = result.output.splitlines()[0].split()[-1]

    result = runner.invoke(cli, ["view", struct_id])
    assert result.exit_code != 0
    assert "pymol" in result.output


def test_view_unknown_structure_fails_cleanly(project_dir: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["view", "nope"])
    assert result.exit_code != 0
