from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def test_replay_cli_dry_run(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    result = runner.invoke(cli, ["replay", "--dry-run"])
    assert result.exit_code == 0
    assert "[1]" in result.output


def test_replay_cli_reconstructs_project(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    runner.invoke(cli, ["mutate", "tiny", "--chain", "A", "--resid", "1", "--to", "VAL"])

    result = runner.invoke(cli, ["replay"])
    assert result.exit_code == 0
    assert "Backed up previous project state" in result.output
    assert "2 step(s), 0 failed" in result.output

    status = runner.invoke(cli, ["status"])
    assert "tiny_A1ALAVAL" in status.output


def test_replay_cli_no_logged_commands_fails_cleanly(project_dir: Path):
    result = CliRunner().invoke(cli, ["replay"])
    assert result.exit_code != 0
