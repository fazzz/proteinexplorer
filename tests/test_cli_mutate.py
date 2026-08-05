from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str = "tiny") -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_mutate_creates_new_structure(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)

    result = runner.invoke(cli, ["mutate", struct_id, "--chain", "A", "--resid", "1", "--to", "VAL"])
    assert result.exit_code == 0
    assert "ALA1 -> VAL" in result.output
    assert "Saved as" in result.output

    status = runner.invoke(cli, ["status"])
    assert "tiny_A1ALAVAL" in status.output


def test_mutate_original_structure_unchanged(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    runner.invoke(cli, ["mutate", struct_id, "--chain", "A", "--resid", "1", "--to", "VAL"])

    info = runner.invoke(cli, ["info", struct_id])
    assert "protein: 4" in info.output  # original still has its 4 protein residues intact


def test_mutate_custom_name(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(
        cli, ["mutate", struct_id, "--chain", "A", "--resid", "1", "--to", "VAL", "--name", "my_mutant"]
    )
    assert result.exit_code == 0
    assert "Saved as 'my_mutant'" in result.output


def test_mutate_explicit_scwrl4_fails_cleanly(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(
        cli, ["mutate", struct_id, "--chain", "A", "--resid", "1", "--to", "VAL", "--method", "scwrl4"]
    )
    assert result.exit_code != 0
    assert "Scwrl4" in result.output


def test_mutate_unknown_residue_fails_cleanly(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["mutate", struct_id, "--chain", "A", "--resid", "999", "--to", "VAL"])
    assert result.exit_code != 0


def test_mutate_unknown_target_fails_cleanly(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    result = runner.invoke(cli, ["mutate", struct_id, "--chain", "A", "--resid", "1", "--to", "XYZ"])
    assert result.exit_code != 0
