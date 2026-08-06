from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str = "tiny") -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_plot_ramachandran(project_dir: Path, tiny_pdb: Path, tmp_path: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    out = tmp_path / "rama.png"
    result = runner.invoke(cli, ["plot", "ramachandran", struct_id, str(out)])
    assert result.exit_code == 0
    assert out.exists()


def test_plot_contact_map(project_dir: Path, tiny_pdb: Path, tmp_path: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    out = tmp_path / "contact.png"
    result = runner.invoke(cli, ["plot", "contact-map", struct_id, str(out)])
    assert result.exit_code == 0
    assert out.exists()


def test_plot_secondary(project_dir: Path, tiny_pdb: Path, tmp_path: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb)
    out = tmp_path / "ss.png"
    result = runner.invoke(cli, ["plot", "secondary", struct_id, str(out), "--method", "geometric"])
    assert result.exit_code == 0
    assert out.exists()
