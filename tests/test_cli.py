from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def test_status_before_import(project_dir: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["status"])
    assert result.exit_code == 0
    assert "No structures imported yet" in result.output


def test_import_then_status(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    assert result.exit_code == 0
    assert "Imported 'tiny'" in result.output

    result = runner.invoke(cli, ["status"])
    assert result.exit_code == 0
    assert "tiny" in result.output
    assert "chains=2" in result.output


def test_export_roundtrip(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    struct_id = result.output.splitlines()[0].split()[-1].rstrip(".")

    out_path = project_dir / "out.cif"
    result = runner.invoke(cli, ["export", struct_id, str(out_path)])
    assert result.exit_code == 0
    assert out_path.exists()


def test_export_unknown_id_is_clean_error(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    result = runner.invoke(cli, ["export", "nonexistent", "out.pdb"])
    assert result.exit_code != 0
    assert "No structure with id or name" in result.output


def test_import_nonexistent_file(project_dir: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["import", "/no/such/file.pdb"])
    assert result.exit_code != 0


def test_info_command(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    struct_id = result.output.splitlines()[0].split()[-1]

    result = runner.invoke(cli, ["info", struct_id])
    assert result.exit_code == 0
    assert "tiny" in result.output
    assert "chains: 2" in result.output
    assert "hetero groups: LIG, ZN" in result.output


def test_descriptor_command(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    struct_id = result.output.splitlines()[0].split()[-1]

    result = runner.invoke(cli, ["descriptor", struct_id])
    assert result.exit_code == 0
    assert "molecular weight" in result.output
    assert "SASA" in result.output
    assert "secondary structure: unavailable" in result.output


def test_info_unknown_id(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    runner.invoke(cli, ["import", str(tiny_pdb), "--name", "tiny"])
    result = runner.invoke(cli, ["info", "nope"])
    assert result.exit_code != 0
