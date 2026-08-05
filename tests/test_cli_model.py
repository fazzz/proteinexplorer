from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str = "gapped") -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def _gapped_pdb() -> Path:
    return Path(__file__).parent / "data" / "gapped.pdb"


def test_model_gaps_detects_gap(project_dir: Path):
    runner = CliRunner()
    struct_id = _import(runner, _gapped_pdb())
    result = runner.invoke(cli, ["model", "gaps", struct_id])
    assert result.exit_code == 0
    assert "2 .. 6" in result.output
    assert "3 residue(s) missing" in result.output


def test_model_gaps_none_for_contiguous_structure(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, name="tiny")
    result = runner.invoke(cli, ["model", "gaps", struct_id])
    assert result.exit_code == 0
    assert "No gaps detected" in result.output


def test_model_loop_fills_gap_and_saves_new_structure(project_dir: Path):
    runner = CliRunner()
    struct_id = _import(runner, _gapped_pdb())
    result = runner.invoke(
        cli, ["model", "loop", struct_id, "--chain", "A", "--start", "3", "--end", "5", "--sequence", "GVL"]
    )
    assert result.exit_code == 0
    assert "GLY3, VAL4, LEU5" in result.output
    assert "Saved as" in result.output

    status = runner.invoke(cli, ["status"])
    assert "gapped_loopA3-5" in status.output


def test_model_loop_default_sequence_is_poly_ala(project_dir: Path):
    runner = CliRunner()
    struct_id = _import(runner, _gapped_pdb())
    result = runner.invoke(cli, ["model", "loop", struct_id, "--chain", "A", "--start", "3", "--end", "5"])
    assert result.exit_code == 0
    assert "ALA3, ALA4, ALA5" in result.output


def test_model_loop_missing_anchor_fails_cleanly(project_dir: Path):
    runner = CliRunner()
    struct_id = _import(runner, _gapped_pdb())
    result = runner.invoke(cli, ["model", "loop", struct_id, "--chain", "A", "--start", "20", "--end", "22"])
    assert result.exit_code != 0


def test_model_homology_fails_cleanly_without_modeller(project_dir: Path, tmp_path: Path):
    runner = CliRunner()
    alignment = tmp_path / "align.pir"
    alignment.write_text(">P1;dummy\n")
    result = runner.invoke(
        cli,
        [
            "model", "homology",
            "--alignment", str(alignment),
            "--template", "tmpl",
            "--target", "tgt",
            "--template-dir", str(tmp_path),
            "--output-dir", str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
    assert "MODELLER" in result.output
