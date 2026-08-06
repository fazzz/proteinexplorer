from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def test_predict_colabfold_fails_cleanly_without_binary(project_dir: Path, tmp_path: Path):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["predict", "colabfold", "MKVLTA", "--output-dir", str(tmp_path / "out")],
    )
    assert result.exit_code != 0
    assert "ColabFold" in result.output


def test_predict_alphafold_fails_cleanly_without_script(project_dir: Path, tmp_path: Path):
    runner = CliRunner()
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">q\nMKVLTA\n")
    result = runner.invoke(
        cli,
        [
            "predict", "alphafold",
            "--fasta", str(fasta),
            "--output-dir", str(tmp_path / "out"),
            "--alphafold-script", str(tmp_path / "nope.sh"),
            "--data-dir", str(tmp_path / "data"),
            "--max-template-date", "2020-01-01",
        ],
    )
    assert result.exit_code != 0
    assert "AlphaFold" in result.output
