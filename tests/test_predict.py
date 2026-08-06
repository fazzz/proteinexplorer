from pathlib import Path

import pytest

from proteinexplorer import predict as pred


# --- write_fasta -----------------------------------------------------

def test_write_fasta_format(tmp_path: Path):
    path = pred.write_fasta("MKV LTA", "query1", tmp_path / "out.fasta")
    content = path.read_text()
    assert content == ">query1\nMKVLTA\n"


def test_write_fasta_strips_whitespace_and_newlines(tmp_path: Path):
    path = pred.write_fasta("MK\nVL\tTA ", "q", tmp_path / "out.fasta")
    assert path.read_text() == ">q\nMKVLTA\n"


# --- ColabFold (binary not installed in this environment) ----------------

def test_colabfold_binary_not_found():
    assert pred.colabfold_binary() is None


def test_colabfold_predict_raises_clean_error_without_binary(tmp_path: Path):
    with pytest.raises(pred.PredictionToolNotAvailableError):
        pred.colabfold_predict("MKVLTA", tmp_path)


def test_colabfold_predict_error_mentions_install_pointer(tmp_path: Path):
    with pytest.raises(pred.PredictionToolNotAvailableError, match="ColabFold"):
        pred.colabfold_predict("MKVLTA", tmp_path)


# --- AlphaFold (script not installed in this environment) ---------------

def test_alphafold_predict_raises_clean_error_without_script(tmp_path: Path):
    with pytest.raises(pred.PredictionToolNotAvailableError):
        pred.alphafold_predict(
            fasta_path=tmp_path / "in.fasta",
            output_dir=tmp_path / "out",
            alphafold_script=tmp_path / "does_not_exist" / "run_alphafold.sh",
            data_dir=tmp_path / "data",
            max_template_date="2020-01-01",
        )


def test_alphafold_predict_error_mentions_install_pointer(tmp_path: Path):
    with pytest.raises(pred.PredictionToolNotAvailableError, match="AlphaFold"):
        pred.alphafold_predict(
            fasta_path=tmp_path / "in.fasta",
            output_dir=tmp_path / "out",
            alphafold_script=tmp_path / "nope.sh",
            data_dir=tmp_path / "data",
            max_template_date="2020-01-01",
        )
