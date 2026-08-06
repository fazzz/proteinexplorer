from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str) -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_cluster_ensemble_identical_structures_form_one_cluster(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    id1 = _import(runner, tiny_pdb, "s1")
    id2 = _import(runner, tiny_pdb, "s2")
    result = runner.invoke(cli, ["cluster", "ensemble", id1, id2])
    assert result.exit_code == 0
    assert "1 cluster(s)" in result.output


def test_cluster_ensemble_hierarchical(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    id1 = _import(runner, tiny_pdb, "s1")
    id2 = _import(runner, tiny_pdb, "s2")
    result = runner.invoke(cli, ["cluster", "ensemble", id1, id2, "--method", "hierarchical", "--n-clusters", "1"])
    assert result.exit_code == 0
    assert "method=hierarchical" in result.output


def test_cluster_ensemble_needs_two_ids(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    id1 = _import(runner, tiny_pdb, "s1")
    result = runner.invoke(cli, ["cluster", "ensemble", id1])
    assert result.exit_code != 0


def test_cluster_models_multimodel_structure(project_dir: Path):
    runner = CliRunner()
    path = Path(__file__).parent / "data" / "tiny.pdb"
    # tiny.pdb is single-model; just verify the clean "only 1 model" error path
    struct_id = _import(runner, path, "tiny")
    result = runner.invoke(cli, ["cluster", "models", struct_id])
    assert result.exit_code != 0
    assert "only 1 model" in result.output


def test_cluster_models_separates_folds(project_dir: Path):
    runner = CliRunner()
    path = Path(__file__).parent / "data" / "multimodel.pdb"
    struct_id = _import(runner, path, "ensemble")
    result = runner.invoke(cli, ["cluster", "models", struct_id, "--threshold", "1.0"])
    assert result.exit_code == 0
    assert "2 cluster(s)" in result.output
