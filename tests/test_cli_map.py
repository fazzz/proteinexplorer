from pathlib import Path

from click.testing import CliRunner

from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str) -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_map_pocket(project_dir: Path, tmp_path: Path):
    runner = CliRunner()
    cage_path = Path(__file__).parent / "data" / "cage.pdb"
    struct_id = _import(runner, cage_path, "cage")
    out = tmp_path / "pocket.pml"
    result = runner.invoke(
        cli, ["map", "pocket", struct_id, str(out), "--spacing", "1.0", "--padding", "1.0",
              "--min-pocket-points", "1"]
    )
    assert result.exit_code == 0
    assert out.exists()
    assert "color" in out.read_text()


def test_map_mutation(project_dir: Path, tiny_pdb: Path, tmp_path: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    out = tmp_path / "mut.pml"
    result = runner.invoke(
        cli, ["map", "mutation", struct_id, str(out), "--residue", "A/1", "--residue", "A/2"]
    )
    assert result.exit_code == 0
    assert out.exists()
    assert "color red" in out.read_text()


def test_map_domain(project_dir: Path, tiny_pdb: Path, tmp_path: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    out = tmp_path / "domain.pml"
    result = runner.invoke(
        cli, ["map", "domain", struct_id, str(out), "--range", "A:1-2:testdomain"]
    )
    assert result.exit_code == 0
    assert out.exists()


def test_map_domain_invalid_range_fails_cleanly(project_dir: Path, tiny_pdb: Path, tmp_path: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    out = tmp_path / "domain.pml"
    result = runner.invoke(cli, ["map", "domain", struct_id, str(out), "--range", "bogus"])
    assert result.exit_code != 0


def test_map_conservation(project_dir: Path, tiny_pdb: Path, tmp_path: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    csv_path = tmp_path / "values.csv"
    csv_path.write_text("A/1,0.9\nA/2,0.1\nA/3,0.5\n")
    out = tmp_path / "cons.pml"
    result = runner.invoke(cli, ["map", "conservation", struct_id, str(csv_path), str(out)])
    assert result.exit_code == 0
    assert out.exists()
    assert "spectrum b" in out.read_text()

    status = runner.invoke(cli, ["status"])
    assert "tiny_bfactor" in status.output
