from pathlib import Path

from click.testing import CliRunner

from proteinexplorer import annotate as ann
from proteinexplorer.cli import cli


def _import(runner: CliRunner, path: Path, name: str) -> str:
    result = runner.invoke(cli, ["import", str(path), "--name", name])
    assert result.exit_code == 0, result.output
    return result.output.splitlines()[0].split()[-1]


def test_annotate_metal_sites(project_dir: Path):
    runner = CliRunner()
    path = Path(__file__).parent / "data" / "metal_site.pdb"
    struct_id = _import(runner, path, "metal")
    result = runner.invoke(cli, ["annotate", "metal-sites", struct_id])
    assert result.exit_code == 0
    assert "A/ZN100" in result.output
    assert "A/HIS1" in result.output


def test_annotate_metal_sites_none_found(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["annotate", "metal-sites", struct_id, "--cutoff", "0.01"])
    assert result.exit_code == 0
    assert "No metal-binding sites found" in result.output


def test_annotate_metadata(project_dir: Path, tiny_pdb: Path):
    runner = CliRunner()
    struct_id = _import(runner, tiny_pdb, "tiny")
    result = runner.invoke(cli, ["annotate", "metadata", struct_id])
    assert result.exit_code == 0
    assert "method:" in result.output
    assert "resolution:" in result.output


def test_annotate_uniprot_success(project_dir: Path, monkeypatch):
    def fake_lookup(accession):
        return ann.UniProtAnnotation(
            accession=accession, gene_names=["lysC"], organism="Gallus gallus",
            taxonomy_id=9031, ec_numbers=["3.2.1.17"], go_terms=["GO:0003796"],
        )

    monkeypatch.setattr(ann, "uniprot_lookup", fake_lookup)
    runner = CliRunner()
    result = runner.invoke(cli, ["annotate", "uniprot", "P00720"])
    assert result.exit_code == 0
    assert "lysC" in result.output
    assert "Gallus gallus" in result.output


def test_annotate_uniprot_fetch_error_is_clean(project_dir: Path, monkeypatch):
    def fake_lookup(accession):
        raise ann.AnnotationFetchError("could not reach host")

    monkeypatch.setattr(ann, "uniprot_lookup", fake_lookup)
    runner = CliRunner()
    result = runner.invoke(cli, ["annotate", "uniprot", "P00720"])
    assert result.exit_code != 0
    assert "could not reach host" in result.output


def test_annotate_pfam_success(project_dir: Path, monkeypatch):
    def fake_pfam(accession):
        return [ann.PfamDomain(accession="PF00062", name="Lysozyme_C")]

    monkeypatch.setattr(ann, "pfam_domains", fake_pfam)
    runner = CliRunner()
    result = runner.invoke(cli, ["annotate", "pfam", "P00720"])
    assert result.exit_code == 0
    assert "PF00062" in result.output
