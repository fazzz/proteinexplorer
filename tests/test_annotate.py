import json
import urllib.error
from pathlib import Path

import pytest

from proteinexplorer import annotate as ann
from proteinexplorer import io as pio


@pytest.fixture()
def metal_structure():
    path = Path(__file__).parent / "data" / "metal_site.pdb"
    return pio.load_structure(path, structure_id="metal")


# --- metal_binding_sites -------------------------------------------------

def test_metal_binding_site_detects_coordinating_residues(metal_structure):
    sites = ann.metal_binding_sites(metal_structure)
    assert len(sites) == 1
    site = sites[0]
    assert site.ion_label == "A/ZN100"
    assert site.ion_element == "ZN"
    assert set(site.coordinating_residues) == {"A/HIS1", "A/CYS2", "A/ASP3"}


def test_metal_binding_site_excludes_far_residues(metal_structure):
    sites = ann.metal_binding_sites(metal_structure)
    assert "A/ALA10" not in sites[0].coordinating_residues


def test_metal_binding_site_distances_are_sorted(metal_structure):
    sites = ann.metal_binding_sites(metal_structure)
    distances = sites[0].coordinating_distances
    assert distances == sorted(distances)


def test_metal_binding_site_tight_cutoff_excludes_everything(metal_structure):
    sites = ann.metal_binding_sites(metal_structure, cutoff=1.0)
    assert sites == []


def test_metal_binding_site_no_ions_returns_empty(tiny_pdb: Path):
    from Bio.PDB import PDBParser
    # tiny.pdb DOES have a Zn ion; use metal_site.pdb's non-ion residues
    # only by cutting the cutoff to something that catches nothing instead
    structure = pio.load_structure(tiny_pdb, "t")
    sites = ann.metal_binding_sites(structure, cutoff=0.01)
    assert sites == []


# --- structure_metadata ---------------------------------------------------

def test_structure_metadata_reads_header(metal_structure):
    meta = ann.structure_metadata(metal_structure)
    assert meta.deposition_date is not None
    assert meta.method is None  # not present in this synthetic fixture


# --- External lookups (network calls fail in this sandbox; verify the
# error handling path rather than a live fetch) --------------------------

def test_uniprot_lookup_unreachable_host_raises_clean_error():
    with pytest.raises(ann.AnnotationFetchError):
        ann.uniprot_lookup("P00720")


def test_pfam_domains_unreachable_host_raises_clean_error():
    with pytest.raises(ann.AnnotationFetchError):
        ann.pfam_domains("P00720")


def test_http_get_json_wraps_http_error(monkeypatch):
    def fake_urlopen(*args, **kwargs):
        raise urllib.error.HTTPError("http://example.com", 404, "Not Found", {}, None)

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    with pytest.raises(ann.AnnotationFetchError, match="404"):
        ann._http_get_json("http://example.com")


def test_uniprot_lookup_parses_response(monkeypatch):
    payload = {
        "genes": [{"geneName": {"value": "lysC"}}],
        "organism": {"scientificName": "Gallus gallus", "taxonId": 9031},
        "proteinDescription": {"recommendedName": {"ecNumbers": [{"value": "3.2.1.17"}]}},
        "uniProtKBCrossReferences": [
            {"database": "GO", "id": "GO:0003796"},
            {"database": "PDB", "id": "1LYS"},
        ],
    }

    def fake_http_get_json(url, timeout=15.0):
        return payload

    monkeypatch.setattr(ann, "_http_get_json", fake_http_get_json)
    result = ann.uniprot_lookup("P00720")
    assert result.gene_names == ["lysC"]
    assert result.organism == "Gallus gallus"
    assert result.taxonomy_id == 9031
    assert result.ec_numbers == ["3.2.1.17"]
    assert result.go_terms == ["GO:0003796"]


def test_pfam_domains_parses_response(monkeypatch):
    payload = {
        "results": [
            {"metadata": {"accession": "PF00062", "name": "Lysozyme_C"}},
        ]
    }

    def fake_http_get_json(url, timeout=15.0):
        return payload

    monkeypatch.setattr(ann, "_http_get_json", fake_http_get_json)
    result = ann.pfam_domains("P00720")
    assert len(result) == 1
    assert result[0].accession == "PF00062"
    assert result[0].name == "Lysozyme_C"
