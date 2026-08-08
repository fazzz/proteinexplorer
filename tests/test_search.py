from pathlib import Path

import pytest

from proteinexplorer import search


# --- foldseek availability -------------------------------------------

def test_foldseek_binary_not_found():
    assert search.foldseek_binary() is None


def test_easy_search_raises_clean_error_without_binary(tmp_path: Path):
    with pytest.raises(search.FoldseekNotAvailableError):
        search.easy_search(tmp_path / "query.pdb", tmp_path / "targets")


def test_easy_search_error_mentions_install_pointer(tmp_path: Path):
    with pytest.raises(search.FoldseekNotAvailableError, match="foldseek"):
        search.easy_search(tmp_path / "query.pdb", tmp_path / "targets")


def test_createdb_raises_clean_error_without_binary(tmp_path: Path):
    with pytest.raises(search.FoldseekNotAvailableError):
        search.createdb(tmp_path / "structures", tmp_path / "db" / "targetDB")


# --- tabular parsing (pure, no binary needed) ---------------------------

def test_parse_tabular_basic():
    columns = ["query", "target", "evalue", "bits"]
    text = "q1\tt1\t1e-10\t120.5\nq1\tt2\t1e-5\t80.0\n"
    hits = search._parse_tabular(text, columns)
    assert len(hits) == 2
    assert hits[0].query == "q1"
    assert hits[0].target == "t1"
    assert hits[0].evalue == pytest.approx(1e-10)
    assert hits[0].bits == pytest.approx(120.5)
    assert hits[1].target == "t2"


def test_parse_tabular_skips_blank_lines():
    columns = ["query", "target"]
    text = "q1\tt1\n\n\nq1\tt2\n"
    hits = search._parse_tabular(text, columns)
    assert len(hits) == 2


def test_parse_tabular_empty_text():
    assert search._parse_tabular("", ["query", "target"]) == []


def test_search_hit_alntmscore_property():
    hit = search.SearchHit(fields={"alntmscore": "0.87"})
    assert hit.alntmscore == pytest.approx(0.87)


def test_search_hit_missing_field_is_none():
    hit = search.SearchHit(fields={"query": "q1"})
    assert hit.evalue is None
    assert hit.bits is None
    assert hit.alntmscore is None


def test_search_hit_fields_preserved_verbatim():
    hit = search.SearchHit(fields={"query": "q1", "custom_col": "xyz"})
    assert hit.fields["custom_col"] == "xyz"
