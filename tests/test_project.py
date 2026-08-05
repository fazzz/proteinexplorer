from pathlib import Path

import pytest

from proteinexplorer import project as proj
from proteinexplorer.project import ProjectError, StructureNotFoundError


def test_import_creates_project_and_record(project_dir: Path, tiny_pdb: Path):
    record = proj.import_structure(project_dir, tiny_pdb, name="tiny")

    assert record.name == "tiny"
    assert record.format == "pdb"
    assert record.n_chains == 2
    assert (project_dir / ".proteinexplorer" / "structures" / f"{record.id}.pdb").exists()
    assert (project_dir / ".proteinexplorer" / "index.json").exists()


def test_import_default_name_is_filename_stem(project_dir: Path, tiny_pdb: Path):
    record = proj.import_structure(project_dir, tiny_pdb)
    assert record.name == "tiny"


def test_list_records_after_two_imports(project_dir: Path, tiny_pdb: Path):
    proj.import_structure(project_dir, tiny_pdb, name="a")
    proj.import_structure(project_dir, tiny_pdb, name="b")
    records = proj.list_records(project_dir)
    assert {r.name for r in records} == {"a", "b"}


def test_get_record_by_id_or_name(project_dir: Path, tiny_pdb: Path):
    record = proj.import_structure(project_dir, tiny_pdb, name="tiny")
    by_id = proj.get_record(project_dir, record.id)
    by_name = proj.get_record(project_dir, "tiny")
    assert by_id.id == by_name.id == record.id


def test_get_record_missing_raises(project_dir: Path, tiny_pdb: Path):
    proj.import_structure(project_dir, tiny_pdb, name="tiny")
    with pytest.raises(StructureNotFoundError):
        proj.get_record(project_dir, "nope")


def test_get_record_prefers_most_recent_when_name_is_duplicated(project_dir: Path, tiny_pdb: Path):
    first = proj.import_structure(project_dir, tiny_pdb, name="dup")
    second = proj.import_structure(project_dir, tiny_pdb, name="dup")
    assert first.id != second.id
    resolved = proj.get_record(project_dir, "dup")
    assert resolved.id == second.id


def test_export_same_format_is_byte_copy(project_dir: Path, tiny_pdb: Path):
    record = proj.import_structure(project_dir, tiny_pdb, name="tiny")
    dest = project_dir / "out.pdb"
    proj.export_structure(project_dir, record.id, dest)
    assert dest.read_text() == tiny_pdb.read_text()


def test_export_converts_format(project_dir: Path, tiny_pdb: Path):
    record = proj.import_structure(project_dir, tiny_pdb, name="tiny")
    dest = project_dir / "out.cif"
    proj.export_structure(project_dir, record.id, dest)
    content = dest.read_text()
    assert content.startswith("data_")


def test_find_project_root_without_project_raises(project_dir: Path):
    with pytest.raises(ProjectError):
        proj.find_project_root(project_dir)


def test_log_command_appends_entries(project_dir: Path, tiny_pdb: Path):
    proj.import_structure(project_dir, tiny_pdb, name="tiny")
    root = proj.find_project_root(project_dir)
    proj.log_command(root, ["import", "tiny.pdb", "--name", "tiny"])
    log_path = root / ".proteinexplorer" / "log.json"
    import json

    entries = json.loads(log_path.read_text())
    assert len(entries) == 1
    assert entries[0]["argv"] == ["import", "tiny.pdb", "--name", "tiny"]
