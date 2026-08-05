"""Project management under a .proteinexplorer/ directory.

Mirrors the ChemExplorer / BioExplorer project layout:

    .proteinexplorer/
        index.json      # metadata for every imported structure
        log.json        # recorded CLI invocations, for future `prot replay`
        structures/
            <id>.pdb | <id>.cif   # untouched copy of the imported file
"""

from __future__ import annotations

import json
import secrets
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from proteinexplorer import io as pio

PROJECT_DIRNAME = ".proteinexplorer"
INDEX_FILENAME = "index.json"
LOG_FILENAME = "log.json"
STRUCTURES_DIRNAME = "structures"


class ProjectError(Exception):
    pass


class StructureNotFoundError(ProjectError):
    pass


@dataclass
class StructureRecord:
    id: str
    name: str
    source_path: str
    format: str
    imported_at: str
    n_models: int
    n_chains: int
    n_residues: int
    n_atoms: int
    category_totals: dict[str, int]
    hetero_resnames: list[str]
    has_altloc: bool
    tags: list[str] = field(default_factory=list)


def _project_dir(root: Path) -> Path:
    return root / PROJECT_DIRNAME


def find_project_root(start: str | Path = ".") -> Path:
    """Walk upward from `start` looking for a .proteinexplorer/ directory."""
    current = Path(start).resolve()
    for candidate in [current, *current.parents]:
        if (candidate / PROJECT_DIRNAME).is_dir():
            return candidate
    raise ProjectError(
        "No .proteinexplorer project found in this directory or any parent. "
        "Run `prot import <file>` first (it initializes a project automatically)."
    )


def init_project(root: str | Path = ".") -> Path:
    """Create a new .proteinexplorer/ project if one doesn't already exist
    at `root`. Returns the project root path."""
    root = Path(root).resolve()
    pdir = _project_dir(root)
    (pdir / STRUCTURES_DIRNAME).mkdir(parents=True, exist_ok=True)
    index_path = pdir / INDEX_FILENAME
    if not index_path.exists():
        index_path.write_text(json.dumps([], indent=2))
    log_path = pdir / LOG_FILENAME
    if not log_path.exists():
        log_path.write_text(json.dumps([], indent=2))
    return root


def _load_index(root: Path) -> list[dict]:
    index_path = _project_dir(root) / INDEX_FILENAME
    if not index_path.exists():
        return []
    return json.loads(index_path.read_text())


def _save_index(root: Path, records: list[dict]) -> None:
    index_path = _project_dir(root) / INDEX_FILENAME
    index_path.write_text(json.dumps(records, indent=2))


def _new_structure_id() -> str:
    return "p_" + secrets.token_hex(4)


def log_command(root: Path, argv: list[str] | None = None) -> None:
    """Append the current CLI invocation to log.json for future `prot replay`."""
    argv = argv if argv is not None else sys.argv[1:]
    log_path = _project_dir(root) / LOG_FILENAME
    entries = json.loads(log_path.read_text()) if log_path.exists() else []
    entries.append(
        {
            "argv": argv,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
    )
    log_path.write_text(json.dumps(entries, indent=2))


def import_structure(
    root: str | Path,
    src_path: str | Path,
    name: str | None = None,
    fmt: str | None = None,
) -> StructureRecord:
    """Import a structure file into the project: copy the raw file,
    parse it to build a summary, and record it in index.json."""
    root = init_project(root)
    src = Path(src_path)
    if not src.exists():
        raise ProjectError(f"File not found: {src}")

    fmt = fmt or pio.infer_format(src)
    ext = ".pdb" if fmt == "pdb" else ".cif"

    struct_id = _new_structure_id()
    stored_path = _project_dir(root) / STRUCTURES_DIRNAME / f"{struct_id}{ext}"
    pio.copy_raw(src, stored_path)

    structure = pio.load_structure(stored_path, structure_id=struct_id, fmt=fmt)
    summary = pio.summarize(structure)

    record = StructureRecord(
        id=struct_id,
        name=name or src.stem,
        source_path=str(src),
        format=fmt,
        imported_at=datetime.now(timezone.utc).isoformat(),
        n_models=summary.n_models,
        n_chains=summary.n_chains,
        n_residues=summary.n_residues,
        n_atoms=summary.n_atoms,
        category_totals=summary.category_totals(),
        hetero_resnames=summary.hetero_resnames,
        has_altloc=summary.has_altloc,
    )

    records = _load_index(root)
    records.append(asdict(record))
    _save_index(root, records)
    return record


def get_record(root: str | Path, struct_id: str) -> StructureRecord:
    root = find_project_root(root)
    entries = _load_index(root)
    # ids are always unique; check those first.
    for entry in entries:
        if entry["id"] == struct_id:
            return StructureRecord(**dict(entry))
    # names are not enforced unique (re-importing the same file, or two
    # different files under the same --name, is allowed) -- prefer the
    # most recently imported match so `prot info <name>` etc. follow the
    # latest import rather than silently pinning to the first one.
    for entry in reversed(entries):
        if entry["name"] == struct_id:
            return StructureRecord(**dict(entry))
    raise StructureNotFoundError(f"No structure with id or name '{struct_id}' in project")


def structure_path(root: str | Path, struct_id: str) -> Path:
    root = find_project_root(root)
    record = get_record(root, struct_id)
    ext = ".pdb" if record.format == "pdb" else ".cif"
    return _project_dir(root) / STRUCTURES_DIRNAME / f"{record.id}{ext}"


def export_structure(
    root: str | Path,
    struct_id: str,
    dest_path: str | Path,
    fmt: str | None = None,
) -> Path:
    root = find_project_root(root)
    src = structure_path(root, struct_id)
    dest = Path(dest_path)
    dest_fmt = fmt or pio.infer_format(dest)

    record = get_record(root, struct_id)
    if dest_fmt == record.format:
        pio.copy_raw(src, dest)
    else:
        structure = pio.load_structure(src, structure_id=record.id, fmt=record.format)
        pio.save_structure(structure, dest, fmt=dest_fmt)
    return dest


def list_records(root: str | Path) -> list[StructureRecord]:
    root = find_project_root(root)
    return [StructureRecord(**entry) for entry in _load_index(root)]
