"""Workflow replay (spec section "Replay").

Re-executes the CLI invocations recorded in `.proteinexplorer/log.json`
(every command that mutates the project has been logging there since
`project.py` was first written) via Click's CliRunner, in-process --
the same approach BioExplorer's `bio replay` uses.

Structure IDs are randomly generated on every `import` (see
project._new_structure_id), so a replayed command that references a
structure by its old literal ID (rather than by --name) would otherwise
break once the project is reset and re-populated with fresh IDs. This is
handled here: before resetting, the current index.json is read to build
an old-id -> name map; as each replayed command creates a new structure
record, its name is matched back to that map to learn the corresponding
old-id -> new-id substitution, which is then applied to every later
command's argv before it runs. This isn't perfect (ambiguous when two
records share a name -- resolved oldest-first) but covers the common
case of copy-pasting an id from one command's output into the next.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from proteinexplorer import project as proj

_ID_PREFIX = "p_"
DEFAULT_SKIP = {"view", "predict", "annotate", "replay"}


class ReplayError(Exception):
    pass


@dataclass
class ReplayStepResult:
    index: int
    argv: list[str]
    rewritten_argv: list[str]
    skipped: bool
    exit_code: int | None
    output: str


@dataclass
class ReplayResult:
    steps: list[ReplayStepResult] = field(default_factory=list)
    backup_dir: Path | None = None

    @property
    def n_failed(self) -> int:
        return sum(1 for s in self.steps if not s.skipped and s.exit_code not in (0, None))


def _load_log(root: Path) -> list[dict]:
    log_path = proj._project_dir(root) / proj.LOG_FILENAME
    if not log_path.exists():
        return []
    return json.loads(log_path.read_text())


def _load_index(root: Path) -> list[dict]:
    index_path = proj._project_dir(root) / proj.INDEX_FILENAME
    if not index_path.exists():
        return []
    return json.loads(index_path.read_text())


def backup_project(root: Path) -> Path:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    backup_dir = root / f".proteinexplorer_prereplay_{timestamp}"
    shutil.copytree(proj._project_dir(root), backup_dir)
    return backup_dir


def _reset_project(root: Path) -> None:
    pdir = proj._project_dir(root)
    structures_dir = pdir / proj.STRUCTURES_DIRNAME
    if structures_dir.exists():
        shutil.rmtree(structures_dir)
    structures_dir.mkdir(parents=True, exist_ok=True)
    (pdir / proj.INDEX_FILENAME).write_text(json.dumps([], indent=2))
    # log.json is intentionally left untouched -- we're replaying from it.


def _rewrite_argv(argv: list[str], id_remap: dict[str, str]) -> list[str]:
    return [id_remap.get(token, token) for token in argv]


def _command_key(argv: list[str]) -> str:
    """The top-level command/group name an argv starts with, used for
    --skip matching (e.g. "predict" matches both "predict colabfold" and
    "predict alphafold")."""
    return argv[0] if argv else ""


def plan(root: str | Path, start: int = 1, end: int | None = None) -> list[dict]:
    """Return the slice of logged entries that would be replayed, without
    running anything (for --dry-run)."""
    root = proj.find_project_root(root)
    entries = _load_log(root)
    end = end if end is not None else len(entries)
    return entries[start - 1:end]


def replay(
    root: str | Path,
    start: int = 1,
    end: int | None = None,
    skip: set[str] | None = None,
    continue_on_error: bool = False,
    reset: bool = True,
    dry_run: bool = False,
) -> ReplayResult:
    root = proj.find_project_root(root)
    skip = DEFAULT_SKIP if skip is None else (skip | {"replay"})
    entries = _load_log(root)
    end = end if end is not None else len(entries)
    selected = entries[start - 1:end]
    if not selected:
        raise ReplayError("No logged commands in the selected range")

    old_index = _load_index(root)
    id_to_name = {e["id"]: e["name"] for e in old_index}
    id_remap: dict[str, str] = {}

    result = ReplayResult()

    if dry_run:
        for i, entry in enumerate(selected, start=start):
            argv = entry["argv"]
            skipped = _command_key(argv) in skip
            result.steps.append(
                ReplayStepResult(index=i, argv=argv, rewritten_argv=_rewrite_argv(argv, id_remap),
                                  skipped=skipped, exit_code=None, output="(dry run)")
            )
        return result

    if reset:
        result.backup_dir = backup_project(root)
        _reset_project(root)

    from click.testing import CliRunner
    from proteinexplorer.cli import cli as cli_group

    runner = CliRunner()
    original_cwd = Path.cwd()
    import os
    os.chdir(root)
    try:
        for i, entry in enumerate(selected, start=start):
            argv = entry["argv"]
            if _command_key(argv) in skip:
                result.steps.append(
                    ReplayStepResult(index=i, argv=argv, rewritten_argv=argv, skipped=True, exit_code=None, output="")
                )
                continue

            rewritten = _rewrite_argv(argv, id_remap)
            index_before = _load_index(root)
            invocation = runner.invoke(cli_group, rewritten)
            step = ReplayStepResult(
                index=i, argv=argv, rewritten_argv=rewritten, skipped=False,
                exit_code=invocation.exit_code, output=invocation.output,
            )
            result.steps.append(step)

            if invocation.exit_code != 0 and not continue_on_error:
                raise ReplayError(
                    f"Step {i} failed (exit {invocation.exit_code}): {' '.join(rewritten)}\n{invocation.output}"
                )

            index_after = _load_index(root)
            if len(index_after) > len(index_before):
                new_entry = index_after[-1]
                old_id = next(
                    (oid for oid, name in id_to_name.items() if name == new_entry["name"] and oid not in id_remap),
                    None,
                )
                if old_id is not None:
                    id_remap[old_id] = new_entry["id"]
    finally:
        os.chdir(original_cwd)

    return result
