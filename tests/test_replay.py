import json
from pathlib import Path

import pytest

from proteinexplorer import project as proj
from proteinexplorer import replay


def _run(project_dir: Path, argv: list[str]):
    from click.testing import CliRunner
    from proteinexplorer.cli import cli

    runner = CliRunner()
    result = runner.invoke(cli, argv)
    assert result.exit_code == 0, result.output
    return result


@pytest.fixture()
def workflow_project(project_dir: Path, tiny_pdb: Path):
    """A project with import -> mutate -> mutate(using the intermediate's
    literal old id) already recorded in log.json."""
    _run(project_dir, ["import", str(tiny_pdb), "--name", "tiny"])
    result = _run(project_dir, ["mutate", "tiny", "--chain", "A", "--resid", "1", "--to", "VAL"])
    mutant_id = result.output.splitlines()[-1].split("(")[-1].rstrip(")")
    _run(project_dir, ["mutate", mutant_id, "--chain", "A", "--resid", "2", "--to", "TRP"])
    return project_dir, mutant_id


# --- plan / dry run ------------------------------------------------------

def test_plan_returns_full_log_by_default(workflow_project):
    project_dir, _ = workflow_project
    entries = replay.plan(project_dir)
    assert len(entries) == 3
    assert entries[0]["argv"][0] == "import"


def test_plan_respects_start_end(workflow_project):
    project_dir, _ = workflow_project
    entries = replay.plan(project_dir, start=2, end=2)
    assert len(entries) == 1
    assert entries[0]["argv"][0] == "mutate"


def test_replay_dry_run_does_not_touch_project(workflow_project):
    project_dir, _ = workflow_project
    before = json.loads((project_dir / ".proteinexplorer" / "index.json").read_text())
    result = replay.replay(project_dir, dry_run=True)
    after = json.loads((project_dir / ".proteinexplorer" / "index.json").read_text())
    assert before == after
    assert result.backup_dir is None
    assert len(result.steps) == 3
    assert all(s.output == "(dry run)" for s in result.steps)


# --- replay with reset (the default) --------------------------------

def test_replay_reconstructs_project(workflow_project):
    project_dir, _ = workflow_project
    result = replay.replay(project_dir, reset=True)
    assert result.backup_dir is not None
    assert result.backup_dir.exists()
    assert result.n_failed == 0

    records = proj.list_records(project_dir)
    names = {r.name for r in records}
    assert names == {"tiny", "tiny_A1ALAVAL", "tiny_A1ALAVAL_A2GLYTRP"}


def test_replay_rewrites_old_literal_id(workflow_project):
    project_dir, old_mutant_id = workflow_project
    result = replay.replay(project_dir, reset=True)
    step3 = result.steps[2]
    assert step3.argv[1] == old_mutant_id  # original argv unchanged
    assert step3.rewritten_argv[1] != old_mutant_id  # rewritten to the new id
    assert step3.rewritten_argv[1].startswith("p_")


def test_replay_backup_contains_pre_replay_state(workflow_project):
    project_dir, _ = workflow_project
    result = replay.replay(project_dir, reset=True)
    backup_index = json.loads((result.backup_dir / "index.json").read_text())
    assert len(backup_index) == 3  # tiny + 2 mutants, before reset


def test_replay_skips_default_skip_list(project_dir: Path, tiny_pdb: Path):
    _run(project_dir, ["import", str(tiny_pdb), "--name", "tiny"])
    root = proj.find_project_root(project_dir)
    # simulate a logged `view` invocation (normally skipped by default)
    proj.log_command(root, ["view", "tiny"])

    result = replay.replay(project_dir, reset=True)
    assert result.steps[-1].skipped is True


def test_replay_custom_skip_set(workflow_project):
    project_dir, _ = workflow_project
    result = replay.replay(project_dir, reset=True, skip={"mutate"})
    assert result.steps[1].skipped is True
    assert result.steps[2].skipped is True
    records = proj.list_records(project_dir)
    assert {r.name for r in records} == {"tiny"}


def test_replay_empty_range_raises(workflow_project):
    project_dir, _ = workflow_project
    with pytest.raises(replay.ReplayError):
        replay.replay(project_dir, start=10, end=10)


def test_replay_stops_on_first_failure_by_default(project_dir: Path, tiny_pdb: Path):
    _run(project_dir, ["import", str(tiny_pdb), "--name", "tiny"])
    root = proj.find_project_root(project_dir)
    # log a command that will fail on replay (nonexistent residue)
    proj.log_command(root, ["mutate", "tiny", "--chain", "A", "--resid", "999", "--to", "VAL"])

    with pytest.raises(replay.ReplayError):
        replay.replay(project_dir, reset=True)


def test_replay_continue_on_error(project_dir: Path, tiny_pdb: Path):
    _run(project_dir, ["import", str(tiny_pdb), "--name", "tiny"])
    root = proj.find_project_root(project_dir)
    proj.log_command(root, ["mutate", "tiny", "--chain", "A", "--resid", "999", "--to", "VAL"])
    proj.log_command(root, ["mutate", "tiny", "--chain", "A", "--resid", "1", "--to", "VAL"])

    result = replay.replay(project_dir, reset=True, continue_on_error=True)
    assert result.n_failed == 1
    assert result.steps[-1].exit_code == 0


def test_replay_no_reset_keeps_existing_state(workflow_project):
    project_dir, _ = workflow_project
    before = len(proj.list_records(project_dir))
    result = replay.replay(project_dir, reset=False, start=1, end=1)
    assert result.backup_dir is None
    after = len(proj.list_records(project_dir))
    assert after == before + 1  # the import ran again, adding one more record
