import pytest

from proteinexplorer import view


# --- viewer_binary -----------------------------------------------------

def test_viewer_binary_not_found_for_any_tool():
    assert view.viewer_binary("pymol") is None
    assert view.viewer_binary("chimerax") is None
    assert view.viewer_binary("vmd") is None


def test_viewer_binary_unknown_tool_raises():
    with pytest.raises(view.ViewError):
        view.viewer_binary("bogus")


# --- build_command -------------------------------------------------------

def test_build_command_raises_without_binary():
    with pytest.raises(view.ViewerNotAvailableError):
        view.build_command("pymol", "structure.pdb")


def test_build_command_error_mentions_tool_and_install_hint():
    with pytest.raises(view.ViewerNotAvailableError, match="pymol"):
        view.build_command("pymol", "structure.pdb")


def test_build_command_pymol_with_binary(monkeypatch):
    monkeypatch.setattr(view, "viewer_binary", lambda tool: "/usr/bin/pymol")
    cmd = view.build_command("pymol", "structure.pdb")
    assert cmd == ["/usr/bin/pymol", "structure.pdb"]


def test_build_command_pymol_with_script(monkeypatch):
    monkeypatch.setattr(view, "viewer_binary", lambda tool: "/usr/bin/pymol")
    cmd = view.build_command("pymol", "structure.pdb", script_path="script.pml")
    assert cmd == ["/usr/bin/pymol", "structure.pdb", "script.pml"]


def test_build_command_vmd_uses_dash_e_for_script(monkeypatch):
    monkeypatch.setattr(view, "viewer_binary", lambda tool: "/usr/bin/vmd")
    cmd = view.build_command("vmd", "structure.pdb", script_path="script.vmd")
    assert cmd == ["/usr/bin/vmd", "structure.pdb", "-e", "script.vmd"]


def test_build_command_chimerax_with_binary(monkeypatch):
    monkeypatch.setattr(view, "viewer_binary", lambda tool: "/usr/bin/chimerax")
    cmd = view.build_command("chimerax", "structure.pdb")
    assert cmd == ["/usr/bin/chimerax", "structure.pdb"]


# --- launch ----------------------------------------------------------

def test_launch_raises_without_binary():
    with pytest.raises(view.ViewerNotAvailableError):
        view.launch("pymol", "structure.pdb")


def test_launch_calls_popen_with_built_command(monkeypatch):
    monkeypatch.setattr(view, "viewer_binary", lambda tool: "/usr/bin/pymol")

    captured = {}

    class FakePopen:
        def __init__(self, cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs

    monkeypatch.setattr(view.subprocess, "Popen", FakePopen)
    view.launch("pymol", "structure.pdb", script_path="s.pml")
    assert captured["cmd"] == ["/usr/bin/pymol", "structure.pdb", "s.pml"]
    assert captured["kwargs"]["start_new_session"] is True
