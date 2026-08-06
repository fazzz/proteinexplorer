from pathlib import Path

import pytest

from proteinexplorer import io as pio
from proteinexplorer import plot as plt_mod


@pytest.fixture()
def tiny_structure(tiny_pdb: Path):
    return pio.load_structure(tiny_pdb, "t")


def test_ramachandran_plot_saves_file(tiny_structure, tmp_path: Path):
    out = tmp_path / "rama.png"
    result = plt_mod.ramachandran_plot(tiny_structure, out)
    assert result == out
    assert out.exists()
    assert out.stat().st_size > 0


def test_ramachandran_plot_returns_figure_without_output_path(tiny_structure):
    fig = plt_mod.ramachandran_plot(tiny_structure)
    assert fig is not None
    assert len(fig.axes) == 1


def test_ramachandran_plot_chain_filter(tiny_structure, tmp_path: Path):
    out = tmp_path / "rama_a.png"
    plt_mod.ramachandran_plot(tiny_structure, out, chain_id="A")
    assert out.exists()


def test_contact_map_plot_saves_file(tiny_structure, tmp_path: Path):
    out = tmp_path / "contact.png"
    result = plt_mod.contact_map_plot(tiny_structure, out)
    assert result == out
    assert out.exists()
    assert out.stat().st_size > 0


def test_contact_map_plot_heavy_mode(tiny_structure, tmp_path: Path):
    out = tmp_path / "contact_heavy.png"
    plt_mod.contact_map_plot(tiny_structure, out, mode="heavy", cutoff=5.0)
    assert out.exists()


def test_secondary_structure_plot_saves_file(tiny_structure, tmp_path: Path):
    out = tmp_path / "ss.png"
    result = plt_mod.secondary_structure_plot(tiny_structure, out, method="geometric")
    assert result == out
    assert out.exists()
    assert out.stat().st_size > 0


def test_secondary_structure_plot_returns_figure(tiny_structure):
    fig = plt_mod.secondary_structure_plot(tiny_structure, method="geometric")
    assert fig is not None


def test_plot_extra_not_available_error_message(tiny_structure, monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "matplotlib":
            raise ImportError("mocked: matplotlib not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(plt_mod.PlotExtraNotAvailableError):
        plt_mod.ramachandran_plot(tiny_structure)
