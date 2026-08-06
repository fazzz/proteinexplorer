"""Static plots (spec section "Plot").

Uses matplotlib, same as ChemExplorer's/BioExplorer's viz.py (the
`[viz]` optional extra: `pip install -e ".[viz]"`). Every function raises
a clear PlotExtraNotAvailableError if matplotlib isn't installed, and
either saves to a given output path or returns the Figure for further
use (e.g. from a notebook).
"""

from __future__ import annotations

from pathlib import Path

from Bio.PDB.Structure import Structure


class PlotExtraNotAvailableError(RuntimeError):
    pass


def _require_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")  # headless-safe backend; caller can re-select before importing pyplot elsewhere
        import matplotlib.pyplot as plt
        return plt
    except ImportError as exc:
        raise PlotExtraNotAvailableError(
            "Plotting needs matplotlib. Install it with `pip install -e \".[viz]\"` "
            "(or `pip install matplotlib`)."
        ) from exc


# --- Ramachandran plot ---------------------------------------------------

def ramachandran_plot(structure: Structure, output_path: str | Path | None = None, chain_id: str | None = None):
    """Phi/psi scatter plot, with the core alpha-helix and beta-strand
    Ramachandran regions used by secondary.py's geometric classifier
    shaded for reference."""
    plt = _require_matplotlib()
    from proteinexplorer import geometry as geom
    from proteinexplorer.models import ResidueCategory, classify_residue
    from proteinexplorer.secondary import _ALPHA_PHI, _ALPHA_PSI, _BETA_PHI

    model = next(iter(structure))
    points = []
    for chain in model:
        if chain_id is not None and chain.id != chain_id:
            continue
        residues = [r for r in chain if classify_residue(r.resname, r.id[0]) is ResidueCategory.PROTEIN]
        for residue in residues:
            try:
                torsions = geom.backbone_torsions(chain, residue.id[1])
            except geom.GeometryError:
                continue
            if torsions.phi is not None and torsions.psi is not None:
                points.append((torsions.phi, torsions.psi))

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.add_patch(
        plt.Rectangle(
            (_ALPHA_PHI[0], _ALPHA_PSI[0]), _ALPHA_PHI[1] - _ALPHA_PHI[0], _ALPHA_PSI[1] - _ALPHA_PSI[0],
            facecolor="tab:red", alpha=0.15, label="alpha region (geometric classifier)",
        )
    )
    ax.add_patch(
        plt.Rectangle(
            (_BETA_PHI[0], 90), _BETA_PHI[1] - _BETA_PHI[0], 90,
            facecolor="tab:blue", alpha=0.15, label="beta region (geometric classifier)",
        )
    )
    ax.add_patch(
        plt.Rectangle(
            (_BETA_PHI[0], -180), _BETA_PHI[1] - _BETA_PHI[0], 30,
            facecolor="tab:blue", alpha=0.15,
        )
    )
    if points:
        phis, psis = zip(*points)
        ax.scatter(phis, psis, s=12, c="black", alpha=0.7)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("phi (deg)")
    ax.set_ylabel("psi (deg)")
    ax.set_title(f"Ramachandran plot ({len(points)} residues)")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.legend(loc="upper right", fontsize=7)
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path)
        plt.close(fig)
        return Path(output_path)
    return fig


# --- Contact map heatmap ---------------------------------------------------

def contact_map_plot(
    structure: Structure,
    output_path: str | Path | None = None,
    atoms=None,
    mode: str = "ca",
    cutoff: float = 8.0,
):
    plt = _require_matplotlib()
    from proteinexplorer import contact as ct

    cm = ct.contact_map(structure, atoms=atoms, mode=mode, cutoff=cutoff)
    n = len(cm.labels)

    fig, ax = plt.subplots(figsize=(max(4, n * 0.15), max(4, n * 0.15)))
    ax.imshow(cm.matrix, cmap="Greys", origin="lower", interpolation="none")
    ax.set_title(f"Contact map ({mode}, cutoff={cutoff} A, {n} residues)")
    step = max(1, n // 20)
    ticks = list(range(0, n, step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([cm.labels[i] for i in ticks], rotation=90, fontsize=6)
    ax.set_yticks(ticks)
    ax.set_yticklabels([cm.labels[i] for i in ticks], fontsize=6)
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path)
        plt.close(fig)
        return Path(output_path)
    return fig


# --- Secondary structure diagram --------------------------------------

_SS_COLORS = {"H": "tab:red", "E": "tab:blue", "C": "lightgray", "-": "lightgray"}


def secondary_structure_plot(
    structure: Structure,
    output_path: str | Path | None = None,
    method: str = "auto",
):
    """Linear secondary structure cartoon: one horizontal track per chain,
    colored by SS class along the residue sequence (helix=red,
    strand=blue, coil=gray)."""
    plt = _require_matplotlib()
    from proteinexplorer import secondary as sec

    residues, used_method = sec.secondary_structure(structure, method=method)
    by_chain: dict[str, list] = {}
    for r in residues:
        by_chain.setdefault(r.chain_id, []).append(r)
    for chain_id in by_chain:
        by_chain[chain_id].sort(key=lambda r: r.resseq)

    chain_ids = sorted(by_chain)
    fig, ax = plt.subplots(figsize=(10, 0.6 * len(chain_ids) + 1))
    for row, chain_id in enumerate(chain_ids):
        chain_residues = by_chain[chain_id]
        for r in chain_residues:
            code = r.code if r.code in ("H", "E") else "C"
            ax.barh(row, 1, left=r.resseq, height=0.8, color=_SS_COLORS[code], edgecolor="none")
        # fixed offset just left of the axes frame, independent of each
        # chain's own residue numbering range (avoids labels landing far
        # off-plot when chains have very different resid ranges)
        ax.text(
            -0.01, row, chain_id, ha="right", va="center", fontsize=9,
            transform=ax.get_yaxis_transform(),
        )

    ax.set_yticks([])
    ax.set_xlabel("residue number")
    ax.set_title(f"Secondary structure (method={used_method})")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=_SS_COLORS["H"], label="helix"),
        plt.Rectangle((0, 0), 1, 1, color=_SS_COLORS["E"], label="strand"),
        plt.Rectangle((0, 0), 1, 1, color=_SS_COLORS["C"], label="coil"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, ncol=3)
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path)
        plt.close(fig)
        return Path(output_path)
    return fig
