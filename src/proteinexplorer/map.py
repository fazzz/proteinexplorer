"""Mapping annotations onto structures for external viewers (spec section
"Mapping").

Pure text/script generation, no external tool needed to *generate* the
script (you still need PyMOL/ChimeraX/VMD installed to actually view the
result -- this module just writes the commands for whichever one you use).

Two mechanisms, matching BioExplorer's `bio structure` conservation
mapping:

- Discrete residue groups (pockets, mutation sites, domains): a
  select-and-color script per group, one color per group.
- Continuous per-residue values (conservation scores, B-factors, any
  other float per residue): written into a copy of the structure's
  B-factor column, plus a short "spectrum/color by B-factor" script
  snippet -- the same trick BioExplorer uses for conservation mapping.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio.PDB.Structure import Structure

from proteinexplorer.models import ResidueCategory, classify_residue

DEFAULT_PALETTE = [
    "red", "blue", "green", "yellow", "orange", "purple",
    "cyan", "magenta", "salmon", "teal",
]

Tool = str  # "pymol" | "chimerax" | "vmd"


class MapError(ValueError):
    pass


@dataclass
class ResidueGroup:
    label: str
    residues: list[str]  # "chain/resnum" labels, e.g. "A/41"
    color: str


def _parse_label(label: str) -> tuple[str, int]:
    """Accepts either "A/41" (from contact.py/pocket.py's <chain>/<resname><resnum>
    style labels) or a bare "A/GLY41"-style label; extracts (chain, resnum)."""
    chain, rest = label.split("/", 1)
    digits = "".join(ch for ch in rest if ch.isdigit())
    if not digits:
        raise MapError(f"Could not parse a residue number out of label {label!r}")
    return chain, int(digits)


# --- Discrete residue-group scripts -----------------------------------

def _pymol_group_script(groups: list[ResidueGroup], object_name: str) -> str:
    lines = [f"# PyMOL coloring script ({len(groups)} group(s))", "color gray80, all"]
    for i, group in enumerate(groups):
        sel_name = f"grp_{i + 1}"
        parts = []
        for label in group.residues:
            chain, resnum = _parse_label(label)
            parts.append(f"(chain {chain} and resi {resnum})")
        selection = " or ".join(parts)
        lines.append(f"select {sel_name}, {object_name} and ({selection})")
        lines.append(f"color {group.color}, {sel_name}")
        lines.append(f"# {sel_name} = {group.label}")
    return "\n".join(lines) + "\n"


def _chimerax_group_script(groups: list[ResidueGroup], object_name: str) -> str:
    lines = [f"# ChimeraX coloring script ({len(groups)} group(s))", "color gray"]
    for group in groups:
        parts = []
        for label in group.residues:
            chain, resnum = _parse_label(label)
            parts.append(f"/{chain}:{resnum}")
        selection = "".join(parts)
        lines.append(f"color {selection} {group.color}")
        lines.append(f"# {group.label}")
    return "\n".join(lines) + "\n"


def _vmd_group_script(groups: list[ResidueGroup], object_name: str) -> str:
    lines = [f"# VMD Tcl coloring script ({len(groups)} group(s))"]
    for i, group in enumerate(groups):
        by_chain: dict[str, list[int]] = {}
        for label in group.residues:
            chain, resnum = _parse_label(label)
            by_chain.setdefault(chain, []).append(resnum)
        clauses = " or ".join(
            f"(chain {chain} and resid {' '.join(str(r) for r in resnums)})"
            for chain, resnums in by_chain.items()
        )
        lines.append(f"set sel{i + 1} [atomselect top \"{clauses}\"]")
        lines.append(f"$sel{i + 1} set colorID {i % 32}")
        lines.append(f"# sel{i + 1} = {group.label} (color {group.color})")
    return "\n".join(lines) + "\n"


_GROUP_GENERATORS = {
    "pymol": _pymol_group_script,
    "chimerax": _chimerax_group_script,
    "vmd": _vmd_group_script,
}


def generate_group_script(groups: list[ResidueGroup], tool: str = "pymol", object_name: str = "structure") -> str:
    if tool not in _GROUP_GENERATORS:
        raise MapError(f"Unknown tool: {tool!r} (expected one of {sorted(_GROUP_GENERATORS)})")
    if not groups:
        raise MapError("No residue groups to map")
    return _GROUP_GENERATORS[tool](groups, object_name)


def assign_colors(labeled_residue_lists: dict[str, list[str]], palette: list[str] | None = None) -> list[ResidueGroup]:
    palette = palette or DEFAULT_PALETTE
    groups = []
    for i, (label, residues) in enumerate(labeled_residue_lists.items()):
        groups.append(ResidueGroup(label=label, residues=residues, color=palette[i % len(palette)]))
    return groups


# --- Pocket / mutation / domain convenience wrappers ------------------

def pocket_map_script(pockets, tool: str = "pymol", object_name: str = "structure") -> str:
    """One color per pocket, using each Pocket's lining residues (from
    pocket.find_pockets)."""
    groups = assign_colors({f"pocket_{p.id}": p.residues for p in pockets if p.residues})
    return generate_group_script(groups, tool=tool, object_name=object_name)


def mutation_map_script(
    residue_labels: list[str], tool: str = "pymol", object_name: str = "structure", color: str = "red"
) -> str:
    """Highlight a set of residues (e.g. known/engineered mutation sites)
    in a single color."""
    groups = [ResidueGroup(label="mutations", residues=residue_labels, color=color)]
    return generate_group_script(groups, tool=tool, object_name=object_name)


@dataclass
class DomainRange:
    label: str
    chain_id: str
    start: int
    end: int


def domain_map_script(domains: list[DomainRange], tool: str = "pymol", object_name: str = "structure") -> str:
    """One color per named residue range (e.g. Pfam domains you already
    know the boundaries of -- this module doesn't fetch domain residue
    ranges itself, see annotate.py for the Pfam accession/name lookup)."""
    labeled = {
        d.label: [f"{d.chain_id}/{r}" for r in range(d.start, d.end + 1)]
        for d in domains
    }
    groups = assign_colors(labeled)
    return generate_group_script(groups, tool=tool, object_name=object_name)


# --- Continuous per-residue values (conservation, etc.) ---------------

def write_bfactors(structure: Structure, values_by_label: dict[str, float], default: float = 0.0) -> Structure:
    """Return a structure (the same object, modified in place) with every
    atom's B-factor set from `values_by_label` (keyed by "chain/resnum"),
    defaulting to `default` for anything not in the map. Combine with
    `spectrum_script()` to color by the written values in PyMOL/ChimeraX/VMD."""
    model = next(iter(structure))
    parsed = {_parse_label(label): value for label in values_by_label for value in [values_by_label[label]]}
    for chain in model:
        for residue in chain:
            if classify_residue(residue.resname, residue.id[0]) is not ResidueCategory.PROTEIN:
                continue
            value = parsed.get((chain.id, residue.id[1]), default)
            for atom in residue:
                atom.bfactor = value
    return structure


def spectrum_script(tool: str = "pymol", object_name: str = "structure", low_color: str = "blue", high_color: str = "red") -> str:
    if tool == "pymol":
        return f"spectrum b, {low_color}_white_{high_color}, {object_name}\n"
    if tool == "chimerax":
        return f"color byattribute bfactor {object_name} palette {low_color}:white:{high_color}\n"
    if tool == "vmd":
        return (
            "mol modcolor 0 top Beta\n"
            f"mol scaleminmax top 0 [measure minmax [atomselect top all]]\n"
            "# VMD: set the color scale to e.g. 'BWR' in the Graphics > Colors menu for blue-white-red\n"
        )
    raise MapError(f"Unknown tool: {tool!r} (expected pymol, chimerax, or vmd)")
