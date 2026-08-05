"""Contact / interaction analysis (spec section "Contact Analysis").

Distance-heavy-atom-only criteria are used throughout (no hydrogens
required), which matches the typical crystallographic PDB/mmCIF input this
tool targets. Cutoffs are the commonly used defaults from the structural
biology literature (Arpeggio/PLIP-style) and are all overridable.

Built on top of geometry.py (distance, fit_plane) and models.py
(classification) so the same primitives are reused rather than
reimplemented.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Structure import Structure

from proteinexplorer import geometry as geom
from proteinexplorer.descriptor import HYDROPHOBIC_RESIDUES
from proteinexplorer.models import ResidueCategory, classify_residue

# --- Interaction group definitions ------------------------------------

HBOND_ELEMENTS = frozenset({"N", "O"})

BASIC_GROUP_ATOMS: dict[str, tuple[str, ...]] = {
    "ARG": ("NH1", "NH2", "NE"),
    "LYS": ("NZ",),
    "HIS": ("ND1", "NE2"),
}
ACIDIC_GROUP_ATOMS: dict[str, tuple[str, ...]] = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
}
CATION_GROUP_ATOMS: dict[str, tuple[str, ...]] = {
    "ARG": ("NH1", "NH2", "NE", "CZ"),
    "LYS": ("NZ",),
}
# 6-membered rings for pi-pi/cation-pi geometry (Trp's indole is
# approximated by its benzo ring, His by its imidazole ring).
AROMATIC_RING_ATOMS: dict[str, tuple[str, ...]] = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
}

_DISULFIDE_SG_DISTANCE = 2.5


def _residue_label(residue) -> str:
    chain_id = residue.get_parent().id
    return f"{chain_id}/{residue.resname}{residue.id[1]}"


def _protein_residues(structure: Structure):
    model = next(iter(structure))
    for chain in model:
        for residue in chain:
            if classify_residue(residue.resname, residue.id[0]) is ResidueCategory.PROTEIN:
                yield residue


def _group_atoms(structure: Structure, resname_to_atoms: dict[str, tuple[str, ...]]):
    """Yield (residue, [Atom,...]) for every residue whose resname is a key
    in resname_to_atoms, restricted to the atom names listed."""
    for residue in _protein_residues(structure):
        names = resname_to_atoms.get(residue.resname.strip().upper())
        if not names:
            continue
        atoms = [residue[name] for name in names if name in residue]
        if atoms:
            yield residue, atoms


# --- Hydrogen bonds ----------------------------------------------------

@dataclass
class HydrogenBond:
    donor_residue: str
    donor_atom: str
    acceptor_residue: str
    acceptor_atom: str
    distance: float


def find_hydrogen_bonds(structure: Structure, cutoff: float = 3.5) -> list[HydrogenBond]:
    """Heavy-atom-only H-bond detection: any N/O...N/O pair within `cutoff`
    Angstrom, on different residues. Since no hydrogens are required, the
    "donor"/"acceptor" labeling is nominal (either atom could be either
    role); both are reported so the caller can judge chemical plausibility.
    """
    model = next(iter(structure))
    candidates = [
        atom for atom in model.get_atoms()
        if (atom.element or "").strip().upper() in HBOND_ELEMENTS
    ]
    if not candidates:
        return []
    ns = NeighborSearch(candidates)
    pairs = ns.search_all(cutoff, level="A")

    results: list[HydrogenBond] = []
    seen = set()
    for atom_a, atom_b in pairs:
        res_a, res_b = atom_a.get_parent(), atom_b.get_parent()
        if res_a is res_b:
            continue
        key = tuple(sorted((id(atom_a), id(atom_b))))
        if key in seen:
            continue
        seen.add(key)
        d = geom.distance([atom_a], [atom_b])
        results.append(
            HydrogenBond(
                donor_residue=_residue_label(res_a), donor_atom=atom_a.get_name(),
                acceptor_residue=_residue_label(res_b), acceptor_atom=atom_b.get_name(),
                distance=d,
            )
        )
    return sorted(results, key=lambda h: h.distance)


# --- Salt bridges --------------------------------------------------------

@dataclass
class SaltBridge:
    basic_residue: str
    acidic_residue: str
    distance: float


def find_salt_bridges(structure: Structure, cutoff: float = 4.0) -> list[SaltBridge]:
    basic = list(_group_atoms(structure, BASIC_GROUP_ATOMS))
    acidic = list(_group_atoms(structure, ACIDIC_GROUP_ATOMS))
    results: list[SaltBridge] = []
    for basic_res, basic_atoms in basic:
        for acidic_res, acidic_atoms in acidic:
            if basic_res is acidic_res:
                continue
            d = geom.distance(basic_atoms, acidic_atoms)
            if d <= cutoff:
                results.append(
                    SaltBridge(
                        basic_residue=_residue_label(basic_res),
                        acidic_residue=_residue_label(acidic_res),
                        distance=d,
                    )
                )
    return sorted(results, key=lambda s: s.distance)


# --- Hydrophobic contacts -------------------------------------------------

@dataclass
class HydrophobicContact:
    residue_a: str
    residue_b: str
    min_distance: float


def find_hydrophobic_contacts(structure: Structure, cutoff: float = 4.5) -> list[HydrophobicContact]:
    from proteinexplorer.models import is_backbone_atom

    residues = [
        r for r in _protein_residues(structure)
        if r.resname.strip().upper() in HYDROPHOBIC_RESIDUES
    ]
    residue_carbon_atoms = {}
    for r in residues:
        atoms = [
            a for a in r if (a.element or "").strip().upper() == "C"
            and not is_backbone_atom(a.get_name(), ResidueCategory.PROTEIN)
        ]
        if atoms:
            residue_carbon_atoms[id(r)] = (r, atoms)

    entries = list(residue_carbon_atoms.values())
    results: list[HydrophobicContact] = []
    for i in range(len(entries)):
        res_a, atoms_a = entries[i]
        for j in range(i + 1, len(entries)):
            res_b, atoms_b = entries[j]
            coords_a = np.array([a.coord for a in atoms_a])
            coords_b = np.array([a.coord for a in atoms_b])
            diff = coords_a[:, None, :] - coords_b[None, :, :]
            dists = np.sqrt((diff ** 2).sum(axis=-1))
            min_d = float(dists.min())
            if min_d <= cutoff:
                results.append(
                    HydrophobicContact(
                        residue_a=_residue_label(res_a), residue_b=_residue_label(res_b),
                        min_distance=min_d,
                    )
                )
    return sorted(results, key=lambda c: c.min_distance)


# --- Pi-pi and cation-pi ---------------------------------------------------

@dataclass
class PiPiInteraction:
    residue_a: str
    residue_b: str
    centroid_distance: float
    plane_angle: float
    stack_type: str


def _ring_plane(residue, atom_names: tuple[str, ...]):
    atoms = [residue[name] for name in atom_names if name in residue]
    if len(atoms) < 3:
        return None
    return geom.fit_plane(atoms), atoms


def find_pipi_interactions(structure: Structure, cutoff: float = 7.0) -> list[PiPiInteraction]:
    aromatics = []
    for residue in _protein_residues(structure):
        names = AROMATIC_RING_ATOMS.get(residue.resname.strip().upper())
        if not names:
            continue
        plane_info = _ring_plane(residue, names)
        if plane_info is not None:
            aromatics.append((residue, *plane_info))

    results: list[PiPiInteraction] = []
    for i in range(len(aromatics)):
        res_a, plane_a, atoms_a = aromatics[i]
        for j in range(i + 1, len(aromatics)):
            res_b, plane_b, atoms_b = aromatics[j]
            d = float(np.linalg.norm(plane_a.point - plane_b.point))
            if d > cutoff:
                continue
            cos_angle = np.clip(np.dot(plane_a.normal, plane_b.normal), -1.0, 1.0)
            raw_angle = float(np.degrees(np.arccos(cos_angle)))
            angle_deg = min(raw_angle, 180.0 - raw_angle)  # fold into 0-90
            stack_type = "parallel" if angle_deg < 30 else ("t-shaped" if angle_deg > 60 else "intermediate")
            results.append(
                PiPiInteraction(
                    residue_a=_residue_label(res_a), residue_b=_residue_label(res_b),
                    centroid_distance=d, plane_angle=angle_deg, stack_type=stack_type,
                )
            )
    return sorted(results, key=lambda p: p.centroid_distance)


@dataclass
class CationPiInteraction:
    cation_residue: str
    aromatic_residue: str
    distance: float


def find_cationpi_interactions(structure: Structure, cutoff: float = 6.0) -> list[CationPiInteraction]:
    cations = list(_group_atoms(structure, CATION_GROUP_ATOMS))
    aromatics = []
    for residue in _protein_residues(structure):
        names = AROMATIC_RING_ATOMS.get(residue.resname.strip().upper())
        if not names:
            continue
        plane_info = _ring_plane(residue, names)
        if plane_info is not None:
            aromatics.append((residue, plane_info[0]))

    results: list[CationPiInteraction] = []
    for cation_res, cation_atoms in cations:
        cation_point = geom.centroid(cation_atoms)
        for aromatic_res, plane in aromatics:
            if cation_res is aromatic_res:
                continue
            d = float(np.linalg.norm(cation_point - plane.point))
            if d <= cutoff:
                results.append(
                    CationPiInteraction(
                        cation_residue=_residue_label(cation_res),
                        aromatic_residue=_residue_label(aromatic_res),
                        distance=d,
                    )
                )
    return sorted(results, key=lambda c: c.distance)


# --- Disulfide bonds -------------------------------------------------------

@dataclass
class DisulfideBond:
    residue_a: str
    residue_b: str
    distance: float


def find_disulfide_bonds(structure: Structure, cutoff: float = _DISULFIDE_SG_DISTANCE) -> list[DisulfideBond]:
    cys_sg = [
        (r, r["SG"]) for r in _protein_residues(structure)
        if r.resname.strip().upper() == "CYS" and "SG" in r
    ]
    results: list[DisulfideBond] = []
    for i in range(len(cys_sg)):
        res_a, sg_a = cys_sg[i]
        for j in range(i + 1, len(cys_sg)):
            res_b, sg_b = cys_sg[j]
            d = geom.distance([sg_a], [sg_b])
            if d <= cutoff:
                results.append(
                    DisulfideBond(residue_a=_residue_label(res_a), residue_b=_residue_label(res_b), distance=d)
                )
    return sorted(results, key=lambda b: b.distance)


# --- Contact map -----------------------------------------------------------

@dataclass
class ContactMap:
    labels: list[str]
    matrix: np.ndarray  # boolean, symmetric, diagonal False
    distances: np.ndarray  # actual distances used


def contact_map(structure: Structure, atoms=None, cutoff: float = 8.0, mode: str = "ca") -> ContactMap:
    """Residue-residue contact map.

    mode="ca": CA-CA distance (protein residues only).
    mode="heavy": minimum heavy-atom-to-heavy-atom distance between the two
    residues (any residue type -- protein/nucleic/ligand/ion, not water).
    `atoms` restricts which atoms are considered (e.g. output of
    selection.select); defaults to the whole first model minus water.
    """
    model = next(iter(structure))
    if atoms is None:
        atoms = [
            a for a in model.get_atoms()
            if classify_residue(a.get_parent().resname, a.get_parent().id[0]) is not ResidueCategory.WATER
        ]

    residue_atoms: dict[int, tuple] = {}
    for atom in atoms:
        residue = atom.get_parent()
        residue_atoms.setdefault(id(residue), (residue, []))[1].append(atom)

    if mode == "ca":
        entries = []
        for residue, res_atoms in residue_atoms.values():
            ca = next((a for a in res_atoms if a.get_name() == "CA"), None)
            if ca is not None:
                entries.append((residue, [ca]))
    else:
        entries = list(residue_atoms.values())

    entries.sort(key=lambda e: (e[0].get_parent().id, e[0].id[1]))
    labels = [_residue_label(r) for r, _ in entries]
    n = len(entries)
    distances = np.full((n, n), np.inf)
    for i in range(n):
        for j in range(i + 1, n):
            coords_a = np.array([a.coord for a in entries[i][1]])
            coords_b = np.array([a.coord for a in entries[j][1]])
            diff = coords_a[:, None, :] - coords_b[None, :, :]
            d = float(np.sqrt((diff ** 2).sum(axis=-1)).min())
            distances[i, j] = distances[j, i] = d
    np.fill_diagonal(distances, 0.0)
    matrix = (distances <= cutoff) & (distances > 0)
    return ContactMap(labels=labels, matrix=matrix, distances=distances)


# --- Residue interaction network -------------------------------------------

@dataclass
class InteractionEdge:
    residue_a: str
    residue_b: str
    kind: str
    value: float  # distance (or centroid distance) in Angstrom


def interaction_network(
    structure: Structure,
    hbond_cutoff: float = 3.5,
    salt_bridge_cutoff: float = 4.0,
    hydrophobic_cutoff: float = 4.5,
    pipi_cutoff: float = 7.0,
    cationpi_cutoff: float = 6.0,
    disulfide_cutoff: float = _DISULFIDE_SG_DISTANCE,
) -> list[InteractionEdge]:
    """Aggregate every interaction type into one edge list, suitable for
    building a residue interaction graph/network."""
    edges: list[InteractionEdge] = []
    for h in find_hydrogen_bonds(structure, hbond_cutoff):
        edges.append(InteractionEdge(h.donor_residue, h.acceptor_residue, "hbond", h.distance))
    for s in find_salt_bridges(structure, salt_bridge_cutoff):
        edges.append(InteractionEdge(s.basic_residue, s.acidic_residue, "salt_bridge", s.distance))
    for c in find_hydrophobic_contacts(structure, hydrophobic_cutoff):
        edges.append(InteractionEdge(c.residue_a, c.residue_b, "hydrophobic", c.min_distance))
    for p in find_pipi_interactions(structure, pipi_cutoff):
        edges.append(InteractionEdge(p.residue_a, p.residue_b, f"pi_pi_{p.stack_type}", p.centroid_distance))
    for cp in find_cationpi_interactions(structure, cationpi_cutoff):
        edges.append(InteractionEdge(cp.cation_residue, cp.aromatic_residue, "cation_pi", cp.distance))
    for d in find_disulfide_bonds(structure, disulfide_cutoff):
        edges.append(InteractionEdge(d.residue_a, d.residue_b, "disulfide", d.distance))
    return edges
