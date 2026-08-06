"""Structure comparison (spec section "Compare").

Every function here compares two already-loaded structures. Several of
them need a residue *correspondence* between the two structures; rather
than performing a real sequence-independent structural alignment (what
TM-align actually does via dynamic programming + iterative superposition),
this module uses a much simpler, explicit convention: residues correspond
by (chain_id, residue_seqnum). That's the right assumption for comparing
two states of "the same" numbered structure (e.g. a mutant vs. its parent,
two experimental structures of one protein with consistent numbering) and
the wrong one for true homologs with different numbering -- for that,
install TM-align/US-align and use the external path.

- rmsd(): thin pass-through to geometry.rmsd on the common-label CA atoms.
- tm_score(): external TMalign/USalign if installed (real sequence-
  independent structural alignment); otherwise a fixed-correspondence
  fallback that computes the standard TM-score distance-scoring formula
  over the common (chain, resseq) CA pairs. The fallback's score is NOT
  numerically comparable to real TM-align output (different, usually
  smaller, normalization length and no search for the optimal alignment)
  -- this is called out in the return value and CLI output.
- secondary_structure_similarity(): Q3-style fraction-identical over
  common residues, both sides collapsed to the 3-class H/E/C scheme so
  a DSSP run can be compared against a geometric-fallback run.
- contact_similarity(): Jaccard similarity between two contact maps,
  restricted to residue pairs where both residues exist in both structures.
- pocket_overlap(): Jaccard similarity between two pockets' lining-residue
  label sets (defaults to comparing each structure's largest pocket).
- ligand_comparison(): which ligand resnames are shared, and RMSD for any
  shared ligand present as exactly one instance in each structure with
  matching atom names.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB.Structure import Structure

from proteinexplorer import geometry as geom
from proteinexplorer.models import ResidueCategory, classify_residue


def _residue_label(residue) -> tuple[str, int]:
    return (residue.get_parent().id, residue.id[1])


def _ca_by_label(structure: Structure) -> dict[tuple[str, int], object]:
    model = next(iter(structure))
    result = {}
    for chain in model:
        for residue in chain:
            if classify_residue(residue.resname, residue.id[0]) is ResidueCategory.PROTEIN and "CA" in residue:
                result[(chain.id, residue.id[1])] = residue["CA"]
    return result


def common_ca_atoms(structure_a: Structure, structure_b: Structure) -> tuple[list, list, list[tuple[str, int]]]:
    """CA atoms present at the same (chain_id, resseq) in both structures,
    in a consistent sorted order."""
    ca_a = _ca_by_label(structure_a)
    ca_b = _ca_by_label(structure_b)
    labels = sorted(set(ca_a) & set(ca_b))
    return [ca_a[label] for label in labels], [ca_b[label] for label in labels], labels


class CompareError(ValueError):
    pass


# --- RMSD ------------------------------------------------------------------

def rmsd(structure_a: Structure, structure_b: Structure, fit: bool = True) -> tuple[float, int]:
    """RMSD over CA atoms common to both structures (by chain_id+resseq).
    Returns (rmsd, n_atoms_used)."""
    atoms_a, atoms_b, labels = common_ca_atoms(structure_a, structure_b)
    if len(labels) < 2:
        raise CompareError("Fewer than 2 common (chain, resid) CA atoms between the two structures")
    return geom.rmsd(atoms_a, atoms_b, fit=fit), len(labels)


# --- TM-score ----------------------------------------------------------

@dataclass
class TMScoreResult:
    score: float
    method: str  # "tmalign" or "fallback"
    n_residues: int
    note: str


def _d0(length: int) -> float:
    if length < 15:
        return 0.5
    value = 1.24 * (length - 15) ** (1 / 3) - 1.8
    return max(value, 0.5)


def tm_score_fallback(structure_a: Structure, structure_b: Structure) -> TMScoreResult:
    atoms_a, atoms_b, labels = common_ca_atoms(structure_a, structure_b)
    if len(labels) < 3:
        raise CompareError("Fewer than 3 common (chain, resid) CA atoms between the two structures")

    from Bio.PDB.Superimposer import Superimposer

    sup = Superimposer()
    sup.set_atoms(list(atoms_a), list(atoms_b))
    rot, tran = sup.rotran
    coords_b = np.array([a.coord for a in atoms_b]) @ rot + tran
    coords_a = np.array([a.coord for a in atoms_a])
    distances = np.sqrt(((coords_a - coords_b) ** 2).sum(axis=1))

    L = len(labels)
    d0 = _d0(L)
    score = float(np.mean(1.0 / (1.0 + (distances / d0) ** 2)))
    return TMScoreResult(
        score=score, method="fallback", n_residues=L,
        note=(
            "Fixed-correspondence fallback (matched by chain+resid, not a real "
            "structural alignment): normalized by the common-residue count, not "
            "the standard reference-length convention, so this is NOT numerically "
            "comparable to real TM-align/TM-score output. Install TMalign or "
            "US-align for a real score."
        ),
    )


def tm_score_external(structure_a_path: str | Path, structure_b_path: str | Path) -> TMScoreResult:
    binary = shutil.which("TMalign") or shutil.which("USalign") or shutil.which("tmalign")
    if binary is None:
        raise RuntimeError("TMalign/USalign not found on PATH")

    result = subprocess.run(
        [binary, str(structure_a_path), str(structure_b_path)],
        capture_output=True, text=True, timeout=300,
    )
    if result.returncode != 0:
        raise RuntimeError(f"{Path(binary).name} failed: {result.stderr or result.stdout}")

    score = None
    for line in result.stdout.splitlines():
        if line.startswith("TM-score") and "Chain_1" in line:
            score = float(line.split("=")[1].split()[0])
            break
    if score is None:
        raise RuntimeError(f"Could not parse TM-score from {Path(binary).name} output")

    return TMScoreResult(
        score=score, method="tmalign", n_residues=-1,
        note=f"Real structural alignment via {Path(binary).name}.",
    )


def tm_score(
    structure_a: Structure,
    structure_b: Structure,
    structure_a_path: str | Path | None = None,
    structure_b_path: str | Path | None = None,
    method: str = "auto",
) -> TMScoreResult:
    if method == "fallback":
        return tm_score_fallback(structure_a, structure_b)
    if method == "tmalign":
        if structure_a_path is None or structure_b_path is None:
            raise ValueError("method='tmalign' requires structure_a_path and structure_b_path")
        return tm_score_external(structure_a_path, structure_b_path)
    if method == "auto":
        if structure_a_path is not None and structure_b_path is not None:
            try:
                return tm_score_external(structure_a_path, structure_b_path)
            except RuntimeError:
                pass
        return tm_score_fallback(structure_a, structure_b)
    raise ValueError(f"Unknown method: {method!r}")


# --- Secondary structure similarity -----------------------------------

def _collapse_3class(code: str) -> str:
    if code in ("H", "G", "I"):
        return "H"
    if code in ("E", "B"):
        return "E"
    return "C"


def secondary_structure_similarity(
    structure_a: Structure, structure_b: Structure,
    method: str = "auto",
) -> tuple[float, int]:
    """Q3-style fraction of common (chain, resid) residues whose secondary
    structure class matches, both sides collapsed to H/E/C so a DSSP
    assignment can be compared against a geometric-fallback assignment."""
    from proteinexplorer import secondary as sec

    residues_a, _ = sec.secondary_structure(structure_a, method=method if method != "auto" else "geometric")
    residues_b, _ = sec.secondary_structure(structure_b, method=method if method != "auto" else "geometric")

    map_a = {(r.chain_id, r.resseq): _collapse_3class(r.code) for r in residues_a}
    map_b = {(r.chain_id, r.resseq): _collapse_3class(r.code) for r in residues_b}
    common = sorted(set(map_a) & set(map_b))
    if not common:
        raise CompareError("No common (chain, resid) residues between the two structures")

    matches = sum(1 for label in common if map_a[label] == map_b[label])
    return matches / len(common), len(common)


# --- Contact similarity --------------------------------------------------

def contact_similarity(
    structure_a: Structure, structure_b: Structure,
    mode: str = "ca", cutoff: float = 8.0,
) -> tuple[float, int, int]:
    """Jaccard similarity between the two structures' contact maps,
    restricted to residue pairs where both residues exist in both
    structures. Returns (jaccard, n_shared_contacts, n_union_contacts)."""
    from proteinexplorer import contact as ct

    cm_a = ct.contact_map(structure_a, mode=mode, cutoff=cutoff)
    cm_b = ct.contact_map(structure_b, mode=mode, cutoff=cutoff)

    common_labels = sorted(set(cm_a.labels) & set(cm_b.labels))
    if len(common_labels) < 2:
        raise CompareError("Fewer than 2 common labeled residues between the two structures")

    index_a = {label: i for i, label in enumerate(cm_a.labels)}
    index_b = {label: i for i, label in enumerate(cm_b.labels)}

    contacts_a = set()
    contacts_b = set()
    for i in range(len(common_labels)):
        for j in range(i + 1, len(common_labels)):
            li, lj = common_labels[i], common_labels[j]
            if cm_a.matrix[index_a[li], index_a[lj]]:
                contacts_a.add((li, lj))
            if cm_b.matrix[index_b[li], index_b[lj]]:
                contacts_b.add((li, lj))

    union = contacts_a | contacts_b
    intersection = contacts_a & contacts_b
    jaccard = len(intersection) / len(union) if union else 1.0
    return jaccard, len(intersection), len(union)


# --- Pocket overlap ------------------------------------------------------

def pocket_overlap(
    structure_a: Structure, structure_b: Structure,
    pocket_index_a: int = 1, pocket_index_b: int = 1,
    atoms_a: list | None = None, atoms_b: list | None = None,
    **pocket_kwargs,
) -> tuple[float, list, list]:
    """Jaccard similarity between two pockets' lining-residue labels
    (default: each structure's largest pocket, i.e. #1). `atoms_a`/
    `atoms_b` optionally restrict each structure's own search region
    (e.g. from a selection); `pocket_kwargs` (spacing, padding, ...) are
    shared by both searches. Returns
    (jaccard, pocket_a_residues, pocket_b_residues)."""
    from proteinexplorer import pocket as pk

    pockets_a = pk.find_pockets(structure_a, atoms=atoms_a, **pocket_kwargs)
    pockets_b = pk.find_pockets(structure_b, atoms=atoms_b, **pocket_kwargs)

    pocket_a = next((p for p in pockets_a if p.id == pocket_index_a), None)
    pocket_b = next((p for p in pockets_b if p.id == pocket_index_b), None)
    if pocket_a is None or pocket_b is None:
        raise CompareError("Requested pocket index not found in one or both structures")

    set_a, set_b = set(pocket_a.residues), set(pocket_b.residues)
    union = set_a | set_b
    jaccard = len(set_a & set_b) / len(union) if union else 1.0
    return jaccard, pocket_a.residues, pocket_b.residues


# --- Ligand comparison -----------------------------------------------------

@dataclass
class LigandComparison:
    common_resnames: list[str]
    only_in_a: list[str]
    only_in_b: list[str]
    rmsd_by_resname: dict[str, float]  # only populated where a 1:1 atom-name match was possible


def _ligand_residues(structure: Structure) -> dict[str, list]:
    model = next(iter(structure))
    result: dict[str, list] = {}
    for chain in model:
        for residue in chain:
            if classify_residue(residue.resname, residue.id[0]) is ResidueCategory.LIGAND:
                result.setdefault(residue.resname.strip().upper(), []).append(residue)
    return result


def ligand_comparison(structure_a: Structure, structure_b: Structure, fit: bool = True) -> LigandComparison:
    ligands_a = _ligand_residues(structure_a)
    ligands_b = _ligand_residues(structure_b)
    names_a, names_b = set(ligands_a), set(ligands_b)
    common = sorted(names_a & names_b)

    rmsd_by_resname: dict[str, float] = {}
    for name in common:
        instances_a, instances_b = ligands_a[name], ligands_b[name]
        if len(instances_a) != 1 or len(instances_b) != 1:
            continue
        res_a, res_b = instances_a[0], instances_b[0]
        names_atoms_a = {a.get_name() for a in res_a}
        names_atoms_b = {a.get_name() for a in res_b}
        if names_atoms_a != names_atoms_b:
            continue
        ordered_names = sorted(names_atoms_a)
        atoms_a = [res_a[n] for n in ordered_names]
        atoms_b = [res_b[n] for n in ordered_names]
        rmsd_by_resname[name] = geom.rmsd(atoms_a, atoms_b, fit=fit)

    return LigandComparison(
        common_resnames=common,
        only_in_a=sorted(names_a - names_b),
        only_in_b=sorted(names_b - names_a),
        rmsd_by_resname=rmsd_by_resname,
    )
