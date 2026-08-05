"""Secondary structure analysis (spec section "Secondary Structure Analysis").

Two independent assignment methods are provided:

- dssp(): wraps an external mkdssp/dssp binary via Biopython's DSSP class.
  This is the standard H-bond-pattern-based DSSP algorithm (full 8-class
  codes: H/G/I/E/B/T/S/-) and should be preferred whenever it's available.
- geometric(): a dependency-free fallback that classifies each residue's
  (phi, psi) backbone dihedrals into the core alpha-helix / beta-strand
  Ramachandran regions (computed via geometry.backbone_torsions, which is
  already dependency-free), with short-run smoothing. Coarser than DSSP
  (3-class H/E/C, no explicit turns/bridges, no H-bond pattern) but
  requires no external binary.

secondary_structure() picks dssp when possible and falls back to
geometric automatically, tagging which method actually produced the
result so callers/CLI output can be honest about the source.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB.Structure import Structure

from proteinexplorer import geometry as geom
from proteinexplorer.models import ResidueCategory, classify_residue

# Core Ramachandran regions (degrees), deliberately conservative/central so
# borderline residues fall back to coil rather than being over-called.
_ALPHA_PHI = (-100.0, -30.0)
_ALPHA_PSI = (-80.0, 0.0)
_BETA_PHI = (-180.0, -45.0)


def _in_alpha_region(phi: float, psi: float) -> bool:
    return _ALPHA_PHI[0] <= phi <= _ALPHA_PHI[1] and _ALPHA_PSI[0] <= psi <= _ALPHA_PSI[1]


def _in_beta_region(phi: float, psi: float) -> bool:
    return _BETA_PHI[0] <= phi <= _BETA_PHI[1] and (psi >= 90.0 or psi <= -150.0)


class DSSPNotAvailableError(RuntimeError):
    pass


@dataclass
class ResidueSS:
    chain_id: str
    resseq: int
    resname: str
    code: str
    phi: float | None
    psi: float | None
    method: str  # "dssp" or "geometric"


def _protein_residue_list(chain):
    return [
        r for r in chain
        if classify_residue(r.resname, r.id[0]) is ResidueCategory.PROTEIN
    ]


# --- DSSP-based assignment --------------------------------------------

def dssp(structure: Structure, pdb_path: str | Path) -> list[ResidueSS]:
    """Full 8-class DSSP secondary structure via an external mkdssp/dssp
    binary. Requires the structure's source file on disk."""
    binary = shutil.which("mkdssp") or shutil.which("dssp")
    if binary is None:
        raise DSSPNotAvailableError(
            "DSSP executable (mkdssp/dssp) not found on PATH. Install DSSP "
            "(e.g. via conda: `conda install -c salilab dssp`) to use the "
            "H-bond-pattern-based assignment; falling back to the built-in "
            "phi/psi-based classifier is also available (--method geometric)."
        )

    from Bio.PDB.DSSP import DSSP

    model = next(iter(structure))
    dssp_result = DSSP(model, str(pdb_path), dssp=binary)

    results = []
    for key in dssp_result.keys():
        chain_id, res_id = key
        entry = dssp_result[key]
        # Bio.PDB.DSSP entry layout: (dssp_index, aa, ss, rel_asa, phi, psi, ...)
        aa, ss, phi, psi = entry[1], entry[2], entry[4], entry[5]
        resseq = res_id[1] if isinstance(res_id, tuple) else res_id
        results.append(
            ResidueSS(
                chain_id=chain_id, resseq=resseq, resname=aa,
                code=ss if ss and ss != " " else "-",
                phi=phi if phi != 360.0 else None,
                psi=psi if psi != 360.0 else None,
                method="dssp",
            )
        )
    return results


# --- Geometric (phi/psi) fallback -------------------------------------

def geometric(structure: Structure, min_helix_run: int = 4, min_strand_run: int = 2) -> list[ResidueSS]:
    """Dependency-free phi/psi-region classification with short-run
    smoothing (isolated 1-2 residue "helices" or single-residue "strands"
    are demoted to coil, since real secondary structure elements span
    several residues)."""
    model = next(iter(structure))
    raw: list[ResidueSS] = []
    for chain in model:
        residues = _protein_residue_list(chain)
        for residue in residues:
            resseq = residue.id[1]
            try:
                torsions = geom.backbone_torsions(chain, resseq)
            except geom.GeometryError:
                continue
            phi, psi = torsions.phi, torsions.psi
            if phi is None or psi is None:
                code = "C"
            elif _in_alpha_region(phi, psi):
                code = "H"
            elif _in_beta_region(phi, psi):
                code = "E"
            else:
                code = "C"
            raw.append(
                ResidueSS(
                    chain_id=chain.id, resseq=resseq, resname=residue.resname,
                    code=code, phi=phi, psi=psi, method="geometric",
                )
            )

    return _smooth_runs(raw, min_helix_run=min_helix_run, min_strand_run=min_strand_run)


def _smooth_runs(residues: list[ResidueSS], min_helix_run: int, min_strand_run: int) -> list[ResidueSS]:
    if not residues:
        return residues

    result = list(residues)
    i = 0
    while i < len(result):
        j = i
        while j < len(result) and result[j].code == result[i].code and result[j].chain_id == result[i].chain_id:
            j += 1
        run_len = j - i
        code = result[i].code
        min_len = min_helix_run if code == "H" else (min_strand_run if code == "E" else None)
        if min_len is not None and run_len < min_len:
            for k in range(i, j):
                r = result[k]
                result[k] = ResidueSS(
                    chain_id=r.chain_id, resseq=r.resseq, resname=r.resname,
                    code="C", phi=r.phi, psi=r.psi, method=r.method,
                )
        i = j
    return result


# --- Unified entry point ------------------------------------------------

def secondary_structure(
    structure: Structure,
    pdb_path: str | Path | None = None,
    method: str = "auto",
) -> tuple[list[ResidueSS], str]:
    """Return (per-residue assignments, method actually used).

    method="auto" (default): try DSSP if a pdb_path was given and the
    binary is available, otherwise fall back to the geometric classifier.
    method="dssp": require DSSP, raise DSSPNotAvailableError if missing.
    method="geometric": always use the phi/psi classifier.
    """
    if method == "geometric":
        return geometric(structure), "geometric"

    if method == "dssp":
        if pdb_path is None:
            raise ValueError("method='dssp' requires pdb_path")
        return dssp(structure, pdb_path), "dssp"

    if method == "auto":
        if pdb_path is not None:
            try:
                return dssp(structure, pdb_path), "dssp"
            except DSSPNotAvailableError:
                pass
        return geometric(structure), "geometric"

    raise ValueError(f"Unknown method: {method!r} (expected 'auto', 'dssp', or 'geometric')")


def composition(residues: list[ResidueSS]) -> dict[str, float]:
    if not residues:
        return {}
    total = len(residues)
    counts: dict[str, int] = {}
    for r in residues:
        counts[r.code] = counts.get(r.code, 0) + 1
    return {code: count / total for code, count in counts.items()}
