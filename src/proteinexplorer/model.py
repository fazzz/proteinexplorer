"""Structure modeling (spec section "Modeling").

Covers missing-residue detection, a dependency-free crude loop filler,
and an external-tool wrapper for homology modeling. Point mutation and
sidechain rebuilding are intentionally NOT duplicated here -- both are
already `prot mutate` (rebuild a residue's side chain by "mutating" it to
its own identity, or to a different one).

- find_gaps(): purely reads residue numbering, no external tool needed.
- fill_loop_linear(): a dependency-free placeholder-backbone filler. This
  is deliberately crude (linear CA interpolation between anchors, idealized
  local backbone geometry, no clash checking, no energy minimization) --
  it produces a topologically continuous trace to fill a gap, not a
  scientifically validated loop model. Real loop modeling needs a tool
  like MODELLER or Rosetta loop modeling, which this module does not
  attempt to reimplement.
- homology_model(): wraps the MODELLER Python package (automodel), which
  requires a license from https://salilab.org/modeller/. No dependency-free
  fallback exists for homology modeling -- template selection, alignment,
  and energy-based model building can't be meaningfully approximated
  without it, so this raises a clear error when MODELLER isn't installed
  rather than attempting a fake substitute.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB.Structure import Structure

from proteinexplorer.models import ResidueCategory, classify_residue
from proteinexplorer.mutate import MutationError, normalize_resname


# --- Missing residue detection -------------------------------------------

@dataclass
class Gap:
    chain_id: str
    prev_resseq: int  # last resolved residue before the gap
    next_resseq: int  # first resolved residue after the gap
    length: int  # number of missing residues implied by the numbering gap


def find_gaps(structure: Structure) -> list[Gap]:
    """Detect numbering discontinuities among a chain's protein residues.

    This only sees what's in the coordinate file: it flags jumps in
    residue sequence numbers (e.g. 41 -> 45 implies 3 missing residues),
    it does not cross-reference SEQRES records, so residues absent from
    both ATOM and SEQRES won't be reported as anything other than "the
    chain doesn't cover that range."
    """
    model = next(iter(structure))
    gaps: list[Gap] = []
    for chain in model:
        resseqs = sorted(
            r.id[1] for r in chain
            if classify_residue(r.resname, r.id[0]) is ResidueCategory.PROTEIN
        )
        for prev, nxt in zip(resseqs, resseqs[1:]):
            if nxt - prev > 1:
                gaps.append(Gap(chain_id=chain.id, prev_resseq=prev, next_resseq=nxt, length=nxt - prev - 1))
    return gaps


# --- Crude dependency-free loop filler -------------------------------------

@dataclass
class LoopFillResult:
    chain_id: str
    start_resseq: int
    end_resseq: int
    residues_added: list[str]
    note: str


def fill_loop_linear(
    structure: Structure,
    chain_id: str,
    start: int,
    end: int,
    sequence: str | None = None,
) -> LoopFillResult:
    """Fill resid range [start, end] with placeholder residues on a
    straight-line CA trace between the flanking anchor residues
    (start-1 and end+1), with idealized local backbone geometry around
    each interpolated CA.

    This is a crude filler, not a real loop model: no clash checking, no
    energy minimization, no use of real bond lengths/angles beyond a
    generic idealized local frame. It exists to give a continuous,
    topologically-plausible starting trace when no external loop-modeling
    tool is available -- treat the result as a placeholder to refine
    further (e.g. with MODELLER/Rosetta), not a finished model.
    """
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Residue import Residue

    model = next(iter(structure))
    if chain_id not in model:
        raise MutationError(f"No chain {chain_id!r} in this structure")
    chain = model[chain_id]

    if end < start:
        raise MutationError(f"end ({end}) must be >= start ({start})")
    n_missing = end - start + 1

    if sequence is not None:
        if len(sequence) != n_missing:
            raise MutationError(
                f"--sequence has {len(sequence)} residues but the gap is {n_missing} long"
            )
        target_resnames = [normalize_resname(c) for c in sequence]
    else:
        target_resnames = ["ALA"] * n_missing

    prev_id = (" ", start - 1, " ")
    next_id = (" ", end + 1, " ")
    if prev_id not in chain or next_id not in chain:
        raise MutationError(
            f"Both anchor residues ({chain_id}/{start - 1} and {chain_id}/{end + 1}) "
            f"must already exist in the structure to fill the gap between them"
        )
    prev_res, next_res = chain[prev_id], chain[next_id]
    if "CA" not in prev_res or "CA" not in next_res:
        raise MutationError("Both anchor residues need a CA atom")

    prev_ca = prev_res["CA"].coord.astype(float)
    next_ca = next_res["CA"].coord.astype(float)
    chain_vec = next_ca - prev_ca
    chain_dir = chain_vec / np.linalg.norm(chain_vec)
    # a fixed arbitrary vector not parallel to chain_dir, for a stable local frame
    ref = np.array([0.0, 0.0, 1.0]) if abs(chain_dir[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    side = np.cross(chain_dir, ref)
    side = side / np.linalg.norm(side)

    added_names: list[str] = []
    for i, resname in enumerate(target_resnames, start=1):
        t = i / (n_missing + 1)
        ca = prev_ca + chain_vec * t
        n_coord = ca - chain_dir * 1.0 + side * 0.3
        c_coord = ca + chain_dir * 1.0 - side * 0.3
        o_coord = c_coord + side * 1.0

        resseq = start + i - 1
        residue = Residue((" ", resseq, " "), resname, getattr(chain, "segid", ""))
        residue.add(Atom("N", n_coord, 0.0, 1.0, " ", "N", 0, element="N"))
        residue.add(Atom("CA", ca, 0.0, 1.0, " ", "CA", 0, element="C"))
        residue.add(Atom("C", c_coord, 0.0, 1.0, " ", "C", 0, element="C"))
        residue.add(Atom("O", o_coord, 0.0, 1.0, " ", "O", 0, element="O"))
        chain.add(residue)
        added_names.append(f"{resname}{resseq}")

    # Bio.PDB's internal child_list isn't automatically kept sorted by
    # resseq when adding out of order; re-sort so downstream iteration
    # (backbone_torsions, secondary structure, ...) sees sequence order.
    chain.child_list.sort(key=lambda r: r.id[1])

    return LoopFillResult(
        chain_id=chain_id, start_resseq=start, end_resseq=end,
        residues_added=added_names,
        note=(
            "Crude placeholder trace: linear CA interpolation with idealized "
            "local backbone geometry, no clash checking or energy minimization. "
            "Refine with a real loop-modeling tool (e.g. MODELLER, Rosetta) "
            "before relying on this region."
        ),
    )


# --- Homology modeling (external tool only) --------------------------------

class ModellerNotAvailableError(RuntimeError):
    pass


def homology_model(
    alignment_pir_path: str | Path,
    template_codes: list[str],
    target_code: str,
    template_search_dir: str | Path,
    output_dir: str | Path,
    n_models: int = 1,
) -> list[Path]:
    """Build homology model(s) via MODELLER's automodel.

    Requires the `modeller` Python package (a license from
    https://salilab.org/modeller/ is needed to install it) -- there is no
    dependency-free fallback for homology modeling, so this raises
    ModellerNotAvailableError with installation guidance when the package
    is missing rather than attempting a fake substitute.

    `alignment_pir_path` must be a PIR-format alignment containing both
    the target sequence (as `target_code`) and the template(s)
    (`template_codes`), following MODELLER's standard alignment format.
    Template structure files must be discoverable in
    `template_search_dir` (MODELLER's default lookup behaviour).
    """
    try:
        from modeller import Environ
        from modeller.automodel import AutoModel
    except ImportError as exc:
        raise ModellerNotAvailableError(
            "The `modeller` Python package is not installed. Homology modeling "
            "has no dependency-free fallback -- install MODELLER (license "
            "required from https://salilab.org/modeller/) to use this command."
        ) from exc

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    env = Environ()
    env.io.atom_files_directory = [str(template_search_dir)]
    a = AutoModel(
        env,
        alnfile=str(alignment_pir_path),
        knowns=template_codes,
        sequence=target_code,
    )
    a.starting_model = 1
    a.ending_model = n_models
    a.make()

    return [Path(m["name"]) for m in a.outputs if m["failure"] is None]
