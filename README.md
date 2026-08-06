# ProteinExplorer

CLI-based structural bioinformatics workbench for **static protein 3D
structure analysis and modeling**. MD trajectory analysis is intentionally
out of scope (see the future MDExplorer). Third tool in the series after
ChemExplorer (small molecules) and BioExplorer (sequences).

## Status

Implemented so far:
- Core data model: PDB/mmCIF I/O, protein/nucleic/water/ion/ligand
  classification, backbone-atom detection (`models.py`, `io.py`)
- Project management under `.proteinexplorer/` (`project.py`)
- CLI skeleton: `prot import` / `prot export` / `prot status`
- `prot info` / `prot descriptor`
- **Common selection language** (`selection.py`): shared by every
  analysis command
- `prot geometry`: `distance`, `angle`, `dihedral`, `backbone-torsions`,
  `rmsd`, `coords`, `distmatrix`
- `prot contact`: `hbond`, `saltbridge`, `hydrophobic`, `pipi`, `cationpi`,
  `disulfide`, `map`, `network`
- `prot secondary`: `dssp` (external binary) or a dependency-free
  `geometric` fallback
- `prot pocket detect`: dependency-free grid/ray-casting cavity detector
- `prot mutate`: point mutation -> new structure in the project.
  `scwrl4` (external) or `cb_only` (dependency-free fallback)
- `prot model`: `gaps` (missing-residue detection), `loop` (crude
  dependency-free gap filler), `homology` (external MODELLER wrapper,
  no fallback)
- `prot compare`:
  - `rmsd` -- over CA atoms common to both structures (matched by chain
    ID + residue number)
  - `tmscore` -- external TMalign/US-align if installed (real structural
    alignment), otherwise a fixed-correspondence fallback using the
    standard TM-score distance formula over the common CA pairs. The
    fallback score is explicitly **not** numerically comparable to real
    TM-align output (different normalization, no alignment search)
  - `secondary` -- Q3-style secondary structure similarity (collapsed to
    H/E/C) over common residues
  - `contact` -- Jaccard similarity between two contact maps
  - `pocket` -- Jaccard similarity between two pockets' lining residues
    (default: each structure's largest pocket)
  - `ligand` -- shared ligand resnames, and RMSD for any ligand present
    as one matching-atom-name instance in each structure

  All `compare` residue-correspondence commands assume matching by
  (chain ID, residue number) -- the right assumption for two states of
  "the same" numbered structure (e.g. a mutant vs. its parent), not for
  true homologs with different numbering.

Remaining spec commands (predict/cluster/annotate/map/plot/view/replay)
are not yet implemented. `geometry`'s convex hull and residue-residue
angle matrix were left out as under-specified.

## Quickstart

```bash
uv sync
uv run prot import structure.pdb --name my_protein
uv run prot status
uv run prot mutate my_protein --chain A --resid 50 --to VAL
uv run prot compare rmsd my_protein my_protein_A50XXXVAL
uv run prot compare tmscore my_protein my_protein_A50XXXVAL
uv run prot compare pocket my_protein my_protein_A50XXXVAL --selection "chain A and resid 40:80"
uv run prot export my_protein out.cif
```

Each project is tracked under a `.proteinexplorer/` directory created
automatically on first `import`, mirroring ChemExplorer's `.chemexplorer/`
and BioExplorer's `.bioexplorer/` layout. Every CLI invocation is logged to
`.proteinexplorer/log.json` for a future `prot replay`.

## Tests

```bash
uv run pytest -q
```
