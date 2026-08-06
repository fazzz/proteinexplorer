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
- `prot model`: `gaps`, `loop` (crude dependency-free gap filler),
  `homology` (external MODELLER wrapper, no fallback)
- `prot compare`: `rmsd`, `tmscore` (external TMalign/US-align, or a
  fixed-correspondence fallback), `secondary`, `contact`, `pocket`,
  `ligand`
- `prot cluster`: ensemble clustering by pairwise RMSD, same two-tier
  pattern as ChemExplorer/BioExplorer's clustering commands:
  - `ensemble` -- cluster several structures already in the project
  - `models` -- cluster the MODEL records within one multi-model file
    (e.g. an NMR ensemble)
  - `--method greedy` (default): pure-Python CD-HIT-style incremental
    clustering, no extra dependency
  - `--method hierarchical`: scipy-based agglomerative clustering,
    needs `pip install -e ".[cluster]"`
  - Both report a medoid as each cluster's representative

Remaining spec commands (predict/annotate/map/plot/view/replay) are not
yet implemented. `geometry`'s convex hull and residue-residue angle
matrix were left out as under-specified.

## Quickstart

```bash
uv sync
uv run prot import structure.pdb --name my_protein
uv run prot status
uv run prot mutate my_protein --chain A --resid 50 --to VAL
uv run prot compare rmsd my_protein my_protein_A50XXXVAL
uv run prot cluster ensemble my_protein my_protein_A50XXXVAL --threshold 2.0
uv run prot export my_protein out.cif
```

For an NMR-style multi-model file:

```bash
uv run prot import ensemble.pdb --name ensemble
uv run prot cluster models ensemble --threshold 2.0
```

Each project is tracked under a `.proteinexplorer/` directory created
automatically on first `import`, mirroring ChemExplorer's `.chemexplorer/`
and BioExplorer's `.bioexplorer/` layout. Every CLI invocation is logged to
`.proteinexplorer/log.json` for a future `prot replay`.

## Tests

```bash
uv run pytest -q
# hierarchical clustering tests need scipy:
uv run pip install -e ".[cluster]" && uv run pytest -q
```
