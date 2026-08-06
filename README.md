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
- `prot compare`: `rmsd`, `tmscore`, `secondary`, `contact`, `pocket`,
  `ligand`
- `prot cluster`: `ensemble` / `models`, `--method greedy` (dependency-free)
  or `hierarchical` (needs `[cluster]` extra)
- `prot plot` (needs `[viz]` extra): `ramachandran`, `contact-map`,
  `secondary`
- `prot predict`: structure prediction via external tools only --
  `colabfold` (wraps `colabfold_batch`) and `alphafold` (wraps
  `run_alphafold.sh`). No dependency-free fallback exists for ab initio
  prediction, so this is a thin, honest wrapper: it errors out clearly
  when the tool isn't installed rather than attempting a fake substitute.
  `--import-name` imports the top-ranked model into the project on success.

Remaining spec commands (annotate/map/view/replay) are not yet
implemented. `geometry`'s convex hull and residue-residue angle matrix
were left out as under-specified.

## Quickstart

```bash
uv sync
uv run prot import structure.pdb --name my_protein
uv run prot status
uv run prot predict colabfold "MKVLTA..." --output-dir out/ --import-name predicted
```

Each project is tracked under a `.proteinexplorer/` directory created
automatically on first `import`, mirroring ChemExplorer's `.chemexplorer/`
and BioExplorer's `.bioexplorer/` layout. Every CLI invocation is logged to
`.proteinexplorer/log.json` for a future `prot replay`.

## Tests

```bash
uv run pytest -q
# hierarchical clustering and plotting need extras:
uv run pip install -e ".[cluster,viz]" && uv run pytest -q
```
