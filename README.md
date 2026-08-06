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
- `prot predict`: `colabfold` / `alphafold`, external tools only
- `prot annotate`: `metal-sites` / `metadata` (built-in) and `uniprot` /
  `pfam` (external REST lookups)
- `prot map`: `pocket` / `mutation` / `domain` / `conservation` coloring
  scripts for PyMOL/ChimeraX/VMD
- `prot view`: launch an external 3D viewer (PyMOL/ChimeraX/VMD) on a
  structure, optionally running a script (e.g. from `prot map`)
  afterward. Starts the viewer as a background process; no
  dependency-free substitute exists for this command.

Remaining spec command: `replay` (workflow re-execution from
`.proteinexplorer/log.json`, which every command already writes to) is
not yet implemented. `geometry`'s convex hull and residue-residue angle
matrix were left out as under-specified.

## Quickstart

```bash
uv sync
uv run prot import structure.pdb --name my_protein
uv run prot status
uv run prot map mutation my_protein mut.pml --residue A/50
uv run prot view my_protein --tool pymol --script mut.pml
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
