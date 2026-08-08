# ProteinExplorer

CLI-based structural bioinformatics workbench for **static protein 3D
structure analysis and modeling**. MD trajectory analysis is intentionally
out of scope (see the future MDExplorer). Third tool in the series after
ChemExplorer (small molecules) and BioExplorer (sequences).

## Status

**All spec sections are implemented.** Summary:

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
- `prot view`: launch an external 3D viewer on a structure
- `prot replay`: re-run the commands recorded in
  `.proteinexplorer/log.json`. Backs up the current project state to
  `.proteinexplorer_prereplay_<timestamp>`, resets, and replays every
  logged command in-process (via Click's CliRunner). Since structure IDs
  are regenerated on every import, any later command whose argv
  literally referenced an old ID gets that ID automatically rewritten to
  the newly-generated one (matched by structure name). `--from`/`--to`
  select a range, `--skip` overrides the default skip list
  (`view,predict,annotate`), `--continue-on-error` keeps going past a
  failed step, `--no-reset` replays onto the current state instead of
  resetting, `--dry-run` shows the plan without running anything.

`geometry`'s convex hull and residue-residue angle matrix were left out
of `distance_matrix`/`bounding_box` as under-specified in the original
spec discussion; happy to add them on request.

## Quickstart

```bash
uv sync
uv run prot import structure.pdb --name my_protein
uv run prot mutate my_protein --chain A --resid 50 --to VAL
uv run prot mutate my_protein_A50XXXVAL --chain A --resid 87 --to TRP
uv run prot replay --dry-run
uv run prot replay
```

Each project is tracked under a `.proteinexplorer/` directory created
automatically on first `import`, mirroring ChemExplorer's `.chemexplorer/`
and BioExplorer's `.bioexplorer/` layout.

## Tutorial

See [`docs/TUTORIAL.md`](docs/TUTORIAL.md) for a full walkthrough of
every command using real data (1A8O, the HIV-1 capsid C-terminal
domain), including the actual output of each command and three
generated plots. `examples/1a8o/` has the real structure file and
generated outputs; `examples/illustrative/` has small synthetic
structures used for demos that need something 1A8O doesn't have
(a bound metal ion, an enclosed cavity).

## Tests

```bash
uv run pytest -q
# hierarchical clustering and plotting need extras:
uv run pip install -e ".[cluster,viz]" && uv run pytest -q
```
