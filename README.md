# ProteinExplorer

CLI-based structural bioinformatics workbench for **static protein 3D
structure analysis and modeling**. MD trajectory analysis is intentionally
out of scope (see the future MDExplorer). Third tool in the series after
ChemExplorer (small molecules) and BioExplorer (sequences).

## Status

**All spec sections are implemented, plus a `search` command added
after the fact (Foldseek integration).** Summary:

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
- `prot search`: structural similarity search via Foldseek
  (external-tool-only, not on PyPI -- no dependency-free fallback exists
  for large-scale structural database search):
  - `foldseek` -- search one structure against a Foldseek database
    (`--target-db`), a directory of structure files (`--target-dir`), or
    every other structure already in the project (`--against-project`)
  - `createdb` -- build a persistent Foldseek database from a directory
    of structure files, for repeated searches
- `prot fix`: structure fixing/cleanup via PDBFixer (`pip install -e
  ".[fix]"` -- free/open-source and pip-installable, unlike
  Scwrl4/MODELLER/Foldseek, so this is treated as a normal optional
  extra rather than an external-tool-only wrapper):
  - `report` -- shows what PDBFixer would find (missing atoms within
    existing residues, missing whole residues from SEQRES, nonstandard
    residues) without changing anything
  - `apply` -- runs the requested repair steps and saves the result as a
    new structure in the project. Adds missing atoms within existing
    residues by default (a real gap `prot mutate`/`prot model` don't
    cover); replacing nonstandard residues (e.g. MSE -> MET), whole
    missing-residue insertion, heterogen removal, and hydrogen addition
    are all available via flags

  **Overlap with `prot model`, resolved on purpose:** PDBFixer's own
  missing-residue detection only sees gaps recorded in the file's SEQRES
  header -- a numbering-only gap with no SEQRES is invisible to it
  (verified directly: it reports nothing for this project's own
  `gapped.pdb` fixture, while `prot model gaps`, which works from
  residue numbering alone, correctly finds the 3-residue gap). Use
  `prot model gaps` for detection; `prot fix apply
  --add-missing-residues` can fill gaps PDBFixer's own detection *does*
  find, with more realistic template-based geometry than `prot model
  loop`'s straight-line placeholder trace.

  One more thing worth knowing: PDBFixer's PDB writer puts water into
  its own chain on output, so a single-chain input can come back as two
  chains purely because of solvent -- not a bug in this wrapper, just
  PDBFixer's own convention.
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
uv run pip install -e ".[cluster,viz,fix]" && uv run pytest -q
```
