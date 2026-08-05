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
- `prot info`: detailed per-structure summary
- `prot descriptor`: molecular weight, atom/residue/chain/ligand/water
  counts, SASA, radius of gyration, CA-CA contact density, hydrophobic
  ratio, disulfide bond count, secondary structure composition (via
  `secondary.py`, DSSP if available else the geometric fallback)
- **Common selection language** (`selection.py`): `protein`/`nucleic`/
  `water`/`ion`/`ligand`, `backbone`/`sidechain`, `chain X`, `resid a:b`,
  `resname X`, `atom X`, `within N <sel>`, boolean `and`/`or`/`not` with
  parentheses -- shared by every analysis command
- `prot geometry`: `distance`, `angle`, `dihedral`, `backbone-torsions`
  (phi/psi/omega + side chain chi1-4), `rmsd` (Kabsch-fit or raw), `coords`
  (centroid/COM/bounding box/radius of gyration/plane fit/principal
  axes/moment of inertia), `distmatrix` (pairwise centroid distance matrix)
- `prot contact`: `hbond`, `saltbridge`, `hydrophobic`, `pipi`, `cationpi`,
  `disulfide`, `map` (residue-residue contact map), `network` (all
  interaction types combined into one edge list)
- `prot secondary`: per-residue secondary structure + composition, via
  `secondary.py`. Two methods: `dssp` (external mkdssp/dssp binary, full
  8-class H/G/I/E/B/T/S/- codes) and `geometric` (dependency-free phi/psi
  Ramachandran-region classifier with short-run smoothing, 3-class H/E/C).
  `--method auto` (default) prefers DSSP and falls back to geometric when
  the binary isn't installed.

Remaining spec commands (pocket/mutate/model/predict/compare/cluster/
annotate/map/plot/view/replay) are not yet implemented. `geometry`'s
convex hull and residue-residue angle matrix were left out as
under-specified; happy to add them on request.

## Quickstart

```bash
uv sync
uv run prot import structure.pdb --name my_protein
uv run prot status
uv run prot info my_protein
uv run prot descriptor my_protein
uv run prot geometry distance my_protein "chain A and resid 10 and atom CA" "chain A and resid 50 and atom CA"
uv run prot geometry backbone-torsions my_protein --chain A --resid 25
uv run prot geometry coords my_protein "chain A and backbone"
uv run prot contact hbond my_protein
uv run prot contact map my_protein --mode heavy --cutoff 5
uv run prot secondary my_protein
uv run prot secondary my_protein --method dssp   # requires mkdssp/dssp on PATH
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
