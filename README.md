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
  ratio, disulfide bond count, secondary structure composition (needs
  external `mkdssp`/`dssp`, degrades gracefully if absent)
- **Common selection language** (`selection.py`): `protein`/`nucleic`/
  `water`/`ion`/`ligand`, `backbone`/`sidechain`, `chain X`, `resid a:b`,
  `resname X`, `atom X`, `within N <sel>`, boolean `and`/`or`/`not` with
  parentheses -- shared by every analysis command
- `prot geometry`: `distance`, `angle`, `dihedral`, `backbone-torsions`
  (phi/psi/omega + side chain chi1-4), `rmsd` (Kabsch-fit or raw), `coords`
  (centroid/COM/bounding box/radius of gyration/plane fit/principal
  axes/moment of inertia), `distmatrix` (pairwise centroid distance matrix)
- `prot contact`: `hbond` (heavy-atom N/O...N/O, no hydrogens needed),
  `saltbridge` (Arg/Lys/His vs Asp/Glu charged-group distance), `hydrophobic`
  (sidechain carbon-carbon), `pipi` (aromatic ring stacking, with
  parallel/t-shaped/intermediate classification via ring-plane angle),
  `cationpi` (Arg/Lys vs aromatic ring), `disulfide` (Cys SG-SG),
  `map` (residue-residue contact map, CA-CA or min-heavy-atom mode),
  `network` (all interaction types combined into one edge list)

Remaining spec commands (secondary/pocket/mutate/model/predict/compare/
cluster/annotate/map/plot/view/replay) are not yet implemented.
`geometry`'s convex hull and residue-residue angle matrix were left out as
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
uv run prot contact saltbridge my_protein
uv run prot contact map my_protein --mode heavy --cutoff 5
uv run prot contact network my_protein
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
