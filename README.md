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
  ratio, disulfide bond count, secondary structure composition
- **Common selection language** (`selection.py`): `protein`/`nucleic`/
  `water`/`ion`/`ligand`, `backbone`/`sidechain`, `chain X`, `resid a:b`,
  `resname X`, `atom X`, `within N <sel>`, boolean `and`/`or`/`not` with
  parentheses -- shared by every analysis command
- `prot geometry`: `distance`, `angle`, `dihedral`, `backbone-torsions`
  (phi/psi/omega + side chain chi1-4), `rmsd`, `coords` (centroid/COM/
  bounding box/radius of gyration/plane fit/principal axes/moment of
  inertia), `distmatrix`
- `prot contact`: `hbond`, `saltbridge`, `hydrophobic`, `pipi`, `cationpi`,
  `disulfide`, `map`, `network`
- `prot secondary`: per-residue secondary structure + composition
  (`dssp` external binary, or a dependency-free `geometric` phi/psi
  fallback; `auto` picks whichever is available)
- `prot pocket detect`: dependency-free grid/ray-casting cavity detector
  (LIGSITE-style approximation -- no fpocket/P2Rank available in this
  environment). Reports volume, a rough surface-area estimate, lining
  residues, hydrophobic fraction, and a simple **heuristic** (not a
  trained model) druggability score. Large structures need `--selection`
  to keep the search grid a manageable size; the command errors out with
  a clear message (and the point count) instead of hanging.

Remaining spec commands (mutate/model/predict/compare/cluster/annotate/
map/plot/view/replay) are not yet implemented. `geometry`'s convex hull
and residue-residue angle matrix were left out as under-specified.

## Quickstart

```bash
uv sync
uv run prot import structure.pdb --name my_protein
uv run prot status
uv run prot info my_protein
uv run prot descriptor my_protein
uv run prot geometry distance my_protein "chain A and resid 10 and atom CA" "chain A and resid 50 and atom CA"
uv run prot contact hbond my_protein
uv run prot secondary my_protein
uv run prot pocket detect my_protein --selection "chain A and resid 40:80"
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
