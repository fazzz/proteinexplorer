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
  (LIGSITE-style approximation)
- `prot mutate`: point mutation -> new structure in the project. `scwrl4`
  (external, full rotamer optimization) or `cb_only` (dependency-free:
  backbone kept, idealized virtual C-beta only, honestly reports it does
  not build the rest of the side chain)
- `prot model`:
  - `gaps` -- detect numbering discontinuities (missing residues), no
    external tool needed
  - `loop` -- crude dependency-free gap filler: linear CA interpolation
    with idealized local backbone geometry between the flanking anchor
    residues. No clash checking or energy minimization -- a placeholder
    trace to refine further, not a real loop model. Saves the result as
    a new structure in the project
  - `homology` -- wraps an external MODELLER installation (license
    required from https://salilab.org/modeller/). No dependency-free
    fallback exists for homology modeling, so this errors out clearly
    when the `modeller` package isn't installed

Sidechain rebuilding/repacking is intentionally not a separate command --
use `prot mutate --to <same residue>` for that (it's the same operation as
a point mutation with an unchanged target identity).

Remaining spec commands (predict/compare/cluster/annotate/map/plot/view/
replay) are not yet implemented. `geometry`'s convex hull and
residue-residue angle matrix were left out as under-specified.

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
uv run prot mutate my_protein --chain A --resid 50 --to VAL
uv run prot model gaps my_protein
uv run prot model loop my_protein --chain A --start 41 --end 43 --sequence GLY
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
