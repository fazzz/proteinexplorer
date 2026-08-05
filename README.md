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
  counts, SASA (Biopython Shrake-Rupley), radius of gyration, CA-CA contact
  density, hydrophobic ratio, disulfide bond count, secondary structure
  composition (needs external `mkdssp`/`dssp`, degrades gracefully if absent)
- **Common selection language** (`selection.py`): `protein`/`nucleic`/
  `water`/`ion`/`ligand`, `backbone`/`sidechain`, `chain X`, `resid a:b`,
  `resname X`, `atom X`, `within N <sel>`, boolean `and`/`or`/`not` with
  parentheses -- shared by every analysis command
- `prot geometry`: `distance` (atom-atom or centroid-based, covers atom/
  residue/chain pairs), `angle` (3-point, covers bond angles), `dihedral`
  (arbitrary 4-point torsion), `backbone-torsions` (phi/psi/omega + side
  chain chi1-4 for one residue), `rmsd` (Kabsch-fit or raw, between two
  structures in the project), `coords` (centroid, center of mass, bounding
  box, radius of gyration, plane fit, principal axes, moment of inertia for
  a selection), `distmatrix` (pairwise centroid distance matrix across
  named selections)

Remaining spec commands (secondary/contact/pocket/mutate/model/predict/
compare/cluster/annotate/map/plot/view/replay) are not yet implemented.
`geometry`'s convex hull and residue-residue angle matrix were left out of
this pass as under-specified; happy to add them on request.

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
