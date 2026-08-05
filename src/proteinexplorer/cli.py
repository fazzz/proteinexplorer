"""ProteinExplorer CLI skeleton.

Implemented so far: `prot import`, `prot export`, `prot status`.
Remaining spec commands (info/clean/descriptor/geometry/...) will be added
on top of this same project/io layer.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import click

from proteinexplorer import io as pio
from proteinexplorer import project as proj
from proteinexplorer.project import ProjectError, StructureNotFoundError


@click.group()
@click.version_option(package_name="proteinexplorer")
def cli() -> None:
    """prot: CLI-based structural bioinformatics workbench for static
    protein 3D structures."""


@cli.command("import")
@click.argument("path", type=click.Path(exists=True))
@click.option("--name", default=None, help="Display name for the structure (default: filename stem).")
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["pdb", "mmcif"]),
    default=None,
    help="Force input format instead of inferring it from the file extension.",
)
def import_cmd(path: str, name: str | None, fmt: str | None) -> None:
    """Import a PDB/mmCIF structure file into the current project."""
    try:
        record = proj.import_structure(".", path, name=name, fmt=fmt)
        proj.log_command(proj.find_project_root("."), sys.argv[1:])
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Imported '{record.name}' as {record.id}")
    click.echo(
        f"  {record.n_chains} chain(s), {record.n_residues} residue(s), "
        f"{record.n_atoms} atom(s)"
    )
    if record.hetero_resnames:
        click.echo(f"  hetero groups: {', '.join(record.hetero_resnames)}")


@cli.command("export")
@click.argument("structure_id")
@click.argument("dest", type=click.Path())
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["pdb", "mmcif"]),
    default=None,
    help="Force output format instead of inferring it from the destination extension.",
)
def export_cmd(structure_id: str, dest: str, fmt: str | None) -> None:
    """Export a structure from the project to a PDB/mmCIF file."""
    try:
        out_path = proj.export_structure(".", structure_id, dest, fmt=fmt)
    except StructureNotFoundError as exc:
        raise click.ClickException(str(exc)) from exc
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Exported {structure_id} -> {out_path}")


@cli.command("status")
def status_cmd() -> None:
    """List structures currently in the project."""
    try:
        records = proj.list_records(".")
    except ProjectError:
        records = []

    if not records:
        click.echo("No structures imported yet. Use `prot import <file>`.")
        return

    click.echo(f"{len(records)} structure(s):")
    for r in records:
        cats = ", ".join(f"{k}={v}" for k, v in r.category_totals.items())
        click.echo(f"  {r.id}  {r.name:<20}  chains={r.n_chains}  atoms={r.n_atoms}  ({cats})")


@cli.command("info")
@click.argument("structure_id")
def info_cmd(structure_id: str) -> None:
    """Show a detailed summary of one structure in the project."""
    try:
        root = proj.find_project_root(".")
        record = proj.get_record(root, structure_id)
        path = proj.structure_path(root, structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    structure = pio.load_structure(path, structure_id=record.id, fmt=record.format)
    header = pio.header_info(structure)

    click.echo(f"{record.id}  ({record.name})")
    click.echo(f"  source: {record.source_path}")
    click.echo(f"  format: {record.format}   imported: {record.imported_at}")
    if header["structure_method"]:
        click.echo(f"  method: {header['structure_method']}")
    if header["resolution"]:
        click.echo(f"  resolution: {header['resolution']} A")
    click.echo(f"  models: {record.n_models}   chains: {record.n_chains}")
    click.echo(f"  residues: {record.n_residues}   atoms: {record.n_atoms}")
    for category, count in record.category_totals.items():
        click.echo(f"    {category}: {count}")
    if record.hetero_resnames:
        click.echo(f"  hetero groups: {', '.join(record.hetero_resnames)}")
    if record.has_altloc:
        click.echo("  note: contains alternate conformations (altloc)")


@cli.command("descriptor")
@click.argument("structure_id")
def descriptor_cmd(structure_id: str) -> None:
    """Compute structure descriptors (MW, SASA, Rg, contact density, ...)."""
    from proteinexplorer import descriptor as desc

    try:
        root = proj.find_project_root(".")
        record = proj.get_record(root, structure_id)
        path = proj.structure_path(root, structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    structure = pio.load_structure(path, structure_id=record.id, fmt=record.format)
    d = desc.compute_descriptors(structure, path, record.category_totals)

    click.echo(f"{record.id}  ({record.name})")
    click.echo(f"  molecular weight: {d.molecular_weight:.1f} Da (heavy-atom)")
    click.echo(f"  atoms: {d.n_atoms}   residues: {d.n_residues}   chains: {d.n_chains}")
    click.echo(f"  ligands: {d.n_ligands}   waters: {d.n_waters}")
    if d.sasa_total is not None:
        click.echo(f"  SASA: {d.sasa_total:.1f} A^2")
    if d.radius_of_gyration is not None:
        click.echo(f"  radius of gyration: {d.radius_of_gyration:.2f} A")
    if d.contact_density is not None:
        click.echo(f"  contact density (CA-CA < 8A / residue): {d.contact_density:.2f}")
    if d.hydrophobic_ratio is not None:
        click.echo(f"  hydrophobic ratio: {d.hydrophobic_ratio:.2f}")
    click.echo(f"  disulfide bonds: {d.disulfide_count}")
    if d.secondary_structure is not None:
        composition = ", ".join(f"{k}={v:.2f}" for k, v in sorted(d.secondary_structure.items()))
        click.echo(f"  secondary structure ({d.secondary_structure_method}): {composition}")
    elif d.secondary_structure_error:
        click.echo(f"  secondary structure: unavailable ({d.secondary_structure_error})")


def _load(structure_id: str):
    """Resolve a structure_id/name to (record, Bio.PDB Structure) using the
    current project."""
    root = proj.find_project_root(".")
    record = proj.get_record(root, structure_id)
    path = proj.structure_path(root, structure_id)
    structure = pio.load_structure(path, structure_id=record.id, fmt=record.format)
    return record, structure


def _select_or_fail(structure, expr: str):
    from proteinexplorer import selection as sel

    try:
        atoms = sel.select(structure, expr)
    except sel.SelectionSyntaxError as exc:
        raise click.ClickException(f"Invalid selection {expr!r}: {exc}") from exc
    if not atoms:
        raise click.ClickException(f"Selection {expr!r} matched no atoms")
    return atoms


@cli.group("geometry")
def geometry_group() -> None:
    """Distance, angle, dihedral, and coordinate-based geometric analysis."""


@geometry_group.command("distance")
@click.argument("structure_id")
@click.argument("selection_a")
@click.argument("selection_b")
def geometry_distance_cmd(structure_id: str, selection_a: str, selection_b: str) -> None:
    """Distance between two selections (atom-atom, or centroid-based for
    multi-atom selections -- covers atom-residue/residue-residue/chain-chain)."""
    from proteinexplorer import geometry as geom

    try:
        _, structure = _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc
    atoms_a = _select_or_fail(structure, selection_a)
    atoms_b = _select_or_fail(structure, selection_b)
    click.echo(f"{geom.distance(atoms_a, atoms_b):.3f} A")


@geometry_group.command("angle")
@click.argument("structure_id")
@click.argument("selection_a")
@click.argument("selection_b")
@click.argument("selection_c")
def geometry_angle_cmd(structure_id: str, selection_a: str, selection_b: str, selection_c: str) -> None:
    """Angle (degrees) at vertex B for selections A-B-C (bond angle when
    each selection is a single atom; arbitrary 3-point angle otherwise)."""
    from proteinexplorer import geometry as geom

    try:
        _, structure = _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc
    a = _select_or_fail(structure, selection_a)
    b = _select_or_fail(structure, selection_b)
    c = _select_or_fail(structure, selection_c)
    click.echo(f"{geom.angle(a, b, c):.2f} deg")


@geometry_group.command("dihedral")
@click.argument("structure_id")
@click.argument("selection_a")
@click.argument("selection_b")
@click.argument("selection_c")
@click.argument("selection_d")
def geometry_dihedral_cmd(
    structure_id: str, selection_a: str, selection_b: str, selection_c: str, selection_d: str
) -> None:
    """Torsion angle (degrees, -180..180) for selections A-B-C-D."""
    from proteinexplorer import geometry as geom

    try:
        _, structure = _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc
    a = _select_or_fail(structure, selection_a)
    b = _select_or_fail(structure, selection_b)
    c = _select_or_fail(structure, selection_c)
    d = _select_or_fail(structure, selection_d)
    click.echo(f"{geom.dihedral(a, b, c, d):.2f} deg")


@geometry_group.command("backbone-torsions")
@click.argument("structure_id")
@click.option("--chain", "chain_id", required=True, help="Chain ID.")
@click.option("--resid", required=True, type=int, help="Residue sequence number.")
def geometry_backbone_torsions_cmd(structure_id: str, chain_id: str, resid: int) -> None:
    """phi/psi/omega and side-chain chi angles for one residue."""
    from proteinexplorer import geometry as geom

    try:
        _, structure = _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    model = next(iter(structure))
    if chain_id not in model:
        raise click.ClickException(f"No chain {chain_id!r} in this structure")
    try:
        result = geom.backbone_torsions(model[chain_id], resid)
    except geom.GeometryError as exc:
        raise click.ClickException(str(exc)) from exc

    def fmt(value):
        return f"{value:.1f}" if value is not None else "n/a"

    click.echo(f"{result.chain_id}/{result.resname}{result.resseq}")
    click.echo(f"  phi={fmt(result.phi)}  psi={fmt(result.psi)}  omega={fmt(result.omega)}")
    if result.chi:
        chi_str = "  ".join(f"chi{i + 1}={v:.1f}" for i, v in enumerate(result.chi))
        click.echo(f"  {chi_str}")


@geometry_group.command("rmsd")
@click.argument("structure_id_a")
@click.argument("structure_id_b")
@click.option("--selection", "selection", default="protein and atom CA",
              help="Selection applied to both structures (default: protein CA atoms).")
@click.option("--no-fit", is_flag=True, help="Skip superposition; compute raw coordinate RMSD.")
def geometry_rmsd_cmd(structure_id_a: str, structure_id_b: str, selection: str, no_fit: bool) -> None:
    """RMSD between matching selections of two structures in the project."""
    from proteinexplorer import geometry as geom

    try:
        _, structure_a = _load(structure_id_a)
        _, structure_b = _load(structure_id_b)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    atoms_a = _select_or_fail(structure_a, selection)
    atoms_b = _select_or_fail(structure_b, selection)
    try:
        value = geom.rmsd(atoms_a, atoms_b, fit=not no_fit)
    except geom.GeometryError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"RMSD ({'no-fit' if no_fit else 'fit'}, {len(atoms_a)} atoms): {value:.3f} A")


@geometry_group.command("coords")
@click.argument("structure_id")
@click.argument("selection", default="all")
def geometry_coords_cmd(structure_id: str, selection: str) -> None:
    """Coordinate analysis for a selection: centroid, center of mass,
    bounding box, radius of gyration, plane fit, principal axes, and
    moment of inertia."""
    from proteinexplorer import geometry as geom

    try:
        _, structure = _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc
    atoms = _select_or_fail(structure, selection)

    click.echo(f"{len(atoms)} atom(s) selected")
    click.echo(f"  centroid: {geom.centroid(atoms).round(3).tolist()}")
    try:
        click.echo(f"  center of mass: {geom.center_of_mass(atoms).round(3).tolist()}")
    except geom.GeometryError as exc:
        click.echo(f"  center of mass: unavailable ({exc})")

    bbox = geom.bounding_box(atoms)
    click.echo(f"  bounding box: min={bbox.min.round(3).tolist()} max={bbox.max.round(3).tolist()}")

    try:
        rg = geom.radius_of_gyration(atoms)
        click.echo(f"  radius of gyration: {rg:.3f} A")
    except geom.GeometryError as exc:
        click.echo(f"  radius of gyration: unavailable ({exc})")

    try:
        plane = geom.fit_plane(atoms)
        click.echo(
            f"  plane fit: normal={plane.normal.round(3).tolist()} "
            f"rms_deviation={plane.rms_deviation:.3f} A"
        )
    except geom.GeometryError as exc:
        click.echo(f"  plane fit: unavailable ({exc})")

    try:
        axes = geom.principal_axes(atoms)
        click.echo(f"  principal axes eigenvalues: {axes.eigenvalues.round(3).tolist()}")
    except geom.GeometryError as exc:
        click.echo(f"  principal axes: unavailable ({exc})")

    try:
        inertia = geom.moment_of_inertia(atoms)
        click.echo(f"  moment of inertia eigenvalues: {inertia.eigenvalues.round(3).tolist()}")
    except geom.GeometryError as exc:
        click.echo(f"  moment of inertia: unavailable ({exc})")


@geometry_group.command("distmatrix")
@click.argument("structure_id")
@click.argument("selections", nargs=-1, required=True)
def geometry_distmatrix_cmd(structure_id: str, selections: tuple[str, ...]) -> None:
    """Pairwise centroid distance matrix across two or more selections
    (e.g. one --selection-style argument per residue/chain/group)."""
    from proteinexplorer import geometry as geom

    if len(selections) < 2:
        raise click.ClickException("Provide at least 2 selections")

    try:
        _, structure = _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    groups = [_select_or_fail(structure, expr) for expr in selections]
    try:
        matrix = geom.distance_matrix(groups)
    except geom.GeometryError as exc:
        raise click.ClickException(str(exc)) from exc

    header = "        " + "  ".join(f"{i + 1:>8}" for i in range(len(selections)))
    click.echo(header)
    for i, row in enumerate(matrix):
        click.echo(f"{i + 1:>6}  " + "  ".join(f"{v:8.3f}" for v in row))
    for i, expr in enumerate(selections):
        click.echo(f"  [{i + 1}] {expr}")


@cli.group("contact")
def contact_group() -> None:
    """Hydrogen bonds, salt bridges, hydrophobic/pi contacts, disulfides,
    contact maps, and residue interaction networks."""


def _contact_load(structure_id: str):
    try:
        return _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc


@contact_group.command("hbond")
@click.argument("structure_id")
@click.option("--cutoff", default=3.5, show_default=True, help="Heavy-atom N/O...N/O distance cutoff (A).")
def contact_hbond_cmd(structure_id: str, cutoff: float) -> None:
    """Heavy-atom hydrogen bond candidates (no hydrogens required)."""
    from proteinexplorer import contact as ct

    _, structure = _contact_load(structure_id)
    bonds = ct.find_hydrogen_bonds(structure, cutoff=cutoff)
    if not bonds:
        click.echo("No hydrogen bond candidates found.")
        return
    for b in bonds:
        click.echo(f"{b.donor_residue}:{b.donor_atom}  --  {b.acceptor_residue}:{b.acceptor_atom}   {b.distance:.2f} A")


@contact_group.command("saltbridge")
@click.argument("structure_id")
@click.option("--cutoff", default=4.0, show_default=True, help="Charged-group centroid distance cutoff (A).")
def contact_saltbridge_cmd(structure_id: str, cutoff: float) -> None:
    """Salt bridges between basic (Arg/Lys/His) and acidic (Asp/Glu) residues."""
    from proteinexplorer import contact as ct

    _, structure = _contact_load(structure_id)
    bridges = ct.find_salt_bridges(structure, cutoff=cutoff)
    if not bridges:
        click.echo("No salt bridges found.")
        return
    for b in bridges:
        click.echo(f"{b.basic_residue}  --  {b.acidic_residue}   {b.distance:.2f} A")


@contact_group.command("hydrophobic")
@click.argument("structure_id")
@click.option("--cutoff", default=4.5, show_default=True, help="Sidechain carbon-carbon distance cutoff (A).")
def contact_hydrophobic_cmd(structure_id: str, cutoff: float) -> None:
    """Hydrophobic sidechain-sidechain contacts."""
    from proteinexplorer import contact as ct

    _, structure = _contact_load(structure_id)
    contacts = ct.find_hydrophobic_contacts(structure, cutoff=cutoff)
    if not contacts:
        click.echo("No hydrophobic contacts found.")
        return
    for c in contacts:
        click.echo(f"{c.residue_a}  --  {c.residue_b}   {c.min_distance:.2f} A")


@contact_group.command("pipi")
@click.argument("structure_id")
@click.option("--cutoff", default=7.0, show_default=True, help="Aromatic ring centroid distance cutoff (A).")
def contact_pipi_cmd(structure_id: str, cutoff: float) -> None:
    """Pi-pi stacking between aromatic residues (Phe/Tyr/Trp/His)."""
    from proteinexplorer import contact as ct

    _, structure = _contact_load(structure_id)
    interactions = ct.find_pipi_interactions(structure, cutoff=cutoff)
    if not interactions:
        click.echo("No pi-pi interactions found.")
        return
    for p in interactions:
        click.echo(
            f"{p.residue_a}  --  {p.residue_b}   {p.centroid_distance:.2f} A   "
            f"angle={p.plane_angle:.1f} deg   ({p.stack_type})"
        )


@contact_group.command("cationpi")
@click.argument("structure_id")
@click.option("--cutoff", default=6.0, show_default=True, help="Cation-ring-centroid distance cutoff (A).")
def contact_cationpi_cmd(structure_id: str, cutoff: float) -> None:
    """Cation-pi interactions between Arg/Lys and aromatic residues."""
    from proteinexplorer import contact as ct

    _, structure = _contact_load(structure_id)
    interactions = ct.find_cationpi_interactions(structure, cutoff=cutoff)
    if not interactions:
        click.echo("No cation-pi interactions found.")
        return
    for c in interactions:
        click.echo(f"{c.cation_residue}  --  {c.aromatic_residue}   {c.distance:.2f} A")


@contact_group.command("disulfide")
@click.argument("structure_id")
@click.option("--cutoff", default=2.5, show_default=True, help="SG-SG distance cutoff (A).")
def contact_disulfide_cmd(structure_id: str, cutoff: float) -> None:
    """Disulfide (Cys-Cys) bonds."""
    from proteinexplorer import contact as ct

    _, structure = _contact_load(structure_id)
    bonds = ct.find_disulfide_bonds(structure, cutoff=cutoff)
    if not bonds:
        click.echo("No disulfide bonds found.")
        return
    for b in bonds:
        click.echo(f"{b.residue_a}  --  {b.residue_b}   {b.distance:.2f} A")


@contact_group.command("map")
@click.argument("structure_id")
@click.option("--selection", default=None, help="Restrict to a selection (default: everything but water).")
@click.option("--mode", type=click.Choice(["ca", "heavy"]), default="ca", show_default=True,
              help="ca: CA-CA distance. heavy: minimum heavy-atom distance.")
@click.option("--cutoff", default=8.0, show_default=True, help="Contact distance cutoff (A).")
def contact_map_cmd(structure_id: str, selection: str | None, mode: str, cutoff: float) -> None:
    """Residue-residue contact map."""
    from proteinexplorer import contact as ct

    _, structure = _contact_load(structure_id)
    atoms = _select_or_fail(structure, selection) if selection else None
    cm = ct.contact_map(structure, atoms=atoms, cutoff=cutoff, mode=mode)

    n = len(cm.labels)
    click.echo(f"{n} residues, {int(cm.matrix.sum() / 2)} contact(s) (mode={mode}, cutoff={cutoff} A)")
    for i in range(n):
        contacts_i = [cm.labels[j] for j in range(n) if cm.matrix[i, j]]
        if contacts_i:
            click.echo(f"  {cm.labels[i]}: {', '.join(contacts_i)}")


@contact_group.command("network")
@click.argument("structure_id")
def contact_network_cmd(structure_id: str) -> None:
    """All interaction types combined into one residue interaction network edge list."""
    from proteinexplorer import contact as ct

    _, structure = _contact_load(structure_id)
    edges = ct.interaction_network(structure)
    if not edges:
        click.echo("No interactions found.")
        return
    click.echo(f"{len(edges)} edge(s):")
    for e in edges:
        click.echo(f"  {e.residue_a} -[{e.kind}]- {e.residue_b}   {e.value:.2f} A")


@cli.command("secondary")
@click.argument("structure_id")
@click.option("--method", type=click.Choice(["auto", "dssp", "geometric"]), default="auto", show_default=True,
              help="auto: DSSP if available, else the built-in phi/psi classifier.")
@click.option("--chain", "chain_id", default=None, help="Restrict output to one chain.")
def secondary_cmd(structure_id: str, method: str, chain_id: str | None) -> None:
    """Per-residue secondary structure assignment and composition."""
    from proteinexplorer import secondary as sec

    try:
        record, structure = _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    path = proj.structure_path(proj.find_project_root("."), structure_id)
    try:
        residues, used_method = sec.secondary_structure(structure, pdb_path=path, method=method)
    except (sec.DSSPNotAvailableError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc

    if chain_id:
        residues = [r for r in residues if r.chain_id == chain_id]
        if not residues:
            raise click.ClickException(f"No residues found for chain {chain_id!r}")

    click.echo(f"{record.id}  ({record.name})   method={used_method}")

    current_chain = None
    codes: list[str] = []
    for r in residues:
        if r.chain_id != current_chain:
            if codes:
                click.echo(f"  {current_chain}: {''.join(codes)}")
            current_chain = r.chain_id
            codes = []
        codes.append(r.code if r.code != "-" else "C")
    if codes:
        click.echo(f"  {current_chain}: {''.join(codes)}")

    comp = sec.composition(residues)
    comp_str = ", ".join(f"{k}={v:.2f}" for k, v in sorted(comp.items()))
    click.echo(f"  composition: {comp_str}")


@cli.group("pocket")
def pocket_group() -> None:
    """Cavity/binding-pocket detection (grid-based, LIGSITE-style
    approximation -- see `prot pocket detect --help` for caveats)."""


@pocket_group.command("detect")
@click.argument("structure_id")
@click.option("--selection", default=None,
              help="Restrict the search region to a selection (default: the whole structure). "
                   "Strongly recommended for large structures to keep the grid small.")
@click.option("--spacing", default=1.5, show_default=True, help="Grid spacing (A).")
@click.option("--padding", default=3.0, show_default=True, help="Padding around the search region (A).")
@click.option("--max-ray-length", default=6.0, show_default=True, help="Max ray-cast distance per direction (A).")
@click.option("--min-enclosed-axes", default=5, show_default=True, help="Axes (of 7) that must be enclosed.")
@click.option("--min-pocket-points", default=3, show_default=True, help="Minimum grid points to report a pocket.")
def pocket_detect_cmd(
    structure_id: str, selection: str | None, spacing: float, padding: float,
    max_ray_length: float, min_enclosed_axes: int, min_pocket_points: int,
) -> None:
    """Detect cavities via a dependency-free grid/ray-casting method
    (LIGSITE-style approximation -- not fpocket; see module docs).
    Reports volume, a rough surface-area estimate, lining residues,
    hydrophobicity, and a simple (non-ML) druggability heuristic."""
    from proteinexplorer import pocket as pk

    _, structure = _contact_load(structure_id)
    atoms = _select_or_fail(structure, selection) if selection else None

    try:
        pockets = pk.find_pockets(
            structure, atoms=atoms, spacing=spacing, padding=padding,
            max_ray_length=max_ray_length, min_enclosed_axes=min_enclosed_axes,
            min_pocket_points=min_pocket_points,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    if not pockets:
        click.echo("No pockets found with the current parameters.")
        return

    click.echo(f"{len(pockets)} pocket(s) found:")
    for p in pockets:
        click.echo(
            f"  #{p.id}  volume={p.volume:.0f} A^3  surface~={p.surface_area:.0f} A^2  "
            f"hydrophobicity={p.hydrophobicity:.2f}  druggability~={p.druggability_score:.2f}"
        )
        click.echo(f"      centroid={p.centroid.round(2).tolist()}")
        click.echo(f"      residues: {', '.join(p.residues) if p.residues else '(none within lining distance)'}")


@cli.command("mutate")
@click.argument("structure_id")
@click.option("--chain", "chain_id", required=True, help="Chain ID of the residue to mutate.")
@click.option("--resid", required=True, type=int, help="Residue sequence number to mutate.")
@click.option("--to", "target", required=True, help="Target amino acid (1- or 3-letter code).")
@click.option("--method", type=click.Choice(["auto", "scwrl4", "cb_only"]), default="auto", show_default=True,
              help="auto: Scwrl4 if installed, else the built-in backbone+C-beta fallback.")
@click.option("--name", default=None, help="Name for the new mutant structure (default: derived automatically).")
def mutate_cmd(structure_id: str, chain_id: str, resid: int, target: str, method: str, name: str | None) -> None:
    """Point-mutate one residue, saving the result as a new structure in
    the project (the original is left untouched)."""
    from proteinexplorer import mutate as mut

    try:
        record, structure = _load(structure_id)
        root = proj.find_project_root(".")
        src_path = proj.structure_path(root, structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    try:
        result = mut.mutate_residue(structure, src_path, chain_id, resid, target, method=method)
    except (mut.MutationError, mut.Scwrl4NotAvailableError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc

    mutant_name = name or (
        f"{record.name}_{chain_id}{resid}{result.original_resname}{result.new_resname}"
    )
    ext = ".pdb" if record.format == "pdb" else ".cif"
    tmp_dir = Path(tempfile.mkdtemp())
    tmp_path = tmp_dir / f"{mutant_name}{ext}"
    pio.save_structure(structure, tmp_path, fmt=record.format)
    try:
        new_record = proj.import_structure(root, tmp_path, name=mutant_name, fmt=record.format)
    finally:
        tmp_path.unlink(missing_ok=True)
        tmp_dir.rmdir()

    click.echo(f"{result.chain_id}/{result.original_resname}{result.resseq} -> {result.new_resname}")
    click.echo(f"  method: {result.method}")
    click.echo(f"  atoms placed: {', '.join(result.atoms_placed)}")
    click.echo(f"  {result.note}")
    click.echo(f"Saved as '{new_record.name}' ({new_record.id})")


@cli.group("model")
def model_group() -> None:
    """Missing-residue detection, a crude loop filler, and homology
    modeling (external MODELLER wrapper)."""


@model_group.command("gaps")
@click.argument("structure_id")
def model_gaps_cmd(structure_id: str) -> None:
    """Detect numbering gaps (missing residues) per chain."""
    from proteinexplorer import model as mdl

    _, structure = _contact_load(structure_id)
    gaps = mdl.find_gaps(structure)
    if not gaps:
        click.echo("No gaps detected.")
        return
    for g in gaps:
        click.echo(f"{g.chain_id}: {g.prev_resseq} .. {g.next_resseq}  ({g.length} residue(s) missing)")


@model_group.command("loop")
@click.argument("structure_id")
@click.option("--chain", "chain_id", required=True, help="Chain ID.")
@click.option("--start", required=True, type=int, help="First missing residue number.")
@click.option("--end", required=True, type=int, help="Last missing residue number.")
@click.option("--sequence", default=None,
              help="1-letter sequence for the gap (default: poly-ALA placeholder).")
@click.option("--name", default=None, help="Name for the new filled structure.")
def model_loop_cmd(
    structure_id: str, chain_id: str, start: int, end: int, sequence: str | None, name: str | None
) -> None:
    """Fill a gap with a crude placeholder backbone trace (dependency-free;
    NOT a real loop model -- see `prot model loop --help` caveats in the
    docs). Saves the result as a new structure in the project."""
    from proteinexplorer import model as mdl

    try:
        record, structure = _load(structure_id)
        root = proj.find_project_root(".")
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    try:
        result = mdl.fill_loop_linear(structure, chain_id, start, end, sequence=sequence)
    except mdl.MutationError as exc:
        raise click.ClickException(str(exc)) from exc

    mutant_name = name or f"{record.name}_loop{chain_id}{start}-{end}"
    ext = ".pdb" if record.format == "pdb" else ".cif"
    tmp_dir = Path(tempfile.mkdtemp())
    tmp_path = tmp_dir / f"{mutant_name}{ext}"
    pio.save_structure(structure, tmp_path, fmt=record.format)
    try:
        new_record = proj.import_structure(root, tmp_path, name=mutant_name, fmt=record.format)
    finally:
        tmp_path.unlink(missing_ok=True)
        tmp_dir.rmdir()

    click.echo(f"Filled {result.chain_id}/{result.start_resseq}-{result.end_resseq}: {', '.join(result.residues_added)}")
    click.echo(f"  {result.note}")
    click.echo(f"Saved as '{new_record.name}' ({new_record.id})")


@model_group.command("homology")
@click.option("--alignment", "alignment_path", required=True, type=click.Path(exists=True),
              help="PIR-format alignment file containing target and template sequences.")
@click.option("--template", "template_codes", required=True, multiple=True,
              help="Template code(s) as named in the alignment (repeatable).")
@click.option("--target", "target_code", required=True, help="Target sequence code as named in the alignment.")
@click.option("--template-dir", required=True, type=click.Path(exists=True),
              help="Directory containing the template structure file(s).")
@click.option("--output-dir", required=True, type=click.Path(), help="Directory to write the model(s) to.")
@click.option("--n-models", default=1, show_default=True, help="Number of models to generate.")
def model_homology_cmd(
    alignment_path: str, template_codes: tuple[str, ...], target_code: str,
    template_dir: str, output_dir: str, n_models: int,
) -> None:
    """Homology modeling via an external MODELLER installation (license
    required from https://salilab.org/modeller/). No dependency-free
    fallback exists for this command."""
    from proteinexplorer import model as mdl

    try:
        outputs = mdl.homology_model(
            alignment_pir_path=alignment_path,
            template_codes=list(template_codes),
            target_code=target_code,
            template_search_dir=template_dir,
            output_dir=output_dir,
            n_models=n_models,
        )
    except mdl.ModellerNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"{len(outputs)} model(s) written:")
    for path in outputs:
        click.echo(f"  {path}")


def main() -> None:
    cli()


if __name__ == "__main__":
    main()
