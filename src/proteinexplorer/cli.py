"""ProteinExplorer CLI skeleton.

Implemented so far: `prot import`, `prot export`, `prot status`.
Remaining spec commands (info/clean/descriptor/geometry/...) will be added
on top of this same project/io layer.
"""

from __future__ import annotations

import sys

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
        click.echo(f"  secondary structure: {composition}")
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


def main() -> None:
    cli()


if __name__ == "__main__":
    main()
