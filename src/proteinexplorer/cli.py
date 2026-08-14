"""ProteinExplorer CLI skeleton.

Implemented so far: `prot import`, `prot export`, `prot status`.
Remaining spec commands (info/clean/descriptor/geometry/...) will be added
on top of this same project/io layer.
"""

from __future__ import annotations

import sys
import tempfile
import shutil
from pathlib import Path

import click

from proteinexplorer import io as pio
from proteinexplorer import project as proj
from proteinexplorer.project import ProjectError, StructureNotFoundError

# The actual argv a command was invoked with. sys.argv[1:] only reflects
# the true invocation for the real `prot ...` entry point -- when invoked
# in-process via Click's CliRunner (as every CLI test, and `prot replay`
# itself, do), sys.argv still belongs to the outer process. ProtGroup.main
# below captures the real args Click received either way, so log_command()
# calls read from here instead of sys.argv directly.
_current_argv: list[str] = []


def current_argv() -> list[str]:
    return list(_current_argv)


class ProtGroup(click.Group):
    def main(self, args=None, **kwargs):
        global _current_argv
        _current_argv = list(args) if args is not None else list(sys.argv[1:])
        return super().main(args=args, **kwargs)


@click.group(cls=ProtGroup)
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
        proj.log_command(proj.find_project_root("."), current_argv())
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
    proj.log_command(root, current_argv())

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
    proj.log_command(root, current_argv())

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


@cli.group("compare")
def compare_group() -> None:
    """Compare two structures already in the project: RMSD, TM-score,
    secondary structure similarity, contact similarity, pocket overlap,
    and ligand comparison."""


def _compare_load(structure_id_a: str, structure_id_b: str):
    try:
        root = proj.find_project_root(".")
        record_a = proj.get_record(root, structure_id_a)
        record_b = proj.get_record(root, structure_id_b)
        path_a = proj.structure_path(root, structure_id_a)
        path_b = proj.structure_path(root, structure_id_b)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc
    structure_a = pio.load_structure(path_a, structure_id=record_a.id, fmt=record_a.format)
    structure_b = pio.load_structure(path_b, structure_id=record_b.id, fmt=record_b.format)
    return (record_a, structure_a, path_a), (record_b, structure_b, path_b)


@compare_group.command("rmsd")
@click.argument("structure_id_a")
@click.argument("structure_id_b")
@click.option("--no-fit", is_flag=True, help="Skip superposition; compute raw coordinate RMSD.")
def compare_rmsd_cmd(structure_id_a: str, structure_id_b: str, no_fit: bool) -> None:
    """RMSD over CA atoms common to both structures (matched by chain ID + residue number)."""
    from proteinexplorer import compare as cmp

    (_, structure_a, _), (_, structure_b, _) = _compare_load(structure_id_a, structure_id_b)
    try:
        value, n = cmp.rmsd(structure_a, structure_b, fit=not no_fit)
    except cmp.CompareError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"RMSD ({'no-fit' if no_fit else 'fit'}, {n} common CA atoms): {value:.3f} A")


@compare_group.command("tmscore")
@click.argument("structure_id_a")
@click.argument("structure_id_b")
@click.option("--method", type=click.Choice(["auto", "tmalign", "fallback"]), default="auto", show_default=True)
def compare_tmscore_cmd(structure_id_a: str, structure_id_b: str, method: str) -> None:
    """TM-score. Uses external TMalign/US-align if installed (real
    structural alignment); otherwise a fixed-correspondence fallback that
    is NOT numerically comparable to real TM-align output (see notes)."""
    from proteinexplorer import compare as cmp

    (_, structure_a, path_a), (_, structure_b, path_b) = _compare_load(structure_id_a, structure_id_b)
    try:
        result = cmp.tm_score(structure_a, structure_b, structure_a_path=path_a, structure_b_path=path_b, method=method)
    except cmp.CompareError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"TM-score: {result.score:.3f}  (method={result.method}, n_residues={result.n_residues})")
    click.echo(f"  {result.note}")


@compare_group.command("secondary")
@click.argument("structure_id_a")
@click.argument("structure_id_b")
def compare_secondary_cmd(structure_id_a: str, structure_id_b: str) -> None:
    """Secondary structure similarity (Q3-style, collapsed to H/E/C) over
    residues common to both structures."""
    from proteinexplorer import compare as cmp

    (_, structure_a, _), (_, structure_b, _) = _compare_load(structure_id_a, structure_id_b)
    try:
        score, n = cmp.secondary_structure_similarity(structure_a, structure_b)
    except cmp.CompareError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Secondary structure similarity: {score:.2f}  ({n} common residues)")


@compare_group.command("contact")
@click.argument("structure_id_a")
@click.argument("structure_id_b")
@click.option("--mode", type=click.Choice(["ca", "heavy"]), default="ca", show_default=True)
@click.option("--cutoff", default=8.0, show_default=True)
def compare_contact_cmd(structure_id_a: str, structure_id_b: str, mode: str, cutoff: float) -> None:
    """Jaccard similarity between the two structures' contact maps
    (restricted to residues present in both)."""
    from proteinexplorer import compare as cmp

    (_, structure_a, _), (_, structure_b, _) = _compare_load(structure_id_a, structure_id_b)
    try:
        jaccard, shared, union = cmp.contact_similarity(structure_a, structure_b, mode=mode, cutoff=cutoff)
    except cmp.CompareError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Contact similarity (Jaccard): {jaccard:.2f}  ({shared}/{union} contacts shared)")


@compare_group.command("pocket")
@click.argument("structure_id_a")
@click.argument("structure_id_b")
@click.option("--pocket-a", default=1, show_default=True, help="Pocket # in structure A (1 = largest).")
@click.option("--pocket-b", default=1, show_default=True, help="Pocket # in structure B (1 = largest).")
@click.option("--selection", default=None, help="Restrict the pocket search region in both structures.")
@click.option("--spacing", default=1.5, show_default=True, help="Grid spacing (A), passed to pocket detection.")
@click.option("--padding", default=3.0, show_default=True, help="Padding around the search region (A).")
@click.option("--min-pocket-points", default=3, show_default=True, help="Minimum grid points to count as a pocket.")
def compare_pocket_cmd(
    structure_id_a: str, structure_id_b: str, pocket_a: int, pocket_b: int,
    selection: str | None, spacing: float, padding: float, min_pocket_points: int,
) -> None:
    """Jaccard similarity between two pockets' lining residues."""
    from proteinexplorer import compare as cmp

    (_, structure_a, _), (_, structure_b, _) = _compare_load(structure_id_a, structure_id_b)
    atoms_a = _select_or_fail(structure_a, selection) if selection else None
    atoms_b = _select_or_fail(structure_b, selection) if selection else None
    try:
        jaccard, res_a, res_b = cmp.pocket_overlap(
            structure_a, structure_b, pocket_index_a=pocket_a, pocket_index_b=pocket_b,
            atoms_a=atoms_a, atoms_b=atoms_b, spacing=spacing, padding=padding,
            min_pocket_points=min_pocket_points,
        )
    except cmp.CompareError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Pocket overlap (Jaccard): {jaccard:.2f}")
    click.echo(f"  A pocket #{pocket_a}: {', '.join(res_a)}")
    click.echo(f"  B pocket #{pocket_b}: {', '.join(res_b)}")


@compare_group.command("ligand")
@click.argument("structure_id_a")
@click.argument("structure_id_b")
@click.option("--no-fit", is_flag=True, help="Skip superposition for ligand RMSD; use raw coordinates.")
def compare_ligand_cmd(structure_id_a: str, structure_id_b: str, no_fit: bool) -> None:
    """Compare bound ligands: which resnames are shared, and RMSD for any
    shared ligand present as a single matching-atom instance in each."""
    from proteinexplorer import compare as cmp

    (_, structure_a, _), (_, structure_b, _) = _compare_load(structure_id_a, structure_id_b)
    result = cmp.ligand_comparison(structure_a, structure_b, fit=not no_fit)

    click.echo(f"Common ligands: {', '.join(result.common_resnames) or '(none)'}")
    if result.only_in_a:
        click.echo(f"  only in A: {', '.join(result.only_in_a)}")
    if result.only_in_b:
        click.echo(f"  only in B: {', '.join(result.only_in_b)}")
    for name, value in result.rmsd_by_resname.items():
        click.echo(f"  {name} RMSD: {value:.3f} A")


@cli.group("cluster")
def cluster_group() -> None:
    """Ensemble clustering by pairwise RMSD -- either several structures
    already in the project, or several MODEL records within one
    multi-model file (e.g. an NMR ensemble)."""


def _cluster_report(result, matrix_note: str = "") -> None:
    click.echo(f"{len(result.clusters)} cluster(s) (method={result.method})")
    for c in result.clusters:
        others = [m for m in c.member_labels if m != c.representative_label]
        others_str = f" + {', '.join(others)}" if others else ""
        click.echo(f"  #{c.id}: representative={c.representative_label}{others_str}")


@cluster_group.command("ensemble")
@click.argument("structure_ids", nargs=-1, required=True)
@click.option("--selection", default="protein and atom CA", show_default=True,
              help="Selection applied to every structure (must match atom count across all of them).")
@click.option("--method", type=click.Choice(["greedy", "hierarchical"]), default="greedy", show_default=True)
@click.option("--threshold", default=2.0, show_default=True, help="RMSD threshold (A) for --method greedy.")
@click.option("--n-clusters", default=None, type=int, help="Target cluster count for --method hierarchical.")
@click.option("--distance-threshold", default=None, type=float,
              help="RMSD cutoff (A) for --method hierarchical (alternative to --n-clusters).")
@click.option("--no-fit", is_flag=True, help="Skip superposition; use raw coordinate RMSD.")
def cluster_ensemble_cmd(
    structure_ids: tuple[str, ...], selection: str, method: str, threshold: float,
    n_clusters: int | None, distance_threshold: float | None, no_fit: bool,
) -> None:
    """Cluster several structures already in the project by pairwise RMSD."""
    from proteinexplorer import cluster as clu

    if len(structure_ids) < 2:
        raise click.ClickException("Provide at least 2 structure IDs")

    atom_groups = []
    for sid in structure_ids:
        _, structure = _contact_load(sid)
        atom_groups.append(_select_or_fail(structure, selection))

    try:
        matrix = clu.pairwise_rmsd_matrix(atom_groups, fit=not no_fit)
        if method == "greedy":
            result = clu.greedy(list(structure_ids), matrix, threshold=threshold)
        else:
            result = clu.hierarchical(
                list(structure_ids), matrix, n_clusters=n_clusters, distance_threshold=distance_threshold
            )
    except (clu.ClusterError, clu.ClusterExtraNotAvailableError) as exc:
        raise click.ClickException(str(exc)) from exc

    _cluster_report(result)


@cluster_group.command("models")
@click.argument("structure_id")
@click.option("--selection", default="protein and atom CA", show_default=True)
@click.option("--method", type=click.Choice(["greedy", "hierarchical"]), default="greedy", show_default=True)
@click.option("--threshold", default=2.0, show_default=True, help="RMSD threshold (A) for --method greedy.")
@click.option("--n-clusters", default=None, type=int, help="Target cluster count for --method hierarchical.")
@click.option("--distance-threshold", default=None, type=float)
@click.option("--no-fit", is_flag=True)
def cluster_models_cmd(
    structure_id: str, selection: str, method: str, threshold: float,
    n_clusters: int | None, distance_threshold: float | None, no_fit: bool,
) -> None:
    """Cluster the MODEL records within one multi-model structure (e.g. an
    NMR ensemble) by pairwise RMSD."""
    from proteinexplorer import cluster as clu
    from proteinexplorer import selection as sel
    import copy

    try:
        record, full_structure = _load(structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    n_models = len(full_structure)
    if n_models < 2:
        raise click.ClickException(
            f"{structure_id} has only 1 model. Use `prot cluster ensemble` to "
            f"compare multiple imported structures instead."
        )

    atom_groups = []
    labels = []
    for model in full_structure:
        model_copy = copy.deepcopy(model)
        model_copy.detach_parent()
        single = full_structure.__class__(full_structure.id)
        single.add(model_copy)
        atoms = sel.select(single, selection)
        if not atoms:
            raise click.ClickException(f"Selection {selection!r} matched no atoms in model {model.id}")
        atom_groups.append(atoms)
        labels.append(f"model_{model.id}")

    try:
        matrix = clu.pairwise_rmsd_matrix(atom_groups, fit=not no_fit)
        if method == "greedy":
            result = clu.greedy(labels, matrix, threshold=threshold)
        else:
            result = clu.hierarchical(labels, matrix, n_clusters=n_clusters, distance_threshold=distance_threshold)
    except (clu.ClusterError, clu.ClusterExtraNotAvailableError) as exc:
        raise click.ClickException(str(exc)) from exc

    _cluster_report(result)


@cli.group("plot")
def plot_group() -> None:
    """Static plots (matplotlib, needs `pip install -e ".[viz]"`):
    Ramachandran, contact map, secondary structure diagram."""


@plot_group.command("ramachandran")
@click.argument("structure_id")
@click.argument("output", type=click.Path())
@click.option("--chain", "chain_id", default=None, help="Restrict to one chain.")
def plot_ramachandran_cmd(structure_id: str, output: str, chain_id: str | None) -> None:
    """Phi/psi scatter plot with the geometric classifier's alpha/beta
    regions shaded for reference."""
    from proteinexplorer import plot as plt_mod

    _, structure = _contact_load(structure_id)
    try:
        path = plt_mod.ramachandran_plot(structure, output, chain_id=chain_id)
    except plt_mod.PlotExtraNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Saved {path}")


@plot_group.command("contact-map")
@click.argument("structure_id")
@click.argument("output", type=click.Path())
@click.option("--selection", default=None, help="Restrict to a selection (default: everything but water).")
@click.option("--mode", type=click.Choice(["ca", "heavy"]), default="ca", show_default=True)
@click.option("--cutoff", default=8.0, show_default=True)
def plot_contact_map_cmd(structure_id: str, output: str, selection: str | None, mode: str, cutoff: float) -> None:
    """Contact map heatmap."""
    from proteinexplorer import plot as plt_mod

    _, structure = _contact_load(structure_id)
    atoms = _select_or_fail(structure, selection) if selection else None
    try:
        path = plt_mod.contact_map_plot(structure, output, atoms=atoms, mode=mode, cutoff=cutoff)
    except plt_mod.PlotExtraNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Saved {path}")


@plot_group.command("secondary")
@click.argument("structure_id")
@click.argument("output", type=click.Path())
@click.option("--method", type=click.Choice(["auto", "dssp", "geometric"]), default="auto", show_default=True)
def plot_secondary_cmd(structure_id: str, output: str, method: str) -> None:
    """Linear secondary structure diagram, one track per chain."""
    from proteinexplorer import plot as plt_mod

    _, structure = _contact_load(structure_id)
    try:
        path = plt_mod.secondary_structure_plot(structure, output, method=method)
    except plt_mod.PlotExtraNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Saved {path}")


@cli.group("predict")
def predict_group() -> None:
    """Structure prediction via external tools only (ColabFold/AlphaFold)
    -- no dependency-free fallback exists for this command."""


@predict_group.command("colabfold")
@click.argument("sequence")
@click.option("--name", default="query", show_default=True, help="Name for the FASTA entry and output files.")
@click.option("--output-dir", required=True, type=click.Path(), help="Directory for ColabFold's output.")
@click.option("--import-name", default=None, help="If set, import the top model into the project under this name.")
def predict_colabfold_cmd(sequence: str, name: str, output_dir: str, import_name: str | None) -> None:
    """Predict a structure from a sequence via a local ColabFold
    installation (`colabfold_batch`)."""
    from proteinexplorer import predict as pred

    try:
        result = pred.colabfold_predict(sequence, output_dir, name=name)
    except pred.PredictionToolNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc
    except pred.PredictionError as exc:
        raise click.ClickException(str(exc)) from exc

    if not result.models:
        click.echo(f"ColabFold ran but no output PDB was found under {result.output_dir}")
        return

    click.echo(f"{len(result.models)} model(s) written to {result.output_dir}:")
    for m in result.models:
        click.echo(f"  {m}")

    if import_name:
        try:
            root = proj.find_project_root(".")
        except ProjectError as exc:
            raise click.ClickException(str(exc)) from exc
        record = proj.import_structure(root, result.models[0], name=import_name, fmt="pdb")
        proj.log_command(root, current_argv())
        click.echo(f"Imported top model as '{record.name}' ({record.id})")


@predict_group.command("alphafold")
@click.option("--fasta", "fasta_path", required=True, type=click.Path(exists=True))
@click.option("--output-dir", required=True, type=click.Path())
@click.option("--alphafold-script", required=True, type=click.Path(), help="Path to run_alphafold.sh.")
@click.option("--data-dir", required=True, type=click.Path(), help="AlphaFold parameter/sequence database directory.")
@click.option("--max-template-date", required=True, help="YYYY-MM-DD cutoff for template search.")
@click.option("--model-preset", default="monomer", show_default=True)
@click.option("--db-preset", default="full_dbs", show_default=True)
@click.option("--import-name", default=None, help="If set, import the top model into the project under this name.")
def predict_alphafold_cmd(
    fasta_path: str, output_dir: str, alphafold_script: str, data_dir: str,
    max_template_date: str, model_preset: str, db_preset: str, import_name: str | None,
) -> None:
    """Predict a structure via a local AlphaFold installation
    (`run_alphafold.sh`). Requires the AlphaFold parameter/sequence
    databases to already be set up locally."""
    from proteinexplorer import predict as pred

    try:
        result = pred.alphafold_predict(
            fasta_path=fasta_path, output_dir=output_dir, alphafold_script=alphafold_script,
            data_dir=data_dir, max_template_date=max_template_date,
            model_preset=model_preset, db_preset=db_preset,
        )
    except pred.PredictionToolNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc
    except pred.PredictionError as exc:
        raise click.ClickException(str(exc)) from exc

    if not result.models:
        click.echo(f"AlphaFold ran but no ranked_*.pdb output was found under {result.output_dir}")
        return

    click.echo(f"{len(result.models)} model(s) written to {result.output_dir}:")
    for m in result.models:
        click.echo(f"  {m}")

    if import_name:
        try:
            root = proj.find_project_root(".")
        except ProjectError as exc:
            raise click.ClickException(str(exc)) from exc
        record = proj.import_structure(root, result.models[0], name=import_name, fmt="pdb")
        proj.log_command(root, current_argv())
        click.echo(f"Imported top model as '{record.name}' ({record.id})")


@cli.group("annotate")
def annotate_group() -> None:
    """Structure annotation: built-in metal-binding site detection and
    experimental metadata, plus external UniProt/Pfam lookups."""


@annotate_group.command("metal-sites")
@click.argument("structure_id")
@click.option("--cutoff", default=3.0, show_default=True, help="Ion-to-coordinating-atom distance cutoff (A).")
def annotate_metal_sites_cmd(structure_id: str, cutoff: float) -> None:
    """Detect metal ions and their protein-atom coordinating residues
    (purely geometric, no external database)."""
    from proteinexplorer import annotate as ann

    _, structure = _contact_load(structure_id)
    sites = ann.metal_binding_sites(structure, cutoff=cutoff)
    if not sites:
        click.echo("No metal-binding sites found.")
        return
    for site in sites:
        pairs = ", ".join(f"{r} ({d:.2f} A)" for r, d in zip(site.coordinating_residues, site.coordinating_distances))
        click.echo(f"{site.ion_label} ({site.ion_element}): {pairs}")


@annotate_group.command("metadata")
@click.argument("structure_id")
def annotate_metadata_cmd(structure_id: str) -> None:
    """Experimental metadata (method, resolution, deposition date) from
    the structure's own header."""
    from proteinexplorer import annotate as ann

    _, structure = _contact_load(structure_id)
    meta = ann.structure_metadata(structure)
    click.echo(f"method: {meta.method or 'unknown'}")
    click.echo(f"resolution: {meta.resolution if meta.resolution is not None else 'unknown'}")
    click.echo(f"deposition date: {meta.deposition_date or 'unknown'}")


@annotate_group.command("uniprot")
@click.argument("accession")
def annotate_uniprot_cmd(accession: str) -> None:
    """Gene name, organism, taxonomy, EC numbers, and GO terms for a
    UniProt accession (external REST lookup)."""
    from proteinexplorer import annotate as ann

    try:
        result = ann.uniprot_lookup(accession)
    except ann.AnnotationFetchError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"{result.accession}")
    click.echo(f"  gene(s): {', '.join(result.gene_names) or 'unknown'}")
    click.echo(f"  organism: {result.organism or 'unknown'} (taxon {result.taxonomy_id or 'unknown'})")
    click.echo(f"  EC number(s): {', '.join(result.ec_numbers) or 'none'}")
    click.echo(f"  GO terms: {', '.join(result.go_terms) or 'none'}")


@annotate_group.command("pfam")
@click.argument("accession")
def annotate_pfam_cmd(accession: str) -> None:
    """Pfam domain hits for a UniProt accession (external REST lookup
    via InterPro)."""
    from proteinexplorer import annotate as ann

    try:
        domains = ann.pfam_domains(accession)
    except ann.AnnotationFetchError as exc:
        raise click.ClickException(str(exc)) from exc

    if not domains:
        click.echo("No Pfam domains found.")
        return
    for d in domains:
        click.echo(f"  {d.accession}: {d.name}")


@cli.group("map")
def map_group() -> None:
    """Generate coloring/selection scripts for external viewers (PyMOL,
    ChimeraX, VMD) from pockets, mutation sites, domain ranges, or
    per-residue values (conservation, etc.)."""


@map_group.command("pocket")
@click.argument("structure_id")
@click.argument("output", type=click.Path())
@click.option("--tool", type=click.Choice(["pymol", "chimerax", "vmd"]), default="pymol", show_default=True)
@click.option("--selection", default=None, help="Restrict the pocket search region.")
@click.option("--spacing", default=1.5, show_default=True)
@click.option("--padding", default=3.0, show_default=True)
@click.option("--min-pocket-points", default=3, show_default=True)
def map_pocket_cmd(
    structure_id: str, output: str, tool: str, selection: str | None,
    spacing: float, padding: float, min_pocket_points: int,
) -> None:
    """One color per detected pocket."""
    from proteinexplorer import map as map_mod
    from proteinexplorer import pocket as pk

    record, structure = _contact_load(structure_id)
    atoms = _select_or_fail(structure, selection) if selection else None
    pockets = pk.find_pockets(
        structure, atoms=atoms, spacing=spacing, padding=padding, min_pocket_points=min_pocket_points
    )
    try:
        script = map_mod.pocket_map_script(pockets, tool=tool, object_name=record.name)
    except map_mod.MapError as exc:
        raise click.ClickException(str(exc)) from exc
    Path(output).write_text(script)
    click.echo(f"Saved {output} ({len(pockets)} pocket(s), tool={tool})")


@map_group.command("mutation")
@click.argument("structure_id")
@click.argument("output", type=click.Path())
@click.option("--residue", "residues", multiple=True, required=True,
              help="Residue label chain/resnum (e.g. A/50), repeatable.")
@click.option("--tool", type=click.Choice(["pymol", "chimerax", "vmd"]), default="pymol", show_default=True)
@click.option("--color", default="red", show_default=True)
def map_mutation_cmd(structure_id: str, output: str, residues: tuple[str, ...], tool: str, color: str) -> None:
    """Highlight a set of residues (e.g. mutation sites) in one color."""
    from proteinexplorer import map as map_mod

    record, _ = _contact_load(structure_id)
    try:
        script = map_mod.mutation_map_script(list(residues), tool=tool, object_name=record.name, color=color)
    except map_mod.MapError as exc:
        raise click.ClickException(str(exc)) from exc
    Path(output).write_text(script)
    click.echo(f"Saved {output} ({len(residues)} residue(s), tool={tool})")


@map_group.command("domain")
@click.argument("structure_id")
@click.argument("output", type=click.Path())
@click.option("--range", "ranges", multiple=True, required=True,
              help="chain:start-end:label, e.g. A:10-45:PF00062, repeatable.")
@click.option("--tool", type=click.Choice(["pymol", "chimerax", "vmd"]), default="pymol", show_default=True)
def map_domain_cmd(structure_id: str, output: str, ranges: tuple[str, ...], tool: str) -> None:
    """One color per named residue range (e.g. domain boundaries)."""
    from proteinexplorer import map as map_mod

    record, _ = _contact_load(structure_id)
    domains = []
    for spec in ranges:
        try:
            chain_part, rest = spec.split(":", 1)
            span, label = rest.split(":", 1)
            start_s, end_s = span.split("-", 1)
            domains.append(map_mod.DomainRange(label=label, chain_id=chain_part, start=int(start_s), end=int(end_s)))
        except ValueError as exc:
            raise click.ClickException(
                f"Invalid --range {spec!r}, expected chain:start-end:label (e.g. A:10-45:PF00062)"
            ) from exc

    try:
        script = map_mod.domain_map_script(domains, tool=tool, object_name=record.name)
    except map_mod.MapError as exc:
        raise click.ClickException(str(exc)) from exc
    Path(output).write_text(script)
    click.echo(f"Saved {output} ({len(domains)} domain(s), tool={tool})")


@map_group.command("conservation")
@click.argument("structure_id")
@click.argument("values_csv", type=click.Path(exists=True))
@click.argument("output_script", type=click.Path())
@click.option("--tool", type=click.Choice(["pymol", "chimerax", "vmd"]), default="pymol", show_default=True)
@click.option("--structure-name", default=None, help="Name for the B-factor-annotated structure saved to the project.")
def map_conservation_cmd(
    structure_id: str, values_csv: str, output_script: str, tool: str, structure_name: str | None
) -> None:
    """Write per-residue values (e.g. conservation scores) into a copy of
    the structure's B-factor column, plus a spectrum/color script.
    `values_csv` has one "chain/resnum,value" pair per line, e.g. `A/41,0.87`."""
    from proteinexplorer import map as map_mod

    try:
        record, structure = _load(structure_id)
        root = proj.find_project_root(".")
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    values: dict[str, float] = {}
    for line in Path(values_csv).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        label, value_s = line.rsplit(",", 1)
        values[label.strip()] = float(value_s)

    map_mod.write_bfactors(structure, values)

    out_name = structure_name or f"{record.name}_bfactor"
    ext = ".pdb" if record.format == "pdb" else ".cif"
    tmp_dir = Path(tempfile.mkdtemp())
    tmp_path = tmp_dir / f"{out_name}{ext}"
    pio.save_structure(structure, tmp_path, fmt=record.format)
    try:
        new_record = proj.import_structure(root, tmp_path, name=out_name, fmt=record.format)
    finally:
        tmp_path.unlink(missing_ok=True)
        tmp_dir.rmdir()
    proj.log_command(root, current_argv())

    Path(output_script).write_text(map_mod.spectrum_script(tool=tool, object_name=out_name))
    click.echo(f"Saved B-factor-annotated structure as '{new_record.name}' ({new_record.id})")
    click.echo(f"Saved {output_script} (tool={tool})")


@cli.command("view")
@click.argument("structure_id")
@click.option("--tool", type=click.Choice(["pymol", "chimerax", "vmd"]), default="pymol", show_default=True)
@click.option("--script", "script_path", default=None, type=click.Path(exists=True),
              help="Script to run in the viewer after loading (e.g. from `prot map`).")
def view_cmd(structure_id: str, tool: str, script_path: str | None) -> None:
    """Launch an external 3D viewer (PyMOL/ChimeraX/VMD) on a structure.
    Starts the viewer as a background process and returns immediately --
    there is no dependency-free substitute for this command."""
    from proteinexplorer import view as view_mod

    try:
        root = proj.find_project_root(".")
        path = proj.structure_path(root, structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    try:
        process = view_mod.launch(tool, path, script_path=script_path)
    except (view_mod.ViewerNotAvailableError, view_mod.ViewError) as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Launched {tool} (pid {process.pid}) on {path}" + (f" with {script_path}" if script_path else ""))


@cli.command("replay")
@click.option("--from", "start", default=1, show_default=True, type=int, help="First log entry to replay (1-based).")
@click.option("--to", "end", default=None, type=int, help="Last log entry to replay (default: the end of the log).")
@click.option("--skip", "skip_csv", default=None,
              help="Comma-separated command names to skip (default: view,predict,annotate).")
@click.option("--continue-on-error", is_flag=True, help="Keep going after a step fails instead of stopping.")
@click.option("--no-reset", is_flag=True, help="Replay onto the current project state instead of backing it up and resetting first.")
@click.option("--dry-run", is_flag=True, help="Show what would be replayed without running anything.")
def replay_cmd(start: int, end: int | None, skip_csv: str | None, continue_on_error: bool, no_reset: bool, dry_run: bool) -> None:
    """Re-run the commands recorded in .proteinexplorer/log.json.

    By default, backs up the current project to
    .proteinexplorer_prereplay_<timestamp>, resets it, and replays every
    logged command from scratch, in-process. Structure IDs are
    regenerated on each import; any later command's argv that literally
    referenced an old ID is rewritten to the new one automatically.
    """
    from proteinexplorer import replay as replay_mod

    skip = set(skip_csv.split(",")) if skip_csv else None
    try:
        result = replay_mod.replay(
            ".", start=start, end=end, skip=skip, continue_on_error=continue_on_error,
            reset=not no_reset, dry_run=dry_run,
        )
    except (replay_mod.ReplayError, ProjectError) as exc:
        raise click.ClickException(str(exc)) from exc

    if result.backup_dir is not None:
        click.echo(f"Backed up previous project state to {result.backup_dir}")

    for step in result.steps:
        if step.skipped:
            click.echo(f"  [{step.index}] skip: {' '.join(step.argv)}")
            continue
        if dry_run:
            marker = "plan"
        else:
            marker = "ok" if step.exit_code in (0, None) else f"FAILED (exit {step.exit_code})"
        rewritten_note = "" if step.rewritten_argv == step.argv else f"  (rewritten: {' '.join(step.rewritten_argv)})"
        click.echo(f"  [{step.index}] {marker}: {' '.join(step.argv)}{rewritten_note}")

    if not dry_run:
        click.echo(f"{len(result.steps)} step(s), {result.n_failed} failed")


@cli.group("search")
def search_group() -> None:
    """Structural similarity search via Foldseek. External-tool-only --
    there is no dependency-free fallback for large-scale structural
    database search."""


@search_group.command("foldseek")
@click.argument("structure_id")
@click.option("--target-db", default=None, help="A pre-built Foldseek database path.")
@click.option("--target-dir", default=None, type=click.Path(exists=True),
              help="A directory of structure files to search against (Foldseek builds a temporary index).")
@click.option("--against-project", is_flag=True,
              help="Search against every other structure currently in this project.")
@click.option("--sensitivity", default=None, type=float, help="Foldseek -s sensitivity (higher = slower, more sensitive).")
@click.option("--max-seqs", default=None, type=int, help="Max results per query.")
def search_foldseek_cmd(
    structure_id: str, target_db: str | None, target_dir: str | None,
    against_project: bool, sensitivity: float | None, max_seqs: int | None,
) -> None:
    """Search one structure against a Foldseek database, a directory of
    structures, or every other structure already in this project."""
    from proteinexplorer import search as search_mod

    n_targets = sum(x is not None for x in (target_db, target_dir)) + (1 if against_project else 0)
    if n_targets != 1:
        raise click.ClickException("Provide exactly one of --target-db, --target-dir, --against-project")

    try:
        root = proj.find_project_root(".")
        query_path = proj.structure_path(root, structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    tmp_dir = None
    try:
        if against_project:
            query_record = proj.get_record(root, structure_id)
            tmp_dir = Path(tempfile.mkdtemp())
            n_copied = 0
            for record in proj.list_records(root):
                if record.id == query_record.id:
                    continue
                src = proj.structure_path(root, record.id)
                ext = ".pdb" if record.format == "pdb" else ".cif"
                shutil.copyfile(src, tmp_dir / f"{record.name}{ext}")
                n_copied += 1
            if n_copied == 0:
                raise click.ClickException("No other structures in this project to search against")
            target = tmp_dir
        else:
            target = target_db or target_dir

        try:
            hits = search_mod.easy_search(
                query_path, target, sensitivity=sensitivity, max_seqs=max_seqs,
            )
        except search_mod.FoldseekNotAvailableError as exc:
            raise click.ClickException(str(exc)) from exc
        except search_mod.SearchError as exc:
            raise click.ClickException(str(exc)) from exc
    finally:
        if tmp_dir is not None:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    if not hits:
        click.echo("No hits found.")
        return
    click.echo(f"{len(hits)} hit(s):")
    for h in hits:
        extras = []
        if h.alntmscore is not None:
            extras.append(f"TM-score={h.alntmscore:.3f}")
        if h.evalue is not None:
            extras.append(f"e-value={h.evalue:.2e}")
        click.echo(f"  {h.target}" + (f"  ({', '.join(extras)})" if extras else ""))


@search_group.command("createdb")
@click.argument("structure_dir", type=click.Path(exists=True))
@click.argument("db_path", type=click.Path())
def search_createdb_cmd(structure_dir: str, db_path: str) -> None:
    """Build a persistent Foldseek database from a directory of structure
    files, for repeated searches against the same target set."""
    from proteinexplorer import search as search_mod

    try:
        path = search_mod.createdb(structure_dir, db_path)
    except search_mod.FoldseekNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc
    except search_mod.SearchError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Database created at {path}")


@cli.group("fix")
def fix_group() -> None:
    """Structure fixing/cleanup via PDBFixer (`pip install -e ".[fix]"` --
    free and pip-installable, unlike Scwrl4/MODELLER/Foldseek). See
    `prot fix apply --help` for how this relates to `prot model gaps`/
    `prot model loop`."""


@fix_group.command("report")
@click.argument("structure_id")
def fix_report_cmd(structure_id: str) -> None:
    """Show what PDBFixer would find (missing atoms/residues, nonstandard
    residues) without changing anything."""
    from proteinexplorer import fix as fix_mod

    try:
        record, _ = _load(structure_id)
        path = proj.structure_path(proj.find_project_root("."), structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    try:
        analysis = fix_mod.analyze(path)
    except fix_mod.PDBFixerNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"{record.id}  ({record.name})")
    if analysis.missing_residues:
        click.echo("  missing residues (from SEQRES):")
        for label, resnames in analysis.missing_residues.items():
            click.echo(f"    {label}: {', '.join(resnames)}")
    else:
        click.echo("  missing residues: none found (PDBFixer needs SEQRES for this -- "
                    "try `prot model gaps` too, which works from numbering alone)")
    if analysis.missing_atoms:
        click.echo("  incomplete residues (missing atoms):")
        for label, atoms in analysis.missing_atoms.items():
            click.echo(f"    {label}: {', '.join(atoms)}")
    else:
        click.echo("  incomplete residues: none found")
    if analysis.nonstandard_residues:
        click.echo("  nonstandard residues:")
        for label, new_name in analysis.nonstandard_residues:
            click.echo(f"    {label} -> {new_name}")
    else:
        click.echo("  nonstandard residues: none found")


@fix_group.command("apply")
@click.argument("structure_id")
@click.option("--add-missing-residues", is_flag=True,
              help="Also insert whole missing residues PDBFixer's own detection found (needs SEQRES; see `prot fix report`).")
@click.option("--no-missing-atoms", is_flag=True, help="Skip completing residues that have some atoms missing.")
@click.option("--no-replace-nonstandard", is_flag=True, help="Skip normalizing nonstandard residues (e.g. MSE -> MET).")
@click.option("--remove-heterogens", type=click.Choice(["water", "all"]), default=None,
              help="water: remove ions/ligands but keep water. all: remove everything including water.")
@click.option("--add-hydrogens", "ph", default=None, type=float,
              help="Add hydrogens at the given pH (e.g. 7.0). Omit to leave the structure heavy-atom-only.")
@click.option("--name", default=None, help="Name for the fixed structure (default: derived automatically).")
def fix_apply_cmd(
    structure_id: str, add_missing_residues: bool, no_missing_atoms: bool,
    no_replace_nonstandard: bool, remove_heterogens: str | None, ph: float | None, name: str | None,
) -> None:
    """Run PDBFixer and save the result as a new structure in the project
    (the original is left untouched, same as `prot mutate`/`prot model
    loop`). Every step defaults to a conservative choice -- see the
    options above to opt in/out of each one."""
    from proteinexplorer import fix as fix_mod

    try:
        record = proj.get_record(proj.find_project_root("."), structure_id)
        root = proj.find_project_root(".")
        src_path = proj.structure_path(root, structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    fixed_name = name or f"{record.name}_fixed"
    tmp_dir = Path(tempfile.mkdtemp())
    tmp_path = tmp_dir / f"{fixed_name}.pdb"
    try:
        try:
            report = fix_mod.fix(
                src_path, tmp_path,
                add_missing_residues=add_missing_residues,
                add_missing_atoms=not no_missing_atoms,
                replace_nonstandard=not no_replace_nonstandard,
                remove_heterogens=remove_heterogens,
                add_hydrogens_ph=ph,
            )
        except fix_mod.PDBFixerNotAvailableError as exc:
            raise click.ClickException(str(exc)) from exc

        new_record = proj.import_structure(root, tmp_path, name=fixed_name, fmt="pdb")
        proj.log_command(root, current_argv())
    finally:
        tmp_path.unlink(missing_ok=True)
        tmp_dir.rmdir()

    if report.residues_added:
        click.echo(f"  residues added: {sum(len(v) for v in report.residues_added.values())}")
    if report.atoms_added:
        click.echo(f"  residues completed (missing atoms added): {len(report.atoms_added)}")
    if report.nonstandard_replaced:
        click.echo(f"  nonstandard residues replaced: {len(report.nonstandard_replaced)}")
    if report.heterogens_removed:
        click.echo("  heterogens removed")
    if report.hydrogens_added_at_ph is not None:
        click.echo(f"  hydrogens added at pH {report.hydrogens_added_at_ph}")
    click.echo(f"Saved as '{new_record.name}' ({new_record.id})")


@cli.group("valid")
def valid_group() -> None:
    """Structure validation: dependency-free steric clash and bond
    geometry checks, plus an external MolProbity wrapper for the real
    calibrated Ramachandran/rotamer/clashscore analysis (no dependency-
    free substitute exists for that -- see `prot valid molprobity --help`)."""


@valid_group.command("clashes")
@click.argument("structure_id")
@click.option("--selection", default=None, help="Restrict to a selection (default: everything but water).")
@click.option("--tolerance", default=0.4, show_default=True, help="Allowed vdW overlap tolerance (A).")
def valid_clashes_cmd(structure_id: str, selection: str | None, tolerance: float) -> None:
    """Steric (van der Waals) clashes between non-bonded atoms."""
    from proteinexplorer import valid as valid_mod

    _, structure = _contact_load(structure_id)
    atoms = _select_or_fail(structure, selection) if selection else None
    result = valid_mod.clashes(structure, atoms=atoms, tolerance=tolerance)
    if not result:
        click.echo("No clashes found.")
        return
    click.echo(f"{len(result)} clash(es):")
    for c in result:
        click.echo(f"  {c.atom_a}  --  {c.atom_b}   dist={c.distance:.2f} A   overlap={c.overlap:.2f} A")


@valid_group.command("geometry")
@click.argument("structure_id")
def valid_geometry_cmd(structure_id: str) -> None:
    """Backbone bond length/angle outliers vs. standard idealized
    covalent geometry."""
    from proteinexplorer import valid as valid_mod

    _, structure = _contact_load(structure_id)
    result = valid_mod.bond_geometry(structure)
    if not result:
        click.echo("No bond geometry outliers found.")
        return
    click.echo(f"{len(result)} outlier(s):")
    for o in result:
        click.echo(f"  {o.residue}  {o.kind}: {o.value:.2f} (ideal {o.ideal:.2f}, deviation {o.deviation:+.2f})")


@valid_group.command("molprobity")
@click.argument("structure_id")
def valid_molprobity_cmd(structure_id: str) -> None:
    """Run an external MolProbity installation for the full calibrated
    validation report (Ramachandran/rotamer outliers, clashscore,
    CaBLAM, ...). No dependency-free substitute exists for this --
    see `prot valid clashes`/`prot valid geometry` for what's available
    without it."""
    from proteinexplorer import valid as valid_mod

    try:
        path = proj.structure_path(proj.find_project_root("."), structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    try:
        result = valid_mod.molprobity(path)
    except valid_mod.MolProbityNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc

    if result.summary:
        click.echo("Summary (best-effort extraction):")
        for key, value in result.summary.items():
            click.echo(f"  {key}: {value}")
    click.echo("--- full report ---")
    click.echo(result.raw_output)


@cli.group("assembly")
def assembly_group() -> None:
    """Biological assembly generation via gemmi (`pip install -e
    ".[assembly]"` -- free/open-source, pip-installable). Bio.PDB (the
    backend behind every other command) only ever gives you the
    asymmetric unit as deposited; this fills that gap."""


@assembly_group.command("list")
@click.argument("structure_id")
def assembly_list_cmd(structure_id: str) -> None:
    """List the biological assemblies documented in a structure's source
    file (PDB REMARK 350 / mmCIF _pdbx_struct_assembly), without
    generating any of them."""
    from proteinexplorer import assembly as asm_mod

    try:
        path = proj.structure_path(proj.find_project_root("."), structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    try:
        infos = asm_mod.list_assemblies(path)
    except asm_mod.GemmiNotAvailableError as exc:
        raise click.ClickException(str(exc)) from exc

    if not infos:
        click.echo("No biological assembly documented in this file.")
        return
    for info in infos:
        details = info.oligomeric_details or "unspecified"
        click.echo(f"  #{info.name}  {details}  chains={','.join(info.chains_involved)}  operators={info.n_operators}")


@assembly_group.command("generate")
@click.argument("structure_id")
@click.option("--assembly-name", default=None, help="Which documented assembly to expand (default: the first one).")
@click.option("--how", type=click.Choice(["number", "short", "duplicate"]), default="number", show_default=True,
              help="How to name copied chains. 'duplicate' can produce colliding chain IDs.")
@click.option("--name", default=None, help="Name for the expanded structure (default: derived automatically).")
def assembly_generate_cmd(structure_id: str, assembly_name: str | None, how: str, name: str | None) -> None:
    """Expand a documented biological assembly and save it as a new
    structure in the project (the original is left untouched)."""
    from proteinexplorer import assembly as asm_mod

    try:
        record = proj.get_record(proj.find_project_root("."), structure_id)
        root = proj.find_project_root(".")
        src_path = proj.structure_path(root, structure_id)
    except ProjectError as exc:
        raise click.ClickException(str(exc)) from exc

    out_name = name or f"{record.name}_assembly"
    tmp_dir = Path(tempfile.mkdtemp())
    tmp_path = tmp_dir / f"{out_name}.pdb"
    try:
        try:
            result = asm_mod.generate_assembly(src_path, tmp_path, assembly_name=assembly_name, how=how)
        except asm_mod.GemmiNotAvailableError as exc:
            raise click.ClickException(str(exc)) from exc
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc

        new_record = proj.import_structure(root, tmp_path, name=out_name, fmt="pdb")
        proj.log_command(root, current_argv())
    finally:
        tmp_path.unlink(missing_ok=True)
        tmp_dir.rmdir()

    click.echo(f"Assembly '{result.assembly_name}': {result.chains_before} -> {result.chains_after}")
    click.echo(f"Saved as '{new_record.name}' ({new_record.id})")


def main() -> None:
    cli()


if __name__ == "__main__":
    main()
