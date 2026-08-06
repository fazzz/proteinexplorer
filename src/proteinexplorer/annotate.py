"""Annotation (spec section "Annotation").

Two layers, same split as BioExplorer's annotate.py/annotate_external.py:

- Built-in, no network needed: metal-binding site detection (purely
  geometric, from the structure's own ion + protein-atom coordinates) and
  experimental metadata (already embedded in the PDB/mmCIF header).
- External DB/API lookups: UniProt (protein-level: gene name, organism,
  taxonomy, EC numbers, GO terms) and Pfam domains via the InterPro REST
  API. Both take an accession the caller already knows (there's no local
  PDB-ID -> UniProt mapping step here) and make a single REST GET each --
  no job-polling APIs (like InterProScan sequence submission) are wired
  up in this pass. Network errors (including an execution sandbox with no
  route to these hosts) surface as a clear AnnotationFetchError rather
  than a stack trace.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from dataclasses import dataclass, field

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Structure import Structure

from proteinexplorer import io as pio
from proteinexplorer.models import ResidueCategory, classify_residue

UNIPROT_REST_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
INTERPRO_PFAM_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/"

_METAL_COORDINATING_ELEMENTS = frozenset({"O", "N", "S"})


def _residue_label(residue) -> str:
    chain_id = residue.get_parent().id
    return f"{chain_id}/{residue.resname}{residue.id[1]}"


# --- Built-in: metal-binding sites ----------------------------------------

@dataclass
class MetalBindingSite:
    ion_label: str
    ion_element: str
    coordinating_residues: list[str]
    coordinating_distances: list[float]


def metal_binding_sites(structure: Structure, cutoff: float = 3.0) -> list[MetalBindingSite]:
    """Ions in the structure plus the protein residues coordinating them
    (any O/N/S atom within `cutoff` Angstrom of the ion) -- purely
    geometric, no external database needed."""
    model = next(iter(structure))
    all_atoms = list(model.get_atoms())
    if not all_atoms:
        return []
    ns = NeighborSearch(all_atoms)

    sites: list[MetalBindingSite] = []
    for chain in model:
        for residue in chain:
            if classify_residue(residue.resname, residue.id[0]) is not ResidueCategory.ION:
                continue
            for ion_atom in residue:
                nearby = ns.search(ion_atom.coord, cutoff, level="A")
                coordinating: dict[str, tuple[str, float]] = {}
                for atom in nearby:
                    other_residue = atom.get_parent()
                    if classify_residue(other_residue.resname, other_residue.id[0]) is not ResidueCategory.PROTEIN:
                        continue
                    if (atom.element or "").strip().upper() not in _METAL_COORDINATING_ELEMENTS:
                        continue
                    d = float(((atom.coord - ion_atom.coord) ** 2).sum() ** 0.5)
                    label = _residue_label(other_residue)
                    if label not in coordinating or d < coordinating[label][1]:
                        coordinating[label] = (atom.get_name(), d)

                if coordinating:
                    ordered = sorted(coordinating.items(), key=lambda kv: kv[1][1])
                    sites.append(
                        MetalBindingSite(
                            ion_label=_residue_label(residue),
                            ion_element=(ion_atom.element or "").strip(),
                            coordinating_residues=[label for label, _ in ordered],
                            coordinating_distances=[d for _, (_, d) in ordered],
                        )
                    )
    return sites


# --- Built-in: experimental metadata ---------------------------------------

@dataclass
class StructureMetadata:
    method: str | None
    resolution: float | None
    deposition_date: str | None


def structure_metadata(structure: Structure) -> StructureMetadata:
    header = pio.header_info(structure)
    return StructureMetadata(
        method=header["structure_method"],
        resolution=header["resolution"],
        deposition_date=header["deposition_date"],
    )


# --- External: UniProt -----------------------------------------------------

class AnnotationFetchError(RuntimeError):
    pass


@dataclass
class UniProtAnnotation:
    accession: str
    gene_names: list[str] = field(default_factory=list)
    organism: str | None = None
    taxonomy_id: int | None = None
    ec_numbers: list[str] = field(default_factory=list)
    go_terms: list[str] = field(default_factory=list)


def _http_get_json(url: str, timeout: float = 15.0) -> dict:
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise AnnotationFetchError(f"HTTP {exc.code} fetching {url}: {exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise AnnotationFetchError(f"Could not reach {url}: {exc.reason}") from exc
    except (TimeoutError, ConnectionError) as exc:
        raise AnnotationFetchError(f"Timed out fetching {url}: {exc}") from exc


def uniprot_lookup(accession: str) -> UniProtAnnotation:
    """Gene name, organism, taxonomy, EC numbers, and GO terms for a
    UniProt accession, via the UniProt REST API."""
    data = _http_get_json(UNIPROT_REST_URL.format(accession=accession))

    gene_names = [
        g.get("geneName", {}).get("value")
        for g in data.get("genes", [])
        if g.get("geneName", {}).get("value")
    ]

    organism_data = data.get("organism", {})
    organism = organism_data.get("scientificName")
    taxonomy_id = organism_data.get("taxonId")

    ec_numbers = []
    for desc in data.get("proteinDescription", {}).get("recommendedName", {}).get("ecNumbers", []):
        if desc.get("value"):
            ec_numbers.append(desc["value"])

    go_terms = [
        xref.get("id")
        for xref in data.get("uniProtKBCrossReferences", [])
        if xref.get("database") == "GO" and xref.get("id")
    ]

    return UniProtAnnotation(
        accession=accession, gene_names=gene_names, organism=organism,
        taxonomy_id=taxonomy_id, ec_numbers=ec_numbers, go_terms=go_terms,
    )


# --- External: Pfam domains (via InterPro) ---------------------------------

@dataclass
class PfamDomain:
    accession: str
    name: str


def pfam_domains(uniprot_accession: str) -> list[PfamDomain]:
    """Pfam domain hits for a UniProt accession, via the InterPro REST API
    (a single direct query, not the job-polling InterProScan sequence
    submission API)."""
    data = _http_get_json(INTERPRO_PFAM_URL.format(accession=uniprot_accession))
    results = []
    for entry in data.get("results", []):
        metadata = entry.get("metadata", {})
        acc = metadata.get("accession")
        name = metadata.get("name")
        if acc:
            results.append(PfamDomain(accession=acc, name=name or ""))
    return results
