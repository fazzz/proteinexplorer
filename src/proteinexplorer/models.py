"""Domain model layer on top of Bio.PDB.

Bio.PDB already provides the Structure/Model/Chain/Residue/Atom hierarchy,
so this module intentionally does not reimplement it. What it adds is the
classification Bio.PDB does not give out of the box (protein / nucleic acid
/ water / ion / ligand), plus small helpers used by the selection language
and by higher-level analysis commands.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


# --- Residue name tables ----------------------------------------------

# Standard 20 amino acids + common non-standard residues that still count
# as "protein" for classification purposes (selenomethionine, etc).
STANDARD_AMINO_ACIDS: frozenset[str] = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    }
)
NONSTANDARD_PROTEIN_RESIDUES: frozenset[str] = frozenset(
    {
        "MSE",  # selenomethionine
        "SEC",  # selenocysteine
        "PYL",  # pyrrolysine
        "HYP",  # hydroxyproline
        "CSO", "CSD",  # oxidized cysteine variants
        "PTR", "SEP", "TPO",  # phosphotyrosine/serine/threonine
    }
)
PROTEIN_RESIDUES: frozenset[str] = STANDARD_AMINO_ACIDS | NONSTANDARD_PROTEIN_RESIDUES

NUCLEIC_RESIDUES: frozenset[str] = frozenset(
    {
        "DA", "DC", "DG", "DT", "DU",  # DNA
        "A", "C", "G", "U",  # RNA (single-letter mmCIF names)
        "RA", "RC", "RG", "RU",  # RNA (some PDB variants)
    }
)

WATER_RESIDUES: frozenset[str] = frozenset({"HOH", "WAT", "H2O", "DOD"})

# Common monatomic/small polyatomic ions seen as crystallographic additives.
ION_RESIDUES: frozenset[str] = frozenset(
    {
        "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "FE2", "CU", "CU1",
        "NI", "CO", "CD", "HG", "BA", "CS", "LI", "AL", "PB", "SR",
        "SO4", "PO4", "NO3", "NH4", "BR", "IOD",
    }
)

# Standard protein backbone atom names.
PROTEIN_BACKBONE_ATOMS: frozenset[str] = frozenset({"N", "CA", "C", "O", "OXT"})
# Standard nucleic acid backbone (sugar-phosphate) atom names.
NUCLEIC_BACKBONE_ATOMS: frozenset[str] = frozenset(
    {
        "P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'",
        "C2'", "O2'", "C1'",
    }
)


class ResidueCategory(str, Enum):
    PROTEIN = "protein"
    NUCLEIC = "nucleic"
    WATER = "water"
    ION = "ion"
    LIGAND = "ligand"  # anything else with HETATM records (small molecules)


def classify_residue(resname: str, hetero_flag: str) -> ResidueCategory:
    """Classify a residue using its 3/1-letter name and Bio.PDB hetero flag.

    ``hetero_flag`` is the first element of Bio.PDB's residue id tuple:
    ``" "`` for standard ATOM records, ``"W"`` for water, or ``"H_XXX"``
    for other HETATM groups.
    """
    name = resname.strip().upper()

    if hetero_flag == "W" or name in WATER_RESIDUES:
        return ResidueCategory.WATER
    if name in PROTEIN_RESIDUES:
        return ResidueCategory.PROTEIN
    if name in NUCLEIC_RESIDUES:
        return ResidueCategory.NUCLEIC
    if name in ION_RESIDUES:
        return ResidueCategory.ION
    return ResidueCategory.LIGAND


def is_backbone_atom(atom_name: str, category: ResidueCategory) -> bool:
    name = atom_name.strip().upper()
    if category is ResidueCategory.PROTEIN:
        return name in PROTEIN_BACKBONE_ATOMS
    if category is ResidueCategory.NUCLEIC:
        return name in NUCLEIC_BACKBONE_ATOMS
    return False


@dataclass
class ChainSummary:
    chain_id: str
    n_residues: int
    categories: dict[str, int] = field(default_factory=dict)  # category -> count


@dataclass
class StructureSummary:
    """Lightweight, serializable summary of a parsed structure.

    Produced by io.summarize() and used by `prot status` / `prot info`
    without needing to keep the full Bio.PDB object tree in memory.
    """

    n_models: int
    chains: list[ChainSummary]
    n_atoms: int
    hetero_resnames: list[str]
    has_altloc: bool

    @property
    def n_chains(self) -> int:
        return len(self.chains)

    @property
    def n_residues(self) -> int:
        return sum(c.n_residues for c in self.chains)

    def category_totals(self) -> dict[str, int]:
        totals: dict[str, int] = {}
        for chain in self.chains:
            for category, count in chain.categories.items():
                totals[category] = totals.get(category, 0) + count
        return totals
