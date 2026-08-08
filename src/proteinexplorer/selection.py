"""Common selection language shared across all `prot` commands.

Grammar (case-insensitive keywords):

    expr        := or_expr
    or_expr     := and_expr ("or" and_expr)*
    and_expr    := not_expr ("and" not_expr)*
    not_expr    := "not" not_expr | atom_expr
    atom_expr   := "(" expr ")"
                 | "within" NUMBER atom_expr
                 | "protein" | "nucleic" | "water" | "ion" | "ligand"
                 | "backbone" | "sidechain"
                 | "chain" ID
                 | "resid" INT [":" INT]
                 | "resname" NAME
                 | "atom" NAME
                 | "all"

Examples:
    "chain A"
    "protein and not backbone"
    "resid 30:50"
    "within 5 ligand"
    "(chain A and backbone) within 6 ligand"
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Structure import Structure

from proteinexplorer.models import ResidueCategory, classify_residue, is_backbone_atom

_TOKEN_RE = re.compile(
    r"""
    (?P<NUMBER>\d+\.\d+|\d+)
  | (?P<COLON>:)
  | (?P<LPAREN>\()
  | (?P<RPAREN>\))
  | (?P<WORD>[A-Za-z0-9_'\-]+)
  | (?P<WS>\s+)
    """,
    re.VERBOSE,
)

_KEYWORDS = {
    "and", "or", "not", "within", "chain", "resid", "resname", "atom",
    "protein", "nucleic", "water", "ion", "ligand", "backbone", "sidechain", "all",
}


class SelectionSyntaxError(ValueError):
    pass


@dataclass
class Token:
    kind: str
    value: str


def tokenize(text: str) -> list[Token]:
    tokens: list[Token] = []
    pos = 0
    while pos < len(text):
        match = _TOKEN_RE.match(text, pos)
        if not match:
            raise SelectionSyntaxError(f"Unexpected character at position {pos}: {text[pos:pos + 10]!r}")
        pos = match.end()
        kind = match.lastgroup
        value = match.group()
        if kind == "WS":
            continue
        tokens.append(Token(kind, value))
    return tokens


# --- AST nodes -----------------------------------------------------------

class Node:
    def matches(self, ctx: "_EvalContext", atom) -> bool:  # pragma: no cover - interface
        raise NotImplementedError


@dataclass
class All(Node):
    def matches(self, ctx, atom) -> bool:
        return True


@dataclass
class Category(Node):
    category: ResidueCategory

    def matches(self, ctx, atom) -> bool:
        return ctx.category_of(atom) is self.category


@dataclass
class Backbone(Node):
    def matches(self, ctx, atom) -> bool:
        return is_backbone_atom(atom.get_name(), ctx.category_of(atom))


@dataclass
class Sidechain(Node):
    def matches(self, ctx, atom) -> bool:
        category = ctx.category_of(atom)
        if category not in (ResidueCategory.PROTEIN, ResidueCategory.NUCLEIC):
            return False
        return not is_backbone_atom(atom.get_name(), category)


@dataclass
class ChainSel(Node):
    chain_id: str

    def matches(self, ctx, atom) -> bool:
        return atom.get_parent().get_parent().id == self.chain_id


@dataclass
class ResidRange(Node):
    low: int
    high: int

    def matches(self, ctx, atom) -> bool:
        resseq = atom.get_parent().id[1]
        return self.low <= resseq <= self.high


@dataclass
class ResName(Node):
    names: frozenset[str]

    def matches(self, ctx, atom) -> bool:
        return atom.get_parent().resname.strip().upper() in self.names


@dataclass
class AtomName(Node):
    names: frozenset[str]

    def matches(self, ctx, atom) -> bool:
        return atom.get_name().strip().upper() in self.names


@dataclass
class Within(Node):
    distance: float
    inner: Node

    def matches(self, ctx, atom) -> bool:
        inner_atoms = ctx.eval_node(self.inner)
        if not inner_atoms:
            return False
        nearby = ctx.neighbor_search.search(atom.coord, self.distance, level="A")
        return any(a in inner_atoms for a in nearby)


@dataclass
class Not(Node):
    inner: Node

    def matches(self, ctx, atom) -> bool:
        return not self.inner.matches(ctx, atom)


@dataclass
class And(Node):
    left: Node
    right: Node

    def matches(self, ctx, atom) -> bool:
        return self.left.matches(ctx, atom) and self.right.matches(ctx, atom)


@dataclass
class Or(Node):
    left: Node
    right: Node

    def matches(self, ctx, atom) -> bool:
        return self.left.matches(ctx, atom) or self.right.matches(ctx, atom)


# --- Parser ----------------------------------------------------------------

class _Parser:
    def __init__(self, tokens: list[Token]):
        self.tokens = tokens
        self.pos = 0

    def _peek(self) -> Token | None:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def _advance(self) -> Token:
        tok = self._peek()
        if tok is None:
            raise SelectionSyntaxError("Unexpected end of selection expression")
        self.pos += 1
        return tok

    def _expect_word(self, *expected: str) -> str:
        tok = self._advance()
        if tok.kind != "WORD" or tok.value.lower() not in expected:
            raise SelectionSyntaxError(f"Expected one of {expected}, got {tok.value!r}")
        return tok.value.lower()

    def parse(self) -> Node:
        node = self._or_expr()
        if self._peek() is not None:
            raise SelectionSyntaxError(f"Unexpected trailing token: {self._peek().value!r}")
        return node

    def _or_expr(self) -> Node:
        node = self._and_expr()
        while self._peek() and self._peek().kind == "WORD" and self._peek().value.lower() == "or":
            self._advance()
            node = Or(node, self._and_expr())
        return node

    def _and_expr(self) -> Node:
        node = self._not_expr()
        while self._peek() and self._peek().kind == "WORD" and self._peek().value.lower() == "and":
            self._advance()
            node = And(node, self._not_expr())
        return node

    def _not_expr(self) -> Node:
        if self._peek() and self._peek().kind == "WORD" and self._peek().value.lower() == "not":
            self._advance()
            return Not(self._not_expr())
        return self._atom_expr()

    def _atom_expr(self) -> Node:
        tok = self._peek()
        if tok is None:
            raise SelectionSyntaxError("Unexpected end of selection expression")

        if tok.kind == "LPAREN":
            self._advance()
            node = self._or_expr()
            if not (self._peek() and self._peek().kind == "RPAREN"):
                raise SelectionSyntaxError("Missing closing parenthesis")
            self._advance()
            return node

        if tok.kind != "WORD":
            raise SelectionSyntaxError(f"Unexpected token: {tok.value!r}")

        word = tok.value.lower()

        if word == "within":
            self._advance()
            num_tok = self._advance()
            if num_tok.kind != "NUMBER":
                raise SelectionSyntaxError("Expected a number after 'within'")
            inner = self._atom_expr()
            return Within(float(num_tok.value), inner)

        if word == "all":
            self._advance()
            return All()

        if word == "protein":
            self._advance()
            return Category(ResidueCategory.PROTEIN)
        if word == "nucleic":
            self._advance()
            return Category(ResidueCategory.NUCLEIC)
        if word == "water":
            self._advance()
            return Category(ResidueCategory.WATER)
        if word == "ion":
            self._advance()
            return Category(ResidueCategory.ION)
        if word == "ligand":
            self._advance()
            return Category(ResidueCategory.LIGAND)
        if word == "backbone":
            self._advance()
            return Backbone()
        if word == "sidechain":
            self._advance()
            return Sidechain()

        if word == "chain":
            self._advance()
            id_tok = self._advance()
            return ChainSel(id_tok.value)

        if word == "resid":
            self._advance()
            low_tok = self._advance()
            if low_tok.kind != "NUMBER":
                raise SelectionSyntaxError("Expected an integer after 'resid'")
            low = int(float(low_tok.value))
            high = low
            if self._peek() and self._peek().kind == "COLON":
                self._advance()
                high_tok = self._advance()
                if high_tok.kind != "NUMBER":
                    raise SelectionSyntaxError("Expected an integer after 'resid N:'")
                high = int(float(high_tok.value))
            return ResidRange(low, high)

        if word == "resname":
            self._advance()
            names = {self._advance().value.upper()}
            while self._peek() and self._peek().kind == "WORD" and self._peek().value == ",":
                self._advance()
                names.add(self._advance().value.upper())
            return ResName(frozenset(names))

        if word == "atom":
            self._advance()
            names = {self._advance().value.upper()}
            return AtomName(frozenset(names))

        raise SelectionSyntaxError(f"Unknown selection keyword: {word!r}")


def parse_selection(text: str) -> Node:
    tokens = tokenize(text)
    if not tokens:
        raise SelectionSyntaxError("Empty selection expression")
    return _Parser(tokens).parse()


class _EvalContext:
    def __init__(self, structure: Structure):
        model = next(iter(structure))
        self.model = model
        self._category_cache: dict[int, ResidueCategory] = {}
        atoms = list(model.get_atoms())
        self.neighbor_search = NeighborSearch(atoms)
        # Canonical (chain_id, resseq, icode) order, NOT raw file/parse
        # order. Bio.PDB's writer groups all HETATM records after ATOM
        # records regardless of their original position, so after a
        # structure has been saved and reloaded (e.g. by `prot mutate`),
        # a residue like MSE (selenomethionine, hetero-flagged) can end
        # up relocated to the end of the atom list even though its
        # residue number sits in the middle of the sequence. Selections
        # from two independently-parsed structures must still line up
        # atom-for-atom (geometry.rmsd, cluster's pairwise_rmsd_matrix,
        # ...), so we sort canonically here rather than trusting parse
        # order. Sort is stable, so atom order within one residue is
        # preserved.
        atoms.sort(key=lambda a: (a.get_parent().get_parent().id, a.get_parent().id[1], a.get_parent().id[2]))
        self._all_atoms = atoms

    def category_of(self, atom) -> ResidueCategory:
        residue = atom.get_parent()
        key = id(residue)
        cached = self._category_cache.get(key)
        if cached is None:
            cached = classify_residue(residue.resname, residue.id[0])
            self._category_cache[key] = cached
        return cached

    def eval_node(self, node: Node) -> set:
        return {atom for atom in self._all_atoms if node.matches(self, atom)}


def select(structure: Structure, expr: str) -> list:
    """Parse and evaluate a selection expression against a structure's
    first model. Returns a list of Bio.PDB Atom objects in canonical
    (chain_id, residue seqnum, insertion code) order -- this is NOT
    necessarily the raw file/parse order (see _EvalContext for why:
    Bio.PDB's writer relocates HETATM-flagged residues, so file order
    isn't stable across a save/reload round-trip)."""
    node = parse_selection(expr)
    ctx = _EvalContext(structure)
    matched = ctx.eval_node(node)
    return [atom for atom in ctx._all_atoms if atom in matched]
