"""Structural similarity search (spec section "Search" -- planned in the
original ProteinExplorer spec draft but never implemented until now).

Wraps Foldseek (https://github.com/steineggerlab/foldseek), which is the
right tool for this: fast structural search over large databases via its
3Di structural alphabet. There is no meaningful dependency-free
approximation for "search a large structural database" -- like
predict.py and model.py's homology command, this is an honest,
external-tool-only wrapper. Foldseek isn't pip-installable (it ships as a
compiled binary), so this raises a clear FoldseekNotAvailableError with
install pointers when it's missing, same pattern as Scwrl4/MODELLER/
TMalign.

`easy_search()` wraps `foldseek easy-search`, whose target can be either
a pre-built Foldseek database (see `createdb()`) or a plain directory of
PDB/mmCIF files -- Foldseek builds a temporary database from a directory
on the fly, so a persistent createdb step is only needed for repeated
searches against the same target set.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_FORMAT_OUTPUT = (
    "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
    "evalue,bits,alntmscore,lddt,rmsd"
)


class FoldseekNotAvailableError(RuntimeError):
    pass


class SearchError(RuntimeError):
    pass


def foldseek_binary() -> str | None:
    return shutil.which("foldseek")


@dataclass
class SearchHit:
    """One search hit. `fields` holds every column from the run's
    --format-output, keyed by column name, as strings -- exactly what
    Foldseek printed, so nothing here silently assumes a column exists
    or guesses at its type beyond the few convenience properties below.
    """

    fields: dict[str, str] = field(default_factory=dict)

    @property
    def query(self) -> str | None:
        return self.fields.get("query")

    @property
    def target(self) -> str | None:
        return self.fields.get("target")

    @property
    def evalue(self) -> float | None:
        value = self.fields.get("evalue")
        return float(value) if value is not None else None

    @property
    def bits(self) -> float | None:
        value = self.fields.get("bits")
        return float(value) if value is not None else None

    @property
    def alntmscore(self) -> float | None:
        value = self.fields.get("alntmscore")
        return float(value) if value is not None else None


def _require_binary() -> str:
    binary = foldseek_binary()
    if binary is None:
        raise FoldseekNotAvailableError(
            "foldseek executable not found on PATH. Install Foldseek "
            "(https://github.com/steineggerlab/foldseek#installation -- "
            "conda/homebrew/static binary, not on PyPI) to search structural "
            "databases -- there is no dependency-free fallback for large-scale "
            "structural similarity search."
        )
    return binary


def _parse_tabular(output_text: str, columns: list[str]) -> list[SearchHit]:
    hits = []
    for line in output_text.splitlines():
        line = line.strip()
        if not line:
            continue
        values = line.split("\t")
        hits.append(SearchHit(fields=dict(zip(columns, values))))
    return hits


def easy_search(
    query_path: str | Path,
    target: str | Path,
    format_output: str = DEFAULT_FORMAT_OUTPUT,
    sensitivity: float | None = None,
    max_seqs: int | None = None,
    extra_args: list[str] | None = None,
    timeout: int = 1800,
) -> list[SearchHit]:
    """Search `query_path` (a single structure file) against `target` (a
    Foldseek database, or a directory of structure files) via
    `foldseek easy-search`. Returns hits sorted the way Foldseek printed
    them (best-first for the default e-value/bit-score ranking)."""
    binary = _require_binary()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        output_path = tmp_path / "result.tsv"
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        cmd = [
            binary, "easy-search", str(query_path), str(target), str(output_path), str(scratch),
            "--format-output", format_output,
        ]
        if sensitivity is not None:
            cmd += ["-s", str(sensitivity)]
        if max_seqs is not None:
            cmd += ["--max-seqs", str(max_seqs)]
        cmd += extra_args or []

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            raise SearchError(f"foldseek easy-search failed (exit {result.returncode}): {result.stderr[-2000:]}")

        if not output_path.exists():
            return []
        text = output_path.read_text()

    return _parse_tabular(text, format_output.split(","))


def createdb(structure_dir: str | Path, db_path: str | Path, timeout: int = 1800) -> Path:
    """Build a persistent Foldseek database from a directory of structure
    files, for repeated searches against the same target set (faster than
    letting `easy_search` rebuild a temporary database every call)."""
    binary = _require_binary()
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [binary, "createdb", str(structure_dir), str(db_path)],
        capture_output=True, text=True, timeout=timeout,
    )
    if result.returncode != 0:
        raise SearchError(f"foldseek createdb failed (exit {result.returncode}): {result.stderr[-2000:]}")
    return db_path
