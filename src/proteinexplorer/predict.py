"""Structure prediction (spec section "Predict").

Ab initio / deep-learning structure prediction has no meaningful
dependency-free approximation (unlike, say, pocket detection or virtual
C-beta placement) -- there is no reasonable geometric heuristic that
substitutes for a trained model. So unlike most other modules in this
package, this one does NOT offer a built-in fallback: it is a thin,
honest wrapper around external tools only, matching the project's stated
philosophy of providing "a unified interface to external structural
biology tools" rather than reimplementing them.

- colabfold_predict(): wraps a local `colabfold_batch` installation.
- alphafold_predict(): wraps a local AlphaFold `run_alphafold.sh`
  installation (requires the multi-hundred-GB AlphaFold parameter/
  sequence databases to be set up separately -- this only invokes the
  script, it does not set up or validate the databases).

Both raise a clear PredictionToolNotAvailableError with installation
pointers when the relevant tool isn't found, rather than attempting any
kind of fake substitute.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


class PredictionToolNotAvailableError(RuntimeError):
    pass


class PredictionError(RuntimeError):
    pass


@dataclass
class PredictionResult:
    tool: str
    output_dir: Path
    models: list[Path]  # sorted best-first when the tool reports a rank
    log_tail: str


def write_fasta(sequence: str, name: str, path: str | Path) -> Path:
    path = Path(path)
    sequence = "".join(sequence.split())  # strip any whitespace/newlines
    path.write_text(f">{name}\n{sequence}\n")
    return path


# --- ColabFold -----------------------------------------------------------

def colabfold_binary() -> str | None:
    return shutil.which("colabfold_batch")


def colabfold_predict(
    sequence: str,
    output_dir: str | Path,
    name: str = "query",
    extra_args: list[str] | None = None,
    timeout: int = 3600,
) -> PredictionResult:
    """Predict a structure via a local ColabFold installation
    (`colabfold_batch`). Requires ColabFold to already be installed
    (https://github.com/sokrypton/ColabFold) -- this does not download
    weights or set anything up, it only invokes the binary.
    """
    binary = colabfold_binary()
    if binary is None:
        raise PredictionToolNotAvailableError(
            "colabfold_batch not found on PATH. Install ColabFold "
            "(https://github.com/sokrypton/ColabFold) to predict structures "
            "locally -- there is no dependency-free fallback for structure "
            "prediction."
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / f"{name}.fasta"
    write_fasta(sequence, name, fasta_path)

    cmd = [binary, str(fasta_path), str(output_dir)] + (extra_args or [])
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        raise PredictionError(f"colabfold_batch failed (exit {result.returncode}): {result.stderr[-2000:]}")

    models = sorted(output_dir.glob(f"{name}*rank_001*.pdb"))
    if not models:
        models = sorted(output_dir.glob(f"{name}*.pdb"))
    return PredictionResult(
        tool="colabfold", output_dir=output_dir, models=models,
        log_tail=result.stdout[-2000:],
    )


# --- AlphaFold -------------------------------------------------------------

def alphafold_predict(
    fasta_path: str | Path,
    output_dir: str | Path,
    alphafold_script: str | Path,
    data_dir: str | Path,
    max_template_date: str,
    model_preset: str = "monomer",
    db_preset: str = "full_dbs",
    timeout: int = 21600,
) -> PredictionResult:
    """Predict a structure via a local AlphaFold installation's
    `run_alphafold.sh`. Requires AlphaFold itself plus its (very large)
    parameter and sequence databases to already be set up
    (https://github.com/google-deepmind/alphafold) -- this only invokes
    the script with the given paths, it does not download or validate
    the databases.
    """
    script = Path(alphafold_script)
    if not script.exists():
        raise PredictionToolNotAvailableError(
            f"AlphaFold script not found at {script}. Install AlphaFold "
            "(https://github.com/google-deepmind/alphafold) and its parameter/"
            "sequence databases to predict structures locally -- there is no "
            "dependency-free fallback for structure prediction."
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(script),
        f"--fasta_paths={fasta_path}",
        f"--output_dir={output_dir}",
        f"--data_dir={data_dir}",
        f"--max_template_date={max_template_date}",
        f"--model_preset={model_preset}",
        f"--db_preset={db_preset}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        raise PredictionError(f"AlphaFold failed (exit {result.returncode}): {result.stderr[-2000:]}")

    models = sorted(output_dir.rglob("ranked_*.pdb"))
    return PredictionResult(
        tool="alphafold", output_dir=output_dir, models=models,
        log_tail=result.stdout[-2000:],
    )
