from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from datetime import date
from pathlib import Path
from typing import Any, Optional


@dataclass(frozen=True)
class RunParams:
    """
    Core parameters that define a run.

    Notes:
    - nr: CD-HIT identity (0.90 -> "nr90")
    - bitscore: SSN edge threshold
    - inflation: MCL inflation
    - max_targets: DIAMOND --max-target-seqs
    """
    nr: float
    bitscore: float
    inflation: float
    max_targets: int


def make_run_id(p: RunParams, prefix: Optional[str] = None) -> str:
    """
    Create a filesystem-friendly run id that encodes key parameters.
    Example: gh20__2026-01-31__nr90__bs100__I2.0__t5000
    """
    nr_tag = f"nr{int(round(p.nr * 100))}"
    bs_tag = f"bs{int(round(p.bitscore))}"
    I_tag = f"I{p.inflation}"
    t_tag = f"t{int(p.max_targets)}"

    base = f"{date.today().isoformat()}__{nr_tag}__{bs_tag}__{I_tag}__{t_tag}"
    return f"{prefix}__{base}" if prefix else base


def make_run_dirs(root: str | Path, run_id: str) -> dict[str, Path]:
    """
    Create and return the standard directory layout for a run.

    returns dict keys:
      run, acquire, clean, ssn, cluster, consensus, figures, logs
    """
    root = Path(root)
    run = root / "runs" / run_id

    d = {
        "run": run,
        "acquire": run / "01_acquire",
        "clean": run / "02_clean",
        "ssn": run / "03_ssn",
        "cluster": run / "04_cluster",
        "consensus": run / "05_consensus",
        "figures": run / "figures",
        "logs": run / "logs",
    }

    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)

    return d


def write_run_config(run_dir: str | Path, params: RunParams, extra: Optional[dict[str, Any]] = None) -> Path:
    """
    Write a JSON file with run parameters for reproducibility.
    Output: <run_dir>/run_config.json
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "params": asdict(params),
        "extra": extra or {},
    }

    out = run_dir / "run_config.json"
    with out.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)

    return out

def load_run_config(run_dir: str | Path) -> RunParams:
    """
    Load RunParams from <run_dir>/run_config.json.
    This is the single source of truth for parameters.
    """
    run_dir = Path(run_dir)
    cfg_path = run_dir / "run_config.json"

    if not cfg_path.exists():
        raise FileNotFoundError(f"run_config.json not found in {run_dir}")

    with cfg_path.open("r", encoding="utf-8") as f:
        payload = json.load(f)

    params = payload.get("params")
    if params is None:
        raise ValueError("run_config.json missing 'params' field")

    return RunParams(**params)
