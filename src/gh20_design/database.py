from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Tuple, Optional, List
import shutil
import subprocess


# -------------------------
# FASTA I/O (simple, fast)
# -------------------------

def iter_fasta(fasta_path: str | Path) -> Iterator[Tuple[str, str]]:
    """Yield (header_without_>, sequence) from a FASTA file."""
    fasta_path = Path(fasta_path)
    header: Optional[str] = None
    seq_chunks: List[str] = []

    with fasta_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            yield header, "".join(seq_chunks)


def write_fasta(records: Iterator[Tuple[str, str]], out_path: str | Path, wrap: int = 80) -> int:
    """
    Write FASTA records and return how many sequences were written.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def _wrap(seq: str) -> str:
        if wrap <= 0:
            return seq + "\n"
        return "\n".join(seq[i:i+wrap] for i in range(0, len(seq), wrap)) + "\n"

    n = 0
    with out_path.open("w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            f.write(_wrap(seq))
            n += 1
    return n


# -------------------------
# Length filtering
# -------------------------

@dataclass
class LengthFilterConfig:
    min_len: int = 200
    max_len: Optional[int] = None
    drop_ambiguous: bool = False  # drop sequences containing non-standard AA chars
    wrap: int = 80


def _is_standard_protein(seq: str) -> bool:
    """
    Very simple protein alphabet check.
    Allows X by default; adjust if needed.
    """
    allowed = set("ACDEFGHIKLMNPQRSTVWYXBZU")  # includes X, B, Z, U
    return all((c in allowed) for c in seq.upper())


def length_filter_fasta(
    in_fasta: str | Path,
    out_fasta: str | Path,
    cfg: LengthFilterConfig,
) -> dict:
    """
    Filter sequences by length (and optionally by allowed alphabet).

    Returns summary dict.
    """
    in_fasta = Path(in_fasta)
    out_fasta = Path(out_fasta)

    kept = 0
    dropped_short = 0
    dropped_long = 0
    dropped_amb = 0
    total = 0

    def gen():
        nonlocal kept, dropped_short, dropped_long, dropped_amb, total
        for header, seq in iter_fasta(in_fasta):
            total += 1
            L = len(seq)
            if L < cfg.min_len:
                dropped_short += 1
                continue
            if cfg.max_len is not None and L > cfg.max_len:
                dropped_long += 1
                continue
            if cfg.drop_ambiguous and not _is_standard_protein(seq):
                dropped_amb += 1
                continue
            kept += 1
            yield header, seq

    write_fasta(gen(), out_fasta, wrap=cfg.wrap)

    return {
        "input": str(in_fasta),
        "output": str(out_fasta),
        "total": total,
        "kept": kept,
        "dropped_short": dropped_short,
        "dropped_long": dropped_long,
        "dropped_ambiguous": dropped_amb,
        "min_len": cfg.min_len,
        "max_len": cfg.max_len,
        "drop_ambiguous": cfg.drop_ambiguous,
    }


# -------------------------
# CD-HIT wrapper
# -------------------------

@dataclass
class CDHitConfig:
    identity: float = 0.90  # -c
    word_length: Optional[int] = None  # -n (auto if None)
    threads: int = 8  # -T
    memory_mb: int = 16000  # -M (MB)
    description_len: int = 0  # -d 0 keeps full header (recommended)
    extra_args: Optional[list[str]] = None


def _auto_word_length(identity: float) -> int:
    """
    Reasonable defaults for cd-hit -n given -c.
    cd-hit requires specific combos; these are common choices:
      -c 0.7-1.0 with -n 2..5
    """
    if identity >= 0.9:
        return 5
    if identity >= 0.88:
        return 5
    if identity >= 0.85:
        return 4
    if identity >= 0.8:
        return 4
    return 3


def run_cdhit(
    in_fasta: str | Path,
    out_fasta: str | Path,
    cfg: CDHitConfig,
    executable: str = "cd-hit",
) -> dict:
    """
    Run cd-hit on a protein FASTA.

    Produces:
      out_fasta
      out_fasta + ".clstr"
    """
    in_fasta = Path(in_fasta)
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    exe = shutil.which(executable)
    if exe is None:
        raise FileNotFoundError(
            f"'{executable}' not found on PATH. Install cd-hit or load the module, then retry."
        )

    n = cfg.word_length if cfg.word_length is not None else _auto_word_length(cfg.identity)

    cmd = [
        exe,
        "-i", str(in_fasta),
        "-o", str(out_fasta),
        "-c", str(cfg.identity),
        "-n", str(n),
        "-T", str(cfg.threads),
        "-M", str(cfg.memory_mb),
        "-d", str(cfg.description_len),
    ]
    if cfg.extra_args:
        cmd.extend(cfg.extra_args)

    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        raise RuntimeError(
            "cd-hit failed.\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"STDOUT:\n{proc.stdout}\n\n"
            f"STDERR:\n{proc.stderr}\n"
        )

    return {
        "input": str(in_fasta),
        "output": str(out_fasta),
        "clusters_file": str(out_fasta) + ".clstr",
        "identity": cfg.identity,
        "word_length": n,
        "threads": cfg.threads,
        "memory_mb": cfg.memory_mb,
        "stdout_tail": proc.stdout[-1000:],
        "stderr_tail": proc.stderr[-1000:],
    }
