from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import shutil
import subprocess


@dataclass
class DiamondConfig:
    diamond_exe: str = "diamond"
    threads: int = 8
    tmpdir: Optional[str] = None  # e.g. "/tmp" or a fast scratch
    block_size: Optional[float] = None  # GB, diamond --block-size
    index_chunks: Optional[int] = None  # diamond --index-chunks


@dataclass
class AllVsAllConfig:
    evalue: float = 1e-5
    max_target_seqs: int = 5000  # higher = denser SSN; can be big
    query_cover: Optional[float] = None  # percent, e.g. 50
    subject_cover: Optional[float] = None  # percent, e.g. 50
    min_identity: Optional[float] = None  # percent, e.g. 30
    mask_low_complexity: bool = True  # diamond default is masking; we expose it


def _check_exe(name: str) -> str:
    exe = shutil.which(name)
    if exe is None:
        raise FileNotFoundError(f"'{name}' not found on PATH. Install it (e.g. micromamba) and retry.")
    return exe


def make_diamond_db(fasta_in: str | Path, db_prefix: str | Path, cfg: DiamondConfig) -> dict:
    """
    diamond makedb --in FASTA --db PREFIX
    Produces PREFIX.dmnd
    """
    fasta_in = Path(fasta_in)
    db_prefix = Path(db_prefix)
    db_prefix.parent.mkdir(parents=True, exist_ok=True)

    exe = _check_exe(cfg.diamond_exe)

    cmd = [exe, "makedb", "--in", str(fasta_in), "--db", str(db_prefix)]
    if cfg.threads:
        cmd += ["--threads", str(cfg.threads)]
    if cfg.tmpdir:
        cmd += ["--tmpdir", str(cfg.tmpdir)]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "diamond makedb failed.\n"
            f"Command: {' '.join(cmd)}\n\nSTDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}\n"
        )

    return {
        "fasta_in": str(fasta_in),
        "db_prefix": str(db_prefix),
        "db_file": str(db_prefix) + ".dmnd",
        "stdout_tail": proc.stdout[-1000:],
        "stderr_tail": proc.stderr[-1000:],
    }


def run_all_vs_all(
    fasta_in: str | Path,
    db_prefix: str | Path,
    out_tsv: str | Path,
    diamond_cfg: DiamondConfig,
    av_cfg: AllVsAllConfig,
) -> dict:
    """
    diamond blastp query vs itself (db built from same FASTA).
    Output format: tabular with fields needed for edge filtering.
    """
    fasta_in = Path(fasta_in)
    db_prefix = Path(db_prefix)
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    exe = _check_exe(diamond_cfg.diamond_exe)

    # IMPORTANT: this DIAMOND expects BLAST-style tokens:
    # --outfmt 6 qseqid sseqid ...
    outfmt_fields = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
    ]

    cmd = [
        exe, "blastp",
        "--query", str(fasta_in),
        "--db", str(db_prefix),
        "--out", str(out_tsv),
        "--outfmt", "6", *outfmt_fields,
        "--evalue", str(av_cfg.evalue),
        "--max-target-seqs", str(av_cfg.max_target_seqs),
        "--threads", str(diamond_cfg.threads),
    ]

    if diamond_cfg.tmpdir:
        cmd += ["--tmpdir", str(diamond_cfg.tmpdir)]
    if diamond_cfg.block_size is not None:
        cmd += ["--block-size", str(diamond_cfg.block_size)]
    if diamond_cfg.index_chunks is not None:
        cmd += ["--index-chunks", str(diamond_cfg.index_chunks)]

    if av_cfg.query_cover is not None:
        cmd += ["--query-cover", str(av_cfg.query_cover)]
    if av_cfg.subject_cover is not None:
        cmd += ["--subject-cover", str(av_cfg.subject_cover)]
    if av_cfg.min_identity is not None:
        cmd += ["--id", str(av_cfg.min_identity)]

    if not av_cfg.mask_low_complexity:
        cmd += ["--masking", "0"]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "diamond blastp failed.\n"
            f"Command: {' '.join(cmd)}\n\nSTDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}\n"
        )

    return {
        "query_fasta": str(fasta_in),
        "db_prefix": str(db_prefix),
        "out_tsv": str(out_tsv),
        "evalue": av_cfg.evalue,
        "max_target_seqs": av_cfg.max_target_seqs,
        "stdout_tail": proc.stdout[-1000:],
        "stderr_tail": proc.stderr[-1000:],
    }



@dataclass
class EdgeFilterConfig:
    bitscore_min: float = 100.0
    evalue_max: float = 1e-5
    min_align_len: int = 0
    min_identity: float = 0.0  # percent
    min_qcov: float = 0.0      # percent
    min_scov: float = 0.0      # percent
    drop_self_hits: bool = True


def filter_edges(
    hits_tsv: str | Path,
    out_edges_tsv: str | Path,
    cfg: EdgeFilterConfig,
) -> dict:
    """
    Convert DIAMOND hits TSV -> SSN edge list TSV.

    Input columns (from run_all_vs_all outfmt):
      qseqid sseqid pident length ... evalue bitscore qlen slen

    Output columns:
      source  target  bitscore  evalue  pident  alnlen  qcov  scov
    """
    hits_tsv = Path(hits_tsv)
    out_edges_tsv = Path(out_edges_tsv)
    out_edges_tsv.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    total = 0

    with hits_tsv.open("r", encoding="utf-8") as fin, out_edges_tsv.open("w", encoding="utf-8") as fout:
        fout.write("source\ttarget\tbitscore\tevalue\tpident\talnlen\tqcov\tscov\n")

        for line in fin:
            line = line.strip()
            if not line:
                continue
            total += 1

            parts = line.split("\t")
            if len(parts) < 14:
                continue

            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            alnlen = int(parts[3])
            evalue = float(parts[10])
            bitscore = float(parts[11])
            qlen = float(parts[12])
            slen = float(parts[13])

            if cfg.drop_self_hits and qseqid == sseqid:
                continue

            if bitscore < cfg.bitscore_min:
                continue
            if evalue > cfg.evalue_max:
                continue
            if alnlen < cfg.min_align_len:
                continue
            if pident < cfg.min_identity:
                continue

            qcov = 100.0 * alnlen / qlen if qlen > 0 else 0.0
            scov = 100.0 * alnlen / slen if slen > 0 else 0.0

            if qcov < cfg.min_qcov:
                continue
            if scov < cfg.min_scov:
                continue

            fout.write(
                f"{qseqid}\t{sseqid}\t{bitscore}\t{evalue}\t{pident}\t{alnlen}\t{qcov:.2f}\t{scov:.2f}\n"
            )
            kept += 1

    return {
        "hits_tsv": str(hits_tsv),
        "edges_tsv": str(out_edges_tsv),
        "total_hits": total,
        "kept_edges": kept,
        "bitscore_min": cfg.bitscore_min,
        "evalue_max": cfg.evalue_max,
        "min_align_len": cfg.min_align_len,
        "min_identity": cfg.min_identity,
        "min_qcov": cfg.min_qcov,
        "min_scov": cfg.min_scov,
        "drop_self_hits": cfg.drop_self_hits,
    }
