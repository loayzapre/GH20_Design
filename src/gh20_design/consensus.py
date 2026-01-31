from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional
import shutil
import subprocess


# -------------------------
# FASTA I/O
# -------------------------

def iter_fasta(fasta_path: str | Path) -> Iterable[Tuple[str, str]]:
    fasta_path = Path(fasta_path)
    header = None
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


def write_fasta(records: Iterable[Tuple[str, str]], out_path: str | Path, wrap: int = 80) -> int:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def _wrap(seq: str) -> str:
        if wrap <= 0:
            return seq + "\n"
        return "\n".join(seq[i:i+wrap] for i in range(0, len(seq), wrap)) + "\n"

    n = 0
    with out_path.open("w", encoding="utf-8") as f:
        for h, s in records:
            f.write(f">{h}\n")
            f.write(_wrap(s))
            n += 1
    return n


# -------------------------
# Split FASTA by clusters
# -------------------------

def load_clusters_tsv(clusters_tsv: str | Path) -> Dict[str, List[str]]:
    """
    Expects columns: cluster_id, size, members (comma-separated)
    """
    clusters_tsv = Path(clusters_tsv)
    clusters: Dict[str, List[str]] = {}

    with clusters_tsv.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for req in ("cluster_id", "members"):
            if req not in idx:
                raise ValueError(f"clusters TSV missing column '{req}'")

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            cid = parts[idx["cluster_id"]]
            members = parts[idx["members"]].split(",") if parts[idx["members"]] else []
            clusters[cid] = members

    return clusters


def split_fasta_by_cluster(
    fasta_in: str | Path,
    clusters_tsv: str | Path,
    out_dir: str | Path,
    min_cluster_size: int = 2,
    wrap: int = 80,
) -> dict:
    """
    Writes one FASTA per cluster: out_dir/<cluster_id>.fasta
    Returns summary stats.
    """
    fasta_in = Path(fasta_in)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    clusters = load_clusters_tsv(clusters_tsv)
    clusters = {cid: mem for cid, mem in clusters.items() if len(mem) >= min_cluster_size}

    # Load sequences into dict (NR90 is not too crazy; OK in memory)
    seqs = {}
    for h, s in iter_fasta(fasta_in):
        # node IDs in edges are qseqid = header up to first space (often). Our pipeline wrote full headers.
        # safest: use full header string as ID
        seqs[h] = s

    written_clusters = 0
    missing_members = 0
    total_members = 0

    for cid, members in clusters.items():
        total_members += len(members)
        recs = []
        for m in members:
            # Members come from SSN nodes; in our edge writer we used qseqid/sseqid exactly.
            # That typically equals the FASTA header token up to whitespace.
            # Try exact header first; fallback to first-token matching.
            if m in seqs:
                recs.append((m, seqs[m]))
            else:
                # fallback: match by first token
                hit = None
                for h in seqs.keys():
                    if h.split()[0] == m:
                        hit = h
                        break
                if hit is not None:
                    recs.append((m, seqs[hit]))
                else:
                    missing_members += 1

        if len(recs) >= min_cluster_size:
            write_fasta(recs, out_dir / f"{cid}.fasta", wrap=wrap)
            written_clusters += 1

    return {
        "fasta_in": str(fasta_in),
        "clusters_total": len(clusters),
        "clusters_written": written_clusters,
        "total_members": total_members,
        "missing_members": missing_members,
        "out_dir": str(out_dir),
    }


# -------------------------
# MAFFT + consensus
# -------------------------

@dataclass
class MafftConfig:
    mafft_exe: str = "mafft"
    threads: int = 8
    mode: str = "auto"  # "auto" | "linsi" | "ginsi" | "einsi"
    max_iterate: Optional[int] = None  # e.g. 1000
    localpair: bool = False
    globalpair: bool = False


def _check_exe(name: str) -> str:
    exe = shutil.which(name)
    if exe is None:
        raise FileNotFoundError(f"'{name}' not found on PATH. Install it and retry.")
    return exe


def run_mafft(in_fasta: str | Path, out_aln: str | Path, cfg: MafftConfig) -> dict:
    in_fasta = Path(in_fasta)
    out_aln = Path(out_aln)
    out_aln.parent.mkdir(parents=True, exist_ok=True)

    exe = _check_exe(cfg.mafft_exe)

    cmd = [exe]

    if cfg.mode == "linsi":
        cmd += ["--localpair", "--maxiterate", str(cfg.max_iterate or 1000)]
    elif cfg.mode == "ginsi":
        cmd += ["--globalpair", "--maxiterate", str(cfg.max_iterate or 1000)]
    elif cfg.mode == "einsi":
        cmd += ["--genafpair", "--maxiterate", str(cfg.max_iterate or 1000)]
    else:
        # auto: let mafft decide
        cmd += ["--auto"]

    cmd += ["--thread", str(cfg.threads), str(in_fasta)]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "mafft failed.\n"
            f"Command: {' '.join(cmd)}\n\nSTDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}\n"
        )

    out_aln.write_text(proc.stdout, encoding="utf-8")
    return {"in_fasta": str(in_fasta), "out_aln": str(out_aln)}


def consensus_from_alignment(
    aln_fasta: str | Path,
    out_fasta: str | Path,
    consensus_id: str,
    gap_char: str = "-",
    min_fraction: float = 0.5,
    wrap: int = 80,
) -> dict:
    """
    Simple majority-rule consensus per column.
    - If the most frequent non-gap character frequency < min_fraction, emit 'X'
    - Ignores gaps when counting
    """
    aln_fasta = Path(aln_fasta)
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    seqs = [seq for _, seq in iter_fasta(aln_fasta)]
    if not seqs:
        raise ValueError(f"No sequences in alignment: {aln_fasta}")

    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("Alignment sequences have inconsistent lengths")

    cons = []
    for i in range(L):
        col = [s[i] for s in seqs]
        # count non-gaps
        counts = {}
        non_gap = 0
        for c in col:
            if c == gap_char:
                continue
            non_gap += 1
            c = c.upper()
            counts[c] = counts.get(c, 0) + 1

        if non_gap == 0:
            cons.append(gap_char)
            continue

        aa, cnt = max(counts.items(), key=lambda kv: kv[1])
        if cnt / non_gap >= min_fraction:
            cons.append(aa)
        else:
            cons.append("X")

    cons_seq = "".join(cons)

    # Write
    def _wrap(seq: str) -> str:
        if wrap <= 0:
            return seq + "\n"
        return "\n".join(seq[i:i+wrap] for i in range(0, len(seq), wrap)) + "\n"

    with out_fasta.open("w", encoding="utf-8") as f:
        f.write(f">{consensus_id}\n")
        f.write(_wrap(cons_seq))

    return {"aln": str(aln_fasta), "out": str(out_fasta), "length": len(cons_seq), "min_fraction": min_fraction}
