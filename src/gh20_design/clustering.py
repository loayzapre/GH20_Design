from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import math
import shutil
import subprocess


@dataclass
class EdgeToMCLConfig:
    weight_mode: str = "bitscore"   # "bitscore" only (for now)
    transform: str = "log10"        # "none" | "log10" | "sqrt"
    undirected: bool = True
    drop_self: bool = True


def edges_tsv_to_abc_stream(
    edges_tsv: str | Path,
    out_abc: str | Path,
    cfg: EdgeToMCLConfig,
) -> dict:
    """
    Convert edges TSV produced by filter_edges() into MCL abc format:
      nodeA <tab> nodeB <tab> weight

    Expected input columns:
      source target bitscore evalue pident alnlen qcov scov
    """
    edges_tsv = Path(edges_tsv)
    out_abc = Path(out_abc)
    out_abc.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    written = 0

    with edges_tsv.open("r", encoding="utf-8") as fin, out_abc.open("w", encoding="utf-8") as fout:
        header = fin.readline().strip().split("\t")
        idx = {name: i for i, name in enumerate(header)}

        req = ["source", "target", "bitscore"]
        for r in req:
            if r not in idx:
                raise ValueError(f"Missing required column '{r}' in {edges_tsv}")

        for line in fin:
            line = line.strip()
            if not line:
                continue
            total += 1
            parts = line.split("\t")

            s = parts[idx["source"]]
            t = parts[idx["target"]]

            if cfg.drop_self and s == t:
                continue

            w = float(parts[idx["bitscore"]])

            if cfg.transform == "log10":
                w = math.log10(max(w, 1e-12))
            elif cfg.transform == "sqrt":
                w = math.sqrt(max(w, 0.0))
            elif cfg.transform == "none":
                pass
            else:
                raise ValueError("transform must be one of: none, log10, sqrt")

            fout.write(f"{s}\t{t}\t{w}\n")
            written += 1
            if cfg.undirected:
                fout.write(f"{t}\t{s}\t{w}\n")
                written += 1

    return {"total_edges_in": total, "lines_written": written, "out_abc": str(out_abc), "transform": cfg.transform}


@dataclass
class MCLConfig:
    inflation: float = 2.0
    threads: int = 8
    mcl_exe: str = "mcl"


def run_mcl(abc_path: str | Path, out_raw: str | Path, cfg: MCLConfig) -> dict:
    abc_path = Path(abc_path)
    out_raw = Path(out_raw)
    out_raw.parent.mkdir(parents=True, exist_ok=True)

    exe = shutil.which(cfg.mcl_exe)
    if exe is None:
        raise FileNotFoundError("mcl not found on PATH. Install it (micromamba install -c conda-forge -c bioconda mcl).")

    cmd = [
        exe, str(abc_path),
        "--abc",
        "-I", str(cfg.inflation),
        "-te", str(cfg.threads),
        "-o", str(out_raw),
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "mcl failed.\n"
            f"Command: {' '.join(cmd)}\n\nSTDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}\n"
        )

    return {"out_raw": str(out_raw), "inflation": cfg.inflation, "threads": cfg.threads}


def clusters_raw_to_tsv(clusters_raw: str | Path, out_tsv: str | Path, min_size: int = 2) -> dict:
    clusters_raw = Path(clusters_raw)
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    kept = 0
    with clusters_raw.open("r", encoding="utf-8") as fin, out_tsv.open("w", encoding="utf-8") as fout:
        fout.write("cluster_id\tsize\tmembers\n")
        for i, line in enumerate(fin, start=1):
            line = line.strip()
            if not line:
                continue
            total += 1
            members = line.split("\t") if "\t" in line else line.split()
            if len(members) < min_size:
                continue
            kept += 1
            fout.write(f"C{i}\t{len(members)}\t{','.join(members)}\n")

    return {"total_clusters": total, "kept_clusters": kept, "out_tsv": str(out_tsv), "min_size": min_size}
