from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple, Optional


# -------------------------
# FASTA I/O
# -------------------------

def iter_fasta(path: str | Path) -> Iterable[Tuple[str, str]]:
    path = Path(path)
    h = None
    seq = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq)
                h = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if h is not None:
            yield h, "".join(seq)


# -------------------------
# Cluster + consensus indexes
# -------------------------

def load_member_to_cluster(clusters_tsv: str | Path) -> Dict[str, str]:
    """
    clusters TSV must have columns: cluster_id, members
    where members is comma-separated.
    Returns: member_id -> cluster_id
    """
    clusters_tsv = Path(clusters_tsv)
    m2c: Dict[str, str] = {}

    with clusters_tsv.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        col = {name: i for i, name in enumerate(header)}
        if "cluster_id" not in col or "members" not in col:
            raise ValueError("clusters TSV must have columns: cluster_id, members")

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            cid = parts[col["cluster_id"]]
            members = parts[col["members"]].split(",") if parts[col["members"]] else []
            for m in members:
                if m:
                    m2c[m] = cid

    return m2c


def load_cluster_to_consensus(consensus_all_fasta: str | Path) -> Dict[str, str]:
    """
    Expects headers like: C12|consensus
    Returns: cluster_id -> consensus_sequence
    """
    consensus_all_fasta = Path(consensus_all_fasta)
    c2seq: Dict[str, str] = {}
    for h, s in iter_fasta(consensus_all_fasta):
        cid = h.split("|", 1)[0].strip()
        c2seq[cid] = s
    return c2seq


def _prefix_match(member_to_cluster: Dict[str, str], query: str) -> Optional[str]:
    qtok = query.split()[0]
    for m, cid in member_to_cluster.items():
        if m.split()[0] == qtok:
            return cid
    return None


# -------------------------
# Main query API
# -------------------------

def consensus_for_protein(
    protein_id: str,
    clusters_tsv: str | Path,
    consensus_all_fasta: str | Path,
    *,
    allow_prefix_match: bool = True,
) -> dict:
    """
    Given a protein id (node id used in SSN), returns its cluster and consensus.
    """
    m2c = load_member_to_cluster(clusters_tsv)
    c2seq = load_cluster_to_consensus(consensus_all_fasta)

    cid = m2c.get(protein_id)
    if cid is None and allow_prefix_match:
        cid = _prefix_match(m2c, protein_id)

    if cid is None:
        return {"found": False, "query": protein_id, "reason": "protein not found in clusters TSV"}

    cons = c2seq.get(cid)
    if cons is None:
        return {"found": False, "query": protein_id, "cluster_id": cid, "reason": "consensus not found"}

    return {
        "found": True,
        "query": protein_id,
        "cluster_id": cid,
        "consensus_len": len(cons),
        "consensus_seq": cons,
    }
