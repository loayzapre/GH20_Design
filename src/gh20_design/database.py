# src/gh20_design/database.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Dict, Tuple

import pandas as pd


@dataclass
class LengthFilterConfig:
    """Filter out extreme sequence lengths using percentiles."""
    low_q: float = 0.05
    high_q: float = 0.95
    min_len: Optional[int] = None  # optional hard cutoff
    max_len: Optional[int] = None  # optional hard cutoff


@dataclass
class DatabaseBuildConfig:
    family: str = "GH20"
    length_filter: LengthFilterConfig = LengthFilterConfig()
    out_dir: Path = Path("outputs/database")
    cache_dir: Path = Path("outputs/cache")
    add_taxonomy: bool = True
    # CD-HIT identity cutoffs you mentioned (100, 95, 90); keep last as downstream
    cdhit_cutoffs: Tuple[float, ...] = (1.00, 0.95, 0.90)


def clean_species_name(name: str) -> str:
    """Normalize species names (keep your current rules)."""
    # TODO: paste your rules from the notebook version
    return " ".join(str(name).strip().split())


def is_valid_protein(seq: str) -> bool:
    """Basic protein sanity checks."""
    if not isinstance(seq, str):
        return False
    seq = seq.strip().upper()
    if len(seq) == 0:
        return False
    # keep only typical AA + X (unknown); tune if you used other rules
    allowed = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ")
    return all(c in allowed for c in seq)


def filter_length(df: pd.DataFrame, cfg: LengthFilterConfig) -> pd.DataFrame:
    """Percentile-based length filter + optional hard min/max."""
    if "length" not in df.columns:
        raise ValueError("Expected a 'length' column.")

    low = int(df["length"].quantile(cfg.low_q))
    high = int(df["length"].quantile(cfg.high_q))

    lo = cfg.min_len if cfg.min_len is not None else low
    hi = cfg.max_len if cfg.max_len is not None else high

    return df[(df["length"] >= lo) & (df["length"] <= hi)].copy()


def fasta_from_dataframe(df: pd.DataFrame, out_fasta: Path, id_col: str = "id", seq_col: str = "sequence") -> None:
    """Write FASTA from a dataframe."""
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            f.write(f">{row[id_col]}\n{row[seq_col]}\n")


# ---------- HIGH-LEVEL PIPELINE ----------

def build_gh20_database(
    identifiers_df: pd.DataFrame,
    cfg: DatabaseBuildConfig,
) -> pd.DataFrame:
    """
    Build curated GH20 database dataframe and write key outputs.

    Expected identifiers_df columns depend on your CAZy scraping step, but typically:
      - 'id' (unique protein id / accession)
      - 'sequence' (if already fetched) OR fields to fetch sequences
      - optional taxonomy fields
    """
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    cfg.cache_dir.mkdir(parents=True, exist_ok=True)

    df = identifiers_df.copy()

    # If your pipeline fills sequences later, call your fill_sequence function here.
    # df = fill_sequence_indatabase(df, ...)

    # Clean + validate
    if "sequence" not in df.columns:
        raise ValueError("identifiers_df must contain a 'sequence' column at this stage.")

    df["sequence"] = df["sequence"].astype(str).str.strip()
    df = df[df["sequence"].map(is_valid_protein)].copy()
    df["length"] = df["sequence"].str.len()

    # Length filter (percentiles like your Chapter III description) :contentReference[oaicite:1]{index=1}
    df = filter_length(df, cfg.length_filter)

    # Optional: taxonomy step (keep your own implementation; call it here)
    if cfg.add_taxonomy:
        # df = add_taxonomy_to_db(df, cache_dir=cfg.cache_dir)
        pass

    # Save curated table + FASTA
    out_tsv = cfg.out_dir / f"{cfg.family}_curated.tsv"
    out_fasta = cfg.out_dir / f"{cfg.family}_curated.fasta"
    df.to_csv(out_tsv, sep="\t", index=False)
    fasta_from_dataframe(df, out_fasta, id_col="id", seq_col="sequence")

    return df
