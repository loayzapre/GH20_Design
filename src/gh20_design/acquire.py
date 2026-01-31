from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple
import re
from datetime import datetime


# =========================
# 1) TXT parsing (generic)
# =========================

@dataclass(frozen=True)
class ProteinIdRecord:
    """
    A single record parsed from a tab-separated text file.

    Expected format (>= 5 tab-separated columns):
        family   kingdom   organism   protein_id   source

    Example:
        GH20    Bacteria   Abiotrophia defectiva ATCC 49176   WWT19304.1   ncbi
        GH18    Eukaryota  Some fungus v1.0                   423730      jgi
    """
    family: str
    kingdom: str
    organism: str
    protein_id: str
    source: str  # e.g. "ncbi", "jgi" (we keep whatever is provided)


def load_records_from_txt(txt_path: str | Path) -> List[ProteinIdRecord]:
    """
    Load records from a tab-separated TXT file.

    - Skips empty lines
    - Skips an optional header line starting with 'family'
    - Requires >= 5 columns
    """
    txt_path = Path(txt_path)
    records: List[ProteinIdRecord] = []

    with txt_path.open("r", encoding="utf-8") as f:
        for ln, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue

            if line.lower().startswith("family") and "\t" in line:
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                raise ValueError(
                    f"[load_records_from_txt] Line {ln} has < 5 columns: {line!r}"
                )

            family, kingdom, organism, pid, source = parts[:5]
            records.append(
                ProteinIdRecord(
                    family=family.strip(),
                    kingdom=kingdom.strip(),
                    organism=organism.strip(),
                    protein_id=pid.strip(),
                    source=source.strip().lower(),
                )
            )

    return records


# =========================
# 2) FASTA reading + indexing
# =========================

# Split header tokens by whitespace or pipe
_HEADER_SPLIT_RE = re.compile(r"[\s\|]+")

# Broad NCBI-like accession pattern, e.g. BCM93776.1, WWT19304.1, WP_012345678.1
_ACCESSION_RE = re.compile(r"\b([A-Z]{2,6}_?\d{3,12}\.\d+)\b")

# Numeric-only tokens often represent JGI protein IDs in exports (e.g. 423730)
_DIGITS_ONLY_RE = re.compile(r"^\d{4,14}$")


def iter_fasta(fasta_path: str | Path) -> Iterator[Tuple[str, str]]:
    """
    Yield (header, sequence) from a FASTA file.
    Header is returned without the leading '>'.
    """
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


def _tokens_from_header(header: str) -> List[str]:
    """
    Extract multiple useful lookup tokens from a FASTA header.

    We index:
    - first token
    - all tokens split by whitespace and '|'
    - accession-like substrings detected by regex
    - pure digit tokens (useful for numeric IDs like JGI)
    """
    header = header.strip()
    tokens = [t for t in _HEADER_SPLIT_RE.split(header) if t]

    collected: List[str] = []
    if tokens:
        collected.append(tokens[0])

    collected.extend(tokens)

    for m in _ACCESSION_RE.finditer(header):
        collected.append(m.group(1))

    for t in tokens:
        if _DIGITS_ONLY_RE.match(t):
            collected.append(t)

    # Deduplicate while preserving order
    seen: set[str] = set()
    uniq: List[str] = []
    for t in collected:
        if t not in seen:
            seen.add(t)
            uniq.append(t)

    return uniq


def build_fasta_index(fasta_path: str | Path) -> Dict[str, Tuple[str, str]]:
    """
    Build a token -> (full_header, sequence) index for a FASTA DB.

    If multiple sequences share the same token, we keep the first occurrence.
    """
    index: Dict[str, Tuple[str, str]] = {}
    for header, seq in iter_fasta(fasta_path):
        for tok in _tokens_from_header(header):
            if tok not in index:
                index[tok] = (header, seq)
    return index


# =========================
# 3) Matching + reporting
# =========================

@dataclass
class MatchReport:
    total: int
    matched: int
    missing: int
    missing_records: List[ProteinIdRecord]


def match_records_to_fasta(
    records: Iterable[ProteinIdRecord],
    fasta_index: Dict[str, Tuple[str, str]],
) -> Tuple[Dict[str, Tuple[str, str, ProteinIdRecord]], MatchReport]:
    """
    Match records to sequences using a pre-built FASTA index.

    Returns:
    - hits: protein_id -> (full_header, sequence, original_record)
    - report: MatchReport
    """
    recs = list(records)
    hits: Dict[str, Tuple[str, str, ProteinIdRecord]] = {}
    missing: List[ProteinIdRecord] = []

    for r in recs:
        pid = r.protein_id.strip()
        if not pid:
            missing.append(r)
            continue

        # 1) Direct lookup
        if pid in fasta_index:
            full_header, seq = fasta_index[pid]
            hits[pid] = (full_header, seq, r)
            continue

        # 2) Common fallback: accession without version (e.g. BCM93776 from BCM93776.1)
        if "." in pid:
            pid_nov = pid.split(".", 1)[0]
            if pid_nov in fasta_index:
                full_header, seq = fasta_index[pid_nov]
                hits[pid] = (full_header, seq, r)
                continue

        # 3) Common fallback: numeric IDs sometimes appear as pid.t1, pid.t2, etc.
        # We try a few typical variants.
        if pid.isdigit():
            for suffix in (".t1", ".t2", ".t3"):
                cand = pid + suffix
                if cand in fasta_index:
                    full_header, seq = fasta_index[cand]
                    hits[pid] = (full_header, seq, r)
                    break
            else:
                missing.append(r)
        else:
            missing.append(r)

    report = MatchReport(
        total=len(recs),
        matched=len(hits),
        missing=len(missing),
        missing_records=missing,
    )
    return hits, report


# =========================
# 4) Writing outputs
# =========================

def write_fasta(
    hits: Dict[str, Tuple[str, str, ProteinIdRecord]],
    out_fasta: str | Path,
    header_style: str = "id_plus_source",
    wrap: int = 80,
) -> None:
    """
    Write matched sequences to a curated FASTA.

    header_style:
      - "keep_db_header": use the original FASTA DB header
      - "id_plus_source": >{protein_id}|{source}
      - "rich": >{protein_id}|{source}|{family}|{kingdom}|{organism}

    wrap:
      - sequence line width; set <=0 to disable wrapping
    """
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    def _wrap_seq(seq: str) -> str:
        if wrap <= 0:
            return seq + "\n"
        return "\n".join(seq[i:i + wrap] for i in range(0, len(seq), wrap)) + "\n"

    with out_fasta.open("w", encoding="utf-8") as f:
        for pid, (full_header, seq, rec) in hits.items():
            if header_style == "keep_db_header":
                header = full_header
            elif header_style == "rich":
                organism = rec.organism.replace("\t", " ").replace("|", "/").strip()
                header = f"{pid}|{rec.source}|{rec.family}|{rec.kingdom}|{organism}"
            else:
                header = f"{pid}|{rec.source}"

            f.write(f">{header}\n")
            f.write(_wrap_seq(seq))


def write_missing_tsv(missing: List[ProteinIdRecord], out_tsv: str | Path) -> None:
    """Write missing records to a TSV file for debugging/auditing."""
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    with out_tsv.open("w", encoding="utf-8") as f:
        f.write("family\tkingdom\torganism\tprotein_id\tsource\n")
        for r in missing:
            f.write(f"{r.family}\t{r.kingdom}\t{r.organism}\t{r.protein_id}\t{r.source}\n")


# =========================
# 5) Convenience: TXT + FASTA -> curated FASTA
# =========================

@dataclass
class ResolveConfig:
    txt_path: str | Path
    db_fasta: str | Path
    out_fasta: str | Path
    out_missing_tsv: Optional[str | Path] = None
    header_style: str = "id_plus_source"
    fetch_missing_from_ncbi: bool = True
    ncbi_batch_size: int = 200
    ncbi_sleep_s: float = 0.34
    ncbi_timeout: int = 30
    ncbi_retries: int = 2


def resolve_sequences(cfg: ResolveConfig) -> MatchReport:
    """
    End-to-end:
      1) parse TXT
      2) build FASTA index from local DB
      3) match records to DB sequences
      4) (optional) fetch missing NCBI records from NCBI and add them
      5) write curated FASTA
      6) optionally write missing TSV
    """
    records = load_records_from_txt(cfg.txt_path)
    fasta_index = build_fasta_index(cfg.db_fasta)

    hits, report = match_records_to_fasta(records, fasta_index)

    # ---- NEW: fetch missing from NCBI (only those labeled as ncbi) ----
    if cfg.fetch_missing_from_ncbi and report.missing_records:
        ncbi_missing = [r for r in report.missing_records if r.source == "ncbi"]

        if ncbi_missing:
            ids = [r.protein_id for r in ncbi_missing]
            fetched = fetch_ncbi_fasta_batch(
                ids,
                timeout=cfg.ncbi_timeout,
                batch_size=cfg.ncbi_batch_size,
                sleep_s=cfg.ncbi_sleep_s,
                retries=cfg.ncbi_retries,
            )

            # Add fetched sequences into hits
            for r in ncbi_missing:
                pid = r.protein_id.strip().split()[0]
                if pid in fetched:
                    # use a minimal "full_header" placeholder
                    full_header = pid
                    hits[r.protein_id] = (full_header, fetched[pid], r)

            # Recompute missing list: keep those still not in hits
            still_missing: List[ProteinIdRecord] = []
            for r in report.missing_records:
                if r.protein_id in hits:
                    continue
                # also consider versionless match
                if "." in r.protein_id and r.protein_id.split(".", 1)[0] in fetched:
                    continue
                still_missing.append(r)

            report = MatchReport(
                total=report.total,
                matched=len(hits),
                missing=len(still_missing),
                missing_records=still_missing,
            )

    # Write outputs
    write_fasta(hits, cfg.out_fasta, header_style=cfg.header_style)

    if cfg.out_missing_tsv is not None:
        write_missing_tsv(report.missing_records, cfg.out_missing_tsv)

    return report

import time
import requests

_NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _parse_fasta_text_to_dict(fasta_text: str) -> Dict[str, str]:
    """
    Parse FASTA returned by NCBI into {id_token -> sequence}.
    We store by first header token (e.g., 'WP_012345678.1').
    """
    out: Dict[str, str] = {}
    header = None
    seq_chunks: List[str] = []

    for raw in fasta_text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                key = header.split()[0]
                out[key] = "".join(seq_chunks)
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line)

    if header is not None:
        key = header.split()[0]
        out[key] = "".join(seq_chunks)

    return out


def fetch_ncbi_fasta_batch(
    protein_ids: List[str],
    *,
    timeout: int = 30,
    batch_size: int = 200,
    sleep_s: float = 0.34,   # ~3 req/s (safe default)
    retries: int = 2,
) -> Dict[str, str]:
    """
    Fetch protein sequences from NCBI in batches via EFetch.

    Returns: {protein_id -> sequence} for IDs successfully fetched.
    """
    results: Dict[str, str] = {}

    # Normalize IDs (first token only)
    ids = [pid.strip().split()[0] for pid in protein_ids if pid and pid.strip()]
    if not ids:
        return results

    for i in range(0, len(ids), batch_size):
        chunk = ids[i:i + batch_size]
        params = {
            "db": "protein",
            "id": ",".join(chunk),
            "rettype": "fasta",
            "retmode": "text",
        }

        last_err = None
        for attempt in range(retries + 1):
            try:
                r = requests.get(_NCBI_EFETCH, params=params, timeout=timeout)
                if r.ok and r.text.startswith(">"):
                    parsed = _parse_fasta_text_to_dict(r.text)
                    # keep only ones we asked for (some may fail silently)
                    for pid in chunk:
                        if pid in parsed and parsed[pid]:
                            results[pid] = parsed[pid]
                    break
                else:
                    last_err = f"HTTP {r.status_code} or empty FASTA"
            except Exception as e:
                last_err = repr(e)

            # backoff
            time.sleep(1.0 + attempt)

        # be gentle to NCBI
        time.sleep(sleep_s)

        # If you want debugging, you can print last_err here (optional)
        # if last_err: print("[NCBI] batch failed:", last_err)

    return results


def fetch_ncbi_fasta_one(protein_id: str, timeout: int = 30) -> Optional[str]:
    """
    Fetch one protein sequence from NCBI EFetch (FASTA).
    Returns the sequence (string) or None.
    """
    pid = protein_id.strip().split()[0]
    if not pid:
        return None

    params = {"db": "protein", "id": pid, "rettype": "fasta", "retmode": "text"}
    r = requests.get(_NCBI_EFETCH, params=params, timeout=timeout)
    if not r.ok or not r.text.startswith(">"):
        return None

    seq = "".join(line.strip() for line in r.text.splitlines() if line and not line.startswith(">"))
    return seq if seq else None


def fasta_has_id(fasta_path: str | Path, protein_id: str) -> bool:
    """
    Check if curated FASTA already contains this protein_id as first token in the header.
    We match headers like: >WWT19304.1|ncbi|...
    """
    fasta_path = Path(fasta_path)
    pid = protein_id.strip().split()[0]
    if not fasta_path.exists():
        return False

    with fasta_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip()
                first = h.split()[0].split("|")[0]
                if first == pid:
                    return True
    return False


def append_fasta_entry_if_new(
    out_fasta: str | Path,
    protein_id: str,
    seq: str,
    *,
    header: str,
    wrap: int = 80,
) -> bool:
    """
    Append a single FASTA entry if protein_id is not already present.
    Returns True if written, False if already existed.
    """
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    if fasta_has_id(out_fasta, protein_id):
        return False

    def _wrap(seq_: str) -> str:
        if wrap <= 0:
            return seq_ + "\n"
        return "\n".join(seq_[i:i + wrap] for i in range(0, len(seq_), wrap)) + "\n"

    with out_fasta.open("a", encoding="utf-8") as f:
        f.write(f">{header}\n")
        f.write(_wrap(seq))

    return True


def append_missing_record(
    missing_tsv: str | Path,
    r: ProteinIdRecord,
) -> None:
    """
    Append one missing record to missing.tsv (creates header if file doesn't exist).
    """
    missing_tsv = Path(missing_tsv)
    missing_tsv.parent.mkdir(parents=True, exist_ok=True)

    write_header = not missing_tsv.exists()
    with missing_tsv.open("a", encoding="utf-8") as f:
        if write_header:
            f.write("family\tkingdom\torganism\tprotein_id\tsource\n")
        f.write(f"{r.family}\t{r.kingdom}\t{r.organism}\t{r.protein_id}\t{r.source}\n")


def append_record_from_sources(
    r: ProteinIdRecord,
    *,
    db_fasta: str | Path,
    curated_fasta: str | Path,
    missing_tsv: str | Path,
    header_style: str = "rich",
    try_ncbi: bool = True,
    ncbi_timeout: int = 30,
) -> dict:
    """
    Resolve ONE ProteinIdRecord and append to curated_fasta if found.

    Strategy:
      - Try local FASTA DB first (fast, offline).
      - If not found and r.source == 'ncbi' and try_ncbi=True: fetch from NCBI.
      - If still not found: append to missing.tsv.

    Returns a dict with status and details.
    """
    db_fasta = Path(db_fasta)
    curated_fasta = Path(curated_fasta)
    missing_tsv = Path(missing_tsv)

    # Already present?
    if fasta_has_id(curated_fasta, r.protein_id):
        return {"status": "already_present", "protein_id": r.protein_id}

    # Local DB lookup (use your existing index logic)
    fasta_index = build_fasta_index(db_fasta)

    seq = None
    full_header = None

    pid = r.protein_id.strip()
    if pid in fasta_index:
        full_header, seq = fasta_index[pid]
        source_used = "local_db"
    else:
        # try versionless
        source_used = None
        if "." in pid:
            pid_nov = pid.split(".", 1)[0]
            if pid_nov in fasta_index:
                full_header, seq = fasta_index[pid_nov]
                source_used = "local_db"

    # NCBI fetch if missing and allowed
    if seq is None and try_ncbi and r.source == "ncbi":
        seq = fetch_ncbi_fasta_one(r.protein_id, timeout=ncbi_timeout)
        if seq:
            full_header = r.protein_id
            source_used = "ncbi"

    if seq is None:
        append_missing_record(missing_tsv, r)
        return {"status": "missing", "protein_id": r.protein_id}

    # Build header in the same style as write_fasta()
    if header_style == "keep_db_header" and full_header:
        header = full_header
    elif header_style == "rich":
        organism = r.organism.replace("\t", " ").replace("|", "/").strip()
        header = f"{r.protein_id}|{r.source}|{r.family}|{r.kingdom}|{organism}"
    else:
        header = f"{r.protein_id}|{r.source}"

    written = append_fasta_entry_if_new(curated_fasta, r.protein_id, seq, header=header, wrap=80)

    return {
        "status": "added" if written else "already_present",
        "protein_id": r.protein_id,
        "source_used": source_used,
        "seq_len": len(seq) if seq else 0,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
