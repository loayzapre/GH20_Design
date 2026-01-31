"""
gh20_design

Pipeline for:
CAZy → curated FASTA → filtering → SSN → clustering → consensus design

During development, keep __init__ minimal.
Import submodules explicitly, e.g.:

  from gh20_design.acquire import resolve_sequences
  from gh20_design.database import length_filter_fasta
  from gh20_design.ssn import run_all_vs_all
  from gh20_design.clustering import run_mcl
"""

__all__ = []
