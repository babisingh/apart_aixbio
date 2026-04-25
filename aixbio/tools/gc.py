from __future__ import annotations

from Bio.SeqUtils import gc_fraction


def compute_gc(dna: str) -> float:
    if not dna:
        return 0.0
    return gc_fraction(dna)
