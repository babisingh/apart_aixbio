from __future__ import annotations

from Bio.Data.CodonTable import standard_dna_table
from Bio.SeqUtils import CodonAdaptationIndex

# ---------------------------------------------------------------------------
# E. coli K-12 codon-usage frequencies (source: Kazusa / REBASE).
# Embedded directly to avoid the ``python-codon-tables`` dependency.
# ---------------------------------------------------------------------------

ECOLI_CODON_USAGE: dict[str, dict[str, float]] = {
    "*": {"TAA": 0.64, "TAG": 0.07, "TGA": 0.29},
    "A": {"GCA": 0.21, "GCC": 0.27, "GCG": 0.36, "GCT": 0.16},
    "C": {"TGC": 0.56, "TGT": 0.44},
    "D": {"GAC": 0.37, "GAT": 0.63},
    "E": {"GAA": 0.69, "GAG": 0.31},
    "F": {"TTC": 0.43, "TTT": 0.57},
    "G": {"GGA": 0.11, "GGC": 0.41, "GGG": 0.15, "GGT": 0.34},
    "H": {"CAC": 0.43, "CAT": 0.57},
    "I": {"ATA": 0.07, "ATC": 0.42, "ATT": 0.51},
    "K": {"AAA": 0.76, "AAG": 0.24},
    "L": {"CTA": 0.04, "CTC": 0.10, "CTG": 0.50, "CTT": 0.10, "TTA": 0.13, "TTG": 0.13},
    "M": {"ATG": 1.00},
    "N": {"AAC": 0.55, "AAT": 0.45},
    "P": {"CCA": 0.19, "CCC": 0.12, "CCG": 0.53, "CCT": 0.16},
    "Q": {"CAA": 0.35, "CAG": 0.65},
    "R": {"AGA": 0.04, "AGG": 0.02, "CGA": 0.06, "CGC": 0.40, "CGG": 0.10, "CGT": 0.38},
    "S": {"AGC": 0.28, "AGT": 0.15, "TCA": 0.12, "TCC": 0.15, "TCG": 0.15, "TCT": 0.15},
    "T": {"ACA": 0.13, "ACC": 0.44, "ACG": 0.27, "ACT": 0.16},
    "V": {"GTA": 0.15, "GTC": 0.22, "GTG": 0.37, "GTT": 0.26},
    "W": {"TGG": 1.00},
    "Y": {"TAC": 0.43, "TAT": 0.57},
}

# ---------------------------------------------------------------------------
# CAI index (cached singleton)
# ---------------------------------------------------------------------------

_CAI_INDEX: CodonAdaptationIndex | None = None


def _get_ecoli_cai_index() -> CodonAdaptationIndex:
    """Build a Biopython CAI index from the E. coli codon usage table.

    The relative adaptiveness weights (w_ij) are computed as
    freq(codon) / max_freq(synonymous codons), matching Sharp & Li 1987.
    The index is cached after first construction.
    """
    global _CAI_INDEX
    if _CAI_INDEX is not None:
        return _CAI_INDEX

    weights: dict[str, float] = {}
    for aa, codon_freqs in ECOLI_CODON_USAGE.items():
        max_freq = max(codon_freqs.values())
        for codon, freq in codon_freqs.items():
            weights[codon] = freq / max_freq if max_freq > 0 else 0.5

    cai = CodonAdaptationIndex.__new__(CodonAdaptationIndex)
    dict.__init__(cai)
    cai._table = standard_dna_table
    cai.update(weights)

    _CAI_INDEX = cai
    return _CAI_INDEX


def compute_cai(dna: str) -> float:
    """Compute the Codon Adaptation Index for E. coli using Biopython.

    Follows Sharp & Li (Nucleic Acids Res. 1987) — single-codon amino
    acids (Met, Trp) and stop codons are excluded from the geometric mean.
    """
    dna = dna.upper()
    if len(dna) < 3:
        return 0.0

    cai_index = _get_ecoli_cai_index()
    try:
        return cai_index.calculate(dna)
    except (KeyError, TypeError, ZeroDivisionError):
        return 0.0
