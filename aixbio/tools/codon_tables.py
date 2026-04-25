from __future__ import annotations

from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Forward table: codon → amino acid
# ---------------------------------------------------------------------------

_FWD = standard_dna_table.forward_table  # excludes stops
_STOPS = set(standard_dna_table.stop_codons)


def translate_codon(codon: str) -> str:
    """Translate a single DNA codon to its amino acid (or ``*`` for stop)."""
    codon = codon.upper()
    if codon in _STOPS:
        return "*"
    return _FWD[codon]


def translate_dna(dna: str) -> str:
    """Translate a full DNA coding sequence to a protein string.

    Raises ``ValueError`` when the length is not divisible by 3 (frameshift).
    Stop codons appear as ``*``.
    """
    if len(dna) % 3 != 0:
        raise ValueError(
            f"DNA length {len(dna)} is not divisible by 3 — possible frameshift. "
            f"Trailing nucleotides: {dna[-(len(dna) % 3):]!r}"
        )
    return str(Seq(dna.upper()).translate())


# ---------------------------------------------------------------------------
# Reverse table: amino acid → list of synonymous codons
# ---------------------------------------------------------------------------

AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in _FWD.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)


def synonymous_alternatives(codon: str) -> list[str]:
    """Return all synonymous codons for *codon*, excluding itself."""
    aa = translate_codon(codon)
    return [c for c in AA_TO_CODONS.get(aa, []) if c != codon.upper()]


# ---------------------------------------------------------------------------
# E. coli codon usage (organism-specific)
# ---------------------------------------------------------------------------

# Rare codons in E. coli K-12 — organism-specific domain knowledge with no
# Biopython equivalent.  Must be DNA codons (T, not U).
RARE_CODONS_ECOLI = frozenset({"AGG", "AGA", "CGA", "CTA", "ATA", "CCC"})


def get_ecoli_table() -> dict[str, dict[str, float]]:
    """Return E. coli K-12 codon-usage frequencies.

    Sourced from the embedded table in ``aixbio.tools.cai``.
    """
    from aixbio.tools.cai import ECOLI_CODON_USAGE

    return ECOLI_CODON_USAGE


def best_ecoli_codon(aa: str) -> str:
    """Return the optimal E. coli codon for *aa* (highest usage frequency)."""
    from aixbio.tools.cai import ECOLI_CODON_USAGE

    codons = ECOLI_CODON_USAGE.get(aa)
    if not codons:
        raise ValueError(f"Unknown amino acid: {aa}")
    return max(codons, key=lambda c: codons[c])


def split_codons(dna: str) -> list[str]:
    """Split a DNA string into a list of triplets."""
    return [dna[i : i + 3] for i in range(0, len(dna), 3)]
