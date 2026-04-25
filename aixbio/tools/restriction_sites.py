from __future__ import annotations

from functools import lru_cache

from Bio.Restriction import AllEnzymes, RestrictionBatch
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# REBASE enzyme-name prefix convention:
#   Eco* → Escherichia coli, Bsu* → Bacillus subtilis, Hin* → Haemophilus
#   influenzae, etc.  This mapping lets the pipeline derive native restriction
#   enzymes from the host organism without hardcoding enzyme lists.
# ---------------------------------------------------------------------------

ORGANISM_ENZYME_PREFIX: dict[str, str] = {
    "Escherichia coli": "Eco",
    "Bacillus subtilis": "Bsu",
    "Haemophilus influenzae": "Hin",
    "Staphylococcus aureus": "Sau",
    "Streptomyces albus": "Sal",
    "Thermus aquaticus": "Taq",
    "Deinococcus radiodurans": "Dra",
}


@lru_cache(maxsize=16)
def get_native_enzymes(host_organism: str) -> tuple[str, ...]:
    """Return enzyme names native to *host_organism*.

    Matches using the REBASE naming prefix convention (e.g.
    ``"Escherichia coli"`` → all ``Eco*`` enzymes).
    """
    prefix = ORGANISM_ENZYME_PREFIX.get(host_organism)
    if prefix is None:
        return ()
    return tuple(
        str(e) for e in AllEnzymes if str(e).startswith(prefix)
    )


def get_recognition_site(enzyme_name: str) -> str:
    """Look up the recognition sequence for a named enzyme.

    Uses the Biopython ``Bio.Restriction`` database (REBASE-complete).
    """
    from Bio import Restriction

    enzyme = getattr(Restriction, enzyme_name, None)
    if enzyme is None:
        raise ValueError(f"Unknown enzyme: {enzyme_name}")
    return enzyme.site


def find_restriction_sites(
    dna: str, enzymes: tuple[str, ...] | list[str],
) -> list[tuple[str, int]]:
    """Find all recognition-site occurrences for the given enzymes.

    Returns ``(enzyme_name, position)`` tuples where *position* is the
    **0-based start** of the recognition site (matching the legacy API).
    """
    if not enzymes:
        return []

    seq = Seq(dna.upper())
    batch = RestrictionBatch(first=[], suppliers=[])
    for name in enzymes:
        from Bio import Restriction

        enzyme = getattr(Restriction, name, None)
        if enzyme is not None:
            batch.add(enzyme)

    results = batch.search(seq, linear=True)
    hits: list[tuple[str, int]] = []
    for enzyme, positions in results.items():
        for cut_pos in positions:
            # Convert Biopython 1-based cut position to 0-based site start
            # If the cleavage site is undefined (fst5 is None), Biopython returns 
            # the 1-based start of the recognition site, so replacing it with 0 works.
            fst5 = enzyme.fst5 if enzyme.fst5 is not None else 0
            site_start = (cut_pos - 1) - fst5
            hits.append((str(enzyme), site_start))
    return hits


def has_restriction_sites(
    dna: str, enzymes: tuple[str, ...] | list[str],
) -> bool:
    """Return ``True`` if *dna* contains any of the given restriction sites."""
    return len(find_restriction_sites(dna, enzymes)) > 0
