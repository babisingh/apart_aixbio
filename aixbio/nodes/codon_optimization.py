from __future__ import annotations

from aixbio.models.dna import DNAChain
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.tools.cai import compute_cai
from aixbio.tools.codon_tables import best_ecoli_codon, split_codons, synonymous_alternatives, translate_codon
from aixbio.tools.gc import compute_gc
from aixbio.tools.restriction_sites import find_restriction_sites, get_recognition_site, has_restriction_sites


def codon_optimization(state: ChainSubgraphState) -> dict:
    chain = state["chain"]
    avoid = state["avoid_sites"]

    codons = [best_ecoli_codon(aa) for aa in chain.aa_sequence]
    dna = "".join(codons)

    for _pass in range(5):
        hits = find_restriction_sites(dna, avoid)
        if not hits:
            break
        for enzyme, pos in hits:
            site_len = len(get_recognition_site(enzyme))
            codon_start = pos // 3
            codon_end = (pos + site_len - 1) // 3 + 1
            site_removed = False
            for ci in range(codon_start, min(codon_end, len(codons))):
                original_codon = codons[ci]
                alts = synonymous_alternatives(codons[ci])
                for alt in alts:
                    codons[ci] = alt
                    test_dna = "".join(codons)
                    # Check if this substitution eliminated the specific site
                    local_dna = test_dna[max(0, pos - site_len):pos + site_len * 2]
                    if not has_restriction_sites(local_dna, [enzyme]):
                        site_removed = True
                        break
                if site_removed:
                    break
                codons[ci] = original_codon  # revert if no alt worked
        dna = "".join(codons)

    cai_score = compute_cai(dna)
    gc_content = compute_gc(dna)

    dna_chain = DNAChain(
        id=chain.id,
        dna_sequence=dna,
        cai_score=cai_score,
        gc_content=gc_content,
    )
    return {"optimized_dna": dna_chain}
