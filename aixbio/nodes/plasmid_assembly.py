from __future__ import annotations

from aixbio.models.plasmid import PlasmidChain
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.tools.genbank import build_plasmid_record
from aixbio.tools.restriction_sites import get_recognition_site


def plasmid_assembly(state: ChainSubgraphState) -> dict:
    cassette = state["cassette"]
    vector = state["vector"]
    cloning_sites = state["cloning_sites"]

    flank_5 = get_recognition_site(cloning_sites[0]) if cloning_sites else ""
    flank_3 = get_recognition_site(cloning_sites[1]) if len(cloning_sites) > 1 else ""

    flanked_dna = flank_5 + cassette.full_dna + flank_3

    genbank_str, insert_size = build_plasmid_record(
        chain_id=cassette.id,
        cassette_dna=flanked_dna,
        vector=vector,
        cloning_sites=cloning_sites,
    )

    plasmid = PlasmidChain(
        id=cassette.id,
        genbank_file=genbank_str,
        vector=vector,
        insert_size=insert_size,
    )
    return {"plasmid": plasmid}
