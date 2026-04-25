from __future__ import annotations

import asyncio

from aixbio.models.structure import StructureReport
from aixbio.state.pipeline_state import PipelineState
from aixbio.tools.alphafold import predict_structure


REFERENCE_PDBS = {
    "Insulin": "4INS",
    "hGH": "1HGU",
    "EPO": "1BUY",
    "Interferon": "1AU1",
    "Chymosin": "4CMS",
    "tPA": "1A5H",
}


def structural_validation(state: PipelineState) -> dict:
    protein = state["protein_record"]
    if protein is None:
        return {"structure_report": None}

    ref_pdb = None
    for name, pdb_id in REFERENCE_PDBS.items():
        if name.lower() in protein.name.lower():
            ref_pdb = pdb_id
            break

    results = []
    for chain in protein.chains:
        result = asyncio.run(predict_structure(
            chain_id=chain.id,
            uniprot_id=protein.uniprot_id,
            reference_pdb=ref_pdb,
        ))
        results.append(result)

    return {"structure_report": StructureReport(chains=tuple(results))}
