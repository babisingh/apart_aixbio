from __future__ import annotations

import logging
import os
import tempfile

import httpx
from Bio.PDB import MMCIFParser, Superimposer

from aixbio.models.structure import StructureResult

logger = logging.getLogger(__name__)

ALPHAFOLD_API = "https://alphafold.com/api"

_TIMEOUT = 60
_MAX_RETRIES = 3
_BACKOFF_BASE = 1.0


def _extract_ca_atoms(structure) -> list:
    return [
        atom
        for model in structure
        for chain in model
        for residue in chain
        if residue.id[0] == " " and "CA" in residue
        for atom in [residue["CA"]]
    ]


def _download_reference_cif(pdb_id: str, dest_dir: str) -> str | None:
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    dest_path = os.path.join(dest_dir, f"{pdb_id.upper()}.cif")
    try:
        with httpx.Client(timeout=_TIMEOUT) as client:
            resp = client.get(url)
            resp.raise_for_status()
            with open(dest_path, "wb") as f:
                f.write(resp.content)
        return dest_path
    except (httpx.HTTPStatusError, httpx.ConnectError, httpx.TimeoutException) as exc:
        logger.warning(f"Could not download reference PDB {pdb_id}: {exc}")
        return None


def _compute_rmsd(
    alphafold_cif_path: str,
    reference_pdb_id: str,
    af_entry_id: str,
) -> float | None:
    try:
        tmpdir = tempfile.mkdtemp(prefix="pdb_ref_")
        ref_path = _download_reference_cif(reference_pdb_id, tmpdir)
        if ref_path is None:
            return None

        parser = MMCIFParser(QUIET=True)
        af_struct = parser.get_structure(af_entry_id, alphafold_cif_path)
        ref_struct = parser.get_structure(reference_pdb_id, ref_path)

        af_ca = _extract_ca_atoms(af_struct)
        ref_ca = _extract_ca_atoms(ref_struct)

        n = min(len(af_ca), len(ref_ca))
        if n < 3:
            logger.warning(
                f"Too few CA atoms to superimpose ({len(af_ca)} vs {len(ref_ca)})"
            )
            return None

        sup = Superimposer()
        sup.set_atoms(ref_ca[:n], af_ca[:n])
        return round(float(sup.rms), 3)
    except Exception:
        logger.warning(
            f"RMSD computation failed for {af_entry_id} vs {reference_pdb_id}",
            exc_info=True,
        )
        return None


async def _fetch_predictions(uniprot_id: str) -> list[dict]:
    import asyncio

    url = f"{ALPHAFOLD_API}/prediction/{uniprot_id}"
    last_error: Exception | None = None
    for attempt in range(_MAX_RETRIES):
        try:
            async with httpx.AsyncClient(timeout=_TIMEOUT) as client:
                resp = await client.get(url)
                resp.raise_for_status()
                return resp.json()
        except (httpx.TimeoutException, httpx.ConnectError, httpx.HTTPStatusError) as exc:
            last_error = exc
            if isinstance(exc, httpx.HTTPStatusError) and exc.response.status_code < 500:
                raise
            wait = _BACKOFF_BASE * (2**attempt)
            logger.warning(
                f"AlphaFold API attempt {attempt + 1}/{_MAX_RETRIES} failed for "
                f"{uniprot_id}: {exc}. Retrying in {wait:.1f}s..."
            )
            await asyncio.sleep(wait)

    raise RuntimeError(
        f"Failed to fetch AlphaFold predictions for {uniprot_id} "
        f"after {_MAX_RETRIES} attempts"
    ) from last_error


async def _download_cif(cif_url: str, dest_path: str) -> None:
    async with httpx.AsyncClient(timeout=_TIMEOUT) as client:
        resp = await client.get(cif_url)
        resp.raise_for_status()
        with open(dest_path, "wb") as f:
            f.write(resp.content)


async def predict_structure(
    chain_id: str,
    uniprot_id: str,
    reference_pdb: str | None = None,
) -> StructureResult:
    """Look up an AlphaFold predicted structure from the AlphaFold DB.

    Uses the AlphaFold Protein Structure Database API at alphafold.com
    (Biopython-compatible endpoints) to retrieve pre-computed predictions
    and optionally compute RMSD against a reference PDB.
    """
    predictions = await _fetch_predictions(uniprot_id)

    if not predictions:
        logger.warning(f"No AlphaFold predictions found for {uniprot_id}")
        return StructureResult(
            id=chain_id,
            plddt_mean=0.0,
            rmsd_to_ref=None,
            perplexity=None,
            structure_file="",
        )

    pred = predictions[0]
    plddt_mean = pred.get("globalMetricValue", 0.0)
    cif_url = pred.get("cifUrl", "")
    entry_id = pred.get("entryId", f"AF-{uniprot_id}-F1")

    cif_path = ""
    rmsd = None

    if cif_url:
        output_dir = os.path.join("output", "structures")
        os.makedirs(output_dir, exist_ok=True)
        cif_path = os.path.join(output_dir, f"{entry_id}.cif")

        if not os.path.exists(cif_path):
            await _download_cif(cif_url, cif_path)
            logger.info(f"Downloaded AlphaFold structure to {cif_path}")

        if reference_pdb:
            rmsd = _compute_rmsd(cif_path, reference_pdb, entry_id)

    logger.info(
        f"AlphaFold prediction for {uniprot_id} chain '{chain_id}': "
        f"pLDDT={plddt_mean:.1f}, RMSD={rmsd}"
    )

    return StructureResult(
        id=chain_id,
        plddt_mean=plddt_mean,
        rmsd_to_ref=rmsd,
        perplexity=None,
        structure_file=cif_path,
    )
