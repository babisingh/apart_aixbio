from __future__ import annotations

import logging
import os
from io import StringIO
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

PET28A_BACKBONE_SIZE = 5369

_PET28A_FEATURES = [
    ("regulatory", 370, 386, {"note": ["T7 promoter"], "regulatory_class": ["promoter"]}),
    ("regulatory", 386, 404, {"note": ["lac operator"]}),
    ("regulatory", 404, 411, {"note": ["ribosome binding site"]}),
    ("CDS", 2970, 3785, {"gene": ["kan"], "note": ["Kanamycin resistance (KanR)"]}),
    ("rep_origin", 4338, 4927, {"note": ["pBR322 origin of replication"]}),
    ("rep_origin", 4990, 5446, {"note": ["f1 origin of replication"]}),
    ("terminator", 450, 497, {"note": ["T7 terminator"]}),
]

_BACKBONE_SEARCH_PATHS = [
    Path(__file__).parent.parent.parent / "data" / "pET-28a.fasta",
    Path(__file__).parent.parent.parent / "data" / "pET-28a.fa",
    Path(__file__).parent.parent.parent / "data" / "pET-28a.gb",
]

_cached_backbone: str | None = None


def _load_backbone() -> str | None:
    global _cached_backbone
    if _cached_backbone is not None:
        return _cached_backbone

    env_path = os.getenv("PET28A_BACKBONE_PATH")
    search = [Path(env_path)] + _BACKBONE_SEARCH_PATHS if env_path else _BACKBONE_SEARCH_PATHS

    for path in search:
        if not path.exists():
            continue
        text = path.read_text()
        if path.suffix in (".fasta", ".fa"):
            lines = [l.strip() for l in text.splitlines() if l.strip() and not l.startswith(">")]
            seq = "".join(lines).upper()
        elif path.suffix == ".gb":
            from Bio import SeqIO
            record = SeqIO.read(StringIO(text), "genbank")
            seq = str(record.seq).upper()
        else:
            continue
        if len(seq) > 0:
            logger.info(f"Loaded pET-28a(+) backbone from {path} ({len(seq)} bp)")
            _cached_backbone = seq
            return seq

    return None


def build_plasmid_record(
    chain_id: str,
    cassette_dna: str,
    vector: str = "pET-28a(+)",
    cloning_sites: tuple[str, ...] = ("BamHI", "XhoI"),
) -> tuple[str, int]:
    insert_size = len(cassette_dna)

    backbone_seq = _load_backbone()
    if backbone_seq:
        backbone_size = len(backbone_seq)
        plasmid_seq = backbone_seq + cassette_dna
        description_note = ""
    else:
        backbone_size = PET28A_BACKBONE_SIZE
        logger.warning(
            "No pET-28a(+) backbone file found. Using placeholder N's. "
            "Place the real sequence at data/pET-28a.fasta or set PET28A_BACKBONE_PATH."
        )
        plasmid_seq = "N" * backbone_size + cassette_dna
        description_note = " (WARNING: backbone is placeholder N's)"

    total_size = backbone_size + insert_size
    record = SeqRecord(
        Seq(plasmid_seq),
        id=f"{vector}_{chain_id}",
        name=f"{vector}_{chain_id}",
        description=f"{vector} with {chain_id} insert{description_note}",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )

    for feat_type, start, end, qualifiers in _PET28A_FEATURES:
        if start < backbone_size:
            record.features.append(SeqFeature(
                FeatureLocation(start, min(end, backbone_size)),
                type=feat_type,
                qualifiers=qualifiers,
            ))

    record.features.append(SeqFeature(
        FeatureLocation(0, backbone_size),
        type="source",
        qualifiers={"note": [f"{vector} backbone"]},
    ))

    record.features.append(SeqFeature(
        FeatureLocation(backbone_size, total_size),
        type="CDS",
        qualifiers={
            "gene": [chain_id],
            "note": [f"Optimized {chain_id} expression cassette"],
            "codon_start": [1],
        },
    ))

    buf = StringIO()
    from Bio import SeqIO
    SeqIO.write(record, buf, "genbank")
    return buf.getvalue(), insert_size
