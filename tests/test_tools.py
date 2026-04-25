from aixbio.tools.codon_tables import (
    best_ecoli_codon,
    split_codons,
    synonymous_alternatives,
    translate_codon,
    translate_dna,
    RARE_CODONS_ECOLI,
)
from aixbio.tools.gc import compute_gc
from aixbio.tools.cai import compute_cai
from aixbio.tools.restriction_sites import (
    find_restriction_sites,
    get_native_enzymes,
    get_recognition_site,
    has_restriction_sites,
)


def test_translate_codon():
    assert translate_codon("ATG") == "M"
    assert translate_codon("TAA") == "*"
    assert translate_codon("CGT") == "R"


def test_best_ecoli_codon():
    codon = best_ecoli_codon("R")
    assert translate_codon(codon) == "R"
    assert codon not in RARE_CODONS_ECOLI


def test_split_codons():
    assert split_codons("ATGCGT") == ["ATG", "CGT"]
    assert split_codons("ATGCG") == ["ATG", "CG"]


def test_translate_dna():
    assert translate_dna("ATGCGT") == "MR"


def test_synonymous_alternatives():
    alts = synonymous_alternatives("CGT")
    assert "CGT" not in alts
    assert all(translate_codon(a) == "R" for a in alts)


def test_gc_content():
    assert compute_gc("GCGC") == 1.0
    assert compute_gc("ATAT") == 0.0
    assert abs(compute_gc("GCATAT") - 1 / 3) < 0.01


def test_cai_perfect():
    aa = "MRK"
    dna = "".join(best_ecoli_codon(a) for a in aa)
    cai = compute_cai(dna)
    assert cai > 0.99


def test_restriction_sites():
    dna = "AAAGGATCCAAA"  # contains BamHI (GGATCC)
    sites = find_restriction_sites(dna, ("BamHI", "XhoI"))
    assert len(sites) == 1
    assert sites[0][0] == "BamHI"
    assert has_restriction_sites(dna, ("BamHI",))
    assert not has_restriction_sites("AAAAAA", ("BamHI",))


def test_rna_fold():
    from aixbio.tools.rna_fold import estimate_five_prime_dg
    dg = estimate_five_prime_dg("ATGATGATGATGATG")
    assert isinstance(dg, float)


def test_native_enzymes_ecoli():
    enzymes = get_native_enzymes("Escherichia coli")
    assert len(enzymes) > 0
    assert all(e.startswith("Eco") for e in enzymes)


def test_native_enzymes_unknown_organism():
    enzymes = get_native_enzymes("Unknown organism")
    assert enzymes == ()


def test_recognition_site():
    assert get_recognition_site("EcoRI") == "GAATTC"
    assert get_recognition_site("BamHI") == "GGATCC"
    assert get_recognition_site("NotI") == "GCGGCCGC"
