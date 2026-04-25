# Handover: Migrate remaining utilities to Biopython

> **Context**: `compute_cai` in `aixbio/tools/cai.py` was already migrated to
> use `Bio.SeqUtils.CodonAdaptationIndex`. This document covers the remaining
> hand-rolled functions that have Biopython equivalents.

---

## 1. `codon_tables.py` — Translation helpers

### 1a. `CODON_TO_AA` dict + `translate_codon(codon)`

**Current**: Hand-maintained 64-entry dict mapping DNA triplets → amino acids.

**Biopython replacement**: `Bio.Data.CodonTable.standard_dna_table.forward_table`
(plus `standard_dna_table.stop_codons` for `*`).

```python
from Bio.Data.CodonTable import standard_dna_table

def translate_codon(codon: str) -> str:
    codon = codon.upper()
    if codon in standard_dna_table.stop_codons:
        return "*"
    return standard_dna_table.forward_table[codon]
```

**Callers** (8 files):
`codon_optimization.py`, `escalation_agent.py`, `remediation_agent.py`,
`cai.py` (internal), `codon_tables.py` (internal), `test_tools.py`,
`test_escalation.py`, `test_deterministic_nodes.py`

---

### 1b. `AA_TO_CODONS` dict + `synonymous_alternatives(codon)`

**Current**: Reverse mapping built at import time from `CODON_TO_AA`.

**Biopython replacement**: Build from `standard_dna_table.forward_table`.

```python
from Bio.Data.CodonTable import standard_dna_table

AA_TO_CODONS: dict[str, list[str]] = {}
for codon, aa in standard_dna_table.forward_table.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

def synonymous_alternatives(codon: str) -> list[str]:
    aa = translate_codon(codon)
    return [c for c in AA_TO_CODONS[aa] if c != codon.upper()]
```

**Callers**: `codon_optimization.py`, `test_tools.py`, `test_escalation.py`

---

### 1c. `translate_dna(dna)`

**Current**: Splits into codons, maps each with `translate_codon`.

**Biopython replacement**: `Bio.Seq.Seq.translate()`.

```python
from Bio.Seq import Seq

def translate_dna(dna: str) -> str:
    if len(dna) % 3 != 0:
        raise ValueError(
            f"DNA length {len(dna)} is not divisible by 3 — possible frameshift. "
            f"Trailing nucleotides: {dna[-(len(dna) % 3):]!r}"
        )
    return str(Seq(dna.upper()).translate())
```

> [!NOTE]
> `Seq.translate()` includes the stop codon as `*` by default. If callers
> expect no trailing `*`, use `.translate(to_stop=True)` or `.rstrip("*")`.

**Callers**: `sequence_validation.py`, `test_deterministic_nodes.py`

---

### 1d. `best_ecoli_codon(aa)` + `get_ecoli_table()`

**Current**: Looks up the most-frequent codon for an amino acid from the
`python_codon_tables` E. coli table.

**Biopython replacement**: Use `CodonAdaptationIndex.optimize()` for full
sequence optimization, or keep the per-codon lookup but source data from
Biopython's codon adaptation index (already built in `cai.py`).

```python
# Option A — keep per-codon function, source from the cached CAI index
from aixbio.tools.cai import _get_ecoli_cai_index
from Bio.Data.CodonTable import standard_dna_table

def best_ecoli_codon(aa: str) -> str:
    cai = _get_ecoli_cai_index()
    candidates = [
        c for c, a in standard_dna_table.forward_table.items() if a == aa
    ]
    if not candidates:
        raise ValueError(f"Unknown amino acid: {aa}")
    return max(candidates, key=lambda c: cai[c])
```

```python
# Option B — use CAI.optimize() for whole-sequence optimization
from aixbio.tools.cai import _get_ecoli_cai_index

def optimize_sequence(dna: str) -> str:
    cai = _get_ecoli_cai_index()
    return str(cai.optimize(dna))
```

> [!IMPORTANT]
> `best_ecoli_codon` is called per-codon in inner loops of `remediation_agent.py`
> and `codon_optimization.py`. Option A preserves that usage pattern. Option B
> would require restructuring those callers to do whole-sequence optimization.

**Callers**: `codon_optimization.py`, `remediation_agent.py`, `test_tools.py`,
`test_escalation.py`

---

### 1e. `split_codons(dna)`

**No Biopython equivalent** — this is a trivial utility
(`[dna[i:i+3] for i in range(0, len(dna), 3)]`). Keep as-is.

---

### 1f. `RARE_CODONS_ECOLI`

**No Biopython equivalent** — organism-specific domain knowledge.
Keep as-is.

---

## 2. `gc.py` — `compute_gc(dna)`

**Current**: Hand-counts G+C characters.

**Biopython replacement**: `Bio.SeqUtils.gc_fraction()`.

```python
from Bio.SeqUtils import gc_fraction

def compute_gc(dna: str) -> float:
    if not dna:
        return 0.0
    return gc_fraction(dna)
```

**Callers** (4 files): `sequence_validation.py`, `codon_optimization.py`,
`remediation_agent.py`, `test_tools.py`, `test_escalation.py`

---

## 3. `restriction_sites.py` — `find_restriction_sites()` / `has_restriction_sites()`

**Current**: Hand-maintained `ENZYME_SITES` dict + regex search.

**Biopython replacement**: `Bio.Restriction` module (REBASE-complete enzyme
database with `RestrictionBatch` and `Analysis`).

```python
from Bio.Restriction import RestrictionBatch, Analysis, AllEnzymes
from Bio.Seq import Seq

# Map enzyme name strings to Bio.Restriction enzyme objects
def _get_enzymes(names: list[str]) -> RestrictionBatch:
    rb = RestrictionBatch()
    for name in names:
        rb.add(AllEnzymes.get(name))
    return rb

def find_restriction_sites(
    dna: str, enzymes: tuple[str, ...] | list[str]
) -> list[tuple[str, int]]:
    seq = Seq(dna.upper())
    rb = _get_enzymes(list(enzymes))
    ana = Analysis(rb, seq, linear=True)
    hits = []
    for enzyme, positions in ana.full().items():
        for pos in positions:
            hits.append((str(enzyme), pos - 1))  # Biopython uses 1-indexed
    return hits
```

> [!WARNING]
> **Breaking change potential**: The current `ENZYME_SITES` dict is also
> imported directly by `plasmid_assembly.py`, `codon_optimization.py`, and
> `remediation_agent.py` to get the raw recognition sequence string and its
> length. The Biopython `Restriction` enzyme objects expose `.site` for this,
> e.g. `EcoRI.site` → `"GAATTC"`. Callers using `ENZYME_SITES[enzyme]` must be
> updated to use `getattr(Restriction, enzyme).site`.

**Callers** (5 files): `sequence_validation.py`, `codon_optimization.py`,
`remediation_agent.py`, `plasmid_assembly.py`, `test_tools.py`

---

## 4. Dependency to remove: `python-codon-tables`

After migrating `get_ecoli_table()` and `best_ecoli_codon()`, the
`python_codon_tables` package is no longer needed. It is currently only
imported in `codon_tables.py` line 3. Remove from `pyproject.toml` after
migration.

---

## Suggested migration order

| Step | File(s) | Function(s) | Risk |
|------|---------|-------------|------|
| 1 | `gc.py` | `compute_gc` | 🟢 Low — drop-in, identical semantics |
| 2 | `codon_tables.py` | `translate_codon`, `translate_dna`, `CODON_TO_AA` | 🟢 Low — standard genetic code |
| 3 | `codon_tables.py` | `AA_TO_CODONS`, `synonymous_alternatives` | 🟢 Low — derived from step 2 |
| 4 | `codon_tables.py` | `best_ecoli_codon`, `get_ecoli_table` | 🟡 Medium — need to decide Option A vs B |
| 5 | `restriction_sites.py` | `find_restriction_sites`, `ENZYME_SITES` | 🟡 Medium — `ENZYME_SITES` used for raw site strings |
| 6 | `pyproject.toml` | Remove `python-codon-tables` dep | 🟢 Low — after step 4 |

### After each step
- Run `uv run python -m pytest tests/ -v` to verify no regressions.
- Check that all callers listed above still function correctly.
