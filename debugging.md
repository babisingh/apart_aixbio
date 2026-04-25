# Adversarial Review: aixbio biocompound pipeline

The code is a **hybrid leaning toward demo-ware**: two of six steps are real engineering, the rest are LLM scaffolding, stubs, or fake data wearing valid file formats. Architecture is also dramatically over-spec'd vs handoff.md.

## Architecture: spec rejected DAGs; code uses LangGraph

handoff.md:472 explicitly says **"Why not a DAG? The pipeline is mostly linear"** and prescribes a ~50-line functional framework with `pipeline()`/`retry()`/`fork()`. Implementation is `aixbio/graph/main_graph.py:18-72` — a full LangGraph StateGraph with subgraphs, conditional_edges, `Send()` fan-out, and human-in-the-loop checkpoints (`human_checkpoint_chains`, `human_checkpoint_plasmid`) **not in the spec at all**. Functionally equivalent, ~10× the code, harder to debug.

## Step-by-step science verdict

| Step | Verdict | Evidence |
|---|---|---|
| 1. Sequence retrieval | 🔴 **LLM hallucination risk** | `nodes/sequence_retrieval.py:19-93` hands UniProt JSON to an OpenRouter LLM and asks it to extract mature chain coordinates. On parse failure (lines 51-61) it **falls back to the full precursor** — for insulin that means signal peptide + C-peptide + A + B as one blob. No deterministic feature-table parser. No insulin A+B extraction test. |
| 2. Codon optimization | 🟢 **Real** | `tools/codon_tables.py:61-66` uses real `python-codon-tables`; `tools/cai.py:8-31` is a correct Sharp & Li CAI. One small bug in `nodes/codon_optimization.py:35` — the post-swap window is `pos + site_len * 2` instead of `pos + site_len`, which can mask whether the site was actually removed. |
| 3. Cassette assembly | 🔴 **6xHis frameshift bug** | `nodes/cassette_assembly.py:12` defines `"6xHis": "CACCACCACCACCACCACC"` — **19 nt**, one trailing `C` too many (correct is 18 nt = `CAC`×6). The comment on line 11 even claims "CAC×6 = 18 nt" but the literal disagrees. Cassette frame = `ATG (3) + tag (19) + Enterokinase (15) = 37 nt` before the gene → `37 mod 3 = 1` → **every emitted cassette reads the gene in frame +1, encoding a nonfunctional protein**. Step 5 back-translation will now raise `ValueError` (`tools/codon_tables.py:52-56`) instead of silently mistranslating, so the pipeline halts rather than producing wrong-but-plausible output. Glycosylation N-X-S/T scan is otherwise solid. **One-char fix: drop the trailing `C`.** |
| 4. Plasmid assembly | 🟡 **Valid GenBank, fake backbone** | `tools/genbank.py:45` literally does `plasmid_seq = "N" * 5369 + cassette_dna`. The pET-28a(+) backbone is 5369 N's. The file parses, the features annotate, but **no lab can synthesize this**. The test only asserts "LOCUS" appears in the output. |
| 5. Validation | 🟡 **Checks real, remediation LLM-driven** | GC/CAI/restriction/back-translation/rare-codon are all deterministic and correct. RNA secondary structure (`tools/rna_fold.py`) is a Nussinov heuristic with hardcoded −1.5 kcal/pair — file itself admits this isn't production. **Remediation loop** (`nodes/remediation_agent.py:16-139`) asks the LLM to suggest codon swaps, then applies them — non-deterministic, can introduce new restriction sites, no test coverage. |
| 6. Structural validation | 🔴 **Stub** | `tools/alphafold.py:20-27` self-documents as STUBBED, returns `plddt_mean=0.0`. |

## Other significant issues

- **No file output.** Spec deliverables list FASTA/GenBank/JSON files. There is **no `open()`/`write()`** in the pipeline nodes — `PlasmidChain.genbank_file` holds the string content, not a path. `__main__.py` runs the graph in memory and exits.
- **UniProt error handling.** `tools/uniprot.py` retries are fine, but if the entry lacks a `sequence` key (or returns 404), downstream `extract_sequence()` KeyErrors with no user-facing message.
- **Test coverage gaps.** `test_deterministic_nodes.py` hardcodes the insulin B-chain AA string, so the LLM-driven Step 1 is **never tested**. No remediation-agent test. No end-to-end test from `compound_id="P01308"`.
- **Data models match spec.** Frozen dataclasses across `models/` are faithful — one of the few clean wins.

## Hottest flaws (in priority order)

1. **6xHis tag is 19 nt, frameshifts the gene** (`cassette_assembly.py:12`) — every cassette currently emitted is biologically broken. Trip-wire: Step 5 back-translation now raises rather than silently mistranslating. One-character fix.
2. **Step 1 is an LLM black box** with a fallback to the full precursor — the silent-corruption mode is the worst kind.
3. **Plasmid backbone is 5369 N's** — output is structurally valid and scientifically useless.
4. **Step 5 remediation uses an LLM** instead of deterministic hill-climbing on GC/CAI; can suggest fixes that re-introduce sites it just removed.
5. **No artifacts written to disk** — the spec's deliverables list isn't actually delivered.
6. **LangGraph over-engineering** — human checkpoints and dual-graph compilation for what the spec described as a linear chain with one retry boundary.

## Bottom line

Steps 2 and 3 are publishable-quality. Steps 1, 4, 5-remediation, and 6 are scaffolding that *looks* like science. For a hackathon demo this is fine; for anything claiming "production-ready DNA" (handoff.md:5) it would lose money. The cheapest fixes with the largest payoff: replace Step 1's LLM with a 30-line deterministic UniProt feature parser, drop a real pET-28a(+) FASTA into the repo, and write the JSON/GenBank artifacts to disk in `__main__.py`.
