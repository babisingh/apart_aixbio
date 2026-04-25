# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

aixbio is a digital-to-biological pipeline that takes a UniProt protein accession (e.g. `P01308` for insulin), reconstructs a production-ready DNA sequence optimized for E. coli expression, assembles it into an expression cassette and plasmid, and validates the result. Built as a hackathon project using LangGraph.

## Commands

```bash
# Install dependencies (uses uv, Python 3.11)
uv sync

# Run pipeline
uv run aixbio P01308                          # insulin, with human checkpoints
uv run aixbio P01308 --auto-approve           # skip human review prompts
uv run aixbio P01308 --structural             # include AlphaFold step (stubbed)
uv run aixbio P01308 --escalation             # enable LLM escalation on remediation failure

# Run tests
uv run pytest                                 # all tests
uv run pytest tests/test_deterministic_nodes.py  # steps 2-5 without LLM
uv run pytest tests/test_tools.py             # unit tests for bio tools
uv run pytest tests/test_chain_subgraph.py    # LangGraph subgraph integration
uv run pytest tests/test_escalation.py        # escalation agent (mocked LLM)
```

## Environment

Requires `OPENROUTER_API_KEY` in `.env`. The LLM model defaults to `deepseek/deepseek-v4-flash` via OpenRouter (`LLM_MODEL` env var). Step 1 (sequence retrieval) always makes an LLM call. The escalation agent (`--escalation` flag) makes one LLM call per chain only when the deterministic remediation loop exhausts its retry budget. All other steps are deterministic.

## Architecture

The pipeline is a LangGraph StateGraph with two levels:

**Main graph** (`graph/main_graph.py`): sequence_retrieval -> human_checkpoint -> fan-out to per-chain processing -> merge_results -> human_checkpoint -> optional structural_validation

**Chain subgraph** (`graph/chain_subgraph.py`): codon_optimization -> cassette_assembly -> plasmid_assembly -> sequence_validation. On validation failure, routes to a deterministic remediation loop that applies synonymous codon swaps and revalidates, up to `max_remediation_attempts`. If remediation is exhausted and `--escalation` is enabled, an LLM escalation agent fires once and chooses one of four outcomes: `apply_plan` (try specific swaps the greedy loop missed), `incompatible` (host cannot express this protein), `change_strategy` (change tag/protease/vector/cloning sites and re-run), or `give_up` (pipeline bug or unfixable). The `escalation_used` flag prevents re-entry.

### Key directories

- `models/` - Frozen dataclasses for inter-step data (ProteinRecord, DNAChain, CassetteChain, PlasmidChain, ValidationReport, EscalationDecision, etc.). New escalation models are a discriminated union on the `kind` field; routers branch on `kind`, never on free text.
- `nodes/` - LangGraph node functions. Each takes and returns a state dict. Deterministic nodes (codon_optimization, cassette_assembly, plasmid_assembly, sequence_validation, remediation_agent) vs LLM-backed nodes (sequence_retrieval, escalation_agent).
- `tools/` - Pure utility functions: codon tables, CAI calculation, GC content, restriction site scanning, RNA fold estimation, GenBank file generation, UniProt/AlphaFold API clients.
- `state/` - TypedDict state definitions. `PipelineState` for the main graph, `ChainSubgraphState` for per-chain processing. Uses LangGraph `Annotated` reducers for list accumulation fields.
- `prompts/` - LLM prompt templates for sequence retrieval and remediation.

### State flow

Nodes return partial state dicts (only the keys they update). List fields (`chain_results`, `decision_log`, `warnings`, `remediation_history`) use an `append_log` reducer -- nodes must return only the **delta**, not the full accumulated list.

### Config

`config.py` defines pipeline defaults (host organism, restriction sites, tag type, vector) and wraps `ChatOpenAI` as `ChatOpenRouter` to fix a `max_tokens` vs `max_completion_tokens` incompatibility with OpenRouter.
