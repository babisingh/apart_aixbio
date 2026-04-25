from __future__ import annotations

import logging
from datetime import datetime, timezone

from aixbio.models.audit import AgentDecision
from aixbio.models.dna import DNAChain
from aixbio.models.remediation import PlannedFix, RemediationAction, RemediationPlan
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.tools.cai import compute_cai
from aixbio.tools.codon_tables import (
    RARE_CODONS_ECOLI,
    best_ecoli_codon,
    get_ecoli_table,
    split_codons,
    synonymous_alternatives,
    translate_codon,
)
from aixbio.tools.gc import compute_gc
from aixbio.tools.restriction_sites import find_restriction_sites, get_recognition_site, has_restriction_sites

logger = logging.getLogger(__name__)


def remediation_agent(state: ChainSubgraphState) -> dict:
    dna_chain = state["optimized_dna"]
    failed_checks = state["failed_checks"]
    attempt = state["remediation_attempt"]

    codons = split_codons(dna_chain.dna_sequence)
    all_fixes: list[PlannedFix] = []
    check_names = {c.name for c in failed_checks}

    if "restriction_sites" in check_names:
        all_fixes.extend(_plan_restriction_site_fixes(codons, state["avoid_sites"]))

    if "rare_codons" in check_names:
        all_fixes.extend(_plan_rare_codon_fixes(codons))

    if "cai_score" in check_names:
        all_fixes.extend(_plan_cai_fixes(codons))

    if "gc_content" in check_names:
        gc_check = next(c for c in failed_checks if c.name == "gc_content")
        gc_val = float(gc_check.value)
        all_fixes.extend(_plan_gc_fixes(codons, gc_val))

    plan = RemediationPlan(
        actions=tuple(all_fixes),
        reasoning=f"Deterministic fixes for: {', '.join(sorted(check_names))}",
        priority_order=tuple(sorted(check_names)),
    )

    decision = AgentDecision(
        node="remediation_agent",
        reasoning=plan.reasoning,
        action=f"Planned {len(plan.actions)} fixes for {len(failed_checks)} failed checks",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_summary=f"Chain {dna_chain.id}, attempt {attempt + 1}, {len(failed_checks)} failures",
        output_summary=f"Fixes: {', '.join(a.check_name for a in plan.actions)}" if plan.actions else "No fixes available",
    )

    return {
        "remediation_plan": plan,
        "decision_log": [decision],
    }


def _plan_restriction_site_fixes(
    codons: list[str], avoid_sites: tuple[str, ...]
) -> list[PlannedFix]:
    dna = "".join(codons)
    hits = find_restriction_sites(dna, avoid_sites)
    fixes = []
    for enzyme, pos in hits:
        site_len = len(get_recognition_site(enzyme))
        codon_start = pos // 3
        codon_end = (pos + site_len - 1) // 3 + 1
        fixed = False
        for ci in range(codon_start, min(codon_end, len(codons))):
            for alt in synonymous_alternatives(codons[ci]):
                old = codons[ci]
                codons[ci] = alt
                local_dna = "".join(codons)[max(0, pos - site_len):pos + site_len * 2]
                if not has_restriction_sites(local_dna, [enzyme]):
                    fixes.append(PlannedFix(
                        check_name="restriction_sites",
                        strategy="synonymous_codon_swap",
                        target_positions=(ci,),
                        replacement_codons=(alt,),
                    ))
                    fixed = True
                    break
                codons[ci] = old
            if fixed:
                break
    return fixes


def _plan_rare_codon_fixes(codons: list[str]) -> list[PlannedFix]:
    fixes = []
    for i, codon in enumerate(codons):
        if codon.upper() in RARE_CODONS_ECOLI:
            best = best_ecoli_codon(translate_codon(codon))
            if best != codon.upper():
                fixes.append(PlannedFix(
                    check_name="rare_codons",
                    strategy="replace_rare_with_optimal",
                    target_positions=(i,),
                    replacement_codons=(best,),
                ))
                codons[i] = best
    return fixes


def _plan_cai_fixes(codons: list[str]) -> list[PlannedFix]:
    table = get_ecoli_table()
    scored = []
    for i, codon in enumerate(codons):
        aa = translate_codon(codon)
        if aa == "*":
            continue
        best = best_ecoli_codon(aa)
        if best == codon.upper():
            continue
        aa_codons = table.get(aa, {})
        max_freq = max(aa_codons.values()) if aa_codons else 1.0
        codon_freq = aa_codons.get(codon.upper(), 0.0)
        relative = codon_freq / max_freq if max_freq > 0 else 0.0
        if relative < 0.5:
            scored.append((i, best, relative))

    scored.sort(key=lambda x: x[2])
    fixes = []
    for i, best, _ in scored:
        fixes.append(PlannedFix(
            check_name="cai_score",
            strategy="replace_suboptimal_codon",
            target_positions=(i,),
            replacement_codons=(best,),
        ))
        codons[i] = best
    return fixes


def _plan_gc_fixes(codons: list[str], current_gc: float) -> list[PlannedFix]:
    need_higher = current_gc < 0.50
    fixes = []

    candidates = []
    for i, codon in enumerate(codons):
        aa = translate_codon(codon)
        if aa == "*":
            continue
        alts = synonymous_alternatives(codon)
        if not alts:
            continue
        codon_gc = compute_gc(codon)
        if need_higher:
            best_alt = max(alts, key=lambda c: compute_gc(c))
            improvement = compute_gc(best_alt) - codon_gc
        else:
            best_alt = min(alts, key=lambda c: compute_gc(c))
            improvement = codon_gc - compute_gc(best_alt)
        if improvement > 0:
            candidates.append((i, best_alt, improvement))

    candidates.sort(key=lambda x: -x[2])

    for i, alt, _ in candidates:
        fixes.append(PlannedFix(
            check_name="gc_content",
            strategy="increase_gc" if need_higher else "decrease_gc",
            target_positions=(i,),
            replacement_codons=(alt,),
        ))
        codons[i] = alt
        new_gc = compute_gc("".join(codons))
        if 0.50 <= new_gc <= 0.60:
            break

    return fixes


def apply_fixes(state: ChainSubgraphState) -> dict:
    dna_chain = state["optimized_dna"]
    plan = state["remediation_plan"]
    codons = split_codons(dna_chain.dna_sequence)

    actions: list[RemediationAction] = []

    for fix in plan.actions:
        for pos, new_codon in zip(fix.target_positions, fix.replacement_codons):
            if pos >= len(codons):
                continue
            old_codon = codons[pos]
            old_aa = translate_codon(old_codon)
            new_aa = translate_codon(new_codon.upper())
            if old_aa != new_aa:
                continue

            actions.append(RemediationAction(
                check_name=fix.check_name,
                fix_type=fix.strategy,
                positions_affected=(pos,),
                codons_before=(old_codon,),
                codons_after=(new_codon.upper(),),
                reasoning=plan.reasoning,
            ))
            codons[pos] = new_codon.upper()

    new_dna = "".join(codons)
    new_dna_chain = DNAChain(
        id=dna_chain.id,
        dna_sequence=new_dna,
        cai_score=compute_cai(new_dna),
        gc_content=compute_gc(new_dna),
    )

    return {
        "optimized_dna": new_dna_chain,
        "remediation_history": actions,
        "remediation_attempt": state["remediation_attempt"] + 1,
    }
