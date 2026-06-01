
# -*- coding: utf-8 -*-
"""Adaptive Simulated Annealing auto tuner."""

import math
import random
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

from .filters import (
    apply_filter_params,
    collect_matched_chr_ids,
    filter_non_dominant_noise_by_dominant_pairs,
    filter_blocks_by_pair_aln,
    filter_blocks_by_pair_count,
    prune_chr_partners_by_block_count,
)

INITIAL_IDENTITY_CEILING = 0.99
REPRESENTATIVE_MIN_ALN_LEN = 1000


def _clip(v, lo, hi, cast=float):
    v = max(lo, min(hi, v))
    return cast(v)


def _dominant_alignment_initial_state(raw_blocks):
    """
    Build data-driven initial filter thresholds from chromosome-dominant pairs.

    For every chromosome on both rows, the dominant partner is the chromosome pair
    with the largest cumulative aligned length. One best alignment block is chosen
    from each dominant pair, and global thresholds are set to the least stringent
    values still preserving these representative dominant relationships.
    """
    pair_stats = {}
    for b in raw_blocks:
        key = (int(b["pair"]), b["upper"], b["lower"])
        stat = pair_stats.setdefault(key, {
            "aln_sum": 0,
            "block_count": 0,
            "blocks": [],
        })
        stat["aln_sum"] += int(b.get("alnLen", 0) or 0)
        stat["block_count"] += 1
        stat["blocks"].append(b)

    if not pair_stats:
        return {
            "state": {
                "min_block": 1000,
                "min_identity": 0.10,
                "min_mapq": 0,
                "min_pair_aln": 0,
                "min_pair_blocks": 0,
                "max_bottom_per_top": 0,
            },
            "summary": {
                "init_method": "dominant_alignment_empty_fallback",
                "init_dominant_pairs": 0,
                "init_representative_blocks": 0,
            },
            "dominant_chr_pair_keys": [],
            "dominant_chr_pair_records": [],
        }

    dominant_keys = set()
    best_upper = {}
    best_lower = {}
    for key, stat in pair_stats.items():
        pair, upper, lower = key
        upper_chr_key = (pair, upper)
        lower_chr_key = (pair + 1, lower)
        if upper_chr_key not in best_upper or stat["aln_sum"] > pair_stats[best_upper[upper_chr_key]]["aln_sum"]:
            best_upper[upper_chr_key] = key
        if lower_chr_key not in best_lower or stat["aln_sum"] > pair_stats[best_lower[lower_chr_key]]["aln_sum"]:
            best_lower[lower_chr_key] = key

    dominant_keys.update(best_upper.values())
    dominant_keys.update(best_lower.values())

    representatives = []
    representative_pair_aln = []
    dominant_records = []
    for key in dominant_keys:
        stat = pair_stats[key]
        representative_candidates = [
            b for b in stat["blocks"]
            if int(b.get("alnLen", 0) or 0) >= REPRESENTATIVE_MIN_ALN_LEN
        ] or stat["blocks"]
        best_block = max(
            representative_candidates,
            key=lambda b: (
                float(b.get("identity", 0.0) or 0.0),
                int(b.get("alnLen", 0) or 0),
                int(b.get("mapq", 0) or 0),
            )
        )
        representatives.append(best_block)
        representative_pair_aln.append(int(best_block.get("alnLen", 0) or 0))
        pair, upper, lower = key
        dominant_records.append({
            "pair": pair,
            "upper": upper,
            "lower": lower,
            "dominant_aln_sum": int(stat.get("aln_sum", 0) or 0),
            "dominant_block_count": int(stat.get("block_count", 0) or 0),
            "representative_upper_start": int(best_block.get("upperStart", 0) or 0),
            "representative_upper_end": int(best_block.get("upperEnd", 0) or 0),
            "representative_lower_start": int(best_block.get("lowerStart", 0) or 0),
            "representative_lower_end": int(best_block.get("lowerEnd", 0) or 0),
            "representative_aln_len": int(best_block.get("alnLen", 0) or 0),
            "representative_identity": float(best_block.get("identity", 0.0) or 0.0),
            "representative_mapq": int(best_block.get("mapq", 0) or 0),
        })

    # These thresholds are inferred from the representative blocks, so the
    # representatives themselves must be able to pass the initial filter.
    min_block = max(1, min(int(b.get("alnLen", 0) or 0) for b in representatives))
    min_identity = min(
        INITIAL_IDENTITY_CEILING,
        max(0.10, min(float(b.get("identity", 0.0) or 0.0) for b in representatives)),
    )
    min_mapq = max(0, min(int(b.get("mapq", 0) or 0) for b in representatives))
    min_pair_aln = max(0, min(int(x) for x in representative_pair_aln))

    state = {
        "min_block": int(min_block),
        "min_identity": float(min_identity),
        "min_mapq": int(min_mapq),
        "min_pair_aln": int(min_pair_aln),
        "min_pair_blocks": 1,
        "max_bottom_per_top": 0,
    }
    summary = {
        "init_method": "chromosome_dominant_alignment",
        "init_dominant_pairs": len(dominant_keys),
        "init_representative_blocks": len(representatives),
        "init_min_block": state["min_block"],
        "init_min_identity": state["min_identity"],
        "init_min_mapq": state["min_mapq"],
        "init_min_pair_aln": state["min_pair_aln"],
        "init_min_pair_blocks": state["min_pair_blocks"],
    }
    return {
        "state": state,
        "summary": summary,
        "dominant_chr_pair_keys": [list(x) for x in sorted(dominant_keys)],
        "dominant_chr_pair_records": sorted(
            dominant_records,
            key=lambda x: (x["pair"], x["upper"], x["lower"]),
        ),
    }


def _apply(raw_blocks, params, max_partners, partner_cap_mode):
    before = apply_filter_params(raw_blocks, params)

    after, dropped, report = prune_chr_partners_by_block_count(
        before,
        max_partners=max_partners,
        mode=partner_cap_mode,
    )

    return before, after, dropped, report


def _chr_fill_metrics(blocks, baseline_chr_ids, target_chr_fill):
    baseline_chr_ids = set(baseline_chr_ids or [])
    fill_sum = defaultdict(int)
    chr_len = {}

    for b in blocks:
        upper_key = (int(b["pair"]), b["upper"])
        lower_key = (int(b["pair"]) + 1, b["lower"])
        upper_span = max(0, int(b.get("upperEnd", 0) or 0) - int(b.get("upperStart", 0) or 0))
        lower_span = max(0, int(b.get("lowerEnd", 0) or 0) - int(b.get("lowerStart", 0) or 0))
        fill_sum[upper_key] += upper_span
        fill_sum[lower_key] += lower_span
        chr_len[upper_key] = max(int(chr_len.get(upper_key, 0) or 0), int(b.get("upperLen", 0) or 0))
        chr_len[lower_key] = max(int(chr_len.get(lower_key, 0) or 0), int(b.get("lowerLen", 0) or 0))

    rates = []
    passed = 0
    for key in baseline_chr_ids:
        rate = min(1.0, fill_sum.get(key, 0) / max(1, int(chr_len.get(key, 0) or 0)))
        rates.append(rate)
        if rate >= target_chr_fill:
            passed += 1

    if not rates:
        return {
            "chr_fill_target": target_chr_fill,
            "chr_fill_min": 0.0,
            "chr_fill_mean": 0.0,
            "chr_fill_pass_fraction": 0.0,
            "chr_fill_pass_count": 0,
            "chr_fill_total_chr": 0,
        }

    return {
        "chr_fill_target": target_chr_fill,
        "chr_fill_min": min(rates),
        "chr_fill_mean": sum(rates) / len(rates),
        "chr_fill_pass_fraction": passed / len(rates),
        "chr_fill_pass_count": passed,
        "chr_fill_total_chr": len(rates),
    }


def _score(
    blocks,
    target_blocks,
    coverage,
    min_coverage,
    params,
    fill_metrics=None,
    target_chr_fill=0.90,
    fill_mode_block_threshold=300,
    noise_ratio=0.0,
):
    count = len(blocks)

    fill_metrics = fill_metrics or {}
    chr_fill_mean = float(fill_metrics.get("chr_fill_mean", 0.0) or 0.0)
    chr_fill_pass_fraction = float(fill_metrics.get("chr_fill_pass_fraction", 0.0) or 0.0)
    chr_fill_min = float(fill_metrics.get("chr_fill_min", 0.0) or 0.0)
    fill_shortage = max(0.0, target_chr_fill - chr_fill_min)
    pass_shortage = max(0.0, min_coverage - chr_fill_pass_fraction)
    block_reward = math.log1p(count)

    # ------------------------------------------------------------------
    # chromosome retention protection
    # prevent ASA from deleting too many chromosomes
    # ------------------------------------------------------------------
    coverage_reward = coverage * 8.0

    coverage_penalty = max(0, min_coverage - coverage) * 12.0

    if coverage < 0.80:
        coverage_penalty += 20.0

    # ------------------------------------------------------------------
    # noise penalty
    # weaker than before to avoid over-cleaning
    # ------------------------------------------------------------------
    aln_lengths = [int(b.get("alnLen", 0) or 0) for b in blocks]
    small_blocks = sum(1 for x in aln_lengths if x < 100000)

    noise_penalty = (
        small_blocks / max(1, count)
    ) * 0.80

    # ------------------------------------------------------------------
    # continuity reward
    # encourage long synteny structures
    # ------------------------------------------------------------------
    total_aln = sum(aln_lengths)
    avg_aln = total_aln / max(1, count)
    continuity_bonus = avg_aln / 200000.0
    long_block_bonus = sum(min(x, 5000000) for x in aln_lengths) / 1000000.0

    # ------------------------------------------------------------------
    # reward chromosome fusion / one-to-many structure
    # ------------------------------------------------------------------
    chr_pairs = set()

    for b in blocks:
        chr_pairs.add((b["upper"], b["lower"]))

    fusion_bonus = len(chr_pairs) * 0.015

    # Blocks that would be removed by the final drawing noise filter should not
    # drive the optimizer. Keep this penalty light so real fusion/one-to-many
    # structures can still survive when the final filter classifies them as kept.
    post_filter_noise_penalty = max(0.0, min(1.0, float(noise_ratio or 0.0))) * 2.0

    # ------------------------------------------------------------------
    # softer complexity penalty
    # ------------------------------------------------------------------
    complexity_penalty = (
        params["min_identity"] * 0.03
        + params["min_mapq"] * 0.003
    )

    score = (
        + block_reward * 3.0
        + chr_fill_mean * 4.0
        + chr_fill_pass_fraction * 3.0
        - fill_shortage * 4.0
        - pass_shortage * 5.0
        + coverage_reward
        - coverage_penalty
        - noise_penalty
        + continuity_bonus
        + long_block_bonus
        + fusion_bonus
        - post_filter_noise_penalty
        - complexity_penalty
    )

    return score


def _score_blocks_after_noise_filter(blocks, dominant_chr_pair_keys, params):
    if not blocks:
        return [], 0.0, {"status": "empty_input"}

    pair_params = {
        int(b["pair"]): params
        for b in blocks
    }
    scored_blocks, _, info = filter_non_dominant_noise_by_dominant_pairs(
        blocks,
        dominant_chr_pair_keys=dominant_chr_pair_keys,
        pair_params=pair_params,
        never_filter_to_zero=True,
    )
    noise_ratio = 1.0 - (len(scored_blocks) / max(1, len(blocks)))
    info = dict(info or {})
    info["noise_ratio"] = noise_ratio
    info["scored_blocks"] = len(scored_blocks)
    info["provisional_blocks"] = len(blocks)
    return scored_blocks, noise_ratio, info


def _block_key(block):
    return (
        int(block.get("pair", -1)),
        str(block.get("upper", "")),
        str(block.get("lower", "")),
        int(block.get("upperStart", 0) or 0),
        int(block.get("upperEnd", 0) or 0),
        int(block.get("lowerStart", 0) or 0),
        int(block.get("lowerEnd", 0) or 0),
        str(block.get("strand", "")),
    )


def _chr_pair_key(value):
    if isinstance(value, (list, tuple)) and len(value) >= 3:
        return (int(value[0]), value[1], value[2])
    if isinstance(value, (list, tuple)) and len(value) >= 2:
        return (value[0], value[1])
    return value


def _block_chr_pair_keys(block):
    pair = int(block.get("pair", -1))
    upper = block.get("upper")
    lower = block.get("lower")
    return (pair, upper, lower), (upper, lower)


def _partition_scored_blocks(blocks, dominant_chr_pair_keys):
    dominant_keys = {
        _chr_pair_key(key)
        for key in (dominant_chr_pair_keys or [])
    }
    dominant_blocks = []
    non_dominant_blocks = []
    for block in blocks:
        pair_key, legacy_pair_key = _block_chr_pair_keys(block)
        if pair_key in dominant_keys or legacy_pair_key in dominant_keys:
            dominant_blocks.append(block)
        else:
            non_dominant_blocks.append(block)
    return dominant_blocks, non_dominant_blocks


def _block_key_set(blocks):
    return {_block_key(block) for block in blocks}


def _alignment_sum(blocks):
    return sum(int(block.get("alnLen", 0) or 0) for block in blocks)


def _gain_metrics(scored_blocks, dominant_chr_pair_keys, previous_dominant_keys=None, previous_noise_keys=None):
    previous_dominant_keys = previous_dominant_keys or set()
    previous_noise_keys = previous_noise_keys or set()

    dominant_blocks, noise_blocks = _partition_scored_blocks(scored_blocks, dominant_chr_pair_keys)
    dominant_keys = _block_key_set(dominant_blocks)
    noise_keys = _block_key_set(noise_blocks)
    new_dominant_keys = dominant_keys - previous_dominant_keys
    new_noise_keys = noise_keys - previous_noise_keys

    new_dominant_blocks = [block for block in dominant_blocks if _block_key(block) in new_dominant_keys]
    new_noise_blocks = [block for block in noise_blocks if _block_key(block) in new_noise_keys]

    return {
        "scored_blocks": len(scored_blocks),
        "previous_scored_blocks": len(previous_dominant_keys) + len(previous_noise_keys),
        "dominant_blocks": dominant_blocks,
        "non_dominant_blocks": noise_blocks,
        "dominant_keys": dominant_keys,
        "noise_keys": noise_keys,
        "dominant_block_count": len(dominant_blocks),
        "non_dominant_block_count": len(noise_blocks),
        "delta_dominant_blocks": len(new_dominant_blocks),
        "delta_noise_blocks": len(new_noise_blocks),
        "delta_dominant_aln": _alignment_sum(new_dominant_blocks),
        "delta_noise_aln": _alignment_sum(new_noise_blocks),
    }


def _should_stop_relaxing(gain_metrics):
    """
    Parameter-free stop rule: stop when the accepted relaxation adds more
    non-dominant kept blocks than dominant backbone blocks, and those non-dominant
    gains also carry more aligned sequence. Blocks already removable by the final
    noise filter are not counted here.
    """
    return (
        int(gain_metrics.get("scored_blocks", 0) or 0)
        >= int(gain_metrics.get("previous_scored_blocks", 0) or 0)
        and
        int(gain_metrics.get("delta_noise_blocks", 0) or 0)
        > int(gain_metrics.get("delta_dominant_blocks", 0) or 0)
        and int(gain_metrics.get("delta_noise_aln", 0) or 0)
        > int(gain_metrics.get("delta_dominant_aln", 0) or 0)
    )


def _candidate_net_gain(gain_metrics):
    return (
        int(gain_metrics.get("delta_dominant_blocks", 0) or 0)
        - int(gain_metrics.get("delta_noise_blocks", 0) or 0)
        + (
            int(gain_metrics.get("delta_dominant_aln", 0) or 0)
            - int(gain_metrics.get("delta_noise_aln", 0) or 0)
        ) / 1000000.0
    )


def _state_snapshot(
    round_id,
    state,
    blocks,
    scored_blocks,
    before,
    dropped,
    report,
    score,
    coverage,
    fill,
    reached,
    noise_ratio,
    dominant_keys,
    noise_keys,
):
    return {
        "round": round_id,
        "state": dict(state),
        "blocks": list(blocks),
        "scored_blocks": list(scored_blocks),
        "before": list(before),
        "dropped": list(dropped),
        "report": report,
        "score": score,
        "coverage": coverage,
        "fill": dict(fill),
        "reached": reached,
        "noise_ratio": noise_ratio,
        "dominant_keys": set(dominant_keys),
        "noise_keys": set(noise_keys),
    }


def _single_param_candidates(state, param_name):
    candidates = []

    def add(value):
        candidate = dict(state)
        candidate[param_name] = value
        candidates.append(candidate)

    if param_name == "min_block":
        current = int(state.get(param_name, 1000) or 1000)
        for factor in (0.50, 0.70, 0.85, 1.15):
            add(max(1, int(current * factor)))
    elif param_name == "min_identity":
        current = float(state.get(param_name, 0.10) or 0.10)
        for step in (-0.050, -0.025, -0.010, 0.010):
            add(_clip(current + step, 0.10, 0.99, float))
    elif param_name == "min_mapq":
        current = int(state.get(param_name, 0) or 0)
        for step in (-10, -5, -2, 2):
            add(_clip(current + step, 0, 60, int))
    elif param_name == "min_pair_aln":
        current = int(state.get(param_name, 0) or 0)
        if current <= 0:
            for value in (0, 1000, 5000):
                add(value)
        else:
            for factor in (0.50, 0.70, 0.85, 1.15):
                add(max(0, int(current * factor)))
    elif param_name == "min_pair_blocks":
        current = int(state.get(param_name, 0) or 0)
        for step in (-2, -1, 1):
            add(_clip(current + step, 0, 200, int))

    unique = []
    seen = set()
    for candidate in candidates:
        key = tuple(sorted(candidate.items()))
        if key not in seen and candidate != state:
            seen.add(key)
            unique.append(candidate)
    return unique


def _random_single_param_candidate(state, param_name):
    candidates = _single_param_candidates(state, param_name)
    if not candidates:
        return None
    return random.choice(candidates)


def _select_rollback_snapshot(snapshots, noise_events):
    if not snapshots:
        return None, None
    if not noise_events:
        return snapshots[-1], None

    worst_event = max(
        noise_events,
        key=lambda event: (
            int(event["gain"].get("delta_noise_blocks", 0) or 0),
            int(event["gain"].get("delta_noise_aln", 0) or 0),
        )
    )
    rollback_round = max(0, int(worst_event.get("round", 0) or 0) - 1)
    rollback_snapshot = snapshots[0]
    for snapshot in snapshots:
        if int(snapshot.get("round", 0) or 0) <= rollback_round:
            rollback_snapshot = snapshot
        else:
            break
    return rollback_snapshot, worst_event


def _perturb(state):
    new_state = dict(state)

    new_state["min_block"] = int(
        max(
            1000,
            state["min_block"] * random.uniform(0.6, 1.4)
        )
    )

    new_state["min_identity"] = _clip(
        state["min_identity"] + random.uniform(-0.05, 0.05),
        0.10,
        0.99,
        float,
    )

    new_state["min_mapq"] = _clip(
        state["min_mapq"] + random.randint(-5, 5),
        0,
        60,
        int,
    )

    new_state["min_pair_aln"] = int(
        max(
            0,
            state["min_pair_aln"] * random.uniform(0.5, 1.5)
        )
    )

    new_state["min_pair_blocks"] = _clip(
        state["min_pair_blocks"] + random.randint(-2, 2),
        0,
        200,
        int,
    )

    return new_state


def _halve_strict_params(state):
    relaxed = dict(state)
    relaxed["min_block"] = max(1, int(int(relaxed.get("min_block", 1000)) * 0.50))
    relaxed["min_identity"] = max(0.10, float(relaxed.get("min_identity", 0.60)) * 0.50)
    relaxed["min_mapq"] = max(0, int(int(relaxed.get("min_mapq", 0)) * 0.50))
    relaxed["min_pair_aln"] = max(0, int(int(relaxed.get("min_pair_aln", 0)) * 0.50))
    relaxed["min_pair_blocks"] = max(0, int(int(relaxed.get("min_pair_blocks", 0)) * 0.50))
    relaxed["max_bottom_per_top"] = int(relaxed.get("max_bottom_per_top", 0) or 0)
    return relaxed


def _format_halving_action(reason, old_state, new_state, before_count, after_count, guard_round=None):
    prefix = f"{reason}_drop_guard"
    if guard_round is not None:
        prefix = f"{prefix}_{guard_round}"
    return (
        f"{prefix}:all_params_halved count {before_count}->{after_count};"
        f"min_block {old_state.get('min_block')}->{new_state.get('min_block')};"
        f"min_identity {old_state.get('min_identity')}->{new_state.get('min_identity')};"
        f"min_mapq {old_state.get('min_mapq')}->{new_state.get('min_mapq')};"
        f"min_pair_aln {old_state.get('min_pair_aln')}->{new_state.get('min_pair_aln')};"
        f"min_pair_blocks {old_state.get('min_pair_blocks')}->{new_state.get('min_pair_blocks')}"
    )


def _relax_empty_state(state):
    """Immediately loosen all strict filters after an empty result."""
    return _halve_strict_params(state)


def _best_single_block_state(raw_blocks):
    best = max(
        raw_blocks,
        key=lambda b: (
            int(b.get("alnLen", 0) or 0),
            float(b.get("identity", 0.0) or 0.0),
            int(b.get("mapq", 0) or 0),
        )
    )
    return {
        "min_block": max(1, int(best.get("alnLen", 0) or 0)),
        "min_identity": min(INITIAL_IDENTITY_CEILING, max(0.10, float(best.get("identity", 0.0) or 0.0))),
        "min_mapq": max(0, int(best.get("mapq", 0) or 0)),
        "min_pair_aln": 0,
        "min_pair_blocks": 0,
        "max_bottom_per_top": 0,
    }


def _guard_steep_drop_state(raw_blocks, state, max_drop_fraction=0.50, max_guard_rounds=6):
    """
    If a single parameter removes more than max_drop_fraction of current blocks,
    immediately halve all strict parameters and restart the step check.
    """
    guarded = dict(state)
    actions = []

    def too_steep(before_count, after_count):
        return before_count > 0 and after_count < before_count * (1.0 - max_drop_fraction)

    for guard_round in range(1, max_guard_rounds + 1):
        current = list(raw_blocks)
        triggered = False

        filter_steps = (
            ("min_block", lambda bs: [b for b in bs if int(b.get("alnLen", 0) or 0) >= int(guarded.get("min_block", 0) or 0)]),
            ("min_identity", lambda bs: [b for b in bs if float(b.get("identity", 0.0) or 0.0) >= float(guarded.get("min_identity", 0.0) or 0.0)]),
            ("min_mapq", lambda bs: [b for b in bs if int(b.get("mapq", 0) or 0) >= int(guarded.get("min_mapq", 0) or 0)]),
            ("min_pair_aln", lambda bs: filter_blocks_by_pair_aln(bs, int(guarded.get("min_pair_aln", 0) or 0))),
            ("min_pair_blocks", lambda bs: filter_blocks_by_pair_count(bs, int(guarded.get("min_pair_blocks", 0) or 0))),
        )

        for param_name, filter_func in filter_steps:
            before_count = len(current)
            after = filter_func(current)
            if too_steep(before_count, len(after)):
                old_state = dict(guarded)
                guarded = _halve_strict_params(guarded)
                actions.append(
                    _format_halving_action(param_name, old_state, guarded, before_count, len(after), guard_round)
                )
                triggered = True
                break
            current = after

        if not triggered:
            break

    return guarded, ";".join(actions)


def _apply_with_empty_rescue(raw_blocks, params, max_partners, partner_cap_mode, max_rescue_rounds=8):
    """
    Apply filters and never let a non-empty PAF collapse into an empty seed state.

    The dominant-pair initialization is supposed to provide a chromosome-scale
    backbone. If the inferred thresholds are still too strict, relax them before
    ASA scoring so the optimizer can expand from at least one real alignment.
    """
    state = dict(params)
    before, blocks, dropped, report = _apply(raw_blocks, state, max_partners, partner_cap_mode)
    rescue_actions = []
    rescue_rounds = 0

    while raw_blocks and not blocks and rescue_rounds < max_rescue_rounds:
        old_state = dict(state)
        state = _relax_empty_state(state)
        rescue_rounds += 1
        before, blocks, dropped, report = _apply(raw_blocks, state, max_partners, partner_cap_mode)
        rescue_actions.append(
            _format_halving_action("empty_result", old_state, state, 0, len(blocks), rescue_rounds)
        )

    if raw_blocks and not blocks:
        old_state = dict(state)
        state = _best_single_block_state(raw_blocks)
        before, blocks, dropped, report = _apply(raw_blocks, state, max_partners, partner_cap_mode)
        rescue_rounds += 1
        rescue_actions.append(
            "empty_result_single_block_fallback:"
            f"min_block {old_state.get('min_block')}->{state.get('min_block')};"
            f"min_identity {old_state.get('min_identity')}->{state.get('min_identity')};"
            f"min_mapq {old_state.get('min_mapq')}->{state.get('min_mapq')};"
            f"min_pair_aln {old_state.get('min_pair_aln')}->{state.get('min_pair_aln')};"
            f"min_pair_blocks {old_state.get('min_pair_blocks')}->{state.get('min_pair_blocks')}"
        )

    return before, blocks, dropped, report, state, rescue_rounds, ";".join(rescue_actions), ""


def auto_tune_paf_filters(raw_blocks, args):

    target_blocks = int(getattr(args, "target_blocks", 1000))
    max_partners = int(getattr(args, "max_partners_per_chr", 4))
    partner_cap_mode = getattr(args, "partner_cap_mode", "both")
    min_coverage = float(getattr(args, "auto_min_coverage", 0.65))
    target_chr_fill = float(getattr(args, "target_chr_fill", 0.90))
    fill_mode_block_threshold = int(getattr(args, "fill_mode_block_threshold", 300))

    asa_rounds = int(getattr(args, "asa_rounds", 200))
    start_temp = float(getattr(args, "asa_start_temp", 5.0))
    cooling = float(getattr(args, "asa_cooling", 0.97))

    dominant_init = _dominant_alignment_initial_state(raw_blocks)
    initial_seed_state = dict(dominant_init["state"])
    initial_seed_state["min_pair_blocks"] = max(
        int(initial_seed_state.get("min_pair_blocks", 1)),
        int(getattr(args, "min_pair_blocks", 0)),
    )
    initial_summary = dict(dominant_init["summary"])
    dominant_chr_pair_keys = list(dominant_init.get("dominant_chr_pair_keys", []))
    dominant_chr_pair_records = list(dominant_init.get("dominant_chr_pair_records", []))

    baseline = apply_filter_params(raw_blocks, {
        "min_block": 0,
        "min_identity": 0,
        "min_mapq": 0,
        "min_pair_aln": 0,
        "min_pair_blocks": 0,
        "max_bottom_per_top": 0,
    })

    baseline_chr = collect_matched_chr_ids(baseline)

    history = []

    before, current_blocks, dropped, report, initial_state, initial_rescue_rounds, initial_rescue_action, initial_guard_action = _apply_with_empty_rescue(
        raw_blocks,
        initial_seed_state,
        max_partners,
        partner_cap_mode,
    )

    current_scored_blocks, current_noise_ratio, current_noise_info = _score_blocks_after_noise_filter(
        current_blocks,
        dominant_chr_pair_keys,
        initial_state,
    )
    current_gain = _gain_metrics(current_scored_blocks, dominant_chr_pair_keys)
    current_dominant_keys = current_gain["dominant_keys"]
    current_noise_keys = current_gain["noise_keys"]
    current_cov = len(collect_matched_chr_ids(current_scored_blocks)) / max(1, len(baseline_chr))
    current_fill = _chr_fill_metrics(current_scored_blocks, baseline_chr, target_chr_fill)
    current_reached = bool(current_scored_blocks)

    if not current_scored_blocks:
        current_score = -1e18
        current_reached = False
    else:
        current_score = _score(
            current_scored_blocks,
            target_blocks,
            current_cov,
            min_coverage,
            initial_state,
            fill_metrics=current_fill,
            target_chr_fill=target_chr_fill,
            fill_mode_block_threshold=fill_mode_block_threshold,
            noise_ratio=current_noise_ratio,
        )

    best_state = dict(initial_state)
    best_blocks = current_blocks
    best_scored_blocks = current_scored_blocks
    best_score = current_score
    best_before = before
    best_drop = dropped
    best_report = report
    best_cov = current_cov
    best_fill = current_fill
    best_reached = current_reached
    best_noise_ratio = current_noise_ratio
    stop_reason = ""
    stopped_round = 0
    noise_events = []
    accepted_snapshots = [
        _state_snapshot(
            0,
            initial_state,
            current_blocks,
            current_scored_blocks,
            before,
            dropped,
            report,
            current_score,
            current_cov,
            current_fill,
            current_reached,
            current_noise_ratio,
            current_dominant_keys,
            current_noise_keys,
        )
    ]

    temperature = start_temp
    history.append({
        "round": 0,
        "stage": "initial",
        "temperature": temperature,
        "score": current_score,
        "accepted": 1,
        "is_best": 1,
        "coverage": current_cov,
        "raw_blocks": len(raw_blocks),
        "blocks_before_partner_prune": len(before),
        "dropped_by_partner_cap": len(dropped),
        "blocks": len(current_blocks),
        "scored_blocks": len(current_scored_blocks),
        "noise_ratio": current_noise_ratio,
        "noise_filter_status": current_noise_info.get("status", ""),
        "dominant_blocks": current_gain["dominant_block_count"],
        "non_dominant_blocks": current_gain["non_dominant_block_count"],
        "delta_dominant_blocks": 0,
        "delta_noise_blocks": 0,
        "delta_dominant_aln": 0,
        "delta_noise_aln": 0,
        "best_blocks": len(best_blocks),
        "best_scored_blocks": len(best_scored_blocks),
        "target_reached": int(current_reached),
        "fill_mode_block_threshold": fill_mode_block_threshold,
        "score_mode": "noise_aware_self_stopping",
        "stopped": 0,
        "stop_reason": "",
        "rescue_rounds": initial_rescue_rounds,
        "rescue_action": initial_rescue_action,
        "guard_action": initial_guard_action,
        **initial_summary,
        **current_fill,
        "delta": 0.0,
        **initial_state,
    })

    for round_id in range(1, asa_rounds + 1):

        candidate = _perturb(best_state)

        before, blocks, dropped, report, candidate, rescue_rounds, rescue_action, guard_action = _apply_with_empty_rescue(
            raw_blocks,
            candidate,
            max_partners,
            partner_cap_mode,
        )

        scored_blocks, noise_ratio, noise_info = _score_blocks_after_noise_filter(
            blocks,
            dominant_chr_pair_keys,
            candidate,
        )
        gain = _gain_metrics(
            scored_blocks,
            dominant_chr_pair_keys,
            current_dominant_keys,
            current_noise_keys,
        )
        cov = len(collect_matched_chr_ids(scored_blocks)) / max(1, len(baseline_chr))
        fill = _chr_fill_metrics(scored_blocks, baseline_chr, target_chr_fill)
        reached = bool(scored_blocks)

        if not scored_blocks:
            score = -1e18
        else:
            score = _score(
                scored_blocks,
                target_blocks,
                cov,
                min_coverage,
                candidate,
                fill_metrics=fill,
                target_chr_fill=target_chr_fill,
                fill_mode_block_threshold=fill_mode_block_threshold,
                noise_ratio=noise_ratio,
            )

        delta = score - best_score

        accept = False
        is_best = False

        stop_relaxing = bool(scored_blocks and _should_stop_relaxing(gain))

        if stop_relaxing:
            accept = False
            stop_reason = (
                "non_dominant_gain_exceeds_dominant_gain:"
                f"blocks {gain['delta_noise_blocks']}>{gain['delta_dominant_blocks']};"
                f"aln {gain['delta_noise_aln']}>{gain['delta_dominant_aln']}"
            )
            stopped_round = round_id
            noise_events.append({
                "round": round_id,
                "state": dict(candidate),
                "gain": dict(gain),
                "score": score,
                "noise_ratio": noise_ratio,
            })
        elif not scored_blocks:
            accept = False
        elif delta > 0:
            accept = True
        else:
            prob = math.exp(delta / max(temperature, 1e-6))
            if random.random() < prob:
                accept = True

        if accept:
            best_state = dict(candidate)
            best_blocks = blocks
            best_scored_blocks = scored_blocks
            best_score = score
            best_before = before
            best_drop = dropped
            best_report = report
            best_cov = cov
            best_fill = fill
            best_reached = reached
            best_noise_ratio = noise_ratio
            current_dominant_keys = gain["dominant_keys"]
            current_noise_keys = gain["noise_keys"]
            accepted_snapshots.append(
                _state_snapshot(
                    round_id,
                    candidate,
                    blocks,
                    scored_blocks,
                    before,
                    dropped,
                    report,
                    score,
                    cov,
                    fill,
                    reached,
                    noise_ratio,
                    current_dominant_keys,
                    current_noise_keys,
                )
            )
            if delta > 0:
                is_best = True

        history.append({
            "round": round_id,
            "stage": "asa",
            "temperature": temperature,
            "score": score,
            "accepted": int(accept),
            "is_best": int(is_best),
            "coverage": cov,
            "raw_blocks": len(raw_blocks),
            "blocks_before_partner_prune": len(before),
            "dropped_by_partner_cap": len(dropped),
            "blocks": len(blocks),
            "scored_blocks": len(scored_blocks),
            "noise_ratio": noise_ratio,
            "noise_filter_status": noise_info.get("status", ""),
            "dominant_blocks": gain["dominant_block_count"],
            "non_dominant_blocks": gain["non_dominant_block_count"],
            "delta_dominant_blocks": gain["delta_dominant_blocks"],
            "delta_noise_blocks": gain["delta_noise_blocks"],
            "delta_dominant_aln": gain["delta_dominant_aln"],
            "delta_noise_aln": gain["delta_noise_aln"],
            "best_blocks": len(best_blocks),
            "best_scored_blocks": len(best_scored_blocks),
            "target_reached": int(best_reached),
            "fill_mode_block_threshold": fill_mode_block_threshold,
            "score_mode": "noise_aware_self_stopping",
            "stopped": int(stop_relaxing),
            "stop_reason": stop_reason if stop_relaxing else "",
            "rescue_rounds": rescue_rounds,
            "rescue_action": rescue_action,
            "guard_action": guard_action,
            **initial_summary,
            **fill,
            "delta": delta,
            **candidate,
        })

        if stop_relaxing:
            break

        temperature *= cooling

    if noise_events:
        rollback_snapshot, worst_event = _select_rollback_snapshot(accepted_snapshots, noise_events)
        if rollback_snapshot is not None:
            best_state = dict(rollback_snapshot["state"])
            best_blocks = rollback_snapshot["blocks"]
            best_scored_blocks = rollback_snapshot["scored_blocks"]
            best_score = rollback_snapshot["score"]
            best_before = rollback_snapshot["before"]
            best_drop = rollback_snapshot["dropped"]
            best_report = rollback_snapshot["report"]
            best_cov = rollback_snapshot["coverage"]
            best_fill = rollback_snapshot["fill"]
            best_reached = rollback_snapshot["reached"]
            best_noise_ratio = rollback_snapshot["noise_ratio"]
            current_dominant_keys = set(rollback_snapshot["dominant_keys"])
            current_noise_keys = set(rollback_snapshot["noise_keys"])
            stop_reason = (
                "rollback_to_before_max_noise_gain:"
                f"noise_round={worst_event['round'] if worst_event else 0};"
                f"rollback_round={rollback_snapshot['round']}"
            )

            single_param_rounds = max(1, asa_rounds - int(stopped_round or 0))
            tunable_params = (
                "min_block",
                "min_identity",
                "min_mapq",
                "min_pair_aln",
                "min_pair_blocks",
            )
            previous_param = None
            no_positive_params = set()

            for local_round in range(1, single_param_rounds + 1):
                available_params = [p for p in tunable_params if p != previous_param]
                if not available_params:
                    available_params = list(tunable_params)
                param_name = random.choice(available_params)
                candidate = _random_single_param_candidate(best_state, param_name)

                if candidate is None:
                    no_positive_params.add(param_name)
                    if len(no_positive_params) >= len(tunable_params):
                        stop_message = "single_param_all_params_no_positive_net_gain"
                    else:
                        continue
                else:
                    before, blocks, dropped, report, candidate, rescue_rounds, rescue_action, guard_action = _apply_with_empty_rescue(
                        raw_blocks,
                        candidate,
                        max_partners,
                        partner_cap_mode,
                    )
                    scored_blocks, noise_ratio, noise_info = _score_blocks_after_noise_filter(
                        blocks,
                        dominant_chr_pair_keys,
                        candidate,
                    )
                    gain = _gain_metrics(
                        scored_blocks,
                        dominant_chr_pair_keys,
                        current_dominant_keys,
                        current_noise_keys,
                    )
                    cov = len(collect_matched_chr_ids(scored_blocks)) / max(1, len(baseline_chr))
                    fill = _chr_fill_metrics(scored_blocks, baseline_chr, target_chr_fill)
                    reached = bool(scored_blocks)
                    if not scored_blocks:
                        score = -1e18
                    else:
                        score = _score(
                            scored_blocks,
                            target_blocks,
                            cov,
                            min_coverage,
                            candidate,
                            fill_metrics=fill,
                            target_chr_fill=target_chr_fill,
                            fill_mode_block_threshold=fill_mode_block_threshold,
                            noise_ratio=noise_ratio,
                        )
                    net_gain = _candidate_net_gain(gain)
                    delta = score - best_score
                    accept = False
                    if scored_blocks and net_gain > 0:
                        accept = True
                    elif scored_blocks:
                        prob = math.exp(min(0.0, net_gain) / max(temperature, 1e-6))
                        accept = random.random() < prob

                    selected = {
                        "param_name": param_name,
                        "candidate": candidate,
                        "before": before,
                        "blocks": blocks,
                        "dropped": dropped,
                        "report": report,
                        "scored_blocks": scored_blocks,
                        "noise_ratio": noise_ratio,
                        "noise_info": noise_info,
                        "gain": gain,
                        "coverage": cov,
                        "fill": fill,
                        "reached": reached,
                        "score": score,
                        "net_gain": net_gain,
                        "rescue_rounds": rescue_rounds,
                        "rescue_action": rescue_action,
                        "guard_action": guard_action,
                    }

                    previous_param = param_name
                    if net_gain > 0:
                        no_positive_params.clear()
                    else:
                        no_positive_params.add(param_name)

                    if accept:
                        gain = selected["gain"]
                        best_state = dict(selected["candidate"])
                        best_blocks = selected["blocks"]
                        best_scored_blocks = selected["scored_blocks"]
                        best_score = selected["score"]
                        best_before = selected["before"]
                        best_drop = selected["dropped"]
                        best_report = selected["report"]
                        best_cov = selected["coverage"]
                        best_fill = selected["fill"]
                        best_reached = selected["reached"]
                        best_noise_ratio = selected["noise_ratio"]
                        current_dominant_keys = gain["dominant_keys"]
                        current_noise_keys = gain["noise_keys"]

                    history.append({
                        "round": int(stopped_round or 0) + local_round,
                        "stage": "single_param_random",
                        "changed_param": selected["param_name"],
                        "temperature": temperature,
                        "score": selected["score"],
                        "accepted": int(accept),
                        "is_best": int(accept and delta > 0),
                        "coverage": selected["coverage"],
                        "raw_blocks": len(raw_blocks),
                        "blocks_before_partner_prune": len(selected["before"]),
                        "dropped_by_partner_cap": len(selected["dropped"]),
                        "blocks": len(selected["blocks"]),
                        "scored_blocks": len(selected["scored_blocks"]),
                        "noise_ratio": selected["noise_ratio"],
                        "noise_filter_status": selected["noise_info"].get("status", ""),
                        "dominant_blocks": gain["dominant_block_count"],
                        "non_dominant_blocks": gain["non_dominant_block_count"],
                        "delta_dominant_blocks": gain["delta_dominant_blocks"],
                        "delta_noise_blocks": gain["delta_noise_blocks"],
                        "delta_dominant_aln": gain["delta_dominant_aln"],
                        "delta_noise_aln": gain["delta_noise_aln"],
                        "net_gain": selected["net_gain"],
                        "best_blocks": len(best_blocks),
                        "best_scored_blocks": len(best_scored_blocks),
                        "target_reached": int(best_reached),
                        "fill_mode_block_threshold": fill_mode_block_threshold,
                        "score_mode": "noise_aware_rollback_random_single_param",
                        "stopped": 0,
                        "stop_reason": stop_reason,
                        "rescue_rounds": selected["rescue_rounds"],
                        "rescue_action": selected["rescue_action"],
                        "guard_action": selected["guard_action"],
                        **initial_summary,
                        **selected["fill"],
                        "delta": delta,
                        **candidate,
                    })

                    temperature *= cooling

                    if len(no_positive_params) < len(tunable_params):
                        continue
                    stop_message = "single_param_all_params_no_positive_net_gain"

                if len(no_positive_params) >= len(tunable_params):
                    history.append({
                        "round": int(stopped_round or 0) + local_round,
                        "stage": "single_param_stop",
                        "changed_param": "",
                        "temperature": temperature,
                        "score": best_score,
                        "accepted": 0,
                        "is_best": 0,
                        "coverage": best_cov,
                        "raw_blocks": len(raw_blocks),
                        "blocks_before_partner_prune": len(best_before),
                        "dropped_by_partner_cap": len(best_drop),
                        "blocks": len(best_blocks),
                        "scored_blocks": len(best_scored_blocks),
                        "noise_ratio": best_noise_ratio,
                        "noise_filter_status": "",
                        "dominant_blocks": len(current_dominant_keys),
                        "non_dominant_blocks": len(current_noise_keys),
                        "delta_dominant_blocks": 0,
                        "delta_noise_blocks": 0,
                        "delta_dominant_aln": 0,
                        "delta_noise_aln": 0,
                        "net_gain": 0,
                        "best_blocks": len(best_blocks),
                        "best_scored_blocks": len(best_scored_blocks),
                        "target_reached": int(best_reached),
                        "fill_mode_block_threshold": fill_mode_block_threshold,
                        "score_mode": "noise_aware_rollback_random_single_param",
                        "stopped": 1,
                        "stop_reason": stop_message,
                        "rescue_rounds": 0,
                        "rescue_action": "",
                        "guard_action": "",
                        **initial_summary,
                        **best_fill,
                        "delta": 0,
                        **best_state,
                    })
                    break

    best_state["max_partners_per_chr"] = max_partners
    best_state["partner_cap_mode"] = partner_cap_mode
    best_state["target_chr_fill"] = target_chr_fill
    best_state["target_reached"] = int(best_reached)
    best_state["fill_mode_block_threshold"] = fill_mode_block_threshold
    best_state["score_mode"] = (
        "noise_aware_rollback_random_single_param"
        if noise_events else
        "noise_aware_self_stopping"
    )
    best_state["noise_ratio"] = best_noise_ratio
    best_state["scored_blocks"] = len(best_scored_blocks)
    best_state["auto_stop_reason"] = stop_reason
    best_state["auto_stopped_round"] = stopped_round
    best_state["dominant_chr_pair_keys"] = dominant_chr_pair_keys
    best_state["dominant_chr_pair_records"] = dominant_chr_pair_records
    for key, value in initial_summary.items():
        best_state[key] = value
    for key, value in best_fill.items():
        best_state[key] = value

    return (
        best_scored_blocks,
        best_state,
        history,
        best_score,
        best_cov,
        best_drop,
        best_report,
        best_before,
    )
