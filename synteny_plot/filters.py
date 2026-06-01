# -*- coding: utf-8 -*-
"""Block-level filtering, capping, and chromosome-pair statistics."""

import argparse
from collections import defaultdict

from .utils import natural_key

def filter_blocks_basic(blocks, min_block=0, min_mapq=0, min_identity=0.0):
    """Apply only single-alignment filters to already parsed PAF blocks."""
    return [
        b for b in blocks
        if b["alnLen"] >= min_block
        and b["mapq"] >= min_mapq
        and b["identity"] >= min_identity
    ]


def collect_matched_chr_ids(blocks):
    matched = set()
    for b in blocks:
        matched.add((b["pair"], b["upper"]))
        matched.add((b["pair"] + 1, b["lower"]))
    return matched


def count_blocks_by_adjacent_pair(blocks, n_genomes):
    """Return block counts for each adjacent genome pair: 0=(G1-G2), 1=(G2-G3), ..."""
    counts = {i: 0 for i in range(max(0, n_genomes - 1))}
    for b in blocks:
        counts[int(b["pair"])] += 1
    return counts


def adjacent_pair_coverage(blocks, n_genomes):
    """Fraction of adjacent genome pairs that still have at least one drawable block."""
    counts = count_blocks_by_adjacent_pair(blocks, n_genomes)
    total = max(1, n_genomes - 1)
    covered = sum(1 for i in range(total) if counts.get(i, 0) > 0)
    return covered / total, counts


def clone_args_with_target(args, target_blocks):
    """Clone argparse Namespace but replace target_blocks for per-pair tuning."""
    d = vars(args).copy()
    d["target_blocks"] = max(1, int(target_blocks))
    return argparse.Namespace(**d)


def limit_blocks_per_adjacent_pair(blocks, max_blocks_per_pair):
    """Cap each adjacent genome-pair independently, keeping the longest blocks within each pair."""
    max_blocks_per_pair = int(max_blocks_per_pair)
    if max_blocks_per_pair <= 0:
        return blocks
    grouped = defaultdict(list)
    for b in blocks:
        grouped[int(b["pair"])].append(b)
    out = []
    for pair in sorted(grouped):
        vals = grouped[pair]
        if len(vals) > max_blocks_per_pair:
            vals = sorted(vals, key=lambda x: x["alnLen"], reverse=True)[:max_blocks_per_pair]
        out.extend(vals)
    return out


def limit_blocks_balanced_by_pair(blocks, max_total, n_genomes):
    """
    Cap total blocks while preserving all non-empty adjacent genome pairs.

    The old global cap could keep only G1-G2 if those blocks were longest. This function allocates
    quota across non-empty adjacent pairs and then takes the longest blocks inside each pair.
    """
    max_total = int(max_total)
    if max_total <= 0 or len(blocks) <= max_total:
        return blocks

    grouped = defaultdict(list)
    for b in blocks:
        grouped[int(b["pair"])].append(b)
    grouped = {p: sorted(v, key=lambda x: x["alnLen"], reverse=True) for p, v in grouped.items() if v}
    if not grouped:
        return []

    active = sorted(grouped)
    if max_total < len(active):
        # Degenerate case: target smaller than number of pair rows. Keep one longest block per earliest rows.
        out = []
        for p in active[:max_total]:
            out.append(grouped[p][0])
        return out

    quota = {p: 1 for p in active}
    remaining = max_total - len(active)

    while remaining > 0:
        candidates = [p for p in active if quota[p] < len(grouped[p])]
        if not candidates:
            break
        # Give the next slot to the pair that still has the most remaining candidate blocks.
        p = max(candidates, key=lambda k: (len(grouped[k]) - quota[k], len(grouped[k]), -k))
        quota[p] += 1
        remaining -= 1

    out = []
    for p in active:
        out.extend(grouped[p][:quota[p]])
    return out


def count_chr_pair_blocks(blocks):
    """Summarize block count and aligned length for each adjacent chromosome pair."""
    stat = {}
    for b in blocks:
        key = (b["pair"], b["upper"], b["lower"])
        if key not in stat:
            stat[key] = {
                "pair": b["pair"],
                "upperGenome": b["upperGenome"],
                "lowerGenome": b["lowerGenome"],
                "upper": b["upper"],
                "lower": b["lower"],
                "block_count": 0,
                "aln_sum": 0,
            }
        stat[key]["block_count"] += 1
        stat[key]["aln_sum"] += b["alnLen"]
    return stat


def filter_blocks_by_chr_pair_count_with_report(blocks, min_blocks):
    """
    Drop weak chromosome-pair links by block count.

    This is stricter and more explicit than only dropping unmatched chromosomes:
    if chrA has one good partner and one weak partner, only the weak partner is removed.
    If all partners of chrA are weak, chrA disappears from the plot naturally.
    """
    if min_blocks <= 0:
        return blocks, []

    stat = count_chr_pair_blocks(blocks)
    keep = {key for key, x in stat.items() if x["block_count"] >= min_blocks}
    dropped = [x for key, x in stat.items() if key not in keep]
    kept_blocks = [b for b in blocks if (b["pair"], b["upper"], b["lower"]) in keep]
    dropped.sort(key=lambda x: (x["pair"], natural_key(x["upper"]), natural_key(x["lower"])))
    return kept_blocks, dropped


def filter_final_chr_pair_noise_by_count(blocks,
                                         min_pair_blocks=20,
                                         keep_top_partners=2,
                                         min_block_ratio=0.05,
                                         never_filter_to_zero=True):
    """Final drawing-only cleanup: drop sparse chr-pairs by block count, not alignment quality."""
    if not blocks:
        return blocks, [], {"status": "empty_input", "kept_pairs": 0, "dropped_pairs": 0}

    min_pair_blocks = max(1, int(min_pair_blocks))
    keep_top_partners = max(1, int(keep_top_partners))
    min_block_ratio = max(0.0, float(min_block_ratio))

    stat = count_chr_pair_blocks(blocks)
    by_upper = defaultdict(list)
    by_lower = defaultdict(list)
    for key, x in stat.items():
        by_upper[(x["pair"], x["upper"])].append((key, x["block_count"], x["aln_sum"]))
        by_lower[(x["pair"], x["lower"])].append((key, x["block_count"], x["aln_sum"]))

    for group in by_upper.values():
        upper_partner_count = len(group)
        group.sort(key=lambda z: (z[1], z[2]), reverse=True)
        max_count = max(1, group[0][1])
        for rank, (key, cnt, _aln) in enumerate(group, start=1):
            stat[key]["upper_block_rank"] = rank
            stat[key]["upper_max_block_count"] = max_count
            stat[key]["upper_block_ratio"] = cnt / max_count
            stat[key]["upper_partner_count"] = upper_partner_count

    for group in by_lower.values():
        lower_partner_count = len(group)
        group.sort(key=lambda z: (z[1], z[2]), reverse=True)
        max_count = max(1, group[0][1])
        for rank, (key, cnt, _aln) in enumerate(group, start=1):
            stat[key]["lower_block_rank"] = rank
            stat[key]["lower_max_block_count"] = max_count
            stat[key]["lower_block_ratio"] = cnt / max_count
            stat[key]["lower_partner_count"] = lower_partner_count

    keep = set()
    report = []
    for key, x in stat.items():
        x.setdefault("upper_block_rank", 999999)
        x.setdefault("lower_block_rank", 999999)
        x.setdefault("upper_max_block_count", 0)
        x.setdefault("lower_max_block_count", 0)
        x.setdefault("upper_block_ratio", 0.0)
        x.setdefault("lower_block_ratio", 0.0)
        x.setdefault("upper_partner_count", 0)
        x.setdefault("lower_partner_count", 0)

        reasons = []
        protected_single_partner = (
            x["upper_partner_count"] == 1
            or x["lower_partner_count"] == 1
        )
        protected_by_rank = (
            x["upper_block_rank"] <= keep_top_partners
            or x["lower_block_rank"] <= keep_top_partners
        )
        protected_by_ratio = (
            x["upper_block_ratio"] >= min_block_ratio
            or x["lower_block_ratio"] >= min_block_ratio
        )

        if protected_single_partner:
            pass
        elif x["block_count"] < min_pair_blocks:
            reasons.append(f"block_count<{min_pair_blocks}")
        elif not protected_by_rank and not protected_by_ratio:
            reasons.append(
                f"not_top{keep_top_partners}_and_block_ratio<{min_block_ratio}"
            )

        row = dict(x)
        row["min_pair_blocks"] = min_pair_blocks
        row["keep_top_partners"] = keep_top_partners
        row["min_block_ratio"] = min_block_ratio
        row["status"] = "dropped" if reasons else "kept"
        row["reason"] = ";".join(reasons)
        report.append(row)
        if not reasons:
            keep.add(key)

    kept_blocks = [b for b in blocks if (b["pair"], b["upper"], b["lower"]) in keep]

    if not kept_blocks and never_filter_to_zero:
        for row in report:
            if row["status"] == "dropped":
                row["reason"] = "rollback_filter_would_remove_all;" + row.get("reason", "")
            row["status"] = "kept"
        return blocks, report, {
            "status": "rollback_filter_would_remove_all",
            "kept_pairs": len(stat),
            "dropped_pairs": 0,
            "min_pair_blocks": min_pair_blocks,
            "keep_top_partners": keep_top_partners,
            "min_block_ratio": min_block_ratio,
        }

    dropped_count = sum(1 for x in report if x["status"] == "dropped")
    report.sort(key=lambda x: (x["pair"], natural_key(x["upper"]), natural_key(x["lower"])))
    return kept_blocks, report, {
        "status": "applied",
        "kept_pairs": len(stat) - dropped_count,
        "dropped_pairs": dropped_count,
        "min_pair_blocks": min_pair_blocks,
        "keep_top_partners": keep_top_partners,
        "min_block_ratio": min_block_ratio,
    }


def _is_non_continuous_chr_pair(blocks, continuity_gap_bp):
    if len(blocks) <= 1:
        return True
    continuity_gap_bp = max(0, int(continuity_gap_bp or 0))
    vals = sorted(
        blocks,
        key=lambda b: (
            int(b.get("upperStart", 0) or 0),
            int(b.get("lowerStart", 0) or 0),
        ),
    )
    continuous_links = 0
    for prev, cur in zip(vals, vals[1:]):
        upper_gap = max(0, int(cur.get("upperStart", 0) or 0) - int(prev.get("upperEnd", 0) or 0))
        lower_gap_forward = abs(int(cur.get("lowerStart", 0) or 0) - int(prev.get("lowerEnd", 0) or 0))
        lower_gap_reverse = abs(int(cur.get("lowerEnd", 0) or 0) - int(prev.get("lowerStart", 0) or 0))
        lower_gap = min(lower_gap_forward, lower_gap_reverse)
        if upper_gap <= continuity_gap_bp and lower_gap <= continuity_gap_bp:
            continuous_links += 1
    return continuous_links == 0


def filter_non_dominant_noise_by_dominant_pairs(blocks,
                                                dominant_chr_pair_keys=None,
                                                pair_params=None,
                                                continuity_gap_bp=5_000_000,
                                                borderline_ratio=1.25,
                                                never_filter_to_zero=True):
    """
    Final drawing-only noise cleanup using initialization-time dominant pairs.

    A chromosome-pair is protected if it was selected as a dominant relationship
    during initialization. Non-dominant single-block pairs are removed directly.
    Non-dominant multi-block pairs are removed only when they are non-continuous
    and their support is close to the loosest retained filter thresholds.
    """
    if not blocks:
        return blocks, [], {"status": "empty_input", "kept_pairs": 0, "dropped_pairs": 0}

    dominant_keys = {
        (int(k[0]), str(k[1]), str(k[2]))
        for k in (dominant_chr_pair_keys or [])
        if len(k) >= 3
    }
    if not dominant_keys:
        report = []
        stat = count_chr_pair_blocks(blocks)
        for x in stat.values():
            row = dict(x)
            row["status"] = "kept"
            row["reason"] = "no_dominant_init_pairs_available"
            report.append(row)
        report.sort(key=lambda x: (x["pair"], natural_key(x["upper"]), natural_key(x["lower"])))
        return blocks, report, {
            "status": "no_dominant_init_pairs_available",
            "kept_pairs": len(stat),
            "dropped_pairs": 0,
        }
    pair_params = pair_params or {}
    borderline_ratio = max(1.0, float(borderline_ratio or 1.0))

    grouped = defaultdict(list)
    for b in blocks:
        grouped[(int(b["pair"]), b["upper"], b["lower"])].append(b)

    keep = set()
    report = []
    for key, vals in grouped.items():
        pair, upper, lower = key
        first = vals[0]
        params = pair_params.get(pair, {})
        min_block = max(0, int(params.get("min_block", 0) or 0))
        min_pair_aln = max(0, int(params.get("min_pair_aln", 0) or 0))
        min_identity = max(0.0, float(params.get("min_identity", 0.0) or 0.0))
        min_mapq = max(0, int(params.get("min_mapq", 0) or 0))

        aln_lens = [int(b.get("alnLen", 0) or 0) for b in vals]
        identities = [float(b.get("identity", 0.0) or 0.0) for b in vals]
        mapqs = [int(b.get("mapq", 0) or 0) for b in vals]
        aln_sum = sum(aln_lens)
        max_aln = max(aln_lens) if aln_lens else 0
        median_aln = sorted(aln_lens)[len(aln_lens) // 2] if aln_lens else 0
        min_obs_identity = min(identities) if identities else 0.0
        min_obs_mapq = min(mapqs) if mapqs else 0

        protected = key in dominant_keys
        single_block_noise = (len(vals) == 1 and not protected)
        non_continuous = _is_non_continuous_chr_pair(vals, continuity_gap_bp)
        borderline_length = (
            (min_block > 0 and median_aln <= int(min_block * borderline_ratio))
            or (min_pair_aln > 0 and aln_sum <= int(min_pair_aln * borderline_ratio))
        )
        borderline_quality = (
            (min_identity > 0 and min_obs_identity <= min(1.0, min_identity + 0.02))
            or (min_mapq > 0 and min_obs_mapq <= min_mapq + 5)
        )

        reasons = []
        if protected:
            reasons.append("dominant_init_pair")
        elif single_block_noise:
            reasons.append("drop_non_dominant_single_block")
        elif non_continuous and (borderline_length or borderline_quality):
            reasons.append("drop_non_dominant_noncontinuous_borderline")

        status = "kept"
        if reasons and reasons[0].startswith("drop_"):
            status = "dropped"
        else:
            keep.add(key)

        report.append({
            "pair": pair,
            "upperGenome": first.get("upperGenome", ""),
            "lowerGenome": first.get("lowerGenome", ""),
            "upper": upper,
            "lower": lower,
            "block_count": len(vals),
            "aln_sum": aln_sum,
            "dominant_partner": int(protected),
            "non_continuous": int(non_continuous),
            "continuity_gap_bp": int(continuity_gap_bp),
            "borderline_ratio": borderline_ratio,
            "min_block": min_block,
            "min_pair_aln": min_pair_aln,
            "min_identity": min_identity,
            "min_mapq": min_mapq,
            "median_aln": median_aln,
            "max_aln": max_aln,
            "min_obs_identity": min_obs_identity,
            "min_obs_mapq": min_obs_mapq,
            "status": status,
            "reason": ";".join(reasons),
        })

    kept_blocks = [b for b in blocks if (int(b["pair"]), b["upper"], b["lower"]) in keep]
    dropped_count = sum(1 for x in report if x["status"] == "dropped")

    if not kept_blocks and never_filter_to_zero:
        for row in report:
            if row["status"] == "dropped":
                row["reason"] = "rollback_filter_would_remove_all;" + row.get("reason", "")
            row["status"] = "kept"
        return blocks, report, {
            "status": "rollback_filter_would_remove_all",
            "kept_pairs": len(grouped),
            "dropped_pairs": 0,
        }

    report.sort(key=lambda x: (x["pair"], natural_key(x["upper"]), natural_key(x["lower"])))
    return kept_blocks, report, {
        "status": "applied",
        "kept_pairs": len(grouped) - dropped_count,
        "dropped_pairs": dropped_count,
    }


def _block_signature(b):
    return (
        b.get("pair"), b.get("upper"), b.get("upperStart"), b.get("upperEnd"),
        b.get("lower"), b.get("lowerStart"), b.get("lowerEnd"), b.get("strand"),
    )


def _fill_chr_pair_stat(blocks):
    stat = count_chr_pair_blocks(blocks)
    for x in stat.values():
        x["upper_span_sum"] = 0
        x["lower_span_sum"] = 0
        x["upper_fill_rate"] = 0.0
        x["lower_fill_rate"] = 0.0
        x["max_fill_rate"] = 0.0
    for b in blocks:
        key = (b["pair"], b["upper"], b["lower"])
        x = stat[key]
        upper_span = max(0, int(b.get("upperEnd", 0) or 0) - int(b.get("upperStart", 0) or 0))
        lower_span = max(0, int(b.get("lowerEnd", 0) or 0) - int(b.get("lowerStart", 0) or 0))
        x["upper_span_sum"] += upper_span
        x["lower_span_sum"] += lower_span
        x["upperLen"] = max(int(x.get("upperLen", 0) or 0), int(b.get("upperLen", 0) or 0))
        x["lowerLen"] = max(int(x.get("lowerLen", 0) or 0), int(b.get("lowerLen", 0) or 0))
    for x in stat.values():
        x["upper_fill_rate"] = min(1.0, x["upper_span_sum"] / max(1, int(x.get("upperLen", 0) or 0)))
        x["lower_fill_rate"] = min(1.0, x["lower_span_sum"] / max(1, int(x.get("lowerLen", 0) or 0)))
        x["max_fill_rate"] = max(x["upper_fill_rate"], x["lower_fill_rate"])
    return stat


def filter_dominant_chr_pairs_by_fill_with_rescue(blocks,
                                                  candidate_blocks=None,
                                                  min_fill_rate=0.99,
                                                  min_dominant_blocks=10,
                                                  never_filter_to_zero=True):
    """
    Final drawing filter for clean figures.

    Keep only each chromosome's dominant partner. A dominant pair is drawable if it
    fills at least min_fill_rate of either chromosome, or if it can retain/rescue at
    least min_dominant_blocks from the same chr-pair.
    """
    if not blocks:
        return blocks, [], {"status": "empty_input", "kept_pairs": 0, "dropped_pairs": 0}

    candidate_blocks = list(candidate_blocks or blocks)
    min_fill_rate = max(0.0, min(1.0, float(min_fill_rate)))
    min_dominant_blocks = max(0, int(min_dominant_blocks or 0))

    candidate_stat = _fill_chr_pair_stat(candidate_blocks)
    current_stat = _fill_chr_pair_stat(blocks)
    if not candidate_stat:
        return blocks, [], {"status": "empty_candidates", "kept_pairs": 0, "dropped_pairs": 0}

    by_upper = defaultdict(list)
    by_lower = defaultdict(list)
    for key, x in candidate_stat.items():
        by_upper[(x["pair"], x["upper"])].append(key)
        by_lower[(x["pair"], x["lower"])].append(key)

    dominant_keys = set()
    for keys in by_upper.values():
        dominant_keys.add(max(
            keys,
            key=lambda k: (
                candidate_stat[k]["upper_fill_rate"],
                candidate_stat[k]["block_count"],
                candidate_stat[k]["aln_sum"],
            ),
        ))
    for keys in by_lower.values():
        dominant_keys.add(max(
            keys,
            key=lambda k: (
                candidate_stat[k]["lower_fill_rate"],
                candidate_stat[k]["block_count"],
                candidate_stat[k]["aln_sum"],
            ),
        ))

    block_protected_keys = {
        key for key, x in candidate_stat.items()
        if int(x.get("block_count", 0) or 0) >= min_dominant_blocks
    }
    fill_protected_keys = {
        key for key, x in candidate_stat.items()
        if float(x.get("max_fill_rate", 0.0) or 0.0) >= min_fill_rate
    }
    protected_keys = dominant_keys | block_protected_keys | fill_protected_keys

    current = [b for b in blocks if (b["pair"], b["upper"], b["lower"]) in protected_keys]
    current_sigs = {_block_signature(b) for b in current}
    current_count = defaultdict(int)
    for b in current:
        current_count[(b["pair"], b["upper"], b["lower"])] += 1

    candidates_by_key = defaultdict(list)
    for b in candidate_blocks:
        key = (b["pair"], b["upper"], b["lower"])
        if key in dominant_keys:
            candidates_by_key[key].append(b)
    for vals in candidates_by_key.values():
        vals.sort(key=lambda b: (int(b.get("alnLen", 0) or 0), float(b.get("identity", 0) or 0)), reverse=True)

    report = []
    keep_keys = set()
    for key in sorted(candidate_stat, key=lambda k: (k[0], natural_key(k[1]), natural_key(k[2]))):
        cand = candidate_stat[key]
        cur = current_stat.get(key, {})
        reasons = []
        status = "kept"
        rescued = 0
        before_count = current_count.get(key, 0)

        block_ok = cand["block_count"] >= min_dominant_blocks
        fill_ok = cand["max_fill_rate"] >= min_fill_rate

        if key not in protected_keys:
            status = "dropped"
            reasons.append(f"non_dominant_partner_and_fill_rate<{min_fill_rate}_and_block_count<{min_dominant_blocks}")
        else:
            need = max(0, min_dominant_blocks - before_count)
            if key in dominant_keys and need > 0 and not fill_ok:
                for b in candidates_by_key.get(key, []):
                    sig = _block_signature(b)
                    if sig in current_sigs:
                        continue
                    current.append(b)
                    current_sigs.add(sig)
                    before_count += 1
                    rescued += 1
                    if before_count >= min_dominant_blocks:
                        break
            if fill_ok or block_ok or before_count >= min_dominant_blocks:
                keep_keys.add(key)
            else:
                status = "dropped"
                reasons.append(f"fill_rate<{min_fill_rate}_and_block_count<{min_dominant_blocks}")

        row = dict(cand)
        row["current_block_count"] = int(cur.get("block_count", 0) or 0)
        row["rescued_blocks"] = rescued
        row["blocks_after_rescue"] = before_count if key in dominant_keys else int(cur.get("block_count", 0) or 0)
        row["min_fill_rate"] = min_fill_rate
        row["min_dominant_blocks"] = min_dominant_blocks
        row["dominant_partner"] = int(key in dominant_keys)
        row["block_protected"] = int(key in block_protected_keys)
        row["fill_protected"] = int(key in fill_protected_keys)
        row["status"] = status
        row["reason"] = ";".join(reasons)
        report.append(row)

    kept_blocks = [b for b in current if (b["pair"], b["upper"], b["lower"]) in keep_keys]
    if not kept_blocks and never_filter_to_zero:
        for row in report:
            if row["status"] == "dropped":
                row["reason"] = "rollback_filter_would_remove_all;" + row.get("reason", "")
            row["status"] = "kept"
        return blocks, report, {
            "status": "rollback_filter_would_remove_all",
            "kept_pairs": len(current_stat),
            "dropped_pairs": 0,
            "min_fill_rate": min_fill_rate,
            "min_dominant_blocks": min_dominant_blocks,
        }

    dropped_count = sum(1 for row in report if row["status"] == "dropped")
    return kept_blocks, report, {
        "status": "applied",
        "kept_pairs": len(keep_keys),
        "dropped_pairs": dropped_count,
        "min_fill_rate": min_fill_rate,
        "min_dominant_blocks": min_dominant_blocks,
    }


def add_chr_pair_dominance_metrics(stat):
    """
    Add relative dominance metrics to chromosome-pair summary.

    For each adjacent genome pair:
      - upper_ratio = pair aligned length / strongest partner aligned length for this upper chromosome
      - lower_ratio = pair aligned length / strongest partner aligned length for this lower chromosome
      - upper_rank/lower_rank are ranks by aligned length within each chromosome's partners
    """
    if not stat:
        return stat

    by_upper = defaultdict(list)
    by_lower = defaultdict(list)
    for key, x in stat.items():
        by_upper[(x["pair"], x["upper"])].append((key, x["aln_sum"], x["block_count"]))
        by_lower[(x["pair"], x["lower"])].append((key, x["aln_sum"], x["block_count"]))

    for group in by_upper.values():
        group.sort(key=lambda z: (z[1], z[2]), reverse=True)
        max_aln = max(1, group[0][1])
        max_count = max(1, group[0][2])
        for rank, (key, aln, cnt) in enumerate(group, start=1):
            stat[key]["upper_rank"] = rank
            stat[key]["upper_ratio"] = aln / max_aln
            stat[key]["upper_count_ratio"] = cnt / max_count

    for group in by_lower.values():
        group.sort(key=lambda z: (z[1], z[2]), reverse=True)
        max_aln = max(1, group[0][1])
        max_count = max(1, group[0][2])
        for rank, (key, aln, cnt) in enumerate(group, start=1):
            stat[key]["lower_rank"] = rank
            stat[key]["lower_ratio"] = aln / max_aln
            stat[key]["lower_count_ratio"] = cnt / max_count

    for x in stat.values():
        x.setdefault("upper_rank", 999999)
        x.setdefault("lower_rank", 999999)
        x.setdefault("upper_ratio", 0.0)
        x.setdefault("lower_ratio", 0.0)
        x.setdefault("upper_count_ratio", 0.0)
        x.setdefault("lower_count_ratio", 0.0)
    return stat


def summarize_chr_pair_dominance(blocks):
    stat = count_chr_pair_blocks(blocks)
    return add_chr_pair_dominance_metrics(stat)


def filter_blocks_by_dominant_chr_pairs_with_report(blocks,
                                                    weak_pair_ratio=0.05,
                                                    max_partners_per_chr=2,
                                                    min_pair_blocks=2,
                                                    never_filter_to_zero=True):
    """
    Remove tiny off-target chromosome-pair links using relative dominance, not a hard absolute cutoff.

    Keep a chr-pair if:
      1. it is within the top N partners for both chromosomes, and
      2. its aligned length is at least weak_pair_ratio of the dominant partner
         from both the upper and lower chromosome perspectives, and
      3. it has at least min_pair_blocks blocks.

    If this would remove everything, return the original blocks and mark the filter as rolled back.
    """
    if not blocks:
        return blocks, [], {"status": "empty_input", "kept_pairs": 0, "dropped_pairs": 0}

    weak_pair_ratio = max(0.0, float(weak_pair_ratio))
    max_partners_per_chr = max(1, int(max_partners_per_chr))
    min_pair_blocks = max(1, int(min_pair_blocks))

    stat = summarize_chr_pair_dominance(blocks)
    keep = set()
    dropped = []

    for key, x in stat.items():
        reasons = []
        if x["block_count"] < min_pair_blocks:
            reasons.append(f"block_count<{min_pair_blocks}")
        if x["upper_rank"] > max_partners_per_chr:
            reasons.append(f"upper_rank>{max_partners_per_chr}")
        if x["lower_rank"] > max_partners_per_chr:
            reasons.append(f"lower_rank>{max_partners_per_chr}")
        if x["upper_ratio"] < weak_pair_ratio:
            reasons.append(f"upper_ratio<{weak_pair_ratio}")
        if x["lower_ratio"] < weak_pair_ratio:
            reasons.append(f"lower_ratio<{weak_pair_ratio}")

        if reasons:
            y = dict(x)
            y["reason"] = ";".join(reasons)
            dropped.append(y)
        else:
            keep.add(key)

    kept_blocks = [b for b in blocks if (b["pair"], b["upper"], b["lower"]) in keep]

    # Guardrail: a clean-mode filter is allowed to improve the plot, not destroy it.
    if not kept_blocks and never_filter_to_zero:
        for x in dropped:
            x["reason"] = "rollback_filter_would_remove_all;" + x.get("reason", "")
        info = {
            "status": "rollback_filter_would_remove_all",
            "kept_pairs": len(stat),
            "dropped_pairs": 0,
            "weak_pair_ratio": weak_pair_ratio,
            "max_partners_per_chr": max_partners_per_chr,
            "min_pair_blocks": min_pair_blocks,
        }
        return blocks, dropped, info

    dropped.sort(key=lambda x: (x["pair"], natural_key(x["upper"]), natural_key(x["lower"])))
    info = {
        "status": "applied",
        "kept_pairs": len(keep),
        "dropped_pairs": len(dropped),
        "weak_pair_ratio": weak_pair_ratio,
        "max_partners_per_chr": max_partners_per_chr,
        "min_pair_blocks": min_pair_blocks,
    }
    return kept_blocks, dropped, info


def apply_filter_params(blocks, params):
    out = filter_blocks_basic(
        blocks,
        min_block=params.get("min_block", 0),
        min_mapq=params.get("min_mapq", 0),
        min_identity=params.get("min_identity", 0.0),
    )
    out = filter_blocks_by_pair_aln(out, params.get("min_pair_aln", 0))
    out = filter_blocks_by_pair_count(out, params.get("min_pair_blocks", 0))
    out = keep_top_blocks_by_max_bottom(out, params.get("max_bottom_per_top", 0))
    return out


def filter_blocks_by_pair_aln(blocks, min_pair_aln):
    """
    Keep chromosome-pair connections only if total aligned length reaches threshold.
    Pair key includes genome-pair index.
    """
    if min_pair_aln <= 0:
        return blocks

    pair_sum = defaultdict(int)
    for b in blocks:
        pair_sum[(b["pair"], b["upper"], b["lower"])] += b["alnLen"]

    return [
        b for b in blocks
        if pair_sum[(b["pair"], b["upper"], b["lower"])] >= min_pair_aln
    ]


def filter_blocks_by_pair_count(blocks, min_pair_blocks):
    if min_pair_blocks <= 0:
        return blocks

    pair_count = defaultdict(int)
    for b in blocks:
        pair_count[(b["pair"], b["upper"], b["lower"])] += 1

    return [
        b for b in blocks
        if pair_count[(b["pair"], b["upper"], b["lower"])] >= min_pair_blocks
    ]


def keep_top_blocks_by_max_bottom(blocks, max_bottom_per_top):
    """
    Per adjacent genome-pair and per upper chromosome, keep top N lower chromosomes.
    """
    if max_bottom_per_top <= 0:
        return blocks

    per_pair_upper = defaultdict(lambda: defaultdict(int))
    for b in blocks:
        per_pair_upper[(b["pair"], b["upper"])][b["lower"]] += b["alnLen"]

    allowed = set()
    for key, lower_map in per_pair_upper.items():
        ranked = sorted(lower_map.items(), key=lambda x: x[1], reverse=True)
        for lower, total in ranked[:max_bottom_per_top]:
            allowed.add((key[0], key[1], lower))

    return [b for b in blocks if (b["pair"], b["upper"], b["lower"]) in allowed]


def limit_blocks(blocks, max_blocks):
    if max_blocks <= 0 or len(blocks) <= max_blocks:
        return blocks
    return sorted(blocks, key=lambda b: b["alnLen"], reverse=True)[:max_blocks]


def _rank_chr_partner_pairs(blocks):
    """Return per chromosome-pair statistics ranked by block count, then aligned length."""
    stat = count_chr_pair_blocks(blocks)
    by_upper = defaultdict(list)
    by_lower = defaultdict(list)
    for key, x in stat.items():
        by_upper[(x["pair"], x["upper"])].append(key)
        by_lower[(x["pair"], x["lower"])].append(key)

    upper_rank = {}
    lower_rank = {}

    for group_key, keys in by_upper.items():
        keys = sorted(
            keys,
            key=lambda k: (-stat[k]["block_count"], -stat[k]["aln_sum"], natural_key(stat[k]["lower"])),
        )
        for idx, key in enumerate(keys, start=1):
            upper_rank[key] = idx

    for group_key, keys in by_lower.items():
        keys = sorted(
            keys,
            key=lambda k: (-stat[k]["block_count"], -stat[k]["aln_sum"], natural_key(stat[k]["upper"])),
        )
        for idx, key in enumerate(keys, start=1):
            lower_rank[key] = idx

    for key, x in stat.items():
        x["upper_block_rank"] = upper_rank.get(key, 999999)
        x["lower_block_rank"] = lower_rank.get(key, 999999)

    return stat


def prune_chr_partners_by_block_count(blocks, max_partners=4, mode="both"):
    """Prune chromosome-pair links by partner count.

    Ranking rule: block_count descending, then aln_sum descending.

    mode:
      upper: each upper chromosome keeps at most N lower partners.
      lower: each lower chromosome keeps at most N upper partners.
      both:  keep only links that pass both upper-side and lower-side top-N limits.

    Returns:
      kept_blocks, dropped_blocks, pair_report
    """
    max_partners = int(max_partners or 0)
    if max_partners <= 0 or not blocks:
        return blocks, [], []

    mode = (mode or "both").lower()
    if mode not in {"upper", "lower", "both"}:
        raise ValueError(f"Unknown partner-cap mode: {mode}; expected upper/lower/both")

    stat = _rank_chr_partner_pairs(blocks)
    keep = set()
    pair_report = []

    for key, x in stat.items():
        reasons = []
        if mode in {"upper", "both"} and x["upper_block_rank"] > max_partners:
            reasons.append(f"upper_partner_rank>{max_partners}")
        if mode in {"lower", "both"} and x["lower_block_rank"] > max_partners:
            reasons.append(f"lower_partner_rank>{max_partners}")
        status = "dropped" if reasons else "kept"
        if not reasons:
            keep.add(key)
        row = dict(x)
        row["max_partners_per_chr"] = max_partners
        row["partner_cap_mode"] = mode
        row["status"] = status
        row["reason"] = ";".join(reasons)
        pair_report.append(row)

    kept_blocks = []
    dropped_blocks = []
    for idx, b in enumerate(blocks):
        key = (b["pair"], b["upper"], b["lower"])
        if key in keep:
            kept_blocks.append(b)
        else:
            y = dict(b)
            y["source_block_index_after_filters"] = idx
            sx = stat.get(key, {})
            y["upper_block_rank"] = sx.get("upper_block_rank", "")
            y["lower_block_rank"] = sx.get("lower_block_rank", "")
            y["chr_pair_block_count"] = sx.get("block_count", "")
            y["chr_pair_aln_sum"] = sx.get("aln_sum", "")
            y["reason"] = next((r["reason"] for r in pair_report if (r["pair"], r["upper"], r["lower"]) == key), "partner_cap")
            dropped_blocks.append(y)

    pair_report.sort(key=lambda x: (x["pair"], natural_key(x["upper"]), x["upper_block_rank"], natural_key(x["lower"])))
    return kept_blocks, dropped_blocks, pair_report


def chr_block_count_rows(blocks, stage=""):
    """Summarize block counts per chromosome per adjacent pair and row role."""
    rows = []
    upper_counts = defaultdict(lambda: {"block_count": 0, "partners": set(), "aln_sum": 0, "genome": ""})
    lower_counts = defaultdict(lambda: {"block_count": 0, "partners": set(), "aln_sum": 0, "genome": ""})
    for b in blocks:
        uk = (b["pair"], "upper", b["upper"])
        lk = (b["pair"], "lower", b["lower"])
        upper_counts[uk]["block_count"] += 1
        upper_counts[uk]["partners"].add(b["lower"])
        upper_counts[uk]["aln_sum"] += b["alnLen"]
        upper_counts[uk]["genome"] = b["upperGenome"]
        lower_counts[lk]["block_count"] += 1
        lower_counts[lk]["partners"].add(b["upper"])
        lower_counts[lk]["aln_sum"] += b["alnLen"]
        lower_counts[lk]["genome"] = b["lowerGenome"]

    for (pair, role, chrom), x in list(upper_counts.items()) + list(lower_counts.items()):
        rows.append({
            "stage": stage,
            "adjacent_pair": pair,
            "row_role": role,
            "genome": x["genome"],
            "chr": chrom,
            "block_count": x["block_count"],
            "partner_count": len(x["partners"]),
            "aln_sum": x["aln_sum"],
        })
    rows.sort(key=lambda r: (r["stage"], r["adjacent_pair"], r["row_role"], natural_key(r["chr"])))
    return rows
