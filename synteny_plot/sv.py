# -*- coding: utf-8 -*-
"""Structural variant detection from synteny blocks.

v3: Strand-based inversion detection.
- Determine dominant strand per chr-pair.
- Runs of non-dominant strand blocks are inversion candidates.
- Two-tier span filtering (coordinate span or alignment length sum).
"""

from collections import defaultdict

from .sv_models import SVEvent
from .utils import natural_key

SV_PRIORITY = {
    "SYN": 0,
    "INS_DEL": 1,
    "INS": 2,
    "DEL": 2,
    "FRAG_INV": 2,
    "INV": 3,
    "TRANS": 4,
    "FUSION": 4,
}


def _sort_blocks_by_chr_pair(blocks):
    """Group blocks by (pair, upper, lower) and sort each group by upperStart."""
    grouped = defaultdict(list)
    for b in blocks:
        grouped[(b["pair"], b["upper"], b["lower"])].append(b)
    for key in grouped:
        grouped[key].sort(key=lambda x: (int(x["upperStart"]), int(x["lowerStart"])))
    return grouped


def detect_inversions(blocks, min_blocks=5, min_size=1_000_000):
    """
    Detect inversions using strand-based approach.

    For each chr-pair:
      1. Determine the dominant strand (the one with more blocks).
      2. Find runs of consecutive blocks on the non-dominant strand.
      3. A run with >= min_blocks and span >= min_size (or alnLen sum >= min_size)
         is marked as an inversion candidate.
    """
    events = []
    grouped = _sort_blocks_by_chr_pair(blocks)

    for key, blist in grouped.items():
        if len(blist) < 2:
            continue

        # Determine dominant strand
        plus_count = sum(1 for b in blist if b["strand"] == "+")
        minus_count = sum(1 for b in blist if b["strand"] == "-")

        # If all same strand, no inversions
        if plus_count == 0 or minus_count == 0:
            continue

        dominant = "+" if plus_count >= minus_count else "-"
        # Inversion strand = the minority strand
        inv_strand = "-" if dominant == "+" else "+"

        # Find runs of consecutive inversion-strand blocks
        i = 0
        while i < len(blist):
            if blist[i]["strand"] != inv_strand:
                i += 1
                continue

            run_start = i
            while i < len(blist) and blist[i]["strand"] == inv_strand:
                i += 1

            run = blist[run_start:i]
            if len(run) < min_blocks:
                continue

            ref_start = int(run[0]["upperStart"])
            ref_end = int(run[-1]["upperEnd"])
            coord_span = ref_end - ref_start

            # Tier-1: coordinate span >= L
            if coord_span >= min_size:
                pass  # candidate
            else:
                # Tier-2: sum of alnLen >= L
                inv_aln_sum = sum(int(b.get("alnLen", 0) or 0) for b in run)
                if inv_aln_sum >= min_size:
                    pass  # candidate
                else:
                    continue  # skip

            confidence = min(1.0, len(run) / 10.0 + 0.3)

            events.append(SVEvent(
                sv_type="INV",
                pair=key[0],
                ref_chr=key[1],
                ref_pos=ref_start,
                qry_chr=key[2],
                qry_pos=int(run[0]["lowerStart"]),
                size=coord_span,
                confidence=confidence,
                support_blocks=list(run),
                description=(
                    f"INV: {key[1]}→{key[2]}, {len(run)} blocks, "
                    f"coord_span={coord_span:,}bp, strand={inv_strand}"
                ),
            ))

    return events


def detect_fragment_inversions(blocks, min_blocks=5, min_size=1_000_000):
    """
    Detect fragment-level inversion signals.

    These are non-dominant-strand blocks/runs within a chromosome pair that prove
    a local reverse orientation, but do not satisfy the large inversion criteria
    because they are too short or too discontinuous.
    """
    events = []
    grouped = _sort_blocks_by_chr_pair(blocks)

    for key, blist in grouped.items():
        if len(blist) < 2:
            continue

        plus_count = sum(1 for b in blist if b["strand"] == "+")
        minus_count = sum(1 for b in blist if b["strand"] == "-")
        if plus_count == 0 or minus_count == 0:
            continue

        dominant = "+" if plus_count >= minus_count else "-"
        frag_strand = "-" if dominant == "+" else "+"

        i = 0
        while i < len(blist):
            if blist[i]["strand"] != frag_strand:
                i += 1
                continue

            run_start = i
            while i < len(blist) and blist[i]["strand"] == frag_strand:
                i += 1

            run = blist[run_start:i]
            ref_start = int(run[0]["upperStart"])
            ref_end = int(run[-1]["upperEnd"])
            coord_span = max(0, ref_end - ref_start)
            aln_sum = sum(int(b.get("alnLen", 0) or 0) for b in run)

            is_large_inv = (
                len(run) >= min_blocks
                and (coord_span >= min_size or aln_sum >= min_size)
            )
            if is_large_inv:
                continue

            confidence = min(0.85, max(0.25, len(run) / max(1, min_blocks) * 0.45 + 0.25))
            events.append(SVEvent(
                sv_type="FRAG_INV",
                pair=key[0],
                ref_chr=key[1],
                ref_pos=ref_start,
                qry_chr=key[2],
                qry_pos=int(run[0]["lowerStart"]),
                size=max(coord_span, aln_sum),
                confidence=confidence,
                support_blocks=list(run),
                description=(
                    f"FRAG_INV: {key[1]}→{key[2]}, {len(run)} reverse-orientation blocks, "
                    f"coord_span={coord_span:,}bp, aln_sum={aln_sum:,}bp, strand={frag_strand}"
                ),
            ))

    return events


def _summarize_blocks(vals):
    if not vals:
        return {
            "block_count": 0,
            "aln_sum": 0,
            "upper_span": 0,
            "lower_span": 0,
            "upper_start": 0,
            "upper_end": 0,
            "lower_start": 0,
            "lower_end": 0,
        }
    upper_start = min(int(b.get("upperStart", 0) or 0) for b in vals)
    upper_end = max(int(b.get("upperEnd", 0) or 0) for b in vals)
    lower_start = min(int(b.get("lowerStart", 0) or 0) for b in vals)
    lower_end = max(int(b.get("lowerEnd", 0) or 0) for b in vals)
    return {
        "block_count": len(vals),
        "aln_sum": sum(int(b.get("alnLen", 0) or 0) for b in vals),
        "upper_span": max(0, upper_end - upper_start),
        "lower_span": max(0, lower_end - lower_start),
        "upper_start": upper_start,
        "upper_end": upper_end,
        "lower_start": lower_start,
        "lower_end": lower_end,
    }


def detect_fusions(blocks, min_support=10, min_size=1_000_000,
                   min_partner_ratio=0.20, min_upper_span_fraction=0.05):
    """
    Detect strong chromosome fusion/translocation-like partner switches.

    Earlier versions treated any persistent lower-chromosome switch along an
    upper chromosome as FUSION. That produces many false positives from weak
    secondary hits, scaffold fragmentation, repeats, and residual noise.

    Current criteria are stricter:
      1. Adjacent runs on the upper chromosome must switch lower partners.
      2. Both sides of the switch must have enough blocks and aligned length.
      3. Both partners must be major partners relative to the dominant partner.
      4. Both partner segments must span a substantial part of the upper chr.
      5. The two partners together must explain most of the upper-chr support.
    """
    events = []
    grouped = defaultdict(list)
    for b in blocks:
        grouped[(b["pair"], b["upper"])].append(b)

    for key, blist in grouped.items():
        if len(blist) < min_support:
            continue

        blist.sort(key=lambda x: int(x["upperStart"]))

        partner_blocks = defaultdict(list)
        for b in blist:
            partner_blocks[b["lower"]].append(b)
        partner_stats = {
            lower: _summarize_blocks(vals)
            for lower, vals in partner_blocks.items()
        }
        if len(partner_stats) < 2:
            continue

        total_aln_all = sum(x["aln_sum"] for x in partner_stats.values())
        dominant_aln = max(x["aln_sum"] for x in partner_stats.values())
        upper_len = max(int(b.get("upperLen", 0) or 0) for b in blist)
        min_upper_span = max(min_size, int(upper_len * float(min_upper_span_fraction or 0.0)))

        major_partners = {
            lower for lower, stat in partner_stats.items()
            if stat["block_count"] >= min_support
            and stat["aln_sum"] >= min_size
            and stat["aln_sum"] >= dominant_aln * float(min_partner_ratio or 0.0)
            and stat["upper_span"] >= min_upper_span
        }
        if len(major_partners) < 2:
            continue

        runs = []
        run_start = 0
        for i in range(1, len(blist) + 1):
            if i == len(blist) or blist[i]["lower"] != blist[run_start]["lower"]:
                run = blist[run_start:i]
                runs.append({
                    "lower": blist[run_start]["lower"],
                    "start_index": run_start,
                    "blocks": run,
                    "stat": _summarize_blocks(run),
                })
                run_start = i

        for prev_run, new_run in zip(runs, runs[1:]):
            prev_chr = prev_run["lower"]
            new_chr = new_run["lower"]

            if prev_chr not in major_partners or new_chr not in major_partners:
                continue

            prev_stat = prev_run["stat"]
            new_stat = new_run["stat"]
            if prev_stat["block_count"] < min_support or new_stat["block_count"] < min_support:
                continue
            if prev_stat["aln_sum"] < min_size or new_stat["aln_sum"] < min_size:
                continue
            if prev_stat["upper_span"] < min_upper_span or new_stat["upper_span"] < min_upper_span:
                continue

            partner_pair_aln = (
                partner_stats[prev_chr]["aln_sum"]
                + partner_stats[new_chr]["aln_sum"]
            )
            if total_aln_all > 0 and partner_pair_aln / total_aln_all < 0.60:
                continue

            support_blocks = list(prev_run["blocks"]) + list(new_run["blocks"])
            support = min(prev_stat["block_count"], new_stat["block_count"])
            total_aln = prev_stat["aln_sum"] + new_stat["aln_sum"]
            ref_pos = int(new_run["blocks"][0]["upperStart"])
            confidence = min(
                1.0,
                0.35
                + min(0.30, support / 40.0)
                + min(0.20, total_aln / max(1, 10 * min_size))
                + min(0.15, partner_pair_aln / max(1, total_aln_all) * 0.15),
            )

            events.append(SVEvent(
                sv_type="FUSION",
                pair=key[0],
                ref_chr=key[1],
                ref_pos=ref_pos,
                qry_chr=new_chr,
                qry_pos=int(new_run["blocks"][0]["lowerStart"]),
                size=total_aln,
                confidence=confidence,
                support_blocks=support_blocks,
                description=(
                    f"FUSION: {key[1]} switches {prev_chr}→{new_chr} at {ref_pos:,}; "
                    f"support={prev_stat['block_count']}+{new_stat['block_count']} blocks, "
                    f"aln_sum={total_aln:,}bp, "
                    f"partner_ratio={partner_stats[new_chr]['aln_sum'] / max(1, dominant_aln):.3f}"
                ),
            ))

    return events


def detect_large_gaps(blocks, min_gap=1_000_000):
    """
    Detect large insertion/deletion events from alignment gaps.

    Strand-aware: for "-" strand blocks, query coordinates are reversed,
    so qry_gap is computed as prev.lowerStart - curr.lowerEnd instead.
    """
    events = []
    grouped = _sort_blocks_by_chr_pair(blocks)

    for key, blist in grouped.items():
        if len(blist) < 2:
            continue

        for i in range(1, len(blist)):
            prev = blist[i - 1]
            curr = blist[i]

            ref_gap = int(curr["upperStart"]) - int(prev["upperEnd"])

            # Strand-aware query gap
            strand = curr.get("strand", "+")
            if strand == "-":
                qry_gap = int(prev["lowerStart"]) - int(curr["lowerEnd"])
            else:
                qry_gap = int(curr["lowerStart"]) - int(prev["lowerEnd"])

            if ref_gap <= 0 and qry_gap <= 0:
                continue

            max_gap = max(ref_gap, qry_gap)

            if max_gap < min_gap:
                continue

            # The larger gap side indicates which genome has the insertion/deletion
            # ref_gap > qry_gap: ref has extra sequence -> DEL (deleted in query)
            # qry_gap > ref_gap: query has extra sequence -> INS (inserted in query)
            if ref_gap > qry_gap:
                sv_type = "DEL"
                size = ref_gap
            else:
                sv_type = "INS"
                size = qry_gap

            ref_pos = int(prev["upperEnd"])
            qry_pos = int(prev["lowerEnd"]) if strand != "-" else int(prev["lowerStart"])
            confidence = min(1.0, size / min_gap * 0.3 + 0.4)

            events.append(SVEvent(
                sv_type=sv_type,
                pair=key[0],
                ref_chr=key[1],
                ref_pos=ref_pos,
                qry_chr=key[2],
                qry_pos=qry_pos,
                size=size,
                confidence=confidence,
                support_blocks=[prev, curr],
                description=f"{sv_type}: {key[1]}→{key[2]}, gap size {size:,}bp at {ref_pos:,}",
            ))

    return events


def detect_all_sv(blocks, min_inv_blocks=5, min_inv_size=1_000_000,
                  min_fusion_support=10, min_fusion_size=1_000_000,
                  min_gap_size=1_000_000):
    """Run all SV detectors and return combined results."""
    events = []
    events.extend(detect_inversions(blocks, min_blocks=min_inv_blocks, min_size=min_inv_size))
    events.extend(detect_fragment_inversions(blocks, min_blocks=min_inv_blocks, min_size=min_inv_size))
    events.extend(detect_fusions(blocks, min_support=min_fusion_support, min_size=min_fusion_size))
    events.extend(detect_large_gaps(blocks, min_gap=min_gap_size))

    events.sort(key=lambda e: (e.pair, natural_key(e.ref_chr), e.ref_pos))

    # Deduplicate overlapping events of same type
    deduped = []
    for e in events:
        is_dup = False
        for d in deduped:
            if (e.sv_type == d.sv_type and e.pair == d.pair
                    and e.ref_chr == d.ref_chr and e.qry_chr == d.qry_chr
                    and abs(e.ref_pos - d.ref_pos) < 500_000):
                is_dup = True
                break
        if not is_dup:
            deduped.append(e)

    return deduped


def annotate_blocks_with_sv_events(blocks, sv_events):
    """Attach block-level SV labels for HTML display."""
    for b in blocks:
        b["svType"] = "SYN"
        b["svDescription"] = ""
        b["svConfidence"] = ""

    for event in sv_events or []:
        sv_type = getattr(event, "sv_type", "SYN") or "SYN"
        priority = SV_PRIORITY.get(sv_type, 1)
        for b in getattr(event, "support_blocks", []) or []:
            old_type = b.get("svType", "SYN") or "SYN"
            if priority < SV_PRIORITY.get(old_type, 0):
                continue
            b["svType"] = sv_type
            b["svDescription"] = getattr(event, "description", "") or sv_type
            b["svConfidence"] = round(float(getattr(event, "confidence", 0.0) or 0.0), 4)

    return blocks
