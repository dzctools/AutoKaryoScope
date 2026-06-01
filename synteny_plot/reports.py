# -*- coding: utf-8 -*-
"""TSV report writers."""

from pathlib import Path

from .filters import count_chr_pair_blocks, summarize_chr_pair_dominance
from .utils import natural_key

def write_blocks_tsv(blocks, path):
    with Path(path).open("w") as out:
        out.write(
            "pair\tupper_genome\tlower_genome\t"
            "upper_chr\tupper_start\tupper_end\tupper_len\t"
            "lower_chr\tlower_start\tlower_end\tlower_len\t"
            "strand\taln_len\tidentity\tmapq\tsv_type\tsv_confidence\tsv_description\n"
        )
        for b in blocks:
            out.write(
                f"{b['pair']}\t{b['upperGenome']}\t{b['lowerGenome']}\t"
                f"{b['upper']}\t{b['upperStart']}\t{b['upperEnd']}\t{b['upperLen']}\t"
                f"{b['lower']}\t{b['lowerStart']}\t{b['lowerEnd']}\t{b['lowerLen']}\t"
                f"{b['strand']}\t{b['alnLen']}\t{b['identity']:.5f}\t{b['mapq']}\t"
                f"{b.get('svType', 'SYN')}\t{b.get('svConfidence', '')}\t"
                f"{str(b.get('svDescription', '')).replace(chr(9), ' ')}\n"
            )


def write_chr_pair_dominance_tsv(blocks, path):
    stat = summarize_chr_pair_dominance(blocks)
    rows = sorted(
        stat.values(),
        key=lambda x: (x["pair"], natural_key(x["upper"]), x["upper_rank"], natural_key(x["lower"]))
    )
    fields = [
        "pair", "upperGenome", "lowerGenome", "upper", "lower",
        "block_count", "aln_sum",
        "upper_rank", "upper_ratio", "upper_count_ratio",
        "lower_rank", "lower_ratio", "lower_count_ratio",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in rows:
            row = []
            for f in fields:
                v = x.get(f, "")
                if isinstance(v, float):
                    v = f"{v:.6f}"
                row.append(str(v))
            out.write("\t".join(row) + "\n")


def write_chr_pair_count_tsv(blocks, path):
    stat = count_chr_pair_blocks(blocks)
    rows = sorted(
        stat.values(),
        key=lambda x: (x["pair"], natural_key(x["upper"]), -x["block_count"], natural_key(x["lower"]))
    )
    with Path(path).open("w") as out:
        out.write("pair\tupper_genome\tlower_genome\tupper_chr\tlower_chr\tblock_count\taln_sum\n")
        for x in rows:
            out.write(
                f"{x['pair']}\t{x['upperGenome']}\t{x['lowerGenome']}\t"
                f"{x['upper']}\t{x['lower']}\t{x['block_count']}\t{x['aln_sum']}\n"
            )


def write_dropped_chr_pair_tsv(dropped, path):
    fields = [
        "pair", "upperGenome", "lowerGenome", "upper", "lower",
        "block_count", "aln_sum",
        "upper_rank", "upper_ratio", "lower_rank", "lower_ratio",
        "reason",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in dropped:
            row = []
            for f in fields:
                v = x.get(f, "")
                if isinstance(v, float):
                    v = f"{v:.6f}"
                row.append(str(v))
            out.write("\t".join(row) + "\n")


def write_final_noise_filter_report(rows, path):
    fields = [
        "pair", "upperGenome", "lowerGenome", "upper", "lower",
        "block_count", "aln_sum",
        "current_block_count", "rescued_blocks", "blocks_after_rescue",
        "upper_fill_rate", "lower_fill_rate", "max_fill_rate",
        "min_fill_rate", "min_dominant_blocks", "dominant_partner",
        "block_protected", "fill_protected",
        "upper_block_rank", "lower_block_rank",
        "upper_partner_count", "lower_partner_count",
        "upper_max_block_count", "lower_max_block_count",
        "upper_block_ratio", "lower_block_ratio",
        "min_pair_blocks", "keep_top_partners", "min_block_ratio",
        "non_continuous", "continuity_gap_bp", "borderline_ratio",
        "min_block", "min_pair_aln", "min_identity", "min_mapq",
        "median_aln", "max_aln", "min_obs_identity", "min_obs_mapq",
        "status", "reason",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in rows:
            row = []
            for f in fields:
                v = x.get(f, "")
                if isinstance(v, float):
                    v = f"{v:.6f}"
                row.append(str(v))
            out.write("\t".join(row) + "\n")


def write_dominant_chr_pair_init_tsv(rows, path):
    fields = [
        "adjacent_pair", "pair_label", "preset",
        "pair", "upper", "lower",
        "dominant_aln_sum", "dominant_block_count",
        "representative_upper_start", "representative_upper_end",
        "representative_lower_start", "representative_lower_end",
        "representative_aln_len", "representative_identity", "representative_mapq",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in rows:
            row = []
            for f in fields:
                v = x.get(f, "")
                if isinstance(v, float):
                    v = f"{v:.6f}"
                row.append(str(v))
            out.write("\t".join(row) + "\n")


def write_auto_tune_report(history, path):
    fields = [
        "adjacent_pair", "preset", "round", "changed_param", "min_block", "min_identity", "min_mapq",
        "min_pair_aln", "min_pair_blocks", "max_bottom_per_top",
        "blocks", "score", "coverage"
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in history:
            out.write("\t".join(str(x.get(f, "")) for f in fields) + "\n")


def write_preset_selection_report(preset_results, selected_preset, path):
    fields = [
        "adjacent_pair", "pair_label", "preset", "raw_blocks", "tuned_blocks", "target_blocks", "reached_target",
        "score", "coverage", "min_block", "min_identity", "min_mapq",
        "min_pair_aln", "min_pair_blocks", "max_bottom_per_top", "selected"
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in preset_results:
            bp = x.get("best_params", {}) or {}
            target = int(x.get("target_blocks", 0) or 0)
            selected_flag = int(bool(x.get("selected", False)) or x.get("preset", "") == selected_preset)
            row = {
                "adjacent_pair": x.get("adjacent_pair", ""),
                "pair_label": x.get("pair_label", ""),
                "preset": x.get("preset", ""),
                "raw_blocks": len(x.get("raw_blocks", [])),
                "tuned_blocks": len(x.get("blocks", [])),
                "target_blocks": target,
                "reached_target": int(len(x.get("blocks", [])) >= target if target > 0 else 0),
                "score": f"{float(x.get('score', 0.0)):.6f}",
                "coverage": f"{float(x.get('coverage', 0.0)):.6f}",
                "min_block": bp.get("min_block", ""),
                "min_identity": bp.get("min_identity", ""),
                "min_mapq": bp.get("min_mapq", ""),
                "min_pair_aln": bp.get("min_pair_aln", ""),
                "min_pair_blocks": bp.get("min_pair_blocks", ""),
                "max_bottom_per_top": bp.get("max_bottom_per_top", ""),
                "selected": selected_flag,
            }
            out.write("\t".join(str(row.get(f, "")) for f in fields) + "\n")


def write_partner_pruning_report(rows, path):
    fields = [
        "adjacent_pair", "preset", "pair_label", "pair", "upperGenome", "lowerGenome",
        "upper", "lower", "block_count", "aln_sum",
        "upper_block_rank", "lower_block_rank",
        "max_partners_per_chr", "partner_cap_mode", "status", "reason",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in rows:
            row = []
            for f in fields:
                v = x.get(f, "")
                if isinstance(v, float):
                    v = f"{v:.6f}"
                row.append(str(v))
            out.write("\t".join(row) + "\n")


def write_dropped_partner_blocks_tsv(blocks, path):
    fields = [
        "adjacent_pair", "preset", "pair_label", "pair",
        "upperGenome", "lowerGenome",
        "upper", "upperStart", "upperEnd", "upperLen",
        "lower", "lowerStart", "lowerEnd", "lowerLen",
        "strand", "alnLen", "identity", "mapq",
        "upper_block_rank", "lower_block_rank", "chr_pair_block_count", "chr_pair_aln_sum",
        "source_block_index_after_filters", "reason",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for b in blocks:
            row = []
            for f in fields:
                v = b.get(f, "")
                if isinstance(v, float):
                    v = f"{v:.6f}"
                row.append(str(v))
            out.write("\t".join(row) + "\n")


def write_chr_block_counts_tsv(rows, path):
    fields = [
        "adjacent_pair", "preset", "pair_label", "stage", "row_role",
        "genome", "chr", "block_count", "partner_count", "aln_sum",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in rows:
            out.write("\t".join(str(x.get(f, "")) for f in fields) + "\n")

# Override older writer with an expanded schema for deterministic auto logic.
def write_auto_tune_report(history, path):
    fields = [
        "adjacent_pair", "pair_label", "preset", "round", "stage", "changed_param",
        "min_block", "min_identity", "min_mapq", "min_pair_aln", "min_pair_blocks",
        "max_bottom_per_top", "max_partners_per_chr", "partner_cap_mode",
        "target_blocks", "raw_blocks", "blocks_before_partner_prune",
        "dropped_by_partner_cap", "blocks", "scored_blocks", "noise_ratio",
        "noise_filter_status", "dominant_blocks", "non_dominant_blocks",
        "delta_dominant_blocks", "delta_noise_blocks",
        "delta_dominant_aln", "delta_noise_aln", "net_gain",
        "best_blocks", "best_scored_blocks", "accepted", "is_best",
        "delta", "temperature", "target_reached", "reached_target",
        "stopped", "stop_reason",
        "rescue_rounds", "rescue_action", "guard_action",
        "init_method", "init_dominant_pairs", "init_representative_blocks",
        "init_min_block", "init_min_identity", "init_min_mapq",
        "init_min_pair_aln", "init_min_pair_blocks",
        "fill_mode_block_threshold", "score_mode",
        "chr_fill_target", "chr_fill_min", "chr_fill_mean",
        "chr_fill_pass_fraction", "chr_fill_pass_count", "chr_fill_total_chr",
        "coverage", "score",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in history:
            row = []
            for f in fields:
                v = x.get(f, "")
                if isinstance(v, float):
                    v = f"{v:.6f}"
                row.append(str(v))
            out.write("\t".join(row) + "\n")


def write_asa_selected_steps_tsv(history, path):
    rows = [x for x in history if int(x.get("accepted", 0) or 0) == 1]
    return write_auto_tune_report(rows, path)


def write_asa_trace_pdf(history, path):
    path = Path(path)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except Exception as exc:
        path.with_suffix(".plot_error.txt").write_text(
            f"Cannot import matplotlib for PDF plotting: {exc}\n",
            encoding="utf-8",
        )
        return False

    if not history:
        path.with_suffix(".plot_error.txt").write_text("No ASA history rows to plot.\n", encoding="utf-8")
        return False

    grouped = {}
    for row in history:
        key = (row.get("adjacent_pair", ""), row.get("pair_label", ""), row.get("preset", ""))
        grouped.setdefault(key, []).append(row)

    with PdfPages(path) as pdf:
        for (pair, pair_label, preset), rows in sorted(grouped.items(), key=lambda x: x[0]):
            rows = sorted(rows, key=lambda x: int(x.get("round", 0) or 0))
            x = [int(r.get("round", 0) or 0) for r in rows]
            blocks = [float(r.get("blocks", 0) or 0) for r in rows]
            best_blocks = [float(r.get("best_blocks", 0) or 0) for r in rows]
            before_prune = [float(r.get("blocks_before_partner_prune", 0) or 0) for r in rows]
            accepted_x = [int(r.get("round", 0) or 0) for r in rows if int(r.get("accepted", 0) or 0) == 1]
            accepted_y = [float(r.get("blocks", 0) or 0) for r in rows if int(r.get("accepted", 0) or 0) == 1]
            target = float(rows[-1].get("target_blocks", 0) or 0)

            fig, ax = plt.subplots(figsize=(11, 7))
            ax.plot(x, before_prune, color="#9aa0a6", linewidth=1.2, label="blocks before partner prune")
            ax.plot(x, blocks, color="#1f77b4", linewidth=1.4, label="trial blocks")
            ax.plot(x, best_blocks, color="#2ca02c", linewidth=2.0, label="best/accepted blocks")
            if accepted_x:
                ax.scatter(accepted_x, accepted_y, color="#d62728", s=14, label="accepted steps")
            if target > 0:
                ax.axhline(target, color="black", linestyle="--", linewidth=1, label=f"target {int(target)}")
            ax.set_title(f"ASA optimizer block trace: pair {pair} {pair_label} {preset}")
            ax.set_xlabel("round")
            ax.set_ylabel("block count")
            ax.grid(True, alpha=0.25)
            ax.legend()
            pdf.savefig(fig)
            plt.close(fig)

            fig, ax = plt.subplots(figsize=(11, 5))
            score = [float(r.get("score", 0) or 0) for r in rows]
            coverage = [float(r.get("coverage", 0) or 0) for r in rows]
            ax.plot(x, score, color="#9467bd", label="score")
            ax.set_xlabel("round")
            ax.set_ylabel("score")
            ax.grid(True, alpha=0.25)
            ax2 = ax.twinx()
            ax2.plot(x, coverage, color="#ff7f0e", label="coverage")
            ax2.set_ylabel("coverage")
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax.legend(lines + lines2, labels + labels2, loc="best")
            ax.set_title(f"ASA score/coverage trace: pair {pair} {pair_label} {preset}")
            pdf.savefig(fig)
            plt.close(fig)

    return True


# Override older writer with additional partner-cap and before-prune columns.
def write_preset_selection_report(preset_results, selected_preset, path):
    fields = [
        "adjacent_pair", "pair_label", "preset", "raw_blocks", "blocks_before_partner_prune",
        "dropped_partner_blocks", "tuned_blocks", "target_blocks", "reached_target",
        "fill_mode_block_threshold", "score_mode",
        "chr_fill_target", "chr_fill_min", "chr_fill_mean",
        "chr_fill_pass_fraction", "chr_fill_pass_count", "chr_fill_total_chr",
        "score", "coverage", "min_block", "min_identity", "min_mapq",
        "min_pair_aln", "min_pair_blocks", "max_partners_per_chr", "partner_cap_mode", "selected",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for x in preset_results:
            bp = x.get("best_params", {}) or {}
            target = int(x.get("target_blocks", 0) or 0)
            selected_flag = int(bool(x.get("selected", False)))
            row = {
                "adjacent_pair": x.get("adjacent_pair", ""),
                "pair_label": x.get("pair_label", ""),
                "preset": x.get("preset", ""),
                "raw_blocks": len(x.get("raw_blocks", [])),
                "blocks_before_partner_prune": len(x.get("before_prune_blocks", [])),
                "dropped_partner_blocks": len(x.get("dropped_partner_blocks", [])),
                "tuned_blocks": len(x.get("blocks", [])),
                "target_blocks": target,
                "reached_target": int(x.get("target_reached", len(x.get("blocks", [])) >= target if target > 0 else 0)),
                "fill_mode_block_threshold": bp.get("fill_mode_block_threshold", ""),
                "score_mode": bp.get("score_mode", ""),
                "chr_fill_target": bp.get("target_chr_fill", bp.get("chr_fill_target", "")),
                "chr_fill_min": bp.get("chr_fill_min", ""),
                "chr_fill_mean": bp.get("chr_fill_mean", ""),
                "chr_fill_pass_fraction": bp.get("chr_fill_pass_fraction", ""),
                "chr_fill_pass_count": bp.get("chr_fill_pass_count", ""),
                "chr_fill_total_chr": bp.get("chr_fill_total_chr", ""),
                "score": f"{float(x.get('score', 0.0)):.6f}",
                "coverage": f"{float(x.get('coverage', 0.0)):.6f}",
                "min_block": bp.get("min_block", ""),
                "min_identity": bp.get("min_identity", ""),
                "min_mapq": bp.get("min_mapq", ""),
                "min_pair_aln": bp.get("min_pair_aln", ""),
                "min_pair_blocks": bp.get("min_pair_blocks", ""),
                "max_partners_per_chr": bp.get("max_partners_per_chr", ""),
                "partner_cap_mode": bp.get("partner_cap_mode", ""),
                "selected": selected_flag,
            }
            out.write("\t".join(str(row.get(f, "")) for f in fields) + "\n")


def write_sv_summary_tsv(sv_events, path):
    fields = [
        "sv_type", "pair", "ref_chr", "ref_pos", "qry_chr", "qry_pos",
        "size", "confidence", "block_count", "description",
    ]
    with Path(path).open("w") as out:
        out.write("\t".join(fields) + "\n")
        for e in sv_events:
            d = e.to_dict()
            row = [str(d.get(f, "")) for f in fields]
            out.write("\t".join(row) + "\n")
