# -*- coding: utf-8 -*-
"""Command-line entry point and end-to-end workflow orchestration."""

import argparse
import sys
import unicodedata
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .autotune import auto_tune_paf_filters
from .config import load_user_config, write_config_template
from .fasta import auto_min_chr_len_for_records, filter_records_by_length, read_fasta_records
from .filters import (
    adjacent_pair_coverage,
    apply_filter_params,
    chr_block_count_rows,
    clone_args_with_target,
    filter_blocks_by_dominant_chr_pairs_with_report,
    filter_blocks_by_pair_aln,
    filter_blocks_by_pair_count,
    filter_dominant_chr_pairs_by_fill_with_rescue,
    filter_final_chr_pair_noise_by_count,
    filter_non_dominant_noise_by_dominant_pairs,
    keep_top_blocks_by_max_bottom,
    prune_chr_partners_by_block_count,
    limit_blocks,
    limit_blocks_balanced_by_pair,
    limit_blocks_per_adjacent_pair,
)
from .html_builder import build_html
from .sv import annotate_blocks_with_sv_events, detect_all_sv
from .ordering import choose_order, choose_synteny_orders, parse_turn_spec_multi
from .paf import parse_paf_for_pair, run_minimap2
from .reports import (
    write_auto_tune_report,
    write_asa_selected_steps_tsv,
    write_asa_trace_pdf,
    write_blocks_tsv,
    write_chr_pair_count_tsv,
    write_chr_pair_dominance_tsv,
    write_dominant_chr_pair_init_tsv,
    write_dropped_chr_pair_tsv,
    write_final_noise_filter_report,
    write_preset_selection_report,
    write_partner_pruning_report,
    write_dropped_partner_blocks_tsv,
    write_chr_block_counts_tsv,
    write_sv_summary_tsv,
)
from .utils import label_from_path, safe_name, split_csv


def display_label_width(label):
    width = 0
    for ch in str(label):
        width += 2 if unicodedata.east_asian_width(ch) in {"F", "W"} else 1
    return width


def auto_margin_left_for_labels(labels, minimum=220):
    longest = max((display_label_width(x) for x in labels), default=0)
    return max(minimum, int(longest * 10 + 70))


def process_adjacent_pair(
    i,
    genome_paths,
    labels,
    records_kept,
    args,
    outdir,
    prefix,
    preset_candidates_base,
    pair_target_blocks,
    use_existing_paf,
    fixed_paf_files,
    threads_per_pair,
):
    pair_label = f"{labels[i]}_vs_{labels[i + 1]}"
    pair_results = []
    preset_results = []
    initial_chr_block_count_rows = []

    print(f"[adjacent pair try] pair={i} {labels[i]} vs {labels[i + 1]}", file=sys.stderr)

    for preset_now in preset_candidates_base:
        print(f"[preset try] pair={i} preset={preset_now}", file=sys.stderr)

        if use_existing_paf:
            paf_file = Path(fixed_paf_files[i])
            if not paf_file.exists():
                raise FileNotFoundError(f"PAF not found: {paf_file}")
            print(f"[use existing PAF] pair {i}: {paf_file}", file=sys.stderr)
        else:
            # Preset must be part of the filename. Otherwise asm10/asm20 can accidentally reuse asm5 PAF.
            paf_file = outdir / (
                f"{prefix}.pair{i}_"
                f"{safe_name(labels[i])}_vs_{safe_name(labels[i + 1])}."
                f"{preset_now}.paf"
            )
            run_minimap2(
                target_fa=genome_paths[i],
                query_fa=genome_paths[i + 1],
                paf_out=paf_file,
                minimap2=args.minimap2,
                threads=threads_per_pair,
                preset=preset_now,
                force=args.force,
            )

        upper_ids = [r["seq_id"] for r in records_kept[i]]
        lower_ids = [r["seq_id"] for r in records_kept[i + 1]]

        pair_raw_blocks = parse_paf_for_pair(
            paf_file=paf_file,
            upper_ids=upper_ids,
            lower_ids=lower_ids,
            pair_index=i,
            upper_label=labels[i],
            lower_label=labels[i + 1],
            min_block=0 if args.auto_tune_paf else args.min_block,
            min_mapq=0 if args.auto_tune_paf else args.min_mapq,
            min_identity=0.0 if args.auto_tune_paf else args.min_identity,
        )

        tag = "raw blocks" if args.auto_tune_paf else "basic blocks"
        print(
            f"[pair {i} {tag}] preset={preset_now} "
            f"{labels[i]} vs {labels[i + 1]}: {len(pair_raw_blocks):,}",
            file=sys.stderr,
        )

        if args.auto_tune_paf:
            initial_params_for_report = {
                "min_block": args.min_block,
                "min_identity": args.min_identity,
                "min_mapq": args.min_mapq,
                "min_pair_aln": args.min_pair_aln,
                "min_pair_blocks": args.min_pair_blocks,
                "max_bottom_per_top": 0,
            }
            initial_before_prune = apply_filter_params(pair_raw_blocks, initial_params_for_report)
            initial_after_prune, initial_dropped_tmp, initial_partner_report_tmp = prune_chr_partners_by_block_count(
                initial_before_prune,
                max_partners=args.max_partners_per_chr,
                mode=args.partner_cap_mode,
            )
            for stage_name, block_set in [
                ("initial_before_partner_prune", initial_before_prune),
                ("initial_after_partner_prune", initial_after_prune),
            ]:
                for r in chr_block_count_rows(block_set, stage=stage_name):
                    r["preset"] = preset_now
                    r["adjacent_pair"] = i
                    r["pair_label"] = pair_label
                    initial_chr_block_count_rows.append(r)

            tune_args = clone_args_with_target(args, pair_target_blocks)
            (
                blocks_this_preset,
                best_params_this,
                tune_history_this,
                best_score_this,
                best_cov_this,
                dropped_partner_blocks_this,
                partner_report_this,
                before_prune_blocks_this,
            ) = auto_tune_paf_filters(pair_raw_blocks, tune_args)
            for row in tune_history_this:
                row["preset"] = preset_now
                row["adjacent_pair"] = i
                row["pair_label"] = pair_label
        else:
            blocks_this_preset = pair_raw_blocks
            blocks_this_preset = filter_blocks_by_pair_aln(blocks_this_preset, args.min_pair_aln)
            blocks_this_preset = filter_blocks_by_pair_count(blocks_this_preset, args.min_pair_blocks)
            blocks_this_preset = keep_top_blocks_by_max_bottom(blocks_this_preset, args.max_bottom_per_top)

            best_params_this = {
                "min_block": args.min_block,
                "min_identity": args.min_identity,
                "min_mapq": args.min_mapq,
                "min_pair_aln": args.min_pair_aln,
                "min_pair_blocks": args.min_pair_blocks,
                "max_bottom_per_top": args.max_bottom_per_top,
            }
            tune_history_this = []
            best_score_this = abs(len(blocks_this_preset) - pair_target_blocks) / max(1, pair_target_blocks)
            best_cov_this = 1.0 if blocks_this_preset else 0.0
            dropped_partner_blocks_this = []
            partner_report_this = []
            before_prune_blocks_this = blocks_this_preset

        result = {
            "adjacent_pair": i,
            "pair_label": pair_label,
            "preset": preset_now,
            "raw_blocks": pair_raw_blocks,
            "blocks": blocks_this_preset,
            "before_prune_blocks": before_prune_blocks_this,
            "dropped_partner_blocks": dropped_partner_blocks_this,
            "partner_report": partner_report_this,
            "best_params": best_params_this,
            "tune_history": tune_history_this,
            "score": best_score_this,
            "coverage": best_cov_this,
            "paf_file": str(paf_file),
            "paf_files": [str(paf_file)],
            "dominant_chr_pair_keys": best_params_this.get("dominant_chr_pair_keys", []),
            "dominant_chr_pair_records": best_params_this.get("dominant_chr_pair_records", []),
            "target_blocks": pair_target_blocks,
            "target_reached": int(best_params_this.get("target_reached", len(blocks_this_preset) >= pair_target_blocks)),
            "selected": False,
        }
        preset_results.append(result)
        pair_results.append(result)

        print(
            "[preset result] "
            f"pair={i} preset={preset_now} "
            f"blocks={len(blocks_this_preset):,} "
            f"pair_target={pair_target_blocks:,} "
            f"target_reached={int(best_params_this.get('target_reached', len(blocks_this_preset) >= pair_target_blocks))} "
            f"chr_fill_mean={best_params_this.get('chr_fill_mean', '')} "
            f"chr_fill_pass_fraction={best_params_this.get('chr_fill_pass_fraction', '')} "
            f"score={best_score_this:.6f} "
            f"coverage={best_cov_this:.4f} "
            f"min_block={best_params_this['min_block']} "
            f"min_identity={best_params_this['min_identity']} "
            f"min_mapq={best_params_this['min_mapq']} "
            f"min_pair_aln={best_params_this['min_pair_aln']} "
            f"min_pair_blocks={best_params_this['min_pair_blocks']} "
            f"max_partners_per_chr={best_params_this.get('max_partners_per_chr', args.max_partners_per_chr)} "
            f"partner_cap_mode={best_params_this.get('partner_cap_mode', args.partner_cap_mode)}",
            file=sys.stderr,
        )

        # Keep presets serial within a pair. Stop once this pair reaches target.
        if not use_existing_paf and args.preset == "auto" and int(best_params_this.get("target_reached", 0) or 0) == 1:
            print(
                f"[preset stop] pair={i} {preset_now} reached target "
                f"(chr_fill_mean={best_params_this.get('chr_fill_mean', '')})",
                file=sys.stderr,
            )
            break

    if not pair_results:
        print(f"[warning] pair={i} produced no preset results", file=sys.stderr)
        return {
            "adjacent_pair": i,
            "selected_result": None,
            "pair_results": preset_results,
            "initial_chr_block_count_rows": initial_chr_block_count_rows,
        }

    reached = [x for x in pair_results if int(x.get("target_reached", 0) or 0) == 1]
    if reached:
        chosen_pair = reached[0]
    else:
        # Preserve the existing selection behavior.
        chosen_pair = sorted(pair_results, key=lambda x: x["score"])[0]

    chosen_pair["selected"] = True

    print(
        "[pair selected] "
        f"pair={i} preset={chosen_pair['preset']} "
        f"blocks={len(chosen_pair['blocks']):,} "
        f"pair_target={pair_target_blocks:,} "
        f"target_reached={chosen_pair.get('target_reached', 0)} "
        f"score={chosen_pair['score']:.6f}",
        file=sys.stderr,
    )
    if len(chosen_pair["blocks"]) == 0:
        print(
            f"[warning] adjacent pair {i} ({labels[i]} vs {labels[i + 1]}) still has 0 blocks after all presets/filters. "
            "This pair will remain empty unless the input FASTA/PAF or filtering candidates are changed.",
            file=sys.stderr,
        )

    return {
        "adjacent_pair": i,
        "selected_result": chosen_pair,
        "pair_results": preset_results,
        "initial_chr_block_count_rows": initial_chr_block_count_rows,
    }


def prepare_genome_inputs(args):
    """
    Return genome paths and labels.

    Multi mode:
      --genomes g1.fa g2.fa g3.fa
      --genome-labels A,B,C

    Two-genome compatibility mode:
      --top-genome top.fa --bottom-genome bottom.fa
    """
    if args.genomes:
        genome_paths = [str(Path(x).expanduser().resolve()) for x in args.genomes]
        if len(genome_paths) < 2:
            raise ValueError("--genomes requires at least two FASTA files")

        if args.genome_labels:
            labels = split_csv(args.genome_labels)
            if len(labels) != len(genome_paths):
                raise ValueError("--genome-labels must have the same number of labels as --genomes")
        else:
            labels = [label_from_path(x) for x in genome_paths]

        return genome_paths, labels

    if not args.top_genome or not args.bottom_genome:
        raise ValueError("Use either --genomes g1.fa g2.fa ... or --top-genome + --bottom-genome")

    genome_paths = [
        str(Path(args.top_genome).expanduser().resolve()),
        str(Path(args.bottom_genome).expanduser().resolve()),
    ]

    labels = [
        args.top_label if args.top_label else label_from_path(args.top_genome),
        args.bottom_label if args.bottom_label else label_from_path(args.bottom_genome),
    ]
    return genome_paths, labels

def main():
    parser = argparse.ArgumentParser(
        description="Generate pure-SVG interactive synteny HTML for two or multiple genomes."
    )
    parser.add_argument("--config", default=None,
                        help="JSON config file. 'parameters' override command defaults; 'colors' override initial HTML colors.")
    parser.add_argument("--write-config-template", default=None,
                        help="Write an example JSON config file and exit.")

    # Multi-genome mode
    parser.add_argument("--genomes", nargs="+", default=None,
                        help="Multiple genome FASTA files. Adjacent genomes are aligned and drawn as rows.")
    parser.add_argument("--genome-labels", default=None,
                        help="Comma-separated labels for --genomes, same order and same count.")
    parser.add_argument("--pafs", default=None,
                        help="Comma-separated existing PAF files for adjacent genome pairs. Count must be N-1.")

    # Backward-compatible two-genome mode
    parser.add_argument("--top-genome", default=None)
    parser.add_argument("--bottom-genome", default=None)
    parser.add_argument("--paf", default=None, help="Existing PAF for two-genome mode")
    parser.add_argument("--top-label", default="Top genome")
    parser.add_argument("--bottom-label", default="Bottom genome")

    # Output
    parser.add_argument("-o", "--outdir", default="wg_html_multi_synteny")
    parser.add_argument("--prefix", default="multi_genome_synteny")
    parser.add_argument("--title", default="Multi-genome synteny")

    # Minimap2
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--preset", choices=["auto", "asm5", "asm10", "asm20"], default="auto",
                        help="Minimap2 preset. Default auto tries asm5 -> asm10 -> asm20 and keeps the first preset that reaches target blocks.")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--pair-workers", type=int, default=4,
                        help="Parallel adjacent genome pair workers")
    parser.add_argument("--tune-workers", type=int, default=8,
                        help="Parallel auto tune workers")

    # Filters
    parser.add_argument("--min-chr-len", type=int, default=3_000_000,
                        help="Manual minimum chromosome/contig length. New auto default: 3Mb. Use 0 to enable auto threshold.")
    parser.add_argument("--auto-min-chr-len-floor", type=int, default=2_000_000,
                        help="Auto contig length floor. Used only when --min-chr-len 0. Default: 2Mb.")
    parser.add_argument("--auto-min-chr-len-ratio", type=float, default=0.01,
                        help="Auto contig length ratio of the longest sequence. Used only when --min-chr-len 0. Default: 0.01.")
    parser.add_argument("--min-block", type=int, default=50_000,
                        help="Initial ASA minimum alignment block length. Default: 50000.")
    parser.add_argument("--min-mapq", type=int, default=60,
                        help="Initial ASA minimum mapping quality. Default: 60.")
    parser.add_argument("--min-identity", type=float, default=0.90,
                        help="Initial ASA minimum identity. Default: 0.90.")
    parser.add_argument("--min-pair-aln", type=int, default=10_000_000,
                        help="Initial ASA minimum total aligned length per chromosome pair. Default: 10000000.")
    parser.add_argument("--min-pair-blocks", type=int, default=0)
    parser.add_argument("--max-bottom-per-top", type=int, default=0)
    parser.add_argument("--max-blocks", type=int, default=0)

    # Auto-tuning for PAF post-processing
    parser.add_argument("--auto-tune-paf", dest="auto_tune_paf", action="store_true", default=True,
                        help="Use ASA auto logic with chromosome fill-rate target plus partner-cap pruning. Enabled by default.")
    parser.add_argument("--no-auto-tune-paf", dest="auto_tune_paf", action="store_false",
                        help="Disable automatic PAF post-filter tuning and use manual filter values.")
    parser.add_argument("--target-blocks", type=int, default=0,
                        help="Legacy total target. In new auto mode, --target-blocks-per-pair is preferred.")
    parser.add_argument("--target-blocks-per-pair", type=int, default=5000,
                        help="Legacy/reporting block target per adjacent genome pair. ASA target is now --target-chr-fill. Default: 5000.")
    parser.add_argument("--target-chr-fill", type=float, default=0.90,
                        help="Auto target: chromosome fill rate needed for a chromosome to be considered drawable. Default: 0.99.")
    parser.add_argument("--auto-rounds", type=int, default=4,
                        help="Legacy option kept for compatibility; deterministic v6 auto uses the candidate value lists instead.")
    parser.add_argument("--auto-min-block-values", default="1000,2000,5000,10000,20000,50000",
                        help="Comma-separated candidate --min-block values for auto-tuning.")
    parser.add_argument("--auto-min-identity-values", default="0.55,0.58,0.60,0.62,0.65,0.70,0.75,0.80",
                        help="Comma-separated candidate --min-identity values for auto-tuning.")
    parser.add_argument("--auto-min-mapq-values", default="0,5,10,20,30,40,60",
                        help="Comma-separated candidate --min-mapq values for auto-tuning.")
    parser.add_argument("--auto-min-pair-blocks-values", default="0,2,5,10,20,50,100",
                        help="Comma-separated candidate --min-pair-blocks values for auto-tuning.")
    parser.add_argument("--auto-max-bottom-per-top-values", default="0",
                        help="Legacy option. New auto mode uses --max-partners-per-chr instead.")
    parser.add_argument("--auto-min-pair-aln-values", default="0,50000,100000,200000,500000,1000000,2000000",
                        help="Comma-separated candidate --min-pair-aln values for deterministic auto mode.")
    parser.add_argument("--partner-cap-mode", choices=["upper", "lower", "both"], default="both",
                        help="How to apply --max-partners-per-chr. Default both means each chromosome on both sides keeps at most N partners.")
    parser.add_argument("--auto-min-coverage", type=float, default=0.65,
                        help="Soft lower bound on chromosome coverage during auto parameter selection. Default: 0.65.")
    parser.add_argument("--clean-weak-pairs", dest="clean_weak_pairs", action="store_true", default=False,
                        help="Clean tiny off-target chromosome links using dominant-pair relative filtering. Disabled by default in v6 because partner-cap pruning is the primary cleaner.")
    parser.add_argument("--no-clean-weak-pairs", dest="clean_weak_pairs", action="store_false",
                        help="Disable dominant-pair weak-link cleaning.")
    parser.add_argument("--weak-pair-ratio", type=float, default=0.05,
                        help="Drop chr-pairs whose aligned length is < this fraction of the dominant partner. Default: 0.05.")
    parser.add_argument("--max-partners-per-chr", type=int, default=4,
                        help="New auto partner cap: each chromosome keeps at most N partner chromosomes, ranked by block count. Default: 4.")
    parser.add_argument("--min-weak-pair-blocks", type=int, default=2,
                        help="Minimum block count retained for a chr-pair during clean-mode filtering. Default: 2.")
    parser.add_argument("--cap-to-target-blocks", dest="cap_to_target_blocks", action="store_true", default=False,
                        help="After auto-tuning and weak-pair pruning, keep only the longest --target-blocks blocks if the count is still higher. Disabled by default in v6.")
    parser.add_argument("--no-cap-to-target-blocks", dest="cap_to_target_blocks", action="store_false",
                        help="Disable final capping to --target-blocks.")
    parser.add_argument("--final-noise-filter", dest="final_noise_filter", action="store_true", default=True,
                        help="Before drawing, remove sparse chr-pair links by block count. Enabled by default.")
    parser.add_argument("--no-final-noise-filter", dest="final_noise_filter", action="store_false",
                        help="Disable final block-count chr-pair noise filtering.")
    parser.add_argument("--final-min-pair-blocks", type=int, default=60,
                        help="Final cleanup: drop chr-pairs with fewer than this many blocks. Default: 60.")
    parser.add_argument("--final-keep-top-partners", type=int, default=2,
                        help="Final cleanup: protect top N block-count partners for each chromosome. Default: 2.")
    parser.add_argument("--final-min-block-ratio", type=float, default=0.05,
                        help="Final cleanup: keep non-top chr-pairs only if block count reaches this fraction of the dominant partner. Default: 0.05.")
    parser.add_argument("--final-min-fill-rate", type=float, default=0.90,
                        help="Final cleanup: keep dominant chr-pairs when fill rate reaches this value. Default: 0.99.")
    parser.add_argument("--final-noise-continuity-gap", type=int, default=5_000_000,
                        help="Final dominant-pair noise filter continuity gap in bp. Default: 5000000.")
    parser.add_argument("--final-noise-borderline-ratio", type=float, default=1.25,
                        help="Final dominant-pair noise filter threshold-borderline ratio. Default: 1.25.")
    # SV detection
    parser.add_argument("--detect-sv", dest="detect_sv", action="store_true", default=True,
                        help="Detect structural variants (INV/FUSION/INS_DEL) from synteny blocks. Default: True.")
    parser.add_argument("--no-detect-sv", dest="detect_sv", action="store_false",
                        help="Disable SV detection.")
    parser.add_argument("--min-inv-blocks", type=int, default=5,
                        help="Minimum consecutive blocks for inversion detection. Default: 5.")
    parser.add_argument("--min-inv-size", type=int, default=1000000,
                        help="Minimum inversion size in bp. Default: 1000000.")
    parser.add_argument("--min-fusion-support", type=int, default=10,
                        help="Minimum support blocks for fusion detection. Default: 10.")
    parser.add_argument("--min-gap-size", type=int, default=1000000,
                        help="Minimum gap size for INS/DEL detection. Default: 1000000.")
    parser.add_argument("--final-dominant-min-blocks", type=int, default=10,
                        help="Final cleanup: rescue at least this many blocks for each dominant chr-pair if fill rate is below target. Default: 10.")
    parser.add_argument("--fill-mode-block-threshold", type=int, default=300,
                        help="Auto logic: use block-count mode at/above this many blocks; below it, use fill-rate mode. Default: 300.")

    # Display
    parser.add_argument("--order", choices=["synteny", "input", "name", "length"], default="synteny",
                        help="Chromosome order. Default synteny orders lower rows by dominant partners in the row above.")
    parser.add_argument("--turn", default="",
                        help="Display-only flipping. Examples: 'chr1,chr2' flips row 0; '0:chr1;2:chr5'; 'T2T:chr1;Pub:chr2'.")

    parser.add_argument("--width", type=int, default=3500)
    parser.add_argument("--height", type=int, default=0,
                        help="SVG height. 0 means auto based on number of genomes.")
    parser.add_argument("--margin-left", type=int, default=None,
                        help="Left margin in pixels. Default auto-expands from the longest genome label.")
    parser.add_argument("--margin-right", type=int, default=70)
    parser.add_argument("--top-y", type=int, default=150)
    parser.add_argument("--row-gap", type=int, default=260)
    parser.add_argument("--bar-height", type=int, default=24)
    parser.add_argument("--gap-ratio", type=float, default=0.012)
    parser.add_argument("--curvature", type=float, default=80.0)
    parser.add_argument("--link-opacity", type=float, default=0.58)
    parser.add_argument("--min-line-width", type=float, default=0.45)
    parser.add_argument("--max-line-width", type=float, default=3.0)
    parser.add_argument("--wide-ribbon-threshold", type=int, default=300,
                        help="Draw scaled wide ribbons only when the current visible block count is below this value. Default: 300.")

    pre_args, _ = parser.parse_known_args()
    if pre_args.write_config_template:
        template_path = write_config_template(pre_args.write_config_template)
        print(f"[config template] {template_path}", file=sys.stderr)
        return

    user_config = load_user_config(pre_args.config)
    config_parameters = dict(user_config.get("parameters", {}) or {})
    valid_dests = {
        action.dest
        for action in parser._actions
        if action.dest and action.dest != argparse.SUPPRESS
    }
    unknown_config_params = sorted(k for k in config_parameters if k not in valid_dests)
    if unknown_config_params:
        raise ValueError(
            "Unknown config parameter(s): "
            + ", ".join(unknown_config_params)
            + ". Use argparse dest names, for example min_block, target_chr_fill, width."
        )
    if config_parameters:
        parser.set_defaults(**config_parameters)
    args = parser.parse_args()
    args.color_config = user_config.get("colors", {})

    print(
        "[auto defaults] "
        f"config={args.config or 'none'} "
        f"auto_tune_paf={args.auto_tune_paf} "
        f"target_blocks={args.target_blocks} "
        f"target_blocks_per_pair={args.target_blocks_per_pair} "
        f"target_chr_fill={args.target_chr_fill} "
        f"fill_mode_block_threshold={args.fill_mode_block_threshold} "
        f"clean_weak_pairs={args.clean_weak_pairs} "
        f"weak_pair_ratio={args.weak_pair_ratio} "
        f"max_partners_per_chr={args.max_partners_per_chr} "
        f"partner_cap_mode={getattr(args, 'partner_cap_mode', 'both')} "
        f"preset={args.preset} "
        f"min_chr_len={args.min_chr_len if args.min_chr_len > 0 else 'auto'} "
        f"cap_to_target_blocks={args.cap_to_target_blocks}",
        file=sys.stderr,
    )

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = safe_name(args.prefix)
    html_file = outdir / f"{prefix}.multi_synteny.html"
    blocks_tsv = outdir / f"{prefix}.filtered_blocks.tsv"
    chr_pair_count_tsv = outdir / f"{prefix}.chr_pair_block_counts.tsv"
    chr_pair_dominance_tsv = outdir / f"{prefix}.chr_pair_dominance.tsv"
    dropped_chr_pair_tsv = outdir / f"{prefix}.dropped_chr_pairs.tsv"
    final_noise_filter_tsv = outdir / f"{prefix}.final_noise_filter_report.tsv"
    dominant_chr_pair_init_tsv = outdir / f"{prefix}.dominant_chr_pair_init.tsv"
    auto_tune_tsv = outdir / f"{prefix}.auto_tune_report.tsv"
    asa_selected_steps_tsv = outdir / f"{prefix}.asa_selected_steps.tsv"
    asa_trace_pdf = outdir / f"{prefix}.asa_trace.pdf"
    partner_pruning_tsv = outdir / f"{prefix}.partner_pruning_report.tsv"
    dropped_partner_blocks_tsv = outdir / f"{prefix}.dropped_partner_blocks.tsv"
    initial_chr_block_counts_tsv = outdir / f"{prefix}.initial_chr_block_counts.tsv"

    genome_paths, labels = prepare_genome_inputs(args)
    n = len(genome_paths)
    if args.margin_left is None:
        args.margin_left = auto_margin_left_for_labels(labels)
        print(f"[display] auto margin-left={args.margin_left}", file=sys.stderr)

    print(f"[genomes] {n}", file=sys.stderr)
    for i, (path, label) in enumerate(zip(genome_paths, labels)):
        print(f"  [{i}] {label}: {path}", file=sys.stderr)

    records_all = [read_fasta_records(path) for path in genome_paths]
    min_chr_lens = []
    records_kept = []
    for i, rs in enumerate(records_all):
        if args.min_chr_len and args.min_chr_len > 0:
            cutoff = int(args.min_chr_len)
            mode = "manual"
        else:
            cutoff = auto_min_chr_len_for_records(
                rs,
                floor=args.auto_min_chr_len_floor,
                ratio=args.auto_min_chr_len_ratio,
            )
            mode = "auto"
        kept = filter_records_by_length(rs, cutoff)
        if not kept:
            # Guardrail for fragmented assemblies: keep the longest sequence rather than stopping immediately.
            longest = max(rs, key=lambda x: x["length"])
            kept = [longest]
            print(
                f"[warning] {labels[i]}: auto/manual min_chr_len={cutoff:,} kept nothing; "
                f"fallback to longest sequence {longest['seq_id']} ({longest['length']:,} bp)",
                file=sys.stderr,
            )
        min_chr_lens.append(cutoff)
        records_kept.append(kept)
        print(
            f"[records kept] {labels[i]}: {len(kept)} / {len(rs)} "
            f"using {mode} min_chr_len={cutoff:,}",
            file=sys.stderr,
        )

    # ============================================================
    # PAF inputs / automatic minimap2 preset ladder
    # ============================================================
    # Default behavior:
    #   asm5  -> tune post-filter parameters -> if target reached, stop
    #   asm10 -> tune post-filter parameters -> if target reached, stop
    #   asm20 -> tune post-filter parameters -> final fallback
    #
    # Existing PAF mode cannot switch preset, because the PAF has already been generated.

    preset_selection_tsv = outdir / f"{prefix}.preset_selection_report.tsv"
    use_existing_paf = bool(args.pafs or (args.paf and n == 2))

    if args.pafs:
        fixed_paf_files = [str(Path(x).expanduser().resolve()) for x in split_csv(args.pafs)]
        if len(fixed_paf_files) != n - 1:
            raise ValueError("--pafs must contain exactly N-1 PAF files for --genomes")
    elif args.paf and n == 2:
        fixed_paf_files = [str(Path(args.paf).expanduser().resolve())]
    else:
        fixed_paf_files = None

    if use_existing_paf:
        preset_candidates_base = ["existing_paf"]
    else:
        preset_candidates_base = ["asm5", "asm10", "asm20"] if args.preset == "auto" else [args.preset]

    n_adjacent_pairs = max(1, n - 1)
    pair_target_blocks = (
        int(args.target_blocks_per_pair)
        if getattr(args, "target_blocks_per_pair", 0) and int(args.target_blocks_per_pair) > 0
        else max(1, (int(args.target_blocks) + n_adjacent_pairs - 1) // n_adjacent_pairs)
    )

    print(
        f"[pairwise tuning] total_target={args.target_blocks:,} "
        f"adjacent_pairs={n_adjacent_pairs} pair_target={pair_target_blocks:,}",
        file=sys.stderr,
    )

    preset_results = []
    selected_presets = []
    raw_blocks = []
    all_blocks = []
    all_before_prune_blocks = []
    tune_history = []
    paf_files = []
    all_dropped_partner_blocks = []
    all_partner_pruning_rows = []
    all_initial_chr_block_count_rows = []
    all_dominant_chr_pair_keys = set()
    all_dominant_chr_pair_records = []
    pair_filter_params = {}

    # Tune and select minimap2 preset independently for each adjacent genome pair.
    # This prevents one easy pair, e.g. G1-G2, from consuming all target blocks and leaving G3/G4 empty.
    pair_workers = max(1, min(int(args.pair_workers), n_adjacent_pairs))
    threads_per_pair = max(1, int(args.threads) // pair_workers)
    print(
        f"[pair parallel] workers={pair_workers} total_threads={args.threads} "
        f"threads_per_pair={threads_per_pair}",
        file=sys.stderr,
    )

    pair_outputs = []
    with ProcessPoolExecutor(max_workers=pair_workers) as executor:
        futures = [
            executor.submit(
                process_adjacent_pair,
                i,
                genome_paths,
                labels,
                records_kept,
                args,
                outdir,
                prefix,
                preset_candidates_base,
                pair_target_blocks,
                use_existing_paf,
                fixed_paf_files,
                threads_per_pair,
            )
            for i in range(n - 1)
        ]
        for future in as_completed(futures):
            pair_outputs.append(future.result())

    pair_outputs = sorted(pair_outputs, key=lambda x: x["adjacent_pair"])

    for pair_output in pair_outputs:
        preset_results.extend(pair_output.get("pair_results", []))
        all_initial_chr_block_count_rows.extend(pair_output.get("initial_chr_block_count_rows", []))

        chosen_pair = pair_output.get("selected_result")
        if not chosen_pair:
            continue

        i = int(chosen_pair["adjacent_pair"])
        pair_label = chosen_pair["pair_label"]
        selected_presets.append(f"pair{i}:{chosen_pair['preset']}")
        raw_blocks.extend(chosen_pair["raw_blocks"])
        all_blocks.extend(chosen_pair["blocks"])
        all_before_prune_blocks.extend(chosen_pair.get("before_prune_blocks", chosen_pair["blocks"]))
        tune_history.extend(chosen_pair["tune_history"])
        pair_filter_params[i] = dict(chosen_pair.get("best_params", {}) or {})
        for key in chosen_pair.get("dominant_chr_pair_keys", []):
            if len(key) >= 3:
                all_dominant_chr_pair_keys.add((int(key[0]), str(key[1]), str(key[2])))
        for r in chosen_pair.get("dominant_chr_pair_records", []):
            rr = dict(r)
            rr["adjacent_pair"] = i
            rr["pair_label"] = pair_label
            rr["preset"] = chosen_pair["preset"]
            all_dominant_chr_pair_records.append(rr)

        for b in chosen_pair.get("dropped_partner_blocks", []):
            b["preset"] = chosen_pair["preset"]
            b["adjacent_pair"] = i
            b["pair_label"] = pair_label
            all_dropped_partner_blocks.append(b)
        for r in chosen_pair.get("partner_report", []):
            r["preset"] = chosen_pair["preset"]
            r["adjacent_pair"] = i
            r["pair_label"] = pair_label
            all_partner_pruning_rows.append(r)
        paf_files.append(chosen_pair["paf_file"])

    selected_preset = ";".join(selected_presets) if selected_presets else "none"
    best_params = {}
    best_score = 0.0
    best_cov = 0.0

    write_preset_selection_report(preset_results, selected_preset, preset_selection_tsv)
    pair_cov, pair_counts = adjacent_pair_coverage(all_blocks, n)
    print(
        "[preset selected] "
        f"selected={selected_preset} "
        f"blocks={len(all_blocks):,} "
        f"total_target={args.target_blocks:,} "
        f"pair_target={pair_target_blocks:,} "
        f"adjacent_pair_coverage={pair_cov:.4f} "
        f"pair_counts={pair_counts}",
        file=sys.stderr,
    )
    print(f"[preset selection report] {preset_selection_tsv}", file=sys.stderr)

    if args.auto_tune_paf:
        write_auto_tune_report(tune_history, auto_tune_tsv)
        write_asa_selected_steps_tsv(tune_history, asa_selected_steps_tsv)
        if write_asa_trace_pdf(tune_history, asa_trace_pdf):
            print(f"[asa trace pdf] {asa_trace_pdf}", file=sys.stderr)
        else:
            print(f"[asa trace pdf] failed; see {asa_trace_pdf.with_suffix('.plot_error.txt')}", file=sys.stderr)
        print(
            "[auto best] "
            f"mode=pairwise "
            f"selected={selected_preset} "
            f"blocks={len(all_blocks):,} "
            f"total_target={args.target_blocks:,} "
            f"pair_target={pair_target_blocks:,} "
            f"adjacent_pair_coverage={pair_cov:.4f}",
            file=sys.stderr,
        )
        print(f"[auto report] {auto_tune_tsv}", file=sys.stderr)
        print(f"[asa selected steps] {asa_selected_steps_tsv}", file=sys.stderr)
        write_partner_pruning_report(all_partner_pruning_rows, partner_pruning_tsv)
        write_dropped_partner_blocks_tsv(all_dropped_partner_blocks, dropped_partner_blocks_tsv)
        write_chr_block_counts_tsv(all_initial_chr_block_count_rows, initial_chr_block_counts_tsv)
        write_dominant_chr_pair_init_tsv(all_dominant_chr_pair_records, dominant_chr_pair_init_tsv)
        print(f"[partner pruning report] {partner_pruning_tsv}", file=sys.stderr)
        print(f"[dropped partner blocks] {dropped_partner_blocks_tsv}", file=sys.stderr)
        print(f"[initial chr block counts] {initial_chr_block_counts_tsv}", file=sys.stderr)
        print(f"[dominant chr pair init] {dominant_chr_pair_init_tsv}", file=sys.stderr)
    else:
        print(f"[all blocks after manual filters] {len(all_blocks):,}", file=sys.stderr)

    dropped_chr_pairs = []
    if args.clean_weak_pairs:
        before_clean = len(all_blocks)
        all_blocks, dropped_chr_pairs, clean_info = filter_blocks_by_dominant_chr_pairs_with_report(
            all_blocks,
            weak_pair_ratio=args.weak_pair_ratio,
            max_partners_per_chr=args.max_partners_per_chr,
            min_pair_blocks=args.min_weak_pair_blocks,
            never_filter_to_zero=True,
        )
        print(
            f"[blocks after dominant-pair clean] {len(all_blocks):,} / {before_clean:,}; "
            f"status={clean_info.get('status')} "
            f"kept_pairs={clean_info.get('kept_pairs')} "
            f"dropped_pairs={clean_info.get('dropped_pairs')} "
            f"weak_pair_ratio={args.weak_pair_ratio} "
            f"max_partners_per_chr={args.max_partners_per_chr}",
            file=sys.stderr,
        )
        write_dropped_chr_pair_tsv(dropped_chr_pairs, dropped_chr_pair_tsv)
        print(f"[dropped weak chr-pairs] {dropped_chr_pair_tsv}", file=sys.stderr)
    else:
        print("[dominant-pair clean] disabled", file=sys.stderr)

    if args.final_noise_filter:
        before_final_noise = len(all_blocks)
        all_blocks, final_noise_rows, final_noise_info = filter_non_dominant_noise_by_dominant_pairs(
            all_blocks,
            dominant_chr_pair_keys=all_dominant_chr_pair_keys,
            pair_params=pair_filter_params,
            continuity_gap_bp=args.final_noise_continuity_gap,
            borderline_ratio=args.final_noise_borderline_ratio,
            never_filter_to_zero=True,
        )
        write_final_noise_filter_report(final_noise_rows, final_noise_filter_tsv)
        print(
            f"[final noise filter] {len(all_blocks):,} / {before_final_noise:,}; "
            f"status={final_noise_info.get('status')} "
            f"kept_pairs={final_noise_info.get('kept_pairs')} "
            f"dropped_pairs={final_noise_info.get('dropped_pairs')} "
            f"dominant_pairs={len(all_dominant_chr_pair_keys)} "
            f"continuity_gap={args.final_noise_continuity_gap} "
            f"borderline_ratio={args.final_noise_borderline_ratio}",
            file=sys.stderr,
        )
        print(f"[final noise filter report] {final_noise_filter_tsv}", file=sys.stderr)
    else:
        print("[final noise filter] disabled", file=sys.stderr)

    if args.auto_tune_paf and args.cap_to_target_blocks and len(all_blocks) > args.target_blocks:
        before_cap = len(all_blocks)
        # First cap each adjacent pair to the pair target, then enforce the total target while preserving pair coverage.
        all_blocks = limit_blocks_per_adjacent_pair(all_blocks, pair_target_blocks)
        all_blocks = limit_blocks_balanced_by_pair(all_blocks, args.target_blocks, n)
        pair_cov_after_cap, pair_counts_after_cap = adjacent_pair_coverage(all_blocks, n)
        print(
            f"[blocks after balanced cap-to-target] {len(all_blocks):,} / {before_cap:,}; "
            f"target={args.target_blocks:,} pair_target={pair_target_blocks:,} "
            f"adjacent_pair_coverage={pair_cov_after_cap:.4f} pair_counts={pair_counts_after_cap}",
            file=sys.stderr,
        )
    else:
        all_blocks = limit_blocks(all_blocks, args.max_blocks)
        pair_cov_drawn, pair_counts_drawn = adjacent_pair_coverage(all_blocks, n)
        print(
            f"[blocks drawn] {len(all_blocks):,}; "
            f"adjacent_pair_coverage={pair_cov_drawn:.4f} pair_counts={pair_counts_drawn}",
            file=sys.stderr,
        )

    write_chr_pair_count_tsv(all_blocks, chr_pair_count_tsv)
    write_chr_pair_dominance_tsv(all_blocks, chr_pair_dominance_tsv)
    print(f"[chr-pair block counts] {chr_pair_count_tsv}", file=sys.stderr)
    print(f"[chr-pair dominance] {chr_pair_dominance_tsv}", file=sys.stderr)

    if not all_blocks:
        raise RuntimeError("No blocks passed filters.")

    # Matched IDs per row
    matched = [set() for _ in range(n)]
    for b in all_blocks:
        matched[b["pair"]].add(b["upper"])
        matched[b["pair"] + 1].add(b["lower"])

    if args.order == "synteny":
        row_orders = choose_synteny_orders(records_kept, matched, all_blocks)
        print("[order] synteny-aware row ordering enabled", file=sys.stderr)
    else:
        row_orders = [choose_order(records_kept[i], matched[i], args.order) for i in range(n)]
        print(f"[order] {args.order}", file=sys.stderr)

    genomes = []
    for i in range(n):
        rec_map = {r["seq_id"]: r for r in records_kept[i]}
        row_records = [rec_map[x] for x in row_orders[i] if x in rec_map]
        genomes.append({"label": labels[i], "records": row_records})

    init_turns = parse_turn_spec_multi(args.turn, labels, n)
    # Remove unresolved chromosome names, but warn.
    clean_turns = []
    for i in range(n):
        valid = {r["seq_id"] for r in genomes[i]["records"]}
        requested = init_turns.get(i, set())
        missing = sorted(requested - valid)
        if missing:
            print(f"[warning] --turn names not found in row {i} ({labels[i]}): {','.join(missing)}", file=sys.stderr)
        clean_turns.append(requested & valid)

    if args.height and args.height > 0:
        height = args.height
    else:
        height = args.top_y + args.row_gap * (n - 1) + 170

    row_ys = [args.top_y + i * args.row_gap for i in range(n)]
    # SV detection
    sv_events = []
    if args.detect_sv:
        sv_events = detect_all_sv(
            all_blocks,
            min_inv_blocks=args.min_inv_blocks,
            min_inv_size=args.min_inv_size,
            min_fusion_support=args.min_fusion_support,
            min_gap_size=args.min_gap_size,
        )
        sv_tsv = outdir / f"{prefix}.sv_summary.tsv"
        write_sv_summary_tsv(sv_events, sv_tsv)
        print(
            f"[sv detection] {len(sv_events)} SV events detected",
            file=sys.stderr,
        )
        print(f"[sv summary] {sv_tsv}", file=sys.stderr)
    else:
        print("[sv detection] disabled", file=sys.stderr)

    all_blocks = annotate_blocks_with_sv_events(all_blocks, sv_events)
    write_blocks_tsv(all_blocks, blocks_tsv)

    html = build_html(
        genomes=genomes,
        blocks=all_blocks,
        init_turns=clean_turns,
        title=args.title,
        width=args.width,
        height=height,
        margin_left=args.margin_left,
        margin_right=args.margin_right,
        row_ys=row_ys,
        bar_height=args.bar_height,
        gap_ratio=args.gap_ratio,
        curvature=args.curvature,
        link_opacity=args.link_opacity,
        min_line_width=args.min_line_width,
        max_line_width=args.max_line_width,
        wide_ribbon_threshold=args.wide_ribbon_threshold,
        sv_events=sv_events,
        colors=args.color_config,
    )

    html_file.write_text(html, encoding="utf-8")

    print("[done]", file=sys.stderr)
    print(f"HTML: {html_file}", file=sys.stderr)
    print(f"Blocks: {blocks_tsv}", file=sys.stderr)
    print(f"Chr-pair counts: {chr_pair_count_tsv}", file=sys.stderr)
    print(f"Selected preset: {selected_preset}", file=sys.stderr)
    print(f"Preset selection report: {preset_selection_tsv}", file=sys.stderr)
    if args.auto_tune_paf:
        print(f"Auto-tune report: {auto_tune_tsv}", file=sys.stderr)
        print(f"ASA selected steps: {asa_selected_steps_tsv}", file=sys.stderr)
        print(f"ASA trace PDF: {asa_trace_pdf}", file=sys.stderr)
        print(f"Partner pruning report: {partner_pruning_tsv}", file=sys.stderr)
        print(f"Dropped partner blocks: {dropped_partner_blocks_tsv}", file=sys.stderr)
        print(f"Initial chr block counts: {initial_chr_block_counts_tsv}", file=sys.stderr)
    if args.clean_weak_pairs:
        print(f"Dropped weak chr-pairs: {dropped_chr_pair_tsv}", file=sys.stderr)
        print(f"Chr-pair dominance: {chr_pair_dominance_tsv}", file=sys.stderr)
    if args.final_noise_filter:
        print(f"Final noise filter report: {final_noise_filter_tsv}", file=sys.stderr)
    for i, paf in enumerate(paf_files):
        print(f"PAF pair {i}: {paf}", file=sys.stderr)


if __name__ == "__main__":
    main()
