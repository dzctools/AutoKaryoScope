# draw_picture_html_auto_v6: deterministic auto logic

This version redesigns the automatic logic for synteny plotting.

## New default auto behavior

1. Try minimap2 presets in this order:
   - `asm5`
   - `asm10`
   - `asm20`

2. The initial post-PAF filter is now:

```bash
--min-chr-len 3000000 \
--min-block 5000 \
--min-mapq 0 \
--min-identity 0.60 \
--min-pair-aln 100000
```

3. Each chromosome keeps at most 4 partner chromosomes by default:

```bash
--max-partners-per-chr 4
```

Partner ranking is based on:

```text
block_count descending, then aligned_length_sum descending
```

The default mode is bidirectional:

```bash
--partner-cap-mode both
```

Available modes:

```text
upper  each upper-row chromosome keeps top N lower partners
lower  each lower-row chromosome keeps top N upper partners
both   a chromosome pair must pass both upper-side and lower-side top-N limits
```

4. The auto target is per adjacent genome pair:

```bash
--target-blocks-per-pair 5000
```

The automatic selector chooses a parameter configuration whose final block count after partner pruning is at least this value, and is as close to the target as possible while preserving chromosome coverage.

## New reports

The following files are written in addition to the old outputs:

```text
<prefix>.auto_tune_report.tsv
<prefix>.preset_selection_report.tsv
<prefix>.partner_pruning_report.tsv
<prefix>.dropped_partner_blocks.tsv
<prefix>.initial_chr_block_counts.tsv
```

Important files:

- `auto_tune_report.tsv`: every automatic parameter attempt and its final block count.
- `initial_chr_block_counts.tsv`: per-chromosome block counts after the initial parameters, before and after partner pruning.
- `partner_pruning_report.tsv`: chromosome-pair level keep/drop decision from the top-N partner cap.
- `dropped_partner_blocks.tsv`: individual PAF-derived blocks removed by the partner cap.

## Example

```bash
python draw_picture_html_auto_v6.py \
  --genomes g1.fa g2.fa g3.fa g4.fa \
  --genome-labels G1,G2,G3,G4 \
  --title "Multi-genome synteny" \
  -o ./wg_html_multi_synteny_v6 \
  --prefix multi_genome_synteny
```

To allow more partners:

```bash
--max-partners-per-chr 6
```

To make each pair keep more blocks:

```bash
--target-blocks-per-pair 8000
```

To make the graph cleaner:

```bash
--target-blocks-per-pair 5000 --max-partners-per-chr 4 --partner-cap-mode both
```
