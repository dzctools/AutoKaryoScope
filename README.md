# AutoKaryoScope

AutoKaryoScope is a Python tool for drawing interactive multi-genome synteny maps from genome FASTA files and PAF alignments. It was designed for chromosome-scale comparative genomics, especially for visualizing T2T/non-T2T genome relationships, chromosome-level collinearity, inversions, fusions, translocations, and insertion/deletion-like structural differences.

The main output is a self-contained HTML file with an interactive SVG canvas. The HTML interface supports chromosome focusing, chromosome flipping, SV coloring, manual chromosome color adjustment, and export of the current view to SVG, PNG, and PDF.

This project is actively used across our research group and will be referenced in multiple papers from our team; it is genuinely convenient for daily comparative genomics work, and we warmly welcome you to try it.

## Quick Start

Clone only the current software snapshot from GitHub. This avoids downloading old large PAF files that may exist in repository history:

```bash
git clone --depth 1 https://github.com/dzctools/AutoKaryoScope.git
cd AutoKaryoScope
```

Large demo data and example PAF files are provided separately through GitHub Releases, not through the recommended software clone.

Create the recommended conda environment:

```bash
conda env create -f environment.yml
conda activate synteny_plot
```

Install the package in editable mode:

```bash
pip install -e .
```

Run a simple two-genome example from existing PAF:

```bash
python AutoKaryoScope.py \
  --top-genome genomeA.fa \
  --bottom-genome genomeB.fa \
  --paf A_vs_B.paf \
  --top-label Genome_A \
  --bottom-label Genome_B \
  -o ./autokaryoscope_out \
  --prefix A_vs_B
```

Run a multi-genome example from existing PAF files:

```bash
python AutoKaryoScope.py \
  --genomes genome1.fa genome2.fa genome3.fa genome4.fa \
  --genome-labels G1,G2,G3,G4 \
  --pafs G1_vs_G2.paf,G2_vs_G3.paf,G3_vs_G4.paf \
  -o ./autokaryoscope_out \
  --prefix four_genomes \
  --title "Four-genome synteny"
```

For large datasets or ordinary desktop computers, use block-limited mode:

```bash
python AutoKaryoScope.py \
  --genomes genome1.fa genome2.fa genome3.fa genome4.fa \
  --genome-labels G1,G2,G3,G4 \
  --pafs G1_vs_G2.paf,G2_vs_G3.paf,G3_vs_G4.paf \
  -o ./autokaryoscope_out \
  --prefix four_genomes_block8000 \
  -block \
  --block-limit 8000
```

Open the output HTML file in a browser:

```text
autokaryoscope_out/<prefix>.multi_synteny.html
```

## Demo

Two lightweight interactive HTML demos are included directly in this repository, so users can open them without downloading the large demo data package.

### Demo 1: seal, sea lion, and walrus

Repository path:

```text
demo/seal_hailion_waruls/seal_hai_lion_walrus.multi_synteny.html
```

Online preview:

[Open demo: seal_hai_lion_walrus.multi_synteny.html](https://htmlpreview.github.io/?https://github.com/dzctools/AutoKaryoScope/blob/main/demo/seal_hailion_waruls/seal_hai_lion_walrus.multi_synteny.html)

### Demo 2: 10 animals

Repository path:

```text
demo/10animals/animal_10species_test9_block10000.multi_synteny.html
```

Online preview:

[Open demo: animal_10species_test9_block10000.multi_synteny.html](https://htmlpreview.github.io/?https://github.com/dzctools/AutoKaryoScope/blob/main/demo/10animals/animal_10species_test9_block10000.multi_synteny.html)

Large demo resources are distributed separately from the source code so that the software package remains lightweight. Download the complete demo package, including compressed PAF files and generated result reports, from GitHub Releases:

```text
https://github.com/dzctools/AutoKaryoScope/releases/tag/demo-data-v1.0.1
```

Direct demo package asset:

```text
https://github.com/dzctools/AutoKaryoScope/releases/download/demo-data-v1.0.1/AutoKaryoScope_demo.zip
```

The public assembly accessions used for the 10-animal demonstration are listed below. The order from first to last is also the adjacent PAF order used by `--pafs`.

```text
GCA_030062865.2
GCF_951394435.1
GCF_949774975.1
GCF_964270905.1
GCF_964374135.1
GCF_963455315.2
GCF_039906515.1
GCF_937001465.1
GCF_963924675.1
GCF_008692025.1
GCA_005190385.4
```

`Neophocaena_sunameri` is our newly assembled T2T genome and has not yet been published, so the genome assembly itself is not included in this repository. The seal, sea lion, and walrus genomes used for the seal/sea-lion/walrus demo are also T2T genomes from our ongoing manuscript and will be released after the manuscript is accepted.

## Recommended Workflow Before Running

### 1. For multi-genome analysis, run minimap2 yourself first

For multi-species datasets, it is recommended to generate PAF files before running AutoKaryoScope. This makes the plotting step reproducible and easier to parallelize.

Recommended minimap2 command:

```bash
minimap2 -x <preset> --secondary=no --cs -c --eqx -t <threads> target.fa query.fa > output.paf
```

Common presets:

```text
asm5   closely related assemblies
asm10  moderately diverged assemblies
asm20  more diverged assemblies
```

For `N` genomes, provide `N - 1` PAF files. The PAF files must match adjacent genome pairs in the same order:

```text
genome1.fa genome2.fa genome3.fa genome4.fa
G1_vs_G2.paf,G2_vs_G3.paf,G3_vs_G4.paf
```

The PAF file names themselves are not used to determine the genome relationship. They can have any name, but the order passed to `--pafs` must follow the adjacent genome order exactly:

```text
--genomes G1.fa G2.fa G3.fa G4.fa
--pafs   first_pair.paf,second_pair.paf,third_pair.paf

first_pair.paf  = G1 vs G2
second_pair.paf = G2 vs G3
third_pair.paf  = G3 vs G4
```

### 2. Limit block number for large HTML output

By default, AutoKaryoScope does not force a fixed final block limit. It tries to keep the best syntenic matches after optimization. For large or noisy datasets, this can produce very large HTML files.

If your computer has limited memory or the browser is slow, use:

```bash
-block --block-limit 8000
```

In block-limited mode, the algorithm first searches for a good syntenic solution, then uses rollback and cleanup steps to reduce noisy links and cap the result near the requested block count.

Recommended starting value:

```text
8000 blocks per adjacent genome pair
```

Increase this value for publication-quality dense figures, or decrease it for faster interactive browsing.

### 3. Try the demo HTML

This repository can include a demo HTML result file. Open it directly in a browser to test the interface and interactive operations:

```text
demo/*.html
```

If the demo file is hosted online, open the provided GitHub Pages or release link in your browser.

### 4. Test with example PAF files

Example PAF files used in our tests can be provided through a release, cloud storage link, or project data directory. After downloading them, run AutoKaryoScope with `--paf` or `--pafs`.

Expected layout for a multi-genome test:

```text
data/
  genome1.fa
  genome2.fa
  genome3.fa
  G1_vs_G2.paf
  G2_vs_G3.paf
```

Command:

```bash
python AutoKaryoScope.py \
  --genomes data/genome1.fa data/genome2.fa data/genome3.fa \
  --genome-labels G1,G2,G3 \
  --pafs data/G1_vs_G2.paf,data/G2_vs_G3.paf \
  -o ./test_out \
  --prefix test
```

## What This Software Can Do

AutoKaryoScope can:

- Draw two-genome and multi-genome synteny maps.
- Use existing PAF files or automatically run minimap2.
- Compare adjacent genomes in a user-defined order.
- Build chromosome-scale dominant synteny skeletons.
- Automatically optimize PAF filtering parameters with an optimizer-inspired algorithm.
- Reduce non-dominant noisy chromosome-pair links.
- Annotate inversion, fusion, translocation, insertion, and deletion-like events.
- Focus on selected chromosomes in the HTML viewer.
- Flip chromosome direction interactively.
- Adjust chromosome colors manually in the browser.
- Export the current browser view as SVG, PNG, or PDF.
- Configure default colors and parameters through JSON.

The tool is useful for chromosome-scale comparative genomics, T2T/non-T2T comparison, assembly inspection, multi-species chromosome conservation, and structural variation visualization.

## Environment

The recommended environment is provided in `environment.yml`.

```yaml
name: synteny_plot
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python>=3.8
  - minimap2
  - pip
  - pip:
    - -e .
```

Install with:

```bash
conda env create -f environment.yml
conda activate synteny_plot
```

If you already have a Python environment:

```bash
pip install -e .
```

`minimap2` is required only when the software is asked to generate PAF files automatically. If you provide `--paf` or `--pafs`, minimap2 is not required during plotting.

## Output Files

The output directory contains files such as:

```text
<prefix>.multi_synteny.html
<prefix>.filtered_blocks.tsv
<prefix>.chr_pair_block_counts.tsv
<prefix>.chromosome_correspondence.tsv
<prefix>.preset_selection_report.tsv
<prefix>.auto_tune_report.tsv
<prefix>.optimizer_selected_steps.tsv
<prefix>.optimizer_trace.pdf
<prefix>.final_noise_filter_report.tsv
<prefix>.sv_summary.tsv
```

Important files:

- `multi_synteny.html`: main interactive result.
- `filtered_blocks.tsv`: final blocks used for drawing.
- `chromosome_correspondence.tsv`: dominant chromosome correspondence table.
- `auto_tune_report.tsv`: full automatic filtering history.
- `optimizer_selected_steps.tsv`: selected optimizer and rollback steps.
- `final_noise_filter_report.tsv`: non-dominant noise filtering report.
- `sv_summary.tsv`: structural variation summary.

## HTML Viewer

The HTML viewer supports:

- Zooming and scrolling in the browser.
- Selecting chromosomes for focused display.
- Showing only the selected chromosome and its corresponding adjacent chromosome.
- Flipping chromosome direction.
- Coloring chromosomes manually.
- Showing SV-colored synteny links.
- Exporting the current view to SVG, PNG, or PDF.

The output is self-contained and can usually be opened without a web server.

## Configuration File

AutoKaryoScope supports a JSON configuration file for default parameters and colors.

Generate a template:

```bash
python AutoKaryoScope.py --write-config-template configure.example.json
```

Run with a configuration file:

```bash
python AutoKaryoScope.py \
  --config configure.json \
  --genomes genome1.fa genome2.fa genome3.fa \
  --pafs G1_vs_G2.paf,G2_vs_G3.paf \
  -o ./out
```

Example structure:

```json
{
  "parameters": {
    "width": 3500,
    "height": 0,
    "margin_left": null,
    "margin_right": 70,
    "min_chr_len": 3000000,
    "min_block": 50000,
    "min_mapq": 60,
    "min_identity": 0.9,
    "max_partners_per_chr": 4,
    "final_noise_filter": true
  },
  "colors": {
    "page_background": "#ffffff",
    "text_color": "#111827",
    "chromosome": "#1a2847",
    "link_forward": "rgba(180, 180, 180, 0.50)",
    "link_reverse": "rgba(214, 120, 95, 0.62)",
    "sv_colors": {
      "SYN": "#b0b0b0",
      "INV": "#dc2626",
      "FUSION": "#dc2626",
      "TRANS": "#dc2626",
      "INS": "#16a34a",
      "DEL": "#f59e0b"
    }
  }
}
```

Command-line arguments override values in the configuration file.

For example:

```bash
python AutoKaryoScope.py --config configure.json --width 5000 ...
```

Here `--width 5000` overrides the width stored in `configure.json`.

## Automatic Filtering Algorithm

The automatic mode is enabled by default:

```bash
--auto-tune-paf
```

The algorithm:

1. Parses raw PAF blocks.
2. Removes low-quality chromosomes before skeleton construction.
3. Builds dominant chromosome-pair skeletons.
4. Initializes filtering parameters from dominant representative blocks.
5. Uses optimizer-inspired parameter search.
6. Rolls back when noisy links increase too much.
7. Runs single-parameter cleanup near the block limit.
8. Applies final non-dominant chromosome-pair noise filtering.
9. Annotates SV events and writes the HTML viewer.

The initialization step is designed to start from reliable chromosome-scale matches, not from arbitrary local blocks.

## Manual Parameter Control

Disable automatic tuning:

```bash
--no-auto-tune-paf
```

Then manually set filtering parameters:

```bash
python AutoKaryoScope.py \
  --genomes g1.fa g2.fa g3.fa \
  --genome-labels G1,G2,G3 \
  --pafs g1_g2.paf,g2_g3.paf \
  --no-auto-tune-paf \
  --min-identity 0.90 \
  --min-block 50000 \
  --min-mapq 30 \
  --min-pair-aln 1000000 \
  --min-pair-blocks 2 \
  -o ./manual_out
```

Example strict setting:

```bash
--no-auto-tune-paf \
--min-identity 0.95 \
--min-block 500000 \
--min-mapq 50
```

Example loose setting:

```bash
--no-auto-tune-paf \
--min-identity 0.80 \
--min-block 10000 \
--min-mapq 10
```

## Parameters

### Input and output

| Parameter | Description |
|---|---|
| `--genomes` | FASTA files for multi-genome mode. |
| `--genome-labels` | Comma-separated labels for `--genomes`. |
| `--pafs` | Comma-separated PAF files for adjacent genome pairs. |
| `--top-genome` | First genome in two-genome mode. |
| `--bottom-genome` | Second genome in two-genome mode. |
| `--paf` | Existing PAF file for two-genome mode. |
| `--top-label` | Label for top genome. |
| `--bottom-label` | Label for bottom genome. |
| `-o`, `--outdir` | Output directory. |
| `--prefix` | Output file prefix. |
| `--title` | HTML title. |
| `--config` | JSON configuration file. |
| `--write-config-template` | Write a configuration template and exit. |

### Alignment

| Parameter | Description |
|---|---|
| `--threads` | Total minimap2 threads. |
| `--minimap2` | minimap2 executable path. |
| `--preset` | `auto`, `asm5`, `asm10`, or `asm20`. |
| `--force` | Regenerate minimap2 PAF files. |
| `--pair-workers` | Number of adjacent genome pairs processed in parallel. |
| `--tune-workers` | Number of workers used in parameter tuning. |

### Basic filtering

| Parameter | Description |
|---|---|
| `--min-chr-len` | Minimum sequence length retained as chromosome/scaffold. |
| `--auto-min-chr-len-floor` | Floor for automatic chromosome length filtering. |
| `--auto-min-chr-len-ratio` | Ratio for automatic chromosome length filtering. |
| `--min-block` | Minimum alignment block length. |
| `--min-mapq` | Minimum PAF MAPQ. |
| `--min-identity` | Minimum alignment identity. |
| `--min-pair-aln` | Minimum cumulative aligned length for chromosome pair. |
| `--min-pair-blocks` | Minimum block count for chromosome pair. |
| `--max-bottom-per-top` | Limit lower partners per upper chromosome. |
| `--max-blocks` | Maximum total blocks after filtering. |

### Automatic tuning and block limit

| Parameter | Description |
|---|---|
| `--auto-tune-paf` | Enable automatic PAF filtering. Default. |
| `--no-auto-tune-paf` | Disable automatic tuning and use manual thresholds. |
| `-block`, `--block` | Enable block-limited mode. |
| `--block-limit` | Block limit per adjacent genome pair in block-limited mode. |
| `--target-blocks` | Legacy compatibility option. New block-limited runs should use `-block --block-limit`. |
| `--target-blocks-per-pair` | Legacy compatibility option. New block-limited runs should use `-block --block-limit`. |
| `--target-chr-fill` | Target chromosome fill rate. |
| `--auto-rounds` | Number of automatic search rounds. |
| `--auto-min-block-values` | Candidate minimum block lengths. |
| `--auto-min-identity-values` | Candidate identities. |
| `--auto-min-mapq-values` | Candidate MAPQ values. |
| `--auto-min-pair-blocks-values` | Candidate chromosome-pair block counts. |
| `--auto-min-pair-aln-values` | Candidate chromosome-pair aligned lengths. |
| `--auto-min-coverage` | Minimum adjacent-pair coverage target. |

### Partner and noise filtering

| Parameter | Description |
|---|---|
| `--partner-cap-mode` | `upper`, `lower`, or `both`. |
| `--max-partners-per-chr` | Maximum partners retained per chromosome. |
| `--clean-weak-pairs` | Enable weak chromosome-pair cleaning. |
| `--no-clean-weak-pairs` | Disable weak chromosome-pair cleaning. |
| `--weak-pair-ratio` | Ratio threshold for weak pair removal. |
| `--min-weak-pair-blocks` | Minimum weak-pair block count. |
| `--final-noise-filter` | Enable final non-dominant noise filtering. Default. |
| `--no-final-noise-filter` | Disable final noise filtering. |
| `--final-min-pair-blocks` | Minimum pair blocks in final filtering. |
| `--final-keep-top-partners` | Number of top partners kept during final filtering. |
| `--final-min-block-ratio` | Minimum block ratio in final filtering. |
| `--final-min-fill-rate` | Minimum fill rate in final filtering. |
| `--final-noise-continuity-gap` | Continuity gap for noise filtering. |
| `--final-noise-borderline-ratio` | Borderline ratio for final noise filtering. |

### SV detection

| Parameter | Description |
|---|---|
| `--detect-sv` | Enable SV detection. Default. |
| `--no-detect-sv` | Disable SV detection. |
| `--min-inv-blocks` | Minimum blocks for inversion call. |
| `--min-inv-size` | Minimum inversion size. |
| `--min-fusion-support` | Minimum support for fusion call. |
| `--min-gap-size` | Minimum gap size for insertion/deletion-like event. |
| `--final-dominant-min-blocks` | Minimum dominant blocks for final SV-aware logic. |
| `--fill-mode-block-threshold` | Threshold for fill-mode behavior. |

### Layout and display

| Parameter | Description |
|---|---|
| `--order` | Genome/chromosome order: `synteny`, `input`, `name`, or `length`. |
| `--turn` | Manual chromosome turning/flipping specification. |
| `--width` | SVG canvas width. |
| `--height` | SVG canvas height. `0` means automatic. |
| `--margin-left` | Left margin. Useful for long species names. |
| `--margin-right` | Right margin. |
| `--top-y` | Top row y position. |
| `--row-gap` | Vertical distance between genome rows. |
| `--bar-height` | Chromosome bar height. |
| `--gap-ratio` | Gap ratio between chromosomes. |
| `--curvature` | Link curvature. |
| `--link-opacity` | Link opacity. |
| `--min-line-width` | Minimum link width. |
| `--max-line-width` | Maximum link width. |
| `--wide-ribbon-threshold` | Threshold for drawing wide ribbons. |

## Performance Tips

- Use existing PAF files for multi-genome projects.
- Use `-block --block-limit 8000` for large datasets.
- Increase `--min-chr-len` to remove small scaffolds.
- Use `--order input` when you want the plot to follow your input species order.
- Increase `--margin-left` if species names are clipped.
- For very large HTML files, reduce block limit or disable dense SV display.

## Citation

If this tool is useful for your work, please cite or acknowledge AutoKaryoScope and include the GitHub repository link.
