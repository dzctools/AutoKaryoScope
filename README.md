# DZC Synteny Plot

DZC Synteny Plot is a Python tool for drawing interactive multi-genome synteny maps from genome FASTA files and PAF alignments. It was designed for chromosome-scale comparative genomics, especially for visualizing T2T/non-T2T genome relationships, chromosome-level collinearity, inversions, fusions, translocations, and insertion/deletion-like structural differences.

The main output is a self-contained HTML file with an interactive SVG canvas. The HTML interface supports chromosome focusing, chromosome flipping, SV coloring, manual chromosome color adjustment, and export to SVG, PNG, and PDF.

## Features

- Multi-genome synteny visualization.
- Two-genome compatibility mode.
- Existing PAF input or automatic minimap2 alignment.
- Adaptive simulated annealing filtering for PAF blocks.
- Dominant chromosome-pair initialization.
- Noise-aware rollback and single-parameter optimization.
- Final non-dominant noise filtering.
- SV annotation for inversion, fusion, translocation, insertion/deletion-like events.
- Interactive HTML/SVG viewer.
- Export current view as SVG, PNG, or PDF.
- JSON configuration file for default parameters and initial colors.

## Installation

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate synteny_plot
```

Or install into an existing Python environment:

```bash
pip install -e .
```

The automatic alignment mode requires `minimap2`. If you already provide PAF files with `--pafs` or `--paf`, minimap2 is not required during plotting.

## Quick Start

### Multi-genome mode with existing PAF files

For `N` genomes, provide `N - 1` PAF files. Adjacent rows are connected in the same order as the genome list.

```bash
python draw_picture_html_auto_v6.py \
  --genomes genome1.fa genome2.fa genome3.fa \
  --genome-labels G1,G2,G3 \
  --pafs G1-G2.paf,G2-G3.paf \
  -o ./html \
  --prefix multi_genome_synteny \
  --title "Multi-genome synteny"
```

### Multi-genome mode with automatic minimap2 alignment

```bash
python draw_picture_html_auto_v6.py \
  --genomes genome1.fa genome2.fa genome3.fa \
  --genome-labels G1,G2,G3 \
  -o ./html \
  --prefix multi_genome_synteny \
  --threads 32
```

By default, the program tries minimap2 presets `asm5`, `asm10`, and `asm20` when PAF files are not provided.

### Two-genome mode

```bash
python draw_picture_html_auto_v6.py \
  --top-genome genomeA.fa \
  --bottom-genome genomeB.fa \
  --paf A-B.paf \
  --top-label A \
  --bottom-label B \
  -o ./html
```

## Configuration File

The software supports a JSON configuration file. You can generate a template:

```bash
python draw_picture_html_auto_v6.py --write-config-template configure.json
```

Then run with:

```bash
python draw_picture_html_auto_v6.py \
  --config configure.json \
  --genomes genome1.fa genome2.fa genome3.fa \
  --pafs G1-G2.paf,G2-G3.paf \
  -o ./html
```

The configuration file contains two major sections:

```json
{
  "parameters": {
    "width": 3500,
    "min_block": 50000,
    "min_identity": 0.9,
    "max_partners_per_chr": 4
  },
  "colors": {
    "chromosome": "#1a2847",
    "link_forward": "rgba(180, 180, 180, 0.50)",
    "link_reverse": "rgba(214, 120, 95, 0.62)",
    "sv_colors": {
      "INV": "#dc2626",
      "FUSION": "#dc2626",
      "INS_DEL": "#0d9488"
    }
  }
}
```

Values in `parameters` override program defaults. Command-line arguments have the highest priority and override the configuration file.

For example:

```bash
python draw_picture_html_auto_v6.py --config configure.json --width 5000 ...
```

In this case, `width` will be `5000`, even if `configure.json` contains another value.

## Important Parameters

Input and output:

```text
--genomes                 Multiple genome FASTA files.
--genome-labels           Comma-separated genome labels.
--pafs                    Comma-separated PAF files for adjacent genome pairs.
--outdir / -o             Output directory.
--prefix                  Output file prefix.
--title                   HTML plot title.
```

Alignment:

```text
--threads                 Total minimap2 threads.
--pair-workers            Number of adjacent genome pairs processed in parallel.
--preset                  minimap2 preset: auto, asm5, asm10, or asm20.
--force                   Regenerate minimap2 PAF files.
```

Initial filtering:

```text
--min-chr-len             Minimum chromosome/contig length.
--min-block               Initial minimum alignment block length.
--min-mapq                Initial minimum mapping quality.
--min-identity            Initial minimum alignment identity.
--min-pair-aln            Initial minimum total aligned length per chromosome pair.
```

Automatic filtering:

```text
--auto-tune-paf           Enable adaptive PAF filtering. Enabled by default.
--no-auto-tune-paf        Disable adaptive filtering.
--max-partners-per-chr    Maximum partner chromosomes kept per chromosome.
--partner-cap-mode        upper, lower, or both.
--target-chr-fill         Chromosome fill-rate target used by ASA scoring.
--final-noise-filter      Enable final non-dominant noise filtering.
```

Display:

```text
--width                   SVG canvas width.
--height                  SVG canvas height. 0 means automatic.
--margin-left             Left margin. Default automatically expands for long labels.
--row-gap                 Vertical distance between genome rows.
--bar-height              Chromosome bar height.
--link-opacity            Synteny link opacity.
--min-line-width          Minimum thin-link width.
--max-line-width          Maximum thin-link width.
--wide-ribbon-threshold   Use wide ribbons when visible block count is below this threshold.
--turn                    Initial chromosome flipping.
```

SV detection:

```text
--detect-sv               Enable SV detection. Enabled by default.
--no-detect-sv            Disable SV detection.
--min-inv-blocks          Minimum consecutive blocks for inversion detection.
--min-inv-size            Minimum inversion size.
--min-fusion-support      Minimum support blocks for fusion detection.
--min-gap-size            Minimum gap size for INS/DEL detection.
```

## HTML Viewer

The generated HTML file contains an interactive SVG viewer.

Main functions:

- Focus on one or more chromosomes.
- Flip chromosome orientation in the display.
- Toggle SV colors.
- Click chromosome bars to change their colors.
- Export the current view as SVG, PNG, or PDF.

The HTML file is self-contained and can usually be opened directly in a browser.

## Output Files

Typical output files include:

```text
<prefix>.multi_synteny.html
<prefix>.filtered_blocks.tsv
<prefix>.chr_pair_block_counts.tsv
<prefix>.chr_pair_dominance.tsv
<prefix>.auto_tune_report.tsv
<prefix>.asa_selected_steps.tsv
<prefix>.asa_trace.pdf
<prefix>.preset_selection_report.tsv
<prefix>.partner_pruning_report.tsv
<prefix>.dropped_partner_blocks.tsv
<prefix>.final_noise_filter_report.tsv
<prefix>.dominant_chr_pair_init.tsv
<prefix>.sv_summary.tsv
```

Important reports:

- `filtered_blocks.tsv`: final synteny blocks used for drawing.
- `auto_tune_report.tsv`: ASA parameter search history.
- `asa_selected_steps.tsv`: selected ASA steps.
- `asa_trace.pdf`: visual trace of ASA optimization.
- `dominant_chr_pair_init.tsv`: dominant chromosome-pair initialization records.
- `final_noise_filter_report.tsv`: final noise filtering decisions.
- `sv_summary.tsv`: detected SV events.

## Algorithm Overview

The main workflow is:

1. Read genome FASTA files and keep chromosome/contig records above the length threshold.
2. Use existing PAF files or run minimap2 for adjacent genome pairs.
3. Parse PAF alignments into synteny blocks.
4. Identify dominant chromosome-pair relationships from raw PAF alignments.
5. Build an initial chromosome-scale synteny backbone.
6. Run adaptive simulated annealing to relax filtering parameters.
7. Penalize blocks that would be removed by the final noise filter.
8. If non-dominant noise expands too much, roll back before the largest noise increase.
9. Continue with random single-parameter ASA optimization.
10. Apply final non-dominant noise filtering.
11. Detect SV events and annotate blocks.
12. Build an interactive HTML/SVG plot.

The optimization is designed to avoid relying on a manually fixed target block number. Instead, it favors recovery of dominant chromosome-scale matches while suppressing non-dominant noise.

## Example: Large Multi-Genome Plot

```bash
python draw_picture_html_auto_v6.py \
  --config configure.json \
  --genomes G1.fa G2.fa G3.fa G4.fa G5.fa \
  --genome-labels G1,G2,G3,G4,G5 \
  --pafs G1-G2.paf,G2-G3.paf,G3-G4.paf,G4-G5.paf \
  -o ./html \
  --prefix multi_genome_synteny \
  --width 3500 \
  --threads 32 \
  --pair-workers 4
```

## Recommended Git Ignore

Large genome and alignment files should not be committed to GitHub. A typical `.gitignore` should include:

```gitignore
__pycache__/
*.pyc
*.paf
*.fa
*.fasta
*.fna
*.gz
*.html
*.pdf
*.png
*.svg
*.tsv
asa_step_html/
backups/
logs/
```

## Citation

If you use this tool in a project or publication, please cite the GitHub repository and describe the version or commit hash used for analysis.

