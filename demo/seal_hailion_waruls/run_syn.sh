#!/bin/bash
#PBS -q high1
#PBS -l nodes=1:ppn=16
#PBS -N seal_hai_lion_walrus
#PBS -j oe

cd /public/home/wangwen_lab/duzecheng/projects/sea/genome/syn

WORKDIR="/public/home/wangwen_lab/duzecheng/projects/sea/genome/syn"
GENOME_DIR="/public/home/wangwen_lab/duzecheng/projects/sea/genome/hal_extracted_genomes"
PIPELINE_DIR="/public/home/wangwen_lab/duzecheng/dzc_pipeline/synteny/test9_small"
PIPELINE="${PIPELINE_DIR}/draw_picture_html_auto_v6.py"
OUTDIR="${WORKDIR}/seal_hai_lion_walrus"
PREFIX="seal_hai_lion_walrus"
THREADS="16"
PAIR_WORKERS="2"

source /public/home/wangwen_lab/duzecheng/soft/miniconda3/etc/profile.d/conda.sh
conda activate python310
export PATH=/public/home/wangwen_lab/duzecheng/soft/miniconda3/bin:${PATH}
export PYTHONPATH="${PIPELINE_DIR}:${PYTHONPATH:-}"

mkdir -p "${OUTDIR}"

python "${PIPELINE}" \
  --genomes \
    "${GENOME_DIR}/seal.fa" \
    "${GENOME_DIR}/haixiang.fa" \
    "${GENOME_DIR}/hailion.fa" \
  --genome-labels "seal,walrus,hai_lion" \
  -o "${OUTDIR}" \
  --prefix "${PREFIX}" \
  --title "seal walrus hai_lion synteny block10000" \
  --threads "${THREADS}" \
  --pair-workers "${PAIR_WORKERS}" \
  -block \
  --block-limit 10000 \
  --detect-sv \
  --min-inv-blocks 10 \
  --min-inv-size 1000000 \
  --min-fusion-support 10 \
  --min-gap-size 1000000 \
  2>&1 | tee "${OUTDIR}/${PREFIX}.log"
