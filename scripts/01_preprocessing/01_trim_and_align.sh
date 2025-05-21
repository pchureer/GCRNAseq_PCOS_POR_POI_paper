#!/usr/bin/env bash
#===============================================================================
# 01_trim_and_align.sh
#
# PURPOSE:
#   1) Trim adapters and low-quality bases with Trimmomatic
#   2) Align paired-end reads to GRCh38 using STAR
#   3) Produce sorted BAM and per-gene raw counts
#
# USAGE:
#   bash scripts/01_preprocessing/01_trim_and_align.sh \
#     <SAMPLE_ID> \
#     <R1.fastq.gz> \
#     <R2.fastq.gz>
#
# ARGS:
#   SAMPLE_ID       Sample identifier (e.g., SRR10239199)
#   R1.fastq.gz     Path to read-1 FASTQ file
#   R2.fastq.gz     Path to read-2 FASTQ file
#
# OUTPUTS (in results/<SAMPLE_ID>/):
#   <SAMPLE_ID>_R1_paired.fq.gz
#   <SAMPLE_ID>_R2_paired.fq.gz
#   <SAMPLE_ID>.Aligned.sortedByCoord.out.bam
#   <SAMPLE_ID>.ReadsPerGene.out.tab
#===============================================================================

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 SAMPLE_ID R1.fastq.gz R2.fastq.gz" >&2
  exit 1
fi

SAMPLE=$1
R1=$2
R2=$3

# --- Configure these paths before running ---
ADAPTERS="../data/adapters.fa"
GENOME_DIR="/path/to/STAR_hg38_v46"    # <-- update to your STAR index directory
THREADS=8
OUTDIR="../results/${SAMPLE}"
# --------------------------------------------

mkdir -p "${OUTDIR}"

echo "[`date`] Trimming adapters and low-quality bases for ${SAMPLE}"
trimmomatic PE \
  -threads ${THREADS} \
  "${R1}" "${R2}" \
  "${OUTDIR}/${SAMPLE}_R1_paired.fq.gz"  "${OUTDIR}/${SAMPLE}_R1_unpaired.fq.gz" \
  "${OUTDIR}/${SAMPLE}_R2_paired.fq.gz"  "${OUTDIR}/${SAMPLE}_R2_unpaired.fq.gz" \
  ILLUMINACLIP:${ADAPTERS}:2:30:10 \
  SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36

echo "[`date`] Aligning ${SAMPLE} to GRCh38 with STAR"
STAR \
  --runThreadN ${THREADS} \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${OUTDIR}/${SAMPLE}_R1_paired.fq.gz" "${OUTDIR}/${SAMPLE}_R2_paired.fq.gz" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUTDIR}/${SAMPLE}." \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts

echo "[`date`] Finished ${SAMPLE}. Outputs:"
ls -1 "${OUTDIR}/${SAMPLE}."*
