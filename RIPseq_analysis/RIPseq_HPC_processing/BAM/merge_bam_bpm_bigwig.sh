#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Merges BAM replicates for RIP-seq samples of interest and
# generates BPM-normalized bigWigs matching DeepTools parameters
# that produced the expected TES signal profile:
#   --normalizeUsing BPM --binSize 20 --smoothLength 60
#
# Samples:
#   POLR2A (WT P2AB): Rep2 + Rep3 merged
#   EZH2   (WT):      Rep1 + Rep2 merged
#   IgG    (WT):      single sample
#   Input  (WT):      single sample
##################################################################

set -euo pipefail

BAM_DIR="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/BAM"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/coverage_bpm"
THREADS=10

mkdir -p "${OUT_DIR}"

module purge
module load samtools/1.14
module load deeptools/3.5.1

##--------------------------------------------------------------## Helper

merge_and_bpm() {
  local OUT_NAME=$1
  shift
  local BAMS=("$@")

  local MERGED_BAM="${OUT_DIR}/${OUT_NAME}_merged.bam"
  local SORTED_BAM="${OUT_DIR}/${OUT_NAME}_merged.sort.bam"
  local BPM_BW="${OUT_DIR}/${OUT_NAME}.bpm.bigWig"

  if [ "${#BAMS[@]}" -eq 1 ]; then
    echo "$(date): Sorting single BAM for ${OUT_NAME}"
    samtools sort -@ "${THREADS}" "${BAMS[0]}" -o "${SORTED_BAM}"
  else
    echo "$(date): Merging ${#BAMS[@]} BAMs for ${OUT_NAME}"
    samtools merge -f -@ "${THREADS}" "${MERGED_BAM}" "${BAMS[@]}"
    echo "$(date): Sorting merged BAM for ${OUT_NAME}"
    samtools sort -@ "${THREADS}" "${MERGED_BAM}" -o "${SORTED_BAM}"
    rm -f "${MERGED_BAM}"
  fi

  echo "$(date): Indexing ${SORTED_BAM}"
  samtools index "${SORTED_BAM}"

  echo "$(date): Flagstat for ${OUT_NAME}"
  samtools flagstat "${SORTED_BAM}" \
    > "${OUT_DIR}/${OUT_NAME}_merged.sort.bam.flagstat.txt"

  echo "$(date): Generating BPM bigWig for ${OUT_NAME}"
  bamCoverage \
    --bam                "${SORTED_BAM}" \
    --normalizeUsing     BPM \
    --binSize            20 \
    --smoothLength       60 \
    --outFileName        "${BPM_BW}" \
    --numberOfProcessors "${THREADS}" \
    --ignoreDuplicates \
    --extendReads

  echo "$(date): Done — ${BPM_BW}"
  echo ""
}

##--------------------------------------------------------------## Run

echo "========================================================"
echo "$(date): RIP-seq BAM merge + BPM bigWig pipeline"
echo "  BAM dir: ${BAM_DIR}"
echo "  Out dir: ${OUT_DIR}"
echo "========================================================"
echo ""

## POLR2A — merge WT_P2AB Rep2 + Rep3
merge_and_bpm "POLR2A" \
  "${BAM_DIR}/WT_P2AB_Rep2_S19_L001_filt.sort.bam" \
  "${BAM_DIR}/WT_P2AB_Rep3_S20_L001_filt.sort.bam"

## EZH2 — merge WT_EZH2 Rep1 + Rep2
merge_and_bpm "EZH2" \
  "${BAM_DIR}/WT_EZH2_Rep1_S9_L001_filt.sort.bam" \
  "${BAM_DIR}/WT_EZH2_Rep2_S10_L001_filt.sort.bam"

## IgG — single sample
merge_and_bpm "IgG" \
  "${BAM_DIR}/WT_IgG_S14_L001_filt.sort.bam"

## Input — single sample
merge_and_bpm "Input" \
  "${BAM_DIR}/WT_Input_S13_L001_filt.sort.bam"

##--------------------------------------------------------------## Summary

echo "========================================================"
echo "$(date): All done."
echo ""
echo "BigWigs:"
ls -lh "${OUT_DIR}"/*.bpm.bigWig
echo ""
echo "Merged BAMs:"
ls -lh "${OUT_DIR}"/*.sort.bam
echo "========================================================"
