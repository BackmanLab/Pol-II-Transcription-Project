#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Merges per-sample EU-seq BAMs into per-condition merged BAMs.
# Run this after all align_EUseq.sh jobs have completed.
#
# Merge groups:
#   POLR1A_8hrs_aux : BM20 + BM21 + BM22
#   POLR2A_8hrs_aux : BM23 + BM24 + BM25
#   ActD_6hrs       : BM27 + BM28  (BM26 excluded)
#   WT              : BM29 + BM30 + BM31
#   No_EU_ctrl      : No_EU_ctrl_S49 (single sample, sort/index only)
#
# Usage:
#   sbatch -A p32171 -p genhimem -t 04:00:00 -n 10 --mem=64G \
#     --wrap="sh merge_EUseq_bams.sh"
##################################################################

set -euo pipefail

BAM_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/BAM"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/BAM/merged"
LOG_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/logs"
THREADS=10

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

module purge
module load samtools/1.14

##--------------------------------------------------------------## Helper

merge_bams() {
  local OUT_NAME=$1
  shift
  local BAMS=("$@")

  local MERGED="${OUT_DIR}/${OUT_NAME}_merged.bam"
  local SORTED="${OUT_DIR}/${OUT_NAME}_merged.sort.bam"

  ## Verify all input BAMs exist
  for bam in "${BAMS[@]}"; do
    [[ -f "${bam}" ]] || { echo "*** error: BAM not found: ${bam} ***"; exit 1; }
  done

  if [ "${#BAMS[@]}" -eq 1 ]; then
    echo "$(date): Single sample — sorting ${OUT_NAME}"
    samtools sort -@ "${THREADS}" "${BAMS[0]}" -o "${SORTED}"
  else
    echo "$(date): Merging ${#BAMS[@]} BAMs for ${OUT_NAME}"
    samtools merge -f -@ "${THREADS}" "${MERGED}" "${BAMS[@]}"
    echo "$(date): Sorting merged BAM for ${OUT_NAME}"
    samtools sort -@ "${THREADS}" "${MERGED}" -o "${SORTED}"
    rm -f "${MERGED}"
  fi

  samtools index "${SORTED}"

  samtools flagstat "${SORTED}" \
    > "${OUT_DIR}/${OUT_NAME}_merged.sort.bam.flagstat.txt"

  echo "$(date): Done — ${SORTED}"
  echo "  Flagstat:"
  cat "${OUT_DIR}/${OUT_NAME}_merged.sort.bam.flagstat.txt"
  echo ""
}

##--------------------------------------------------------------## Run merges

echo "========================================================"
echo "$(date): EU-seq BAM merge pipeline"
echo "  Input BAM dir : ${BAM_DIR}"
echo "  Output dir    : ${OUT_DIR}"
echo "========================================================"
echo ""

## POLR1A 8hrs auxin — BM20, BM21, BM22
merge_bams "POLR1A_8hrs_aux" \
  "${BAM_DIR}/BM20_S37_filt.sort.bam" \
  "${BAM_DIR}/BM21_S38_filt.sort.bam" \
  "${BAM_DIR}/BM22_S39_filt.sort.bam"

## POLR2A 8hrs auxin — BM23, BM24, BM25
merge_bams "POLR2A_8hrs_aux" \
  "${BAM_DIR}/BM23_S40_filt.sort.bam" \
  "${BAM_DIR}/BM24_S41_filt.sort.bam" \
  "${BAM_DIR}/BM25_S42_filt.sort.bam"

## ActD 6hrs — BM27, BM28 (BM26 excluded: no nascent signal)
merge_bams "ActD_6hrs" \
  "${BAM_DIR}/BM27_S44_filt.sort.bam" \
  "${BAM_DIR}/BM28_S45_filt.sort.bam"

## WT — BM29, BM30, BM31
merge_bams "WT" \
  "${BAM_DIR}/BM29_S46_filt.sort.bam" \
  "${BAM_DIR}/BM30_S47_filt.sort.bam" \
  "${BAM_DIR}/BM31_S48_filt.sort.bam"

## No EU control — single sample
merge_bams "No_EU_ctrl" \
  "${BAM_DIR}/No_EU_ctrl_S49_filt.sort.bam"

##--------------------------------------------------------------## Summary

echo "========================================================"
echo "$(date): All merges complete."
echo ""
echo "Merged BAMs:"
ls -lh "${OUT_DIR}"/*.sort.bam
echo "========================================================"
