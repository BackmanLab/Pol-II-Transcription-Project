#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH -t 12:00:00
#SBATCH --job-name=merge_euseq_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/merge_out.txt
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/merge_err.txt
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Merges per-sample hg38 EU-seq BAMs into per-condition merged
# BAMs (for IGV / merged tracks; the RPKM and scale-factor steps use
# the individual BAMs, not these).
#
# Merge groups:
#   POLR2A_8hrs_aux : BM23_S40 + BM24_S41 + BM25_S42
#   WT              : BM29_S46 + BM30_S47 + BM31_S48
#   No_EU_ctrl      : single sample (S49) -- sort/index only
#
# Run this after all align_euseq_hg38.sh array tasks complete.
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs
#
# Usage:
#   sbatch merge_euseq_bams_hg38.sh
##################################################################

module purge
module load samtools/1.14

BAM_DIR="/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/BAM"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/BAM/merged"
THREADS=10

mkdir -p "${OUT_DIR}"

##--------------------------------------------------------------## Helper

merge_bams() {
  local OUT_NAME=$1
  shift
  local BAMS=("$@")

  local MERGED="${OUT_DIR}/${OUT_NAME}_merged.bam"
  local SORTED="${OUT_DIR}/${OUT_NAME}_merged.sort.bam"

  for bam in "${BAMS[@]}"; do
    [[ -f "${bam}" ]] || { echo "*** error: BAM not found: ${bam} ***" >&2; exit 1; }
  done

  if [ "${#BAMS[@]}" -eq 1 ]; then
    echo "$(date): Single sample -- sorting ${OUT_NAME}"
    samtools sort -@ "${THREADS}" "${BAMS[0]}" -o "${SORTED}"
  else
    echo "$(date): Merging ${#BAMS[@]} BAMs for ${OUT_NAME}"
    samtools merge -f -@ "${THREADS}" "${MERGED}" "${BAMS[@]}"
    echo "$(date): Sorting merged BAM for ${OUT_NAME}"
    samtools sort -@ "${THREADS}" "${MERGED}" -o "${SORTED}"
    rm -f "${MERGED}"
  fi

  samtools index "${SORTED}"

  ## self-check: merged BAM holds exactly the reads of its inputs
  local expect=0 n got
  for bam in "${BAMS[@]}"; do n=$(samtools view -c "${bam}"); expect=$((expect + n)); done
  got=$(samtools view -c "${SORTED}")
  if [[ "${got}" -ne "${expect}" || "${got}" -eq 0 ]]; then
    echo "*** error: ${OUT_NAME} merged reads ${got} != sum of inputs ${expect} ***" >&2; exit 1
  fi
  echo "  CHECK PASS  ${OUT_NAME}: ${got} reads == sum of ${#BAMS[@]} input BAM(s)"

  samtools flagstat "${SORTED}" \
    > "${OUT_DIR}/${OUT_NAME}_merged.sort.bam.flagstat.txt"

  echo "$(date): Done -- ${SORTED}"
  echo "  Flagstat:"
  cat "${OUT_DIR}/${OUT_NAME}_merged.sort.bam.flagstat.txt"
  echo ""
}

##--------------------------------------------------------------## Run merges

echo "========================================================"
echo "$(date): hg38 EU-seq BAM merge pipeline"
echo "  Input BAM dir : ${BAM_DIR}"
echo "  Output dir    : ${OUT_DIR}"
echo "========================================================"
echo ""

## POLR2A 8hrs auxin -- BM23 (S40), BM24 (S41), BM25 (S42)
merge_bams "POLR2A_8hrs_aux" \
  "${BAM_DIR}/BM23_S40_filt.sort.bam" \
  "${BAM_DIR}/BM24_S41_filt.sort.bam" \
  "${BAM_DIR}/BM25_S42_filt.sort.bam"

## WT -- BM29 (S46), BM30 (S47), BM31 (S48)
merge_bams "WT" \
  "${BAM_DIR}/BM29_S46_filt.sort.bam" \
  "${BAM_DIR}/BM30_S47_filt.sort.bam" \
  "${BAM_DIR}/BM31_S48_filt.sort.bam"

## No EU control -- single sample
merge_bams "No_EU_ctrl" \
  "${BAM_DIR}/No_EU_ctrl_S49_filt.sort.bam"

##--------------------------------------------------------------## Summary

echo "========================================================"
echo "$(date): All merges complete."
echo ""
echo "Merged BAMs:"
ls -lh "${OUT_DIR}"/*.sort.bam
echo "========================================================"
