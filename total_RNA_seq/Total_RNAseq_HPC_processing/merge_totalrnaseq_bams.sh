#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Merges per-sample total RNA-seq BAMs into per-condition merged
# BAMs. Run this after all align_totalRNAseq.sh jobs complete.
# Merged BAMs are used to extract total mapped read counts for
# EU-seq normalization scale factor computation.
#
# Merge groups:
#   POLR1A_8hrs_aux : Rep1 + Rep2 + Rep3 + Rep4
#   POLR2A_8hrs_aux : Rep1 + Rep2 + Rep3 + Rep4
#   ActD_6hrs       : Rep1 + Rep2 + Rep3
#   WT              : Rep1 + Rep2 + Rep3
#
# Usage:
#   sbatch -A p32171 -p genhimem -t 04:00:00 -n 10 --mem=64G \
#     --wrap="sh merge_totalRNAseq_bams.sh"
##################################################################

set -euo pipefail

BAM_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/BAM"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/BAM/merged"
LOG_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/logs"
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

  ## Also print total mapped reads for easy scale factor computation
  MAPPED=$(samtools view -c -F 4 "${SORTED}")
  echo "$(date): Done — ${SORTED}"
  echo "  Total mapped reads: ${MAPPED}"
  echo "  Flagstat:"
  cat "${OUT_DIR}/${OUT_NAME}_merged.sort.bam.flagstat.txt"
  echo ""
}

##--------------------------------------------------------------## Run merges

echo "========================================================"
echo "$(date): Total RNA-seq BAM merge pipeline"
echo "  Input BAM dir : ${BAM_DIR}"
echo "  Output dir    : ${OUT_DIR}"
echo "========================================================"
echo ""

## POLR1A 8hrs auxin — Rep1, Rep2, Rep3, Rep4
merge_bams "POLR1A_8hrs_aux" \
  "${BAM_DIR}/Pol1_8hrs_Aux_Rep1_S1_filt.sort.bam" \
  "${BAM_DIR}/Pol1_8hrs_Aux_Rep2_S2_filt.sort.bam" \
  "${BAM_DIR}/Pol1_8hrs_Aux_Rep3_S3_filt.sort.bam" \
  "${BAM_DIR}/Pol1_8hrs_Aux_Rep4_S4_filt.sort.bam"

## POLR2A 8hrs auxin — Rep1, Rep2, Rep3, Rep4
merge_bams "POLR2A_8hrs_aux" \
  "${BAM_DIR}/Pol2_8hrs_Aux_Rep1_S5_filt.sort.bam" \
  "${BAM_DIR}/Pol2_8hrs_Aux_Rep2_S6_filt.sort.bam" \
  "${BAM_DIR}/Pol2_8hrs_Aux_Rep3_S7_filt.sort.bam" \
  "${BAM_DIR}/Pol2_8hrs_Aux_Rep4_S8_filt.sort.bam"

## ActD 6hrs — Rep1, Rep2, Rep3
merge_bams "ActD_6hrs" \
  "${BAM_DIR}/ActD_6hrs_Rep1_S9_filt.sort.bam" \
  "${BAM_DIR}/ActD_6hrs_Rep2_S10_filt.sort.bam" \
  "${BAM_DIR}/ActD_6hrs_Rep3_S11_filt.sort.bam"

## WT — Rep1, Rep2, Rep3
merge_bams "WT" \
  "${BAM_DIR}/WT_HCT116_CTRL_Rep1_S12_filt.sort.bam" \
  "${BAM_DIR}/WT_HCT116_CTRL_Rep2_S13_filt.sort.bam" \
  "${BAM_DIR}/WT_HCT116_CTRL_Rep3_S14_filt.sort.bam"

##--------------------------------------------------------------## Compute scale factors
## Scale factors for EU-seq normalization.
## Convention: normalize all conditions to the smallest library.
## scaleFactor = min_mapped / condition_mapped
## Apply these in bamCoverage --scaleFactor for EU-seq bigWigs.

echo "========================================================"
echo "$(date): Computing EU-seq scale factors from total RNA-seq"
echo ""

declare -A MAPPED_READS

for COND in POLR1A_8hrs_aux POLR2A_8hrs_aux ActD_6hrs WT; do
  BAM="${OUT_DIR}/${COND}_merged.sort.bam"
  MAPPED_READS[$COND]=$(samtools view -c -F 4 "${BAM}")
  echo "  ${COND}: ${MAPPED_READS[$COND]} mapped reads"
done

## Find minimum
MIN_READS=$(echo "${MAPPED_READS[@]}" | tr ' ' '\n' | sort -n | head -1)
echo ""
echo "  Minimum library size: ${MIN_READS}"
echo ""
echo "  Scale factors (min / condition):"
for COND in POLR1A_8hrs_aux POLR2A_8hrs_aux ActD_6hrs WT; do
  SF=$(echo "scale=6; ${MIN_READS} / ${MAPPED_READS[$COND]}" | bc)
  echo "    ${COND}: ${SF}"
done

echo ""
echo "Apply these scale factors in bamCoverage for EU-seq bigWigs:"
echo "  bamCoverage --bam EU_CONDITION_merged.sort.bam \\"
echo "    --scaleFactor <SF> \\"
echo "    --binSize 20 \\"
echo "    --smoothLength 60 \\"
echo "    --outFileName EU_CONDITION.totalRNA_norm.bigWig"
echo "========================================================"

##--------------------------------------------------------------## Summary

echo ""
echo "Merged BAMs:"
ls -lh "${OUT_DIR}"/*.sort.bam
