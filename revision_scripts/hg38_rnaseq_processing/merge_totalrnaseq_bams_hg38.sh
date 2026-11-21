#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=150G
#SBATCH -t 12:00:00
#SBATCH --job-name=merge_rnaseq_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/logs/merge_out.txt
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/logs/merge_err.txt
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Merges per-sample hg38 total RNA-seq BAMs (single-end) into
# per-condition merged BAMs for generate_total_bigwigs_mito_hg38.sh:
#   POLR2A_8hrs_aux : Rep1 (S5) + Rep2 (S6) + Rep3 (S7)   (Rep4 / S8 dropped)
#   WT              : Rep1 (S12) + Rep2 (S13) + Rep3 (S14)
# Scale factors are not computed here: see
# compute_scale_factors_mito_total_hg38.sh.
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/logs
#
# Usage:
#   sbatch merge_totalrnaseq_bams_hg38.sh
##################################################################

module purge
module load samtools/1.14

BAM_DIR="/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/BAM"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/BAM/merged"
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

  echo "$(date): Merging ${#BAMS[@]} BAMs for ${OUT_NAME}"
  samtools merge -f -@ "${THREADS}" "${MERGED}" "${BAMS[@]}"
  samtools sort -@ "${THREADS}" "${MERGED}" -o "${SORTED}"
  rm -f "${MERGED}"
  samtools index "${SORTED}"
  samtools flagstat "${SORTED}" > "${OUT_DIR}/${OUT_NAME}_merged.sort.bam.flagstat.txt"

  ## self-check: merged BAM holds exactly the reads of its inputs
  local expect=0 n got
  for bam in "${BAMS[@]}"; do n=$(samtools view -c "${bam}"); expect=$((expect + n)); done
  got=$(samtools view -c "${SORTED}")
  if [[ "${got}" -ne "${expect}" || "${got}" -eq 0 ]]; then
    echo "*** error: ${OUT_NAME} merged reads ${got} != sum of inputs ${expect} ***" >&2; exit 1
  fi
  echo "  CHECK PASS  ${OUT_NAME}: ${got} reads == sum of ${#BAMS[@]} input BAMs"

  MAPPED=$(samtools view -c -F 4 "${SORTED}")
  echo "$(date): Done -- ${SORTED}"
  echo "  Total mapped reads: ${MAPPED}"
  cat "${OUT_DIR}/${OUT_NAME}_merged.sort.bam.flagstat.txt"
  echo ""
}

##--------------------------------------------------------------## Run merges

echo "========================================================"
echo "$(date): hg38 total RNA-seq BAM merge pipeline"
echo "  Input BAM dir : ${BAM_DIR}"
echo "  Output dir    : ${OUT_DIR}"
echo "========================================================"
echo ""

## POLR2A 8hrs auxin -- Rep1 (S5), Rep2 (S6), Rep3 (S7); Rep4 (S8) dropped
merge_bams "POLR2A_8hrs_aux" \
  "${BAM_DIR}/Pol2_8hrs_Aux_Rep1_S5_filt.sort.bam" \
  "${BAM_DIR}/Pol2_8hrs_Aux_Rep2_S6_filt.sort.bam" \
  "${BAM_DIR}/Pol2_8hrs_Aux_Rep3_S7_filt.sort.bam"

## WT -- Rep1 (S12), Rep2 (S13), Rep3 (S14)
merge_bams "WT" \
  "${BAM_DIR}/WT_HCT116_CTRL_Rep1_S12_filt.sort.bam" \
  "${BAM_DIR}/WT_HCT116_CTRL_Rep2_S13_filt.sort.bam" \
  "${BAM_DIR}/WT_HCT116_CTRL_Rep3_S14_filt.sort.bam"

echo "========================================================"
echo "Merged BAMs:"
ls -lh "${OUT_DIR}"/*.sort.bam
