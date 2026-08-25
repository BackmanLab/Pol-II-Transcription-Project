#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=euseq_subtract
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/090124_euseq/logs/euseq_dog_int_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/090124_euseq/logs/euseq_dog_int_%j.err

##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Computes exonic RPKM directly from NoEU-subtracted bedgraphs
# by mapping coverage over exon-only intervals.
#
# This is the bedgraph approach to computing exon RPKM —
# more accurate than the count subtraction approach since
# exon signal is measured directly rather than derived as
# gene body minus intron.
#
# Reuses:
#   - Subtracted bedgraphs from euseq_pipeline.sh
#   - Merged exon BED from euseq_intron_dog.sh
#   - Library sizes from euseq_pipeline.sh
#
# Output: euseq_exon_rpkm.txt — same format as intron matrix
#
# Usage:
#   sbatch -A p32171 -p genhimem -t 02:00:00 -n 1 \
#     --cpus-per-task=4 --mem=32G \
#     euseq_exon_rpkm.sh
##################################################################

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 64G"

set -euo pipefail

module purge
module load samtools/1.6
module load bedtools/2.29.2

##--------------------------------------------------------------## 1. Paths

PIPELINE_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline"
SUB_DIR="${PIPELINE_DIR}/subtracted"
BED_DIR="${PIPELINE_DIR}/intron_dog/beds"
OUT_DIR="${PIPELINE_DIR}/exon_counts"
LIBSIZE_FILE="${PIPELINE_DIR}/library_sizes.txt"

mkdir -p "${OUT_DIR}"

##--------------------------------------------------------------## 2. Sample definitions

POLR2A_NAMES=("POLR2A_1" "POLR2A_2" "POLR2A_3")
WT_NAMES=("WT_1" "WT_2" "WT_3")
ALL_NAMES=("${POLR2A_NAMES[@]}" "${WT_NAMES[@]}")

## Exon BEDs split by strand — already built by euseq_intron_dog.sh
EXON_BED_PLUS="${BED_DIR}/exons_plus.bed"
EXON_BED_MINUS="${BED_DIR}/exons_minus.bed"

## Verify exon BEDs exist
if [[ ! -s "${EXON_BED_PLUS}" || ! -s "${EXON_BED_MINUS}" ]]; then
  echo "ERROR: Exon BED files not found." >&2
  echo "  Expected: ${EXON_BED_PLUS}" >&2
  echo "  Run euseq_intron_dog.sh first." >&2
  exit 1
fi

echo "Exon BED files found:"
echo "  Plus:  $(wc -l < ${EXON_BED_PLUS}) intervals"
echo "  Minus: $(wc -l < ${EXON_BED_MINUS}) intervals"

##--------------------------------------------------------------## 3. Compute exon RPKM per sample
## Map subtracted RPM-normalized bedgraph coverage over
## exon-only intervals. Sum coverage per exon interval,
## then aggregate per gene and divide by total exon length.
##
## RPKM = sum(RPM_coverage × bp) / total_exon_length_kb
##
## This directly measures nascent signal over exonic sequence
## without deriving it by subtraction from gene body signal.

compute_exon_rpkm() {
  local name="$1"
  local out_file="${OUT_DIR}/${name}_exon_rpkm.txt"

  if [[ -s "${out_file}" ]]; then
    echo "$(date): Exon RPKM already computed for ${name} — skipping."
    return 0
  fi

  echo "$(date): Computing exon RPKM for ${name}..."

  local plus_bgr="${SUB_DIR}/${name}_plus.sub.bgr"
  local minus_bgr="${SUB_DIR}/${name}_minus.sub.bgr"

  if [[ ! -s "${plus_bgr}" || ! -s "${minus_bgr}" ]]; then
    echo "  ERROR: subtracted bedgraphs not found for ${name}" >&2
    exit 1
  fi

  local tmp_plus="${OUT_DIR}/${name}_plus_exon_tmp.txt"
  local tmp_minus="${OUT_DIR}/${name}_minus_exon_tmp.txt"

  ## Map RPM coverage over each exon interval
  ## Output: gene_id, coverage_sum, exon_interval_length
  bedtools map \
    -a "${EXON_BED_PLUS}" \
    -b "${plus_bgr}" \
    -c 4 -o sum \
    -null 0 | \
  awk '{
    gene_id  = $4
    cov_sum  = $7
    exon_len = $3 - $2
    if (exon_len <= 0) next
    print gene_id"\t"cov_sum"\t"exon_len
  }' > "${tmp_plus}"

  bedtools map \
    -a "${EXON_BED_MINUS}" \
    -b "${minus_bgr}" \
    -c 4 -o sum \
    -null 0 | \
  awk '{
    gene_id  = $4
    cov_sum  = $7
    exon_len = $3 - $2
    if (exon_len <= 0) next
    print gene_id"\t"cov_sum"\t"exon_len
  }' > "${tmp_minus}"

  ## Aggregate across all exon intervals per gene
  ## Sum coverage and sum exon lengths, then compute RPKM
  cat "${tmp_plus}" "${tmp_minus}" | \
  sort -k1,1 | \
  awk '{
    gene  = $1
    cov  += $2
    len  += $3
    if (prev != "" && gene != prev) {
      if (len > 0)
        rpkm = cov / (len / 1000.0)
      else
        rpkm = 0
      printf "%s\t%.6f\n", prev, rpkm
      cov = $2; len = $3
    }
    prev = gene
  }
  END {
    if (prev != "" && len > 0)
      printf "%s\t%.6f\n", prev, cov / (len / 1000.0)
  }' | \
  sort -k1,1 > "${out_file}"

  rm -f "${tmp_plus}" "${tmp_minus}"

  echo "  Written: ${out_file}"
  echo "  Genes with exon RPKM > 0: $(awk '$2 > 0' ${out_file} | wc -l)"
}

## Compute exon RPKM for all EU-seq samples
for name in "${ALL_NAMES[@]}"; do
  compute_exon_rpkm "${name}"
done

##--------------------------------------------------------------## 4. Merge into matrix

echo "$(date): Building exon RPKM matrix..."

MATRIX="${OUT_DIR}/euseq_exon_rpkm.txt"

## Header
HEADER="gene_id"
for name in "${ALL_NAMES[@]}"; do
  HEADER="${HEADER}\t${name}"
done

## Gene IDs from first sample
awk '{print $1}' "${OUT_DIR}/${ALL_NAMES[0]}_exon_rpkm.txt" > \
  "${OUT_DIR}/_genes_tmp.txt"

## Extract RPKM columns
col_files=()
for name in "${ALL_NAMES[@]}"; do
  awk '{print $2}' "${OUT_DIR}/${name}_exon_rpkm.txt" > \
    "${OUT_DIR}/_col_${name}_tmp.txt"
  col_files+=("${OUT_DIR}/_col_${name}_tmp.txt")
done

paste "${OUT_DIR}/_genes_tmp.txt" "${col_files[@]}" > "${MATRIX}.tmp"
echo -e "${HEADER}" | cat - "${MATRIX}.tmp" > "${MATRIX}"
rm -f "${OUT_DIR}/_genes_tmp.txt" "${MATRIX}.tmp"
for name in "${ALL_NAMES[@]}"; do
  rm -f "${OUT_DIR}/_col_${name}_tmp.txt"
done

echo "  Written: ${MATRIX}"
echo "  $(wc -l < ${MATRIX}) genes x ${#ALL_NAMES[@]} samples"

echo ""
echo "========================================================"
echo "$(date): Exon RPKM pipeline complete."
echo "Matrix: ${MATRIX}"
echo "Next: run euseq_ei_comparison_bedgraph.R"
echo "========================================================"
