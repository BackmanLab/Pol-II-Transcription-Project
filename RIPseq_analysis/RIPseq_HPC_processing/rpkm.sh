#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=ripseq_rpkm
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/022025_RIPseq/logs/ripseq_rpkm_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/022025_RIPseq/logs/ripseq_rpkm_%j.err

##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# RIP-seq exon/intron RPKM quantification for Pol2 (P2AB) IP,
# EZH2 IP, WT Input, and POLR2A-KO Input — no subtraction (unlike
# the EU-seq NoEU-subtraction pipeline). IP/Input ratios and any
# replicate averaging/handling are left entirely to the downstream
# R analysis — this script only produces per-replicate RPKM
# matrices.
#
# Strandedness: reverse-stranded, CONFIRMED both empirically
# (check_strandedness.sh on WT_Input_S13: 99.7% read1-antisense-
# to-gene, n=107,493 pooled reads) and by kit chemistry (Illumina
# Stranded Total RNA Prep w/ Ribo-Zero Plus uses the dUTP
# second-strand method; Illumina's own documentation confirms
# Read 1 = first-strand cDNA = antisense to the original
# transcript). Because of this, bedgraph strand assignment below
# is FLIPPED relative to a naive read1-strand-equals-gene-strand
# approach: read1 reads aligned to the "-" strand are used to
# build "+" strand GENE coverage, and vice versa.
#
# NOTE: this is the OPPOSITE assignment used in the EU-seq
# pipeline's code (euseq_subtract_count.sh), which does NOT flip.
# That EU-seq behavior has been separately verified as correct for
# that library (check_strandedness.sh on BM29: 98.7% read1-sense-
# to-gene) — the two pipelines are on genuinely different library
# preps and are each handled correctly for their own chemistry.
#
# Exon/intron BED intervals are reused from the EU-seq project
# (090124_euseq/euseq_pipeline/intron_dog/beds), built from
# eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf — confirmed as the
# correct reference for this RIP-seq data too, so E/I interval
# definitions stay consistent across both pipelines.
#
# Chromosome naming (chr-prefixed) manually verified consistent
# across the RIP-seq BAMs, the EU-seq BED intervals, and
# hg38.chrom.sizes before this script was finalized — no runtime
# check needed.
#
# Two ratios to be computed downstream in R:
#   log2(P2AB / WT_Input)
#   log2(P2AB / P2KO_Input)
# EZH2 RPKM is also computed here for future reference/QC, but is
# not part of either ratio above.
#
# Usage:
#   sbatch ripseq_exon_intron_rpkm.sh
##################################################################

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

set -euo pipefail

module purge
module load samtools
module load bedtools

##--------------------------------------------------------------## 1. Paths and parameters

BAM_DIR="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/BAM"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/ripseq_pipeline"
CHROM_SIZES="/projects/b1042/BackmanLab/Lucas/eric_splicing/chroms/hg38.chrom.sizes"

## Reused from the EU-seq project — built from eric_splicing/gtf/
## Homo_sapiens.GRCh38.115.gtf, confirmed as the shared reference
EUSEQ_BED_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline/intron_dog/beds"
EXON_BED_PLUS="${EUSEQ_BED_DIR}/exons_plus.bed"
EXON_BED_MINUS="${EUSEQ_BED_DIR}/exons_minus.bed"
INTRON_BED_PLUS="${EUSEQ_BED_DIR}/hg38_introns_plus.bed"
INTRON_BED_MINUS="${EUSEQ_BED_DIR}/hg38_introns_minus.bed"

for f in "${EXON_BED_PLUS}" "${EXON_BED_MINUS}" "${INTRON_BED_PLUS}" "${INTRON_BED_MINUS}"; do
  if [[ ! -s "${f}" ]]; then
    echo "ERROR: Required BED not found: ${f}" >&2
    echo "       Run the EU-seq intron/exon BED-building steps first." >&2
    exit 1
  fi
done

BGR_DIR="${OUT_DIR}/bedgraphs"
EXON_CNT_DIR="${OUT_DIR}/exon_counts"
INTRON_CNT_DIR="${OUT_DIR}/intron_counts"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${BGR_DIR}" "${EXON_CNT_DIR}" "${INTRON_CNT_DIR}" "${LOG_DIR}"

##--------------------------------------------------------------## 2. Sample definitions
## No replicate averaging — P2AB and EZH2 replicates are kept as
## separate columns all the way through to the final matrices.

P2AB_BAMS=(
  "${BAM_DIR}/WT_P2AB_Rep2_S19_L001_filt.sort.bam"
  "${BAM_DIR}/WT_P2AB_Rep3_S20_L001_filt.sort.bam"
)
P2AB_NAMES=("P2AB_Rep2" "P2AB_Rep3")

EZH2_BAMS=(
  "${BAM_DIR}/WT_EZH2_Rep1_S9_L001_filt.sort.bam"
  "${BAM_DIR}/WT_EZH2_Rep2_S10_L001_filt.sort.bam"
)
EZH2_NAMES=("EZH2_Rep1" "EZH2_Rep2")

WT_INPUT_BAM="${BAM_DIR}/WT_Input_S13_L001_filt.sort.bam"
WT_INPUT_NAME="WT_Input"

P2KO_INPUT_BAM="${BAM_DIR}/P2KO_input_S12_L001_filt.sort.bam"
P2KO_INPUT_NAME="P2KO_Input"

ALL_BAMS=("${P2AB_BAMS[@]}" "${EZH2_BAMS[@]}" "${WT_INPUT_BAM}" "${P2KO_INPUT_BAM}")
ALL_NAMES=("${P2AB_NAMES[@]}" "${EZH2_NAMES[@]}" "${WT_INPUT_NAME}" "${P2KO_INPUT_NAME}")

##--------------------------------------------------------------## 3. Library sizes file

LIBSIZE_FILE="${OUT_DIR}/library_sizes.txt"
> "${LIBSIZE_FILE}"

##--------------------------------------------------------------## 4. BAM -> RPM-normalized stranded bedGraph
## No NoEU-style subtraction here — RIP-seq IP/Input correction is
## a ratio taken downstream in R after RPKM aggregation, not a
## per-base-pair subtraction like the EU-seq pipeline.
##
## Strand assignment is FLIPPED (reverse-stranded, confirmed —
## see docstring): read1 on "-" builds "+"-strand gene coverage,
## read1 on "+" builds "-"-strand gene coverage.

process_sample() {
  local bam="$1"
  local name="$2"

  local bed_file="${BGR_DIR}/${name}.bed"
  local plus_bgr="${BGR_DIR}/${name}_plus.norm.bgr"
  local minus_bgr="${BGR_DIR}/${name}_minus.norm.bgr"

  if [[ -s "${plus_bgr}" && -s "${minus_bgr}" ]]; then
    echo "$(date): ${name} bedgraphs already exist — recording library size only."
    local lib_size
    lib_size=$(samtools view -c -F 4 -f 64 "${bam}")
    echo -e "${name}\t${lib_size}" >> "${LIBSIZE_FILE}"
    return 0
  fi

  echo "$(date): Processing ${name}..."

  local lib_size
  lib_size=$(samtools view -c -F 4 -f 64 "${bam}")
  echo "  ${name} library size (read1 mapped): ${lib_size}"
  echo -e "${name}\t${lib_size}" >> "${LIBSIZE_FILE}"

  local rpm
  rpm=$(awk "BEGIN {printf \"%.10f\", 1000000.0 / ${lib_size}}")
  echo "  ${name} RPM scalar: ${rpm}"

  echo "  BAM -> BED (read1 only)..."
  samtools view -F 4 -f 64 -b "${bam}" | \
    bedtools bamtobed -i stdin | \
    sort -k1,1 -k2,2n > "${bed_file}"

  echo "  Reads in BED: $(wc -l < ${bed_file})"

  echo "  Generating RPM-normalized bedGraphs (strand-flipped, reverse-stranded)..."

  ## "+" strand GENE coverage comes from read1 reads aligned to "-"
  awk '$6 == "-"' "${bed_file}" | \
  bedtools genomecov -i stdin -g "${CHROM_SIZES}" -bg | \
  awk -v rpm="${rpm}" '{printf "%s\t%s\t%s\t%.10f\n", $1,$2,$3,$4*rpm}' | \
  sort -k1,1 -k2,2n > "${plus_bgr}"

  ## "-" strand GENE coverage comes from read1 reads aligned to "+"
  awk '$6 == "+"' "${bed_file}" | \
  bedtools genomecov -i stdin -g "${CHROM_SIZES}" -bg | \
  awk -v rpm="${rpm}" '{printf "%s\t%s\t%s\t%.10f\n", $1,$2,$3,$4*rpm}' | \
  sort -k1,1 -k2,2n > "${minus_bgr}"

  echo "  ${name} complete."
  rm -f "${bed_file}"
}

for i in "${!ALL_BAMS[@]}"; do
  process_sample "${ALL_BAMS[$i]}" "${ALL_NAMES[$i]}"
done

echo ""
echo "Library sizes:"
cat "${LIBSIZE_FILE}"

##--------------------------------------------------------------## 5. Compute exon RPKM and intron RPKM per sample
## RPKM = sum(RPM_coverage x bp) / interval_length_kb
## Aggregated across all exon (or intron) intervals per gene.

compute_interval_rpkm() {
  local name="$1"
  local interval_plus="$2"
  local interval_minus="$3"
  local out_file="$4"
  local interval_type="$5"

  if [[ -s "${out_file}" ]]; then
    echo "$(date): ${interval_type} RPKM already computed for ${name} — skipping."
    return 0
  fi

  echo "$(date): Computing ${interval_type} RPKM for ${name}..."

  local plus_bgr="${BGR_DIR}/${name}_plus.norm.bgr"
  local minus_bgr="${BGR_DIR}/${name}_minus.norm.bgr"

  if [[ ! -s "${plus_bgr}" || ! -s "${minus_bgr}" ]]; then
    echo "  ERROR: bedgraphs not found for ${name}" >&2
    exit 1
  fi

  local tmp_plus="${out_file%.txt}_plus_tmp.txt"
  local tmp_minus="${out_file%.txt}_minus_tmp.txt"

  bedtools map -a "${interval_plus}" -b "${plus_bgr}" -c 4 -o sum -null 0 | \
  awk '{
    gene_id = $4; cov_sum = $7; int_len = $3 - $2
    if (int_len <= 0) next
    print gene_id"\t"cov_sum"\t"int_len
  }' > "${tmp_plus}"

  bedtools map -a "${interval_minus}" -b "${minus_bgr}" -c 4 -o sum -null 0 | \
  awk '{
    gene_id = $4; cov_sum = $7; int_len = $3 - $2
    if (int_len <= 0) next
    print gene_id"\t"cov_sum"\t"int_len
  }' > "${tmp_minus}"

  cat "${tmp_plus}" "${tmp_minus}" | \
  sort -k1,1 | \
  awk '{
    gene = $1; cov += $2; len += $3
    if (prev != "" && gene != prev) {
      rpkm = (len > 0) ? cov / (len / 1000.0) : 0
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

  local n_signal
  n_signal=$(awk '$2 > 0' "${out_file}" | wc -l)
  echo "  Written: ${out_file}"
  echo "  Genes with ${interval_type} RPKM > 0: ${n_signal}"
}

echo ""
echo "$(date): Computing exon RPKM..."
for name in "${ALL_NAMES[@]}"; do
  compute_interval_rpkm "${name}" "${EXON_BED_PLUS}" "${EXON_BED_MINUS}" \
    "${EXON_CNT_DIR}/${name}_exon_rpkm.txt" "exon"
done

echo ""
echo "$(date): Computing intron RPKM..."
for name in "${ALL_NAMES[@]}"; do
  compute_interval_rpkm "${name}" "${INTRON_BED_PLUS}" "${INTRON_BED_MINUS}" \
    "${INTRON_CNT_DIR}/${name}_intron_rpkm.txt" "intron"
done

##--------------------------------------------------------------## 6. Merge into matrices

build_matrix() {
  local out_file="$1"
  local cnt_dir="$2"
  local suffix="$3"
  shift 3
  local names=("$@")

  local header="gene_id"
  for name in "${names[@]}"; do
    header="${header}\t${name}"
  done

  awk '{print $1}' "${cnt_dir}/${names[0]}_${suffix}.txt" > "${cnt_dir}/_genes_tmp.txt"

  local col_files=()
  for name in "${names[@]}"; do
    awk '{print $2}' "${cnt_dir}/${name}_${suffix}.txt" > "${cnt_dir}/_col_${name}_tmp.txt"
    col_files+=("${cnt_dir}/_col_${name}_tmp.txt")
  done

  paste "${cnt_dir}/_genes_tmp.txt" "${col_files[@]}" > "${out_file}.tmp"
  echo -e "${header}" | cat - "${out_file}.tmp" > "${out_file}"
  rm -f "${cnt_dir}/_genes_tmp.txt" "${out_file}.tmp"
  for name in "${names[@]}"; do
    rm -f "${cnt_dir}/_col_${name}_tmp.txt"
  done

  echo "  Written: ${out_file}"
  echo "  $(wc -l < ${out_file}) genes x ${#names[@]} samples"
}

echo ""
echo "$(date): Building RPKM matrices..."

build_matrix "${EXON_CNT_DIR}/ripseq_exon_rpkm.txt" "${EXON_CNT_DIR}" "exon_rpkm" "${ALL_NAMES[@]}"
build_matrix "${INTRON_CNT_DIR}/ripseq_intron_rpkm.txt" "${INTRON_CNT_DIR}" "intron_rpkm" "${ALL_NAMES[@]}"

echo ""
echo "========================================================"
echo "$(date): Pipeline complete."
echo ""
echo "Exon RPKM matrix:   ${EXON_CNT_DIR}/ripseq_exon_rpkm.txt"
echo "Intron RPKM matrix: ${INTRON_CNT_DIR}/ripseq_intron_rpkm.txt"
echo "Library sizes:      ${LIBSIZE_FILE}"
echo ""
echo "Matrices contain per-replicate columns (P2AB_Rep2, P2AB_Rep3,"
echo "EZH2_Rep1, EZH2_Rep2) — no replicate averaging is done here;"
echo "replicate handling (averaging, modeling separately, etc.) is"
echo "left to the R analysis step."
echo "========================================================"
