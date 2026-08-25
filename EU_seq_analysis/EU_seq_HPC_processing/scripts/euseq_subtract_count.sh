#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=euseq_subtract
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/090124_euseq/logs/euseq_subtract_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/090124_euseq/logs/euseq_subtract_%j.err

##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# EU-seq NoEU background subtraction and gene-level counting.
#
# Workflow (per replicate, per strand):
#   1. BAM -> bedgraph (strand-specific, rf-stranded library)
#   2. Subtract NoEU bedgraph on a base-pair basis
#   3. Floor negative values at zero
#   4. Integrate coverage over gene bodies (TSS to TES)
#      using bedtools map to get per-gene counts
#   5. Output integer count matrix for DESeq2
#
# Library: paired-end, reverse stranded (flags 83/163)
# Conditions: POLR2A (BM23-25) vs WT (BM29-31)
# NoEU control: No_EU_ctrl_S49
#
# Usage:
#   sbatch -A p32171 -p genhimem -t 06:00:00 -n 12 --mem=64G \
#     --output=.../logs/euseq_subtract_out.txt \
#     --error=.../logs/euseq_subtract_err.txt \
#     --wrap="sh euseq_subtract_count.sh"
##################################################################

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 64G"

set -euo pipefail

module purge
module load samtools/1.6
module load bedtools/2.29.2

##--------------------------------------------------------------## 1. Paths and parameters

BAM_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/BAM"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline"
GTF="/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
CHROM_SIZES="/projects/b1042/BackmanLab/Lucas/eric_splicing/chroms/hg38.chrom.sizes"

BED_DIR="${OUT_DIR}/beds"
BGR_DIR="${OUT_DIR}/bedgraphs"
SUB_DIR="${OUT_DIR}/subtracted"
CNT_DIR="${OUT_DIR}/counts"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${BED_DIR}" "${BGR_DIR}" "${SUB_DIR}" "${CNT_DIR}" "${LOG_DIR}"

THREADS=12

##--------------------------------------------------------------## 2. Sample definitions

POLR2A_BAMS=(
  "${BAM_DIR}/BM23_S40_filt.sort.bam"
  "${BAM_DIR}/BM24_S41_filt.sort.bam"
  "${BAM_DIR}/BM25_S42_filt.sort.bam"
)
POLR2A_NAMES=("POLR2A_1" "POLR2A_2" "POLR2A_3")

WT_BAMS=(
  "${BAM_DIR}/BM29_S46_filt.sort.bam"
  "${BAM_DIR}/BM30_S47_filt.sort.bam"
  "${BAM_DIR}/BM31_S48_filt.sort.bam"
)
WT_NAMES=("WT_1" "WT_2" "WT_3")

NOEU_BAM="${BAM_DIR}/No_EU_ctrl_S49_filt.sort.bam"
NOEU_NAME="NoEU"

ALL_BAMS=("${POLR2A_BAMS[@]}" "${WT_BAMS[@]}")
ALL_NAMES=("${POLR2A_NAMES[@]}" "${WT_NAMES[@]}")
ALL_NAMES_WITH_NOEU=("${ALL_NAMES[@]}" "${NOEU_NAME}")
ALL_BAMS_WITH_NOEU=("${ALL_BAMS[@]}" "${NOEU_BAM}")

##--------------------------------------------------------------## 3. Build gene body BED from GTF
## One interval per gene: TSS to TES including introns
## Strand-specific for correct read counting

GENE_BED="${OUT_DIR}/hg38_genes.bed"

if [[ ! -s "${GENE_BED}" ]]; then
  echo "$(date): Building gene body BED from GTF..."
  grep -v "^#" "${GTF}" | \
  awk '$3 == "gene" {
    gene_id = ""
    for (i = 1; i <= NF; i++) {
      if ($i == "gene_id") {
        gene_id = $(i+1)
        gsub(/"/, "", gene_id)
        gsub(/;/, "", gene_id)
        break
      }
    }
    if (gene_id == "") next
    print "chr"$1"\t"($4-1)"\t"$5"\t"gene_id"\t0\t"$7
  }' | \
  grep -E "^chr[0-9XY]+" | \
  sort -k1,1 -k2,2n > "${GENE_BED}"
  echo "  Gene BED written: $(wc -l < ${GENE_BED}) genes"
else
  echo "Gene BED already exists: ${GENE_BED}"
fi

## Separate by strand for counting
GENE_BED_PLUS="${OUT_DIR}/hg38_genes_plus.bed"
GENE_BED_MINUS="${OUT_DIR}/hg38_genes_minus.bed"
if [[ ! -s "${GENE_BED_PLUS}" ]]; then
  awk '$6 == "+"' "${GENE_BED}" > "${GENE_BED_PLUS}"
  awk '$6 == "-"' "${GENE_BED}" > "${GENE_BED_MINUS}"
  echo "  Plus genes:  $(wc -l < ${GENE_BED_PLUS})"
  echo "  Minus genes: $(wc -l < ${GENE_BED_MINUS})"
fi

##--------------------------------------------------------------## 4. Library sizes file

LIBSIZE_FILE="${OUT_DIR}/library_sizes.txt"
> "${LIBSIZE_FILE}"

##--------------------------------------------------------------## 5. BAM -> BED -> extend -> bedGraph -> RPM normalize
##
## For rf-stranded paired-end library:
##   Read1 (flag 64) is the sequenced end — its 5' position
##   indicates where RNA synthesis was occurring.
##   Read1 on + strand -> maps to + strand genes (flag 99)
##   Read1 on - strand -> maps to - strand genes (flag 83)
##
## No read extension — reads are used as aligned by HISAT2.
## HISAT2 is splice-aware so spliced reads are correctly placed.
## Intronic reads represent genuine nascent pre-mRNA signal.

process_sample() {
  local bam="$1"
  local name="$2"

  local bed_file="${BED_DIR}/${name}.bed"
  local plus_bgr="${BGR_DIR}/${name}_plus.norm.bgr"
  local minus_bgr="${BGR_DIR}/${name}_minus.norm.bgr"

  ## Skip if already processed
  if [[ -s "${plus_bgr}" && -s "${minus_bgr}" ]]; then
    echo "$(date): ${name} bedgraphs already exist — recording library size only."
    local lib_size
    lib_size=$(samtools view -c -F 4 -f 64 "${bam}")
    echo -e "${name}\t${lib_size}" >> "${LIBSIZE_FILE}"
    return 0
  fi

  echo "$(date): Processing ${name}..."

  ## Count read1 mapped fragments for RPM denominator
  local lib_size
  lib_size=$(samtools view -c -F 4 -f 64 "${bam}")
  echo "  ${name} library size (read1 mapped): ${lib_size}"
  echo -e "${name}\t${lib_size}" >> "${LIBSIZE_FILE}"

  ## RPM scalar = 1,000,000 / lib_size
  local rpm
  rpm=$(awk "BEGIN {printf \"%.10f\", 1000000.0 / ${lib_size}}")
  echo "  ${name} RPM scalar: ${rpm}"

  ## BAM -> BED (read1 only, mapped)
  ## -F 4: mapped; -f 64: read1
  echo "  BAM -> BED (read1 only)..."
  samtools view -F 4 -f 64 -b "${bam}" | \
    bedtools bamtobed -i stdin | \
    sort -k1,1 -k2,2n > "${bed_file}"

  echo "  Reads in BED: $(wc -l < ${bed_file})"

  ## Generate strand-specific bedGraphs — RPM normalized
  ## No extension: reads used exactly as aligned by HISAT2
  echo "  Generating RPM-normalized bedGraphs..."

  ## Plus strand bedGraph
  awk '$6 == "+"' "${bed_file}" | \
  bedtools genomecov \
    -i stdin \
    -g "${CHROM_SIZES}" \
    -bg | \
  awk -v rpm="${rpm}" '{printf "%s\t%s\t%s\t%.10f\n", $1,$2,$3,$4*rpm}' | \
  sort -k1,1 -k2,2n > "${plus_bgr}"

  ## Minus strand bedGraph
  awk '$6 == "-"' "${bed_file}" | \
  bedtools genomecov \
    -i stdin \
    -g "${CHROM_SIZES}" \
    -bg | \
  awk -v rpm="${rpm}" '{printf "%s\t%s\t%s\t%.10f\n", $1,$2,$3,$4*rpm}' | \
  sort -k1,1 -k2,2n > "${minus_bgr}"

  echo "  ${name} complete."

  ## Cleanup intermediate BED file
  rm -f "${bed_file}"
}

## Process NoEU first
process_sample "${NOEU_BAM}" "${NOEU_NAME}"

## Process all EU-seq samples
for i in "${!ALL_BAMS_WITH_NOEU[@]}"; do
  ## Skip NoEU — already processed
  if [[ "${ALL_NAMES_WITH_NOEU[$i]}" == "${NOEU_NAME}" ]]; then continue; fi
  process_sample "${ALL_BAMS_WITH_NOEU[$i]}" "${ALL_NAMES_WITH_NOEU[$i]}"
done

echo ""
echo "Library sizes:"
cat "${LIBSIZE_FILE}"

##--------------------------------------------------------------## 6. Subtract NoEU RPM-normalized bedGraph
## Both bedGraphs are already RPM-normalized so they are on
## the same scale before subtraction — no additional scaling needed.
## Subtract per strand, floor at zero.

subtract_noeu() {
  local name="$1"
  local strand="$2"   ## "plus" or "minus"

  local eu_bgr="${BGR_DIR}/${name}_${strand}.norm.bgr"
  local noeu_bgr="${BGR_DIR}/${NOEU_NAME}_${strand}.norm.bgr"
  local out_bgr="${SUB_DIR}/${name}_${strand}.sub.bgr"

  if [[ -s "${out_bgr}" ]]; then
    echo "$(date): Subtracted bedGraph already exists for ${name} (${strand}) — skipping."
    return 0
  fi

  echo "$(date): Subtracting NoEU from ${name} (${strand})..."

  bedtools unionbedg \
    -i "${eu_bgr}" "${noeu_bgr}" \
    -g "${CHROM_SIZES}" \
    -header \
    -names "${name}" "${NOEU_NAME}" | \
  awk 'NR > 1 {
    diff = $4 - $5
    if (diff < 0) diff = 0
    if (diff > 0) printf "%s\t%s\t%s\t%.10f\n", $1,$2,$3,diff
  }' | \
  sort -k1,1 -k2,2n > "${out_bgr}"

  echo "  Written: ${out_bgr}"
}

for name in "${ALL_NAMES[@]}"; do
  subtract_noeu "${name}" "plus"
  subtract_noeu "${name}" "minus"
done

##--------------------------------------------------------------## 7. Count read starts within gene bodies
## Count the number of read start positions (5' of read1)
## falling within each gene's TSS-TES interval, strand-specific.
##
## We use the subtracted bedGraph coverage values rather than
## raw read counts since the subtraction produces fractional RPM
## values. The gene-level score is the sum of RPM × base_pairs
## within the gene body — then divided by gene length (kb) to
## give RPKM.
##
## This matches the Zaret lab approach:
##   score = sum(RPM_coverage × bp) / gene_length_kb
##   which equals RPKM when RPM normalization is done first

compute_rpkm() {
  local name="$1"
  local out_file="${CNT_DIR}/${name}_rpkm.txt"

  if [[ -s "${out_file}" ]]; then
    echo "$(date): RPKM already computed for ${name} — skipping."
    return 0
  fi

  echo "$(date): Computing RPKM for ${name}..."

  local plus_bgr="${SUB_DIR}/${name}_plus.sub.bgr"
  local minus_bgr="${SUB_DIR}/${name}_minus.sub.bgr"

  ## Plus strand genes: sum RPM coverage over gene body, divide by length kb
  bedtools map \
    -a "${GENE_BED_PLUS}" \
    -b "${plus_bgr}" \
    -c 4 -o sum \
    -null 0 | \
  awk '{
    gene_id  = $4
    cov_sum  = $7
    gene_len = $3 - $2
    if (gene_len <= 0) next
    rpkm = cov_sum / (gene_len / 1000.0)
    printf "%s\t%.6f\n", gene_id, rpkm
  }' > "${CNT_DIR}/${name}_plus_rpkm.txt"

  ## Minus strand genes
  bedtools map \
    -a "${GENE_BED_MINUS}" \
    -b "${minus_bgr}" \
    -c 4 -o sum \
    -null 0 | \
  awk '{
    gene_id  = $4
    cov_sum  = $7
    gene_len = $3 - $2
    if (gene_len <= 0) next
    rpkm = cov_sum / (gene_len / 1000.0)
    printf "%s\t%.6f\n", gene_id, rpkm
  }' > "${CNT_DIR}/${name}_minus_rpkm.txt"

  ## Combine and sort
  cat "${CNT_DIR}/${name}_plus_rpkm.txt" \
      "${CNT_DIR}/${name}_minus_rpkm.txt" | \
  sort -k1,1 > "${out_file}"

  local n_expressed
  n_expressed=$(awk '$2 > 0' "${out_file}" | wc -l)
  echo "  Written: ${out_file}"
  echo "  Genes with RPKM > 0: ${n_expressed}"

  ## Cleanup intermediate files
  rm -f "${CNT_DIR}/${name}_plus_rpkm.txt" \
        "${CNT_DIR}/${name}_minus_rpkm.txt"
}

## Also compute RPKM for NoEU (for QC matrix)
## Use raw (unsubtracted) bedgraphs for NoEU
compute_rpkm_noeu() {
  local out_file="${CNT_DIR}/${NOEU_NAME}_rpkm.txt"
  if [[ -s "${out_file}" ]]; then
    echo "$(date): RPKM already computed for ${NOEU_NAME} — skipping."
    return 0
  fi

  echo "$(date): Computing RPKM for ${NOEU_NAME}..."

  bedtools map \
    -a "${GENE_BED_PLUS}" \
    -b "${BGR_DIR}/${NOEU_NAME}_plus.norm.bgr" \
    -c 4 -o sum -null 0 | \
  awk '{
    gene_len = $3 - $2
    if (gene_len <= 0) next
    printf "%s\t%.6f\n", $4, $7 / (gene_len / 1000.0)
  }' > "${CNT_DIR}/${NOEU_NAME}_plus_rpkm.txt"

  bedtools map \
    -a "${GENE_BED_MINUS}" \
    -b "${BGR_DIR}/${NOEU_NAME}_minus.norm.bgr" \
    -c 4 -o sum -null 0 | \
  awk '{
    gene_len = $3 - $2
    if (gene_len <= 0) next
    printf "%s\t%.6f\n", $4, $7 / (gene_len / 1000.0)
  }' > "${CNT_DIR}/${NOEU_NAME}_minus_rpkm.txt"

  cat "${CNT_DIR}/${NOEU_NAME}_plus_rpkm.txt" \
      "${CNT_DIR}/${NOEU_NAME}_minus_rpkm.txt" | \
  sort -k1,1 > "${out_file}"

  rm -f "${CNT_DIR}/${NOEU_NAME}_plus_rpkm.txt" \
        "${CNT_DIR}/${NOEU_NAME}_minus_rpkm.txt"
  echo "  Written: ${out_file}"
}

for name in "${ALL_NAMES[@]}"; do
  compute_rpkm "${name}"
done
compute_rpkm_noeu

##--------------------------------------------------------------## 8. Merge RPKM files into matrices

build_matrix() {
  local out_file="$1"
  shift
  local names=("$@")

  local header="gene_id"
  for name in "${names[@]}"; do
    header="${header}\t${name}"
  done

  ## Gene IDs from first sample
  awk '{print $1}' "${CNT_DIR}/${names[0]}_rpkm.txt" > \
    "${CNT_DIR}/_genes_tmp.txt"

  local col_files=()
  for name in "${names[@]}"; do
    awk '{print $2}' "${CNT_DIR}/${name}_rpkm.txt" > \
      "${CNT_DIR}/_col_${name}_tmp.txt"
    col_files+=("${CNT_DIR}/_col_${name}_tmp.txt")
  done

  paste "${CNT_DIR}/_genes_tmp.txt" "${col_files[@]}" > "${out_file}.tmp"
  echo -e "${header}" | cat - "${out_file}.tmp" > "${out_file}"
  rm -f "${CNT_DIR}/_genes_tmp.txt" "${out_file}.tmp"
  for name in "${names[@]}"; do
    rm -f "${CNT_DIR}/_col_${name}_tmp.txt"
  done

  echo "  Written: ${out_file}"
  echo "  $(wc -l < ${out_file}) genes x ${#names[@]} samples"
}

echo "$(date): Building RPKM matrices..."

## QC matrix — all 7 samples including NoEU
build_matrix \
  "${CNT_DIR}/euseq_rpkm_qc.txt" \
  "${ALL_NAMES[@]}" "${NOEU_NAME}"

## Analysis matrix — 6 EU-seq samples only
build_matrix \
  "${CNT_DIR}/euseq_rpkm_analysis.txt" \
  "${ALL_NAMES[@]}"

echo ""
echo "========================================================"
echo "$(date): Pipeline complete."
echo ""
echo "QC matrix (7 samples):       ${CNT_DIR}/euseq_rpkm_qc.txt"
echo "Analysis matrix (6 samples): ${CNT_DIR}/euseq_rpkm_analysis.txt"
echo "Library sizes:               ${LIBSIZE_FILE}"
echo ""
echo "Expected RPKM range for expressed genes: 0.1 - 1000"
echo "If values are still inflated, check fragment extension"
echo "and bedGraph coverage units."
echo "========================================================"
