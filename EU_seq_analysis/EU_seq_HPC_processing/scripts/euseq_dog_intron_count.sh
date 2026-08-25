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
# Counts EU-seq signal over two interval types:
#
#   1. INTRONIC intervals — gene body minus exon intervals
#      per gene, collapsed to total intronic RPKM
#      Captures co-transcriptional splicing signal and
#      Pol II elongation through introns
#
#   2. DOWNSTREAM OF GENE (DoG) windows — 5kb downstream
#      of each gene's TES, filtered to exclude windows
#      overlapping any annotated gene on the same strand
#      Captures transcriptional readthrough past poly-A site
#
# Both use the NoEU-subtracted RPM-normalized bedgraphs
# from euseq_pipeline.sh as input.
#
# RPKM formula (same as main pipeline):
#   RPKM = sum(RPM_coverage × bp) / interval_length_kb
#
# Output: RPKM matrices ready for R analysis
#
# Usage:
#   sbatch -A p32171 -p genhimem -t 04:00:00 -n 1 \
#     --cpus-per-task=8 --mem=32G \
#     --output=.../logs/euseq_intron_dog_%j.out \
#     --error=.../logs/euseq_intron_dog_%j.err \
#     euseq_intron_dog.sh
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
OUT_DIR="${PIPELINE_DIR}/intron_dog"
GTF="/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
CHROM_SIZES="/projects/b1042/BackmanLab/Lucas/eric_splicing/chroms/hg38.chrom.sizes"
LIBSIZE_FILE="${PIPELINE_DIR}/library_sizes.txt"

INTRON_DIR="${OUT_DIR}/intron_counts"
DOG_DIR="${OUT_DIR}/dog_counts"
BED_DIR="${OUT_DIR}/beds"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${INTRON_DIR}" "${DOG_DIR}" "${BED_DIR}" "${LOG_DIR}"

DOG_WINDOW=5000   ## 5kb downstream of TES

##--------------------------------------------------------------## 2. Sample definitions

POLR2A_NAMES=("POLR2A_1" "POLR2A_2" "POLR2A_3")
WT_NAMES=("WT_1" "WT_2" "WT_3")
ALL_NAMES=("${POLR2A_NAMES[@]}" "${WT_NAMES[@]}")

##--------------------------------------------------------------## 3. Build interval BED files from GTF

## Gene body BED (reuse from main pipeline if exists)
GENE_BED="${PIPELINE_DIR}/hg38_genes.bed"
if [[ ! -s "${GENE_BED}" ]]; then
  echo "ERROR: Gene BED not found at ${GENE_BED}" >&2
  echo "       Run euseq_pipeline.sh first to generate gene BED" >&2
  exit 1
fi
echo "$(date): Using gene BED: ${GENE_BED}"
echo "  $(wc -l < ${GENE_BED}) genes"

##--------------------------------------------------------------## 3a. Build exon BED from GTF

EXON_BED="${BED_DIR}/hg38_exons_merged.bed"

if [[ ! -s "${EXON_BED}" ]]; then
  echo "$(date): Building exon BED from GTF..."

  ## Extract exon features — one row per exon
  ## Merge overlapping exons per gene to avoid double-counting
  ## alternative isoforms
  grep -v "^#" "${GTF}" | \
  awk '$3 == "exon" {
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
  sort -k1,1 -k2,2n > "${BED_DIR}/hg38_exons_raw.bed"

  echo "  Raw exon intervals: $(wc -l < ${BED_DIR}/hg38_exons_raw.bed)"

  ## Merge overlapping exons per gene using bedtools merge
  ## Sort by gene_id then coordinate, merge within gene
  ## Use groupBy to preserve gene_id after merge
  bedtools merge \
    -i "${BED_DIR}/hg38_exons_raw.bed" \
    -s \
    -c 4,6 \
    -o distinct,distinct | \
  awk '{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5}' | \
  sort -k1,1 -k2,2n > "${EXON_BED}"

  echo "  Merged exon intervals: $(wc -l < ${EXON_BED})"
else
  echo "Exon BED already exists: ${EXON_BED}"
fi

##--------------------------------------------------------------## 3b. Build intronic intervals
## Intron = gene body interval MINUS all exon intervals for that gene
## Use bedtools subtract: for each gene, subtract its exons from
## its gene body interval to get pure intronic regions

INTRON_BED_PLUS="${BED_DIR}/hg38_introns_plus.bed"
INTRON_BED_MINUS="${BED_DIR}/hg38_introns_minus.bed"

if [[ ! -s "${INTRON_BED_PLUS}" ]]; then
  echo "$(date): Building intronic intervals..."

  ## Separate gene and exon BEDs by strand
  awk '$6 == "+"' "${GENE_BED}"  > "${BED_DIR}/genes_plus.bed"
  awk '$6 == "-"' "${GENE_BED}"  > "${BED_DIR}/genes_minus.bed"
  awk '$6 == "+"' "${EXON_BED}"  > "${BED_DIR}/exons_plus.bed"
  awk '$6 == "-"' "${EXON_BED}"  > "${BED_DIR}/exons_minus.bed"

  ## Subtract exons from gene bodies per strand
  ## -s: strand-specific subtraction
  ## Result: intronic intervals with gene_id in col4
  bedtools subtract \
    -a "${BED_DIR}/genes_plus.bed" \
    -b "${BED_DIR}/exons_plus.bed" \
    -s | \
  sort -k1,1 -k2,2n > "${INTRON_BED_PLUS}"

  bedtools subtract \
    -a "${BED_DIR}/genes_minus.bed" \
    -b "${BED_DIR}/exons_minus.bed" \
    -s | \
  sort -k1,1 -k2,2n > "${INTRON_BED_MINUS}"

  echo "  Plus strand intronic intervals:  $(wc -l < ${INTRON_BED_PLUS})"
  echo "  Minus strand intronic intervals: $(wc -l < ${INTRON_BED_MINUS})"
  echo "  Unique genes with introns (plus):  $(awk '{print $4}' ${INTRON_BED_PLUS}  | sort -u | wc -l)"
  echo "  Unique genes with introns (minus): $(awk '{print $4}' ${INTRON_BED_MINUS} | sort -u | wc -l)"
else
  echo "Intronic BEDs already exist."
fi

##--------------------------------------------------------------## 3c. Build DoG windows
## 5kb downstream of TES, filtered to exclude windows that
## overlap any annotated gene on the same strand

DOG_BED_PLUS="${BED_DIR}/dog_windows_plus.bed"
DOG_BED_MINUS="${BED_DIR}/dog_windows_minus.bed"

if [[ ! -s "${DOG_BED_PLUS}" ]]; then
  echo "$(date): Building ${DOG_WINDOW}bp DoG windows..."

  ## For + strand genes: window is TES to TES+DOG_WINDOW
  ##   TES = col3 (BED end)
  ##   Window: chr, TES, TES+window, gene_id, score, strand
  awk -v win="${DOG_WINDOW}" '$6 == "+" {
    print $1"\t"$3"\t"($3+win)"\t"$4"\t0\t"$6
  }' "${GENE_BED}" | \
  ## Clip to chromosome boundaries
  bedtools slop -i stdin -g "${CHROM_SIZES}" -l 0 -r 0 | \
  sort -k1,1 -k2,2n > "${BED_DIR}/dog_raw_plus.bed"

  ## For - strand genes: window is TSS-DOG_WINDOW to TSS
  ##   TSS of - strand gene = col2 (BED start)
  ##   Window: chr, TSS-window, TSS, gene_id, score, strand
  awk -v win="${DOG_WINDOW}" '$6 == "-" {
    start = $2 - win
    if (start < 0) start = 0
    print $1"\t"start"\t"$2"\t"$4"\t0\t"$6
  }' "${GENE_BED}" | \
  sort -k1,1 -k2,2n > "${BED_DIR}/dog_raw_minus.bed"

  echo "  Raw DoG windows (plus):  $(wc -l < ${BED_DIR}/dog_raw_plus.bed)"
  echo "  Raw DoG windows (minus): $(wc -l < ${BED_DIR}/dog_raw_minus.bed)"

  ## Filter: remove windows that overlap any annotated gene
  ## on the same strand — these would contaminate DoG signal
  ## with legitimate transcription from neighboring genes
  ## -v: report only windows with NO overlap
  ## -s: strand-specific overlap check
  bedtools intersect \
    -a "${BED_DIR}/dog_raw_plus.bed" \
    -b "${GENE_BED}" \
    -s -v | \
  sort -k1,1 -k2,2n > "${DOG_BED_PLUS}"

  bedtools intersect \
    -a "${BED_DIR}/dog_raw_minus.bed" \
    -b "${GENE_BED}" \
    -s -v | \
  sort -k1,1 -k2,2n > "${DOG_BED_MINUS}"

  echo "  Clean DoG windows (plus):  $(wc -l < ${DOG_BED_PLUS})"
  echo "  Clean DoG windows (minus): $(wc -l < ${DOG_BED_MINUS})"
  echo "  Genes with clean DoG windows: $(cat ${DOG_BED_PLUS} ${DOG_BED_MINUS} | awk '{print $4}' | sort -u | wc -l)"
else
  echo "DoG window BEDs already exist."
fi

##--------------------------------------------------------------## 4. Compute RPKM over intervals
## For each sample, map subtracted RPM-normalized bedGraph
## coverage over intronic and DoG intervals using bedtools map
##
## RPKM = sum(RPM_coverage × bp) / interval_length_kb
##
## For introns: aggregate all intronic intervals per gene
## using bedtools groupby after mapping

compute_interval_rpkm() {
  local name="$1"
  local interval_plus="$2"
  local interval_minus="$3"
  local out_file="$4"
  local interval_type="$5"   ## "intron" or "dog"

  if [[ -s "${out_file}" ]]; then
    echo "$(date): ${interval_type} RPKM already computed for ${name} — skipping."
    return 0
  fi

  echo "$(date): Computing ${interval_type} RPKM for ${name}..."

  local plus_bgr="${SUB_DIR}/${name}_plus.sub.bgr"
  local minus_bgr="${SUB_DIR}/${name}_minus.sub.bgr"

  if [[ ! -s "${plus_bgr}" || ! -s "${minus_bgr}" ]]; then
    echo "  ERROR: subtracted bedgraphs not found for ${name}" >&2
    echo "         Expected: ${plus_bgr}" >&2
    exit 1
  fi

  local tmp_plus="${out_file%.txt}_plus_tmp.txt"
  local tmp_minus="${out_file%.txt}_minus_tmp.txt"

  ## Map RPM coverage to intervals — sum coverage over each interval
  bedtools map \
    -a "${interval_plus}" \
    -b "${plus_bgr}" \
    -c 4 -o sum \
    -null 0 | \
  awk '{
    gene_id  = $4
    cov_sum  = $7
    int_len  = $3 - $2
    if (int_len <= 0) next
    ## Store: gene_id, coverage_sum, interval_length
    print gene_id"\t"cov_sum"\t"int_len
  }' > "${tmp_plus}"

  bedtools map \
    -a "${interval_minus}" \
    -b "${minus_bgr}" \
    -c 4 -o sum \
    -null 0 | \
  awk '{
    gene_id  = $4
    cov_sum  = $7
    int_len  = $3 - $2
    if (int_len <= 0) next
    print gene_id"\t"cov_sum"\t"int_len
  }' > "${tmp_minus}"

  ## For introns: multiple intervals per gene — sum coverage
  ## and sum lengths across all introns, then compute RPKM
  ## For DoG: one interval per gene — compute RPKM directly
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

  local n_signal
  n_signal=$(awk '$2 > 0' "${out_file}" | wc -l)
  echo "  Written: ${out_file}"
  echo "  Genes with ${interval_type} RPKM > 0: ${n_signal}"
}

## Compute intronic RPKM for all samples
echo ""
echo "$(date): Computing intronic RPKM..."
for name in "${ALL_NAMES[@]}"; do
  compute_interval_rpkm \
    "${name}" \
    "${INTRON_BED_PLUS}" \
    "${INTRON_BED_MINUS}" \
    "${INTRON_DIR}/${name}_intron_rpkm.txt" \
    "intron"
done

## Compute DoG RPKM for all samples
echo ""
echo "$(date): Computing DoG RPKM..."
for name in "${ALL_NAMES[@]}"; do
  compute_interval_rpkm \
    "${name}" \
    "${DOG_BED_PLUS}" \
    "${DOG_BED_MINUS}" \
    "${DOG_DIR}/${name}_dog_rpkm.txt" \
    "dog"
done

##--------------------------------------------------------------## 5. Merge into matrices

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

  awk '{print $1}' "${cnt_dir}/${names[0]}_${suffix}.txt" > \
    "${cnt_dir}/_genes_tmp.txt"

  local col_files=()
  for name in "${names[@]}"; do
    awk '{print $2}' "${cnt_dir}/${name}_${suffix}.txt" > \
      "${cnt_dir}/_col_${name}_tmp.txt"
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

build_matrix \
  "${INTRON_DIR}/euseq_intron_rpkm.txt" \
  "${INTRON_DIR}" \
  "intron_rpkm" \
  "${ALL_NAMES[@]}"

build_matrix \
  "${DOG_DIR}/euseq_dog_rpkm.txt" \
  "${DOG_DIR}" \
  "dog_rpkm" \
  "${ALL_NAMES[@]}"

echo ""
echo "========================================================"
echo "$(date): Pipeline complete."
echo ""
echo "Intronic RPKM matrix: ${INTRON_DIR}/euseq_intron_rpkm.txt"
echo "DoG RPKM matrix:      ${DOG_DIR}/euseq_dog_rpkm.txt"
echo ""
echo "Next steps:"
echo "  1. Load matrices in R"
echo "  2. Apply total RNA-seq scale factors"
echo "  3. Compute intron/exon RPKM ratio per gene"
echo "  4. Compare DoG RPKM POLR2A vs WT"
echo "========================================================"
