#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=euseq_subtract
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/090124_euseq/logs/euseq_term_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/090124_euseq/logs/euseq_term_%j.err

##################################################################
# Description:
# Computes read-2 3' end termination position pileups within
# introns and exons, from filtered BAMs (flag 128 = read2,
# -F 4/-F 256 exclude unmapped/secondary). Output consumed by
# EU_seq_termination.R for plotting.
#
# NOTE: header previously copy-pasted from euseq_dog_intron_count.sh
# and described the wrong script (intron/DoG RPKM counting) -- fixed.
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

BAM_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/BAM"
BED_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline/intron_dog/beds"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline/termination"

mkdir -p "${OUT_DIR}/intron" "${OUT_DIR}/exon" "${OUT_DIR}/logs"

INTRON_PLUS="${BED_DIR}/hg38_introns_plus.bed"
INTRON_MINUS="${BED_DIR}/hg38_introns_minus.bed"
EXON_PLUS="${BED_DIR}/exons_plus.bed"
EXON_MINUS="${BED_DIR}/exons_minus.bed"

THREADS=8

## Verify BED files exist
for f in "${INTRON_PLUS}" "${INTRON_MINUS}" "${EXON_PLUS}" "${EXON_MINUS}"; do
  if [[ ! -s "${f}" ]]; then
    echo "ERROR: BED file not found: ${f}" >&2
    exit 1
  fi
  echo "Found: $(basename ${f}) — $(wc -l < ${f}) intervals"
done

##--------------------------------------------------------------## 2. Sample definitions

declare -A SAMPLES
SAMPLES["POLR2A_1"]="BM23_S40_filt.sort.bam"
SAMPLES["POLR2A_2"]="BM24_S41_filt.sort.bam"
SAMPLES["POLR2A_3"]="BM25_S42_filt.sort.bam"
SAMPLES["WT_1"]="BM29_S46_filt.sort.bam"
SAMPLES["WT_2"]="BM30_S47_filt.sort.bam"
SAMPLES["WT_3"]="BM31_S48_filt.sort.bam"

##--------------------------------------------------------------## 3. Process each sample

process_sample() {
  local name="$1"
  local bam="${BAM_DIR}/${SAMPLES[$name]}"
  local tmp_dir="${OUT_DIR}/tmp_${name}"
  mkdir -p "${tmp_dir}"

  echo "$(date): Processing ${name}..."
  echo "  BAM: ${bam}"

  if [[ ! -s "${bam}" ]]; then
    echo "ERROR: BAM not found: ${bam}" >&2
    exit 1
  fi

  if [[ ! -s "${bam}.bai" ]]; then
    echo "  Indexing BAM..."
    samtools index -@ "${THREADS}" "${bam}"
  fi

  ## --- Step 1: Extract read 2 → BED ---
  ## flag 128 = read2, -F 4 = exclude unmapped, -F 256 = exclude secondary
  echo "  Extracting read 2..."
  samtools view \
    -@ "${THREADS}" \
    -f 128 \
    -F 4 \
    -F 256 \
    -b "${bam}" | \
  bedtools bamtobed -i stdin | \
  sort -k1,1 -k2,2n > "${tmp_dir}/${name}_read2.bed"

  echo "  Read 2 count: $(wc -l < ${tmp_dir}/${name}_read2.bed)"

  ## --- Step 2: Compute 3' end per read ---
  ## + strand: 3' end = end-1 to end (rightmost base, 0-based)
  ## - strand: 3' end = start to start+1 (leftmost base, 0-based)
  echo "  Computing 3' ends..."
  awk 'BEGIN{OFS="\t"} {
    if ($6 == "+") {
      print $1, $3-1, $3, $4, $5, $6
    } else {
      print $1, $2, $2+1, $4, $5, $6
    }
  }' "${tmp_dir}/${name}_read2.bed" | \
  sort -k1,1 -k2,2n > "${tmp_dir}/${name}_3prime.bed"

  awk '$6 == "+"' "${tmp_dir}/${name}_3prime.bed" > \
    "${tmp_dir}/${name}_3prime_plus.bed"
  awk '$6 == "-"' "${tmp_dir}/${name}_3prime.bed" > \
    "${tmp_dir}/${name}_3prime_minus.bed"

  echo "  3' ends plus:  $(wc -l < ${tmp_dir}/${name}_3prime_plus.bed)"
  echo "  3' ends minus: $(wc -l < ${tmp_dir}/${name}_3prime_minus.bed)"

  ## --- Step 3: Intron intersection ---
  ## Absolute position from intron start
  ## + strand: pos = 3p_start - intron_start
  ## - strand: pos = intron_end - 3p_end
  echo "  Intersecting with introns..."

  bedtools intersect \
    -a "${tmp_dir}/${name}_3prime_plus.bed" \
    -b "${INTRON_PLUS}" \
    -wa -wb -s | \
  awk 'BEGIN{OFS="\t"} {
    pos        = $2 - $8
    intron_len = $9 - $8
    gene_id    = $10
    if (pos >= 0 && pos < intron_len)
      print gene_id, intron_len, pos, "+"
  }' > "${tmp_dir}/${name}_intron_pos_plus.txt"

  bedtools intersect \
    -a "${tmp_dir}/${name}_3prime_minus.bed" \
    -b "${INTRON_MINUS}" \
    -wa -wb -s | \
  awk 'BEGIN{OFS="\t"} {
    pos        = $9 - $3
    intron_len = $9 - $8
    gene_id    = $10
    if (pos >= 0 && pos < intron_len)
      print gene_id, intron_len, pos, "-"
  }' > "${tmp_dir}/${name}_intron_pos_minus.txt"

  cat "${tmp_dir}/${name}_intron_pos_plus.txt" \
      "${tmp_dir}/${name}_intron_pos_minus.txt" \
  > "${OUT_DIR}/intron/${name}_intron_positions.txt"

  echo "  Intron positions: $(wc -l < ${OUT_DIR}/intron/${name}_intron_positions.txt)"

  ## --- Step 4: Exon intersection ---
  ## Relative position within exon: pos / exon_len (0 to 1)
  echo "  Intersecting with exons..."

  bedtools intersect \
    -a "${tmp_dir}/${name}_3prime_plus.bed" \
    -b "${EXON_PLUS}" \
    -wa -wb -s | \
  awk 'BEGIN{OFS="\t"} {
    pos      = $2 - $8
    exon_len = $9 - $8
    gene_id  = $10
    if (pos >= 0 && pos < exon_len && exon_len > 0) {
      rel_pos = pos / exon_len
      print gene_id, exon_len, rel_pos, "+"
    }
  }' > "${tmp_dir}/${name}_exon_pos_plus.txt"

  bedtools intersect \
    -a "${tmp_dir}/${name}_3prime_minus.bed" \
    -b "${EXON_MINUS}" \
    -wa -wb -s | \
  awk 'BEGIN{OFS="\t"} {
    pos      = $9 - $3
    exon_len = $9 - $8
    gene_id  = $10
    if (pos >= 0 && pos < exon_len && exon_len > 0) {
      rel_pos = pos / exon_len
      print gene_id, exon_len, rel_pos, "-"
    }
  }' > "${tmp_dir}/${name}_exon_pos_minus.txt"

  cat "${tmp_dir}/${name}_exon_pos_plus.txt" \
      "${tmp_dir}/${name}_exon_pos_minus.txt" \
  > "${OUT_DIR}/exon/${name}_exon_positions.txt"

  echo "  Exon positions: $(wc -l < ${OUT_DIR}/exon/${name}_exon_positions.txt)"

  rm -rf "${tmp_dir}"
  echo "$(date): ${name} done."
  echo ""
}

for name in "${!SAMPLES[@]}"; do
  process_sample "${name}"
done

echo "========================================================"
echo "$(date): All samples complete."
echo ""
echo "Output:"
ls -lh "${OUT_DIR}/intron/"
ls -lh "${OUT_DIR}/exon/"
echo ""
echo "Next: run euseq_termination_pileup.R"
echo "========================================================"
