#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --job-name=euseq_subtract_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/euseq_subtract_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/euseq_subtract_%j.err
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# EU-seq NoEU background subtraction on hg38, per
# nucleotide, in bedGraph space. Based on euseq_subtract_count.sh.
# Mito (chrM) scaling is NOT applied here: this step is unchanged
# (RPM per library, NoEU subtracted per nucleotide, floored at 0).
# The chrM factor (compute_scale_factors_mito_hg38.sh) multiplies the
# subtracted RPKM afterwards (bigWigs; R scripts).
#   1. BAM -> BED (read1 only, primary alignments; spliced reads are
#      split at N so only aligned blocks carry coverage; the library
#      is sense-stranded, read1 strand = gene strand)
#   2. Strand-specific coverage, RPM-normalized (1e6 / read1 library)
#   3. NoEU subtracted per nucleotide, per strand, floored at 0
#   4. library_sizes.txt: sample, read1 library size, mean aligned
#      bases per read (from the RPM bedGraph itself); the latter
#      converts RPM coverage to RPKM in euseq_rpkm_hg38.sh and in the
#      EU-seq bigWigs.
#
# Outputs (OUT_DIR = .../euseq/euseq_pipeline):
#   bedgraphs/{sample}_{plus,minus}.norm.bgr   RPM coverage
#   subtracted/{sample}_{plus,minus}.sub.bgr   NoEU-subtracted RPM
#   library_sizes.txt
#
# Conditions: POLR2A_8hrs_aux (R1-3) vs WT (R1-3); NoEU control.
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs
# Usage:
#   sbatch euseq_subtract_hg38.sh
# Next: sbatch euseq_rpkm_hg38.sh
##################################################################

set -o pipefail

module purge
module load samtools/1.6
module load bedtools/2.29.2

##--------------------------------------------------------------## 1. Paths

BASE="${BASE:-/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq}"
BAM_DIR="${BAM_DIR:-${BASE}/euseq/BAM}"
OUT_DIR="${OUT_DIR:-${BASE}/euseq/euseq_pipeline}"
GENOME_FA="${GENOME_FA:-${BASE}/reference/hg38.fa}"
CHROM_SIZES="${CHROM_SIZES:-${BASE}/reference/hg38.chrom.sizes}"

BGR_DIR="${OUT_DIR}/bedgraphs"
SUB_DIR="${OUT_DIR}/subtracted"
MIN_READS="${MIN_READS:-1000000}"   ## self-check floor; the toy test lowers it
BED_DIR="${OUT_DIR}/beds"
LIBSIZE_FILE="${OUT_DIR}/library_sizes.txt"

mkdir -p "${BED_DIR}" "${BGR_DIR}" "${SUB_DIR}"

die() { echo "*** error: $*" >&2; exit 1; }

## chrom.sizes from the hg38 FASTA index if it does not exist yet (same contigs as the BAM headers)
if [[ ! -s "${CHROM_SIZES}" ]]; then
  [[ -s "${GENOME_FA}" ]] || die "no chrom.sizes (${CHROM_SIZES}) and no FASTA to build it from (${GENOME_FA}); set GENOME_FA"
  samtools faidx "${GENOME_FA}" || die "samtools faidx failed on ${GENOME_FA}"
  cut -f1,2 "${GENOME_FA}.fai" > "${CHROM_SIZES}" || die "could not write ${CHROM_SIZES}"
  echo "Built ${CHROM_SIZES} from ${GENOME_FA}.fai"
fi
[[ -s "${CHROM_SIZES}" ]] || die "chrom.sizes not found: ${CHROM_SIZES}"

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

for bam in "${ALL_BAMS[@]}" "${NOEU_BAM}"; do
  [[ -f "${bam}" ]] || die "BAM not found: ${bam}"
done

##--------------------------------------------------------------## 3. BAM -> RPM bedGraph, library size, mean aligned length
## read1 (flag 64), mapped reads with MAPQ > 0 (-q 1), primary alignments only
## (-F 2308: unmapped, secondary, supplementary). The same read set is the
## library size (the "per million mapped reads" of RPKM) and the coverage. bamtobed -split reports one BED line per aligned
## block, so the intron of a spliced read carries no coverage.

LIB_TMP="${LIBSIZE_FILE}.tmp"
> "${LIB_TMP}"

process_sample() {
  local bam="$1"
  local name="$2"

  local bed_file="${BED_DIR}/${name}.bed"
  local plus_bgr="${BGR_DIR}/${name}_plus.norm.bgr"
  local minus_bgr="${BGR_DIR}/${name}_minus.norm.bgr"
  local done_marker="${BGR_DIR}/${name}.split_v1.done"

  local lib_size
  lib_size=$(samtools view -c -q 1 -F 2308 -f 64 "${bam}") || die "samtools failed on ${bam}"
  [[ "${lib_size}" =~ ^[0-9]+$ && "${lib_size}" -gt 0 ]] || die "bad library size '${lib_size}' for ${name}"
  echo "$(date): ${name} library size (read1 mapped): ${lib_size}"

  local rpm
  rpm=$(awk -v n="${lib_size}" 'BEGIN {printf "%.10f", 1000000.0 / n}')
  echo "  RPM scalar: ${rpm}"

  if [[ -s "${plus_bgr}" && -s "${minus_bgr}" && -f "${done_marker}" ]]; then
    echo "  bedGraphs already built by this version -- reusing."
  else
    rm -f "${done_marker}"
    echo "  BAM -> BED (read1, split at N)..."
    samtools view -q 1 -F 2308 -f 64 -b "${bam}" | \
      bedtools bamtobed -split -i stdin | \
      LC_ALL=C sort -k1,1 -k2,2n > "${bed_file}" || die "BAM -> BED failed for ${name}"
    [[ -s "${bed_file}" ]] || die "no reads in BED for ${name}"
    echo "  Aligned blocks in BED: $(wc -l < "${bed_file}")"

    local strand out_bgr
    for strand in plus minus; do
      if [[ "${strand}" == "plus" ]]; then out_bgr="${plus_bgr}"; sym="+"; else out_bgr="${minus_bgr}"; sym="-"; fi
      awk -v s="${sym}" '$6 == s' "${bed_file}" | \
        bedtools genomecov -i stdin -g "${CHROM_SIZES}" -bg | \
        awk -v rpm="${rpm}" '{printf "%s\t%s\t%s\t%.10f\n", $1, $2, $3, $4 * rpm}' | \
        LC_ALL=C sort -k1,1 -k2,2n > "${out_bgr}" || die "genomecov failed for ${name} ${strand}"
    done
    rm -f "${bed_file}"
    touch "${done_marker}"
  fi

  ## Mean aligned bases per read = (integral of RPM coverage / rpm) / reads
  local mean_len
  mean_len=$(cat "${plus_bgr}" "${minus_bgr}" | \
    awk -v rpm="${rpm}" -v n="${lib_size}" '{bp += $4 * ($3 - $2)} END {printf "%.4f", bp / rpm / n}')
  echo "  Mean aligned bases per read: ${mean_len}"
  [[ -n "${mean_len}" ]] || die "mean read length failed for ${name}"

  printf "%s\t%s\t%s\n" "${name}" "${lib_size}" "${mean_len}" >> "${LIB_TMP}"
}

process_sample "${NOEU_BAM}" "${NOEU_NAME}"
for i in "${!ALL_BAMS[@]}"; do
  process_sample "${ALL_BAMS[$i]}" "${ALL_NAMES[$i]}"
done

mv "${LIB_TMP}" "${LIBSIZE_FILE}"
echo ""
echo "Library sizes (sample, read1 library, mean aligned bases per read):"
cat "${LIBSIZE_FILE}"

##--------------------------------------------------------------## 4. Subtract NoEU per nucleotide, floor at 0
## Both bedGraphs are RPM, so no additional scaling. Per strand,
## EU - NoEU, negative values set to 0; zeros are not written. Differences
## <= 1e-6 RPM are treated as 0: the bedGraphs are printed to 10 decimals,
## so an exact cancellation leaves a ~1e-10 rounding residue (one read at
## a 10M-read library is 0.1 RPM, so 1e-6 RPM is far below any real signal).

subtract_noeu() {
  local name="$1"
  local strand="$2"

  local eu_bgr="${BGR_DIR}/${name}_${strand}.norm.bgr"
  local noeu_bgr="${BGR_DIR}/${NOEU_NAME}_${strand}.norm.bgr"
  local out_bgr="${SUB_DIR}/${name}_${strand}.sub.bgr"

  echo "$(date): Subtracting NoEU from ${name} (${strand})..."

  bedtools unionbedg -i "${eu_bgr}" "${noeu_bgr}" -filler 0 | \
  awk 'BEGIN {OFS = "\t"}
    NF != 5 { print "*** error: unionbedg line has " NF " fields, expected 5" > "/dev/stderr"; exit 1 }
    { diff = $4 - $5
      if (diff > 1e-6) printf "%s\t%s\t%s\t%.10f\n", $1, $2, $3, diff }' | \
  LC_ALL=C sort -k1,1 -k2,2n > "${out_bgr}" || die "subtraction failed for ${name} ${strand}"

  echo "  Written: ${out_bgr} ($(wc -l < "${out_bgr}") intervals)"
}

for name in "${ALL_NAMES[@]}"; do
  subtract_noeu "${name}" "plus"
  subtract_noeu "${name}" "minus"
done

echo ""
echo "$(date): Subtraction complete."
echo "Subtracted bedGraphs: ${SUB_DIR}"
echo "Library sizes:        ${LIBSIZE_FILE}"

##--------------------------------------------------------------## 5. Self-checks on the live outputs
echo "========================================================"
echo "Self-checks"
NFAIL=0
chk() { if [[ "$1" == ok ]]; then echo "  CHECK PASS  $2"; else echo "  CHECK FAIL  $2"; NFAIL=$((NFAIL + 1)); fi; }

## bgr_stats FILE -> "problems integral": problems = unknown chrom, interval past chrom end,
## empty interval, value <= 0, overlap / unsorted / non-contiguous chromosome; integral = sum(value x bp)
bgr_stats() {
  awk -F'\t' 'FNR == NR { size[$1] = $2; next }
    { if (!($1 in size)) bad++
      else if ($3 > size[$1]) bad++
      if ($2 >= $3 || $4 <= 0) bad++
      if ($1 != pc) { if ($1 in seen) bad++; seen[$1] = 1 } else if ($2 < pe) bad++
      pc = $1; pe = $3; sig += $4 * ($3 - $2) }
    END { printf "%d %.6f", bad + 0, sig }' "${CHROM_SIZES}" "$1"
}

chk "$([[ $(wc -l < "${LIBSIZE_FILE}") -eq $(( ${#ALL_NAMES[@]} + 1 )) ]] && echo ok)" "library_sizes.txt has $(( ${#ALL_NAMES[@]} + 1 )) samples"
chk "$(awk -F'\t' -v m="${MIN_READS}" '$2 < m || $3 < 20 || $3 > 1000 {bad++} END {print bad ? "" : "ok"}' "${LIBSIZE_FILE}")" "read1 library >= ${MIN_READS} and mean aligned bases per read in 20-1000"
for name in "${ALL_NAMES[@]}"; do
  for strand in plus minus; do
    read -r prob sub_sig <<< "$(bgr_stats "${SUB_DIR}/${name}_${strand}.sub.bgr")"
    read -r prob_n norm_sig <<< "$(bgr_stats "${BGR_DIR}/${name}_${strand}.norm.bgr")"
    chk "$([[ "${prob}" == 0 && "${prob_n}" == 0 && $(wc -l < "${SUB_DIR}/${name}_${strand}.sub.bgr") -gt 0 ]] && echo ok)" "${name} ${strand}: bedGraphs non-empty, sorted, non-overlapping, in chrom bounds, values > 0"
    chk "$(awk -v a="${sub_sig}" -v b="${norm_sig}" 'BEGIN {print (a <= b * (1 + 1e-6) && a > 0) ? "ok" : ""}')" "${name} ${strand}: subtracted signal <= EU signal"
  done
done
echo "========================================================"
[[ ${NFAIL} -eq 0 ]] || { echo "*** ${NFAIL} self-check(s) failed" >&2; exit 1; }
