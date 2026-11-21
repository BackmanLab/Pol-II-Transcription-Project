#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH -t 12:00:00
#SBATCH --array=0-6
#SBATCH --job-name=dedup_euseq_fq
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/%x_%a_out.txt
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/%x_%a_err.txt
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Removes repeated reads from the EU-seq paired-end FASTQs before
# alignment. In the delivered FASTQs the same sequencer cluster (identical
# header = instrument:run:flowcell:lane:tile:x:y, identical sequence and
# quality) appears more than once, at different places in the file.
# A record is dropped when its read name has already been seen; the first
# occurrence is kept. Only exact repeats of a cluster are removed: reads
# that merely map to the same position (different names) are untouched.
#
# R1 and R2 are processed together so the pairs stay in sync:
#   - a pair is written only if its R1 and R2 headers match; pairs where
#     they do not match (the two files are out of step at that position)
#     are skipped, and the cluster is kept from its other, matched copy
#   - a dropped record must have the same sequence start as the kept
#     record in both mates (else the job stops)
# Input : BASE/euseq/fastq/<sample>_R{1,2}_001.fastq.gz   (not modified)
# Output: BASE/euseq/fastq_dedup/<sample>_R{1,2}_001.fastq.gz
#         BASE/euseq/fastq_dedup/<sample>_dedup_stats.txt
#
# SLURM array index -> sample (same order as align_euseq_hg38.sh):
#   0-2: POLR2A_8hrs_aux BM23, BM24, BM25
#   3-5: WT BM29, BM30, BM31
#   6:   No_EU_ctrl
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/fastq_dedup
#
# Usage:
#   sbatch /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/scripts/dedup_euseq_fastq_hg38.sh
##################################################################

module purge

FASTQ_ARRAY=(
  BM23_S40_R1_001.fastq.gz       ## 0: POLR2A_8hrs_aux Rep1
  BM24_S41_R1_001.fastq.gz       ## 1: POLR2A_8hrs_aux Rep2
  BM25_S42_R1_001.fastq.gz       ## 2: POLR2A_8hrs_aux Rep3
  BM29_S46_R1_001.fastq.gz       ## 3: WT Rep1
  BM30_S47_R1_001.fastq.gz       ## 4: WT Rep2
  BM31_S48_R1_001.fastq.gz       ## 5: WT Rep3
  No_EU_ctrl_S49_R1_001.fastq.gz ## 6: No_EU control
)

R1="${FASTQ_ARRAY[$SLURM_ARRAY_TASK_ID]}"
R2=$(echo "${R1}" | sed 's/R1_001/R2_001/g')
SAMPLE=$(echo "${R1}" | sed 's/_R1_001.fastq.gz//g')

BASE="${BASE:-/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq}"
IN_DIR="${IN_DIR:-${BASE}/euseq/fastq}"
OUT_DIR="${OUT_DIR:-${BASE}/euseq/fastq_dedup}"

die() { echo "*** error: $*" >&2; exit 1; }

[[ -f "${IN_DIR}/${R1}" ]] || die "R1 not found: ${IN_DIR}/${R1}"
[[ -f "${IN_DIR}/${R2}" ]] || die "R2 not found: ${IN_DIR}/${R2}"
mkdir -p "${OUT_DIR}" || die "cannot create ${OUT_DIR}"

OUT1="${OUT_DIR}/${R1}"
OUT2="${OUT_DIR}/${R2}"
STATS="${OUT_DIR}/${SAMPLE}_dedup_stats.txt"
rm -f "${OUT1}.tmp" "${OUT2}.tmp" "${STATS}"

echo "========================================================"
echo "$(date): Removing repeated reads (by read name)"
echo "  Sample : ${SAMPLE}"
echo "  Input  : ${IN_DIR}/${R1}  +  ${R2}"
echo "  Output : ${OUT_DIR}"
echo "========================================================"

## one line per read pair: R1 (4 fields) + R2 (4 fields), tab separated
paste <(zcat "${IN_DIR}/${R1}" | paste - - - -) <(zcat "${IN_DIR}/${R2}" | paste - - - -) | \
awk -F'\t' -v c1="gzip -c > ${OUT1}.tmp" -v c2="gzip -c > ${OUT2}.tmp" -v stats="${STATS}" '
  {
    n++
    split($1, a, " "); split($5, b, " ")
    key = a[1]
    if (key == "") { bad_header++; next }
    ## R1 and R2 out of sync at this position: not written; the cluster is
    ## still registered so a correctly paired copy elsewhere can be kept
    if (key != b[1]) {
      out_of_sync++
      if (!(key in seen)) { seen[key] = "-"; pending++ }
      next
    }
    fp = substr($2, 1, 20) "|" substr($6, 1, 20)
    if (key in seen) {
      if (seen[key] == "-") { seen[key] = fp; pending--; kept++; print $1 "\n" $2 "\n" $3 "\n" $4 | c1; print $5 "\n" $6 "\n" $7 "\n" $8 | c2; next }
      dropped++
      if (seen[key] != fp) content_mismatch++
      next
    }
    seen[key] = fp
    kept++
    print $1 "\n" $2 "\n" $3 "\n" $4 | c1
    print $5 "\n" $6 "\n" $7 "\n" $8 | c2
  }
  END {
    close(c1); close(c2)
    printf "input_pairs\t%d\nkept_pairs\t%d\ndropped_pairs\t%d\nout_of_sync_pairs\t%d\ndistinct_names\t%d\nclusters_without_paired_copy\t%d\ncontent_mismatch\t%d\nempty_header\t%d\n", n, kept, dropped, out_of_sync, length(seen), pending, content_mismatch, bad_header > stats
    close(stats)
    if (content_mismatch > 0 || bad_header > 0) exit 1
  }'
RC=$?
[[ "${RC}" -eq 0 ]] || { rm -f "${OUT1}.tmp" "${OUT2}.tmp"; cat "${STATS}" 2>/dev/null; die "dedup stopped: a repeated name carries different sequence, or a record has an empty header (see ${STATS})"; }

get() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "${STATS}"; }
IN=$(get input_pairs); KEPT=$(get kept_pairs); DROP=$(get dropped_pairs); DIST=$(get distinct_names); SYNC=$(get out_of_sync_pairs); LOST=$(get clusters_without_paired_copy)

echo ""
cat "${STATS}"
echo "  repeated-read rate: $(awk -v d="${DROP}" -v n="${IN}" 'BEGIN {printf "%.2f", 100 * d / n}')% of input pairs"

##--------------------------------------------------------------## Self-checks on the live output
echo ""
echo "========================================================"
echo "Self-checks"
NFAIL=0
ok()  { echo "  PASS  $*"; }
bad() { echo "  FAIL  $*"; NFAIL=$((NFAIL + 1)); }

[[ "${IN}" -eq $((KEPT + DROP + SYNC)) ]]  && ok "input = kept + dropped + out-of-sync (${IN})" || bad "input ${IN} != kept ${KEPT} + dropped ${DROP} + out-of-sync ${SYNC}"
[[ "${KEPT}" -eq $((DIST - LOST)) ]]       && ok "kept pairs = distinct read names - clusters without a paired copy (${KEPT})" || bad "kept ${KEPT} != distinct names ${DIST} - ${LOST}"
[[ "${LOST}" -eq 0 ]]                      && ok "every read name has a correctly paired copy in the output" || echo "  NOTE  ${LOST} read names have no correctly paired copy and are not in the output"
[[ "${KEPT}" -gt 0 ]]               && ok "output is not empty" || bad "no pairs kept"

gzip -t "${OUT1}.tmp" && gzip -t "${OUT2}.tmp" && ok "output gzip files are intact" || bad "gzip integrity test failed"
L1=$(zcat "${OUT1}.tmp" | wc -l); L2=$(zcat "${OUT2}.tmp" | wc -l)
[[ "${L1}" -eq $((KEPT * 4)) && "${L2}" -eq $((KEPT * 4)) ]] && ok "R1 and R2 both have ${KEPT} reads" || bad "read counts differ: R1 $((L1 / 4)), R2 $((L2 / 4)), expected ${KEPT}"

if [[ "${NFAIL}" -gt 0 ]]; then
  rm -f "${OUT1}.tmp" "${OUT2}.tmp"
  die "${NFAIL} self-check(s) failed; no output written"
fi

mv "${OUT1}.tmp" "${OUT1}" && mv "${OUT2}.tmp" "${OUT2}" || die "could not move output into place"

echo "========================================================"
echo "$(date): Done for ${SAMPLE}"
echo "  ${OUT1}"
echo "  ${OUT2}"
echo "========================================================"
