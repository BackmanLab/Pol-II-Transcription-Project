#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH -t 12:00:00
#SBATCH --job-name=sf_mito_total_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/logs/scale_factors_mito_total_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/logs/scale_factors_mito_total_%j.err
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Mitochondrial (chrM) scale factors for the Pol II TOTAL RNA-seq
# (hg38, single-end), the total RNA-seq counterpart of
# compute_scale_factors_mito_hg38.sh. Same definition:
#
#   chrM_rpm(sample) = chrM reads / mapped library * 1e6
#   mito_sf(condition) = mean chrM_rpm (WT) / mean chrM_rpm (condition)
#                        (WT = 1)
#
# Counted on the INDIVIDUAL total RNA-seq BAMs (POLR2A Rep1-3, WT
# Rep1-3; POLR2A Rep4 dropped) with the read set bamCoverage sees:
# mapped, MAPQ > 0, primary (samtools view -c -q 1 -F 2308).
# generate_total_bigwigs_mito_hg38.sh applies mito_sf to the merged
# BAMs.
#
# Outputs (BASE):
#   scale_factors_mito_total_samples.txt : sample, condition, library,
#                                          chrM_reads, chrM_rpm
#   scale_factors_mito_total.txt         : condition, n_reps,
#                                          mean_chrM_rpm, sum_chrM_reads,
#                                          mito_sf (col 5), sum_reads
#
# Usage:
#   sbatch compute_scale_factors_mito_total_hg38.sh
##################################################################

set -o pipefail

module purge
module load samtools/1.14

BASE="${BASE:-/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq}"
BAM_DIR="${BAM_DIR:-${BASE}/total/BAM}"
MITO="${MITO:-chrM}"
OUT_SAMPLES="${BASE}/scale_factors_mito_total_samples.txt"
OUT_SF="${BASE}/scale_factors_mito_total.txt"
EU_SF_FILE="${EU_SF_FILE:-${BASE}/scale_factors_mito.txt}"

MIN_READS="${MIN_READS:-1000000}"   ## self-check floor; the toy test lowers it
MIN_MITO="${MIN_MITO:-1000}"        ## chrM reads needed for a stable factor
REF_COND="WT"
COND_ORDER=("WT" "POLR2A_8hrs_aux")

die() { echo "*** error: $*" >&2; exit 1; }

## sample | condition | BAM        (POLR2A Rep4 / S8 dropped: 3 vs 3)
SAMPLES=(
  "WT_1|WT|${BAM_DIR}/WT_HCT116_CTRL_Rep1_S12_filt.sort.bam"
  "WT_2|WT|${BAM_DIR}/WT_HCT116_CTRL_Rep2_S13_filt.sort.bam"
  "WT_3|WT|${BAM_DIR}/WT_HCT116_CTRL_Rep3_S14_filt.sort.bam"
  "POLR2A_1|POLR2A_8hrs_aux|${BAM_DIR}/Pol2_8hrs_Aux_Rep1_S5_filt.sort.bam"
  "POLR2A_2|POLR2A_8hrs_aux|${BAM_DIR}/Pol2_8hrs_Aux_Rep2_S6_filt.sort.bam"
  "POLR2A_3|POLR2A_8hrs_aux|${BAM_DIR}/Pol2_8hrs_Aux_Rep3_S7_filt.sort.bam"
)

echo "========================================================"
echo "$(date): Counting ${MITO} reads in the total RNA-seq BAMs"
echo "========================================================"

## count_pair BAM -> "library chrM"   (same filters; library = whole BAM, chrM = region)
count_pair() {
  local bam="$1" lib m
  [[ -f "${bam}" ]]     || die "BAM not found: ${bam}"
  [[ -f "${bam}.bai" ]] || [[ -f "${bam%.bam}.bai" ]] || die "BAM index not found for ${bam} (samtools index)"
  samtools view -H "${bam}" | grep -q -P "^@SQ\tSN:${MITO}\t" || die "contig ${MITO} not in header of ${bam}"
  lib=$(samtools view -c -q 1 -F 2308 "${bam}")         || die "samtools failed on ${bam}"
  m=$(samtools view -c -q 1 -F 2308 "${bam}" "${MITO}") || die "samtools failed on ${bam} ${MITO}"
  [[ "${lib}" =~ ^[0-9]+$ && "${lib}" -gt 0 ]] || die "bad library size '${lib}' for ${bam}"
  [[ "${m}" =~ ^[0-9]+$ ]] || die "bad ${MITO} count '${m}' for ${bam}"
  echo "${lib} ${m}"
}

printf "sample\tcondition\tlibrary\tchrM_reads\tchrM_rpm\n" > "${OUT_SAMPLES}"
for entry in "${SAMPLES[@]}"; do
  IFS='|' read -r sample cond bam <<< "${entry}"
  read -r lib m <<< "$(count_pair "${bam}")"
  [[ -n "${m}" ]] || exit 1
  rpm=$(awk -v m="${m}" -v n="${lib}" 'BEGIN {printf "%.6f", m * 1000000.0 / n}')
  printf "%s\t%s\t%s\t%s\t%s\n" "${sample}" "${cond}" "${lib}" "${m}" "${rpm}" >> "${OUT_SAMPLES}"
  echo "  ${sample} (${cond}): library ${lib} | ${MITO} ${m} | ${rpm} RPM"
done

echo ""
echo "$(date): Computing condition-level factors"

awk -F'\t' -v ref="${REF_COND}" -v order="${COND_ORDER[*]}" 'BEGIN {
    OFS = "\t"
    nc = split(order, ord, " ")
  }
  NR > 1 { n[$2]++; s[$2] += $5; r[$2] += $4; l[$2] += $3 }
  END {
    if (!(ref in n)) { print "*** error: reference condition " ref " has no samples" > "/dev/stderr"; exit 1 }
    print "condition", "n_reps", "mean_chrM_rpm", "sum_chrM_reads", "mito_sf", "sum_reads"
    for (i = 1; i <= nc; i++) {
      c = ord[i]
      if (!(c in n)) { print "*** error: condition " c " has no samples" > "/dev/stderr"; exit 1 }
      if (s[c] <= 0) { print "*** error: condition " c " has no chrM reads" > "/dev/stderr"; exit 1 }
      printf "%s\t%d\t%.6f\t%d\t%.6f\t%d\n", c, n[c], s[c] / n[c], r[c], (s[ref] / n[ref]) / (s[c] / n[c]), l[c]
    }
  }' "${OUT_SAMPLES}" > "${OUT_SF}.tmp" || { rm -f "${OUT_SF}.tmp"; exit 1; }
mv "${OUT_SF}.tmp" "${OUT_SF}"

echo ""
cat "${OUT_SF}"

echo ""
echo "Pairwise WT/POLR2A chrM RPM ratios (each WT replicate vs each POLR2A replicate):"
awk -F'\t' 'NR > 1 { if ($2 == "WT") { w[++nw] = $5; wn[nw] = $1 } else { p[++np] = $5; pn[np] = $1 } }
  END { for (i = 1; i <= nw; i++) for (j = 1; j <= np; j++) {
          v = w[i] / p[j]; printf "  %s / %s = %.4f\n", wn[i], pn[j], v; s += v; q += v * v; k++ ; if (min == "" || v < min) min = v; if (v > max) max = v }
        m = s / k; sd = sqrt((q - k * m * m) / (k - 1))
        printf "  mean %.4f | sd %.4f | range %.4f-%.4f (%d pairs)\n", m, sd, min, max, k }' "${OUT_SAMPLES}"

if [[ -s "${EU_SF_FILE}" ]]; then
  echo ""
  echo "For comparison, EU-seq POLR2A mito factor (${EU_SF_FILE}): $(awk -F'\t' '$1 == "POLR2A_8hrs_aux" {print $5}' "${EU_SF_FILE}")"
fi

echo ""
echo "Written: ${OUT_SAMPLES}"
echo "Written: ${OUT_SF}"

##--------------------------------------------------------------## Self-checks on the live output
echo "========================================================"
echo "Self-checks"
NFAIL=0
chk() { if [[ "$1" == ok ]]; then echo "  CHECK PASS  $2"; else echo "  CHECK FAIL  $2"; NFAIL=$((NFAIL + 1)); fi; }

chk "$([[ $(( $(wc -l < "${OUT_SAMPLES}") - 1 )) -eq ${#SAMPLES[@]} ]] && echo ok)" "${#SAMPLES[@]} samples counted"
chk "$(awk -F'\t' -v m="${MIN_READS}" 'NR > 1 && $3 < m {bad++} END {print bad ? "" : "ok"}' "${OUT_SAMPLES}")" "every sample has >= ${MIN_READS} mapped reads"
chk "$(awk -F'\t' -v m="${MIN_MITO}" 'NR > 1 && $4 < m {bad++} END {print bad ? "" : "ok"}' "${OUT_SAMPLES}")" "every sample has >= ${MIN_MITO} ${MITO} reads"
chk "$(awk -F'\t' 'NR > 1 && $4 > $3 {bad++} END {print bad ? "" : "ok"}' "${OUT_SAMPLES}")" "${MITO} reads <= library in every sample"
chk "$(awk -F'\t' '$1 == "WT" {print ($5 == 1) ? "ok" : ""}' "${OUT_SF}")" "WT mito_sf is exactly 1"
chk "$(awk -F'\t' 'NR > 1 && ($5 < 0.1 || $5 > 10) {bad++} END {print bad ? "" : "ok"}' "${OUT_SF}")" "mito_sf within 0.1-10"
chk "$(awk -F'\t' 'NR > 1 {n++; if ($2 != 3) bad++} END {print (n == 2 && !bad) ? "ok" : ""}' "${OUT_SF}")" "2 conditions x 3 replicates"
chk "$(awk -F'\t' 'NR > 1 { m = $4 * 1000000.0 / $3; s[$2] += m; n[$2]++ }
  END { e = (s["WT"] / n["WT"]) / (s["POLR2A_8hrs_aux"] / n["POLR2A_8hrs_aux"]); print e }' "${OUT_SAMPLES}" | \
  awk -v f="$(awk -F'\t' '$1 == "POLR2A_8hrs_aux" {print $5}' "${OUT_SF}")" '{d = $1 - f; if (d < 0) d = -d; print (d < 1e-5) ? "ok" : ""}')" "POLR2A mito_sf recomputed from raw counts matches"
chk "$(awk -F'\t' 'FNR == NR { if (FNR > 1) { l[$2] += $3; r[$2] += $4 } ; next } FNR > 1 { if ($6 != l[$1] || $4 != r[$1]) bad++ } END {print bad ? "" : "ok"}' "${OUT_SAMPLES}" "${OUT_SF}")" "sum_reads / sum_chrM_reads equal the sums of the sample table"
echo "========================================================"
[[ ${NFAIL} -eq 0 ]] || { echo "*** ${NFAIL} self-check(s) failed" >&2; exit 1; }
