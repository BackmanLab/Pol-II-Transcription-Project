#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH -t 12:00:00
#SBATCH --job-name=total_bigwigs_mito_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/logs/total_bigwigs_mito_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/logs/total_bigwigs_mito_%j.err
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Total RNA-seq bigWigs on hg38 from the merged BAMs (POLR2A Rep1-3
# vs WT Rep1-3), scaled with the chrM (mito) factors of
# compute_scale_factors_mito_total_hg38.sh. Based on
# generate_bigwigs.sh (Pol II repo).
# EU-seq bigWigs are made by generate_euseq_bigwigs_hg38.sh.
#
# Scaling: the repo scaled raw coverage by a library-size factor
# (bamCoverage --scaleFactor SF, WT = 1, i.e. read counts at WT
# depth). Here
#
#   scale(cond) = mito_sf(cond) x sum_reads(WT) / sum_reads(cond)
#
# i.e. library-size matching to WT, then the chrM factor: the chrM
# signal of the scaled bigWigs is the same in both conditions (WT = 1).
# Scale factors are NOT computed here: mito_sf and the library sizes are
# read from scale_factors_mito_total.txt, and the script stops if a
# merged BAM's mapped-read count does not equal the sum of its
# individual BAMs listed there.
#
# Outputs (RNA_COV_DIR = .../total/coverage):
#   cpm/{cond}.cpm.bigWig                     CPM, no chrM factor
#   scaled_mito/{cond}.scaled.bigWig          merged BAM x scale(cond)
#   log2ratio_mito/POLR2A_vs_WT.log2ratio.bigWig  (>0: more signal in POLR2A)
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/total/logs
# Usage (after merge_totalrnaseq_bams_hg38.sh and
# compute_scale_factors_mito_total_hg38.sh):
#   sbatch generate_total_bigwigs_mito_hg38.sh
##################################################################

set -o pipefail

module purge
module load samtools/1.14

BASE="${BASE:-/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq}"
RNA_BAM_DIR="${RNA_BAM_DIR:-${BASE}/total/BAM/merged}"
RNA_COV_DIR="${RNA_COV_DIR:-${BASE}/total/coverage}"
SF_FILE="${SF_FILE:-${BASE}/scale_factors_mito_total.txt}"
MITO="${MITO:-chrM}"

mkdir -p "${RNA_COV_DIR}/cpm" "${RNA_COV_DIR}/scaled_mito" "${RNA_COV_DIR}/log2ratio_mito"

THREADS="${THREADS:-16}"
BINSIZE=20
SMOOTH=60

die() { echo "*** error: $*" >&2; exit 1; }

POLR2A_BAM="${RNA_BAM_DIR}/POLR2A_8hrs_aux_merged.sort.bam"
WT_BAM="${RNA_BAM_DIR}/WT_merged.sort.bam"
[[ -f "${POLR2A_BAM}" && -f "${WT_BAM}" ]] || die "merged BAM not found in ${RNA_BAM_DIR} (run merge_totalrnaseq_bams_hg38.sh)"
[[ -s "${SF_FILE}" ]] || die "scale factors not found: ${SF_FILE} (run compute_scale_factors_mito_total_hg38.sh)"

##--------------------------------------------------------------## 1. Scale factors from the file, checked against the merged BAMs

sf_col() {  ## condition column
  awk -F'\t' -v c="$1" -v k="$2" '$1 == c {print $k}' "${SF_FILE}"
}

declare -A SUM_READS MITO_SF BAMS MERGED_N SCALE
BAMS["POLR2A_8hrs_aux"]="${POLR2A_BAM}"
BAMS["WT"]="${WT_BAM}"

echo "========================================================"
echo "$(date): Total RNA-seq mito scale factors (from ${SF_FILE})"
for COND in POLR2A_8hrs_aux WT; do
  SUM_READS[$COND]=$(sf_col "${COND}" 6)
  MITO_SF[$COND]=$(sf_col "${COND}" 5)
  [[ "${SUM_READS[$COND]}" =~ ^[0-9]+$ && "${MITO_SF[$COND]}" =~ ^[0-9.]+$ ]] || die "could not read ${COND} from ${SF_FILE}"
  MERGED=$(samtools view -c -q 1 -F 2308 "${BAMS[$COND]}") || die "samtools failed on ${BAMS[$COND]}"
  MERGED_N[$COND]="${MERGED}"
  echo "  ${COND}: merged BAM ${MERGED} reads | sum of individual BAMs ${SUM_READS[$COND]} | mito_sf ${MITO_SF[$COND]}"
  [[ "${MERGED}" -eq "${SUM_READS[$COND]}" ]] || die "${COND} merged BAM (${MERGED}) != sum of individual BAMs (${SUM_READS[$COND]}); re-run merge_totalrnaseq_bams_hg38.sh"
done
[[ $(awk -v x="${MITO_SF[WT]}" 'BEGIN {print (x == 1) ? "ok" : "bad"}') == "ok" ]] || die "WT mito_sf is not 1 in ${SF_FILE}"

for COND in POLR2A_8hrs_aux WT; do
  SCALE[$COND]=$(awk -v sf="${MITO_SF[$COND]}" -v nw="${SUM_READS[WT]}" -v n="${SUM_READS[$COND]}" 'BEGIN {printf "%.8f", sf * nw / n}')
  echo "  scale(${COND}) = mito_sf x N_WT / N_cond = ${SCALE[$COND]}"
done

## chrM length (for the chrM signal check at the end)
MITO_LEN=$(samtools view -H "${WT_BAM}" | awk -F'\t' -v c="${MITO}" '$1 == "@SQ" && $2 == "SN:" c {sub("LN:", "", $3); print $3}')
[[ "${MITO_LEN}" =~ ^[0-9]+$ ]] || die "contig ${MITO} not in header of ${WT_BAM}"
echo "========================================================"

module purge
module load deeptools/3.5.6

##--------------------------------------------------------------## 2. bigWigs

for COND in POLR2A_8hrs_aux WT; do
  BAM="${BAMS[$COND]}"
  SF="${SCALE[$COND]}"

  echo "$(date): ${COND} CPM"
  bamCoverage --bam "${BAM}" --normalizeUsing CPM --binSize "${BINSIZE}" --smoothLength "${SMOOTH}" \
    --outFileName "${RNA_COV_DIR}/cpm/${COND}.cpm.bigWig" \
    --numberOfProcessors "${THREADS}" --ignoreDuplicates || die "bamCoverage CPM failed: ${COND}"

  echo "$(date): ${COND} mito-scaled (scale=${SF})"
  bamCoverage --bam "${BAM}" --scaleFactor "${SF}" --binSize "${BINSIZE}" --smoothLength "${SMOOTH}" \
    --outFileName "${RNA_COV_DIR}/scaled_mito/${COND}.scaled.bigWig" \
    --numberOfProcessors "${THREADS}" --ignoreDuplicates || die "bamCoverage scaled failed: ${COND}"
done

echo "$(date): log2(POLR2A / WT), mito-scaled"
bamCompare --bamfile1 "${POLR2A_BAM}" --bamfile2 "${WT_BAM}" \
  --scaleFactors "${SCALE[POLR2A_8hrs_aux]}:${SCALE[WT]}" \
  --operation log2 --pseudocount 1 --binSize "${BINSIZE}" --smoothLength "${SMOOTH}" \
  --outFileName "${RNA_COV_DIR}/log2ratio_mito/POLR2A_vs_WT.log2ratio.bigWig" \
  --numberOfProcessors "${THREADS}" --ignoreDuplicates || die "bamCompare failed"

echo ""
echo "$(date): Complete."
echo "CPM:        ${RNA_COV_DIR}/cpm/"
echo "Scaled:     ${RNA_COV_DIR}/scaled_mito/"
echo "Log2 ratio: ${RNA_COV_DIR}/log2ratio_mito/POLR2A_vs_WT.log2ratio.bigWig"

##--------------------------------------------------------------## 3. Self-checks on the live outputs
echo "========================================================"
echo "Self-checks"
NFAIL=0
chk() { if [[ "$1" == ok ]]; then echo "  CHECK PASS  $2"; else echo "  CHECK FAIL  $2"; NFAIL=$((NFAIL + 1)); fi; }

for f in "${RNA_COV_DIR}"/cpm/{POLR2A_8hrs_aux,WT}.cpm.bigWig "${RNA_COV_DIR}"/scaled_mito/{POLR2A_8hrs_aux,WT}.scaled.bigWig \
         "${RNA_COV_DIR}/log2ratio_mito/POLR2A_vs_WT.log2ratio.bigWig"; do
  chk "$([[ -s "${f}" ]] && echo ok)" "exists, non-empty: $(basename "${f}")"
done

## chrM signal: scaled / CPM over chrM must equal scale x N_cond / 1e6 (CPM = coverage x 1e6/N, scaled = coverage x scale),
## i.e. the factor in scale_factors_mito_total.txt is what was applied. The chrM signal of the two scaled
## bigWigs is printed as well (it is equal up to the mean-of-replicate-RPM vs pooled-reads difference).
module load kentUtils/302
if ! command -v bigWigAverageOverBed >/dev/null 2>&1; then
  chk "" "bigWigAverageOverBed not available (kentUtils/302)"
else
  MBED="${RNA_COV_DIR}/.${MITO}.bed"
  printf "%s\t0\t%s\t%s\n" "${MITO}" "${MITO_LEN}" "${MITO}" > "${MBED}"
  bw_mito_sum() {  ## bigWig -> summed signal over chrM
    local tab="${RNA_COV_DIR}/.${MITO}.tab"
    bigWigAverageOverBed "$1" "${MBED}" "${tab}" || { echo NA; return; }
    awk -F'\t' '{printf "%.6f", $4}' "${tab}"
  }
  declare -A MSUM
  for COND in POLR2A_8hrs_aux WT; do
    S_SC=$(bw_mito_sum "${RNA_COV_DIR}/scaled_mito/${COND}.scaled.bigWig")
    S_CPM=$(bw_mito_sum "${RNA_COV_DIR}/cpm/${COND}.cpm.bigWig")
    MSUM[$COND]="${S_SC}"
    RATIO=$(awk -v a="${S_SC}" -v b="${S_CPM}" -v sc="${SCALE[$COND]}" -v n="${SUM_READS[$COND]}" \
      'BEGIN { if (a + 0 > 0 && b + 0 > 0 && sc > 0 && n > 0) printf "%.5f", (a / b) / (sc * n / 1e6); else print "NA" }')
    chk "$([[ "${RATIO}" != "NA" ]] && awk -v r="${RATIO}" 'BEGIN { exit !(r >= 0.99 && r <= 1.01) }' && echo ok)" \
        "${COND}: ${MITO} scaled/CPM signal ratio / (scale x reads/1e6) = ${RATIO} (within 1%)"
  done
  echo "  INFO  ${MITO} signal in the scaled bigWigs, POLR2A / WT = $(awk -v a="${MSUM[POLR2A_8hrs_aux]}" -v b="${MSUM[WT]}" 'BEGIN { if (b > 0) printf "%.4f", a / b; else print "NA" }') (expected ~1)"
  rm -f "${MBED}" "${RNA_COV_DIR}/.${MITO}.tab"
fi
echo "========================================================"
[[ ${NFAIL} -eq 0 ]] || { echo "*** ${NFAIL} self-check(s) failed" >&2; exit 1; }
