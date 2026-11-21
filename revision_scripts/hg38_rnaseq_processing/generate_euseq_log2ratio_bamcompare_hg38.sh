#!/bin/bash
#SBATCH -A p32171
#SBATCH -p genhimem
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G
#SBATCH -t 12:00:00
#SBATCH --job-name=euseq_bamcompare_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/euseq_bamcompare_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/euseq_bamcompare_%j.err
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# EU-seq log2(POLR2A / WT) bigWig on hg38 from the merged EU-seq BAMs,
# scaled with the chrM (mito) factor. Same bamCompare call as the
# EU-seq step of generate_pol2_vs_wt_bigwigs.sh (Pol II repo): raw merged
# BAMs, no NoEU subtraction, --binSize 20 --smoothLength 60
# --ignoreDuplicates --extendReads --pseudocount 1. Only the scale
# factors differ: the total-RNA-seq factor (0.747437) is replaced by
# the chrM factor of compute_scale_factors_mito_hg38.sh.
#
# Scaling, as in generate_total_bigwigs_mito_hg38.sh:
#
#   scale(cond) = eu_sf(cond) x N_read1(WT) / N_read1(cond)
#
# i.e. library-size matching to WT, then the chrM factor (WT = 1), so
# the chrM signal is the same in both conditions. eu_sf and the read1
# libraries are read from scale_factors_mito.txt and
# scale_factors_mito_samples.txt; the script stops if a merged BAM's
# read1 count does not equal the sum of its individual BAMs.
#
# Output (separate from the per-nucleotide tracks of
# generate_euseq_bigwigs_hg38.sh, which are not touched):
#   euseq/coverage/log2ratio_bamcompare/POLR2A_vs_WT.log2ratio.bamcompare.bigWig
#   (>0: more signal in POLR2A-degraded than WT)
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs
# Usage (after merge_euseq_bams_hg38.sh and compute_scale_factors_mito_hg38.sh):
#   sbatch generate_euseq_log2ratio_bamcompare_hg38.sh
##################################################################

set -o pipefail

module purge
module load samtools/1.14

##--------------------------------------------------------------## 1. Paths

BASE="${BASE:-/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq}"
EU_BAM_DIR="${EU_BAM_DIR:-${BASE}/euseq/BAM/merged}"
SF_FILE="${SF_FILE:-${BASE}/scale_factors_mito.txt}"
SAMPLES_FILE="${SAMPLES_FILE:-${BASE}/scale_factors_mito_samples.txt}"
OUT_DIR="${OUT_DIR:-${BASE}/euseq/coverage/log2ratio_bamcompare}"
OUT_BW="${OUT_DIR}/POLR2A_vs_WT.log2ratio.bamcompare.bigWig"
MITO="${MITO:-chrM}"

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-10}}"
BINSIZE=20
SMOOTH=60

die() { echo "*** error: $*" >&2; exit 1; }

POLR2A_BAM="${EU_BAM_DIR}/POLR2A_8hrs_aux_merged.sort.bam"
WT_BAM="${EU_BAM_DIR}/WT_merged.sort.bam"
for bam in "${POLR2A_BAM}" "${WT_BAM}"; do
  [[ -f "${bam}" ]] || die "merged BAM not found: ${bam} (run merge_euseq_bams_hg38.sh)"
  [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]] || die "BAM index not found for ${bam} (samtools index)"
done
[[ -s "${SF_FILE}" ]]      || die "scale factors not found: ${SF_FILE} (run compute_scale_factors_mito_hg38.sh)"
[[ -s "${SAMPLES_FILE}" ]] || die "sample table not found: ${SAMPLES_FILE} (run compute_scale_factors_mito_hg38.sh)"
[[ ! -e "${OUT_BW}" ]]     || die "output already exists, not overwriting: ${OUT_BW}"
mkdir -p "${OUT_DIR}" || die "cannot create ${OUT_DIR}"

##--------------------------------------------------------------## 2. Scale factors from the files, checked against the merged BAMs

declare -A EU_SF SUM_READS BAMS SCALE
BAMS["POLR2A_8hrs_aux"]="${POLR2A_BAM}"
BAMS["WT"]="${WT_BAM}"

echo "========================================================"
echo "$(date): EU-seq chrM scale factors (${SF_FILE}, ${SAMPLES_FILE})"
for COND in POLR2A_8hrs_aux WT; do
  EU_SF[$COND]=$(awk -F'\t' -v c="${COND}" '$1 == c {print $5}' "${SF_FILE}")
  SUM_READS[$COND]=$(awk -F'\t' -v c="${COND}" 'NR > 1 && $2 == c {s += $3} END {printf "%d", s}' "${SAMPLES_FILE}")
  [[ "${EU_SF[$COND]}" =~ ^[0-9.]+$ && "${SUM_READS[$COND]}" =~ ^[0-9]+$ && "${SUM_READS[$COND]}" -gt 0 ]] || die "could not read ${COND} from ${SF_FILE} / ${SAMPLES_FILE}"
  ## same read set as the library in compute_scale_factors_mito_hg38.sh
  MERGED=$(samtools view -c -q 1 -F 2308 -f 64 "${BAMS[$COND]}") || die "samtools failed on ${BAMS[$COND]}"
  echo "  ${COND}: merged BAM ${MERGED} read1 | sum of individual BAMs ${SUM_READS[$COND]} | eu_sf ${EU_SF[$COND]}"
  [[ "${MERGED}" -eq "${SUM_READS[$COND]}" ]] || die "${COND} merged BAM (${MERGED}) != sum of individual BAMs (${SUM_READS[$COND]}); re-run merge_euseq_bams_hg38.sh"
done
[[ $(awk -v x="${EU_SF[WT]}" 'BEGIN {print (x == 1) ? "ok" : "bad"}') == "ok" ]] || die "WT eu_sf is not 1 in ${SF_FILE}"

for COND in POLR2A_8hrs_aux WT; do
  SCALE[$COND]=$(awk -v sf="${EU_SF[$COND]}" -v nw="${SUM_READS[WT]}" -v n="${SUM_READS[$COND]}" 'BEGIN {printf "%.8f", sf * nw / n}')
  echo "  scale(${COND}) = eu_sf x N_WT / N_cond = ${SCALE[$COND]}"
done
echo "========================================================"

##--------------------------------------------------------------## 3. EU-seq log2(POLR2A / WT)

module purge
module load deeptools/3.5.6

echo "$(date): bamCompare log2(POLR2A / WT)"
bamCompare \
  --bamfile1           "${POLR2A_BAM}" \
  --bamfile2           "${WT_BAM}" \
  --scaleFactors       "${SCALE[POLR2A_8hrs_aux]}:${SCALE[WT]}" \
  --operation          log2 \
  --binSize            "${BINSIZE}" \
  --smoothLength       "${SMOOTH}" \
  --outFileName        "${OUT_BW}" \
  --numberOfProcessors "${THREADS}" \
  --ignoreDuplicates \
  --extendReads \
  --pseudocount        1 || die "bamCompare failed"

echo "  Done: ${OUT_BW}"

##--------------------------------------------------------------## 4. Self-checks on the live output

echo "========================================================"
echo "Self-checks"
NFAIL=0
chk() { if [[ "$1" == ok ]]; then echo "  CHECK PASS  $2"; else echo "  CHECK FAIL  $2"; NFAIL=$((NFAIL + 1)); fi; }

chk "$([[ -s "${OUT_BW}" ]] && echo ok)" "exists, non-empty: $(basename "${OUT_BW}")"

## chrM is the spike-in: after scaling its signal is equal in both conditions, so the mean
## log2 ratio over chrM must be ~0 (|mean| < 0.15; the factor is a mean of per-replicate RPM
## while the BAMs are pooled)
module purge
module load samtools/1.14
module load kentUtils/302
MITO_LEN=$(samtools view -H "${WT_BAM}" | awk -F'\t' -v c="${MITO}" '$1 == "@SQ" && $2 == "SN:" c {sub("LN:", "", $3); print $3}')
if [[ ! "${MITO_LEN}" =~ ^[0-9]+$ ]] || ! command -v bigWigAverageOverBed >/dev/null 2>&1; then
  chk "" "chrM check needs ${MITO} in the BAM header and bigWigAverageOverBed (kentUtils/302)"
else
  MBED="${OUT_DIR}/.${MITO}.bed"; MTAB="${OUT_DIR}/.${MITO}.tab"
  printf "%s\t0\t%s\t%s\n" "${MITO}" "${MITO_LEN}" "${MITO}" > "${MBED}"
  bigWigAverageOverBed "${OUT_BW}" "${MBED}" "${MTAB}" || die "bigWigAverageOverBed failed"
  MEAN=$(awk -F'\t' '{printf "%.4f", $5}' "${MTAB}")
  chk "$(awk -v m="${MEAN}" 'BEGIN { exit !(m == m + 0 && m > -0.15 && m < 0.15) }' && echo ok)" \
      "mean log2(POLR2A/WT) over ${MITO} = ${MEAN} (expected ~0)"
  rm -f "${MBED}" "${MTAB}"
fi
echo "========================================================"
[[ ${NFAIL} -eq 0 ]] || { echo "*** ${NFAIL} self-check(s) failed" >&2; exit 1; }
echo "$(date): done -> ${OUT_BW}"
