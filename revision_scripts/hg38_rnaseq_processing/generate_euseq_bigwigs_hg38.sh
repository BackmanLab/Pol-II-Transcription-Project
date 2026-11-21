#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH -t 12:00:00
#SBATCH --job-name=euseq_bw_bgr_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/euseq_bigwigs_bgr_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/euseq_bigwigs_bgr_%j.err
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# EU-seq bigWigs on hg38 built from the SAME per-nucleotide
# NoEU-subtracted RPM bedGraphs the RPKM matrices are built from
# (euseq_subtract_hg38.sh -> euseq_pipeline/subtracted/*.sub.bgr).
#
# Units: RPKM per nucleotide. Each replicate bedGraph is multiplied by
# 1000 / (its mean aligned bases per read, library_sizes.txt), so the
# mean of a bigWig over any interval equals that interval's RPKM as
# computed in euseq_rpkm_hg38.sh (length-weighted).
#
# Per condition and strand:
#   1. Replicate bedGraphs x (1000 / mean aligned bases per read),
#      averaged per nucleotide (bedtools unionbedg, missing = 0)
#   2. x eu_sf for the condition (scale_factors_mito.txt: chrM-derived,
#      mean chrM RPM in WT / mean chrM RPM in the condition; WT = 1, from
#      compute_scale_factors_mito_hg38.sh)
#   3. Values are >= 0 (NoEU subtraction was floored at 0 per
#      nucleotide); positions with no signal are not written
#   4. Combined track = plus + minus; bedGraphToBigWig
# POLR2A vs WT: log2((POLR2A + pc) / (WT + pc)) per nucleotide.
# pc = signal of one read at the mean library size, in the same
# RPKM-per-nucleotide units (printed below). Override: PSEUDO=x.
# The log2 ratio is a ratio and is negative where POLR2A < WT.
#
# Outputs (EU_COV_DIR = .../euseq/coverage):
#   subtracted/{POLR2A_8hrs_aux,WT}.subtracted.{plus,minus}.bigWig
#   subtracted/{POLR2A_8hrs_aux,WT}.subtracted.bigWig          (plus+minus)
#   log2ratio/POLR2A_vs_WT.log2ratio.{plus,minus}.bigWig
#   log2ratio/POLR2A_vs_WT.log2ratio.bigWig                    (combined)
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs
# Usage (after compute_scale_factors_mito_hg38.sh and euseq_subtract_hg38.sh):
#   sbatch generate_euseq_bigwigs_hg38.sh
##################################################################

set -o pipefail

module purge
module load bedtools/2.29.2

##--------------------------------------------------------------## 1. Paths

BASE="${BASE:-/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq}"
SUB_DIR="${SUB_DIR:-${BASE}/euseq/euseq_pipeline/subtracted}"
LIBSIZE_FILE="${LIBSIZE_FILE:-${BASE}/euseq/euseq_pipeline/library_sizes.txt}"
SF_FILE="${SF_FILE:-${BASE}/scale_factors_mito.txt}"
CHROM_SIZES="${CHROM_SIZES:-${BASE}/reference/hg38.chrom.sizes}"
EU_COV_DIR="${EU_COV_DIR:-${BASE}/euseq/coverage}"

BGR_OUT="${EU_COV_DIR}/bedgraphs"
mkdir -p "${BGR_OUT}" "${EU_COV_DIR}/subtracted" "${EU_COV_DIR}/log2ratio"

TMPSORT="${TMPDIR:-/tmp}"
SORT_MEM="${SORT_MEM:-32G}"

die() { echo "*** error: $*" >&2; exit 1; }

POLR2A_NAMES=(POLR2A_1 POLR2A_2 POLR2A_3)
WT_NAMES=(WT_1 WT_2 WT_3)

[[ -s "${CHROM_SIZES}" ]]  || die "chrom.sizes not found: ${CHROM_SIZES}"
[[ -s "${LIBSIZE_FILE}" ]] || die "library sizes not found: ${LIBSIZE_FILE}"
[[ -s "${SF_FILE}" ]]      || die "scale factors not found: ${SF_FILE} (run compute_scale_factors_mito_hg38.sh)"
for n in "${POLR2A_NAMES[@]}" "${WT_NAMES[@]}"; do
  for s in plus minus; do
    [[ -s "${SUB_DIR}/${n}_${s}.sub.bgr" ]] || die "missing ${SUB_DIR}/${n}_${s}.sub.bgr"
  done
done

##--------------------------------------------------------------## 2. Scale factors, per-sample RPKM multipliers, pseudocount

SF_POLR2A=$(awk -F'\t' '$1 == "POLR2A_8hrs_aux" {print $5}' "${SF_FILE}")
SF_WT=$(awk -F'\t' '$1 == "WT" {print $5}' "${SF_FILE}")
[[ "${SF_POLR2A}" =~ ^[0-9.]+$ && "${SF_WT}" =~ ^[0-9.]+$ ]] || die "could not read eu_sf from ${SF_FILE}"
[[ $(awk -v x="${SF_WT}" 'BEGIN {print (x == 1) ? "ok" : "bad"}') == "ok" ]] || die "WT eu_sf is not 1 in ${SF_FILE}"

## k = 1000 / mean aligned bases per read
kmult() {
  local n="$1" ml
  ml=$(awk -v s="${n}" '$1 == s {print $3}' "${LIBSIZE_FILE}")
  [[ -n "${ml}" ]] || die "no mean aligned length for ${n} in ${LIBSIZE_FILE}"
  awk -v ml="${ml}" 'BEGIN {printf "%.8f", 1000.0 / ml}'
}
K_POLR2A=""; K_WT=""
for n in "${POLR2A_NAMES[@]}"; do K_POLR2A="${K_POLR2A} $(kmult "${n}")"; done
for n in "${WT_NAMES[@]}";     do K_WT="${K_WT} $(kmult "${n}")"; done

if [[ -z "${PSEUDO:-}" ]]; then
  PSEUDO=$(awk -v ids="${POLR2A_NAMES[*]} ${WT_NAMES[*]}" '
    BEGIN { n = split(ids, a, " "); for (i = 1; i <= n; i++) want[a[i]] = 1 }
    ($1 in want) { s += (1000000.0 / $2) * (1000.0 / $3); k++ }
    END { if (k != n) { print "NA"; exit }
          printf "%.8f", s / k }' "${LIBSIZE_FILE}")
fi
[[ "${PSEUDO}" != "NA" && -n "${PSEUDO}" ]] || die "pseudocount failed; ${LIBSIZE_FILE} must list all 6 samples"

echo "========================================================"
echo "$(date): Scale factors (from ${SF_FILE})"
echo "  eu_sf  POLR2A_8hrs_aux = ${SF_POLR2A} | WT = ${SF_WT}"
echo "  RPKM multipliers 1000/mean_len  POLR2A:${K_POLR2A} | WT:${K_WT}"
echo "  Pseudocount for log2 ratio (one read, RPKM-per-nt units): ${PSEUDO}"
echo "========================================================"

##--------------------------------------------------------------## 3. bedGraph functions

## mean over replicates of (value x k_i), x sf; zeros are not written
mean_scale() {
  local out="$1" sf="$2" ks="$3"; shift 3
  local files=("$@")
  echo "$(date): mean x ${sf} -> $(basename "${out}")"
  bedtools unionbedg -i "${files[@]}" -filler 0 | \
  awk -v n="${#files[@]}" -v sf="${sf}" -v ks="${ks}" 'BEGIN {OFS = "\t"; split(ks, k, " ")}
    NF != 3 + n { print "*** error: unionbedg line has " NF " fields, expected " 3 + n > "/dev/stderr"; exit 1 }
    { s = 0; for (i = 1; i <= n; i++) s += $(3 + i) * k[i]
      v = s / n * sf
      if (v > 0) printf "%s\t%s\t%s\t%.6g\n", $1, $2, $3, v }' > "${out}" || die "mean_scale failed: ${out}"
}

## sum of two bedGraphs
sum_bgr() {
  local out="$1" a="$2" b="$3"
  echo "$(date): plus + minus -> $(basename "${out}")"
  bedtools unionbedg -i "${a}" "${b}" -filler 0 | \
  awk 'BEGIN {OFS = "\t"}
    NF != 5 { print "*** error: unionbedg line has " NF " fields, expected 5" > "/dev/stderr"; exit 1 }
    { v = $4 + $5; if (v > 0) printf "%s\t%s\t%s\t%.6g\n", $1, $2, $3, v }' > "${out}" || die "sum_bgr failed: ${out}"
}

## log2((a + pc) / (b + pc)), a = POLR2A, b = WT
log2_bgr() {
  local out="$1" a="$2" b="$3"
  echo "$(date): log2 ratio -> $(basename "${out}")"
  bedtools unionbedg -i "${a}" "${b}" -filler 0 | \
  awk -v pc="${PSEUDO}" 'BEGIN {OFS = "\t"}
    NF != 5 { print "*** error: unionbedg line has " NF " fields, expected 5" > "/dev/stderr"; exit 1 }
    { printf "%s\t%s\t%s\t%.6g\n", $1, $2, $3, log(($4 + pc) / ($5 + pc)) / log(2) }' > "${out}" || die "log2_bgr failed: ${out}"
}

##--------------------------------------------------------------## 4. Build bedGraphs

echo "========================================================"
echo "$(date): Building bedGraphs"
echo "========================================================"

for s in plus minus; do
  P_FILES=(); W_FILES=()
  for n in "${POLR2A_NAMES[@]}"; do P_FILES+=("${SUB_DIR}/${n}_${s}.sub.bgr"); done
  for n in "${WT_NAMES[@]}";     do W_FILES+=("${SUB_DIR}/${n}_${s}.sub.bgr"); done
  mean_scale "${BGR_OUT}/POLR2A_8hrs_aux.subtracted.${s}.bgr" "${SF_POLR2A}" "${K_POLR2A}" "${P_FILES[@]}"
  mean_scale "${BGR_OUT}/WT.subtracted.${s}.bgr"              "${SF_WT}"     "${K_WT}"     "${W_FILES[@]}"
done
for c in POLR2A_8hrs_aux WT; do
  sum_bgr "${BGR_OUT}/${c}.subtracted.bgr" \
          "${BGR_OUT}/${c}.subtracted.plus.bgr" "${BGR_OUT}/${c}.subtracted.minus.bgr"
done

log2_bgr "${BGR_OUT}/POLR2A_vs_WT.log2ratio.plus.bgr" \
         "${BGR_OUT}/POLR2A_8hrs_aux.subtracted.plus.bgr"  "${BGR_OUT}/WT.subtracted.plus.bgr"
log2_bgr "${BGR_OUT}/POLR2A_vs_WT.log2ratio.minus.bgr" \
         "${BGR_OUT}/POLR2A_8hrs_aux.subtracted.minus.bgr" "${BGR_OUT}/WT.subtracted.minus.bgr"
log2_bgr "${BGR_OUT}/POLR2A_vs_WT.log2ratio.bgr" \
         "${BGR_OUT}/POLR2A_8hrs_aux.subtracted.bgr"       "${BGR_OUT}/WT.subtracted.bgr"

##--------------------------------------------------------------## 5. bedGraph -> bigWig

module purge
module load kentUtils/302

to_bigwig() {
  local bgr="$1" bw="$2"
  [[ -s "${bgr}" ]] || die "empty bedGraph: ${bgr}"
  echo "$(date): bigWig -> $(basename "${bw}")"
  LC_ALL=C sort -S "${SORT_MEM}" -T "${TMPSORT}" -k1,1 -k2,2n "${bgr}" > "${bgr}.sorted" || die "sort failed: ${bgr}"
  bedGraphToBigWig "${bgr}.sorted" "${CHROM_SIZES}" "${bw}" || die "bedGraphToBigWig failed: ${bw}"
  rm -f "${bgr}.sorted"
}

echo "========================================================"
echo "$(date): Converting to bigWig"
echo "========================================================"

for c in POLR2A_8hrs_aux WT; do
  to_bigwig "${BGR_OUT}/${c}.subtracted.plus.bgr"  "${EU_COV_DIR}/subtracted/${c}.subtracted.plus.bigWig"
  to_bigwig "${BGR_OUT}/${c}.subtracted.minus.bgr" "${EU_COV_DIR}/subtracted/${c}.subtracted.minus.bigWig"
  to_bigwig "${BGR_OUT}/${c}.subtracted.bgr"       "${EU_COV_DIR}/subtracted/${c}.subtracted.bigWig"
done
to_bigwig "${BGR_OUT}/POLR2A_vs_WT.log2ratio.plus.bgr"  "${EU_COV_DIR}/log2ratio/POLR2A_vs_WT.log2ratio.plus.bigWig"
to_bigwig "${BGR_OUT}/POLR2A_vs_WT.log2ratio.minus.bgr" "${EU_COV_DIR}/log2ratio/POLR2A_vs_WT.log2ratio.minus.bigWig"
to_bigwig "${BGR_OUT}/POLR2A_vs_WT.log2ratio.bgr"       "${EU_COV_DIR}/log2ratio/POLR2A_vs_WT.log2ratio.bigWig"

##--------------------------------------------------------------## 6. Summary

echo ""
echo "$(date): bigWigs complete."
echo "SF_POLR2A = ${SF_POLR2A} | SF_WT = ${SF_WT} | pseudocount = ${PSEUDO}"
echo "Subtracted bigWigs: ${EU_COV_DIR}/subtracted/"
echo "Log2 ratio bigWigs: ${EU_COV_DIR}/log2ratio/"
echo "Intermediate bedGraphs (safe to delete): ${BGR_OUT}/"

##--------------------------------------------------------------## 7. Self-checks on the live outputs
echo "========================================================"
echo "Self-checks"
NFAIL=0
chk() { if [[ "$1" == ok ]]; then echo "  CHECK PASS  $2"; else echo "  CHECK FAIL  $2"; NFAIL=$((NFAIL + 1)); fi; }

REG_DIR="${REG_DIR:-${BASE}/euseq/euseq_pipeline/regions}"
GENE_MAT="${GENE_MAT:-${BASE}/euseq/euseq_pipeline/counts/euseq_rpkm_analysis.txt}"

ALL_BW=()
for c in POLR2A_8hrs_aux WT; do
  for x in subtracted.plus subtracted.minus subtracted; do ALL_BW+=("${EU_COV_DIR}/subtracted/${c}.${x}.bigWig"); done
done
for x in log2ratio.plus log2ratio.minus log2ratio; do ALL_BW+=("${EU_COV_DIR}/log2ratio/POLR2A_vs_WT.${x}.bigWig"); done
miss=0; for f in "${ALL_BW[@]}"; do [[ -s "${f}" ]] || { miss=$((miss + 1)); echo "    missing/empty: ${f}"; }; done
chk "$([[ ${miss} -eq 0 ]] && echo ok)" "all ${#ALL_BW[@]} bigWigs exist and are non-empty"

neg=$(cat "${BGR_OUT}"/*.subtracted*.bgr | awk '$4 < 0' | wc -l)
chk "$([[ ${neg} -eq 0 ]] && echo ok)" "subtracted bedGraphs contain no negative values (floored at 0)"

if command -v bigWigAverageOverBed >/dev/null 2>&1 && [[ -s "${GENE_MAT}" && -s "${REG_DIR}/genes_plus.bed" ]]; then
  ## the track must carry the same signal as the RPKM matrix: mean over a gene body == gene RPKM x eu_sf
  for c in POLR2A_8hrs_aux:"${SF_POLR2A}":"${POLR2A_NAMES[*]}" WT:"${SF_WT}":"${WT_NAMES[*]}"; do
    IFS=':' read -r cond sf cols <<< "${c}"
    for s in plus minus; do
      bigWigAverageOverBed "${EU_COV_DIR}/subtracted/${cond}.subtracted.${s}.bigWig" "${REG_DIR}/genes_${s}.bed" "${BGR_OUT}/${cond}.${s}.avg.tab" || die "bigWigAverageOverBed failed"
      res=$(awk -F'\t' -v sf="${sf}" -v cols="${cols}" '
        FILENAME ~ /avg\.tab$/ { avg[$1] = $5; next }
        FNR == 1 { nc = split(cols, cn, " "); for (i = 1; i <= NF; i++) for (j = 1; j <= nc; j++) if ($i == cn[j]) ix[j] = i; next }
        ($1 in avg) { m = 0; for (j = 1; j <= nc; j++) m += $ix[j]; m = m / nc * sf; n++
                      d = avg[$1] - m; if (d < 0) d = -d
                      if (d > 1e-3 * (m > 1e-3 ? m : 1e-3)) bad++ }
        END { print (n > 0 && bad == 0) ? "ok" : "bad " bad + 0 "/" n }' "${BGR_OUT}/${cond}.${s}.avg.tab" "${GENE_MAT}")
      chk "$([[ "${res}" == ok ]] && echo ok)" "${cond} ${s}: bigWig mean over each gene == gene RPKM x eu_sf, all genes (${res})"
      rm -f "${BGR_OUT}/${cond}.${s}.avg.tab"
    done
  done
else
  echo "  (bigWigAverageOverBed or the RPKM matrix not available -- bigWig-vs-RPKM check skipped)"
fi
echo "========================================================"
[[ ${NFAIL} -eq 0 ]] || { echo "*** ${NFAIL} self-check(s) failed" >&2; exit 1; }
