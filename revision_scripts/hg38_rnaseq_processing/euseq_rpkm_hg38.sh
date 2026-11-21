#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --job-name=euseq_rpkm_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/euseq_rpkm_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/euseq_rpkm_%j.err
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# EU-seq RPKM on hg38 from the NoEU-subtracted RPM
# bedGraphs (euseq_subtract_hg38.sh). Replaces the gene-body,
# exon, intron and DoG RPKM steps of euseq_subtract_count_hg38.sh,
# euseq_exon_rpkm_hg38.sh and euseq_dog_intron_count_hg38.sh.
#
# RPKM of a region set (gene body, merged exons, introns, DoG):
#
#   signal   = sum over bedGraph intervals of  value * overlap_bp
#              (length-weighted integral of RPM coverage over the
#              strand-matched regions of the gene)
#   RPKM     = signal / mean_aligned_bases_per_read / region_kb
#
# Because RPM coverage integrates to (reads * aligned bases / library
# * 1e6), dividing by the mean aligned bases per read gives
# reads / library * 1e6, i.e. reads per kb per million.
# The mean aligned bases per read comes from library_sizes.txt.
#
# Regions (built here from the GTF):
#   gene body : min exon start - max exon end per gene_id
#   exons     : exons merged per gene_id
#   introns   : gene body minus all exons on the same strand
#   DoG       : DOG_WINDOW bp downstream of TES, windows overlapping
#               any same-strand gene removed, clipped to chrom ends
#
# Outputs (PIPELINE_DIR = .../euseq/euseq_pipeline):
#   counts/euseq_rpkm_analysis.txt             gene body, 6 samples
#   counts/euseq_rpkm_qc.txt                   gene body, 6 samples + NoEU
#   exon_counts/euseq_exon_rpkm.txt
#   intron_dog/intron_counts/euseq_intron_rpkm.txt
#   intron_dog/dog_counts/euseq_dog_rpkm.txt
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs
# Usage (after euseq_subtract_hg38.sh):
#   sbatch euseq_rpkm_hg38.sh
##################################################################

set -o pipefail

module purge
module load samtools/1.6
module load bedtools/2.29.2

##--------------------------------------------------------------## 1. Paths

BASE="${BASE:-/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq}"
PIPELINE_DIR="${PIPELINE_DIR:-${BASE}/euseq/euseq_pipeline}"
GTF="${GTF:-${BASE}/reference/Homo_sapiens.GRCh38.115.chr.gtf}"
CHROM_SIZES="${CHROM_SIZES:-${BASE}/reference/hg38.chrom.sizes}"
DOG_WINDOW="${DOG_WINDOW:-5000}"

BGR_DIR="${PIPELINE_DIR}/bedgraphs"
SUB_DIR="${PIPELINE_DIR}/subtracted"
LIBSIZE_FILE="${PIPELINE_DIR}/library_sizes.txt"
REG_DIR="${PIPELINE_DIR}/regions"
GENE_DIR="${PIPELINE_DIR}/counts"
EXON_DIR="${PIPELINE_DIR}/exon_counts"
INTRON_DIR="${PIPELINE_DIR}/intron_dog/intron_counts"
DOG_DIR="${PIPELINE_DIR}/intron_dog/dog_counts"
TMP_DIR="${PIPELINE_DIR}/rpkm_tmp"

mkdir -p "${REG_DIR}" "${GENE_DIR}" "${EXON_DIR}" "${INTRON_DIR}" "${DOG_DIR}" "${TMP_DIR}"

die() { echo "*** error: $*" >&2; exit 1; }

[[ -s "${GTF}" ]]          || die "GTF not found: ${GTF}"
[[ -s "${CHROM_SIZES}" ]]  || die "chrom.sizes not found: ${CHROM_SIZES}"
[[ -s "${LIBSIZE_FILE}" ]] || die "library sizes not found: ${LIBSIZE_FILE} (run euseq_subtract_hg38.sh first)"

ALL_NAMES=("POLR2A_1" "POLR2A_2" "POLR2A_3" "WT_1" "WT_2" "WT_3")
NOEU_NAME="NoEU"

for n in "${ALL_NAMES[@]}"; do
  for s in plus minus; do
    [[ -s "${SUB_DIR}/${n}_${s}.sub.bgr" ]] || die "missing ${SUB_DIR}/${n}_${s}.sub.bgr"
  done
done
for s in plus minus; do
  [[ -s "${BGR_DIR}/${NOEU_NAME}_${s}.norm.bgr" ]] || die "missing ${BGR_DIR}/${NOEU_NAME}_${s}.norm.bgr"
done

##--------------------------------------------------------------## 2. Region BEDs from the GTF

GENE_BED="${REG_DIR}/hg38_genes.bed"
EXON_BED="${REG_DIR}/hg38_exons_merged.bed"

## gene body = min exon start .. max exon end per gene_id (Ensembl GTF; only
## exon rows are used, chr1-22/X/Y only: chrM and unplaced contigs are excluded). Exons on a different chromosome/strand than the
## gene's first exon are skipped and counted.
echo "$(date): Building gene body BED..."
grep -v "^#" "${GTF}" | \
awk 'BEGIN {FS = "\t"; OFS = "\t"}
  $3 == "exon" && $1 ~ /^chr[0-9XY]+/ {
    if (!match($9, /gene_id "[^"]+"/)) next
    id = substr($9, RSTART + 9, RLENGTH - 10)
    s = $4 - 1; e = $5
    if (!(id in chr)) { chr[id] = $1; str[id] = $7; st[id] = s; en[id] = e; order[++n] = id; next }
    if ($1 != chr[id] || $7 != str[id]) { skipped++; next }
    if (s < st[id]) st[id] = s
    if (e > en[id]) en[id] = e
  }
  END {
    for (i = 1; i <= n; i++) { g = order[i]; print chr[g], st[g], en[g], g, 0, str[g] }
    if (skipped > 0) print "  exons skipped (different chr/strand within gene_id): " skipped > "/dev/stderr"
  }' | LC_ALL=C sort -k1,1 -k2,2n > "${GENE_BED}" || die "gene BED failed"
[[ -s "${GENE_BED}" ]] || die "gene BED is empty -- check GTF gene_id parsing"
echo "  Genes: $(wc -l < "${GENE_BED}")"
[[ $(cut -f4 "${GENE_BED}" | sort | uniq -d | wc -l) -eq 0 ]] || die "duplicate gene_id in gene BED"

## exons merged within each gene_id (overlapping or book-ended exons of the
## same gene are joined; exons of different genes are never joined)
echo "$(date): Building per-gene merged exon BED..."
grep -v "^#" "${GTF}" | \
awk 'BEGIN {FS = "\t"; OFS = "\t"}
  $3 == "exon" && $1 ~ /^chr[0-9XY]+/ {
    if (!match($9, /gene_id "[^"]+"/)) next
    print $1, $4 - 1, $5, substr($9, RSTART + 9, RLENGTH - 10), 0, $7
  }' | LC_ALL=C sort -k4,4 -k1,1 -k2,2n | \
awk 'BEGIN {OFS = "\t"}
  function flush() { if (g != "") print c, s, e, g, 0, t }
  { if ($4 == g && $1 == c && $6 == t && $2 <= e) { if ($3 > e) e = $3 }
    else { flush(); g = $4; c = $1; s = $2; e = $3; t = $6 } }
  END { flush() }' | LC_ALL=C sort -k1,1 -k2,2n > "${EXON_BED}" || die "exon BED failed"
[[ -s "${EXON_BED}" ]] || die "exon BED is empty"
echo "  Merged exon intervals: $(wc -l < "${EXON_BED}")"

split_strand() {  ## bed prefix  -> prefix_plus.bed / prefix_minus.bed
  awk '$6 == "+"' "$1" > "$2_plus.bed"
  awk '$6 == "-"' "$1" > "$2_minus.bed"
}
split_strand "${GENE_BED}" "${REG_DIR}/genes"
split_strand "${EXON_BED}" "${REG_DIR}/exons"

echo "$(date): Building intron BEDs (gene body minus same-strand exons)..."
for s in plus minus; do
  bedtools subtract -a "${REG_DIR}/genes_${s}.bed" -b "${REG_DIR}/exons_${s}.bed" -s | \
    LC_ALL=C sort -k1,1 -k2,2n > "${REG_DIR}/introns_${s}.bed" || die "intron subtract failed (${s})"
done
echo "  Intron intervals plus/minus: $(wc -l < "${REG_DIR}/introns_plus.bed") / $(wc -l < "${REG_DIR}/introns_minus.bed")"

echo "$(date): Building ${DOG_WINDOW} bp DoG windows..."
awk -v win="${DOG_WINDOW}" 'BEGIN {OFS = "\t"}
  FNR == NR { size[$1] = $2; next }
  $6 == "+" { s = $3; e = $3 + win; if (e > size[$1]) e = size[$1]; if (s < e) print $1, s, e, $4, 0, $6 }
  $6 == "-" { s = $2 - win; e = $2; if (s < 0) s = 0; if (s < e) print $1, s, e, $4, 0, $6 }' \
  "${CHROM_SIZES}" "${GENE_BED}" | LC_ALL=C sort -k1,1 -k2,2n > "${REG_DIR}/dog_raw.bed" || die "DoG windows failed"
bedtools intersect -a "${REG_DIR}/dog_raw.bed" -b "${GENE_BED}" -s -v | \
  LC_ALL=C sort -k1,1 -k2,2n > "${REG_DIR}/dog_clean.bed" || die "DoG filter failed"
split_strand "${REG_DIR}/dog_clean.bed" "${REG_DIR}/dog"
echo "  DoG windows raw / clean: $(wc -l < "${REG_DIR}/dog_raw.bed") / $(wc -l < "${REG_DIR}/dog_clean.bed")"

##--------------------------------------------------------------## 3. RPKM = integral / mean aligned bases / kb

mean_len_of() {
  local n
  n=$(awk -v s="$1" '$1 == s {print $3}' "${LIBSIZE_FILE}")
  [[ -n "${n}" ]] || die "no mean aligned length for $1 in ${LIBSIZE_FILE}"
  echo "${n}"
}

## region_rpkm REGIONS.bed BEDGRAPH MEAN_LEN  ->  "region_id <tab> rpkm" for every
## region id in REGIONS.bed (0 when no signal), summed over all of the id's intervals
region_rpkm() {
  local regions="$1" bgr="$2" mean_len="$3"
  [[ -s "${regions}" ]] || return 0
  bedtools intersect -a "${regions}" -b "${bgr}" -wo | \
  awk -v ml="${mean_len}" 'BEGIN {OFS = "\t"}
    FNR == NR { len[$4] += $3 - $2; if (!($4 in seen)) { seen[$4] = 1; order[++n] = $4 }; next }
    NF != 11 { print "*** error: intersect line has " NF " fields, expected 11" > "/dev/stderr"; exit 1 }
    { sig[$4] += $10 * $11 }
    END { for (i = 1; i <= n; i++) { g = order[i]; printf "%s\t%.10f\n", g, sig[g] / ml / (len[g] / 1000.0) } }' \
    "${regions}" - | LC_ALL=C sort -k1,1
}

## sample_rpkm NAME PLUS_BGR MINUS_BGR PLUS_REGIONS MINUS_REGIONS OUT
sample_rpkm() {
  local name="$1" pbgr="$2" mbgr="$3" preg="$4" mreg="$5" out="$6"
  local ml; ml=$(mean_len_of "${name}")
  echo "$(date): ${name}: $(basename "${out}") (mean aligned bases/read ${ml})"
  { region_rpkm "${preg}" "${pbgr}" "${ml}"; region_rpkm "${mreg}" "${mbgr}" "${ml}"; } | \
    LC_ALL=C sort -k1,1 > "${out}.tmp" || die "RPKM failed: ${name} ${out}"
  [[ -s "${out}.tmp" ]] || die "empty RPKM output for ${name} ${out}"
  mv "${out}.tmp" "${out}"
  echo "  regions: $(wc -l < "${out}") | with RPKM > 0: $(awk '$2 > 0' "${out}" | wc -l)"
}

## build_matrix OUT_FILE  FILE1 FILE2 ...  (column names from -v NAMES)
build_matrix() {
  local out="$1" names="$2"; shift 2
  awk -v names="${names}" 'BEGIN {OFS = "\t"; nf_ = split(names, nm, " ")}
    FNR == 1 { f++ }
    { v[f, $1] = $2; cnt[f]++; if (f == 1) order[++g] = $1 }
    END {
      if (f != nf_) { print "*** error: " f " files but " nf_ " names" > "/dev/stderr"; exit 1 }
      for (k = 2; k <= f; k++) if (cnt[k] != cnt[1]) {
        print "*** error: file " k " has " cnt[k] " rows, file 1 has " cnt[1] > "/dev/stderr"; exit 1 }
      printf "gene_id"; for (k = 1; k <= f; k++) printf "\t%s", nm[k]; printf "\n"
      for (i = 1; i <= g; i++) {
        printf "%s", order[i]
        for (k = 1; k <= f; k++) {
          if (!((k, order[i]) in v)) { print "*** error: " order[i] " missing in file " k > "/dev/stderr"; exit 1 }
          printf "\t%s", v[k, order[i]]
        }
        printf "\n"
      }
    }' "$@" > "${out}.tmp" || die "matrix failed: ${out}"
  mv "${out}.tmp" "${out}"
  echo "  Written: ${out} ($(( $(wc -l < "${out}") - 1 )) rows)"
}

run_region_set() {  ## label  plus_regions  minus_regions  out_dir  matrix_file
  local label="$1" preg="$2" mreg="$3" outdir="$4" matrix="$5"
  echo ""
  echo "========================================================"
  echo "$(date): ${label} RPKM"
  echo "========================================================"
  local files=() n
  for n in "${ALL_NAMES[@]}"; do
    sample_rpkm "${n}" "${SUB_DIR}/${n}_plus.sub.bgr" "${SUB_DIR}/${n}_minus.sub.bgr" \
                "${preg}" "${mreg}" "${TMP_DIR}/${label}_${n}.txt"
    files+=("${TMP_DIR}/${label}_${n}.txt")
  done
  build_matrix "${matrix}" "${ALL_NAMES[*]}" "${files[@]}"
  if [[ "${label}" == "gene" ]]; then
    sample_rpkm "${NOEU_NAME}" "${BGR_DIR}/${NOEU_NAME}_plus.norm.bgr" "${BGR_DIR}/${NOEU_NAME}_minus.norm.bgr" \
                "${preg}" "${mreg}" "${TMP_DIR}/${label}_${NOEU_NAME}.txt"
    build_matrix "${GENE_DIR}/euseq_rpkm_qc.txt" "${ALL_NAMES[*]} ${NOEU_NAME}" "${files[@]}" "${TMP_DIR}/${label}_${NOEU_NAME}.txt"
  fi
}

run_region_set gene   "${REG_DIR}/genes_plus.bed"  "${REG_DIR}/genes_minus.bed"  "${GENE_DIR}"   "${GENE_DIR}/euseq_rpkm_analysis.txt"
run_region_set exon   "${REG_DIR}/exons_plus.bed"  "${REG_DIR}/exons_minus.bed"  "${EXON_DIR}"   "${EXON_DIR}/euseq_exon_rpkm.txt"
run_region_set intron "${REG_DIR}/introns_plus.bed" "${REG_DIR}/introns_minus.bed" "${INTRON_DIR}" "${INTRON_DIR}/euseq_intron_rpkm.txt"
run_region_set dog    "${REG_DIR}/dog_plus.bed"    "${REG_DIR}/dog_minus.bed"    "${DOG_DIR}"    "${DOG_DIR}/euseq_dog_rpkm.txt"

rm -rf "${TMP_DIR}"

echo ""
echo "$(date): RPKM complete."
echo "Gene body:  ${GENE_DIR}/euseq_rpkm_analysis.txt (+ euseq_rpkm_qc.txt)"
echo "Exon:       ${EXON_DIR}/euseq_exon_rpkm.txt"
echo "Intron:     ${INTRON_DIR}/euseq_intron_rpkm.txt"
echo "DoG:        ${DOG_DIR}/euseq_dog_rpkm.txt"

##--------------------------------------------------------------## 4. Self-checks on the live outputs
echo "========================================================"
echo "Self-checks"
NFAIL=0
chk() { if [[ "$1" == ok ]]; then echo "  CHECK PASS  $2"; else echo "  CHECK FAIL  $2"; NFAIL=$((NFAIL + 1)); fi; }

GENE_M="${GENE_DIR}/euseq_rpkm_analysis.txt"; EXON_M="${EXON_DIR}/euseq_exon_rpkm.txt"
INTRON_M="${INTRON_DIR}/euseq_intron_rpkm.txt"; DOG_M="${DOG_DIR}/euseq_dog_rpkm.txt"
HDR=$(printf 'gene_id\t%s' "$(IFS=$'\t'; echo "${ALL_NAMES[*]}")")

## matrix_ok MATRIX REGION_BEDS... : header, unique ids, id set == ids in the region BEDs, numeric and >= 0, some signal in every column
matrix_ok() {
  local m="$1"; shift
  awk -F'\t' -v hdr="${HDR}" 'FNR == NR { want[$4] = 1; nw = length(want); next }
    FNR == 1 { if ($0 != hdr) { print "bad header"; exit } ; nc = NF; next }
    { if ($1 in got) { print "duplicate id " $1; exit }; got[$1] = 1
      for (i = 2; i <= nc; i++) { if ($i !~ /^[0-9]+(\.[0-9]+)?$/) { print "non-numeric/negative value in row " $1; exit }; if ($i > 0) pos[i]++ } }
    END { if (length(got) != nw) { print "rows " length(got) " != region ids " nw; exit }
          for (i = 2; i <= nc; i++) if (!pos[i]) { print "column " i " has no signal"; exit }
          for (g in want) if (!(g in got)) { print "missing id " g; exit }
          print "ok" }' <(cat "$@") "${m}" | grep -q '^ok$' && echo ok
}
chk "$(matrix_ok "${GENE_M}"   "${REG_DIR}/genes_plus.bed" "${REG_DIR}/genes_minus.bed")"     "gene matrix: header, unique ids, ids == gene BED, numeric >= 0, signal in every sample"
chk "$(matrix_ok "${EXON_M}"   "${REG_DIR}/exons_plus.bed" "${REG_DIR}/exons_minus.bed")"     "exon matrix: same checks"
chk "$(matrix_ok "${INTRON_M}" "${REG_DIR}/introns_plus.bed" "${REG_DIR}/introns_minus.bed")" "intron matrix: same checks"
chk "$(matrix_ok "${DOG_M}"    "${REG_DIR}/dog_plus.bed" "${REG_DIR}/dog_minus.bed")"         "DoG matrix: same checks"

## signal = RPKM x kb. Exons and introns of a gene are disjoint and lie inside the gene body, so
## exon signal + intron signal <= gene signal (every gene, every sample)
inv=$(awk -F'\t' 'FILENAME ~ /\.bed$/ { split(FILENAME, p, "/"); k = p[length(p)]; sub(/_(plus|minus)\.bed$/, "", k); L[k, $4] += $3 - $2; next }
  FNR == 1 { f++; nc = NF; next }
  { r[f, $1] = $0 }
  END { bad = 0; n = 0
    for (key in r) { split(key, a, SUBSEP); if (a[1] != 1) continue; g = a[2]
      split(r[1, g], G, "\t"); split(r[2, g], E, "\t"); split(r[3, g], I, "\t")
      for (i = 2; i <= nc; i++) { n++
        sg = G[i] * L["genes", g] / 1000; se = E[i] * L["exons", g] / 1000
        si = 0; if (I[i] != "") si = I[i] * L["introns", g] / 1000
        if (se > sg * (1 + 1e-4) + 1e-6 || se + si > sg * (1 + 1e-4) + 1e-6) bad++ } }
    print (n > 0 && bad == 0) ? "ok" : "bad " bad "/" n }' \
  "${REG_DIR}/genes_plus.bed" "${REG_DIR}/genes_minus.bed" "${REG_DIR}/exons_plus.bed" "${REG_DIR}/exons_minus.bed" \
  "${REG_DIR}/introns_plus.bed" "${REG_DIR}/introns_minus.bed" "${GENE_M}" "${EXON_M}" "${INTRON_M}")
chk "$([[ "${inv}" == ok ]] && echo ok)" "exon signal + intron signal <= gene signal for every gene and sample (${inv})"
echo "========================================================"
[[ ${NFAIL} -eq 0 ]] || { echo "*** ${NFAIL} self-check(s) failed" >&2; exit 1; }
