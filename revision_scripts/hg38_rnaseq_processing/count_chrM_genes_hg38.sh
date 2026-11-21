#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH -t 12:00:00
#SBATCH --job-name=chrM_genes_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/analysis/mito_qc/logs/count_chrM_genes_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/analysis/mito_qc/logs/count_chrM_genes_%j.err
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Per-gene / per-class chrM read counts for the mito scaling QC
# (SI figure; mito_qc.R). Counts the SAME read set that
# the mito scale factors were computed from:
#   EU-seq        (paired-end): read1, mapped, MAPQ>0, primary
#                               samtools view -q 1 -F 2308 -f 64
#   total RNA-seq (single-end): mapped, MAPQ>0, primary
#                               samtools view -q 1 -F 2308
# so that sum(all classes) == chrM_reads in
# scale_factors_mito_samples.txt / scale_factors_mito_total_samples.txt
# (checked at the end).
#
# chrM genes come from the same Ensembl GRCh38.115 GTF used for the
# RPKM script (gene rows on chrM / MT / M; biotype -> class):
#   Mt_rRNA -> rRNA | protein_coding -> mRNA | Mt_tRNA -> tRNA |
#   reads overlapping no chrM gene -> Unannotated
# Each read is assigned to ONE gene (the gene with the largest overlap;
# ties -> first in the BED), so the classes add up to the chrM total
# (ATP8/ATP6 and ND4L/ND4 overlap). Unstranded assignment; the sense
# column records whether the read is on the gene's strand.
#
# Outputs (OUT_DIR = BASE/analysis/mito_qc):
#   chrM_gene_counts.tsv  : assay, sample, condition, gene_id, gene_name,
#                           biotype, class, gene_len, reads, reads_sense
#                           (unannotated reads: gene_id = "Unannotated")
#   chrM_bin_counts.tsv   : assay, sample, condition, bin_start, bin_end,
#                           reads   (read midpoints, BIN bp bins)
#   chrM_genes.bed        : the chrM gene table used
#
# Usage:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/analysis/mito_qc/logs
#   sbatch count_chrM_genes_hg38.sh
##################################################################

set -o pipefail

module purge
module load samtools/1.14
module load bedtools/2.29.2

BASE="${BASE:-/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq}"
EU_BAM_DIR="${EU_BAM_DIR:-${BASE}/euseq/BAM}"
TOT_BAM_DIR="${TOT_BAM_DIR:-${BASE}/total/BAM}"
GTF="${GTF:-${BASE}/reference/Homo_sapiens.GRCh38.115.chr.gtf}"
MITO="${MITO:-chrM}"
BIN="${BIN:-100}"
OUT_DIR="${OUT_DIR:-${BASE}/analysis/mito_qc}"
EU_SAMPLES_FILE="${EU_SAMPLES_FILE:-${BASE}/scale_factors_mito_samples.txt}"
TOT_SAMPLES_FILE="${TOT_SAMPLES_FILE:-${BASE}/scale_factors_mito_total_samples.txt}"

die() { echo "*** error: $*" >&2; exit 1; }

mkdir -p "${OUT_DIR}" || die "cannot create ${OUT_DIR}"
GENES_BED="${OUT_DIR}/chrM_genes.bed"
OUT_GENES="${OUT_DIR}/chrM_gene_counts.tsv"
OUT_BINS="${OUT_DIR}/chrM_bin_counts.tsv"
TMP="$(mktemp -d "${OUT_DIR}/.tmp.XXXXXX")" || die "cannot create temp dir"
trap 'rm -rf "${TMP}"' EXIT

## assay | sample | condition | BAM | extra samtools flags
SAMPLES=(
  "EU-seq|WT_1|WT|${EU_BAM_DIR}/BM29_S46_filt.sort.bam|-f 64"
  "EU-seq|WT_2|WT|${EU_BAM_DIR}/BM30_S47_filt.sort.bam|-f 64"
  "EU-seq|WT_3|WT|${EU_BAM_DIR}/BM31_S48_filt.sort.bam|-f 64"
  "EU-seq|POLR2A_1|POLR2A_8hrs_aux|${EU_BAM_DIR}/BM23_S40_filt.sort.bam|-f 64"
  "EU-seq|POLR2A_2|POLR2A_8hrs_aux|${EU_BAM_DIR}/BM24_S41_filt.sort.bam|-f 64"
  "EU-seq|POLR2A_3|POLR2A_8hrs_aux|${EU_BAM_DIR}/BM25_S42_filt.sort.bam|-f 64"
  "EU-seq|NoEU|NoEU|${EU_BAM_DIR}/No_EU_ctrl_S49_filt.sort.bam|-f 64"
  "Total RNA-seq|WT_1|WT|${TOT_BAM_DIR}/WT_HCT116_CTRL_Rep1_S12_filt.sort.bam|"
  "Total RNA-seq|WT_2|WT|${TOT_BAM_DIR}/WT_HCT116_CTRL_Rep2_S13_filt.sort.bam|"
  "Total RNA-seq|WT_3|WT|${TOT_BAM_DIR}/WT_HCT116_CTRL_Rep3_S14_filt.sort.bam|"
  "Total RNA-seq|POLR2A_1|POLR2A_8hrs_aux|${TOT_BAM_DIR}/Pol2_8hrs_Aux_Rep1_S5_filt.sort.bam|"
  "Total RNA-seq|POLR2A_2|POLR2A_8hrs_aux|${TOT_BAM_DIR}/Pol2_8hrs_Aux_Rep2_S6_filt.sort.bam|"
  "Total RNA-seq|POLR2A_3|POLR2A_8hrs_aux|${TOT_BAM_DIR}/Pol2_8hrs_Aux_Rep3_S7_filt.sort.bam|"
)

[[ -s "${GTF}" ]] || die "GTF not found: ${GTF}"

##--------------------------------------------------------------## 1. chrM gene table
echo "========================================================"
echo "$(date): Building chrM gene table from ${GTF}"
echo "========================================================"

## gene rows on chrM / MT / M (Ensembl names it MT; the hg38 BAMs use ${MITO});
## BED name = gene_id|gene_name|gene_biotype
grep -v "^#" "${GTF}" | \
awk -v mito="${MITO}" 'BEGIN {FS = "\t"; OFS = "\t"}
  ($1 == "chrM" || $1 == "MT" || $1 == "M" || $1 == mito) && $3 == "gene" {
    id = ""; nm = ""; bt = ""
    if (match($9, /gene_id "[^"]+"/))      id = substr($9, RSTART + 9, RLENGTH - 10)
    if (match($9, /gene_name "[^"]+"/))    nm = substr($9, RSTART + 11, RLENGTH - 12)
    if (match($9, /gene_biotype "[^"]+"/)) bt = substr($9, RSTART + 14, RLENGTH - 15)
    if (id == "") next
    if (nm == "") nm = id
    if (bt == "") bt = "NA"
    print mito, $4 - 1, $5, id "|" nm "|" bt, 0, $7
  }' | LC_ALL=C sort -k2,2n > "${GENES_BED}" || die "chrM gene table failed"

NG=$(wc -l < "${GENES_BED}")
[[ "${NG}" -gt 0 ]] || die "no gene rows on chrM/MT/M in ${GTF} (check the chromosome names and that the GTF has 'gene' rows)"
echo "  chrM genes: ${NG}"
awk -F'\t' '{split($4, a, "|"); n[a[3]]++} END {for (b in n) printf "    %-20s %d\n", b, n[b]}' "${GENES_BED}" | sort

printf "assay\tsample\tcondition\tgene_id\tgene_name\tbiotype\tclass\tgene_len\treads\treads_sense\n" > "${OUT_GENES}"
printf "assay\tsample\tcondition\tbin_start\tbin_end\treads\n" > "${OUT_BINS}"

##--------------------------------------------------------------## 2. Per-sample counting
echo ""
echo "========================================================"
echo "$(date): Assigning chrM reads to genes"
echo "========================================================"

for entry in "${SAMPLES[@]}"; do
  IFS='|' read -r assay sample cond bam flags <<< "${entry}"
  [[ -f "${bam}" ]]     || die "BAM not found: ${bam}"
  [[ -f "${bam}.bai" ]] || [[ -f "${bam%.bam}.bai" ]] || die "BAM index not found for ${bam}"
  samtools view -H "${bam}" | grep -q -P "^@SQ\tSN:${MITO}\t" || die "contig ${MITO} not in header of ${bam}"

  ## chrM reads, same filters as the scale factors -> BED6 (one row per read)
  # shellcheck disable=SC2086
  samtools view -b -q 1 -F 2308 ${flags} "${bam}" "${MITO}" | \
    bedtools bamtobed -i stdin > "${TMP}/reads.bed" || die "chrM read extraction failed for ${bam}"
  NR_READS=$(wc -l < "${TMP}/reads.bed")

  ## read -> gene with the largest overlap (-wao keeps non-overlapping reads, gene cols ".")
  LC_ALL=C sort -k1,1 -k2,2n -k3,3n "${TMP}/reads.bed" | \
    bedtools intersect -a stdin -b "${GENES_BED}" -wao | \
  awk -v assay="${assay}" -v sample="${sample}" -v cond="${cond}" -v bin="${BIN}" \
      -v genes="${TMP}/gene_out.tsv" -v bins="${TMP}/bin_out.tsv" '
    BEGIN {FS = OFS = "\t"}
    function flush() {
      if (rn == "") return
      if (best == "") { cnt["Unannotated"]++ }
      else { cnt[best]++; if (bsense) sen[best]++ }
      mid = int((rs + re) / 2 / bin)
      bc[mid]++
    }
    {
      # cols: 1-6 read (chr,start,end,name,score,strand), 7-12 gene, 13 overlap
      if ($4 != rn) {
        flush()
        rn = $4; rs = $2; re = $3; best = ""; bov = 0; bsense = 0
      }
      if ($7 != "." && $13 + 0 > bov) { bov = $13 + 0; best = $10; bsense = ($6 == $12) }
    }
    END {
      flush()
      for (g in cnt) print g, cnt[g], sen[g] + 0 > genes
      for (b in bc) print b * bin, (b + 1) * bin, bc[b] > bins
    }'
  [[ -f "${TMP}/gene_out.tsv" ]] || : > "${TMP}/gene_out.tsv"
  [[ -f "${TMP}/bin_out.tsv" ]]  || : > "${TMP}/bin_out.tsv"

  ## full gene table per sample (zeros included); class from the biotype
  awk -F'\t' -v assay="${assay}" -v sample="${sample}" -v cond="${cond}" \
      -v counts="${TMP}/gene_out.tsv" 'BEGIN {
        OFS = "\t"
        while ((getline line < counts) > 0) { split(line, a, "\t"); c[a[1]] = a[2]; s[a[1]] = a[3] }
      }
      function cls(bt) {
        if (bt == "Mt_rRNA") return "rRNA"
        if (bt == "Mt_tRNA") return "tRNA"
        if (bt == "protein_coding") return "mRNA"
        return "Other"
      }
      {
        split($4, f, "|"); key = $4
        print assay, sample, cond, f[1], f[2], f[3], cls(f[3]), $3 - $2, c[key] + 0, s[key] + 0
      }
      END { print assay, sample, cond, "Unannotated", "Unannotated", "NA", "Unannotated", 0, c["Unannotated"] + 0, 0 }' \
      "${GENES_BED}" >> "${OUT_GENES}"

  sort -k1,1n "${TMP}/bin_out.tsv" | awk -F'\t' -v OFS="\t" -v a="${assay}" -v s="${sample}" -v c="${cond}" \
      '{print a, s, c, $1, $2, $3}' >> "${OUT_BINS}"
  rm -f "${TMP}/gene_out.tsv" "${TMP}/bin_out.tsv" "${TMP}/reads.bed"

  echo "  ${assay} ${sample} (${cond}): ${NR_READS} chrM reads"
done

##--------------------------------------------------------------## 3. Self-checks
echo ""
echo "========================================================"
echo "Self-checks"
NFAIL=0
ok()   { echo "  PASS  $*"; }
bad()  { echo "  FAIL  $*"; NFAIL=$((NFAIL + 1)); }

## (a) classes add up to the chrM reads of the scale-factor scripts
check_totals() {
  local assay="$1" file="$2" label="$3"
  if [[ ! -s "${file}" ]]; then echo "  SKIP  ${label}: ${file} not found"; return; fi
  awk -F'\t' -v assay="${assay}" -v label="${label}" '
    FNR == NR { if (FNR > 1) ref[$1] = $4; next }
    FNR > 1 && $1 == assay { tot[$2] += $9 }
    END {
      bad = 0
      for (s in ref) {
        if (!(s in tot)) { printf "  FAIL  %s %s: no gene counts\n", label, s; bad++; continue }
        if (tot[s] != ref[s]) { printf "  FAIL  %s %s: classes sum to %d, scale-factor file has %d\n", label, s, tot[s], ref[s]; bad++ }
        else printf "  PASS  %s %s: classes sum to chrM_reads (%d)\n", label, s, tot[s]
      }
      exit bad
    }' "${file}" "${OUT_GENES}" || NFAIL=$((NFAIL + $?))
}
check_totals "EU-seq"        "${EU_SAMPLES_FILE}"  "EU-seq"
check_totals "Total RNA-seq" "${TOT_SAMPLES_FILE}" "Total RNA-seq"

## (b) bin counts add up to the same totals
awk -F'\t' 'FNR == NR { if (FNR > 1) g[$1 "|" $2] += $9; next }
            FNR > 1 { b[$1 "|" $2] += $6 }
            END { bad = 0; for (k in g) if (g[k] != b[k] + 0) { printf "  FAIL  bin counts differ from gene counts for %s (%d vs %d)\n", k, b[k], g[k]; bad++ }
                  if (!bad) print "  PASS  bin counts == gene counts for every sample"; exit bad }' \
    "${OUT_GENES}" "${OUT_BINS}" || NFAIL=$((NFAIL + 1))

## (c) every sample has all genes + the Unannotated row
EXPECT=$((NG + 1))
awk -F'\t' -v e="${EXPECT}" 'NR > 1 { n[$1 "|" $2]++ } END { bad = 0; for (k in n) if (n[k] != e) { printf "  FAIL  %s has %d rows, expected %d\n", k, n[k], e; bad++ }
                                      if (!bad) print "  PASS  every sample has " e " rows (genes + Unannotated)"; exit bad }' "${OUT_GENES}" || NFAIL=$((NFAIL + 1))

echo ""
echo "Per-class chrM reads (all samples):"
awk -F'\t' 'NR > 1 { k = $1 "\t" $2 "\t" $7; r[k] += $9; t[$1 "\t" $2] += $9; keys[k] = 1 }
  END { for (k in keys) { split(k, a, "\t"); printf "  %-14s %-10s %-13s %10d  (%.1f%%)\n", a[1], a[2], a[3], r[k], 100 * r[k] / t[a[1] "\t" a[2]] } }' "${OUT_GENES}" | sort

echo ""
echo "Written: ${OUT_GENES}"
echo "Written: ${OUT_BINS}"
echo "Written: ${GENES_BED}"
if [[ "${NFAIL}" -gt 0 ]]; then
  echo "*** ${NFAIL} self-check(s) FAILED"
  exit 1
fi
echo "$(date): done, all self-checks passed"
