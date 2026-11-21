#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH -t 12:00:00
#SBATCH --job-name=hisat2_index_pol2
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101823_rnaseq/reference/logs/hisat2_build_out.txt
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101823_rnaseq/reference/logs/hisat2_build_err.txt
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Downloads the hg38 reference and builds the HISAT2 splice-aware
# index for the Pol II total RNA-seq re-alignment.
#   FASTA : UCSC hg38.fa.gz (chr-prefixed names; md5 checked against
#           UCSC md5sum.txt)
#   GTF   : Ensembl release 115, Homo_sapiens.GRCh38.115.gtf.gz
# Existing downloads are reused. If the compute node has no internet,
# run the four wget/md5 lines of section 1 on a login node first.
#
# Ensembl names chromosomes 1..22, X, Y, MT, but the FASTA uses chr1..
# chrX, chrM, so splice sites extracted from the raw GTF would never
# match the genome. The GTF is first rewritten with chr-prefixed names
# (scaffolds / patches not in the FASTA are dropped); splice sites and
# exons are extracted from that copy.
#
# hisat2_extract_splice_sites.py / hisat2_extract_exons.py are not in
# Quest's hisat2/2.1.0 module; they are taken from the HISAT2 2.2.1
# source tree (plain Python, no compilation).
#
# Outputs (REFDIR = /projects/b1042/BackmanLab/Lucas/101823_rnaseq/reference):
#   hg38.fa, Homo_sapiens.GRCh38.115.gtf, Homo_sapiens.GRCh38.115.chr.gtf
#   hisat/hg38.{ss,exon}, hisat/hg38.[1-8].ht2
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101823_rnaseq/reference/logs
# Usage:
#   sbatch build_hisat2_index_pol2.sh
# Next: sbatch align_totalRNAseq_pol2.sh
##################################################################

set -o pipefail

module purge
module load hisat2/2.1.0

REFDIR="${REFDIR:-/projects/b1042/BackmanLab/Lucas/101823_rnaseq/reference}"
FASTA="${FASTA:-${REFDIR}/hg38.fa}"
GTF="${GTF:-${REFDIR}/Homo_sapiens.GRCh38.115.gtf}"
UCSC_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips"
GTF_URL="https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"
HISAT2_PKG="${HISAT2_PKG:-/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer/references/hisat2-2.2.1}"

OUTDIR="${REFDIR}/hisat"
OUTPREFIX="${OUTDIR}/hg38"
GTF_CHR="${REFDIR}/Homo_sapiens.GRCh38.115.chr.gtf"

mkdir -p "${OUTDIR}" "${REFDIR}/logs"

die() { echo "*** error: $*" >&2; exit 1; }

FASTA_MD5_OK="no"
for p in hisat2_extract_splice_sites.py hisat2_extract_exons.py; do
  [[ -x "${HISAT2_PKG}/${p}" ]] || die "${HISAT2_PKG}/${p} not found or not executable (chmod +x it)"
done

##--------------------------------------------------------------## 1. Download the reference (skipped if already present)
if [[ ! -s "${FASTA}" ]]; then
  echo "$(date): Downloading hg38.fa.gz and md5sum.txt from UCSC"
  wget -c --tries=5 -q -O "${FASTA}.gz" "${UCSC_URL}/hg38.fa.gz"      || die "FASTA download failed (run the wget on a login node)"
  wget --tries=5 -q -O "${REFDIR}/ucsc_md5sum.txt" "${UCSC_URL}/md5sum.txt" || die "md5sum.txt download failed"
  want=$(awk '$2 == "hg38.fa.gz" {print $1}' "${REFDIR}/ucsc_md5sum.txt")
  have=$(md5sum "${FASTA}.gz" | awk '{print $1}')
  [[ -n "${want}" && "${want}" == "${have}" ]] || die "hg38.fa.gz md5 mismatch (expected ${want}, got ${have})"
  FASTA_MD5_OK="yes"
  gunzip -c "${FASTA}.gz" > "${FASTA}" || die "FASTA unzip failed"
fi
if [[ ! -s "${GTF}" ]]; then
  echo "$(date): Downloading ${GTF_URL}"
  wget -c --tries=5 -q -O "${GTF}.gz" "${GTF_URL}" || die "GTF download failed (run the wget on a login node)"
  gzip -t "${GTF}.gz"                              || die "GTF download is not a valid gzip file"
  gunzip -c "${GTF}.gz" > "${GTF}"                 || die "GTF unzip failed"
fi
[[ -s "${FASTA}" ]] || die "hg38 FASTA not found: ${FASTA}"
[[ -s "${GTF}" ]]   || die "GTF not found: ${GTF}"

##--------------------------------------------------------------## 2. chr-prefixed GTF restricted to FASTA contigs
echo "$(date): FASTA contigs"
grep '^>' "${FASTA}" | sed 's/^>//; s/[ \t].*//' > "${REFDIR}/hg38_contigs.txt" || die "contig listing failed"
echo "  contigs: $(wc -l < "${REFDIR}/hg38_contigs.txt")"

echo "$(date): Writing chr-prefixed GTF"
grep -v '^#' "${GTF}" | awk 'BEGIN {FS = OFS = "\t"}
  FNR == NR { ok[$1] = 1; next }
  { c = $1
    if (c ~ /^([0-9]+|X|Y)$/) c = "chr" c
    else if (c == "MT") c = "chrM"
    else { dropped++; next }
    if (!(c in ok)) { missing[c] = 1; dropped++; next }
    $1 = c; print }
  END { for (m in missing) print "  contig absent from FASTA: " m > "/dev/stderr"
        print "  GTF lines dropped (scaffolds / not in FASTA): " dropped + 0 > "/dev/stderr" }' \
  "${REFDIR}/hg38_contigs.txt" - > "${GTF_CHR}" || die "GTF rewrite failed"
[[ -s "${GTF_CHR}" ]] || die "chr-prefixed GTF is empty"

##--------------------------------------------------------------## 3. Splice sites, exons, index
echo "$(date): Extracting splice sites and exons"
"${HISAT2_PKG}/hisat2_extract_splice_sites.py" "${GTF_CHR}" > "${OUTDIR}/hg38.ss"   || die "splice site extraction failed"
"${HISAT2_PKG}/hisat2_extract_exons.py"        "${GTF_CHR}" > "${OUTDIR}/hg38.exon" || die "exon extraction failed"

echo "$(date): Building HISAT2 index"
hisat2-build -p 12 --ss "${OUTDIR}/hg38.ss" --exon "${OUTDIR}/hg38.exon" "${FASTA}" "${OUTPREFIX}" || die "hisat2-build failed"

##--------------------------------------------------------------## 4. Self-checks on the live outputs
echo "========================================================"
echo "Self-checks"
NFAIL=0
chk() { if [[ "$1" == ok ]]; then echo "  CHECK PASS  $2"; else echo "  CHECK FAIL  $2"; NFAIL=$((NFAIL + 1)); fi; }

chk "$([[ $(cut -f1 "${GTF_CHR}" | sort -u | grep -c -E '^chr([0-9]+|X|Y|M)$') -ge 24 ]] && echo ok)" "chr-prefixed GTF covers chr1-22, X, Y, M"
chk "$([[ -z $(cut -f1 "${GTF_CHR}" | sort -u | grep -v -x -F -f "${REFDIR}/hg38_contigs.txt") ]] && echo ok)" "every GTF contig exists in the FASTA"
[[ "${FASTA_MD5_OK}" == "yes" ]] && chk ok "hg38.fa.gz md5 matches UCSC md5sum.txt"
chk "$([[ $(grep -c -E '^chr([0-9]+|X|Y|M)$' "${REFDIR}/hg38_contigs.txt") -eq 25 ]] && echo ok)" "FASTA has chr1-22, X, Y, M"
chk "$(grep -q 'gene_id' "${GTF}" && grep -q $'\texon\t' "${GTF}" && echo ok)" "GTF has exon rows with gene_id"
chk "$([[ $(wc -l < "${OUTDIR}/hg38.ss") -gt 100000 ]] && echo ok)" "splice-site file has > 100,000 sites ($(wc -l < "${OUTDIR}/hg38.ss"))"
chk "$([[ $(wc -l < "${OUTDIR}/hg38.exon") -gt 100000 ]] && echo ok)" "exon file has > 100,000 exons ($(wc -l < "${OUTDIR}/hg38.exon"))"
chk "$([[ $(ls "${OUTPREFIX}".*.ht2 2>/dev/null | wc -l) -eq 8 ]] && echo ok)" "index has 8 .ht2 files"
echo "========================================================"
[[ ${NFAIL} -eq 0 ]] || { echo "*** ${NFAIL} self-check(s) failed" >&2; exit 1; }

echo "$(date): HISAT2 index complete: ${OUTPREFIX}.*.ht2"
