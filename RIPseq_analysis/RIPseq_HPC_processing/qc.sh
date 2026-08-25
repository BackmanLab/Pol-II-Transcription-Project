#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=rip_q

##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Basic RIP-seq QC at the BAM level: mapping rate and duplication
# (whatever samtools flagstat reports — see note below).
#
# PCR Bottleneck Coefficient (PBC1) is deliberately NOT included.
# It's a ChIP-seq metric that assumes independent DNA fragments
# landing at the same start position by chance is rare — for
# RNA-based RIP/CLIP-seq, real reads legitimately stack at exons,
# splice junctions, and TSSs, so positional duplication doesn't
# distinguish PCR artifact from real signal the way it does for
# ChIP-seq. The CLIP-seq field's own standard (eCLIP/ENCODE,
# nf-core/clipseq) uses UMI-based deduplication instead of
# coordinate-based methods for exactly this reason.
#
# Peak-calling-based FRiP is DELIBERATELY NOT included here.
# Already tried (Piranha, MACS2) and found uninformative for this
# dataset: WT_Input reflects transcription from all three RNA
# polymerases, while P2AB is specifically Pol2-associated RNA —
# there's no sparse "peaks vs flat background" structure the way
# there would be for a site-specific DNA-binding factor, since most
# genes are Pol2-transcribed regardless. The R companion script
# (ripseq_qc_analysis.R) computes genomic feature distribution
# (Promoter/Enhancer/Exon/Intron/TES fraction of signal) as the
# more appropriate substitute for this specific assay.
##################################################################

set -euo pipefail

BAM_DIR="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/BAM"
GTF_PATH="/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
OUT_DIR="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/ripseq_pipeline/qc"
mkdir -p "${OUT_DIR}"

## Real filenames on disk (verified against actual directory listing,
## not guessed) — samtools/bedtools naming convention SAMPLE=filename
## does not match this project's simplified sample names.
declare -A BAM_FILES=(
  [P2AB_Rep2]="WT_P2AB_Rep2_S19_L001_filt.sort.bam"
  [P2AB_Rep3]="WT_P2AB_Rep3_S20_L001_filt.sort.bam"
  [EZH2_Rep1]="WT_EZH2_Rep1_S9_L001_filt.sort.bam"
  [EZH2_Rep2]="WT_EZH2_Rep2_S10_L001_filt.sort.bam"
  [WT_Input]="WT_Input_S13_L001_filt.sort.bam"
  [P2KO_input]="P2KO_input_S12_L001_filt.sort.bam"
  [WT_IgG]="WT_IgG_S14_L001_filt.sort.bam"
)
SAMPLES=(P2AB_Rep2 P2AB_Rep3 EZH2_Rep1 EZH2_Rep2 WT_Input P2KO_input WT_IgG)

module load samtools

##--------------------------------------------------------------## 1. Per-sample QC metrics
SUMMARY="${OUT_DIR}/qc_summary.tsv"
echo -e "sample\ttotal_reads\tmapped_reads\tmapping_pct\tflagstat_duplicates" > "${SUMMARY}"

for sample in "${SAMPLES[@]}"; do
  bam_filename="${BAM_FILES[$sample]:-}"
  if [[ -z "${bam_filename}" ]]; then
    echo "  WARNING: no filename mapping defined for ${sample}, skipping."
    continue
  fi
  bam="${BAM_DIR}/${bam_filename}"
  if [[ ! -f "${bam}" ]]; then
    echo "  WARNING: BAM not found for ${sample}, skipping: ${bam}"
    continue
  fi
  echo "Processing ${sample} (${bam_filename})..."

  ## --- Mapping / duplication ---
  ## Reuses the flagstat file already on disk (${bam}.flagstat.txt)
  ## instead of re-running samtools flagstat.
  flagstat_out="${bam}.flagstat.txt"
  if [[ ! -f "${flagstat_out}" ]]; then
    echo "  No existing flagstat file found (${flagstat_out}), running samtools flagstat..."
    flagstat_out="${OUT_DIR}/${sample}_flagstat.txt"
    samtools flagstat "${bam}" > "${flagstat_out}"
  fi

  ## Every grep here needs a "|| echo" fallback: with set -euo pipefail
  ## active, grep returning exit 1 (no match — not itself an error
  ## message, just a silent nonzero exit) kills the ENTIRE script
  ## instantly via pipefail, with nothing printed. This is almost
  ## certainly what "died silently" was — three of these four lines
  ## were unprotected.
  total_reads=$(grep "in total" "${flagstat_out}" | awk '{print $1}' || echo "NA")
  mapped_reads=$(grep " mapped (" "${flagstat_out}" | head -1 | awk '{print $1}' || echo "NA")
  mapping_pct=$(grep " mapped (" "${flagstat_out}" | head -1 | grep -oP '\(\K[0-9.]+(?=%)' || echo "NA")
  dup_reads=$(grep "duplicates" "${flagstat_out}" | awk '{print $1}' || echo "NA")

  if [[ "${total_reads}" == "NA" || "${mapped_reads}" == "NA" ]]; then
    echo "  WARNING: could not parse flagstat output for ${sample} —"
    echo "           samtools flagstat wording may differ from what this"
    echo "           script expects ('in total', ' mapped ('). Raw output:"
    sed 's/^/    /' "${flagstat_out}"
  fi

  if [[ "${dup_reads}" == "0" ]]; then
    echo "  NOTE: ${sample} shows 0 duplicates in flagstat — this BAM was"
    echo "        likely never dup-marked upstream (e.g. via samtools markdup"
    echo "        or Picard MarkDuplicates). This does NOT mean duplication is"
    echo "        actually zero — flagstat only reports what's already flagged"
    echo "        in the BAM, it doesn't detect duplicates itself."
  fi

  ## PCR Bottleneck Coefficient (PBC1) deliberately NOT computed here.
  ## It's a ChIP-seq metric that assumes independent DNA fragments
  ## landing at the same start position by chance is rare — for
  ## RNA-based RIP/CLIP-seq, real reads legitimately stack at exons,
  ## splice junctions, and TSSs, so positional duplication doesn't
  ## distinguish PCR artifact from real signal the way it does for
  ## ChIP-seq. This is exactly why the CLIP-seq field (eCLIP/ENCODE,
  ## nf-core/clipseq) uses UMI-based deduplication instead of
  ## coordinate-based methods. If this library's prep included UMIs,
  ## that's the metric to compute instead — ask before assuming either
  ## way.

  echo -e "${sample}\t${total_reads}\t${mapped_reads}\t${mapping_pct}\t${dup_reads}" >> "${SUMMARY}"
  echo "  Mapping: ${mapping_pct}%"
done

echo ""
echo "QC summary written to: ${SUMMARY}"
