#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Generates bigWig coverage files for EU-seq and total RNA-seq.
#
# EU-seq outputs per condition:
#   1. Scaled bigWig        — EU scaled by total RNA-seq SF
#   2. Subtracted bigWig    — EU - NoEU, floor at 0
#   3. Log2 ratio bigWig    — log2(EU/NoEU) with pseudocount
#
# Total RNA-seq outputs per condition:
#   1. CPM bigWig           — standard CPM normalization
#   2. Scaled bigWig        — scaled by same total RNA-seq SF
#                             for direct co-plotting with EU-seq
#
# Scale factors derived from total RNA-seq mapped read counts
# (min_reads / condition_reads), anchored to WT (smallest library).
#
# No EU control (WT cells) scaled by WT scale factor before
# subtraction — consistent with ChIP-seq IgG handling.
#
# Usage:
#   sbatch -A b1042 -p genomics-himem -t 12:00:00 -n 16 --mem=300G \
#     --wrap="sh generate_bigwigs.sh"
##################################################################

set -euo pipefail

##--------------------------------------------------------------## 1. Scale factors
## Derived from total RNA-seq mapped read counts
## (min_reads / condition_reads), WT is reference (SF = 1.0)

#SF_POLR1A="0.552524"
SF_POLR2A="0.747437"
SF_ActD="0.887481"
SF_WT="1.000000"
SF_NoEU="${SF_WT}"    ## No EU ctrl from WT cells, scaled by WT factor

##--------------------------------------------------------------## 2. Paths

EU_BAM_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/BAM/merged"
RNA_BAM_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/BAM/merged"

EU_COV_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/coverage"
RNA_COV_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/coverage"
LOG_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/logs"

mkdir -p "${EU_COV_DIR}/scaled"
mkdir -p "${EU_COV_DIR}/subtracted"
mkdir -p "${EU_COV_DIR}/log2ratio"
mkdir -p "${RNA_COV_DIR}/cpm"
mkdir -p "${RNA_COV_DIR}/scaled"
mkdir -p "${LOG_DIR}"

THREADS=10
BINSIZE=20
SMOOTH=60

NO_EU_BAM="${EU_BAM_DIR}/No_EU_ctrl_merged.sort.bam"
[[ -f "${NO_EU_BAM}" ]] || {
  echo "*** error: No EU BAM not found: ${NO_EU_BAM} ***"
  exit 1
}

module purge
module load deeptools/3.5.6

##--------------------------------------------------------------## Helper functions

## EU-seq: scaled by total RNA-seq SF
make_scaled_bw() {
  local BAM=$1
  local OUT=$2
  local SF=$3
  echo "$(date): Scaled bigWig (SF=${SF}) — $(basename ${BAM})"
  bamCoverage \
    --bam                "${BAM}" \
    --scaleFactor        "${SF}" \
    --binSize            "${BINSIZE}" \
    --smoothLength       "${SMOOTH}" \
    --outFileName        "${OUT}" \
    --numberOfProcessors "${THREADS}" \
    --ignoreDuplicates \
    --extendReads
  echo "  Done: ${OUT}"
}

## EU-seq: subtracted (EU - NoEU), floor at 0
make_subtracted_bw() {
  local EU_BAM=$1
  local OUT=$2
  local SF_EU=$3
  echo "$(date): Subtracted bigWig (EU - NoEU, floor=0) — $(basename ${EU_BAM})"
  bamCompare \
    --bamfile1           "${EU_BAM}" \
    --bamfile2           "${NO_EU_BAM}" \
    --scaleFactors       "${SF_EU}:${SF_NoEU}" \
    --operation          subtract \
    --binSize            "${BINSIZE}" \
    --smoothLength       "${SMOOTH}" \
    --outFileName        "${OUT}" \
    --numberOfProcessors "${THREADS}" \
    --ignoreDuplicates \
    --extendReads \
    --pseudocount        0
  echo "  Done: ${OUT}"
}

## EU-seq: log2 ratio (log2(EU/NoEU))
make_log2_bw() {
  local EU_BAM=$1
  local OUT=$2
  local SF_EU=$3
  echo "$(date): Log2 ratio bigWig — $(basename ${EU_BAM})"
  bamCompare \
    --bamfile1           "${EU_BAM}" \
    --bamfile2           "${NO_EU_BAM}" \
    --scaleFactors       "${SF_EU}:${SF_NoEU}" \
    --operation          log2 \
    --binSize            "${BINSIZE}" \
    --smoothLength       "${SMOOTH}" \
    --outFileName        "${OUT}" \
    --numberOfProcessors "${THREADS}" \
    --ignoreDuplicates \
    --extendReads \
    --pseudocount        1
  echo "  Done: ${OUT}"
}

## Total RNA-seq: CPM normalization
make_cpm_bw() {
  local BAM=$1
  local OUT=$2
  echo "$(date): CPM bigWig — $(basename ${BAM})"
  bamCoverage \
    --bam                "${BAM}" \
    --normalizeUsing     CPM \
    --binSize            "${BINSIZE}" \
    --smoothLength       "${SMOOTH}" \
    --outFileName        "${OUT}" \
    --numberOfProcessors "${THREADS}" \
    --ignoreDuplicates
  echo "  Done: ${OUT}"
}

## Total RNA-seq: scaled by same SF as EU-seq (for co-plotting)
make_rna_scaled_bw() {
  local BAM=$1
  local OUT=$2
  local SF=$3
  echo "$(date): Scaled bigWig (SF=${SF}) — $(basename ${BAM})"
  bamCoverage \
    --bam                "${BAM}" \
    --scaleFactor        "${SF}" \
    --binSize            "${BINSIZE}" \
    --smoothLength       "${SMOOTH}" \
    --outFileName        "${OUT}" \
    --numberOfProcessors "${THREADS}" \
    --ignoreDuplicates
  echo "  Done: ${OUT}"
}

##--------------------------------------------------------------## 3. EU-seq bigWigs

echo "========================================================"
echo "$(date): Generating EU-seq bigWigs"
echo "========================================================"
echo ""

declare -A EU_BAMS
#EU_BAMS["POLR1A_8hrs_aux"]="${EU_BAM_DIR}/POLR1A_8hrs_aux_merged.sort.bam"
EU_BAMS["POLR2A_8hrs_aux"]="${EU_BAM_DIR}/POLR2A_8hrs_aux_merged.sort.bam"
EU_BAMS["ActD_6hrs"]="${EU_BAM_DIR}/ActD_6hrs_merged.sort.bam"
EU_BAMS["WT"]="${EU_BAM_DIR}/WT_merged.sort.bam"

declare -A SF_MAP
#SF_MAP["POLR1A_8hrs_aux"]="${SF_POLR1A}"
SF_MAP["POLR2A_8hrs_aux"]="${SF_POLR2A}"
SF_MAP["ActD_6hrs"]="${SF_ActD}"
SF_MAP["WT"]="${SF_WT}"

for COND in POLR2A_8hrs_aux ActD_6hrs WT; do

  BAM="${EU_BAMS[$COND]}"
  SF="${SF_MAP[$COND]}"

  [[ -f "${BAM}" ]] || {
    echo "*** error: BAM not found: ${BAM} ***"
    exit 1
  }

  echo "--- EU-seq: ${COND} (SF = ${SF}) ---"

  make_scaled_bw \
    "${BAM}" \
    "${EU_COV_DIR}/scaled/${COND}.scaled.bigWig" \
    "${SF}"

  make_subtracted_bw \
    "${BAM}" \
    "${EU_COV_DIR}/subtracted/${COND}.subtracted.bigWig" \
    "${SF}"

  make_log2_bw \
    "${BAM}" \
    "${EU_COV_DIR}/log2ratio/${COND}.log2ratio.bigWig" \
    "${SF}"

  echo ""
done

##--------------------------------------------------------------## 4. Total RNA-seq bigWigs

echo "========================================================"
echo "$(date): Generating total RNA-seq bigWigs"
echo "========================================================"
echo ""

declare -A RNA_BAMS
#RNA_BAMS["POLR1A_8hrs_aux"]="${RNA_BAM_DIR}/POLR1A_8hrs_aux_merged.sort.bam"
RNA_BAMS["POLR2A_8hrs_aux"]="${RNA_BAM_DIR}/POLR2A_8hrs_aux_merged.sort.bam"
RNA_BAMS["ActD_6hrs"]="${RNA_BAM_DIR}/ActD_6hrs_merged.sort.bam"
RNA_BAMS["WT"]="${RNA_BAM_DIR}/WT_merged.sort.bam"

for COND in POLR2A_8hrs_aux ActD_6hrs WT; do

  BAM="${RNA_BAMS[$COND]}"
  SF="${SF_MAP[$COND]}"

  [[ -f "${BAM}" ]] || {
    echo "*** error: BAM not found: ${BAM} ***"
    exit 1
  }

  echo "--- Total RNA-seq: ${COND} ---"

  ## CPM — for general visualization and signal shape
  make_cpm_bw \
    "${BAM}" \
    "${RNA_COV_DIR}/cpm/${COND}.cpm.bigWig"

  ## Scaled — same SF as EU-seq, for direct co-plotting
  make_rna_scaled_bw \
    "${BAM}" \
    "${RNA_COV_DIR}/scaled/${COND}.scaled.bigWig" \
    "${SF}"

  echo ""
done

##--------------------------------------------------------------## 5. Summary

echo "========================================================"
echo "$(date): All bigWigs complete."
echo ""
echo "EU-seq scaled bigWigs:          ${EU_COV_DIR}/scaled/"
echo "EU-seq subtracted bigWigs:      ${EU_COV_DIR}/subtracted/"
echo "EU-seq log2 ratio bigWigs:      ${EU_COV_DIR}/log2ratio/"
echo "Total RNA-seq CPM bigWigs:      ${RNA_COV_DIR}/cpm/"
echo "Total RNA-seq scaled bigWigs:   ${RNA_COV_DIR}/scaled/"
echo ""
echo "For metagene co-plotting use:"
echo "  EU-seq:       ${EU_COV_DIR}/subtracted/"
echo "  Total RNA-seq: ${RNA_COV_DIR}/scaled/"
echo "========================================================"
