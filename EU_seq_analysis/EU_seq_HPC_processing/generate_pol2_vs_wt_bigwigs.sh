#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Generates log2(POLR2A / WT) bigWigs for both EU-seq and total
# RNA-seq, showing fold-change in signal relative to wildtype.
#
# These differ from the EU/NoEU log2ratio files — those correct
# for technical background. These files show biological change:
#   > 0: more signal in POLR2A-degraded than WT
#   < 0: less signal in POLR2A-degraded than WT
#
# Scale factors applied before ratio (from total RNA-seq):
#   POLR2A: 0.747437
#   WT:     1.000000
#
# Outputs:
#   EU-seq:        POLR2A_vs_WT.log2ratio.bigWig
#   Total RNA-seq: POLR2A_vs_WT.log2ratio.bigWig
#
# Usage:
#   sbatch -A p32171 -p genhimem -t 02:00:00 -n 10 --mem=32G \
#     --output=/projects/b1042/BackmanLab/Lucas/090124_euseq/logs/pol2_vs_wt_out.txt \
#     --error=/projects/b1042/BackmanLab/Lucas/090124_euseq/logs/pol2_vs_wt_err.txt \
#     --wrap="sh generate_pol2_vs_wt_bigwigs.sh"
##################################################################

set -euo pipefail

##--------------------------------------------------------------## 1. Scale factors

SF_POLR2A="0.747437"
SF_WT="1.000000"

##--------------------------------------------------------------## 2. Paths

EU_BAM_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/BAM/merged"
RNA_BAM_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/BAM/merged"

EU_COV_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/coverage/log2ratio"
RNA_COV_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/coverage/log2ratio"
LOG_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/logs"

mkdir -p "${EU_COV_DIR}" "${RNA_COV_DIR}" "${LOG_DIR}"

THREADS=10
BINSIZE=20
SMOOTH=60

## BAMs
EU_POL2_BAM="${EU_BAM_DIR}/POLR2A_8hrs_aux_merged.sort.bam"
EU_WT_BAM="${EU_BAM_DIR}/WT_merged.sort.bam"

RNA_POL2_BAM="${RNA_BAM_DIR}/POLR2A_8hrs_aux_merged.sort.bam"
RNA_WT_BAM="${RNA_BAM_DIR}/WT_merged.sort.bam"

## Validate
for BAM in "${EU_POL2_BAM}" "${EU_WT_BAM}" "${RNA_POL2_BAM}" "${RNA_WT_BAM}"; do
  [[ -f "${BAM}" ]] || { echo "*** error: BAM not found: ${BAM} ***"; exit 1; }
done

module purge
module load deeptools/3.5.6

##--------------------------------------------------------------## 3. EU-seq log2(POLR2A / WT)

echo "========================================================"
echo "$(date): EU-seq log2(POLR2A / WT)"
echo "========================================================"

bamCompare \
  --bamfile1           "${EU_POL2_BAM}" \
  --bamfile2           "${EU_WT_BAM}" \
  --scaleFactors       "${SF_POLR2A}:${SF_WT}" \
  --operation          log2 \
  --binSize            "${BINSIZE}" \
  --smoothLength       "${SMOOTH}" \
  --outFileName        "${EU_COV_DIR}/POLR2A_vs_WT.log2ratio.bigWig" \
  --numberOfProcessors "${THREADS}" \
  --ignoreDuplicates \
  --extendReads \
  --pseudocount        1

echo "  Done: ${EU_COV_DIR}/POLR2A_vs_WT.log2ratio.bigWig"

##--------------------------------------------------------------## 4. Total RNA-seq log2(POLR2A / WT)

echo ""
echo "========================================================"
echo "$(date): Total RNA-seq log2(POLR2A / WT)"
echo "========================================================"

bamCompare \
  --bamfile1           "${RNA_POL2_BAM}" \
  --bamfile2           "${RNA_WT_BAM}" \
  --scaleFactors       "${SF_POLR2A}:${SF_WT}" \
  --operation          log2 \
  --binSize            "${BINSIZE}" \
  --smoothLength       "${SMOOTH}" \
  --outFileName        "${RNA_COV_DIR}/POLR2A_vs_WT.log2ratio.bigWig" \
  --numberOfProcessors "${THREADS}" \
  --ignoreDuplicates \
  --pseudocount        1

echo "  Done: ${RNA_COV_DIR}/POLR2A_vs_WT.log2ratio.bigWig"

##--------------------------------------------------------------## 5. Summary

echo ""
echo "========================================================"
echo "$(date): Complete."
echo ""
echo "EU-seq log2(POLR2A/WT):       ${EU_COV_DIR}/POLR2A_vs_WT.log2ratio.bigWig"
echo "Total RNA-seq log2(POLR2A/WT): ${RNA_COV_DIR}/POLR2A_vs_WT.log2ratio.bigWig"
echo ""
echo "Values > 0: more signal in POLR2A-degraded than WT"
echo "Values < 0: less signal in POLR2A-degraded than WT"
echo "========================================================"
