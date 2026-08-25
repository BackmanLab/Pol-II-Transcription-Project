
#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Submission wrapper for merge_bam_bpm_bigwig.sh
# Submits a single Quest job to merge RIP-seq BAM replicates
# and generate BPM bigWigs for all four conditions.
#
# Inputs (hardcoded in merge_bam_bpm_bigwig.sh):
#   POLR2A: WT_P2AB_Rep2_S19_L001_filt.sort.bam
#           WT_P2AB_Rep3_S20_L001_filt.sort.bam
#   EZH2:   WT_EZH2_Rep1_S9_L001_filt.sort.bam
#           WT_EZH2_Rep2_S10_L001_filt.sort.bam
#   IgG:    WT_IgG_S14_L001_filt.sort.bam
#   Input:  WT_Input_S13_L001_filt.sort.bam
#
# Outputs written to:
#   /projects/b1042/BackmanLab/Lucas/022025_RIPseq/coverage_bpm/
#
# Usage: sh submit_merge_bpm.sh
##################################################################

SCRIPT="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/BAM/merge_bam_bpm_bigwig.sh"
LOG_DIR="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/logs"

mkdir -p "${LOG_DIR}"

## Verify script exists before submitting
if [ ! -f "${SCRIPT}" ]; then
  echo "ERROR: script not found at ${SCRIPT}"
  exit 1
fi

echo "Submitting: ${SCRIPT}"
echo ""
echo "Input BAMs:"
echo "  POLR2A: WT_P2AB_Rep2_S19_L001_filt.sort.bam"
echo "          WT_P2AB_Rep3_S20_L001_filt.sort.bam"
echo "  EZH2:   WT_EZH2_Rep1_S9_L001_filt.sort.bam"
echo "          WT_EZH2_Rep2_S10_L001_filt.sort.bam"
echo "  IgG:    WT_IgG_S14_L001_filt.sort.bam"
echo "  Input:  WT_Input_S13_L001_filt.sort.bam"
echo ""

sbatch \
  -A b1042 \
  -p genomics \
  -t 18:00:00 \
  -n 10 \
  --mem=150G \
  --job-name=rip_merge_bpm \
  --output="${LOG_DIR}/merge_bpm_out.txt" \
  --error="${LOG_DIR}/merge_bpm_err.txt" \
  --mail-type=END,FAIL \
  --wrap="sh ${SCRIPT}"

echo "Job submitted. Monitor with:"
echo "  tail -f ${LOG_DIR}/merge_bpm_out.txt"
