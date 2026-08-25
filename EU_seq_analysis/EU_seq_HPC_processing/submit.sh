#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Submission wrapper for EU-seq alignment pipeline.
# Submits one job per sample via sbatch --wrap. Run from any
# directory on Quest.
# Usage: sh submit_EUseq.sh
##################################################################

SCRIPT="/projects/b1042/BackmanLab/Lucas/090124_euseq/align_euseq.sh"
FASTQ_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/fastq"
GENOME="/projects/b1042/BackmanLab/Lucas/090124_euseq/genome/ref_hisat/hg38"
GTF="/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
LOG_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/logs"
THREADS=10

mkdir -p "${LOG_DIR}"

## R1 files only — R2 derived automatically in the alignment script
array=(
  BM20_S37_R1_001.fastq.gz      ## POLR1A 8hrs aux Rep1
  BM21_S38_R1_001.fastq.gz      ## POLR1A 8hrs aux Rep2
  BM22_S39_R1_001.fastq.gz      ## POLR1A 8hrs aux Rep3
  BM23_S40_R1_001.fastq.gz      ## POLR2A 8hrs aux Rep1
  BM24_S41_R1_001.fastq.gz      ## POLR2A 8hrs aux Rep2
  BM25_S42_R1_001.fastq.gz      ## POLR2A 8hrs aux Rep3
  BM27_S44_R1_001.fastq.gz      ## ActD 6hrs Rep2 (BM26 excluded)
  BM28_S45_R1_001.fastq.gz      ## ActD 6hrs Rep3
  BM29_S46_R1_001.fastq.gz      ## WT Rep1
  BM30_S47_R1_001.fastq.gz      ## WT Rep2
  BM31_S48_R1_001.fastq.gz      ## WT Rep3
  No_EU_ctrl_S49_R1_001.fastq.gz ## No EU control
)

for file in "${array[@]}"; do
  echo "Submitting: ${file}"
  sbatch \
    -A b1042 \
    -p genomics-himem \
    -t 12:00:00 \
    -n ${THREADS} \
    --mem=64G \
    --job-name=euseq_${file} \
    --output="${LOG_DIR}/${file}_out.txt" \
    --error="${LOG_DIR}/${file}_err.txt" \
    --mail-type=END,FAIL \
    --wrap="sh ${SCRIPT} -f ${file} -w ${FASTQ_DIR} -g ${GENOME} -r ${GTF} -t ${THREADS}"
  sleep 1
done

echo ""
echo "All EU-seq jobs submitted."
echo "Monitor with: squeue -u lmc0633"
echo "After all jobs complete, run: sh merge_EUseq_bams.sh"
