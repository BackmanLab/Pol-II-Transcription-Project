#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Submission wrapper for total RNA-seq alignment
# pipeline. Single-end reads. Submits one job per sample via
# sbatch --wrap.
# Usage: sh submit_totalRNAseq.sh
##################################################################

SCRIPT="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/align_totalRNAseq.sh"
FASTQ_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/fastq"
GENOME="/projects/b1042/BackmanLab/Lucas/090124_euseq/genome/ref_hisat/hg38"
GTF="/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
LOG_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/logs"
THREADS=10

mkdir -p "${LOG_DIR}"

array=(
  Pol1_8hrs_Aux_Rep1_S1_R1_001.fastq.gz    ## POLR1A 8hrs aux Rep1
  Pol1_8hrs_Aux_Rep2_S2_R1_001.fastq.gz    ## POLR1A 8hrs aux Rep2
  Pol1_8hrs_Aux_Rep3_S3_R1_001.fastq.gz    ## POLR1A 8hrs aux Rep3
  Pol1_8hrs_Aux_Rep4_S4_R1_001.fastq.gz    ## POLR1A 8hrs aux Rep4
  Pol2_8hrs_Aux_Rep1_S5_R1_001.fastq.gz    ## POLR2A 8hrs aux Rep1
  Pol2_8hrs_Aux_Rep2_S6_R1_001.fastq.gz    ## POLR2A 8hrs aux Rep2
  Pol2_8hrs_Aux_Rep3_S7_R1_001.fastq.gz    ## POLR2A 8hrs aux Rep3
  Pol2_8hrs_Aux_Rep4_S8_R1_001.fastq.gz    ## POLR2A 8hrs aux Rep4
  ActD_6hrs_Rep1_S9_R1_001.fastq.gz        ## ActD 6hrs Rep1
  ActD_6hrs_Rep2_S10_R1_001.fastq.gz       ## ActD 6hrs Rep2
  ActD_6hrs_Rep3_S11_R1_001.fastq.gz       ## ActD 6hrs Rep3
  WT_HCT116_CTRL_Rep1_S12_R1_001.fastq.gz  ## WT Rep1
  WT_HCT116_CTRL_Rep2_S13_R1_001.fastq.gz  ## WT Rep2
  WT_HCT116_CTRL_Rep3_S14_R1_001.fastq.gz  ## WT Rep3
)

for file in "${array[@]}"; do
  echo "Submitting: ${file}"
  sbatch \
    -A b1042 \
    -p genomics \
    -t 12:00:00 \
    -n ${THREADS} \
    --mem=64G \
    --job-name=rnaseq_${file} \
    --output="${LOG_DIR}/${file}_out.txt" \
    --error="${LOG_DIR}/${file}_err.txt" \
    --mail-type=END,FAIL \
    --wrap="sh ${SCRIPT} -f ${file} -w ${FASTQ_DIR} -g ${GENOME} -r ${GTF} -t ${THREADS}"
  sleep 1
done

echo ""
echo "All total RNA-seq jobs submitted."
echo "Monitor with: squeue -u lmc0633"
echo "After all jobs complete, run: sh merge_totalRNAseq_bams.sh"
