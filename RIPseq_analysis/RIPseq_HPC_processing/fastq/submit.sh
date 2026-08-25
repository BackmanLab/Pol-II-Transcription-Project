#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Submission wrapper for align.sh. Loops over the
# full RIP-seq sample list (POLR2A/EZH2 IP, IgG, Input, across WT
# and P1KO backgrounds) and submits align.sh for each R1 fastq.
##################################################################

array=(
    # WT POLR2A (P1AB) replicates
    WT_P1AB_Rep1_S4_L001_R1_001.fastq.gz
    WT_P1AB_Rep2_S5_L001_R1_001.fastq.gz
    WT_P1AB_rRNA_S16_L001_R1_001.fastq.gz

    # WT POLR2A (P2AB) replicates
    WT_P2AB_Rep2_S19_L001_R1_001.fastq.gz
    WT_P2AB_Rep3_S20_L001_R1_001.fastq.gz
    WT_P2AB_rRNA_S17_L001_R1_001.fastq.gz

    # WT EZH2
    WT_EZH2_Rep1_S9_L001_R1_001.fastq.gz
    WT_EZH2_Rep2_S10_L001_R1_001.fastq.gz

    # WT controls
    WT_IgG_S14_L001_R1_001.fastq.gz
    WT_Input_S13_L001_R1_001.fastq.gz

    # P1KO POLR2A (P1AB) replicates
    P1KO_P1AB_Rep1_S1_L001_R1_001.fastq.gz
    P1KO_P1AB_Rep2_S2_L001_R1_001.fastq.gz
    P1KO_P1AB_Rep3_S3_L001_R1_001.fastq.gz
    P1KO_P1AB_rRNA_S15_L001_R1_001.fastq.gz

    # P1KO controls
    P1KO_Input_S11_L001_R1_001.fastq.gz

    # P2KO POLR2A (P1AB) replicates
    P2KO_P1AB_Rep1_S6_L001_R1_001.fastq.gz
    P2KO_P1AB_Rep2_S7_L001_R1_001.fastq.gz
    P2KO_P1AB_Rep3_S8_L001_R1_001.fastq.gz
    P2KO_P1AB_rRNA_S18_L001_R1_001.fastq.gz

    # P2KO controls
    P2KO_input_S12_L001_R1_001.fastq.gz
)

SCRIPT="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/align.sh"
FASTQ="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/fastq"
GENOME="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/genome/ref_hisat/hg38"
GTF="/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
THREADS=10

for file in "${!array[@]}"; do

    echo "Submitting: ${array[file]}"

    sbatch \
        -A p32171 \
        -p genhimem \
        -t 24:00:00 \
        -n ${THREADS} \
        --mem=250G \
        --job-name=rip_${array[file]} \
        --output=${FASTQ}/../logs/${array[file]}_out.txt \
        --error=${FASTQ}/../logs/${array[file]}_err.txt \
        --wrap="sh ${SCRIPT} -f ${array[file]} -w ${FASTQ} -g ${GENOME} -r ${GTF} -t ${THREADS}"

    sleep 1

done
