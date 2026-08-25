#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Submission wrapper for RREA.R -- loops over the
# primary-only BAM sample list and submits RREA.R for each.
#
# Usage:
  # Rscript RREA.R --bam <file> --gtf <file> --out <file> [options]

cd /projects/b1042/BackmanLab/Lucas/eric_splicing/bam

array=(GM19102_chr_rep1a_primary_only.bam  K562_4sU_ONT_rep2_primary_only.bam  K562_DMSO_rep2_primary_only.bam  K562_PlaB_rep2_primary_only.bam  K562_dRNA_rep2_primary_only.bam
K562_4sU_ONT_primary_only.bam K562_DMSO_primary_only.bam K562_PlaB_rep1_primary_only.bam  K562_dRNA_rep1_primary_only.bam)

for file in "${!array[@]}";
do echo "Submitting job for: ${array[file]}" >&2;

ROOT="/projects/b1042/BackmanLab/Lucas/eric_splicing"
GTF="/gtf/gencode.v47.basic.annotation.gtf"
SCRIPT="/scripts/opt_scripts/RREA.R"

OUT=$(echo ${array[file]} | sed 's/_primary_only.bam/_rrea_output.txt/g')
LOG=$(echo ${array[file]} | sed 's/_primary_only.bam//g')

sbatch -A p32171 -p short -t 0-04:00 -n 8 \
--job-name run_rrea_${array[file]} \
--mem=100G \
--output=out_${LOG}.log \
--error=err_${LOG}.err \
--wrap="module load R/4.3.0; Rscript ${ROOT}${SCRIPT} --bam ${ROOT}/bam/${array[file]} --gtf ${ROOT}${GTF} --out ${ROOT}/out/${OUT} --workers 8"

sleep 1	    # wait 1 second between each job submission

done

