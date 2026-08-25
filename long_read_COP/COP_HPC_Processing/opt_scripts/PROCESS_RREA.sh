#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Submission wrapper for PROCESS_RREA.R -- loops
# over the RREA.R output list and submits PROCESS_RREA.R for each.
#
# Usage: Rscript PROCESS_RREA.R --input <rrea_output.txt> --output <file>
##################################################################

cd /projects/b1042/BackmanLab/Lucas/eric_splicing/out

array=(GM19102_chr_rep1a_rrea_output.txt  K562_4sU_ONT_rrea_output.txt    K562_DMSO_rrea_output.txt       K562_PlaB_rep2_rrea_output.txt  K562_dRNA_rep2_rrea_output.txt
K562_4sU_ONT_rep2_rrea_output.txt  K562_DMSO_rep2_rrea_output.txt  K562_PlaB_rep1_rrea_output.txt  K562_dRNA_rep1_rrea_output.txt)

for file in "${!array[@]}";
do echo "Submitting job for: ${array[file]}" >&2;

ROOT="/projects/b1042/BackmanLab/Lucas/eric_splicing"
SCRIPT="/scripts/opt_scripts/PROCESS_RREA.R"

OUT=$(echo ${array[file]} | sed 's/_output.txt//g')
LOG=$(echo ${array[file]} | sed 's/_output.txt//g')

sbatch -A p32171 -p short -t 0-04:00 -n 1 \
--job-name process_rrea_${array[file]} \
--mem=100G \
--output=out_${LOG}.log \
--error=err_${LOG}.err \
--wrap="module load R/4.3.0; Rscript ${ROOT}${SCRIPT} --input ${ROOT}/out/${array[file]} --output ${ROOT}/out/${OUT}"

sleep 1	    # wait 1 second between each job submission

done
