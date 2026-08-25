#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=40G
#SBATCH -t 2880
#SBATCH --array=1-6
#SBATCH --job-name=count_reads
#SBATCH -o logs/count_%a.out
#SBATCH -e logs/count_%a.err
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Count reads in each R1 FASTQ (SLURM array, 6 tasks).
# After all 6 tasks finish, run: bash 00b_make_samplesheet.sh
#
# Submit: sbatch 00_count_reads.sh

mkdir -p logs counts

FASTQ_DIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp

# Map array index to condition and replicate (R1 only — R2 count is identical)
SAMPLES=(
    "WT_HCT116_CTRL 1"
    "WT_HCT116_CTRL 2"
    "Pol2_6hrs_Aux 1"
    "Pol2_6hrs_Aux 2"
    "1hr_ActD 1"
    "1hr_ActD 2"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
COND=$(echo "${SAMPLES[$IDX]}" | awk '{print $1}')
REP=$(echo  "${SAMPLES[$IDX]}" | awk '{print $2}')

R1=$(ls ${FASTQ_DIR}/${COND}/*R1*Rep${REP}*.fastq.gz)
R2=$(ls ${FASTQ_DIR}/${COND}/*R2*Rep${REP}*.fastq.gz)

echo "Task ${SLURM_ARRAY_TASK_ID}: ${COND} Rep${REP}"
echo "  R1: $R1"

# Count reads: wc -l is faster than awk; divide by 4 for FASTQ format
N=$(( $(zcat "$R1" | wc -l) / 4 ))

echo "  Reads: $N"

# Write result to a per-task file (avoids write collisions between array tasks)
echo -e "${COND}\t${REP}\t${R1}\t${R2}\t${N}" \
    > counts/count_${COND}_Rep${REP}.txt

echo "Done: ${COND} Rep${REP} = ${N} reads"
