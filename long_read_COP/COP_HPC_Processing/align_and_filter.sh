#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cop_align
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Long-read alignment -> primary-only BAM. Aligns a long-read
# FASTQ (minimap2, splice-aware ONT preset), restricts to
# uniquely-mapping reads, then filters to primary alignments only
# (-F 0x900: excludes secondary + supplementary). Output
# (<NAME>_primary_only.bam, indexed) is the expected input for
# opt_scripts/RREA.sh.
#
# Usage: sbatch align_and_filter.sh <NAME> <fastq> <reference.fa>
##################################################################
set -euo pipefail

NAME=$1
FASTQ=$2
REF=$3
THREADS=${SLURM_CPUS_PER_TASK:-16}

module load minimap2
module load samtools

# Align (splice-aware, ONT preset)
minimap2 -t ${THREADS} -ax splice -uf -k14 "${REF}" "${FASTQ}" > "${NAME}_minimap2.sam"

# Restrict to uniquely-mapping reads (drop any read name appearing
# more than once among mapped alignments)
grep '^@' "${NAME}_minimap2.sam" > "${NAME}_headers.sam"
grep -v '^@' "${NAME}_minimap2.sam" | sort > "${NAME}_alignment.sam"
awk '$3!="*" {print}' "${NAME}_alignment.sam" | cut -f1 | sort | uniq -c | \
  awk '$1==1 {print $2}' > "${NAME}_uniq_names.txt"
grep -F -f "${NAME}_uniq_names.txt" "${NAME}_alignment.sam" > "${NAME}_uniq_alignment.sam"
cat "${NAME}_headers.sam" "${NAME}_uniq_alignment.sam" > "${NAME}_uniq.sam"

# SAM -> sorted BAM
samtools view -bT "${REF}" "${NAME}_uniq.sam" | samtools sort -@ ${THREADS} -o "${NAME}_uniq_sort.bam"

# Primary alignments only (exclude secondary 0x100 + supplementary 0x800)
samtools view -b -F 0x900 "${NAME}_uniq_sort.bam" > "${NAME}_primary_only.bam"
samtools index "${NAME}_primary_only.bam"

rm -f "${NAME}_minimap2.sam" "${NAME}_headers.sam" "${NAME}_alignment.sam" \
      "${NAME}_uniq_names.txt" "${NAME}_uniq_alignment.sam" "${NAME}_uniq.sam" \
      "${NAME}_uniq_sort.bam"

echo "Done: ${NAME}_primary_only.bam"
