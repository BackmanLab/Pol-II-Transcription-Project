#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=180G
#SBATCH -t 16:00:00
#SBATCH --job-name=bwa_index
#SBATCH -o logs/bwa_index.out
#SBATCH -e logs/bwa_index.err
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Builds the BWA index for hg38, required by Juicer's
# alignment step.

module load bwa/0.7.17
JUICEDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer

bwa index ${JUICEDIR}/references/hg38.fa
echo "BWA index complete"
