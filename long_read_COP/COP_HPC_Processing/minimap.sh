#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=minimap
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Single-sample long-read alignment (splice-aware).
# Bare minimap2 alignment only (no read filtering). See
# align_and_filter.sh for the full alignment -> unique-reads ->
# primary-only-BAM pipeline expected by opt_scripts/RREA.sh.
##################################################################

module load minimap2
minimap2 -ax splice -uf -k14 Homo_sapiens.GRCh38.dna.primary_assembly.fa SRR8268942.fastq > K562_4sU_ONT.sam

