#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --job-name=site_positions
#SBATCH -o logs/site_positions.out
#SBATCH -e logs/site_positions.err
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Generates DpnII (GATC) restriction cut site positions
# for hg38, required by Juicer's alignment/dedup pipeline.

module load python
JUICEDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer

# Generate DpnII (GATC) cut site positions for hg38
python ${JUICEDIR}/scripts/generate_site_positions.py \
    DpnII hg38 ${JUICEDIR}/references/hg38.fa

# Move output to correct location
mv hg38_DpnII.txt ${JUICEDIR}/restriction_sites/
ln -sf ${JUICEDIR}/restriction_sites/hg38_DpnII.txt ${JUICEDIR}/references/

# Derive chrom.sizes from the site positions file
awk 'BEGIN{OFS="\t"}{print $1, $NF}' \
    ${JUICEDIR}/restriction_sites/hg38_DpnII.txt > \
    ${JUICEDIR}/references/hg38.chrom.sizes

echo "Site positions and chrom.sizes generated"
