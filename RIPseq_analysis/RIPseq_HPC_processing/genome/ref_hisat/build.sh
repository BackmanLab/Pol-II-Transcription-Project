#!/bin/bash
#SBATCH -A p32171 ##--## SLURM start
#SBATCH -p short
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=150G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --time 4:00:00
#SBATCH --job-name=hisat_build
#SBATCH --output=build_out
#SBATCH --error=build_err ##--## SLURM end
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Builds the HISAT2 splice-aware genome index (hg38)
# used for RIP-seq alignment.

module load hisat2/2.1.0

ref=/projects/b1042/BackmanLab/Lucas/022025_RIPseq/genome/ref_hisat/hg38.fa
exon=/projects/b1042/BackmanLab/Lucas/022025_RIPseq/genome/ref_hisat/genome.exon
ss=/projects/b1042/BackmanLab/Lucas/022025_RIPseq/genome/ref_hisat/genome.ss
cd /projects/b1042/BackmanLab/Lucas/022025_RIPseq/genome/ref_hisat/

hisat2-build -p 12 --exon $exon --ss $ss $ref hg38
