#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH --mem=250G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=rsem_build_ref
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/090124_euseq/genome/rsem/logs/build_ref_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/090124_euseq/genome/rsem/logs/build_ref_%j.err
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Builds RSEM transcriptome reference from hg38 genome fasta
# and Ensembl GTF. Run once before rsem-calculate-expression.
#
# Output: indexed RSEM reference in REF_DIR under prefix "hg38"
#
# Usage:
#   sbatch build_rsem_ref.sh
##################################################################

set -euo pipefail

##--------------------------------------------------------------## Paths

GENOME_FA="/projects/b1042/BackmanLab/Lucas/090124_euseq/genome/ref_hisat/hg38.fa"
GTF="/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
REF_DIR="/projects/b1042/BackmanLab/Lucas/090124_euseq/genome/rsem"
LOG_DIR="${REF_DIR}/logs"
REF_PREFIX="${REF_DIR}/hg38"
THREADS=12

mkdir -p "${REF_DIR}" "${LOG_DIR}"

echo "========================================================"
echo "$(date): Building RSEM reference"
echo "  Genome FA  : ${GENOME_FA}"
echo "  GTF        : ${GTF}"
echo "  Output dir : ${REF_DIR}"
echo "  Prefix     : ${REF_PREFIX}"
echo "========================================================"

[[ -f "${GENOME_FA}" ]] || { echo "*** error: genome fasta not found: ${GENOME_FA} ***"; exit 1; }
[[ -f "${GTF}"       ]] || { echo "*** error: GTF not found: ${GTF} ***";               exit 1; }

##--------------------------------------------------------------## Build reference

module purge
module load rsem

rsem-prepare-reference \
  --gtf "${GTF}" \
  --num-threads "${THREADS}" \
  "${GENOME_FA}" \
  "${REF_PREFIX}"

echo ""
echo "========================================================"
echo "$(date): RSEM reference build complete"
echo "  Reference files in : ${REF_DIR}"
ls -lh "${REF_DIR}/"
echo "========================================================"
