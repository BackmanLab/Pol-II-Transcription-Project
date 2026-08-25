#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# This script performs QC, trimming, mapping, and coverage
# generation for RIP-seq data (POLR2A, EZH2, Input) in HCT116
# cells. Processes one paired-end sample per call. Designed to
# be submitted via sbatch --wrap on Northwestern Quest HPC.
#
# Usage (via submission wrapper):
#   sbatch -A p32171 -p genhimem -t 1-24:00 -n 10 --mem=300G \
#     --wrap="sh align_HISAT_RIPseq.sh \
#       -f SAMPLE_R1_001.fastq.gz \
#       -w /projects/b1042/BackmanLab/Lucas/022025_RIPseq/fastq \
#       -g /projects/b1042/BackmanLab/Lucas/022025_RIPseq/genome/ref_hisat/hg38 \
#       -r /projects/b1042/BackmanLab/Lucas/022025_RIPseq/genome/gtf/ensembl_gtf/Homo_sapiens.GRCh38.113.gtf \
#       -t 10"
#
# Output directories are created as siblings of the fastq dir:
#   BAM/  SAM/  OUT/  coverage/  fastqc/
##################################################################

set -euo pipefail

##--------------------------------------------------------------## Help

helpFunction()
{
   echo ""
   echo "************** RIP-seq alignment pipeline (Quest) **************"
   echo ""
   echo "Usage: sh $0 -f <R1 fastq> -w <fastq dir> -g <HISAT2 index> -r <GTF> -t <threads>"
   echo ""
   echo "  -f   Forward read filename (R1), e.g. SAMPLE_R1_001.fastq.gz"
   echo "  -w   Full path to directory containing FASTQs"
   echo "  -g   Full path to HISAT2 genome index prefix"
   echo "  -r   Full path to primary GTF annotation file"
   echo "  -t   Number of threads"
   echo ""
   echo "****************************************************************"
   exit 1
}

##--------------------------------------------------------------## Parse arguments

while getopts "f:w:g:r:t:" opt
do
   case $opt in
      f) R1=$OPTARG ;;
      w) wdir=$OPTARG ;;
      g) gendir=$OPTARG ;;
      r) gtfdir=$OPTARG ;;
      t) threads=$OPTARG ;;
      ?) echo "Unknown argument --${OPTARG}"; helpFunction ;;
   esac
done

[[ -z "${R1:-}"      ]] && { echo "*** error: forward read (-f) required ***";      helpFunction; }
[[ -z "${wdir:-}"    ]] && { echo "*** error: fastq directory (-w) required ***";   helpFunction; }
[[ -z "${gendir:-}"  ]] && { echo "*** error: HISAT2 index (-g) required ***";      helpFunction; }
[[ -z "${gtfdir:-}"  ]] && { echo "*** error: GTF file (-r) required ***";          helpFunction; }
[[ -z "${threads:-}" ]] && { echo "*** error: thread count (-t) required ***";      helpFunction; }

[[ -d "${wdir}" ]]       || { echo "*** error: fastq directory not found: ${wdir} ***"; exit 1; }
[[ -f "${wdir}/${R1}" ]] || { echo "*** error: R1 not found: ${wdir}/${R1} ***";        exit 1; }

##--------------------------------------------------------------## Derive R2

R2=$(echo "${R1}" | sed 's/R1_001/R2_001/g')
[[ -f "${wdir}/${R2}" ]] || { echo "*** error: R2 not found: ${wdir}/${R2} ***"; exit 1; }

##--------------------------------------------------------------## Directory structure
## All output dirs are siblings of the fastq directory (one level up from wdir).
## Example: wdir = /projects/.../022025_RIPseq/fastq/
##          root = /projects/.../022025_RIPseq/
##          BAM  = /projects/.../022025_RIPseq/BAM/

root="$(dirname "${wdir}")/"

outfold="${root}OUT";      mkdir -p "${outfold}"
samfold="${root}SAM";      mkdir -p "${samfold}"
qcfold="${root}fastqc";    mkdir -p "${qcfold}"
covfold="${root}coverage"; mkdir -p "${covfold}"
bamfold="${root}BAM";      mkdir -p "${bamfold}"

##--------------------------------------------------------------## Derive file names

SAMPLE=$(echo "${R1}" | sed 's/_R1_001.fastq.gz//g')

T1=$(echo "${R1}" | sed 's/R1_001.fastq.gz/R1_001_val_1.fq/g')
T2=$(echo "${R2}" | sed 's/R2_001.fastq.gz/R2_001_val_2.fq/g')

S1="${SAMPLE}_R1_001.sam"

B_FILT="${SAMPLE}_filt.bam"
B_PSORT="${SAMPLE}_filt.sort.bam"

BW_CPM="${SAMPLE}.cpm.bigWig"
BW_RPGC="${SAMPLE}.rpgc.bigWig"

HISAT_LOG="${samfold}/${SAMPLE}_hisat2.log"

echo "========================================================"
echo "$(date): RIP-seq pipeline starting"
echo "  Sample   : ${SAMPLE}"
echo "  Root     : ${root}"
echo "  Threads  : ${threads}"
echo "========================================================"

##--------------------------------------------------------------## 1. FastQC — raw reads

echo ""
echo "--- 1. FastQC (raw reads) ---"

module purge
module load fastqc/0.12.0

fastqc -t "${threads}" \
  "${wdir}/${R1}" \
  "${wdir}/${R2}" \
  -o "${qcfold}"

echo "FastQC raw reads complete."

##--------------------------------------------------------------## 2. Trim with TrimGalore

echo ""
echo "--- 2. TrimGalore ---"

module purge
module load TrimGalore/0.6.10

trim_galore \
  --paired \
  --retain_unpaired \
  --dont_gzip \
  --cores "${threads}" \
  -o "${outfold}" \
  "${wdir}/${R1}" \
  "${wdir}/${R2}"

[[ -f "${outfold}/${T1}" ]] || { echo "*** error: trimmed R1 not found: ${outfold}/${T1} ***"; exit 1; }
[[ -f "${outfold}/${T2}" ]] || { echo "*** error: trimmed R2 not found: ${outfold}/${T2} ***"; exit 1; }

echo "TrimGalore complete."

##--------------------------------------------------------------## 3. Align with HISAT2

echo ""
echo "--- 3. HISAT2 alignment ---"

module purge
module load hisat2/2.1.0

hisat2 \
  -p "${threads}" \
  --dta \
  -x "${gendir}" \
  -1 "${outfold}/${T1}" \
  -2 "${outfold}/${T2}" \
  -S "${samfold}/${S1}" \
  2> "${HISAT_LOG}"

echo "HISAT2 complete. Alignment summary:"
cat "${HISAT_LOG}"

##--------------------------------------------------------------## 4. Filter, sort, index with samtools

echo ""
echo "--- 4. samtools: filter, sort, index ---"

module purge
module load samtools/1.14

## -F 3852 removes: unmapped (4), mate unmapped (8), not primary (256),
##   fails QC (512), PCR duplicate (1024), supplementary (2048)
## -q 30: retain only uniquely mapping reads (MAPQ >= 30 in HISAT2)
samtools view \
  -bS \
  -F 3852 \
  -q 30 \
  "${samfold}/${S1}" \
  -o "${bamfold}/${B_FILT}"

samtools sort \
  -@ "${threads}" \
  "${bamfold}/${B_FILT}" \
  -o "${bamfold}/${B_PSORT}"

samtools index "${bamfold}/${B_PSORT}"

samtools flagstat "${bamfold}/${B_PSORT}" \
  > "${bamfold}/${B_PSORT}.flagstat.txt"

echo "samtools complete. Flagstat:"
cat "${bamfold}/${B_PSORT}.flagstat.txt"

##--------------------------------------------------------------## 5. BigWig coverage

echo ""
echo "--- 5. BigWig coverage (deeptools) ---"

module purge
module load deeptools/3.5.1

## CPM normalization: comparable enrichment across samples at different depths
bamCoverage \
  --bam "${bamfold}/${B_PSORT}" \
  --normalizeUsing CPM \
  --outFileName "${covfold}/${BW_CPM}" \
  --binSize 10 \
  --numberOfProcessors "${threads}" \
  --ignoreDuplicates \
  --extendReads

## RPGC normalization: 1x genome depth, useful for absolute occupancy comparisons
bamCoverage \
  --bam "${bamfold}/${B_PSORT}" \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 2747877702 \
  --outFileName "${covfold}/${BW_RPGC}" \
  --binSize 10 \
  --numberOfProcessors "${threads}" \
  --ignoreDuplicates \
  --extendReads

echo "BigWig files generated."

##--------------------------------------------------------------## 6. FastQC — trimmed reads

echo ""
echo "--- 6. FastQC (trimmed reads) ---"

module purge
module load fastqc/0.12.0

fastqc -t "${threads}" \
  "${outfold}/${T1}" \
  "${outfold}/${T2}" \
  -o "${qcfold}"

echo "FastQC trimmed reads complete."

##--------------------------------------------------------------## Done

echo ""
echo "========================================================"
echo "$(date): Pipeline complete for ${SAMPLE}"
echo ""
echo "  Position-sorted BAM : ${bamfold}/${B_PSORT}"
echo "  BigWig CPM          : ${covfold}/${BW_CPM}"
echo "  BigWig RPGC         : ${covfold}/${BW_RPGC}"
echo "  HISAT2 log          : ${HISAT_LOG}"
echo "  Flagstat            : ${bamfold}/${B_PSORT}.flagstat.txt"
echo "========================================================"
