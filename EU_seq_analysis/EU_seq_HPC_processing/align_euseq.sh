#!/bin/bash
##################################################################
# Author: Lucas Carter
# Email: lucascarter2025@u.northwestern.edu
# PI: Vadim Backman
# Description:
# EU-seq alignment pipeline (paired-end). Performs trimming,
# HISAT2 alignment, BAM filtering/sorting, and replicate merging.
# Does NOT generate bigWigs — those are generated separately
# using bamCompare with total RNA-seq derived scale factors.
#
# Conditions and replicates:
#   POLR1A_8hrs_aux : BM20, BM21, BM22
#   POLR2A_8hrs_aux : BM23, BM24, BM25
#   ActD_6hrs       : BM27, BM28  (BM26 excluded — no signal)
#   WT              : BM29, BM30, BM31
#   No_EU_ctrl      : No_EU_ctrl_S49 (single sample)
#
# Usage (via sbatch --wrap):
#   sbatch -A p32171 -p genhimem -t 12:00:00 -n 10 --mem=64G \
#     --wrap="sh align_EUseq.sh -f SAMPLE_R1_001.fastq.gz \
#       -w /projects/b1042/BackmanLab/Lucas/090124_euseq/fastq \
#       -g /projects/b1042/BackmanLab/Lucas/090124_euseq/genome/ref_hisat/hg38 \
#       -r /projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf \
#       -t 10"
#
# After all per-sample jobs complete, run the merge step:
#   sh merge_EUseq_bams.sh
##################################################################

set -euo pipefail

##--------------------------------------------------------------## Help

helpFunction()
{
  echo ""
  echo "************** EU-seq alignment pipeline (Quest) **************"
  echo "Usage: sh $0 -f <R1 fastq> -w <fastq dir> -g <HISAT2 index> -r <GTF> -t <threads>"
  echo "  -f   Forward read R1 filename"
  echo "  -w   Full path to fastq directory"
  echo "  -g   Full path to HISAT2 genome index prefix"
  echo "  -r   Full path to GTF annotation file"
  echo "  -t   Number of threads"
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

[[ -z "${R1:-}"      ]] && { echo "*** error: R1 required (-f) ***";          helpFunction; }
[[ -z "${wdir:-}"    ]] && { echo "*** error: fastq directory required (-w) ***"; helpFunction; }
[[ -z "${gendir:-}"  ]] && { echo "*** error: HISAT2 index required (-g) ***";    helpFunction; }
[[ -z "${gtfdir:-}"  ]] && { echo "*** error: GTF required (-r) ***";             helpFunction; }
[[ -z "${threads:-}" ]] && { echo "*** error: threads required (-t) ***";         helpFunction; }

[[ -d "${wdir}" ]]       || { echo "*** error: fastq directory not found: ${wdir} ***"; exit 1; }
[[ -f "${wdir}/${R1}" ]] || { echo "*** error: R1 not found: ${wdir}/${R1} ***";        exit 1; }

##--------------------------------------------------------------## Derive R2 (paired-end)

R2=$(echo "${R1}" | sed 's/R1_001/R2_001/g')
[[ -f "${wdir}/${R2}" ]] || { echo "*** error: R2 not found: ${wdir}/${R2} ***"; exit 1; }

##--------------------------------------------------------------## Directory structure

root="/projects/b1042/BackmanLab/Lucas/090124_euseq/"

outfold="${root}OUT";      mkdir -p "${outfold}"
samfold="${root}SAM";      mkdir -p "${samfold}"
qcfold="${root}fastqc";    mkdir -p "${qcfold}"
bamfold="${root}BAM";      mkdir -p "${bamfold}"
logfold="${root}logs";     mkdir -p "${logfold}"

##--------------------------------------------------------------## Derive file names

SAMPLE=$(echo "${R1}" | sed 's/_R1_001.fastq.gz//g')
T1=$(echo "${R1}" | sed 's/R1_001.fastq.gz/R1_001_val_1.fq/g')
T2=$(echo "${R2}" | sed 's/R2_001.fastq.gz/R2_001_val_2.fq/g')
S1="${SAMPLE}_R1_001.sam"
B_FILT="${SAMPLE}_filt.bam"
B_PSORT="${SAMPLE}_filt.sort.bam"
HISAT_LOG="${logfold}/${SAMPLE}_hisat2.log"

echo "========================================================"
echo "$(date): EU-seq alignment pipeline"
echo "  Sample  : ${SAMPLE}"
echo "  R1      : ${R1}"
echo "  R2      : ${R2}"
echo "  Root    : ${root}"
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

echo "FastQC complete."

##--------------------------------------------------------------## 2. TrimGalore (paired-end)

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

[[ -f "${outfold}/${T1}" ]] || { echo "*** error: trimmed R1 not found ***"; exit 1; }
[[ -f "${outfold}/${T2}" ]] || { echo "*** error: trimmed R2 not found ***"; exit 1; }
echo "TrimGalore complete."

##--------------------------------------------------------------## 3. HISAT2 alignment

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

##--------------------------------------------------------------## 4. Filter, sort, index

echo ""
echo "--- 4. samtools: filter, sort, index ---"
module purge
module load samtools/1.14

## Remove: unmapped, mate unmapped, not primary, fails QC,
##         PCR duplicate, supplementary (-F 3852), MAPQ < 30
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

##--------------------------------------------------------------## 5. FastQC — trimmed reads

echo ""
echo "--- 5. FastQC (trimmed reads) ---"
module purge
module load fastqc/0.12.0

fastqc -t "${threads}" \
  "${outfold}/${T1}" \
  "${outfold}/${T2}" \
  -o "${qcfold}"

echo "FastQC trimmed complete."

##--------------------------------------------------------------## Done

echo ""
echo "========================================================"
echo "$(date): Alignment complete for ${SAMPLE}"
echo "  Sorted BAM : ${bamfold}/${B_PSORT}"
echo "  Flagstat   : ${bamfold}/${B_PSORT}.flagstat.txt"
echo "  HISAT2 log : ${HISAT_LOG}"
echo "========================================================"
