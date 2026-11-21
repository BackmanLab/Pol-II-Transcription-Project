#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=150G
#SBATCH -t 12:00:00
#SBATCH --array=0-6
#SBATCH --job-name=euseq_hg38
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/%x_%a_out.txt
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs/%x_%a_err.txt
#SBATCH --mail-type=END,FAIL
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Pol II EU-seq alignment to hg38 (paired-end). Input FASTQs are the
# read-name deduplicated ones from dedup_euseq_fastq_hg38.sh (the delivered
# FASTQs contain repeated reads). Same steps as
# align_euseq.sh in the Pol-II-Transcription-Project repository:
# FastQC -> TrimGalore (--paired --retain_unpaired) -> HISAT2 (--dta)
# -> samtools filter (-F 3852 -q 30)/sort/index. chrM is NOT removed.
# Ends with a chrM read-count check (idxstats) per sample.
# Does NOT generate bigWigs or run counting.
#
# Scope: POLR2A_8hrs_aux (BM23-25) + WT (BM29-31) + No_EU_ctrl.
# POLR1A and ActD samples are not aligned here.
#
# SLURM array index -> sample:
#   0-2: POLR2A_8hrs_aux BM23, BM24, BM25
#   3-5: WT BM29, BM30, BM31
#   6:   No_EU_ctrl
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/logs
#
# Usage:
#   sbatch /projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/scripts/align_euseq_hg38.sh
##################################################################

module purge

##--------------------------------------------------------------## Sample array (R1 fastq only -- R2 derived below)

FASTQ_ARRAY=(
  BM23_S40_R1_001.fastq.gz       ## 0: POLR2A_8hrs_aux Rep1
  BM24_S41_R1_001.fastq.gz       ## 1: POLR2A_8hrs_aux Rep2
  BM25_S42_R1_001.fastq.gz       ## 2: POLR2A_8hrs_aux Rep3
  BM29_S46_R1_001.fastq.gz       ## 3: WT Rep1
  BM30_S47_R1_001.fastq.gz       ## 4: WT Rep2
  BM31_S48_R1_001.fastq.gz       ## 5: WT Rep3
  No_EU_ctrl_S49_R1_001.fastq.gz ## 6: No_EU control
)

R1="${FASTQ_ARRAY[$SLURM_ARRAY_TASK_ID]}"
THREADS=10

##--------------------------------------------------------------## Paths

BASE="/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq"
FASTQ_DIR="${BASE}/euseq/fastq_dedup"   ## read-name deduplicated (dedup_euseq_fastq_hg38.sh)
GENOME="${BASE}/reference/hisat/hg38"
GTF="${BASE}/reference/Homo_sapiens.GRCh38.115.chr.gtf"
root="${BASE}/euseq/"

outfold="${root}OUT";      mkdir -p "${outfold}"
samfold="${root}SAM";      mkdir -p "${samfold}"
qcfold="${root}fastqc";    mkdir -p "${qcfold}"
bamfold="${root}BAM";      mkdir -p "${bamfold}"
logfold="${root}logs";     mkdir -p "${logfold}"

[[ -f "${FASTQ_DIR}/${R1}" ]] || { echo "*** error: R1 not found: ${FASTQ_DIR}/${R1} ***" >&2; exit 1; }
[[ -f "${GENOME}.1.ht2" ]]    || { echo "*** error: hg38 HISAT2 index not found at ${GENOME} ***" >&2; exit 1; }
[[ -f "${GTF}" ]]             || { echo "*** error: hg38 GTF not found: ${GTF} ***" >&2; exit 1; }

##--------------------------------------------------------------## Derive R2 (paired-end)

R2=$(echo "${R1}" | sed 's/R1_001/R2_001/g')
[[ -f "${FASTQ_DIR}/${R2}" ]] || { echo "*** error: R2 not found: ${FASTQ_DIR}/${R2} ***" >&2; exit 1; }

##--------------------------------------------------------------## Derive file names

SAMPLE=$(echo "${R1}" | sed 's/_R1_001.fastq.gz//g')
T1=$(echo "${R1}" | sed 's/R1_001.fastq.gz/R1_001_val_1.fq/g')
T2=$(echo "${R2}" | sed 's/R2_001.fastq.gz/R2_001_val_2.fq/g')
S1="${SAMPLE}_R1_001.sam"
B_FILT="${SAMPLE}_filt.bam"
B_PSORT="${SAMPLE}_filt.sort.bam"
HISAT_LOG="${logfold}/${SAMPLE}_hisat2.log"

echo "========================================================"
echo "$(date): EU-seq alignment to hg38 (paired-end)"
echo "  Sample  : ${SAMPLE}"
echo "  R1      : ${R1}"
echo "  R2      : ${R2}"
echo "  Root    : ${root}"
echo "========================================================"

##--------------------------------------------------------------## 1. FastQC -- raw reads

module load fastqc/0.12.0
fastqc -t "${THREADS}" "${FASTQ_DIR}/${R1}" "${FASTQ_DIR}/${R2}" -o "${qcfold}"
echo "FastQC (raw) complete."

##--------------------------------------------------------------## 2. TrimGalore (paired-end)

module purge
module load TrimGalore/0.6.10

trim_galore \
  --paired \
  --retain_unpaired \
  --dont_gzip \
  --cores "${THREADS}" \
  -o "${outfold}" \
  "${FASTQ_DIR}/${R1}" \
  "${FASTQ_DIR}/${R2}"

[[ -f "${outfold}/${T1}" ]] || { echo "*** error: trimmed R1 not found: ${outfold}/${T1} ***" >&2; exit 1; }
[[ -f "${outfold}/${T2}" ]] || { echo "*** error: trimmed R2 not found: ${outfold}/${T2} ***" >&2; exit 1; }
echo "TrimGalore complete."

##--------------------------------------------------------------## 3. HISAT2 alignment (paired-end, hg38)

module purge
module load hisat2/2.1.0

hisat2 \
  -p "${THREADS}" \
  --dta \
  -x "${GENOME}" \
  -1 "${outfold}/${T1}" \
  -2 "${outfold}/${T2}" \
  -S "${samfold}/${S1}" \
  2> "${HISAT_LOG}"

echo "HISAT2 complete. Alignment summary:"
cat "${HISAT_LOG}"

##--------------------------------------------------------------## 4. Filter, sort, index

module purge
module load samtools/1.14

## -F 3852: unmapped, mate unmapped, not primary, fails QC, dup, supplementary
## -q 30: MAPQ >= 30
samtools view \
  -bS \
  -F 3852 \
  -q 30 \
  "${samfold}/${S1}" \
  -o "${bamfold}/${B_FILT}"

samtools sort \
  -@ "${THREADS}" \
  "${bamfold}/${B_FILT}" \
  -o "${bamfold}/${B_PSORT}"

samtools index "${bamfold}/${B_PSORT}"

samtools flagstat "${bamfold}/${B_PSORT}" \
  > "${bamfold}/${B_PSORT}.flagstat.txt"

echo "samtools complete. Flagstat:"
cat "${bamfold}/${B_PSORT}.flagstat.txt"

##--------------------------------------------------------------## 5. FastQC -- trimmed reads

module purge
module load fastqc/0.12.0

fastqc -t "${THREADS}" "${outfold}/${T1}" "${outfold}/${T2}" -o "${qcfold}"
echo "FastQC (trimmed) complete."

##--------------------------------------------------------------## 6. chrM check

module purge
module load samtools/1.14

samtools idxstats "${bamfold}/${B_PSORT}" > "${bamfold}/${B_PSORT}.idxstats.txt"
CHRM=$(awk '$1=="chrM"{s+=$3} END{print s+0}' "${bamfold}/${B_PSORT}.idxstats.txt")
TOTAL=$(awk '{s+=$3} END{print s+0}' "${bamfold}/${B_PSORT}.idxstats.txt")
echo "chrM check: ${CHRM} chrM reads of ${TOTAL} mapped (${SAMPLE})"
[[ "${TOTAL}" -gt 0 ]] || { echo "*** error: no mapped reads in ${B_PSORT} ***" >&2; exit 1; }
[[ "${CHRM}" -gt 0 ]]  || echo "*** WARNING: zero chrM reads in ${B_PSORT} ***" >&2

##--------------------------------------------------------------## Done

echo "========================================================"
echo "$(date): Alignment complete for ${SAMPLE}"
echo "  Sorted BAM : ${bamfold}/${B_PSORT}"
echo "  Flagstat   : ${bamfold}/${B_PSORT}.flagstat.txt"
echo "  HISAT2 log : ${HISAT_LOG}"
echo "========================================================"
