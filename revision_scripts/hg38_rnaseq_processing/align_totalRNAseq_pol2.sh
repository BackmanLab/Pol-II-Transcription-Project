#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH -t 18:00:00
#SBATCH --array=0-9
#SBATCH --job-name=rnaseq_pol2
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/101823_rnaseq/total/logs/%x_%a_out.txt
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/101823_rnaseq/total/logs/%x_%a_err.txt
#SBATCH --mail-type=END,FAIL
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Total RNA-seq re-alignment to hg38 for the Pol II data
# (single-end). Same steps as total_RNA_seq/Total_RNAseq_HPC_processing/
# align_totalRNAseq.sh: TrimGalore -> HISAT2 (--dta) -> samtools
# filter/sort/index/flagstat. Reads the original fastq directory
# directly; the index is the one built by build_hisat2_index_pol2.sh;
# index in 101823_rnaseq/reference/hisat/; outputs land in 101823_rnaseq/total/ (BAM/, SAM/, OUT/, fastqc/).
#
# SLURM array index -> sample:
#   0-3: POLR2A_8hrs_aux Rep1-4
#   4-6: ActD_6hrs Rep1-3
#   7-9: WT Rep1-3
#
# Before submitting:
#   mkdir -p /projects/b1042/BackmanLab/Lucas/101823_rnaseq/total/logs
#
# Usage:
#   sbatch align_totalRNAseq_pol2.sh
# After all array tasks complete, run:
#   sbatch merge_totalrnaseq_bams_pol2.sh
##################################################################

module purge

##--------------------------------------------------------------## Sample array

FASTQ_ARRAY=(
  Pol2_8hrs_Aux_Rep1_S5_R1_001.fastq.gz     ## 0: POLR2A Rep1
  Pol2_8hrs_Aux_Rep2_S6_R1_001.fastq.gz     ## 1: POLR2A Rep2
  Pol2_8hrs_Aux_Rep3_S7_R1_001.fastq.gz     ## 2: POLR2A Rep3
  Pol2_8hrs_Aux_Rep4_S8_R1_001.fastq.gz     ## 3: POLR2A Rep4
  ActD_6hrs_Rep1_S9_R1_001.fastq.gz         ## 4: ActD Rep1
  ActD_6hrs_Rep2_S10_R1_001.fastq.gz        ## 5: ActD Rep2
  ActD_6hrs_Rep3_S11_R1_001.fastq.gz        ## 6: ActD Rep3
  WT_HCT116_CTRL_Rep1_S12_R1_001.fastq.gz   ## 7: WT Rep1
  WT_HCT116_CTRL_Rep2_S13_R1_001.fastq.gz   ## 8: WT Rep2
  WT_HCT116_CTRL_Rep3_S14_R1_001.fastq.gz   ## 9: WT Rep3
)

R1="${FASTQ_ARRAY[$SLURM_ARRAY_TASK_ID]}"
THREADS=10

##--------------------------------------------------------------## Paths

FASTQ_DIR="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/fastq"
GENOME="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/reference/hisat/hg38"
GTF="/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
root="/projects/b1042/BackmanLab/Lucas/101823_rnaseq/total/"

outfold="${root}OUT";      mkdir -p "${outfold}"
samfold="${root}SAM";      mkdir -p "${samfold}"
qcfold="${root}fastqc";    mkdir -p "${qcfold}"
bamfold="${root}BAM";      mkdir -p "${bamfold}"
logfold="${root}logs";     mkdir -p "${logfold}"

[[ -f "${FASTQ_DIR}/${R1}" ]]  || { echo "*** error: read file not found: ${FASTQ_DIR}/${R1} ***" >&2; exit 1; }
[[ -f "${GENOME}.1.ht2" ]]     || { echo "*** error: hg38 HISAT2 index not found at ${GENOME} -- run build_hisat2_index_pol2.sh first ***" >&2; exit 1; }
[[ -f "${GTF}" ]]              || { echo "*** error: GTF not found: ${GTF} ***" >&2; exit 1; }

SAMPLE=$(echo "${R1}" | sed 's/_R1_001.fastq.gz//g')
T1=$(echo "${R1}" | sed 's/_R1_001.fastq.gz/_R1_001_trimmed.fq/g')
S1="${SAMPLE}_R1_001.sam"
B_FILT="${SAMPLE}_filt.bam"
B_PSORT="${SAMPLE}_filt.sort.bam"
HISAT_LOG="${logfold}/${SAMPLE}_hisat2.log"

echo "========================================================"
echo "$(date): Total RNA-seq alignment to hg38 (single-end)"
echo "  Sample  : ${SAMPLE}"
echo "  Read    : ${R1}"
echo "  Root    : ${root}"
echo "========================================================"

##--------------------------------------------------------------## 1. FastQC -- raw reads

module load fastqc/0.12.0
fastqc -t "${THREADS}" "${FASTQ_DIR}/${R1}" -o "${qcfold}"
echo "FastQC (raw) complete."

##--------------------------------------------------------------## 2. TrimGalore (single-end)

module purge
module load TrimGalore/0.6.10
trim_galore \
  --dont_gzip \
  --cores "${THREADS}" \
  -o "${outfold}" \
  "${FASTQ_DIR}/${R1}"

[[ -f "${outfold}/${T1}" ]] || { echo "*** error: trimmed file not found: ${outfold}/${T1} ***" >&2; exit 1; }
echo "TrimGalore complete."

##--------------------------------------------------------------## 3. HISAT2 alignment (single-end, hg38)

module purge
module load hisat2/2.1.0

hisat2 \
  -p "${THREADS}" \
  --dta \
  -x "${GENOME}" \
  -U "${outfold}/${T1}" \
  -S "${samfold}/${S1}" \
  2> "${HISAT_LOG}"

echo "HISAT2 complete. Alignment summary:"
cat "${HISAT_LOG}"

##--------------------------------------------------------------## 4. Filter, sort, index

module purge
module load samtools/1.14

## -F 3852 removes unmapped, not primary, fails QC, duplicate, supplementary
## -q 30 retains uniquely mapping reads
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

fastqc -t "${THREADS}" "${outfold}/${T1}" -o "${qcfold}"
echo "FastQC (trimmed) complete."

##--------------------------------------------------------------## Done

echo "========================================================"
echo "$(date): Alignment complete for ${SAMPLE}"
echo "  Sorted BAM : ${bamfold}/${B_PSORT}"
echo "  Flagstat   : ${bamfold}/${B_PSORT}.flagstat.txt"
echo "  HISAT2 log : ${HISAT_LOG}"
echo "========================================================"
