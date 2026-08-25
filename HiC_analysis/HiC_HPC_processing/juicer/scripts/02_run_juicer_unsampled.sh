#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Runs Juicer on the original unsampled FASTQs for all
# 3 conditions. Input: original FASTQs in 060825_hic-samp/<COND>/.
# Output: juicer_unsampled/<COND>_Rep<N>/
#
# Use these .hic files when depth-matching is not required (e.g. visualization,
# compartment analysis, TAD calling at maximum resolution).
#
# IMPORTANT: Run from LOGIN NODE -- juicer.sh submits its own SLURM jobs.
# Usage: bash 02_run_juicer_unsampled.sh

JUICER_SCRIPT=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer/scripts/juicer.sh
JUICER_DIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer
GENOME_FASTA=${JUICER_DIR}/references/hg38.fa
CHROM_SIZES=${JUICER_DIR}/references/hg38.chrom.sizes
SITE_FILE=${JUICER_DIR}/restriction_sites/hg38_DpnII.txt

FASTQ_RAW=/projects/b1042/BackmanLab/Lucas/060825_hic-samp
WORKDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer_unsampled

# Condition directory names and their FASTQ subdirectories
declare -A COND_DIRS=(
    [WT_HCT116_CTRL]="WT_HCT116_CTRL"
    [Pol2_6hrs_Aux]="Pol2_6hrs_Aux"
    [1hr_ActD]="1hr_ActD"
)

CONDITIONS=(WT_HCT116_CTRL Pol2_6hrs_Aux 1hr_ActD)

# ── Verify inputs ─────────────────────────────────────────────────────────
echo "Verifying reference files..."
for f in "$JUICER_SCRIPT" "$GENOME_FASTA" "$CHROM_SIZES" "$SITE_FILE"; do
    [[ -f "$f" ]] && echo "  ✓ $f" || { echo "  ✗ MISSING: $f" >&2; exit 1; }
done

# ── Submit Juicer for each condition × replicate ──────────────────────────
for COND in "${CONDITIONS[@]}"; do
    for REP in 1 2; do

        SAMPLE="${COND}_Rep${REP}"
        SAMPLE_DIR=${WORKDIR}/${SAMPLE}
        FASTQ_DIR=${SAMPLE_DIR}/fastq

        echo ""
        echo "═══════════════════════════════════════════════"
        echo "Submitting: ${SAMPLE}"

        mkdir -p ${FASTQ_DIR}

        # Symlink only Rep-specific FASTQs into fastq/ subdirectory
        # Original files follow pattern: <COND>_R{1,2}_Rep<N>.fastq.gz
        for F in ${FASTQ_RAW}/${COND_DIRS[$COND]}/*Rep${REP}.fastq.gz; do
            [[ -f "$F" ]] || { echo "ERROR: No FASTQs found for ${COND} Rep${REP}" >&2; exit 1; }
            ln -sf $F ${FASTQ_DIR}/$(basename $F)
        done

        echo "  FASTQs:"
        ls ${FASTQ_DIR}/

        bash ${JUICER_SCRIPT} \
            -d ${SAMPLE_DIR} \
            -g hg38 \
            -z ${GENOME_FASTA} \
            -p ${CHROM_SIZES} \
            -s DpnII \
            -y ${SITE_FILE} \
            -A b1042 &

        echo "  Jobs submitted for ${SAMPLE} (pid $!)"
    done
done

wait

echo ""
echo "All jobs submitted."
echo "Monitor: squeue -u $USER"
echo "After completion, run: bash 03_merge_mega_unsampled.sh"
