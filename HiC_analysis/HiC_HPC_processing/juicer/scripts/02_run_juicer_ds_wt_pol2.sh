#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Runs Juicer on WT and Pol2 downsampled conditions only
# (ds_wt_pol2 group). Input: downsampled_wt_pol2/ (depth-matched to
# Pol2 ~971M reads). Output: juicer_ds_wt_pol2/<COND>_Rep<N>/
#
# IMPORTANT: Run from LOGIN NODE -- juicer.sh submits its own SLURM jobs.
# Usage: bash 02_run_juicer_ds_wt_pol2.sh              # full pipeline
#        bash 02_run_juicer_ds_wt_pol2.sh -S merge     # restart from merge stage
#        bash 02_run_juicer_ds_wt_pol2.sh -S afterdedup # restart from post-dedup
# Prerequisites: 01_downsampled_wt_pol2.sh complete and verified

# ── Parse optional -S stage flag ──────────────────────────────────────────
STAGE_FLAG=""
while getopts "S:" opt; do
    case $opt in
        S) STAGE_FLAG="-S $OPTARG" ;;
        *) echo "Usage: bash $0 [-S stage]" >&2; exit 1 ;;
    esac
done

[[ -n "$STAGE_FLAG" ]] && echo "Restarting from stage: $STAGE_FLAG"

JUICER_SCRIPT=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer/scripts/juicer.sh
JUICER_DIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer
GENOME_FASTA=${JUICER_DIR}/references/hg38.fa
CHROM_SIZES=${JUICER_DIR}/references/hg38.chrom.sizes
SITE_FILE=${JUICER_DIR}/restriction_sites/hg38_DpnII.txt

FASTQ_DS=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/downsampled_wt_pol2
WORKDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer_ds_wt_pol2_genomics

CONDITIONS=(WT_HCT116_CTRL Pol2_6hrs_Aux)

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
        echo "Submitting: ${SAMPLE}${STAGE_FLAG:+ ($STAGE_FLAG)}"

        mkdir -p ${FASTQ_DIR}

        # Symlink downsampled FASTQs into fastq/ subdirectory
        for F in ${FASTQ_DS}/${COND}/Rep${REP}/*_ds.fastq.gz; do
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
            -A b1042 \
            ${STAGE_FLAG} &

        echo "  Jobs submitted for ${SAMPLE} (pid $!)"
    done
done

wait

echo ""
echo "All jobs submitted for ds_wt_pol2."
echo "Monitor: squeue -u \$USER"
echo "After completion, run: bash 03_merge_mega_ds_wt_pol2.sh"
