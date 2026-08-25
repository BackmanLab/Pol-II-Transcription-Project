#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Runs Juicer on all 3 downsampled conditions (WT, Pol2,
# ActD; ds_all group), genomics-himem partition variant of
# 02_run_juicer_ds_all.sh. Input: downsampled/ (depth-matched to ActD
# ~712M reads per condition). Output: juicer_ds_all/<COND>_Rep<N>/
#
# IMPORTANT: Run from LOGIN NODE -- juicer.sh submits its own SLURM jobs.
# Usage: bash 02_run_juicer_ds_all_genhm.sh              # full pipeline
#        bash 02_run_juicer_ds_all_genhm.sh -S merge     # restart from merge stage
#        bash 02_run_juicer_ds_all_genhm.sh -S afterdedup # restart from post-dedup
# Prerequisites: 01_downsample_fastqs.sh complete and verified

# ── Parse optional -S stage flag ──────────────────────────────────────────
STAGE_FLAG=""
while getopts "S:" opt; do
    case $opt in
        S) STAGE_FLAG="-S $OPTARG" ;;
        *) echo "Usage: bash $0 [-S stage]" >&2; exit 1 ;;
    esac
done

[[ -n "$STAGE_FLAG" ]] && echo "Restarting from stage: $STAGE_FLAG"

JUICER_SCRIPT=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer/scripts/juicer_genhm.sh
JUICER_DIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer
GENOME_FASTA=${JUICER_DIR}/references/hg38.fa
CHROM_SIZES=${JUICER_DIR}/references/hg38.chrom.sizes
SITE_FILE=${JUICER_DIR}/restriction_sites/hg38_DpnII.txt

FASTQ_DS=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/downsampled
WORKDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer_ds_all

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
echo "All jobs submitted for ds_all."
echo "Monitor: squeue -u \$USER"
echo "After completion, run: bash 03_merge_mega_ds_all.sh"
