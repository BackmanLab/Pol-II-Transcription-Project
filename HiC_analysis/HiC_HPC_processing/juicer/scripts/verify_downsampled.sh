#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=150G
#SBATCH -t 480
#SBATCH --job-name=verify_ds
#SBATCH -o logs/verify_downsampled.out
#SBATCH -e logs/verify_downsampled.err
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Verifies all downsampled FASTQs across both output
# directories (downsampled/ = 3-condition WT/Pol2/ActD;
# downsampled_wt_pol2/ = 2-condition WT/Pol2). Checks gzip integrity
# (pigz -t) and read counts (zcat | wc -l / 4), run in parallel
# across all files.
#
# Usage: sbatch verify_downsampled.sh

mkdir -p logs

module load pigz

BASE=/projects/b1042/BackmanLab/Lucas/060825_hic-samp
OUTFILE=logs/verify_downsampled.out
THREADS=${SLURM_NTASKS_PER_NODE}

# ── Collect all downsampled FASTQs ────────────────────────────────────────
FILES=$(find \
    ${BASE}/downsampled \
    ${BASE}/downsampled_wt_pol2 \
    -name "*_ds.fastq.gz" | sort)

N_FILES=$(echo "$FILES" | wc -l)
echo "Found ${N_FILES} downsampled FASTQ files"
echo ""

# ── Step 1: Gzip integrity check (parallel) ───────────────────────────────
echo "══════════════════════════════════════════"
echo "Step 1 — Gzip integrity (pigz -t)"
echo "══════════════════════════════════════════"

echo "$FILES" | xargs -P ${THREADS} -I{} bash -c \
    'pigz -t {} 2>&1 && echo "  OK:     {}" || echo "  CORRUPT: {}"'

echo ""

# ── Step 2: Read counts (parallel) ────────────────────────────────────────
echo "══════════════════════════════════════════"
echo "Step 2 — Read counts"
echo "══════════════════════════════════════════"

echo "$FILES" | xargs -P ${THREADS} -I{} bash -c \
    'printf "  %-70s %s reads\n" "{}" "$(( $(zcat {} | wc -l) / 4 ))"' \
    | sort

echo ""

# ── Summary: R1/R2 pair matching ──────────────────────────────────────────
echo "══════════════════════════════════════════"
echo "Step 3 — R1/R2 pair count check"
echo "══════════════════════════════════════════"

for DIR in ${BASE}/downsampled ${BASE}/downsampled_wt_pol2; do
    [[ -d "$DIR" ]] || continue
    echo "Directory: $DIR"
    for COND_DIR in ${DIR}/*/; do
        COND=$(basename $COND_DIR)
        for REP_DIR in ${COND_DIR}/Rep*/; do
            REP=$(basename $REP_DIR)
            R1=$(ls ${REP_DIR}/*R1*_ds.fastq.gz 2>/dev/null)
            R2=$(ls ${REP_DIR}/*R2*_ds.fastq.gz 2>/dev/null)
            [[ -z "$R1" || -z "$R2" ]] && \
                echo "  MISSING: ${COND} ${REP}" && continue
            N_R1=$(( $(zcat $R1 | wc -l) / 4 ))
            N_R2=$(( $(zcat $R2 | wc -l) / 4 ))
            if [[ $N_R1 -eq $N_R2 ]]; then
                echo "  ✓ ${COND} ${REP}: ${N_R1} pairs"
            else
                echo "  ✗ MISMATCH ${COND} ${REP}: R1=${N_R1} R2=${N_R2}" >&2
            fi
        done
    done
    echo ""
done

echo "Verification complete."
