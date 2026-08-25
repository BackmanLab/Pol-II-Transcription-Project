#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=hic_downsample
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=150G
#SBATCH -t 2880
#SBATCH --array=1-6
#SBATCH --output=logs/downsample_%a.out
#SBATCH --error=logs/downsample_%a.err
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Downsamples each replicate to matched sequencing depth
# using seqtk (ds_all group -- ActD/Pol2/WT). Array index 1-6:
# 1=ActD_Rep1  2=ActD_Rep2  3=Pol2_Rep1  4=Pol2_Rep2  5=WT_Rep1  6=WT_Rep2
#
# Prerequisites: 00_count_reads.sh + 00b_make_samplesheet.sh must be complete
# Usage:         sbatch 01_downsample_fastqs.sh
#                sbatch --array=1-2 01_downsample_fastqs.sh   (ActD only)
#                sbatch --array=3-5 01_downsample_fastqs.sh   (Pol2+WT only)

mkdir -p logs

module load seqtk/1.4

# Use pigz for parallel compression if available, otherwise gzip
module load pigz 2>/dev/null || true
if command -v pigz &>/dev/null; then
    COMPRESS="pigz -p ${SLURM_NTASKS_PER_NODE}"
    echo "Compressor: pigz (${SLURM_NTASKS_PER_NODE} threads)"
else
    COMPRESS="gzip"
    echo "Compressor: gzip (pigz not available)"
fi

OUTDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/downsampled
SHEET=sample_sheet_fractions.txt
SEED=42

mkdir -p $OUTDIR

# ── Read sample sheet line for this array task ────────────────────────────
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SHEET)
COND=$(echo "$LINE" | awk '{print $1}')
REP=$(echo  "$LINE" | awk '{print $2}')
R1=$(echo   "$LINE" | awk '{print $3}')
R2=$(echo   "$LINE" | awk '{print $4}')
FRAC=$(echo "$LINE" | awk '{print $7}')

echo "Task ${SLURM_ARRAY_TASK_ID}: ${COND} Rep${REP}  fraction=${FRAC}"
echo "  R1: $R1"
echo "  R2: $R2"

OUTDIR_COND=${OUTDIR}/${COND}/Rep${REP}
mkdir -p ${OUTDIR_COND}

R1_OUT=${OUTDIR_COND}/$(basename $R1 .fastq.gz)_ds.fastq.gz
R2_OUT=${OUTDIR_COND}/$(basename $R2 .fastq.gz)_ds.fastq.gz

# ── Helper: fast completion check — reads only first FASTQ record ─────────
is_complete() {
    local f=$1
    [[ -s "$f" ]] || return 1
    zcat "$f" 2>/dev/null | head -4 | grep -q "^@"
}

# ── Process R1 ────────────────────────────────────────────────────────────
if is_complete "$R1_OUT"; then
    echo "R1 already complete — skipping"
else
    if (( $(echo "$FRAC >= 0.9999" | bc -l) )); then
        echo "Fraction=1.0 — copying R1..."
        cp $R1 $R1_OUT
    else
        echo "Downsampling R1..."
        seqtk sample -s${SEED} $R1 $FRAC | ${COMPRESS} > $R1_OUT
    fi
fi

# ── Process R2 ────────────────────────────────────────────────────────────
if is_complete "$R2_OUT"; then
    echo "R2 already complete — skipping"
else
    if (( $(echo "$FRAC >= 0.9999" | bc -l) )); then
        echo "Fraction=1.0 — copying R2..."
        cp $R2 $R2_OUT
    else
        echo "Downsampling R2..."
        seqtk sample -s${SEED} $R2 $FRAC | ${COMPRESS} > $R2_OUT
    fi
fi

# ── Verify ────────────────────────────────────────────────────────────────
for f in "$R1_OUT" "$R2_OUT"; do
    if is_complete "$f"; then
        echo "  ✓ $(basename $f)"
    else
        echo "  ✗ FAILED or empty: $(basename $f)" >&2
        exit 1
    fi
done

echo "Done: ${COND} Rep${REP}"
