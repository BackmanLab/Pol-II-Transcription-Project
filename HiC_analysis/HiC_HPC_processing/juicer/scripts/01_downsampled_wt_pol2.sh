#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=150G
#SBATCH -t 2880
#SBATCH --array=1-4
#SBATCH --job-name=downsample_wt_pol2
#SBATCH -o logs/downsample_wt_pol2_%a.out
#SBATCH -e logs/downsample_wt_pol2_%a.err
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Downsamples WT and Pol2 FASTQs to matched depth
# (Pol2 = target ~971M reads). Array: 1=WT_Rep1  2=WT_Rep2  3=Pol2_Rep1  4=Pol2_Rep2
# Usage: sbatch 01_downsampled_wt_pol2.sh
#        sbatch --array=1-2 01_downsampled_wt_pol2.sh   (WT only)

mkdir -p logs

module load seqtk/1.4
module load pigz 2>/dev/null || true

# ── Set compressor ────────────────────────────────────────────────────────
if command -v pigz &>/dev/null; then
    COMPRESS="pigz -p ${SLURM_NTASKS_PER_NODE}"
    echo "Compressor: pigz (${SLURM_NTASKS_PER_NODE} threads)"
else
    COMPRESS="gzip"
    echo "Compressor: gzip (pigz not available)"
fi

FASTQ_DIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp
OUTDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/downsampled_wt_pol2
SEED=42

# Array index → condition, rep, fraction
SAMPLES=(
    "WT_HCT116_CTRL 1 0.606033"
    "WT_HCT116_CTRL 2 0.606033"
    "Pol2_6hrs_Aux  1 1.000000"
    "Pol2_6hrs_Aux  2 1.000000"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
COND=$(echo "${SAMPLES[$IDX]}" | awk '{print $1}')
REP=$(echo  "${SAMPLES[$IDX]}" | awk '{print $2}')
FRAC=$(echo "${SAMPLES[$IDX]}" | awk '{print $3}')

echo "Task ${SLURM_ARRAY_TASK_ID}: ${COND} Rep${REP}  fraction=${FRAC}"

OUTDIR_COND=${OUTDIR}/${COND}/Rep${REP}
mkdir -p ${OUTDIR_COND}

R1=$(ls ${FASTQ_DIR}/${COND}/*R1*Rep${REP}*.fastq.gz)
R2=$(ls ${FASTQ_DIR}/${COND}/*R2*Rep${REP}*.fastq.gz)
R1_OUT=${OUTDIR_COND}/$(basename $R1 .fastq.gz)_ds.fastq.gz
R2_OUT=${OUTDIR_COND}/$(basename $R2 .fastq.gz)_ds.fastq.gz

echo "  R1: $R1"
echo "  R2: $R2"

# ── Helper: fast completion check ────────────────────────────────────────
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
