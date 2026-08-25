#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Calls A/B compartments from merged mega .hic files
# using juicer_tools eigenvector at 500kb resolution. Eigenvector
# signs are oriented using GC content computed from hg38.fa (A
# compartments are GC-rich). Per-chromosome output is concatenated
# into a single genome-wide file.
#
# Steps per .hic file:
#   1. Generate 500kb bins from chrom.sizes (once, shared across all jobs)
#   2. Compute GC content per bin from hg38.fa (once, shared)
#   3. Run juicer_tools eigenvector KR per chromosome
#   4. Orient eigenvector sign using GC content correlation
#   5. Concatenate all chromosomes -> <label>_compartments.bedgraph
#
# CPU only -- runs on genomics partition.
# Prerequisites: mega .hic files complete.
# Usage: bash run_compartments.sh

BASE=/projects/b1042/BackmanLab/Lucas/060825_hic-samp
JUICER_DIR=${BASE}/juicer
JUICER_TOOLS=${JUICER_DIR}/scripts/juicer_tools
CHROM_SIZES=${JUICER_DIR}/references/hg38.chrom.sizes
GENOME_FA=${JUICER_DIR}/references/hg38.fa
LOG_DIR=${BASE}/logs/compartments
COMP_DIR=${BASE}/compartments
REFDIR=${COMP_DIR}/reference
mkdir -p ${LOG_DIR} ${REFDIR}

# Chromosomes to process (bare number format — matches data chromosome naming)
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

RESOLUTION=500000  # 500kb — juicer_tools 2.x eigenvector does not support below 500kb

# ── .hic files to process ─────────────────────────────────────────────────
declare -A HIC_FILES=(
    ["ds_wt_pol2_WT"]="${BASE}/mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"
    ["ds_wt_pol2_Pol2"]="${BASE}/mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"
    ["ds_all_WT"]="${BASE}/mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"
    ["ds_all_Pol2"]="${BASE}/mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"
    ["ds_all_ActD"]="${BASE}/mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic"
    ["unsampled_WT"]="${BASE}/mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic"
    ["unsampled_Pol2"]="${BASE}/mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"
    ["unsampled_ActD"]="${BASE}/mega_unsampled/1hr_ActD/mega/aligned/inter_30.hic"
)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 & 2: Generate 500kb bins and GC content (shared reference, run once)
# ══════════════════════════════════════════════════════════════════════════════
GC_BED=${REFDIR}/hg38_500kb_gc.bed
BINS_BED=${REFDIR}/hg38_500kb_bins.bed

REF_JID=""

if [[ ! -f "$GC_BED" ]]; then
    echo "Submitting GC content reference job..."
    REF_JID=$(sbatch --account=b1042 --parsable << EOF
#!/bin/bash
#SBATCH -p genomics
#SBATCH -A b1042
#SBATCH --mem=50G
#SBATCH -c 4
#SBATCH -t 2:00:00
#SBATCH -o ${LOG_DIR}/gc_reference-%j.out
#SBATCH -e ${LOG_DIR}/gc_reference-%j.err
#SBATCH -J compartments_gc_ref

module load bedtools/2.30.0
module load samtools/1.6

echo "Creating 500kb bins..."
bedtools makewindows \
    -g ${CHROM_SIZES} \
    -w ${RESOLUTION} \
    > ${BINS_BED}

echo "Computing GC content per bin..."
bedtools nuc \
    -fi ${GENOME_FA} \
    -bed ${BINS_BED} \
    > ${GC_BED}

echo "Done. Bins: \$(wc -l < ${BINS_BED})"
EOF
)
    echo "  GC reference job → ${REF_JID}"
else
    echo "GC content reference already exists, skipping."
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3-5: Eigenvector + orientation + concatenation per .hic file
# ══════════════════════════════════════════════════════════════════════════════

# Inline Python script for sign orientation and concatenation
# Written once and reused by all jobs
ORIENT_SCRIPT=${COMP_DIR}/orient_eigenvector.py
cat > ${ORIENT_SCRIPT} << 'PYEOF'
#!/usr/bin/env python3
"""
Orient eigenvector sign using GC content and concatenate chromosomes.

Usage: python3 orient_eigenvector.py <gc_bed> <eigen_dir> <output_bedgraph> <resolution>

Reads per-chromosome eigenvector files from eigen_dir, correlates each with
GC content from gc_bed, flips sign if correlation is negative (A = GC-rich
= positive eigenvector), and writes a genome-wide bedgraph.
"""
import sys
import os
import numpy as np

gc_bed    = sys.argv[1]   # bedtools nuc output
eigen_dir = sys.argv[2]   # directory with per-chr eigenvector .txt files
outfile   = sys.argv[3]   # output bedgraph
resolution= int(sys.argv[4])

# Parse GC content: cols 0=chr, 1=start, 2=end, 4=%AT, 5=%GC (0-indexed)
gc = {}
with open(gc_bed) as f:
    for line in f:
        if line.startswith('#') or line.startswith('chr'):
            continue
        parts = line.strip().split('\t')
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
        try:
            gc_frac = float(parts[5])  # %GC column (bedtools nuc col 6, 0-indexed 5)
        except (IndexError, ValueError):
            continue
        gc.setdefault(chrom, []).append((start, end, gc_frac))

chroms_written = 0
skipped = []

with open(outfile, 'w') as out:
    for eigfile in sorted(os.listdir(eigen_dir)):
        if not eigfile.endswith('_eigen.txt'):
            continue
        chrom = eigfile.replace('_eigen.txt', '')
        fpath = os.path.join(eigen_dir, eigfile)

        # Read eigenvector values
        values = []
        with open(fpath) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    values.append(float(line))
                except ValueError:
                    values.append(np.nan)

        if not values or all(np.isnan(v) for v in values):
            skipped.append(chrom)
            continue

        # Match GC content for this chromosome
        gc_vals = gc.get(chrom, [])
        if not gc_vals:
            skipped.append(f"{chrom}(no GC)")
            continue

        # Align by bin index
        n = min(len(values), len(gc_vals))
        eig_arr = np.array(values[:n])
        gc_arr  = np.array([g[2] for g in gc_vals[:n]])

        # Mask NaN
        mask = ~np.isnan(eig_arr)
        if mask.sum() < 2:
            skipped.append(f"{chrom}(all NaN)")
            continue

        corr = np.corrcoef(eig_arr[mask], gc_arr[mask])[0, 1]
        sign = 1 if corr >= 0 else -1

        # Write oriented bedgraph
        for i, val in enumerate(values):
            if i < len(gc_vals):
                start = gc_vals[i][0]
                end   = gc_vals[i][1]
            else:
                start = i * resolution
                end   = start + resolution
            oriented = (sign * val) if not np.isnan(val) else 'NA'
            out.write(f"{chrom}\t{start}\t{end}\t{oriented}\n")

        chroms_written += 1

print(f"Written {chroms_written} chromosomes. Skipped: {skipped if skipped else 'none'}")
PYEOF

echo ""

for LABEL in "${!HIC_FILES[@]}"; do
    HIC=${HIC_FILES[$LABEL]}
    OUTDIR=${COMP_DIR}/${LABEL}

    if [[ ! -f "$HIC" ]]; then
        echo "  SKIP ${LABEL}: inter_30.hic not found"
        continue
    fi

    mkdir -p ${OUTDIR}
    echo "Submitting compartments for: ${LABEL}"

    # Set dependency on GC reference job if it was just submitted
    DEP_FLAG=""
    [[ -n "$REF_JID" ]] && DEP_FLAG="--dependency=afterok:${REF_JID}"

    JID=$(sbatch --account=b1042 ${DEP_FLAG} --parsable << EOF
#!/bin/bash
#SBATCH -p genomics
#SBATCH -A b1042
#SBATCH --mem=50G
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH -o ${LOG_DIR}/compartments_${LABEL}-%j.out
#SBATCH -e ${LOG_DIR}/compartments_${LABEL}-%j.err
#SBATCH -J compartments_${LABEL}

module load java/jdk11.0.10

EIGEN_DIR=${OUTDIR}/eigenvectors
mkdir -p \${EIGEN_DIR}

echo "Computing eigenvectors for: ${LABEL}"
echo "HiC file: ${HIC}"
echo ""

# Run eigenvector per chromosome
for CHR in ${CHROMS[@]}; do
    OUTFILE=\${EIGEN_DIR}/\${CHR}_eigen.txt
    echo "  chr\${CHR}..."
    ${JUICER_TOOLS} eigenvector \
        SCALE \
        ${HIC} \
        \${CHR} \
        BP \
        ${RESOLUTION} \
        \${OUTFILE}
    if [[ \$? -ne 0 ]]; then
        echo "  WARNING: eigenvector failed for chr\${CHR}, skipping"
    fi
done

echo ""
echo "Orienting eigenvectors using GC content and concatenating..."

python3 ${ORIENT_SCRIPT} \
    ${GC_BED} \
    \${EIGEN_DIR} \
    ${OUTDIR}/${LABEL}_compartments.bedgraph \
    ${RESOLUTION}

echo ""
echo "Output: ${OUTDIR}/${LABEL}_compartments.bedgraph"
echo "(-: Compartment calling complete for ${LABEL}"
EOF
)
    echo "  job ${JID} → ${OUTDIR}/${LABEL}_compartments.bedgraph"
done

echo ""
echo "All compartment jobs submitted."
echo "Monitor: squeue -u \$USER -j \$(squeue -u \$USER -h -o '%i' | tr '\n' ',')"
echo "Logs:    ${LOG_DIR}/"
echo "Output:  ${COMP_DIR}/<label>/<label>_compartments.bedgraph"
