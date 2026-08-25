#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Diagnostic test for HiCCUPS on a single chromosome
# and resolution. Runs HiCCUPS on chr1 at 10kb using SCALE
# normalization. Checks whether output files are produced and
# reports results.
#
# Usage: sbatch test_hiccups.sh
#SBATCH -p genomics-gpu
#SBATCH -A b1042
#SBATCH --mem=50G
#SBATCH --gres=gpu:a100:1
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH -o /projects/b1042/BackmanLab/Lucas/060825_hic-samp/logs/test_hiccups-%j.out
#SBATCH -e /projects/b1042/BackmanLab/Lucas/060825_hic-samp/logs/test_hiccups-%j.err
#SBATCH -J test_hiccups

module load cuda/11.8.0-gcc-12.3.0
module load java/jdk11.0.10

export IBM_JAVA_OPTIONS="-Xmx40000m -Xgcthreads1"
export _JAVA_OPTIONS="-Xms5000m -Xmx40000m"

BASE=/projects/b1042/BackmanLab/Lucas/060825_hic-samp
JUICER_TOOLS=${BASE}/juicer/scripts/juicer_tools
HIC=${BASE}/mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic
OUTDIR=${BASE}/logs/test_hiccups_output

mkdir -p ${OUTDIR}

echo "══════════════════════════════════════════════════"
echo "HiCCUPS Diagnostic Test"
echo "══════════════════════════════════════════════════"
echo "HiC file: ${HIC}"
echo "Output:   ${OUTDIR}"
echo "CUDA:     $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null)"
echo ""
date

# ── Test 1: Verify SCALE norm is readable at chr1 10kb ─────────────────────
echo "── Test 1: SCALE normalization available at 10kb? ──"
NORM_LINES=$(${JUICER_TOOLS} dump norm SCALE ${HIC} 1 BP 10000 2>/dev/null | wc -l)
if [[ ${NORM_LINES} -gt 0 ]]; then
    echo "  PASS — ${NORM_LINES} norm values at chr1 10kb"
else
    echo "  FAIL — no SCALE norm values at chr1 10kb"
    exit 1
fi

# ── Test 2: Verify contact data is readable at chr1 10kb ───────────────────
echo ""
echo "── Test 2: Contact data available at chr1 10kb? ──"
CONTACT_LINES=$(${JUICER_TOOLS} dump observed SCALE ${HIC} 1 1 BP 10000 2>/dev/null | wc -l)
if [[ ${CONTACT_LINES} -gt 0 ]]; then
    echo "  PASS — ${CONTACT_LINES} contact records at chr1 10kb"
else
    echo "  FAIL — no contact data at chr1 10kb with SCALE"
    exit 1
fi

# ── Test 3: Run HiCCUPS on chr1 only at 10kb ───────────────────────────────
echo ""
echo "── Test 3: HiCCUPS on chr1 at 10kb only ──"
echo "  Running..."

${JUICER_TOOLS} hiccups \
    -k SCALE \
    -r 10000 \
    -c 1 \
    ${HIC} \
    ${OUTDIR}

EXIT_CODE=$?
echo ""
echo "── Results ──"
echo "  Exit code: ${EXIT_CODE}"
echo ""
echo "  Output files:"
ls -lh ${OUTDIR}/ 2>/dev/null || echo "  (no output directory)"
echo ""

if [[ -f "${OUTDIR}/merged_loops.bedpe" ]]; then
    LOOP_COUNT=$(wc -l < ${OUTDIR}/merged_loops.bedpe)
    echo "  PASS — merged_loops.bedpe produced: ${LOOP_COUNT} loops"
elif ls ${OUTDIR}/enriched_pixels_10000*.bedpe 2>/dev/null | head -1 | grep -q bedpe; then
    echo "  PARTIAL — enriched pixels found but no merged file"
    ls -lh ${OUTDIR}/*.bedpe 2>/dev/null
    wc -l ${OUTDIR}/*.bedpe 2>/dev/null
else
    echo "  FAIL — no loop files produced"
    echo ""
    echo "  fdr_thresholds contents:"
    cat ${OUTDIR}/fdr_thresholds_10000 2>/dev/null || echo "  (none)"
fi

echo ""
date
echo "══════════════════════════════════════════════════"
echo "Test complete."
