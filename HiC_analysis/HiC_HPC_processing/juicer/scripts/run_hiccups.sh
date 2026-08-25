#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Submits GPU-accelerated HiCCUPS loop calling for
# merged (mega) .hic files per condition and individual replicate
# .hic files.
#
# Output written to a dedicated top-level loops/ directory:
#   ${BASE}/loops/<dataset>/<condition>/<mega|Rep1|Rep2>/merged_loops.bedpe
#
# Uses SCALE normalization (default in juicer_tools 2.x addNorm).
# Runs on genomics-gpu partition (A100 GPUs) using CUDA 11.8.
#
# Usage: bash run_hiccups.sh

BASE=/projects/b1042/BackmanLab/Lucas/060825_hic-samp
JUICER_TOOLS=${BASE}/juicer/scripts/juicer_tools
LOG_DIR=${BASE}/logs/hiccups
LOOPS_DIR=${BASE}/loops
mkdir -p ${LOG_DIR}

CHROMS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"

# ── Helper: submit one HiCCUPS SBATCH job ──────────────────────────────────
# Usage: submit_hiccups <label> <hic_file> [dependency_jid]
submit_hiccups() {
    local LABEL=$1
    local HIC=$2
    local DEP_JID=$3
    local OUTDIR=${LOOPS_DIR}/${LABEL}
    local SAFE_LABEL=${LABEL//\//_}
    local DEP_FLAG=""

    if [[ ! -f "$HIC" ]]; then
        echo "  SKIP ${LABEL}: inter_30.hic not found" >&2
        echo ""
        return
    fi

    [[ -n "$DEP_JID" ]] && DEP_FLAG="--dependency=afterok:${DEP_JID}"
    mkdir -p ${OUTDIR}

    JID=$(sbatch --account=b1042 ${DEP_FLAG} --parsable << EOF
#!/bin/bash
#SBATCH -p genomics-gpu
#SBATCH -A b1042
#SBATCH --mem=150G
#SBATCH --gres=gpu:a100:1
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -o ${LOG_DIR}/hiccups_${SAFE_LABEL}-%j.out
#SBATCH -e ${LOG_DIR}/hiccups_${SAFE_LABEL}-%j.err
#SBATCH -J hiccups_${SAFE_LABEL}

module load cuda/11.8.0-gcc-12.3.0
module load java/jdk11.0.10

export IBM_JAVA_OPTIONS="-Xmx120000m -Xgcthreads1"
export _JAVA_OPTIONS="-Xms10000m -Xmx120000m"

echo "Running HiCCUPS on: ${HIC}"
echo "Output directory:   ${OUTDIR}"
echo "CUDA device:        \$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
date

${JUICER_TOOLS} hiccups \
    -k SCALE \
    -r 5000,10000,25000 \
    -c ${CHROMS} \
    ${HIC} \
    ${OUTDIR}

EXIT_CODE=\$?
date
if [[ \$EXIT_CODE -eq 0 ]]; then
    echo "(-: HiCCUPS completed for ${LABEL}"
    [[ -f "${OUTDIR}/merged_loops.bedpe" ]] && \
        echo "    Loops: \$(wc -l < ${OUTDIR}/merged_loops.bedpe)" || \
        echo "    WARNING: no merged_loops.bedpe produced"
else
    echo "***! HiCCUPS failed for ${LABEL} (exit code \$EXIT_CODE)"
    exit \$EXIT_CODE
fi
EOF
)
    echo "$JID"
}

# ── Process one condition group ─────────────────────────────────────────────
process_group() {
    local DATASET=$1
    local COND=$2
    local MEGA_HIC=$3
    local REP1_HIC=$4
    local REP2_HIC=$5
    local LABEL="${DATASET}/${COND}"

    echo ""
    echo "═══ ${LABEL} ═══"

    MEGA_JID=$(submit_hiccups "${LABEL}/mega" "${MEGA_HIC}")
    [[ -n "$MEGA_JID" ]] && echo "  mega  → job ${MEGA_JID}"

    REP1_JID=$(submit_hiccups "${LABEL}/Rep1" "${REP1_HIC}" "${MEGA_JID}")
    [[ -n "$REP1_JID" ]] && echo "  Rep1  → job ${REP1_JID}"

    REP2_JID=$(submit_hiccups "${LABEL}/Rep2" "${REP2_HIC}" "${MEGA_JID}")
    [[ -n "$REP2_JID" ]] && echo "  Rep2  → job ${REP2_JID}"
}

# ── ds_wt_pol2 ──────────────────────────────────────────────────────────────
process_group "ds_wt_pol2" "WT_HCT116_CTRL" \
    ${BASE}/mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic \
    ${BASE}/juicer_ds_wt_pol2_genomics/WT_HCT116_CTRL_Rep1/aligned/inter_30.hic \
    ${BASE}/juicer_ds_wt_pol2_genomics/WT_HCT116_CTRL_Rep2/aligned/inter_30.hic

process_group "ds_wt_pol2" "Pol2_6hrs_Aux" \
    ${BASE}/mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic \
    ${BASE}/juicer_ds_wt_pol2_genomics/Pol2_6hrs_Aux_Rep1/aligned/inter_30.hic \
    ${BASE}/juicer_ds_wt_pol2_genomics/Pol2_6hrs_Aux_Rep2/aligned/inter_30.hic

# ── ds_all ───────────────────────────────────────────────────────────────────
process_group "ds_all" "WT_HCT116_CTRL" \
    ${BASE}/mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic \
    ${BASE}/juicer_ds_all/WT_HCT116_CTRL_Rep1/aligned/inter_30.hic \
    ${BASE}/juicer_ds_all/WT_HCT116_CTRL_Rep2/aligned/inter_30.hic

process_group "ds_all" "Pol2_6hrs_Aux" \
    ${BASE}/mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic \
    ${BASE}/juicer_ds_all/Pol2_6hrs_Aux_Rep1/aligned/inter_30.hic \
    ${BASE}/juicer_ds_all/Pol2_6hrs_Aux_Rep2/aligned/inter_30.hic

process_group "ds_all" "1hr_ActD" \
    ${BASE}/mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic \
    ${BASE}/juicer_ds_all/1hr_ActD_Rep1/aligned/inter_30.hic \
    ${BASE}/juicer_ds_all/1hr_ActD_Rep2/aligned/inter_30.hic

# ── unsampled ────────────────────────────────────────────────────────────────
process_group "unsampled" "WT_HCT116_CTRL" \
    ${BASE}/mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic \
    ${BASE}/juicer_unsampled/WT_HCT116_CTRL_Rep1/aligned/inter_30.hic \
    ${BASE}/juicer_unsampled/WT_HCT116_CTRL_Rep2/aligned/inter_30.hic

process_group "unsampled" "Pol2_6hrs_Aux" \
    ${BASE}/mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic \
    ${BASE}/juicer_unsampled/Pol2_6hrs_Aux_Rep1/aligned/inter_30.hic \
    ${BASE}/juicer_unsampled/Pol2_6hrs_Aux_Rep2/aligned/inter_30.hic

process_group "unsampled" "1hr_ActD" \
    ${BASE}/mega_unsampled/1hr_ActD/mega/aligned/inter_30.hic \
    ${BASE}/juicer_unsampled/1hr_ActD_Rep1/aligned/inter_30.hic \
    ${BASE}/juicer_unsampled/1hr_ActD_Rep2/aligned/inter_30.hic

echo ""
echo "All HiCCUPS jobs submitted."
echo "Output: ${LOOPS_DIR}/<dataset>/<condition>/<mega|Rep1|Rep2>/merged_loops.bedpe"
echo "Logs:   ${LOG_DIR}/"
