#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: HiCCUPS loop calling with KR normalization. Mirrors
# run_hiccups.sh but uses -k KR, outputs to loops_kr/. Requires KR
# normalization to be stored in .hic files via:
#   juicer_tools addNorm -k KR <hic_file>
#SBATCH --job-name=hiccups_kr
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/logs/hiccups/hiccups_kr_%A_%a.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/logs/hiccups/hiccups_kr_%A_%a.err
#SBATCH --array=0-7

module load java/jdk11.0.10
module load cuda/11.8.0

BASE="/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
JUICER_TOOLS="${BASE}/juicer/scripts/juicer_tools.jar"
OUTDIR="${BASE}/loops_kr"

# ── Chromosomes ───────────────────────────────────────────────────────────────
CHROMS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"

# ── Job array: dataset × condition ────────────────────────────────────────────
# Index: 0=ds_wt_pol2/WT  1=ds_wt_pol2/Pol2
#        2=ds_all/WT       3=ds_all/Pol2       4=ds_all/ActD
#        5=unsampled/WT    6=unsampled/Pol2    7=unsampled/ActD

declare -a HIC_FILES=(
    "${BASE}/mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"
    "${BASE}/mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"
    "${BASE}/mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"
    "${BASE}/mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"
    "${BASE}/mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic"
    "${BASE}/mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic"
    "${BASE}/mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"
    "${BASE}/mega_unsampled/1hr_ActD/mega/aligned/inter_30.hic"
)

declare -a DATASETS=(
    "ds_wt_pol2" "ds_wt_pol2"
    "ds_all"     "ds_all"     "ds_all"
    "unsampled"  "unsampled"  "unsampled"
)

declare -a CONDITIONS=(
    "WT_HCT116_CTRL" "Pol2_6hrs_Aux"
    "WT_HCT116_CTRL" "Pol2_6hrs_Aux" "1hr_ActD"
    "WT_HCT116_CTRL" "Pol2_6hrs_Aux" "1hr_ActD"
)

HIC="${HIC_FILES[$SLURM_ARRAY_TASK_ID]}"
DS="${DATASETS[$SLURM_ARRAY_TASK_ID]}"
COND="${CONDITIONS[$SLURM_ARRAY_TASK_ID]}"
OUT="${OUTDIR}/${DS}/${COND}/mega"

echo "=== HiCCUPS KR: ${DS} / ${COND} ==="
echo "Input:  ${HIC}"
echo "Output: ${OUT}"
echo "Started: $(date)"

# ── Verify KR normalization is stored ────────────────────────────────────────
echo "Checking KR normalization availability..."
java -jar "${JUICER_TOOLS}" addNorm --check-only "${HIC}" 2>&1 | grep -i "KR" || {
    echo "WARNING: KR normalization may not be stored in ${HIC}"
    echo "Run: java -jar ${JUICER_TOOLS} addNorm -k KR ${HIC}"
}

mkdir -p "${OUT}" "${BASE}/logs/hiccups"

# ── HiCCUPS ───────────────────────────────────────────────────────────────────
java -Xmx60g -jar "${JUICER_TOOLS}" hiccups \
    --cpu \
    -k KR \
    -c "${CHROMS}" \
    -r 5000,10000,25000 \
    -f 0.1,0.1,0.1 \
    -p 4,2,1 \
    -i 7,5,3 \
    -t 0.02,1.5,1.75,2 \
    -d 20000,20000,50000 \
    "${HIC}" \
    "${OUT}"

EXIT=$?
echo "HiCCUPS exit code: ${EXIT}"

if [[ ${EXIT} -eq 0 && -f "${OUT}/merged_loops.bedpe" ]]; then
    N=$(grep -v "^#" "${OUT}/merged_loops.bedpe" | wc -l)
    echo "Loops called: ${N}"
else
    echo "ERROR: merged_loops.bedpe not found or job failed"
    exit 1
fi

echo "Finished: $(date)"
