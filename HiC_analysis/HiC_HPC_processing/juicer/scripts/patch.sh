#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Quest-specific adaptation of juicer.sh -- disables Aiden Lab's
# built-in cluster auto-detection (isRice/isBCM/isVoltron) on the
# working copy of juicer.sh, sets Quest module load commands, and
# patches the else block (queues, juiceDir, hardcoded usernames) so
# the script runs correctly on Quest instead of assuming a
# Rice/BCM/Broad cluster environment. Run once after installing
# Juicer, before the numbered pipeline.
#
# Figure output: none -- infrastructure/patch script.
##################################################################
JUICEDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer
SCRIPT=${JUICEDIR}/scripts/juicer.sh

# ── Disable Aiden Lab cluster detection ───────────────────────────────────
sed -i 's/^isRice=.*/isRice=0/' ${SCRIPT}
sed -i 's/^isBCM=.*/isBCM=0/' ${SCRIPT}
sed -i 's/^isVoltron=.*/isVoltron=0/' ${SCRIPT}

# ── Set Quest module load commands ────────────────────────────────────────
sed -i 's|^load_bwa=.*|load_bwa="module load bwa/0.7.17"|' ${SCRIPT}
sed -i 's|^load_java=.*|load_java="module load java/jdk11.0.10"|' ${SCRIPT}
sed -i 's|^load_samtools=.*|load_samtools="module load samtools/1.6"|' ${SCRIPT}
sed -i 's|^load_gpu=.*|load_gpu=""|' ${SCRIPT}
sed -i 's|^load_awk=.*|load_awk=""|' ${SCRIPT}

# ── Patch the else block: isVoltron=1 → 0, set Quest paths and queues ────
# The else block resets isVoltron=1 — find and patch it
sed -i "/^else$/,/^fi$/ s/isVoltron=1/isVoltron=0/" ${SCRIPT}
sed -i "/^else$/,/^fi$/ s|juiceDir=.*|juiceDir=\"${JUICEDIR}\"|" ${SCRIPT}
sed -i "/^else$/,/^fi$/ s/queue=\".*\"/queue=\"genomics\"/" ${SCRIPT}
sed -i "/^else$/,/^fi$/ s/queue_time=\".*\"/queue_time=\"2880\"/" ${SCRIPT}
sed -i "/^else$/,/^fi$/ s/long_queue=\".*\"/long_queue=\"genomics\"/" ${SCRIPT}
sed -i "/^else$/,/^fi$/ s/long_queue_time=\".*\"/long_queue_time=\"2880\"/" ${SCRIPT}

# ── Replace any remaining genomics-himem or genomics-gpu with genomics ───
sed -i 's/genomics-himem/genomics/g' ${SCRIPT}
sed -i 's/genomics-gpu/genomics/g' ${SCRIPT}

# ── Fix any hardcoded usernames ───────────────────────────────────────────
sed -i 's/squeue -u lmc0633/squeue -u $USER/g' ${SCRIPT}
sed -i 's/squeue -u lmc0633/squeue -u $USER/g' \
    ${JUICEDIR}/scripts/split_rmdups_sam.awk

# ── Verify ────────────────────────────────────────────────────────────────
echo "=== Key variable settings ==="
grep "^isRice\|^isBCM\|^isVoltron\|^load_bwa\|^load_java\|^load_samtools\|^load_gpu\|^load_awk" \
    ${SCRIPT}

echo "=== Queue settings in else block ==="
awk '/^else$/,/^fi$/{print NR": "$0}' ${SCRIPT} | grep "queue\|juiceDir\|isVoltron"

echo "=== Remaining non-genomics partition references ==="
grep "genomics-himem\|genomics-gpu\|lmc0633" ${SCRIPT} \
    ${JUICEDIR}/scripts/split_rmdups_sam.awk 2>/dev/null \
    && echo "WARNING: manual cleanup needed above" \
    || echo "Clean — no non-standard partitions or hardcoded usernames found"

echo "=== juicer_tools call passes -v user correctly ==="
grep "split_rmdups_sam.awk" ${SCRIPT}
