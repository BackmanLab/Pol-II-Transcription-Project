#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Merges per-replicate Juicer outputs into a mega .hic
# per condition. Uses Juicer's mega.sh from the existing installation.
#
# IMPORTANT: Run from the LOGIN NODE after ALL 02_run_juicer.sh jobs complete.
#
# Usage: bash 03_merge_mega.sh

JUICER_DIR=/projects/b1042/BackmanLab/juicer
MEGA_SCRIPT=${JUICER_DIR}/scripts/mega.sh
WORKDIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer_ds
MEGADIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/mega_hic
GENOME_ID=hg38

CONDITIONS=(WT_HCT116_CTRL Pol2_6hrs_Aux 1hr_ActD)

# ── Verify mega.sh exists ─────────────────────────────────────────────────
if [[ ! -f "$MEGA_SCRIPT" ]]; then
    echo "mega.sh not found at $MEGA_SCRIPT" >&2
    echo "Checking for alternative location..."
    MEGA_SCRIPT=$(find /projects/b1042/BackmanLab/juicer -name "mega.sh" 2>/dev/null | head -1)
    if [[ -z "$MEGA_SCRIPT" ]]; then
        echo "ERROR: mega.sh not found in Juicer installation." >&2
        echo "Check: ls /projects/b1042/BackmanLab/juicer/scripts/" >&2
        exit 1
    fi
    echo "Found: $MEGA_SCRIPT"
fi

# ── Verify per-rep outputs ────────────────────────────────────────────────
echo "Verifying per-rep outputs..."
ALL_OK=true
for COND in "${CONDITIONS[@]}"; do
    for REP in 1 2; do
        HIC=${WORKDIR}/${COND}_Rep${REP}/aligned/inter_30.hic
        NODUPS=${WORKDIR}/${COND}_Rep${REP}/aligned/merged_nodups.txt
        if [[ ! -f "$HIC" ]]; then
            echo "  MISSING: $HIC" >&2; ALL_OK=false
        else
            echo "  OK: $(du -sh $HIC | cut -f1)  $HIC"
        fi
        if [[ ! -f "$NODUPS" ]]; then
            echo "  MISSING: $NODUPS" >&2; ALL_OK=false
        fi
    done
done

if [[ "$ALL_OK" != "true" ]]; then
    echo "Some outputs missing — check Juicer job logs before merging." >&2; exit 1
fi

# ── Run mega.sh for each condition ────────────────────────────────────────
for COND in "${CONDITIONS[@]}"; do

    echo ""
    echo "═══════════════════════════════════════════════════"
    echo "Merging: ${COND}"

    MEGA_TOP=${MEGADIR}/${COND}
    mkdir -p ${MEGA_TOP}

    # mega.sh expects subdirectories under topDir each containing aligned/
    for REP in 1 2; do
        ln -sfn ${WORKDIR}/${COND}_Rep${REP} ${MEGA_TOP}/Rep${REP}
    done

    # Run mega.sh — submits SLURM jobs internally, run from login node
    bash ${MEGA_SCRIPT} \
        -d ${MEGA_TOP} \
        -D ${JUICER_DIR} \
        -g ${GENOME_ID}

    echo "  Mega jobs submitted for ${COND}"
done

echo ""
echo "══════════════════════════════════════════════════════════"
echo "All mega merge jobs submitted."
echo ""
echo "Final outputs (update APA pipeline paths to these):"
for COND in "${CONDITIONS[@]}"; do
    echo "  ${MEGADIR}/${COND}/aligned/inter_30.hic"
done
echo ""
echo "In the APA pipeline, set all DEPTH_SCALE values to 1.0:"
echo "  DEPTH_SCALE <- list(wt=1.0, pol2=1.0, actd=1.0)"
