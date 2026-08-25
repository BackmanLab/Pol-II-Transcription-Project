#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Runs mega.sh to merge replicates per condition for the
# 3-condition depth-matched dataset (ds_all: WT, Pol2, ActD --
# depth-matched to ActD).
#
# Required per replicate before running:
#   juicer_ds_all/<COND>_Rep<N>/aligned/merged1.txt
#   juicer_ds_all/<COND>_Rep<N>/aligned/inter.txt
#
# Output per condition:
#   mega_ds_all/<COND>/mega/aligned/inter_30.hic
#
# IMPORTANT: Run from LOGIN NODE.
# Usage: bash run_mega_all.sh

BASE=/projects/b1042/BackmanLab/Lucas/060825_hic-samp
JUICER_SCRIPT=${BASE}/juicer/scripts/mega.sh
JUICER_DIR=${BASE}/juicer
SITE_FILE=${JUICER_DIR}/restriction_sites/hg38_DpnII.txt

JUICER_OUT=${BASE}/juicer_ds_all
MEGA_BASE=${BASE}/mega_ds_all

CONDITIONS=(WT_HCT116_CTRL Pol2_6hrs_Aux 1hr_ActD)

[[ -f "$JUICER_SCRIPT" ]] || { echo "ERROR: mega.sh not found at $JUICER_SCRIPT" >&2; exit 1; }

for COND in "${CONDITIONS[@]}"; do
    echo ""
    echo "═══════════════════════════════════════════════"
    echo "Setting up mega topDir for: ${COND}"

    TOPDIR=${MEGA_BASE}/${COND}

    for REP in 1 2; do
        SAMPLE="${COND}_Rep${REP}"
        SRC_ALIGNED=${JUICER_OUT}/${SAMPLE}/aligned

        for REQFILE in merged1.txt inter.txt; do
            [[ -f "${SRC_ALIGNED}/${REQFILE}" ]] || {
                echo "  ERROR: ${SRC_ALIGNED}/${REQFILE} not found" >&2
                exit 1
            }
        done

        LINK_DIR=${TOPDIR}/${SAMPLE}/aligned
        mkdir -p ${LINK_DIR}

        for FILE in $(ls ${SRC_ALIGNED}/); do
            [[ -f "${SRC_ALIGNED}/${FILE}" ]] && \
                ln -sf ${SRC_ALIGNED}/${FILE} ${LINK_DIR}/${FILE}
        done

        echo "  ✓ Linked ${SAMPLE}/aligned/"
    done

    echo "  Running mega.sh for ${COND}..."
    bash ${JUICER_SCRIPT} \
        -d ${TOPDIR} \
        -g ${JUICER_DIR}/references/hg38.chrom.sizes \
        -s DpnII \
        -y ${SITE_FILE} \
        -D ${JUICER_DIR} \
        -q genhimem \
        -l genhimem \
        -Q 2880 \
        -L 2880 \
        -A p32171 &

    echo "  Jobs submitted for ${COND} (pid $!)"
done

wait
echo ""
echo "All mega jobs submitted for ds_all."
echo "Output: ${MEGA_BASE}/<COND>/mega/aligned/inter_30.hic"
