#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: One-off resubmission of Pol2_6hrs_Aux Rep2 in
# juicer_ds_wt_pol2 -- node failure on qnode3111 killed chunk
# fastq012 alignment.

JUICER_SCRIPT=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer/scripts/juicer.sh
JUICER_DIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer
GENOME_FASTA=${JUICER_DIR}/references/hg38.fa
CHROM_SIZES=${JUICER_DIR}/references/hg38.chrom.sizes
SITE_FILE=${JUICER_DIR}/restriction_sites/hg38_DpnII.txt

SAMPLE_DIR=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer_ds_wt_pol2/Pol2_6hrs_Aux_Rep2

# Wipe just this sample's output
rm -rf ${SAMPLE_DIR}/aligned ${SAMPLE_DIR}/splits \
       ${SAMPLE_DIR}/debug ${SAMPLE_DIR}/HIC_tmp

# Relink FASTQs
FASTQ_DIR=${SAMPLE_DIR}/fastq
mkdir -p ${FASTQ_DIR}
FASTQ_DS=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/downsampled_wt_pol2
for F in ${FASTQ_DS}/Pol2_6hrs_Aux/Rep2/*_ds.fastq.gz; do
    ln -sf $F ${FASTQ_DIR}/$(basename $F)
done

echo "Resubmitting Pol2_6hrs_Aux_Rep2 (wt_pol2)..."
bash ${JUICER_SCRIPT} \
    -d ${SAMPLE_DIR} \
    -g hg38 \
    -z ${GENOME_FASTA} \
    -p ${CHROM_SIZES} \
    -s DpnII \
    -y ${SITE_FILE} \
    -A b1042

echo "Done."
