#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH --job-name=microc_dl
#SBATCH --output=/projects/b1042/BackmanLab/Lucas/microc_encode/logs/microc_dl_%j.out
#SBATCH --error=/projects/b1042/BackmanLab/Lucas/microc_encode/logs/microc_dl_%j.err

##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Downloads the ENCODE POLR2A-degron HCT116 Micro-C contact maps (Table 5):
#   ENCFF534TYD  untreated (ENCSR533VUS)
#   ENCFF302CUI  6 h 5-Ph-IAA (ENCSR218URU)
# Each download is checked against the expected file size and the .hic header.
##################################################################

DIR=/projects/b1042/BackmanLab/Lucas/microc_encode

download_hic() {
    ACC=$1
    EXPECT=$2
    curl -sLf -C - "https://www.encodeproject.org/files/${ACC}/@@download/${ACC}.hic" -o ${DIR}/${ACC}.hic
    echo "${ACC}: curl exit code $?"
    GOT=$(stat -c %s ${DIR}/${ACC}.hic)
    MAGIC=$(head -c 3 ${DIR}/${ACC}.hic)
    if [ "${GOT}" = "${EXPECT}" ] && [ "${MAGIC}" = "HIC" ]; then
        echo "${ACC}: OK (${GOT} bytes, header ${MAGIC})"
    else
        echo "${ACC}: PROBLEM got ${GOT} bytes (expected ${EXPECT}), header '${MAGIC}'"
    fi
}

download_hic ENCFF534TYD 27301466053
download_hic ENCFF302CUI 33785386023
