#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=strand
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Empirically determines library strandedness by comparing read1
# alignment strand against the annotated strand of a handful of
# well-annotated, non-overlapping, robustly-expressed reference
# genes.
#
# Method: for each candidate gene, count read1 alignments that
# land on the SAME strand as the gene (sense) vs the OPPOSITE
# strand (antisense), pooled across genes. A stranded library
# will show a strong majority (~85%+) in one direction; an
# unstranded library will be close to 50/50.
#
# Reuses hg38_genes.bed from the EU-seq pipeline (090124_euseq)
# so gene intervals/strand calls match that pipeline exactly —
# this is a deliberate cross-project reuse, not a stray path.
#
# Usage:
#   sbatch check_strandedness.sh
##################################################################

set -euo pipefail
module purge
module load samtools/1.6

## ── Paths — hardcoded for the 022025_RIPseq project ──────────────────────────
BAM="/projects/b1042/BackmanLab/Lucas/022025_RIPseq/BAM/BM29_S46_filt.sort.bam"
GENE_BED="/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline/hg38_genes.bed"

MIN_GENE_LEN=3000    # bp — lowered from 10000; classic housekeeping genes
                     # (GAPDH, ACTB, etc.) are compact and were all being
                     # filtered out at the old threshold
MIN_READS=50         # minimum read1 reads in-region to trust a gene's contribution

if [[ ! -s "${GENE_BED}" ]]; then
  echo "ERROR: Gene BED not found: ${GENE_BED}" >&2
  exit 1
fi
if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
  echo "ERROR: BAM index not found for ${BAM} — run 'samtools index ${BAM}' first" >&2
  exit 1
fi

echo "=========================================================="
echo "Strandedness check via read1-vs-gene-strand concordance"
echo "BAM:      ${BAM}"
echo "Gene BED: ${GENE_BED}"
echo "=========================================================="

## Candidate reference genes — bare Ensembl gene IDs. Matched against
## hg38_genes.bed allowing for a version suffix, since the GTF-derived
## BED may retain versioned IDs. Mix of compact classics (short, very
## highly expressed) and longer standard qPCR/RNA-seq reference genes
## (more introns, more total sequence for the strand call).
CANDIDATE_IDS=(
  "ENSG00000111640"  # GAPDH   ~4.9kb
  "ENSG00000075624"  # ACTB    ~3.7kb (often skipped — overlaps neighboring genes)
  "ENSG00000196230"  # TUBB    ~8.1kb
  "ENSG00000142541"  # RPL13A  ~5.5kb
  "ENSG00000167658"  # EEF2    ~9.4kb
  "ENSG00000166710"  # B2M     ~7.6kb
  "ENSG00000089157"  # RPLP0   ~4.5kb
  "ENSG00000112592"  # TBP     ~9.1kb
  "ENSG00000023330"  # ALAS1   ~16kb
  "ENSG00000102144"  # PGK1    ~20kb
  "ENSG00000164924"  # YWHAZ   ~30kb
  "ENSG00000073578"  # SDHA    ~38kb
  "ENSG00000165704"  # HPRT1   ~44kb
)

TOTAL_ANTISENSE=0   # read1 opposite strand from gene (matches EU-seq's rf/reverse convention)
TOTAL_SENSE=0       # read1 same strand as gene       (fr/forward convention)
N_GENES_USED=0

for gid in "${CANDIDATE_IDS[@]}"; do
  gene_line=$(awk -v gid="${gid}" '$4==gid || $4 ~ "^"gid"\\."' "${GENE_BED}" | head -1)

  if [[ -z "${gene_line}" ]]; then
    echo "  [skip] ${gid}: not found in gene BED"
    continue
  fi

  chrom=$(echo  "${gene_line}" | cut -f1)
  start=$(echo  "${gene_line}" | cut -f2)
  end=$(echo    "${gene_line}" | cut -f3)
  strand=$(echo "${gene_line}" | cut -f6)
  gene_len=$((end - start))

  if (( gene_len < MIN_GENE_LEN )); then
    echo "  [skip] ${gid}: too short (${gene_len}bp < ${MIN_GENE_LEN}bp)"
    continue
  fi

  ## Require zero overlap with any OTHER annotated gene — an overlapping
  ## antisense neighbor would contaminate the strand call
  overlap_count=$(awk -v c="${chrom}" -v s="${start}" -v e="${end}" -v gid="${gid}" \
    '$1==c && $4!=gid && !($3<=s || $2>=e)' "${GENE_BED}" | wc -l)

  if (( overlap_count > 0 )); then
    echo "  [skip] ${gid}: overlaps ${overlap_count} other annotated gene(s)"
    continue
  fi

  region="${chrom}:$((start+1))-${end}"

  n_fwd=$(samtools view -c -f 64 -F 16 "${BAM}" "${region}")   # read1, mapped to + strand
  n_rev=$(samtools view -c -f 80        "${BAM}" "${region}")   # read1, mapped to - strand (64+16)
  n_total=$((n_fwd + n_rev))

  if (( n_total < MIN_READS )); then
    echo "  [skip] ${gid}: too few read1 reads in region (${n_total})"
    continue
  fi

  if [[ "${strand}" == "+" ]]; then
    antisense=${n_rev}; sense=${n_fwd}
  else
    antisense=${n_fwd}; sense=${n_rev}
  fi

  echo "  [use]  ${gid}  strand=${strand}  len=${gene_len}bp  read1_fwd=${n_fwd}  read1_rev=${n_rev}  n=${n_total}"

  TOTAL_ANTISENSE=$((TOTAL_ANTISENSE + antisense))
  TOTAL_SENSE=$((TOTAL_SENSE + sense))
  N_GENES_USED=$((N_GENES_USED + 1))
done

echo ""
echo "=========================================================="

if (( N_GENES_USED == 0 )); then
  echo "ERROR: no candidate genes passed filtering — cannot determine strandedness." >&2
  echo "Check that gene_id format in ${GENE_BED} matches Ensembl gene IDs." >&2
  exit 1
fi

TOTAL=$((TOTAL_ANTISENSE + TOTAL_SENSE))
PCT_ANTI=$(awk -v a="${TOTAL_ANTISENSE}" -v t="${TOTAL}" 'BEGIN{printf "%.1f", 100*a/t}')
PCT_SENSE=$(awk -v a="${TOTAL_SENSE}"    -v t="${TOTAL}" 'BEGIN{printf "%.1f", 100*a/t}')

echo "Genes used:                ${N_GENES_USED}"
echo "Total read1 reads pooled:  ${TOTAL}"
echo ""
echo "read1 ANTISENSE to gene (reverse-stranded, matches EU-seq convention): ${TOTAL_ANTISENSE}  (${PCT_ANTI}%)"
echo "read1 SENSE to gene     (forward-stranded):                           ${TOTAL_SENSE}  (${PCT_SENSE}%)"
echo ""

if (( $(awk -v p="${PCT_ANTI}"  'BEGIN{print (p>=85)}') )); then
  echo "RESULT: Reverse-stranded library (read1 = antisense to gene)."
  echo "        Matches the EU-seq pipeline's convention — same strand-"
  echo "        assignment logic can be reused as-is."
elif (( $(awk -v p="${PCT_SENSE}" 'BEGIN{print (p>=85)}') )); then
  echo "RESULT: Forward-stranded library (read1 = sense to gene)."
  echo "        OPPOSITE the EU-seq pipeline's convention — strand"
  echo "        assignment must be flipped (swap +/- logic) for this data."
else
  echo "RESULT: No strong majority (~50/50) — library appears UNSTRANDED."
  echo "        Strand-specific counting isn't appropriate; pool both"
  echo "        strands together per interval instead."
fi
echo "=========================================================="
