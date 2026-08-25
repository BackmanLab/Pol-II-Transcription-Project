# RIP-seq Analysis

RIP-seq for POLR2A (P2AB) and EZH2, with IgG and Input controls, in WT and
POLR2A-KO (P1KO) HCT116 backgrounds. `RIPseq_HPC_processing/` handles
alignment through coverage/RPKM generation; the top-level `.R` scripts are
the downstream feature-level and figure-generating analyses.

## Pipeline (`RIPseq_HPC_processing/`)

Run roughly in this order:

1. **`genome/ref_hisat/build.sh`** — builds the HISAT2 genome index (run once).
2. **`fastq/submit.sh`** → **`align.sh`** — per-sample QC (FastQC), trimming
   (TrimGalore), HISAT2 alignment, samtools filtering/sorting, and bigWig
   generation (CPM + RPGC).
3. **`strandedness.sh`** — empirically confirms library strandedness by
   comparing read1 alignment strand against a set of well-annotated,
   non-overlapping reference genes. Confirmed reverse-stranded for this
   library prep (dUTP second-strand method); this result is load-bearing
   for the strand-flip logic in `rpkm.sh` — re-run if the library kit changes.
4. **`BAM/submit.sh`** → **`BAM/merge_bam_bpm_bigwig.sh`** — merges replicate
   BAMs (POLR2A: Rep2+Rep3; EZH2: Rep1+Rep2) and generates BPM-normalized
   bigWigs matching the parameters that reproduce the expected TES signal
   profile (binSize 20, smoothLength 60).
5. **`rpkm.sh`** — per-replicate exon/intron RPKM matrices (strand-specific,
   RPM-normalized bedGraphs intersected with exon/intron BED intervals). No
   IP/Input subtraction or replicate averaging here — that's left to the R
   analysis scripts below. Exon/intron BED intervals are reused from the
   EU-seq pipeline (same GTF, so E/I interval definitions stay consistent
   across both).
6. **`qc.sh`** / **`rip_qc.R`** — mapping/duplication rate, replicate
   correlation, IP enrichment over Input, and genomic feature distribution
   of signal. PCR Bottleneck Coefficient (PBC1) is deliberately excluded —
   see `qc.sh`'s docstring: RIP/CLIP-seq reads legitimately stack at exons,
   splice junctions, and TSSs, so positional duplication doesn't separate
   PCR artifact from real signal the way it does for ChIP-seq. Peak-calling-
   based FRiP was tried (Piranha/MACS2) and found uninformative for this
   IP-vs-multi-polymerase-Input comparison; the genomic feature distribution
   is used instead.

## Downstream analysis (top level)

| Script | Output | Figure |
|---|---|---|
| `rip_qc.R` | Replicate correlation, IP enrichment, feature distribution | **SI Fig. 12** |
| `rip_seq_logfc_analysis.R` | log2(P2AB/Input) by feature type, faceted by expression and length | **Fig. 10D** |
| `rip_seq_intron_dog_heatmap.R` | TES and intron-termination heatmaps (PRO-seq style) | **Fig. 10E-F** |
| `rip_seq_effective_ei.R` | Exon/intron RPKM, log2(P2AB/WT_Input) and log2(P2AB/P2KO_Input) | **Fig. 10G** |

`rip_seq_intron_dog_heatmap.R` and `rip_seq_logfc_analysis.R` independently
build the same five feature-type BED interval sets (Promoter, Enhancer,
Exon, Intron, TES) as a shared first step before diverging — this is real
duplication, not a copy-paste header error. If editing that shared
BED-building logic, update both files.

## Notes on RPKM strand handling

`rpkm.sh`'s strand assignment is **flipped** relative to a naive
read1-strand-equals-gene-strand approach (read1 on "−" builds "+"-strand
gene coverage, and vice versa) — this is the opposite of the EU-seq
pipeline's own strand handling, but both are independently confirmed
correct for their respective library preps via `strandedness.sh`
(RIP-seq: 99.7% read1-antisense-to-gene; EU-seq: 98.7% read1-sense-to-gene).
