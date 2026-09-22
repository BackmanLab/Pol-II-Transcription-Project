# EU-seq Analysis

Nascent transcription (EU-seq) for WT and POLR2A-degron HCT116 cells, with a
NoEU control for background subtraction. Depends on `total_RNA_seq/` for
normalization scale factors (see that folder's README).

## Pipeline (`EU_seq_HPC_processing/`)

1. `genome/ref_hisat/build.sh`, `genome/rsem/build_reference.sh` — build
   references (run once).
2. `submit_EUseq.sh` → `align_euseq.sh` — per-sample QC, trimming,
   HISAT2 alignment, filtering, and sorting.
3. `BAM/strand.sh` — empirically confirms library strandedness (same method
   as `RIPseq_analysis`'s `strandedness.sh`; confirmed **sense**-stranded
   for this library — the opposite convention from RIP-seq's confirmed
   antisense result, and `euseq_subtract_count.sh`'s strand handling does
   **not** flip accordingly. Don't copy either pipeline's strand assumption
   into the other.).
4. `merge_euseq_bams.sh` — merges per-sample BAMs into per-condition merged BAMs.
5. `generate_bigwigs.sh` / `generate_pol2_vs_wt_bigwigs.sh` — coverage and
   log2(POLR2A/WT) bigWigs.
6. `scripts/euseq_subtract_count.sh` — NoEU background subtraction
   (base-pair level) and gene-level integer counts for DESeq2. Reverse-
   stranded, flags 83/163. Conditions: POLR2A (BM23-25) vs. WT (BM29-31),
   NoEU control: No_EU_ctrl_S49.
7. `scripts/euseq_dog_intron_count.sh` — intronic RPKM (gene body minus
   exons) and downstream-of-gene (DoG, 5kb past TES) RPKM, from the
   NoEU-subtracted bedgraphs.
8. `scripts/euseq_exon_rpkm.sh` — exonic RPKM computed directly from
   NoEU-subtracted bedgraphs (more accurate than deriving it as gene-body-
   minus-intron).
9. `scripts/euseq_termination.sh` — read-2 3' end termination position
   pileups within introns/exons.

## Downstream analysis (top level)

| Script | Purpose | Figure |
|---|---|---|
| `EU_seq_QC.R` | Pearson correlation heatmap + PCA (raw, unsubtracted featureCounts) | possibly a superseded/redundant QC check; the correlation matrix panels (**SI Fig. 1I**, **SI Fig. 3A**) come from `total_RNA_seq/total_rnaseq_deseq2_analysis.R` and `EU_seq_zscore_analysis.R` instead |
| `EU_seq_termination.R` | Plots termination pileups from `euseq_termination.sh` output | **SI Fig. 13E** (literal title "Read termination pileup within exons by exon length (control)") |
| `EU_Seq_DOG_Intron_analysis.R` | Intronic/DoG signal analysis (QC/distribution scatter, not a metagene profile) | no figure identified |
| `EU_seq_dog_bedgraph_tes.R` | DoG metagene around TES (−2kb to +10kb), readthrough | no figure identified |
| `EU_seq_log2ratio_coverage_analysis.R` | Metagene of log2(POLR2A/WT), stratified by E/I ratio, expression, and DEG status (Upstream/Body/Downstream/Intron features) | **Fig. 3C**; **SI Fig. 4A-C** (color match EU_log2="#E69F00", RNA_log2="#0072B2"; x_pos/bin_idx structure) |
| `EU_seq_zscore_analysis.R` | Nascent transcription analysis, FPKM from NoEU-subtracted bedgraphs; DEG histogram, DEG overlap, MA plot, PCA, E/I-tertile stratification (`ei_group`) | **Fig. 2I** (EU-seq MA plot); **Fig. 3E-F** (E/I tertile); **SI Fig. 3B-D, 3F**; **Fig. 3A-B**; **SI Fig. 3E** |
| `EU_seq_effective_ei.R` | Genomic vs. nascent effective E/I ratio, stratified by expression | no figure identified |
