# ATAC-seq Analysis

Downstream analysis of public ENCODE ATAC-seq data (POLR2A-degron HCT116,
+auxin 6h vs. untreated — see Table 5 for accessions). No alignment/peak
calling here; both scripts consume bigWig and narrowPeak files directly.

| Script | Purpose | Figure |
|---|---|---|
| `atac_analysis.R` | Gene body / promoter / peak-level accessibility, stratified by expression quartile, genomic E/I ratio, and DEG status (matched total RNA-seq) | **SI Fig. 11E** (accessibility vs. neighborhood transcription by E/I group) |
| `atac_gene_density.R` | Exon/intron accessibility vs. local gene density (gene-count density, expression density, and genomic E/I ratio — each as Q1/Q4 split), POLR2A-KO delta tracks, continuous dose-response scan, correlation against EU-seq nascent transcription change | **SI Fig. 11D** (mean enrichment + log2FC by gene-count density) |

Both scripts use MACS2 fold-enrichment-over-background bigWig values as the
accessibility signal, and resolve narrowPeak duplicates (merged replicates)
by keeping the highest fold-enrichment value per unique coordinate.
