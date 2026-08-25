# Contact Scaling Analysis

| Script | Purpose | Figure |
|---|---|---|
| `scaling_analysis.R` | Gene-level contact scaling (S): genes stratified by E/I ratio (100 bins → tertiles for two-group comparisons) or expression tertile; per-bin contact probability P(s) pooled across genes, S = log10(P(s)) vs. log10(distance) regression slope. Sliding-window S(d) by TPM quartile for the distance-resolved panel. | **Fig. 4F-H** |
| `tandem_gene_analysis.R` | Diagnostic only — genome-wide inter-gene gap distribution, used to justify the 50kb `MAX_GAP` tandem-cluster threshold used by `scaling_analysis.R`'s short-gene (<100kb) cluster branch (exon/intergenic-ratio analysis, separate from the long-gene exon/intron branch used for **Fig. 4F-H**). Not a figure-generating script itself. | — |

See `scaling_analysis.R`'s own header for the full breakdown of the
≥100kb single-gene (exon/intron) vs. <100kb tandem-cluster (exon/intergenic)
branches — these compute different ratios on different gene
populations, matching the two non-coding-packing mechanisms described in the
paper (introns within long genes vs. intergenic space between short genes).
