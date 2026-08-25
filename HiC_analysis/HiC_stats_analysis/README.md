# Hi-C Stats / QC

| Script | Purpose | Figure |
|---|---|---|
| `scc.r` | HiCRep stratum-adjusted correlation coefficient (SCC), pairwise across conditions | **SI Fig. 5B** |
| `hic_stats.R` | General Hi-C QC: contacts relative to WT across downsampling groups, alignment fate/contact classification, library complexity, duplication rate | **SI Fig. 5A, C-F** |

`hic_stats.R` is written to be run interactively/section-by-section in
RStudio (no auto-save) — see its own header before scripting it end-to-end.
