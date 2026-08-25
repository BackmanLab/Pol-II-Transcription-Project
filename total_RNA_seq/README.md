# Total RNA-seq

Conditions: CTRL (0h), ActD (6h), POLR1A-degron (8h aux), POLR2A-degron (8h aux).

## Pipeline (`Total_RNAseq_HPC_processing/`)

1. `submit.sh` → `align_totalRNAseq.sh` — single-end trimming, HISAT2
   alignment, BAM filtering/sorting per sample.
2. `merge_totalrnaseq_bams.sh` — merges per-sample BAMs into per-condition
   merged BAMs (POLR1A: Rep1-4; POLR2A: Rep1-4; ActD: Rep1-3).

**This pipeline is also a dependency for EU-seq**: merged BAM read counts
here are used to compute EU-seq's normalization scale factors. Don't treat
this folder as standalone QC-only — `EU_seq_analysis/` reads its output.

## Downstream analysis (top level)

`total_rnaseq_deseq2_analysis.R` predates this paper (written for an earlier project) —
only a few plots from it were reused here:

| Script | Output | Figure |
|---|---|---|
| `total_rnaseq_deseq2_analysis.R` | ActD and POLR2A MA plots; sample correlation matrix | **Fig. 2G-H**; **SI Fig. 1I**; **SI Fig. 3A** |
| `functions.R` | Shared plotting utilities (`DE_Vol_Plot` volcano-plot generator, etc.) — sourced by other scripts, not standalone | — |
| `damage_markers.R` | Fold-change heatmaps + MA plot for three curated pathway gene panels (Core Splicing, Core Apoptosis, Core UPR), ActD/POLR2A vs. control | **SI Fig. 2A** (literal title match "Core Pathway Gene Expression Is Stable Upon Pol2 Depletion"; exact color scheme) |

`damage_markers.R` sources `functions.R` (`source(here::here("scripts/R/functions.R"))`)

