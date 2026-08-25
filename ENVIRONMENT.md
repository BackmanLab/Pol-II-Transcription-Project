# Software Environment

R packages below are pinned to the exact versions in the `chiapet-hic-env`
conda/mamba environment used on Northwestern's Quest HPC for this project's
R analysis work.

## R and command-line tools

| Tool | Version(s) used | Notes |
|---|---|---|
| R | 4.4.3 (`chiapet-hic-env`, interactive analysis) | HPC submission scripts separately specify `module load R/4.3.0` — these are two different R installations; check the specific script for which applies |
| samtools | 1.6, 1.14 | Version varies by pipeline/script age — check the specific script |
| bedtools | 2.29.2, 2.30.0 | |
| bowtie2 | 2.5.4 | CUT&Tag alignment |
| bwa | 0.7.8, 0.7.17 | |
| hisat2 | 2.1.0 | RNA-seq/RIP-seq/EU-seq alignment |
| minimap2 | (unpinned) | Long-read alignment |
| rsem | 1.3.3 | |
| htseq | 2.0.2 | |
| TrimGalore | 0.6.10 | |
| fastqc | 0.12.0 | |
| deeptools | 3.1.1, 3.5.1, 3.5.6 | Version varies by script — check header |
| MACS2 | 2.2.9.1 | Peak calling |
| picard | 2.21.4 | Duplicate marking |
| seqtk | 1.4 | |
| pigz | (unpinned) | |
| Java | 1.8.0_131/162/191, 8.0.3.22, 11.0.10 | Multiple JDK versions across juicer_tools/HiCCUPS/Picard calls |
| CUDA | 7.0.28, 7.5.18, 8.0.44/54, 11.8.0 | HiCCUPS GPU loop calling |
| Juicer | see `HiC_analysis/README.md` | Excluded from repo — install separately (github.com/aidenlab/juicer) |

## R packages — pinned (confirmed in `chiapet-hic-env`)

Bioconductor:

| Package | Version |
|---|---|
| AnnotationDbi | 1.68.0 |
| BiocGenerics | 0.52.0 |
| BiocParallel | 1.40.0 |
| ChIPseeker | 1.42.0 |
| DESeq2 | 1.46.0 |
| GenomeInfoDb | 1.42.0 |
| GenomicAlignments | 1.42.0 |
| GenomicFeatures | 1.58.0 |
| GenomicRanges | 1.58.0 |
| Gviz | 1.50.0 |
| IRanges | 2.40.0 |
| InteractionSet | 1.34.0 |
| mariner | 1.6.0 |
| Rsamtools | 2.22.0 |
| Rsubread | 2.20.0 |
| TxDb.Hsapiens.UCSC.hg38.knownGene | 3.20.0 |
| biomaRt | 2.62.0 |
| org.Hs.eg.db | 3.20.0 |
| plyranges | 1.26.0 |
| rtracklayer | 1.66.0 |

CRAN:

| Package | Version | | Package | Version |
|---|---|---|---|---|
| RColorBrewer | 1.1_3 | | patchwork | 1.3.2 |
| cowplot | 1.2.0 | | pheatmap | 1.0.13 |
| data.table | 1.17.8 | | plyr | 1.8.9 |
| dplyr | 1.2.1 | | purrr | 1.2.2 |
| forcats | 1.0.1 | | readr | 2.2.0 |
| ggExtra | 0.11.0 | | readxl | 1.5.0 |
| ggdist | 3.3.3 | | rstatix | 0.7.3 |
| ggforce | 0.5.0 | | scales | 1.4.0 |
| ggnewscale | 0.5.2 | | strawr | 0.0.92 |
| ggplot2 | 4.0.2 | | stringr | 1.6.0 |
| ggplotify | 0.1.3 | | tibble | 3.3.1 |
| ggpubr | 0.6.3 | | tidyr | 1.3.2 |
| ggrepel | 0.9.8 | | tidyverse | 2.0.0 |
| ggridges | 0.5.7 | | viridis | 0.6.5 |
| ggside | 0.4.1 | | zoo | 1.8_15 |
| gridExtra | 2.3 | | here | 1.0.2 |
| matrixStats | 1.5.0 | | | |

Install Bioconductor packages via `BiocManager::install()`, CRAN packages
via `install.packages()`, or reproduce the full environment directly:
`mamba env create -f environment.yml` (see that file for the complete,
exported `chiapet-hic-env` specification).

## R packages — used in this repo but not in `chiapet-hic-env`

These packages are called by scripts in this repo but do not appear in the
`chiapet-hic-env` conda environment, meaning they were installed and used locally for analysis
(directly via `install.packages()`/`devtools::install_github()`)

| Package | Used in | Likely explanation |
|---|---|---|
| PupillometryR | `Imaging_Analysis.R`, `total_rnaseq_deseq2_analysis.R` | `Imaging_Analysis.R` |
| arrow | `EDU.K9_STORM.R` |
| beeswarm | `Imaging_Analysis.R` | 
| binom | `Imaging_Analysis.R` | 
| ggbeeswarm | `Imaging_Analysis.R` |
| hrbrthemes | `Imaging_Analysis.R` | 
| dendsort | `total_rnaseq_deseq2_analysis.R` | 
| rcartocolor | `total_rnaseq_deseq2_analysis.R` | 

## Notes

- Tool version spread (e.g., samtools 1.6 vs. 1.14, deepTools 3.1.1 vs.
  3.5.x) reflects differences across pipelines developed at
  different times, not an error — check the specific script's header/module
  load line rather than assuming one version applies repo-wide.
- The pinned R package table above reflects one snapshot of `chiapet-hic-env`
  (exported via `mamba env export --no-builds`) — it will drift from
  whatever's actually installed if the environment is updated later without
  re-exporting.
