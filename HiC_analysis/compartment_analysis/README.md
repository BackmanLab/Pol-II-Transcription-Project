# A/B Compartment Analysis

Eigenvectors are computed upstream
(`HiC_HPC_processing/juicer/scripts/run_compartments.sh` — Juicer eigenvector
tool, SCALE-normalized, 500kb, sign-oriented by GC content). Both scripts
here consume those bedgraphs.

| Script | Purpose | Figure |
|---|---|---|
| `saddle_analysis.R` | Saddle plots (SCALE-normalized O/E, 500kb, 20 quantile bins), compartment strength CS = (AA+BB)/(AB+BA), coding-density stratification | **Fig. 4C-E**; **SI Fig. 7G-H** (log2 coding/noncoding distribution, per-corner O/E) |
| `compartment_analysis.R` | Pipeline 2 — Series 1: E/I boxplot marginals by coding-density tertile; Series 2: PC1 bar marginals by coding-density tertile; expression-stratified PC1 marginals | **SI Fig. 6C-E** |

Both scripts use the same `DS_GROUP` downsampling switch as the loop-analysis
scripts — confirm it's set to the intended group before running.
