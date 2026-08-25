# Imaging Analysis

PWS, EU-STORM, EdU/K9-STORM, and ChromSTEM analysis scripts.

| Script | Purpose | Figure |
|---|---|---|
| `Imaging_Analysis.R` | Large multi-panel script — PWS metrics (D, diffusion, FMM), widefield EU signal, CTCF degron validation, RNAP I/K9 clustering, ChromSTEM Domain Radius/CVC | **Fig. 1C-F**; **Fig. 2B** (`ggtitle("EU Signal Intensity - Widefield")`); **Fig. 9A-C, 9E**; **SI Fig. 1F** (`ggtitle("Pol 2 CTCF")`); **SI Fig. 11F-G** ("Plot quadrants for CVC versus size"); **Fig. 7D-E**; **SI Fig. 1A-E, 1H** |
| `EUstorm.R` | EU-STORM nucleus analysis, POLR2A (Auxin) vs. WT (DMSO): blink density, cluster density, cluster size — boxplot, raincloud, ECDF, and PDF per metric, plus a stats table (ROI- and cell-level, with effect sizes) | **Fig. 2E-F** (cluster size and blinks/µm² cumulative-fraction panels match exactly) |
| `EDU.K9_STORM.R` | STORM domain classification scatter (density vs. domain radius), WT vs. treatment, quadrant thresholds from WT median; both EdU-STORM and K9-STORM, raw and log10-transformed versions | **Fig. 8F**; **SI Fig. 10D-E**; **Fig. 8C-D** |
| `chrom_stem_radial_density.R` | ChromEM radial density, Wildtype vs. POLR2A-degraded — 4-panel figure (Nascent/Small Mature/Decaying/Large Mature domains), density vs. distance | **Fig. 9F** |

## Notes

- `Imaging_Analysis.R` mixes hardcoded personal macOS paths
  (`/Users/lucascarter/Documents/...`) with Quest paths — needs
  genericizing before it's runnable by anyone else. Given its size (1600+
  lines, many distinct panels), this is a strong candidate for splitting
  into separate per-figure scripts if it gets touched again.
- `EUstorm.R`'s own header documents that `col_Data()` (referenced but not
  included in what was originally provided) was reimplemented as
  `read_eu_storm_file()`, based on the transposed structure of the actual
  input CSVs — worth reading that header in full before modifying the
  parsing logic.
