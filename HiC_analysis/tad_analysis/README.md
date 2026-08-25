# TAD / Insulation Analysis

| Script | Purpose | Figure |
|---|---|---|
| `insulation_analysis.R` | 1D insulation score (25kb, 250kb window, SCALE-normalized, z-scored per chromosome), 4-way split: WT/POLR2A-degraded × CTCF-anchored/non-CTCF TADs, expression-stratified (Low/Medium/High) | **Fig. 4A-B**; **SI Fig. 7A-F** (non-CTCF and CTCF-bounded insulation/ATA) |
| `ata_analysis.R` | Aggregate TAD Analysis (GENOVA::ATA, 100×100 bins, 25kb, SCALE-normalized), CTCF classification, gene-linked TPM tertiles | **SI Fig. 6A-B** |

Both use GENOVA/strawr-based contact loading against the same `DS_GROUP`
downsampled Hi-C files as the rest of the Hi-C pipeline.
