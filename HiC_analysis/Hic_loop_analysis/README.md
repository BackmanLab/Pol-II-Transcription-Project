# Hi-C Loop Analysis

Loop calling itself happens upstream (`HiC_HPC_processing/juicer/scripts/run_hiccups.sh`).
Everything here is downstream: APA, differential loop classification, anchor
composition, and locus-level visualization.

| Script | Purpose | Figure |
|---|---|---|
| `final_APA_analysis.R` | Aggregate Peak Analysis — mariner, NORM="NONE" (raw counts; matches mariner's own documented convention of extracting raw contacts and normalizing downstream, rather than pre-balancing the matrix), 15×15 bins (±35kb), `ds_wt_pol2` | **SI Fig. 8B-D**; **Fig. 5B, 5E** |
| `hic_loop_lfc_analysis.R` | Differential loop classification — Gained/Stable/Lost by \|log2FC\| > 0.5, contacts extracted via strawr at fixed loop positions in both conditions | **Fig. 5A, 5C, 5D**; **SI Fig. 8A, 8F** |
| `EP_analysis_functions.R` | **Shared functions, not standalone** — loop anchor cCRE classification (Promoter/Proximal-Enhancer/Distal-Enhancer/Open-Chromatin priority hierarchy). Sourced by other scripts in this folder and reused (adapted) by `RIPseq_analysis/rip_seq_intron_dog_heatmap.R` / `rip_seq_logfc_analysis.R` for their own BED-building step. | — |
| `FINAL_EP_loop_analysis_hiccups.R` | Loop anchor cCRE composition; pure CTCF-CTCF loop classification (no Promoter/Enhancer co-occupancy) | **SI Fig. 8E** (non-gene anchor classification); **SI Fig. 9E** (pure CTCF Gained/Stable/Lost) |
| `final_loop_scatter_analysis.R` | Loop length vs. contact strength scatter, by CTCF pair type | **SI Fig. 9B-D**; **SI Fig. 9A** |
| `hic_general_analysis.R` | Loop strength vs. DEG status (MA-style plot); also contains the anchor-type classification/raincloud logic extracted into `hic_loop_strength_by_anchor_type.R` below | no figure identified for the MA-plot portion specifically |
| `hic_loop_strength_by_anchor_type.R` | Loop strength (SCALE-normalized contacts) by anchor type (Exon-Exon/Exon-cCRE/Exon-Intron/cCRE-cCRE/Intron-Intron), Wildtype Hi-C only | **Fig. 5G** — extracted from the already-active code path in `hic_general_analysis.R`; see `chia_loop_analysis/README.md` for the paired **Fig. 5F** script |
| `hic_loop_track_plotter.R` | Locus-level contact pyramids (GENOVA-style) + gene/cCRE tracks, ported from the ChIA-PET locus plotter (`chia_loop_analysis/updated_CHIA_PET_loci_plotter.R`) | **Fig. 6A-B** (RAD50/MMP7 locus tracks) |

All loop-analysis scripts share a `DS_GROUP` switch to select which
downsampling group (`ds_wt_pol2`, `ds_all`, `unsampled`) to load — check this
is set correctly before running, since it silently determines both the input
`.hic` files and the loop list used.
