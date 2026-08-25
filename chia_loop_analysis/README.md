# ChIA-PET Analysis

| Script | Purpose | Figure |
|---|---|---|
| `updated_CHIA_PET_loci_plotter.R` | Self-contained locus plotter: ChIA-PET arcs, gene annotation + cCREs, POLR2A CUT&Tag signal, Total RNA-seq, EU-seq (5 tracks) | **Fig. 10A-B** |
| `chiapet_loop_strength_by_anchor_type.R` | Loop strength (PET count) by anchor type (Exon-Exon/Exon-cCRE/Exon-Intron/cCRE-cCRE/Intron-Intron), public POLR2A ChIA-PET (ENCSR035PVZ) | **Fig. 5F** |

`chiapet_loop_strength_by_anchor_type.R` is the direct ChIA-PET counterpart
to `HiC_analysis/Hic_loop_analysis/hic_loop_strength_by_anchor_type.R`
(same anchor-type categories, same color palette, same plot style) — the
two are meant to be read side by side as **Fig. 5F-G**.

