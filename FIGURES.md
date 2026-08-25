# Figure Index

One row per individual panel, tying it to the script that produces it (or its
non-script source: collaborator-provided, raw image, schematic, IGV screenshot).

| Figure | Panel | Description | Source |
|---|---|---|---|
| Figure 1 | A | Self-Avoiding Walk model: random walk, no forced returns, no domains | Provided by external collaborator |
| Figure 1 | B | Self-Returning-Excluded Volume (SR-EV) model: forced returns form packing domains | Provided by external collaborator |
| Figure 1 | C | PWS following ActD: average nuclear packing scaling of packing domains | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 1 | D | ChromSTEM following treatment: packing scaling of individual domains | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 1 | E | ChromSTEM following treatment: domain radius of individual domains | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 1 | F | ChromSTEM following treatment: chromatin volume concentration of individual domains | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 1 | G | ChromSTEM tomograms: untreated chromatin state vs. domain expansion with treatment | Provided by in-lab collaborator |
| Figure 1 | H | Representative domain radial density from center to periphery (CVC) | Provided by in-lab collaborator |
| Figure 1 | I | Representative map of D_nucleus chromatin scaling values (wildtype vs. ActD) | Provided by in-lab collaborator |
| Figure 2 | A | Representative widefield fluorescence image per condition (DNA/POLR2A/nascent RNA) | raw microscopy image |
| Figure 2 | B | Cell total corrected fluorescence, widefield, by condition | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 2 | C | Representative SMLM images of EU-labeled RNA per condition | raw microscopy image |
| Figure 2 | D | Representative SMLM images of EU-labeled RNA per condition | raw microscopy image |
| Figure 2 | E | Cumulative distribution of mean EU-STORM cluster size per nucleus (WT-DMSO vs. POLR2A-degraded) | `imaging_analysis_scripts/EUstorm.R` |
| Figure 2 | F | Cumulative distribution of blink density per nucleus (WT-DMSO vs. POLR2A-degraded) | `imaging_analysis_scripts/EUstorm.R` |
| Figure 2 | G | MA plot, Total RNA-seq, POLR2A-degraded cells | `total_RNA_seq/total_rnaseq_deseq2_analysis.R` |
| Figure 2 | H | MA plot, Total RNA-seq, ActD-treated cells | `total_RNA_seq/total_rnaseq_deseq2_analysis.R` |
| Figure 2 | I | MA plot, EU-seq nascent transcription, POLR2A-depleted cells | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| Figure 3 | A | EU-seq log2FC and steady-state RPKM density, by expression tertile (Low/Unexpressed, Average, High) | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| Figure 3 | B | Cumulative distribution of log2FC (POLR2A/WT), EU-seq (solid) vs. Total RNA-seq (dashed), by DEG status | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| Figure 3 | C | Metagene profiles, log2(POLR2A/WT) across gene body and intron length, EU-seq vs. Total RNA-seq | `EU_seq_analysis/EU_seq_log2ratio_coverage_analysis.R` |
| Figure 3 | D | Mean CUT&Tag signal (H3K4me3, POL2), Up-low-expression vs. Down-high-expression gene sets | `Chip_analysis/CT_compute_coverage.R` |
| Figure 3 | E | EU-seq log2FC by DEG status, faceted by E/I ratio category (Intron-rich/Average/Exon-rich) | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| Figure 3 | F | Distribution of EU-seq log2FC by E/I ratio category (Intron-rich/Average/Exon-rich) | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| Figure 4 | A | Mean insulation score vs. TAD-boundary distance, CTCF-anchored WT TADs, by expression | `HiC_analysis/tad_analysis/insulation_analysis.R` |
| Figure 4 | B | Same as A, for TADs called in POLR2A-degraded condition | `HiC_analysis/tad_analysis/insulation_analysis.R` |
| Figure 4 | C | Compartment saddle plot, PC1 eigenvector, Low/High expression bins, Wildtype | `HiC_analysis/compartment_analysis/saddle_analysis.R` |
| Figure 4 | D | Same as C, for POLR2A-degraded condition | `HiC_analysis/compartment_analysis/saddle_analysis.R` |
| Figure 4 | E | Per-corner log2 O/E compartment interaction strength (AA/BB/AB/BA), WT vs. POLR2A-degraded | `HiC_analysis/compartment_analysis/saddle_analysis.R` |
| Figure 4 | F | Contact-probability scaling |S| by E/I ratio percentile, WT vs. POLR2A-degraded; and delta-S | `HiC_analysis/scaling_analysis/scaling_analysis.R` |
| Figure 4 | G | Same as F, genes binned by expression percentile instead of E/I ratio | `HiC_analysis/scaling_analysis/scaling_analysis.R` |
| Figure 4 | H | S(d) as a function of genomic distance, WT vs. POLR2A-degraded, by TPM quartile | `HiC_analysis/scaling_analysis/scaling_analysis.R` |
| Figure 5 | A | Loop contact strength distribution (log10 contacts) by Gained/Stable/Lost, stratified by expression | `HiC_analysis/Hic_loop_analysis/hic_loop_lfc_analysis.R` |
| Figure 5 | B | APA of loops at highly expressed genes: WT-unique, POLR2A-degraded-unique, shared | `HiC_analysis/Hic_loop_analysis/final_APA_analysis.R` |
| Figure 5 | C | Number of Gained/Stable/Lost loops, by expression level and loop-length tertile | `HiC_analysis/Hic_loop_analysis/hic_loop_lfc_analysis.R` |
| Figure 5 | D | Loop contact strength for E-P and P-P loops (ENCODE ChIP-seq/ATAC-seq classified), by length tertile | `HiC_analysis/Hic_loop_analysis/hic_loop_lfc_analysis.R` |
| Figure 5 | E | APA of WT-unique loops, split by P-P and E-P classification | `HiC_analysis/Hic_loop_analysis/final_APA_analysis.R` |
| Figure 5 | F | Loop strength (PET count), public POLR2A ChIA-PET, by anchor element type | `chia_loop_analysis/chiapet_loop_strength_by_anchor_type.R` |
| Figure 5 | G | Loop strength (log10 SCALE contacts), Wildtype Hi-C, same anchor categories as F | `HiC_analysis/Hic_loop_analysis/hic_loop_strength_by_anchor_type.R` |
| Figure 6 | A | Hi-C contact matrix, COL26A1/MYL10/CUX1 locus, WT/POLR2A-degraded/difference, with loop arcs | `HiC_analysis/Hic_loop_analysis/hic_loop_track_plotter.R` |
| Figure 6 | B | Same as A, FBXO7/TIMP3/SYN3 locus | `HiC_analysis/Hic_loop_analysis/hic_loop_track_plotter.R` |
| Figure 6 | C | Schematic: intron-rich vs. exon-rich gene domain organization model | conceptual schematic, no script in repo |
| Figure 6 | D | Schematic: domain density vs. distance from center, core/optimal/surface zones | conceptual schematic, no script in repo |
| Figure 7 | A | Multiplexed SMLM of H3K9me3/H3K27ac/Pol-II in BJ fibroblasts, individual + merged channels | Provided by in-lab collaborator |
| Figure 7 | B | H3K9me3 cluster ROI example and zoomed representative heterochromatin core | Provided by in-lab collaborator |
| Figure 7 | C | Merged multi-color image of core from B, showing spatial relationship to H3K27ac/Pol-II | Provided by in-lab collaborator |
| Figure 7 | D | Fraction of H3K9me3 cores associated with H3K27ac / K27+Pol-II / unassociated | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 7 | E | Radial density of H3K9me3/Pol-II/H3K27ac from domain core, core/intermediate/boundary zones | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 7 | F | Cumulative fraction of public ChIP-seq signal (H3K27ac/H3K9me3/Pol-II) vs. intron length, by expression | `Chip_analysis/peak_scatter.R` |
| Figure 7 | G | KDE of intergenic gap width vs. ChIP mark, by expression level | `Chip_analysis/intergenic_coverage.R` |
| Figure 8 | A | Representative EdU-SMLM nucleus with DBSCAN domain-identification zoom-ins | raw microscopy image |
| Figure 8 | B | Representative EdU-SMLM images, Wildtype vs. POLR2A-degraded | raw microscopy image |
| Figure 8 | C | Radial density of EdU-labeled domains from domain core, WT vs. POLR2A-degraded | `imaging_analysis_scripts/EDU.K9_STORM.R` |
| Figure 8 | D | Domain density vs. domain radius, colored by size category, WT vs. POLR2A-degraded + counts | `imaging_analysis_scripts/EDU.K9_STORM.R` |
| Figure 8 | E | Representative H3K9me3-SMLM images, Wildtype vs. POLR2A-degraded | raw microscopy image |
| Figure 8 | F | Radial density of H3K9me3-defined domain cores from center, WT vs. POLR2A-degraded | `imaging_analysis_scripts/EDU.K9_STORM.R` |
| Figure 9 | A | PWS following treatment: average nuclear packing scaling (Dnucleus) | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 9 | B | PWS following treatment: diffusion (rate of chromatin clutch movement) | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 9 | C | PWS following treatment: fractional moving mass (FMM) | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 9 | D | ChromSTEM tomograms: intact vs. degraded packing domains | Provided by in-lab collaborator |
| Figure 9 | E | Domain classification by size/packing efficiency (Nascent/Small Mature/Large Mature/Decaying) + gain/loss bar chart | `imaging_analysis_scripts/Imaging_Analysis.R` |
| Figure 9 | F | Radial density by domain category (Large Mature/Small Mature/Nascent/Decaying), WT vs. POLR2A-degraded | `imaging_analysis_scripts/chrom_stem_radial_density.R` |
| Figure 10 | A | MMP7/MMP20 locus: ChIA-PET arcs + cCRE/gene tracks, CUT&Tag, Total RNA-seq + EU-seq | `chia_loop_analysis/updated_CHIA_PET_loci_plotter.R` |
| Figure 10 | B | Same as A, RAD50/IL13/IL4 locus | `chia_loop_analysis/updated_CHIA_PET_loci_plotter.R` |
| Figure 10 | C | Normalized log2 gene-body coverage (antibody/input), EZH2 and POLR2A RIP-seq, exonic/intronic | `RIPseq_analysis/rip_seq_logfc_analysis.R` |
| Figure 10 | D | POLR2A RIP-seq enrichment over input (log2 P2AB/WT_Input) by feature type and intron-length category, by expression | `RIPseq_analysis/rip_seq_logfc_analysis.R` |
| Figure 10 | E | Row-normalized density of RIP-seq signal termination position within introns and past TES | `RIPseq_analysis/rip_seq_intron_dog_heatmap.R` |
| Figure 10 | F | Cumulative distribution of EU-seq read termination position, short vs. long introns | `EU_seq_analysis/EU_seq_termination.R` |
| Figure 10 | G | Annotated vs. effective E/I ratio (RIP-seq), short/long introns, by expression | `RIPseq_analysis/rip_seq_effective_ei.R` |
| Figure 11 | A | Schematic: packing domain life cycle (unconstrained walk -> nascent domain -> mature domain) | conceptual schematic, no script in repo |
| Figure 11 | B | Schematic: model of healthy chromatin domain (Pol-II on exons, optimal zone) | conceptual schematic, no script in repo |
| Figure 11 | C | Schematic: model of dysregulated domain following Pol-II loss (intron exposure, readthrough) | conceptual schematic, no script in repo |
| SI Figure 1 | A | Widefield fluorescence, ActD vs. DMSO: DAPI/POLR2A/EU-labeled nascent RNA | raw microscopy image |
| SI Figure 1 | B | Coefficient of variation, DAPI-stained DNA, before/after ActD treatment | `imaging_analysis_scripts/Imaging_Analysis.R` |
| SI Figure 1 | C | Flow cytometry, EU-labeled nascent RNA, before/1hr/6hr treatment | `imaging_analysis_scripts/Imaging_Analysis.R` |
| SI Figure 1 | D | PWS following treatment: average diffusion and fractional moving mass | `imaging_analysis_scripts/Imaging_Analysis.R` |
| SI Figure 1 | E | Flow cytometry, target protein GFP signal, auxin time course | `imaging_analysis_scripts/Imaging_Analysis.R` |
| SI Figure 1 | F | Verification of POLR2A depletion: degron GFP vs. immunofluorescent signal | `imaging_analysis_scripts/Imaging_Analysis.R` |
| SI Figure 1 | G | Western blot, POLR2A-AID2 line, 6hr degradation | Provided by in-lab collaborator |
| SI Figure 1 | H | Nuclear scaling analysis (PWS), A549/HeLa/U2OS following treatment | `imaging_analysis_scripts/Imaging_Analysis.R` |
| SI Figure 1 | I | Sample Pearson correlation between RNA-seq replicates | `total_RNA_seq/total_rnaseq_deseq2_analysis.R` |
| SI Figure 2 | A | Core Pathway Gene Expression Is Stable Upon Pol-II Depletion (MA plot + ridgeline + heatmap: Core Apoptosis/Splicing/UPR) | `total_RNA_seq/damage_markers.R` |
| SI Figure 2 | B | Chromatin fractionation western blot confirming incomplete degradation of Pol-II-bound chromatin, compared to whole-cell degradation (SI Fig. 1G) | Provided by external collaborator |
| SI Figure 2 | C | Lamin B1 / Histone 3 western blot | Provided by external collaborator |
| SI Figure 2 | D | p53 / Histone 3 western blot | Provided by external collaborator |
| SI Figure 2 | E | Gamma H2AX / Histone 3 western blot | Provided by external collaborator |
| SI Figure 2 | F | Quantification of (B-E) as fold change relative to DMSO control | not yet identified — likely a real script (band quantification + stats), not just collaborator-provided; please confirm |
| SI Figure 3 | A | Pearson correlation of scaled EU-seq signal across replicates (NoEU/POLR2A/Wildtype), clustered | `total_RNA_seq/total_rnaseq_deseq2_analysis.R` |
| SI Figure 3 | B | Distribution of log2FC (POLR2A/WT) EU-seq per gene, colored by DE status | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| SI Figure 3 | C | Number of DEGs: EU-seq only, Total RNA-seq only, or both | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| SI Figure 3 | D | MA plot of EU-seq log2FC (POLR2A/WT) vs. log10(mean FPKM+1) | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| SI Figure 3 | E | EU-seq vs. Total RNA-seq log2FC comparison scatter, colored by concordance category | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| SI Figure 3 | F | PCA of EU-seq replicates (NoEU/POLR2A/Wildtype) | `EU_seq_analysis/EU_seq_zscore_analysis.R` |
| SI Figure 3 | G | Mean CUT&Tag signal (POL2, H3K4me3), Up/Down by expression tier, full facet set | `Chip_analysis/CT_compute_coverage.R` |
| SI Figure 4 | A | Mean EU-subtracted log2(POLR2A/WT), EU-seq + Total RNA-seq, gene body + intron, by E/I ratio | `EU_seq_analysis/EU_seq_log2ratio_coverage_analysis.R` |
| SI Figure 4 | B | Same as A, for upregulated genes by expression level | `EU_seq_analysis/EU_seq_log2ratio_coverage_analysis.R` |
| SI Figure 4 | C | Same as A, for downregulated genes by expression level | `EU_seq_analysis/EU_seq_log2ratio_coverage_analysis.R` |
| SI Figure 5 | A | Total valid Hi-C contacts per replicate (% of mean WT unsampled), 3 downsampling strategies | `HiC_analysis/HiC_stats_analysis/hic_stats.R` |
| SI Figure 5 | B | Pairwise SCC (HiCRep) between all 6 replicates, unsampled, 5kb, clustered | `HiC_analysis/HiC_stats_analysis/scc.r` |
| SI Figure 5 | C | Alignment outcome for ds_all-downsampled reads (unmapped/normal/chimeric/multi) | `HiC_analysis/HiC_stats_analysis/hic_stats.R` |
| SI Figure 5 | D | Classification of unique mapped read pairs, ds_all-downsampled (duplicates/intra-fragment/etc) | `HiC_analysis/HiC_stats_analysis/hic_stats.R` |
| SI Figure 5 | E | Library complexity (unique molecules) per replicate, 3 processing groups | `HiC_analysis/HiC_stats_analysis/hic_stats.R` |
| SI Figure 5 | F | Duplicate rate per replicate, 3 processing groups | `HiC_analysis/HiC_stats_analysis/hic_stats.R` |
| SI Figure 6 | A | ATA, CTCF-associated TAD boundaries called from POLR2A-degraded data ("pol2 positions") | `HiC_analysis/tad_analysis/ata_analysis.R` |
| SI Figure 6 | B | ATA, CTCF-associated TAD boundaries called from WT data ("WT positions") | `HiC_analysis/tad_analysis/ata_analysis.R` |
| SI Figure 6 | C | Compartment saddle plots, Medium expression tertile, WT and POLR2A-degraded | `HiC_analysis/compartment_analysis/compartment_analysis.R` |
| SI Figure 6 | D | Compartment saddle plots, POLR2A-degraded, all 3 coding-density tertiles | `HiC_analysis/compartment_analysis/compartment_analysis.R` |
| SI Figure 6 | E | Same as D, for WT cells | `HiC_analysis/compartment_analysis/compartment_analysis.R` |
| SI Figure 7 | A | Mean insulation score, non-CTCF WT TAD boundaries, WT+POLR2A-degraded, by expression | `HiC_analysis/tad_analysis/insulation_analysis.R` |
| SI Figure 7 | B | Per-TAD ATA score, non-CTCF WT TADs, WT+POLR2A-degraded, by expression | `HiC_analysis/tad_analysis/insulation_analysis.R` |
| SI Figure 7 | C | Same as A, non-CTCF POLR2A-degraded TAD boundaries | `HiC_analysis/tad_analysis/insulation_analysis.R` |
| SI Figure 7 | D | Same as B, non-CTCF POLR2A-degraded TADs | `HiC_analysis/tad_analysis/insulation_analysis.R` |
| SI Figure 7 | E | Per-TAD ATA score distributions, CTCF-associated WT TADs, by expression | `HiC_analysis/tad_analysis/insulation_analysis.R` |
| SI Figure 7 | F | Same as E, CTCF-associated POLR2A-degraded TADs | `HiC_analysis/tad_analysis/insulation_analysis.R` |
| SI Figure 7 | G | log2(coding/noncoding) distribution by coding-density tertile and condition, by compartment | `HiC_analysis/compartment_analysis/saddle_analysis.R` |
| SI Figure 7 | H | Per-corner log2(O/E) compartment interaction scores, Medium expression tertile | `HiC_analysis/compartment_analysis/saddle_analysis.R` |
| SI Figure 8 | A | Per-loop log2FC contact strength distribution, by loop classification (WT-only/Pol2-only/shared) | `HiC_analysis/Hic_loop_analysis/hic_loop_lfc_analysis.R` |
| SI Figure 8 | B | APA, P-P loops unique to POLR2A-degraded or shared | `HiC_analysis/Hic_loop_analysis/final_APA_analysis.R` |
| SI Figure 8 | C | Same as B, E-P loops | `HiC_analysis/Hic_loop_analysis/final_APA_analysis.R` |
| SI Figure 8 | D | APA, E-E loops, WT-unique/POLR2A-degraded-unique/shared | `HiC_analysis/Hic_loop_analysis/final_APA_analysis.R` |
| SI Figure 8 | E | Non-gene-body anchor classification (Promoter/Proximal-Enhancer/Distal-Enhancer), up/down genes | `HiC_analysis/Hic_loop_analysis/FINAL_EP_loop_analysis_hiccups.R` |
| SI Figure 8 | F | Loop contact strength, Gained/Lost loops at DEG anchors, by source (WT-only/Pol2-only/shared) | `HiC_analysis/Hic_loop_analysis/hic_loop_lfc_analysis.R` |
| SI Figure 9 | A | HiCCUPS loop fold enrichment vs. loop size, CTCF-anchored, by pair sub-type and exclusivity | `HiC_analysis/Hic_loop_analysis/final_loop_scatter_analysis.R` |
| SI Figure 9 | B | O/E signal distribution, CTCF loops unique to POLR2A-degraded, by pair sub-type | `HiC_analysis/Hic_loop_analysis/final_loop_scatter_analysis.R` |
| SI Figure 9 | C | Same as B, CTCF loops unique to WT | `HiC_analysis/Hic_loop_analysis/final_loop_scatter_analysis.R` |
| SI Figure 9 | D | Total loop counts per CTCF pair sub-type, WT vs. POLR2A-degraded | `HiC_analysis/Hic_loop_analysis/final_loop_scatter_analysis.R` |
| SI Figure 9 | E | Loop contact strength, pure CTCF loops (no Promoter/Enhancer overlap), Gained/Stable/Lost | `HiC_analysis/Hic_loop_analysis/FINAL_EP_loop_analysis_hiccups.R` |
| SI Figure 10 | A | KDE, EdU domain shape properties (r_t, r_min/r_max, r_g, circularity), WT vs. POLR2A-degraded | Provided by in-lab collaborator |
| SI Figure 10 | B | Regression, localizations vs. domain area (log-log), power-law fit, WT vs. POLR2A-degraded | Provided by in-lab collaborator |
| SI Figure 10 | C | Representative EdU domain: localizations, fitted boundary polygon, centroid | Provided by in-lab collaborator |
| SI Figure 10 | D | Density vs. radius, H3K9me3-defined domains, WT/Auxin, by core size category | `imaging_analysis_scripts/EDU.K9_STORM.R` |
| SI Figure 10 | E | KDE domain density + core counts, same 3 H3K9me3 size categories as D | `imaging_analysis_scripts/EDU.K9_STORM.R` |
| SI Figure 10 | F | Representative H3K9me3-SMLM image with core identification zoom-ins | Provided by in-lab collaborator |
| SI Figure 11 | A | STEM HAADF image, representative POLR2A-degraded nucleus | Provided by in-lab collaborator |
| SI Figure 11 | B | STEM HAADF image, representative WT (DMSO) nucleus | Provided by in-lab collaborator |
| SI Figure 11 | C | Example ChromSTEM tomogram, boxed region from A/B, with domain centers/sizes overlaid | Provided by in-lab collaborator |
| SI Figure 11 | D | ATAC-seq mean enrichment + log2(KO/DMSO) accessibility vs. local gene-count density, exon/intron | `atac_analysis/atac_gene_density.R` |
| SI Figure 11 | E | ATAC-seq accessibility change vs. local neighborhood transcription, by genomic E/I tertile | `atac_analysis/atac_analysis.R` |
| SI Figure 11 | F | Domain radius vs. CVC scatter, individual ChromSTEM domains, by size category | `imaging_analysis_scripts/Imaging_Analysis.R` |
| SI Figure 11 | G | CVC distributions by domain size category (Small/Average/Large), Auxin vs. DMSO | `imaging_analysis_scripts/Imaging_Analysis.R` |
| SI Figure 12 | A | HISAT2 concordant alignment breakdown per library | `RIPseq_analysis/rip_qc.R` |
| SI Figure 12 | B | Total and mapped read counts per library | `RIPseq_analysis/rip_qc.R` |
| SI Figure 12 | C | Fraction of bedGraph signal per genomic feature category, P2AB + WT Input | `RIPseq_analysis/rip_qc.R` |
| SI Figure 12 | D | Genome-wide replicate correlation between P2AB replicates | `RIPseq_analysis/rip_qc.R` |
| SI Figure 12 | E | P2AB IP enrichment over Input, exon vs. intron intervals | `RIPseq_analysis/rip_qc.R` |
| SI Figure 12 | F | Genome-wide P2AB IP enrichment over Input (log2, 1000bp bins) | `RIPseq_analysis/rip_qc.R` |
| SI Figure 13 | A | Public Pro-seq, RAD50 locus, POLR2A-degraded (6hr Auxin vs. DMSO) | IGV genome browser screenshot, no script |
| SI Figure 13 | B | Public Pro-seq, additional locus, transcriptional dysregulation | IGV genome browser screenshot, no script |
| SI Figure 13 | C | Public Pro-seq, additional locus, transcriptional dysregulation | IGV genome browser screenshot, no script |
| SI Figure 13 | D | Public Pro-seq, additional locus, transcriptional dysregulation | IGV genome browser screenshot, no script |
| SI Figure 13 | E | Cumulative distribution, RIP-seq read termination within exons, by exon-length tertile | `EU_seq_analysis/EU_seq_termination.R` |
| SI Figure 13 | F | Density of E/I ratio, Average intron-length tertile, genomic reference vs. RIP-seq (P2AB), by expression | `RIPseq_analysis/rip_seq_effective_ei.R` |
| SI Figure 13 | G | Annotated vs. effective E/I ratio, long-read K562 DMSO vs. PlaB, by intron-length tertile | `long_read_COP/RREA_analysis.R` |
| SI Figure 13 | H | Fraction of long reads unspliced/partially retained/retained, K562 DMSO vs. PlaB, by intron length | `long_read_COP/RREA_analysis.R` |

## Open items

- SI Figure 2F (quantification of western blots B-E) has no confirmed source yet --
  likely a real analysis script (band quantification + statistics), not simply
  collaborator-provided like the raw blots themselves. Needs confirmation.