# Revision scripts

Analysis scripts added during revision of the manuscript. Each script is
listed with the figure panels or tables it produces. Scripts were run in R
on Northwestern's Quest HPC cluster; input and output paths are set at the
top of each script and need to be updated for other environments. Public
datasets (Micro-C, Mint-ChIP-seq, ChIP-seq, cCREs, RSEM gene quantification)
are listed with accessions in Table 5 of the manuscript.

## hg38 EU-seq and total RNA-seq processing (`hg38_rnaseq_processing/`)

SLURM scripts that align the EU-seq and total RNA-seq to hg38, compute the
chrM (mitochondrial) scale factors and build the count tables and bigWigs read
by the R scripts below. Listed in run order.

| Script | Description |
|---|---|
| `build_hisat2_index_pol2.sh` | Downloads UCSC hg38 and the Ensembl 115 GTF (chr-prefixed) and builds the HISAT2 splice-aware index |
| `align_totalRNAseq_pol2.sh` | Total RNA-seq (single-end): TrimGalore, HISAT2, samtools filter/sort/index |
| `dedup_euseq_fastq_hg38.sh` | Removes repeated read names from the EU-seq paired-end FASTQs before alignment |
| `align_euseq_hg38.sh` | EU-seq (paired-end): FastQC, TrimGalore, HISAT2, samtools filter (-F 3852 -q 30)/sort/index |
| `euseq_subtract_hg38.sh` | Per-nucleotide NoEU background subtraction of EU-seq RPM bedGraphs, and library sizes |
| `compute_scale_factors_mito_hg38.sh` | EU-seq chrM scale factors (mean Wildtype chrM RPM / mean condition chrM RPM) |
| `compute_scale_factors_mito_total_hg38.sh` | Total RNA-seq chrM scale factors |
| `count_chrM_genes_hg38.sh` | Per-gene and per-bin chrM read counts for `transcriptomics/mito_qc.R` |
| `euseq_rpkm_hg38.sh` | EU-seq RPKM over gene bodies, exons, introns and downstream-of-gene regions (`euseq_rpkm_qc.txt`, `euseq_rpkm_analysis.txt`) |
| `merge_euseq_bams_hg38.sh` | Merges EU-seq BAMs per condition |
| `merge_totalrnaseq_bams_hg38.sh` | Merges total RNA-seq BAMs per condition |
| `generate_euseq_bigwigs_hg38.sh` | NoEU-subtracted, chrM-scaled EU-seq bigWigs (RPKM per nucleotide) and log2(POLR2A-degraded/Wildtype) bigWigs |
| `generate_euseq_log2ratio_bamcompare_hg38.sh` | EU-seq log2(POLR2A-degraded/Wildtype) bigWig from the merged BAMs (deepTools bamCompare, chrM-scaled) |
| `generate_total_bigwigs_mito_hg38.sh` | chrM-scaled total RNA-seq bigWigs and log2(POLR2A-degraded/Wildtype) bigWig from the merged BAMs |

## Public data downloads (`public_data/`)

| Script | Description |
|---|---|
| `download_microc_hic.sh` | ENCODE POLR2A-degron HCT116 Micro-C contact maps (ENCFF534TYD, ENCFF302CUI) |
| `03_download_mintchip.sh` | ENCODE Mint-ChIP-seq fold-change bigWigs and narrowPeaks for H3K27ac, H3K4me3 and H3K9me3 (POLR2A-degron HCT116, 5-Ph-IAA and DMSO, 6 h) |

The ENCODE POLR2A ChIP-seq bigWig (ENCSR000EUU, `POLR2A_ENCODE.bw`) and the
other public files listed in Table 5 were downloaded from the ENCODE portal.

## Transcriptomics (`transcriptomics/`)

| Script | Figure / table | Description |
|---|---|---|
| `EU_seq_zscore_analysis_updated.R` | Fig. 2i, Fig. 3a-b, Fig. 3e-f, SI Fig. 3c-f | EU-seq nascent transcription (hg38, chrM-scaled, NoEU-subtracted): sample correlation and PCA, Z-test differential analysis, comparison with total RNA-seq, and stratification by expression and exon/intron ratio |
| `updated_metagene_analysis.R` | Fig. 3c | Metagene profiles of log2(POLR2A-degraded/Wildtype) EU-seq and total RNA-seq signal over gene bodies and introns, by expression and DEG status |
| `mito_qc.R` | SI Fig. 3a-b | chrM coverage and chrM share of the library for EU-seq and total RNA-seq, with the chrM scale factors |

## Hi-C (`hic/`)

| Script | Figure / table | Description |
|---|---|---|
| `hic_contact_map_eigenvector.R` | SI Fig. 8b | Representative Hi-C contact map (Wildtype lower-left, POLR2A-degraded upper-right) with A/B compartment eigenvector tracks for each condition |

## Micro-C and contact scaling (`microc/`)

| Script | Figure / table | Description |
|---|---|---|
| `scaling_controls_main.R` | Shared | Functions and constants for contact-scaling analysis with length-matched non-genic control windows. Sourced by the three scripts below and by `effect_sizes/scaling_effect_size.R` |
| `scaling_controls_hic.R` | SI Fig. 5h | Hi-C S(d) and ΔS(d) for genes ≥100 kb by expression, with non-genic controls |
| `scaling_controls_microc.R` | SI Fig. 6a-d | Micro-C S(d), ΔS(d) and binned ΔS for genes ≥100 kb and 1–100 kb, with non-genic controls |
| `scaling_controls_microc_AB.R` | SI Fig. 6e-j | The Micro-C analysis run separately within A and B compartments |
| `scaling_chromatin_vs_transcription.R` | Shared | ΔS against the change in nascent transcription and in H3K27ac, H3K4me3 and H3K9me3 (Mint-ChIP-seq). Sourced by the script below |
| `scaling_chromatin_vs_neighborhood_transcription.R` | SI Fig. 12a-d | ΔS by tertile of the change in neighborhood nascent transcription (gene plus genes within 100 kb) and of each histone-mark change |
| `loop_anchor_analysis.R` | Fig. 4d-e, Supplementary Table 2 | Placement of Micro-C loop anchors on exons, cCREs and introns compared with random placement of the same loops, gene relationship of loop anchors, and loop length. Exports the per-loop anchor annotation used by the script below |
| `loop_anchor_contact_strength_analysis.R` | Fig. 4c | Micro-C loop strength by anchor pair for loops with small anchors (≤2 kb) |

`scaling_controls_hic.R`, `scaling_controls_microc.R`, `scaling_controls_microc_AB.R`
and `scaling_chromatin_vs_neighborhood_transcription.R` load their shared script
from the same directory.

## Hi-C loops and transcription (`loop_transcription/`)

| Script | Figure / table | Description |
|---|---|---|
| `loop_regression_main.R` | Shared | Builds the per-loop feature table (Wildtype and POLR2A-degraded HiCCUPS loops, contacts, log2 fold change, anchor Pol-II occupancy, histone marks, anchor-gene expression). Sourced by `fork_pooled_curves.R` |
| `fork_pooled_curves.R` | SI Fig. 10h-i, SI Fig. 12e | Pooled loop contact change across the change in EU-seq, total RNA-seq and histone marks at the anchor gene, and across transcription level in Wildtype and POLR2A-degraded cells |
| `eu_seq_stratified_loop_analysis.R` | SI Fig. 10e-g | Gained and Lost loops, EU-seq fold change, residual EU-seq signal and loop contact change, by tertile of absolute EU-seq loss at the anchor gene |

`fork_pooled_curves.R` loads `loop_regression_main.R` from the same directory.
`eu_seq_stratified_loop_analysis.R` uses the anchor classification functions in
`HiC_analysis/Hic_loop_analysis/EP_analysis_functions.R`.

## Effect sizes (`effect_sizes/`)

Each script writes a CSV of standardized effect sizes (Hedges' g or the
measure named below, with 95% confidence intervals) for POLR2A degradation and
ActD treatment. `compile_effect_size_table.R` combines them into Supplementary
Table 1.

| Script | Figure / table | Description |
|---|---|---|
| `pws_effect_size.R` | Supplementary Table 1 | PWS nuclear packing scaling, diffusion and fractional moving mass |
| `chromstem_domain_loss_effect_size.R` | Fig. 7e (bar chart), SI Fig. 14h, Supplementary Table 1 | ChromSTEM domain count ratios (treated to control) by domain category and size class, and domain packing scaling, radius and chromatin volume concentration |
| `eustorm_effect_size.R` | Supplementary Table 1 | EU-STORM cluster size and blink density |
| `smlm_effect_size.R` | Supplementary Table 1 | EdU and H3K9me3 SMLM domain counts, radius and density, by domain size category |
| `euseq_effect_size.R` | Supplementary Table 1 | EU-seq median log2 fold change, median per-gene Hedges' g and fraction of genes with log2 fold change < −1, overall and by expression tertile |
| `tad_insulation_effect_size.R` | SI Fig. 7c-f, Supplementary Table 1 | Boundary insulation score for TADs called in Wildtype and POLR2A-degraded cells, with and without CTCF-bounded boundaries |
| `compartment_effect_sizes.R` | SI Fig. 8e, Supplementary Table 1 | Genome-wide and per-chromosome compartment strength, and per-corner compartment interaction scores by expression |
| `scaling_effect_size.R` | Supplementary Table 1 | Contact-scaling slope S and the gradient of ΔS across exon/intron ratio and expression percentile |
| `loop_effect_size.R` | Supplementary Table 1 | Loop strength, fraction of loops lost and ratio of lost to gained loops, for all loops and by anchor type (P-P, E-P, E-E, promoter-anchored). Uses `HiC_analysis/Hic_loop_analysis/EP_analysis_functions.R` |
| `compile_effect_size_table.R` | Supplementary Table 1 | Combines the CSVs above into the Supplementary Table 1 file and its display version |
