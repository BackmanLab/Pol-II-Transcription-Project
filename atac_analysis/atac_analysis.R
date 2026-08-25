##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# ATAC-seq chromatin accessibility analysis for POLR2A AID
# degradation (6hr auxin) vs DMSO control in HCT116 cells.
#
# Three analyses:
#   1. Gene body accessibility (bigWig signal TSS-TES)
#   2. Promoter accessibility (bigWig signal TSS +/- 1kb)
#   3. Peak-level accessibility (narrowPeak overlap with genes)
#
# Each analysis stratified by:
#   - Gene expression quartile (TPM from total RNA-seq)
#   - Genomic E/I group (exon/intron length ratio from GTF)
#   - DEG status (DESeq2 results from matched total RNA-seq)
#
# BigWig values = MACS2 fold enrichment over local background
# NarrowPeak duplicates (merged replicates) resolved by keeping
# highest fold enrichment per unique coordinate
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(ggside)
  library(patchwork)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(BiocGenerics)
  library(scales)
})

##--------------------------------------------------------------## 1. Parameters

ATAC_DIR  <- "/projects/b1042/BackmanLab/Lucas/120525_pol2/data/atac"
DMSO_BW   <- file.path(ATAC_DIR, "pol2_dmso.bw")
KO_BW     <- file.path(ATAC_DIR, "pol2_ko.bw")
DMSO_NP   <- file.path(ATAC_DIR, "pol2_dmso.narrowPeak")
KO_NP     <- file.path(ATAC_DIR, "pol2_ko.narrowPeak")

GTF_PATH  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
TPM_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/HCT116_RNAseq_quantification.tsv"
DEG_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"

PROMOTER_WINDOW <- 1000   ## bp each side of TSS
LFC_THRESHOLD   <- 1      ## |log2FC| > 1 for DEG classification
PADJ_THRESHOLD  <- 0.05
TPM_THRESHOLD   <- 1      ## minimum TPM to include gene

## Colors
CONDITION_COLORS <- c(KO = "#7B2D8B", DMSO = "#1B7837")
EI_COLORS <- c("Low E/I" = "#0072B2", "Medium E/I" = "#E69F00", "High E/I" = "#CC79A7")
EXPR_COLORS <- c("Q1 (Low)" = "#0072B2", "Q2" = "#56B4E9",
                 "Q3" = "#E69F00", "Q4 (High)" = "#D55E00")
DEG_COLORS <- c("Up" = "#D55E00", "Down" = "#0072B2",
                "Not DE" = "grey70")

base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    strip.text       = element_text(face = "bold", size = 11)
  )

##--------------------------------------------------------------## 2. Load gene annotations from GTF

message("Building gene annotations from GTF...")

txdb         <- GenomicFeatures::makeTxDbFromGFF(GTF_PATH, format = "gtf",
                                                 taxonomyId = 9606)
genes_gr     <- GenomicFeatures::genes(txdb)

## Gene body GRanges
gene_body_gr <- genes_gr
gene_body_gr$gene_id <- sub("\\..*", "", gene_body_gr$gene_id)

## Promoter GRanges: TSS +/- PROMOTER_WINDOW
promoter_gr  <- GenomicRanges::promoters(genes_gr,
                                         upstream   = PROMOTER_WINDOW,
                                         downstream = PROMOTER_WINDOW)
promoter_gr$gene_id <- sub("\\..*", "", promoter_gr$gene_id)

cat("Gene bodies:", length(gene_body_gr), "\n")
cat("Promoters:  ", length(promoter_gr),  "\n")

##--------------------------------------------------------------## 3. Compute genomic E/I groups

message("Computing genomic E/I groups...")

exons_by_gene <- GenomicRanges::reduce(
  GenomicFeatures::exonsBy(txdb, by = "gene"))
exon_lengths  <- sum(BiocGenerics::width(exons_by_gene))
gene_lengths  <- setNames(BiocGenerics::width(genes_gr), genes_gr$gene_id)
gtf_common    <- intersect(names(exon_lengths), names(gene_lengths))
intron_len    <- pmax(gene_lengths[gtf_common] - exon_lengths[gtf_common], 0)
ei_genomic    <- exon_lengths[gtf_common] / (intron_len + 1)

genomic_ei_df <- data.frame(
  gene_id          = sub("\\..*", "", gtf_common),
  gene_length      = as.numeric(gene_lengths[gtf_common]),
  ei_genomic_ratio = as.numeric(ei_genomic),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(ei_genomic_ratio), ei_genomic_ratio > 0) %>%
  dplyr::mutate(
    genomic_ei_group = dplyr::case_when(
      ei_genomic_ratio < quantile(ei_genomic_ratio, 1/3, na.rm = TRUE) ~ "Low E/I",
      ei_genomic_ratio < quantile(ei_genomic_ratio, 2/3, na.rm = TRUE) ~ "Medium E/I",
      TRUE                                                               ~ "High E/I"
    ),
    genomic_ei_group = factor(genomic_ei_group,
                              levels = c("Low E/I","Medium E/I","High E/I"))
  )

cat("E/I groups computed:", nrow(genomic_ei_df), "genes\n")
print(table(genomic_ei_df$genomic_ei_group))

##--------------------------------------------------------------## 4. Load expression and DEG data

message("Loading expression and DEG data...")

## TPM values — strip Ensembl version suffix
tpm_raw <- read_tsv(TPM_FILE, show_col_types = FALSE)

## Identify TPM column — typically "TPM" or contains "tpm"
tpm_col <- names(tpm_raw)[grepl("tpm|TPM", names(tpm_raw))][1]
id_col  <- names(tpm_raw)[1]

expr_df <- tpm_raw %>%
  dplyr::rename(gene_id = !!id_col, tpm = !!tpm_col) %>%
  dplyr::mutate(
    gene_id = sub(".*\\|(ENSG[0-9]+)\\.[0-9]+\\|.*", "\\1", gene_id)
  ) %>%
  dplyr::select(gene_id, tpm) %>%
  dplyr::filter(!is.na(tpm), tpm >= 0)

cat("Genes with TPM data:", nrow(expr_df), "\n")
cat("Genes with TPM >= ", TPM_THRESHOLD, ":",
    sum(expr_df$tpm >= TPM_THRESHOLD), "\n")

## DESeq2 results
deg_raw <- read_csv(DEG_FILE, show_col_types = FALSE)

## Identify gene ID column
id_col_deg <- names(deg_raw)[grepl("gene|Gene|id|ID", names(deg_raw))][1]

deg_df <- deg_raw %>%
  dplyr::rename(gene_id = !!id_col_deg) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
  dplyr::mutate(
    deg_status = dplyr::case_when(
      !is.na(padj) & padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD ~ "Up",
      !is.na(padj) & padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD ~ "Down",
      TRUE ~ "Not DE"
    ),
    deg_status = factor(deg_status, levels = c("Up","Not DE","Down"))
  ) %>%
  dplyr::select(gene_id, log2FoldChange, padj, deg_status)

cat("DEG counts:\n")
print(table(deg_df$deg_status))

gene_meta <- genomic_ei_df %>%
  dplyr::left_join(expr_df, by = "gene_id") 

## Build master gene annotation table
gene_meta <- genomic_ei_df %>%
  dplyr::left_join(expr_df, by = "gene_id") %>%
  dplyr::left_join(deg_df,  by = "gene_id") %>%
  dplyr::filter(!is.na(tpm)) %>%
  dplyr::mutate(
    expr_group = dplyr::case_when(
      tpm < quantile(tpm[tpm >= TPM_THRESHOLD], 0.25, na.rm = TRUE) ~ "Q1 (Low)",
      tpm < quantile(tpm[tpm >= TPM_THRESHOLD], 0.50, na.rm = TRUE) ~ "Q2",
      tpm < quantile(tpm[tpm >= TPM_THRESHOLD], 0.75, na.rm = TRUE) ~ "Q3",
      tpm >= quantile(tpm[tpm >= TPM_THRESHOLD], 0.75, na.rm = TRUE) ~ "Q4 (High)",
      TRUE ~ NA_character_
    ),
    expr_group = factor(expr_group,
                        levels = c("Q1 (Low)","Q2","Q3","Q4 (High)"))
  ) %>%
  dplyr::filter(!is.na(expr_group))

cat("Genes in master annotation:", nrow(gene_meta), "\n")

##--------------------------------------------------------------## 5. Helper: extract bigWig signal over GRanges

library(GenomeInfoDb) 

STD_CHRS <- paste0("chr", c(1:22, "X", "Y"))

extract_bw_signal <- function(bw_path, gr, label = "") {
  message("  Importing bigWig: ", basename(bw_path), " over ", length(gr), " regions...")
  
  gr_filt <- gr[as.character(seqnames(gr)) %in% STD_CHRS]
  gr_filt <- GenomeInfoDb::keepSeqlevels(gr_filt,
                                         value        = STD_CHRS,
                                         pruning.mode = "coarse")
  
  bw_rle <- rtracklayer::import.bw(bw_path,
                                   which = gr_filt,
                                   as    = "RleList")
  
  common_chrs <- intersect(names(bw_rle), seqlevels(gr_filt))
  bw_rle  <- bw_rle[common_chrs]
  gr_filt <- GenomeInfoDb::keepSeqlevels(gr_filt,
                                         value        = common_chrs,
                                         pruning.mode = "coarse")
  
  signal <- GenomicRanges::binnedAverage(
    bins    = gr_filt,
    numvar  = bw_rle,
    varname = "score"
  )$score
  
  data.frame(
    gene_id = gr_filt$gene_id,
    signal  = signal,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::rename(!!label := signal)
}

##==============================================================##
##              ANALYSIS 1: GENE BODY ACCESSIBILITY             ##
##==============================================================##

message("\n========== ANALYSIS 1: GENE BODY ACCESSIBILITY ==========\n")

## Subset gene body GRanges to expressed genes
expressed_ids <- gene_meta$gene_id
gb_gr <- gene_body_gr[gene_body_gr$gene_id %in% expressed_ids]
cat("Gene bodies to profile:", length(gb_gr), "\n")

seqlevelsStyle(gb_gr) <- "UCSC"

## Extract signal from both bigWigs
gb_dmso <- extract_bw_signal(DMSO_BW, gb_gr, "dmso_signal")
gb_ko   <- extract_bw_signal(KO_BW,   gb_gr, "ko_signal")

## Compute LFC
gb_df <- gb_dmso %>%
  dplyr::inner_join(gb_ko, by = "gene_id") %>%
  dplyr::mutate(
    atac_lfc = log2(ko_signal + 0.01) - log2(dmso_signal + 0.01)
  ) %>%
  dplyr::left_join(gene_meta, by = "gene_id") %>%
  dplyr::filter(!is.na(genomic_ei_group), !is.na(expr_group))

cat("Gene body ATAC LFC summary:\n")
print(summary(gb_df$atac_lfc))

## --- Plot A1-1: Gene body ATAC LFC by expression quartile ---
a1_wilcox_expr <- gb_df %>%
  dplyr::group_by(expr_group) %>%
  dplyr::summarise(
    p_value    = wilcox.test(atac_lfc, mu = 0)$p.value,
    median_lfc = median(atac_lfc, na.rm = TRUE),
    n          = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = max(gb_df$atac_lfc, na.rm = TRUE) * 0.9
  )

library(ggdist)

## --- Plot A1-1: Gene body ATAC LFC by expression quartile --- rainplot

## Pre-compute n per group for labels
expr_n_labels <- gb_df %>%
  dplyr::count(expr_group) %>%
  dplyr::mutate(label = paste0("n=", scales::comma(n)))

## Recompute wilcox labels with y_pos inside plot limits
a1_wilcox_expr <- gb_df %>%
  dplyr::group_by(expr_group) %>%
  dplyr::summarise(
    p_value    = wilcox.test(atac_lfc, mu = 0)$p.value,
    median_lfc = median(atac_lfc, na.rm = TRUE),
    n          = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = quantile(gb_df$atac_lfc, 0.99, na.rm = TRUE) * 1.05
  )

p_gb_expr <- ggplot(gb_df,
                    aes(x = expr_group, y = atac_lfc,
                        fill = expr_group, color = expr_group)) +
  ggdist::stat_halfeye(
    adjust = 0.5, width = 0.5,
    justification = -0.05,
    .width = 0, point_colour = NA, alpha = 0.8
  ) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA,
    color = "grey30", linewidth = 0.5, alpha = 0.8,
    position = position_nudge(x = -0.12)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_text(data = a1_wilcox_expr,
            aes(x = expr_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4.5, color = "grey20") +
  geom_text(data = expr_n_labels,
            aes(x = expr_group,
                y = quantile(gb_df$atac_lfc, 0.99, na.rm=TRUE) * 1.13,
                label = label),
            inherit.aes = FALSE, size = 3.5, color = "grey30") +
  scale_fill_manual(values  = EXPR_COLORS) +
  scale_color_manual(values = EXPR_COLORS) +
  coord_cartesian(ylim = c(
    quantile(gb_df$atac_lfc, 0.01, na.rm = TRUE),
    quantile(gb_df$atac_lfc, 0.99, na.rm = TRUE) * 1.25
  )) +
  labs(x = NULL,
       y = "log2 FC ATAC signal\n(KO / DMSO)",
       title    = "Gene body accessibility by expression level",
       subtitle = "One-sample Wilcoxon vs zero  |  ATAC fold enrichment signal") +
  base_theme + theme(legend.position = "none")

print(p_gb_expr)

## --- Plot A1-2: Gene body ATAC LFC by genomic E/I group ---
a1_wilcox_ei <- gb_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(
    p_value    = wilcox.test(atac_lfc, mu = 0)$p.value,
    median_lfc = median(atac_lfc, na.rm = TRUE),
    n          = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = max(gb_df$atac_lfc, na.rm = TRUE) * 0.9
  )


## --- Plot A1-2: Gene body ATAC LFC by genomic E/I group --- rainplot
ei_n_labels <- gb_df %>%
  dplyr::count(genomic_ei_group) %>%
  dplyr::mutate(label = paste0("n=", scales::comma(n)))

a1_wilcox_ei <- gb_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(
    p_value    = wilcox.test(atac_lfc, mu = 0)$p.value,
    median_lfc = median(atac_lfc, na.rm = TRUE),
    n          = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = quantile(gb_df$atac_lfc, 0.99, na.rm = TRUE) * 1.05
  )

p_gb_ei <- ggplot(gb_df,
                  aes(x = genomic_ei_group, y = atac_lfc,
                      fill = genomic_ei_group, color = genomic_ei_group)) +
  ggdist::stat_halfeye(
    adjust = 0.5, width = 0.5,
    justification = -0.05,
    .width = 0, point_colour = NA, alpha = 0.8
  ) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA,
    color = "grey30", linewidth = 0.5, alpha = 0.8,
    position = position_nudge(x = -0.12)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_text(data = a1_wilcox_ei,
            aes(x = genomic_ei_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4.5, color = "grey20") +
  geom_text(data = ei_n_labels,
            aes(x = genomic_ei_group,
                y = quantile(gb_df$atac_lfc, 0.99, na.rm=TRUE) * 1.13,
                label = label),
            inherit.aes = FALSE, size = 3.5, color = "grey30") +
  scale_fill_manual(values  = EI_COLORS) +
  scale_color_manual(values = EI_COLORS) +
  scale_x_discrete(labels = c("Low E/I"    = "Low E/I\n(intron-rich)",
                              "Medium E/I" = "Medium E/I",
                              "High E/I"   = "High E/I\n(exon-rich)")) +
  coord_cartesian(ylim = c(
    quantile(gb_df$atac_lfc, 0.01, na.rm = TRUE),
    quantile(gb_df$atac_lfc, 0.99, na.rm = TRUE) * 1.25
  )) +
  labs(x = NULL,
       y = "log2 FC ATAC signal\n(KO / DMSO)",
       title    = "Gene body accessibility by genomic E/I group",
       subtitle = "One-sample Wilcoxon vs zero  |  ATAC fold enrichment signal") +
  base_theme + theme(legend.position = "none")

print(p_gb_ei)

## --- Plot A1-3: Gene body ATAC LFC by DEG status ---
p_gb_deg <- ggplot(gb_df %>% dplyr::filter(!is.na(deg_status)),
                   aes(x = deg_status, y = atac_lfc,
                       fill = deg_status)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, width = 0.6,
               alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = DEG_COLORS) +
  coord_cartesian(ylim = c(
    quantile(gb_df$atac_lfc, 0.01, na.rm = TRUE),
    quantile(gb_df$atac_lfc, 0.99, na.rm = TRUE)
  )) +
  labs(x = NULL,
       y = "log2 FC ATAC signal\n(KO / DMSO)",
       title    = "Gene body accessibility by DEG status",
       subtitle = "Total RNA-seq DESeq2 |LFC| > 1, padj < 0.05") +
  base_theme + theme(legend.position = "none")

print(p_gb_deg)

## --- Plot A1-4: ATAC LFC vs RNA-seq LFC scatter with marginal boxplots ---
ATAC_LFC_THRESHOLD=1

gb_rna_df <- gb_df %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(deg_status)) %>%
  dplyr::mutate(
    quadrant = dplyr::case_when(
      atac_lfc >  ATAC_LFC_THRESHOLD & deg_status == "Up"   ~ "Open & more transcribed",
      atac_lfc < -ATAC_LFC_THRESHOLD & deg_status == "Down" ~ "Closed & less transcribed",
      atac_lfc >  ATAC_LFC_THRESHOLD & deg_status == "Down" ~ "Open & less transcribed",
      atac_lfc < -ATAC_LFC_THRESHOLD & deg_status == "Up"   ~ "Closed & more transcribed",
      TRUE                                                   ~ "No change"
    ),
    quadrant = factor(quadrant,
                      levels = c("Open & more transcribed",
                                 "Closed & less transcribed",
                                 "Open & less transcribed",
                                 "Closed & more transcribed",
                                 "No change"))
  )

x_lims <- quantile(gb_rna_df$log2FoldChange, c(0.005, 0.995), na.rm = TRUE)
y_lims <- quantile(gb_rna_df$atac_lfc,       c(0.005, 0.995), na.rm = TRUE)

quad_labels <- gb_rna_df %>%
  dplyr::count(quadrant) %>%
  dplyr::mutate(
    pct   = round(n / sum(n) * 100, 1),
    label = paste0("n=", scales::comma(n), "\n(", pct, "%)"),
    x_pos = dplyr::case_when(
      quadrant == "Open & more transcribed"   ~ x_lims[2] * 0.85,
      quadrant == "Closed & less transcribed" ~ x_lims[1] * 0.85,
      quadrant == "Open & less transcribed"   ~ x_lims[1] * 0.85,
      quadrant == "Closed & more transcribed" ~ x_lims[2] * 0.85,
      TRUE ~ NA_real_
    ),
    y_pos = dplyr::case_when(
      quadrant == "Open & more transcribed"   ~ y_lims[2] * 0.85,
      quadrant == "Closed & less transcribed" ~ y_lims[1] * 0.85,
      quadrant == "Open & less transcribed"   ~ y_lims[2] * 0.85,
      quadrant == "Closed & more transcribed" ~ y_lims[1] * 0.85,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(x_pos))

quad_colors <- c(
  "Open & more transcribed"   = "#D55E00",
  "Closed & less transcribed" = "#0072B2",
  "Open & less transcribed"   = "#E69F00",
  "Closed & more transcribed" = "#56B4E9",
  "No change"                 = "grey80"
)

marginal_df <- gb_rna_df %>% dplyr::filter(quadrant != "No change")
r_atac_rna  <- cor(gb_rna_df$atac_lfc, gb_rna_df$log2FoldChange, use = "complete.obs")

p_gb_scatter <- ggplot(
  gb_rna_df %>% dplyr::arrange(quadrant == "No change"),
  aes(x = log2FoldChange, y = atac_lfc, color = quadrant)
) +
  geom_point(size = 0.6, alpha = 0.5) +
  ggside::geom_xsideboxplot(data = marginal_df,
                            aes(y = log2FoldChange, fill = quadrant),
                            orientation = "y", outlier.shape = NA, linewidth = 0.4, width = 0.6) +
  ggside::geom_ysideboxplot(data = marginal_df,
                            aes(x = atac_lfc, fill = quadrant),
                            orientation = "x", outlier.shape = NA, linewidth = 0.4, width = 0.6) +
  geom_hline(yintercept = c(-ATAC_LFC_THRESHOLD, ATAC_LFC_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_text(data = quad_labels,
            aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3.2, color = "grey20", hjust = 0.5) +
  scale_color_manual(values = quad_colors) +
  scale_fill_manual(values  = quad_colors) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(r_atac_rna, 3)),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "grey40",
           fontface = "italic") +
  ggside::scale_xsidey_continuous(breaks = NULL) +
  ggside::scale_ysidex_continuous(breaks = NULL) +
  ggside::theme_ggside_void() +
  coord_cartesian(xlim = x_lims, ylim = y_lims) +
  labs(x        = "RNA-seq log2 FC (POLR2A / WT)",
       y        = "ATAC log2 FC gene body (KO / DMSO)",
       title    = "Gene body accessibility vs transcriptional change",
       subtitle = paste0("Colored = DEG (padj < ", PADJ_THRESHOLD,
                         ", |LFC| > ", LFC_THRESHOLD,
                         ") AND ATAC |LFC| > ", ATAC_LFC_THRESHOLD)) +
  base_theme +
  theme(legend.position    = "right",
        legend.text        = element_text(size = 9),
        ggside.panel.scale = 0.15)

print(p_gb_scatter)

## --- Plot A1-5: ATAC LFC vs RNA-seq LFC scatter colored by expression level ---

gb_expr_df <- gb_df %>%
  dplyr::filter(!is.na(tpm), tpm > 0) %>%
  dplyr::mutate(
    log10_tpm  = log10(tpm + 0.01),
    atac_group = dplyr::case_when(
      atac_lfc >  ATAC_LFC_THRESHOLD ~ "More accessible",
      atac_lfc < -ATAC_LFC_THRESHOLD ~ "Less accessible",
      TRUE                            ~ "No change"
    ),
    atac_group = factor(atac_group,
                        levels = c("More accessible", "No change",
                                   "Less accessible"))
  )

x_lims2 <- quantile(gb_expr_df$log10_tpm, c(0.005, 0.995), na.rm = TRUE)
y_lims2 <- quantile(gb_expr_df$atac_lfc,  c(0.005, 0.995), na.rm = TRUE)

## Count and position labels for each ATAC group
atac_labels <- gb_expr_df %>%
  dplyr::count(atac_group) %>%
  dplyr::mutate(
    pct   = round(n / sum(n) * 100, 1),
    label = paste0("n=", scales::comma(n), "\n(", pct, "%)"),
    x_pos = mean(x_lims2),
    y_pos = dplyr::case_when(
      atac_group == "More accessible" ~ y_lims2[2] * 0.85,
      atac_group == "Less accessible" ~ y_lims2[1] * 0.85,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(y_pos))

atac_colors <- c(
  "More accessible" = "#D55E00",
  "No change"       = "grey80",
  "Less accessible" = "#0072B2"
)

marginal_expr_df <- gb_expr_df %>% dplyr::filter(atac_group != "No change")
r_atac_expr <- cor(gb_expr_df$atac_lfc, gb_expr_df$log10_tpm, use = "complete.obs")

p_gb_scatter_expr <- ggplot(
  gb_expr_df %>% dplyr::arrange(atac_group == "No change"),
  aes(x = log10_tpm, y = atac_lfc, color = atac_group)
) +
  geom_point(size = 0.6, alpha = 0.5) +
  ggside::geom_xsideboxplot(data = marginal_expr_df,
                            aes(y = log10_tpm, fill = atac_group),
                            orientation = "y", outlier.shape = NA, linewidth = 0.4, width = 0.6) +
  ggside::geom_ysideboxplot(data = marginal_expr_df,
                            aes(x = atac_lfc, fill = atac_group),
                            orientation = "x", outlier.shape = NA, linewidth = 0.4, width = 0.6) +
  geom_hline(yintercept = c(-ATAC_LFC_THRESHOLD, ATAC_LFC_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_text(data = atac_labels,
            aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3.2, color = "grey20", hjust = 0.5) +
  scale_color_manual(values = atac_colors) +
  scale_fill_manual(values  = atac_colors) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(r_atac_expr, 3)),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "grey40",
           fontface = "italic") +
  ggside::scale_xsidey_continuous(breaks = NULL) +
  ggside::scale_ysidex_continuous(breaks = NULL) +
  ggside::theme_ggside_void() +
  coord_cartesian(xlim = x_lims2, ylim = y_lims2) +
  labs(x        = "Gene expression log10(TPM)",
       y        = "ATAC log2 FC gene body (KO / DMSO)",
       title    = "Gene body accessibility vs expression level",
       subtitle = paste0("Colored = ATAC |LFC| > ", ATAC_LFC_THRESHOLD)) +
  base_theme +
  theme(legend.position    = "right",
        legend.text        = element_text(size = 9),
        ggside.panel.scale = 0.15)

print(p_gb_scatter_expr)

##==============================================================##
##              ANALYSIS 2: PROMOTER ACCESSIBILITY              ##
##==============================================================##

message("\n========== ANALYSIS 2: PROMOTER ACCESSIBILITY ==========\n")

## Subset promoter GRanges to expressed genes
pr_gr <- promoter_gr[promoter_gr$gene_id %in% expressed_ids]
cat("Promoters to profile:", length(pr_gr), "\n")

seqlevelsStyle(pr_gr) <- "UCSC"


## Extract signal
pr_dmso <- extract_bw_signal(DMSO_BW, pr_gr, "dmso_signal")
pr_ko   <- extract_bw_signal(KO_BW,   pr_gr, "ko_signal")

## Compute LFC
pr_df <- pr_dmso %>%
  dplyr::inner_join(pr_ko, by = "gene_id") %>%
  dplyr::mutate(
    atac_lfc = log2(ko_signal + 0.01) - log2(dmso_signal + 0.01)
  ) %>%
  dplyr::left_join(gene_meta, by = "gene_id") %>%
  dplyr::filter(!is.na(genomic_ei_group), !is.na(expr_group))

cat("Promoter ATAC LFC summary:\n")
print(summary(pr_df$atac_lfc))

## --- Plot A2-1: Promoter ATAC LFC by expression quartile ---
a2_wilcox_expr <- pr_df %>%
  dplyr::group_by(expr_group) %>%
  dplyr::summarise(
    p_value    = wilcox.test(atac_lfc, mu = 0)$p.value,
    median_lfc = median(atac_lfc, na.rm = TRUE),
    n          = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = max(pr_df$atac_lfc, na.rm = TRUE) * 0.9
  )

p_pr_expr <- ggplot(pr_df,
                    aes(x = expr_group, y = atac_lfc,
                        fill = expr_group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, width = 0.6,
               alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_text(data = a2_wilcox_expr,
            aes(x = expr_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, color = "grey20") +
  scale_fill_manual(values = EXPR_COLORS) +
  coord_cartesian(ylim = c(
    quantile(pr_df$atac_lfc, 0.01, na.rm = TRUE),
    quantile(pr_df$atac_lfc, 0.99, na.rm = TRUE) * 1.15
  )) +
  labs(x = NULL,
       y = "log2 FC ATAC signal\n(KO / DMSO)",
       title    = paste0("Promoter accessibility by expression level (TSS +/- ",
                         PROMOTER_WINDOW/1000, "kb)"),
       subtitle = "One-sample Wilcoxon vs zero  |  ATAC fold enrichment signal") +
  base_theme + theme(legend.position = "none")

print(p_pr_expr)

## --- Plot A2-2: Promoter ATAC LFC by genomic E/I group ---
a2_wilcox_ei <- pr_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(
    p_value    = wilcox.test(atac_lfc, mu = 0)$p.value,
    median_lfc = median(atac_lfc, na.rm = TRUE),
    n          = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = max(pr_df$atac_lfc, na.rm = TRUE) * 0.9
  )

p_pr_ei <- ggplot(pr_df,
                  aes(x = genomic_ei_group, y = atac_lfc,
                      fill = genomic_ei_group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, width = 0.6,
               alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_text(data = a2_wilcox_ei,
            aes(x = genomic_ei_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, color = "grey20") +
  scale_fill_manual(values = EI_COLORS) +
  scale_x_discrete(labels = c("Low E/I"    = "Low E/I\n(intron-rich)",
                              "Medium E/I" = "Medium E/I",
                              "High E/I"   = "High E/I\n(exon-rich)")) +
  coord_cartesian(ylim = c(
    quantile(pr_df$atac_lfc, 0.01, na.rm = TRUE),
    quantile(pr_df$atac_lfc, 0.99, na.rm = TRUE) * 1.15
  )) +
  labs(x = NULL,
       y = "log2 FC ATAC signal\n(KO / DMSO)",
       title    = paste0("Promoter accessibility by genomic E/I group (TSS +/- ",
                         PROMOTER_WINDOW/1000, "kb)"),
       subtitle = "One-sample Wilcoxon vs zero  |  ATAC fold enrichment signal") +
  base_theme + theme(legend.position = "none")

print(p_pr_ei)

## --- Plot A2-3: Promoter vs gene body ATAC LFC scatter ---
pr_gb_df <- pr_df %>%
  dplyr::select(gene_id, pr_lfc = atac_lfc) %>%
  dplyr::inner_join(
    gb_df %>% dplyr::select(gene_id, gb_lfc = atac_lfc,
                            deg_status, expr_group),
    by = "gene_id"
  )

r_pr_gb <- cor(pr_gb_df$pr_lfc, pr_gb_df$gb_lfc, use = "complete.obs")

p_pr_gb <- ggplot(pr_gb_df,
                  aes(x = gb_lfc, y = pr_lfc, color = deg_status)) +
  geom_point(size = 0.6, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = DEG_COLORS) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(r_pr_gb, 3)),
           hjust = 1.1, vjust = -0.5,
           size = 4, color = "grey30", fontface = "italic") +
  labs(x = "Gene body ATAC log2 FC",
       y = "Promoter ATAC log2 FC",
       title    = "Promoter vs gene body accessibility change",
       subtitle = paste0("Diagonal = equal change  |  r = ",
                         round(r_pr_gb, 3))) +
  base_theme

print(p_pr_gb)

##==============================================================##
##           ANALYSIS 3: PEAK-LEVEL ACCESSIBILITY               ##
##==============================================================##

message("\n========== ANALYSIS 3: PEAK-LEVEL ACCESSIBILITY ==========\n")

## Load and deduplicate narrowPeak files
## Columns: chr, start, end, name, score, strand,
##          fold_enrichment, -log10p, -log10q, summit_offset
np_cols <- c("chr","start","end","name","score","strand",
             "fold_enrichment","neg_log10p","neg_log10q","summit_offset")

load_narrowpeak <- function(path, condition) {
  read_tsv(path, col_names = np_cols, show_col_types = FALSE) %>%
    ## Deduplicate: keep highest fold enrichment per unique coordinate
    dplyr::group_by(chr, start, end) %>%
    dplyr::slice_max(fold_enrichment, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(condition = condition,
                  peak_id   = paste0(chr, ":", start, "-", end))
}

message("Loading narrowPeak files...")
dmso_peaks <- load_narrowpeak(DMSO_NP, "DMSO")
ko_peaks   <- load_narrowpeak(KO_NP,   "KO")

cat("DMSO peaks after deduplication:", nrow(dmso_peaks), "\n")
cat("KO peaks after deduplication:  ", nrow(ko_peaks),   "\n")

## Convert to GRanges for intersection
peaks_to_gr <- function(np_df) {
  GenomicRanges::GRanges(
    seqnames = np_df$chr,
    ranges   = IRanges::IRanges(start = np_df$start + 1,
                                end   = np_df$end),
    peak_id         = np_df$peak_id,
    fold_enrichment = np_df$fold_enrichment,
    neg_log10p      = np_df$neg_log10p
  )
}

dmso_gr <- peaks_to_gr(dmso_peaks)
ko_gr   <- peaks_to_gr(ko_peaks)

## Annotate peaks to gene bodies and promoters
annotate_peaks <- function(peaks_gr, genes_gr, promoters_gr) {
  ## Gene body overlaps
  gb_hits  <- GenomicRanges::findOverlaps(peaks_gr, genes_gr,
                                          ignore.strand = TRUE)
  ## Promoter overlaps
  pr_hits  <- GenomicRanges::findOverlaps(peaks_gr, promoters_gr,
                                          ignore.strand = TRUE)
  
  gb_annot <- data.frame(
    peak_idx = queryHits(gb_hits),
    gene_id  = genes_gr$gene_id[subjectHits(gb_hits)],
    region   = "gene_body",
    stringsAsFactors = FALSE
  )
  
  pr_annot <- data.frame(
    peak_idx = queryHits(pr_hits),
    gene_id  = promoters_gr$gene_id[subjectHits(pr_hits)],
    region   = "promoter",
    stringsAsFactors = FALSE
  )
  
  dplyr::bind_rows(gb_annot, pr_annot) %>%
    dplyr::distinct(peak_idx, gene_id, region) %>%
    dplyr::mutate(
      peak_id         = peaks_gr$peak_id[peak_idx],
      fold_enrichment = peaks_gr$fold_enrichment[peak_idx]
    ) %>%
    dplyr::select(-peak_idx)
}

seqlevelsStyle(promoter_gr) <- "UCSC"
seqlevelsStyle(gene_body_gr) <- "UCSC"

message("Annotating DMSO peaks to genes...")
dmso_annot <- annotate_peaks(dmso_gr, gene_body_gr, promoter_gr)
message("Annotating KO peaks to genes...")
ko_annot   <- annotate_peaks(ko_gr,   gene_body_gr, promoter_gr)

cat("DMSO peak-gene pairs:", nrow(dmso_annot), "\n")
cat("KO peak-gene pairs:  ", nrow(ko_annot),   "\n")

## Per-gene peak summary: mean fold enrichment, peak count
summarise_peaks <- function(annot_df, condition) {
  annot_df %>%
    dplyr::group_by(gene_id, region) %>%
    dplyr::summarise(
      mean_fe  = mean(fold_enrichment, na.rm = TRUE),
      n_peaks  = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(condition = condition)
}

dmso_summary <- summarise_peaks(dmso_annot, "DMSO")
ko_summary   <- summarise_peaks(ko_annot,   "KO")

## Join DMSO and KO per gene per region, compute LFC
peak_lfc_df <- dplyr::bind_rows(dmso_summary, ko_summary) %>%
  tidyr::pivot_wider(
    id_cols     = c(gene_id, region),
    names_from  = condition,
    values_from = c(mean_fe, n_peaks),
    values_fill = list(mean_fe = 0, n_peaks = 0)
  ) %>%
  dplyr::mutate(
    peak_lfc     = log2(mean_fe_KO + 0.01) - log2(mean_fe_DMSO + 0.01),
    peak_gained  = mean_fe_KO  > 0 & mean_fe_DMSO == 0,
    peak_lost    = mean_fe_DMSO > 0 & mean_fe_KO  == 0,
    peak_shared  = mean_fe_KO  > 0 & mean_fe_DMSO > 0
  ) %>%
  dplyr::left_join(gene_meta, by = "gene_id") %>%
  dplyr::filter(!is.na(genomic_ei_group), !is.na(expr_group))

cat("\nPeak LFC summary:\n")
print(peak_lfc_df %>%
        dplyr::group_by(region) %>%
        dplyr::summarise(n = dplyr::n(),
                         gained = sum(peak_gained),
                         lost   = sum(peak_lost),
                         shared = sum(peak_shared),
                         .groups = "drop"))

## --- Plot A3-1: Peak fold enrichment LFC by expression — gene body and promoter ---
p_peak_expr <- ggplot(
  peak_lfc_df %>% dplyr::filter(peak_shared),
  aes(x = expr_group, y = peak_lfc, fill = expr_group)
) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, width = 0.6,
               alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = EXPR_COLORS) +
  coord_cartesian(ylim = c(
    quantile(peak_lfc_df$peak_lfc, 0.01, na.rm = TRUE),
    quantile(peak_lfc_df$peak_lfc, 0.99, na.rm = TRUE)
  )) +
  facet_wrap(~ region, nrow = 1,
             labeller = labeller(region = c(gene_body = "Gene body peaks",
                                            promoter  = "Promoter peaks"))) +
  labs(x = NULL,
       y = "log2 FC peak fold enrichment\n(KO / DMSO)",
       title    = "Peak accessibility change by expression level",
       subtitle = "Shared peaks only (present in both conditions)") +
  base_theme + theme(legend.position = "none")

print(p_peak_expr)

## --- Plot A3-2: Peak LFC by genomic E/I group ---
p_peak_ei <- ggplot(
  peak_lfc_df %>% dplyr::filter(peak_shared),
  aes(x = genomic_ei_group, y = peak_lfc, fill = genomic_ei_group)
) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, width = 0.6,
               alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = EI_COLORS) +
  scale_x_discrete(labels = c("Low E/I"    = "Low E/I\n(intron-rich)",
                              "Medium E/I" = "Medium E/I",
                              "High E/I"   = "High E/I\n(exon-rich)")) +
  coord_cartesian(ylim = c(
    quantile(peak_lfc_df$peak_lfc, 0.01, na.rm = TRUE),
    quantile(peak_lfc_df$peak_lfc, 0.99, na.rm = TRUE)
  )) +
  facet_wrap(~ region, nrow = 1,
             labeller = labeller(region = c(gene_body = "Gene body peaks",
                                            promoter  = "Promoter peaks"))) +
  labs(x = NULL,
       y = "log2 FC peak fold enrichment\n(KO / DMSO)",
       title    = "Peak accessibility change by genomic E/I group",
       subtitle = "Shared peaks only (present in both conditions)") +
  base_theme + theme(legend.position = "none")

print(p_peak_ei)

## --- Plot A3-3: Peak count change (gained vs lost) by expression ---
peak_change_df <- peak_lfc_df %>%
  dplyr::filter(region == "promoter") %>%
  dplyr::group_by(expr_group) %>%
  dplyr::summarise(
    gained = sum(peak_gained),
    lost   = sum(peak_lost),
    shared = sum(peak_shared),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(gained, lost, shared),
                      names_to = "status", values_to = "n") %>%
  dplyr::mutate(
    status = factor(status, levels = c("gained","shared","lost")),
    pct    = n / sum(n) * 100
  )

p_peak_change <- ggplot(peak_change_df,
                        aes(x = expr_group, y = n, fill = status)) +
  geom_col(position = "fill", width = 0.7) +
  scale_fill_manual(values = c(gained = "#D55E00",
                               shared = "grey70",
                               lost   = "#0072B2")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL,
       y = "Fraction of genes",
       title    = "Promoter peak gain/loss by expression level",
       subtitle = "Gained = KO only  |  Lost = DMSO only  |  Shared = both") +
  base_theme

print(p_peak_change)

message("\nATAC-seq accessibility analysis complete.")

library(dplyr)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)

##--------------------------------------------------------------## Load EU-seq RPKM

RPKM_FILE <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline/counts/euseq_rpkm_analysis.txt"
WINDOW_BP  <- 500000   ## 500kb neighborhood

rpkm_raw <- read_tsv(RPKM_FILE, show_col_types = FALSE) %>%
  dplyr::mutate(
    mean_polr2a = rowMeans(dplyr::select(., POLR2A_1, POLR2A_2, POLR2A_3)),
    mean_wt     = rowMeans(dplyr::select(., WT_1, WT_2, WT_3))
  ) %>%
  dplyr::select(gene_id, mean_polr2a, mean_wt)

##--------------------------------------------------------------## Compute neighborhood EU-seq signal per gene
## For each gene: sum mean WT RPKM of all genes within 500kb window
## Use gene midpoints for distance calculation

## Get gene midpoints from gene_body_gr (already built in script)
gene_coords <- data.frame(
  gene_id = gene_body_gr$gene_id,
  chr     = as.character(seqnames(gene_body_gr)),
  start   = start(gene_body_gr),
  end     = end(gene_body_gr),
  midpoint = (start(gene_body_gr) + end(gene_body_gr)) / 2,
  stringsAsFactors = FALSE
) %>%
  dplyr::inner_join(rpkm_raw, by = "gene_id") %>%
  dplyr::filter(!is.na(mean_wt))

## Build neighborhood windows
neighbor_gr <- GRanges(
  seqnames = gene_coords$chr,
  ranges   = IRanges(
    start = pmax(1, gene_coords$midpoint - WINDOW_BP / 2),
    end   = gene_coords$midpoint + WINDOW_BP / 2
  ),
  gene_id  = gene_coords$gene_id
)

## Gene midpoint GRanges for counting neighbors
midpoint_gr <- GRanges(
  seqnames = gene_coords$chr,
  ranges   = IRanges(
    start = gene_coords$midpoint,
    end   = gene_coords$midpoint
  ),
  gene_id    = gene_coords$gene_id,
  mean_wt    = gene_coords$mean_wt,
  mean_polr2a = gene_coords$mean_polr2a
)

## For each gene, sum WT RPKM of neighbors within 500kb
message("Computing neighborhood EU-seq signal...")
hits <- findOverlaps(neighbor_gr, midpoint_gr, ignore.strand = TRUE)

neighborhood_signal <- data.frame(
  gene_id              = neighbor_gr$gene_id[queryHits(hits)],
  neighbor_gene_id     = midpoint_gr$gene_id[subjectHits(hits)],
  neighbor_wt_rpkm     = midpoint_gr$mean_wt[subjectHits(hits)],
  neighbor_polr2a_rpkm = midpoint_gr$mean_polr2a[subjectHits(hits)]
) %>%
  dplyr::filter(gene_id != neighbor_gene_id) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    neighborhood_wt_rpkm     = sum(neighbor_wt_rpkm,     na.rm = TRUE),
    neighborhood_polr2a_rpkm = sum(neighbor_polr2a_rpkm, na.rm = TRUE),
    n_neighbors              = n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ## LFC of neighborhood transcription — did the neighborhood
    ## lose transcriptional activity after Pol II degradation?
    neighborhood_lfc = log2(neighborhood_polr2a_rpkm + 0.1) -
      log2(neighborhood_wt_rpkm     + 0.1),
    log10_nbhd_wt    = log10(neighborhood_wt_rpkm     + 0.1),
    log10_nbhd_ko    = log10(neighborhood_polr2a_rpkm + 0.1)
  )

explore_df <- gb_df %>%
  dplyr::select(gene_id, atac_lfc, expr_group, genomic_ei_group) %>%
  dplyr::inner_join(neighborhood_signal, by = "gene_id") %>%
  dplyr::mutate(
    n_neighbor_group = dplyr::case_when(
      n_neighbors < quantile(n_neighbors, 0.33) ~ "Low density",
      n_neighbors < quantile(n_neighbors, 0.67) ~ "Medium density",
      TRUE                                       ~ "High density"
    ),
    n_neighbor_group = factor(n_neighbor_group,
                              levels = c("Low density",
                                         "Medium density",
                                         "High density"))
  )

cat("Genes with neighborhood signal:", nrow(neighborhood_signal), "\n")
print(summary(neighborhood_signal$neighborhood_wt_rpkm))
print(summary(neighborhood_signal$n_neighbors))

##--------------------------------------------------------------## Join with ATAC LFC

explore_df <- gb_df %>%
  dplyr::select(gene_id, atac_lfc, expr_group, genomic_ei_group) %>%
  dplyr::inner_join(neighborhood_signal, by = "gene_id") %>%
  dplyr::mutate(
    log10_neighborhood = log10(neighborhood_wt_rpkm + 0.1),
    n_neighbor_group   = dplyr::case_when(
      n_neighbors < quantile(n_neighbors, 0.33) ~ "Low density",
      n_neighbors < quantile(n_neighbors, 0.67) ~ "Medium density",
      TRUE                                       ~ "High density"
    ),
    n_neighbor_group = factor(n_neighbor_group,
                              levels = c("Low density",
                                         "Medium density",
                                         "High density"))
  )

##--------------------------------------------------------------## Exploratory plots

## 1. Simple scatter: neighborhood transcription vs ATAC LFC
p1 <- ggplot(explore_df,
             aes(x = log10_neighborhood, y = atac_lfc)) +
  geom_point(size = 0.4, alpha = 0.3, color = "grey50") +
  geom_smooth(method = "loess", span = 0.3,
              color = "#D55E00", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("r = ", round(cor(explore_df$log10_neighborhood,
                                            explore_df$atac_lfc,
                                            use = "complete.obs"), 3)),
           hjust = 1.1, vjust = 1.5, size = 4, color = "grey30",
           fontface = "italic") +
  labs(x     = "log10(neighborhood WT EU-seq RPKM, 500kb)",
       y     = "ATAC log2 FC gene body (KO / DMSO)",
       title = "Gene body accessibility vs neighborhood transcriptional activity") +
  base_theme

print(p1)

## 2. Bin neighborhood signal and show mean ATAC LFC per bin
bin_df <- explore_df %>%
  dplyr::mutate(
    nbhd_bin = cut(log10_neighborhood,
                   breaks = 20, include.lowest = TRUE)
  ) %>%
  dplyr::group_by(nbhd_bin) %>%
  dplyr::summarise(
    mean_atac = mean(atac_lfc, na.rm = TRUE),
    se_atac   = sd(atac_lfc, na.rm = TRUE) / sqrt(n()),
    n         = n(),
    mid_x     = mean(log10_neighborhood, na.rm = TRUE),
    .groups   = "drop"
  )

p2 <- ggplot(bin_df, aes(x = mid_x, y = mean_atac)) +
  geom_ribbon(aes(ymin = mean_atac - se_atac,
                  ymax = mean_atac + se_atac),
              fill = "#D55E00", alpha = 0.2) +
  geom_line(color = "#D55E00", linewidth = 1) +
  geom_point(aes(size = n), color = "#D55E00") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_size_continuous(range = c(2, 6)) +
  labs(x     = "log10(neighborhood WT EU-seq RPKM, 500kb)",
       y     = "Mean ATAC log2 FC ± SE",
       size  = "n genes",
       title = "Mean accessibility change by neighborhood transcription level") +
  base_theme

print(p2)

## 3. Gene density (n_neighbors) vs ATAC LFC — rainplot
p3 <- ggplot(explore_df,
             aes(x = n_neighbor_group, y = atac_lfc,
                 fill = n_neighbor_group, color = n_neighbor_group)) +
  ggdist::stat_halfeye(adjust = 0.5, width = 0.5,
                       justification = -0.05, .width = 0,
                       point_colour = NA, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA,
               color = "grey30", linewidth = 0.5, alpha = 0.8,
               position = position_nudge(x = -0.12)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values  = c("Low density"    = "#0072B2",
                                "Medium density" = "#E69F00",
                                "High density"   = "#D55E00")) +
  scale_color_manual(values = c("Low density"    = "#0072B2",
                                "Medium density" = "#E69F00",
                                "High density"   = "#D55E00")) +
  coord_cartesian(ylim = c(
    quantile(explore_df$atac_lfc, 0.01, na.rm = TRUE),
    quantile(explore_df$atac_lfc, 0.99, na.rm = TRUE) * 1.25
  )) +
  labs(x = NULL, y = "ATAC log2 FC gene body (KO / DMSO)",
       title    = "Gene body accessibility by gene density",
       subtitle = "Gene density = n expressed neighbors within 500kb") +
  base_theme + theme(legend.position = "none")

print(p3)

## 4. Neighborhood signal stratified by E/I group
p4 <- ggplot(explore_df,
             aes(x = log10_neighborhood, y = atac_lfc,
                 color = genomic_ei_group)) +
  geom_smooth(method = "loess", span = 0.4,
              linewidth = 1, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_color_manual(values = EI_COLORS) +
  labs(x     = "log10(neighborhood WT EU-seq RPKM, 500kb)",
       y     = "ATAC log2 FC gene body (KO / DMSO)",
       title = "Accessibility vs neighborhood transcription by E/I group") +
  base_theme

print(p4)

## 5. Neighborhood transcription LFC vs ATAC LFC
## Does loss of neighborhood transcription predict chromatin opening?
r_nbhd_lfc <- cor(explore_df$neighborhood_lfc,
                  explore_df$atac_lfc, use = "complete.obs")

p5 <- ggplot(explore_df,
             aes(x = neighborhood_lfc, y = atac_lfc)) +
  geom_point(size = 0.4, alpha = 0.3, color = "grey50") +
  geom_smooth(method = "loess", span = 0.3,
              color = "#7B2D8B", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("r = ", round(r_nbhd_lfc, 3)),
           hjust = 1.1, vjust = 1.5, size = 4,
           color = "grey30", fontface = "italic") +
  labs(x     = "Neighborhood EU-seq log2 FC (POLR2A / WT, 500kb)",
       y     = "ATAC log2 FC gene body (KO / DMSO)",
       title = "Accessibility change vs neighborhood transcription change",
       subtitle = "Does loss of neighborhood Pol II predict chromatin opening?") +
  base_theme

print(p5)

## 6. Both conditions side by side — loess of neighborhood
## signal vs ATAC LFC for WT baseline vs KO baseline
long_nbhd_df <- explore_df %>%
  tidyr::pivot_longer(
    cols      = c(log10_nbhd_wt, log10_nbhd_ko),
    names_to  = "condition",
    values_to = "log10_nbhd"
  ) %>%
  dplyr::mutate(
    condition = dplyr::recode(condition,
                              log10_nbhd_wt = "WT",
                              log10_nbhd_ko = "POLR2A")
  )

p6 <- ggplot(long_nbhd_df,
             aes(x = log10_nbhd, y = atac_lfc,
                 color = condition)) +
  geom_smooth(method = "loess", span = 0.3,
              linewidth = 1, se = TRUE, alpha = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x     = "log10(neighborhood EU-seq RPKM, 500kb)",
       y     = "ATAC log2 FC gene body (KO / DMSO)",
       title = "Accessibility change vs neighborhood transcription",
       subtitle = "WT vs POLR2A neighborhood signal — both conditions") +
  base_theme

print(p6)

ATAC_LFC_THRESHOLD = 1

## Compute EU-seq LFC per gene
euseq_lfc_df <- rpkm_raw %>%
  dplyr::mutate(
    euseq_lfc = log2(mean_polr2a + 0.1) - log2(mean_wt + 0.1)
  ) %>%
  dplyr::select(gene_id, euseq_lfc, mean_wt, mean_polr2a)

## Join with gb_df
gb_eu_df <- gb_df %>%
  dplyr::inner_join(euseq_lfc_df, by = "gene_id") %>%
  dplyr::filter(!is.na(euseq_lfc),
                mean_wt >= 0.5) %>%   ## expressed in WT EU-seq
  dplyr::mutate(
    quadrant = dplyr::case_when(
      atac_lfc >  ATAC_LFC_THRESHOLD & euseq_lfc >  LFC_THRESHOLD ~ "Open & more transcribed",
      atac_lfc < -ATAC_LFC_THRESHOLD & euseq_lfc < -LFC_THRESHOLD ~ "Closed & less transcribed",
      atac_lfc >  ATAC_LFC_THRESHOLD & euseq_lfc < -LFC_THRESHOLD ~ "Open & less transcribed",
      atac_lfc < -ATAC_LFC_THRESHOLD & euseq_lfc >  LFC_THRESHOLD ~ "Closed & more transcribed",
      TRUE                                                          ~ "No change"
    ),
    quadrant = factor(quadrant,
                      levels = c("Open & more transcribed",
                                 "Closed & less transcribed",
                                 "Open & less transcribed",
                                 "Closed & more transcribed",
                                 "No change"))
  )

x_lims_eu <- quantile(gb_eu_df$euseq_lfc, c(0.005, 0.995), na.rm = TRUE)
y_lims_eu <- quantile(gb_eu_df$atac_lfc,  c(0.005, 0.995), na.rm = TRUE)

quad_colors <- c(
  "Open & more transcribed"   = "#D55E00",
  "Closed & less transcribed" = "#0072B2",
  "Open & less transcribed"   = "#E69F00",
  "Closed & more transcribed" = "#56B4E9",
  "No change"                 = "grey80"
)

quad_labels_eu <- gb_eu_df %>%
  dplyr::count(quadrant) %>%
  dplyr::mutate(
    pct   = round(n / sum(n) * 100, 1),
    label = paste0("n=", scales::comma(n), "\n(", pct, "%)"),
    x_pos = dplyr::case_when(
      quadrant == "Open & more transcribed"   ~ x_lims_eu[2] * 0.85,
      quadrant == "Closed & less transcribed" ~ x_lims_eu[1] * 0.85,
      quadrant == "Open & less transcribed"   ~ x_lims_eu[1] * 0.85,
      quadrant == "Closed & more transcribed" ~ x_lims_eu[2] * 0.85,
      TRUE ~ NA_real_
    ),
    y_pos = dplyr::case_when(
      quadrant == "Open & more transcribed"   ~ y_lims_eu[2] * 0.85,
      quadrant == "Closed & less transcribed" ~ y_lims_eu[1] * 0.85,
      quadrant == "Open & less transcribed"   ~ y_lims_eu[2] * 0.85,
      quadrant == "Closed & more transcribed" ~ y_lims_eu[1] * 0.85,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(x_pos))

marginal_eu_df <- gb_eu_df %>% dplyr::filter(quadrant != "No change")
r_atac_eu <- cor(gb_eu_df$atac_lfc, gb_eu_df$euseq_lfc, use = "complete.obs")

p_gb_scatter_eu <- ggplot(
  gb_eu_df %>% dplyr::arrange(quadrant == "No change"),
  aes(x = euseq_lfc, y = atac_lfc, color = quadrant)
) +
  geom_point(size = 0.6, alpha = 0.5) +
  ggside::geom_xsideboxplot(data = marginal_eu_df,
                            aes(y = euseq_lfc, fill = quadrant),
                            orientation = "y", outlier.shape = NA, linewidth = 0.4, width = 0.6) +
  ggside::geom_ysideboxplot(data = marginal_eu_df,
                            aes(x = atac_lfc, fill = quadrant),
                            orientation = "x", outlier.shape = NA, linewidth = 0.4, width = 0.6) +
  geom_hline(yintercept = c(-ATAC_LFC_THRESHOLD, ATAC_LFC_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_text(data = quad_labels_eu,
            aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3.2, color = "grey20", hjust = 0.5) +
  scale_color_manual(values = quad_colors) +
  scale_fill_manual(values  = quad_colors) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(r_atac_eu, 3)),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "grey40",
           fontface = "italic") +
  ggside::scale_xsidey_continuous(breaks = NULL) +
  ggside::scale_ysidex_continuous(breaks = NULL) +
  ggside::theme_ggside_void() +
  coord_cartesian(xlim = x_lims_eu, ylim = y_lims_eu) +
  labs(x        = "EU-seq log2 FC (POLR2A / WT)",
       y        = "ATAC log2 FC gene body (KO / DMSO)",
       title    = "Gene body accessibility vs nascent transcriptional change",
       subtitle = paste0("Colored = both ATAC |LFC| > ", ATAC_LFC_THRESHOLD,
                         " AND EU-seq |LFC| > ", LFC_THRESHOLD,
                         "  |  r = ", round(r_atac_eu, 3))) +
  base_theme +
  theme(legend.position    = "right",
        legend.text        = element_text(size = 9),
        ggside.panel.scale = 0.15)

print(p_gb_scatter_eu)
