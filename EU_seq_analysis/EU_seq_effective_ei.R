##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Compares genomic E/I ratio (reference) to nascent effective
# E/I ratio (exon RPKM / intron RPKM) for WT and POLR2A
# degraded conditions, stratified by gene expression level (WT TPM).
#
# Three lines per panel:
#   Grey dashed — genomic E/I (exon bp / intron bp from GTF)
#   Blue solid  — WT effective E/I (nascent exon/intron RPKM)
#   Rose solid  — POLR2A effective E/I
#
# Expression is binned into tertiles from WT TPM:
#   Low/Unexpressed (bottom TPM tertile + TPM==0/NA)
#   Medium          (middle TPM tertile)
#   High            (top TPM tertile)
# Medium is computed and reported in summary tables, but only
# Low/Unexpressed and High are shown in the plot.
#
# Each panel = density curve on top, outline-only boxplot strip
# (Genomic / WT / POLR2A) underneath, sharing the same x-axis.
#
# X-axis: log10(E/I ratio)
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(BiocGenerics)
})

##--------------------------------------------------------------## 1. Parameters

CNT_DIR    <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline"
INTRON_MAT <- file.path(CNT_DIR, "intron_dog/intron_counts/euseq_intron_rpkm.txt")
GENE_MAT   <- file.path(CNT_DIR, "counts/euseq_rpkm_analysis.txt")
GTF_PATH   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
EXPR_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"

## NOTE: RPKM_THRESHOLD is now UNUSED. It previously drove a gene_wt
## (gene-body RPKM) filter applied uniformly to both WT and POLR2A before
## computing effective E/I. That filter has been replaced by the
## exon-RPKM-is-zero-gets-removed rule (see section 4), applied
## independently per condition, per the pseudocount-methodology sync with
## the RIP-seq script. Left defined here in case anything downstream
## outside this script still expects it.
RPKM_THRESHOLD <- 0.5

POLR2A_SAMPLES <- c("POLR2A_1","POLR2A_2","POLR2A_3")
WT_SAMPLES     <- c("WT_1","WT_2","WT_3")
ALL_SAMPLES    <- c(POLR2A_SAMPLES, WT_SAMPLES)

CONDITION_COLORS <- c(POLR2A = "#CC6677", WT = "#88CCEE")

TOTAL_RNASEQ_SF <- c(
  POLR2A_1 = 0.747437, POLR2A_2 = 0.747437, POLR2A_3 = 0.747437,
  WT_1     = 1.000000, WT_2     = 1.000000, WT_3     = 1.000000
)

base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    strip.text       = element_text(face = "bold", size = 11)
  )

##--------------------------------------------------------------## 2. Load matrices and scale

message("Loading matrices...")

to_matrix <- function(df) {
  mat <- as.matrix(df[, -1])
  rownames(mat) <- sub("\\..*", "", df[[1]])
  mode(mat) <- "numeric"
  mat
}

apply_sf <- function(mat, sf) {
  sweep(mat, 2, sf[colnames(mat)], FUN = "*")
}

intron_mat    <- to_matrix(read_tsv(INTRON_MAT, show_col_types = FALSE))
gene_mat      <- to_matrix(read_tsv(GENE_MAT,   show_col_types = FALSE))

intron_scaled <- apply_sf(intron_mat, TOTAL_RNASEQ_SF)
gene_scaled   <- apply_sf(gene_mat,   TOTAL_RNASEQ_SF)

## Exonic RPKM = gene body - intron, floor at zero
common_gi   <- intersect(rownames(gene_scaled), rownames(intron_scaled))
exon_scaled <- gene_scaled[common_gi, ] - intron_scaled[common_gi, ]
exon_scaled[exon_scaled < 0] <- 0

cat("Matrices loaded and scaled.\n")
cat("Genes with both gene body and intron signal:", length(common_gi), "\n")

##--------------------------------------------------------------## 3. Compute genomic E/I (GTF) and expression bins (WT TPM)

message("Computing genomic E/I reference from GTF...")

txdb          <- GenomicFeatures::makeTxDbFromGFF(GTF_PATH, format = "gtf",
                                                  taxonomyId = 9606)
exons_by_gene <- GenomicRanges::reduce(
  GenomicFeatures::exonsBy(txdb, by = "gene"))
exon_lengths  <- sum(BiocGenerics::width(exons_by_gene))
genes_gr      <- GenomicFeatures::genes(txdb)
gene_lengths  <- setNames(BiocGenerics::width(genes_gr), genes_gr$gene_id)

gtf_common    <- intersect(names(exon_lengths), names(gene_lengths))
intron_len    <- pmax(gene_lengths[gtf_common] - exon_lengths[gtf_common], 0)
ei_genomic    <- exon_lengths[gtf_common] / (intron_len + 1)

genomic_df <- data.frame(
  gene_id       = sub("\\..*", "", gtf_common),
  intron_len_bp = as.numeric(intron_len),
  ei_genomic    = as.numeric(ei_genomic),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(ei_genomic), ei_genomic > 0, intron_len_bp >= 0) %>%
  ## Categorize by total intron length — kept alongside (not instead of)
  ## the WT-TPM expression bin below, so both stratifications are available
  dplyr::mutate(
    intron_cat = dplyr::case_when(
      intron_len_bp <  1000  ~ "Short (<1kb)",
      intron_len_bp <= 10000 ~ "Average (1-10kb)",
      TRUE                   ~ "Long (>10kb)"
    ),
    intron_cat = factor(intron_cat,
                        levels = c("Short (<1kb)",
                                   "Average (1-10kb)",
                                   "Long (>10kb)"))
  )

cat("\nGenomic E/I reference computed for:", nrow(genomic_df), "genes\n")
cat("Genomic E/I summary (overall):\n")
print(summary(genomic_df$ei_genomic))
cat("\nIntron length categories:\n")
print(table(genomic_df$intron_cat))

message("Loading WT TPM and computing expression tertiles...")

tpm_df <- read_tsv(EXPR_FILE, show_col_types = FALSE) %>%
  dplyr::select(gene_id, TPM) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id))

expressed_tpm <- tpm_df %>% dplyr::filter(!is.na(TPM), TPM > 0) %>% dplyr::pull(TPM)
tpm_medium <- quantile(expressed_tpm, 1/3, na.rm = TRUE)
tpm_high   <- quantile(expressed_tpm, 2/3, na.rm = TRUE)

tpm_df <- tpm_df %>%
  dplyr::mutate(
    expr_bin = dplyr::case_when(
      is.na(TPM) | TPM == 0 ~ "Low/Unexpressed",
      TPM <= tpm_medium     ~ "Low/Unexpressed",
      TPM <= tpm_high       ~ "Medium",
      TRUE                  ~ "High"
    ),
    expr_bin = factor(expr_bin, levels = c("Low/Unexpressed","Medium","High"))
  )

cat("\nExpression bins (all genes with TPM data):\n")
print(table(tpm_df$expr_bin))

##--------------------------------------------------------------## 4. Compute effective E/I per gene, anchored on the WT baseline
## ONE gene set, defined by WT (baseline) exon signal: exon_wt > 0.
## That same gene set is used to compute BOTH ei_wt and ei_polr2a — a
## gene qualifies via its WT baseline, since WT is the reference state
## a change is measured against. There's no equivalent reason to also
## require POLR2A exon signal independently.
##
## If POLR2A exon RPKM is exactly 0 for a WT-qualifying gene, the gene
## is KEPT (not dropped) and ei_polr2a is computed as exactly 0 — total
## loss of exonic signal in the treatment condition is a real result,
## not something to discard. These rows can't appear on the log10
## plots below (log10(0) is undefined) and get excluded there via the
## existing ei_ratio>0 filters, but they're counted explicitly here
## first so that exclusion is never silent.
##
## Intron pseudocount is EMPIRICAL: the minimum non-zero intron RPKM
## actually observed in this dataset, pooled across WT and POLR2A —
## not a fixed value. Intron RPKM==0 rows are never removed for being
## zero; "no observed intron signal" is itself a real finding.

message("Computing effective E/I ratios (WT-anchored)...")

ei_common <- intersect(common_gi, genomic_df$gene_id)

raw_ei_data <- data.frame(
  gene_id       = ei_common,
  intron_wt     = rowMeans(intron_scaled[ei_common, WT_SAMPLES]),
  exon_wt       = rowMeans(exon_scaled[ei_common,   WT_SAMPLES]),
  intron_polr2a = rowMeans(intron_scaled[ei_common, POLR2A_SAMPLES]),
  exon_polr2a   = rowMeans(exon_scaled[ei_common,   POLR2A_SAMPLES]),
  stringsAsFactors = FALSE
)

## Empirical intron pseudocount — pooled across WT and POLR2A
intron_values_pooled <- c(raw_ei_data$intron_wt, raw_ei_data$intron_polr2a)
intron_values_pooled <- intron_values_pooled[!is.na(intron_values_pooled) & intron_values_pooled > 0]
INTRON_PSEUDOCOUNT <- min(intron_values_pooled)

cat("\nIntron pseudocount (empirical minimum non-zero intron RPKM,",
    "pooled WT + POLR2A):", INTRON_PSEUDOCOUNT, "\n")
cat("  (", length(intron_values_pooled), "non-zero intron RPKM values contributed to this minimum )\n")

## WT-anchored gene set — the single baseline filter, applied once
ei_df <- raw_ei_data %>%
  dplyr::filter(exon_wt > 0) %>%
  dplyr::mutate(
    ei_wt     = exon_wt     / (intron_wt     + INTRON_PSEUDOCOUNT),
    ei_polr2a = exon_polr2a / (intron_polr2a + INTRON_PSEUDOCOUNT)
  ) %>%
  dplyr::left_join(genomic_df %>% dplyr::select(gene_id, ei_genomic, intron_cat), by = "gene_id") %>%
  dplyr::left_join(tpm_df     %>% dplyr::select(gene_id, expr_bin),               by = "gene_id") %>%
  dplyr::filter(!is.na(ei_genomic), !is.na(expr_bin), !is.na(intron_cat))

cat("\nWT-anchored gene set (exon WT RPKM > 0):", nrow(ei_df), "genes\n")

## CHECK: how many WT-anchored genes have ZERO POLR2A exon signal?
n_polr2a_exon_zero <- sum(ei_df$exon_polr2a == 0)
cat("Of these, genes with ZERO POLR2A exon signal (kept in ei_df with",
    "ei_polr2a=0 exactly; excluded from log10 plots below, not from the data):",
    n_polr2a_exon_zero, " (", round(100 * n_polr2a_exon_zero / nrow(ei_df), 1), "% )\n")

cat("\nGenes per expression bin:\n"); print(table(ei_df$expr_bin))
cat("Genes per intron length category:\n"); print(table(ei_df$intron_cat))

##--------------------------------------------------------------## 5. Build long format for plotting
## Three sources per gene per panel:
##   Reference  — genomic E/I (grey dashed)
##   WT         — effective E/I WT (blue solid)
##   POLR2A     — effective E/I POLR2A (rose solid)
## Medium expression is computed/counted above but intentionally excluded
## here — only Low/Unexpressed and High are plotted.

plot_df <- dplyr::bind_rows(
  ## Genomic reference — same WT-anchored gene set as WT/POLR2A below
  ei_df %>%
    dplyr::select(gene_id, expr_bin, ei_ratio = ei_genomic) %>%
    dplyr::mutate(source = "Genomic\n(reference)"),
  
  ## WT effective
  ei_df %>%
    dplyr::select(gene_id, expr_bin, ei_ratio = ei_wt) %>%
    dplyr::mutate(source = "WT\n(effective)"),
  
  ## POLR2A effective — rows with ei_ratio==0 (zero POLR2A exon signal)
  ## are dropped by the ei_ratio>0 filter below; they're real, kept in
  ## ei_df, just excluded from THIS log10 plot specifically (see the
  ## explicit count printed in section 4)
  ei_df %>%
    dplyr::select(gene_id, expr_bin, ei_ratio = ei_polr2a) %>%
    dplyr::mutate(source = "POLR2A\n(effective)")
) %>%
  dplyr::filter(ei_ratio > 0, !is.na(expr_bin), expr_bin %in% c("Low/Unexpressed","High")) %>%
  dplyr::mutate(
    log10_ei = log10(ei_ratio),
    source   = factor(source,
                      levels = c("Genomic\n(reference)",
                                 "WT\n(effective)",
                                 "POLR2A\n(effective)")),
    expr_bin = factor(expr_bin, levels = c("Low/Unexpressed","High"))
  )

cat("\nPlot data dimensions (Low/Unexpressed + High only):", nrow(plot_df), "rows\n")

##--------------------------------------------------------------## 6. Plot — density + outline-box strip under each expression panel

message("Plotting genomic vs effective E/I ratio by expression group...")

source_colors <- c(
  "Genomic\n(reference)" = "grey50",
  "WT\n(effective)"      = "#88CCEE",
  "POLR2A\n(effective)"  = "#CC6677"
)

source_linetypes <- c(
  "Genomic\n(reference)" = "dashed",
  "WT\n(effective)"      = "solid",
  "POLR2A\n(effective)"  = "solid"
)

source_levels <- c("Genomic\n(reference)", "WT\n(effective)", "POLR2A\n(effective)")

build_ei_panel_generic <- function(sub_df, panel_title, x_range) {
  p_density <- ggplot(sub_df, aes(x = log10_ei, color = source, linetype = source)) +
    geom_density(linewidth = 0.9, adjust = 1.2) +
    geom_vline(xintercept = 0, color = "grey70", linetype = "dotted", linewidth = 0.4) +
    scale_color_manual(values = source_colors, name = NULL) +
    scale_linetype_manual(values = source_linetypes, name = NULL) +
    scale_x_continuous(limits = x_range) +
    labs(x = NULL, y = "Density", title = panel_title) +
    base_theme +
    theme(
      legend.position = "top",
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.margin     = margin(b = 2)
    )
  
  p_strip <- ggplot(sub_df, aes(x = log10_ei, y = source, color = source)) +
    geom_boxplot(fill = NA, width = 0.5, outlier.shape = NA,
                 linewidth = 0.6, show.legend = FALSE) +
    scale_color_manual(values = source_colors) +
    scale_x_continuous(limits = x_range) +
    scale_y_discrete(limits = rev(source_levels)) +
    labs(x = "log10 E/I", y = NULL) +
    base_theme +
    theme(
      legend.position = "none",
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      plot.margin     = margin(t = 2)
    )
  
  p_density / p_strip + patchwork::plot_layout(heights = c(3, 1))
}

## Shared x-axis range across BOTH expression-bin panels (Low/Unexpressed,
## High) so the two facets are directly comparable, matching facet_wrap's
## scales="fixed" behavior rather than each panel autoscaling independently.
plot_df_x_range <- range(plot_df$log10_ei, na.rm = TRUE)

build_ei_panel <- function(bin_name) {
  sub_df <- plot_df %>% dplyr::filter(expr_bin == bin_name)
  build_ei_panel_generic(sub_df, bin_name, x_range = plot_df_x_range)
}

panels <- lapply(levels(plot_df$expr_bin), build_ei_panel)

final_plot <- (wrap_plots(panels, nrow = 1) + patchwork::plot_layout(guides = "collect")) +
  plot_annotation(
    title    = "Genomic vs nascent exon/intron ratio by expression group",
    subtitle = paste0(
      "Grey dashed = genomic E/I (exon bp / intron bp)  |  ",
      "Solid = nascent exon RPKM / intron RPKM  |  ",
      "Low/Unexpressed and High WT TPM tertiles only  |  ",
      "WT-anchored n = ", nrow(ei_df), " genes"
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(color = "grey40", size = 11))
  ) &
  theme(legend.position = "top")

print(final_plot)

## Summary statistics per expression bin and source (Low/Unexpressed + High,
## matching what's plotted — see ei_df above for all three bins)
summary_df <- plot_df %>%
  dplyr::group_by(expr_bin, source) %>%
  dplyr::summarise(
    median_log10_ei = median(log10_ei, na.rm = TRUE),
    mean_log10_ei   = mean(log10_ei,   na.rm = TRUE),
    n               = dplyr::n(),
    .groups         = "drop"
  )

cat("\nMedian log10(E/I) per expression bin and source:\n")
print(summary_df)

##--------------------------------------------------------------## 6b. Gene-length faceted plot — ALL genes, no expression filter
##   This is the original plot from before the script pivoted to
##   expression-based faceting (Genomic/WT/POLR2A by intron length,
##   Short/Average/Long, every gene included regardless of expr_bin).
##   It was dropped when expression faceting was added and never
##   restored — restored here as its own standalone plot, alongside
##   (not replacing) the expression-based plots in sections 6 and 7.

message("\nBuilding plain intron-length plot (all genes, no expression filter)...")

length_only_plot_df <- dplyr::bind_rows(
  ei_df %>%
    dplyr::select(gene_id, intron_cat, ei_ratio = ei_genomic) %>%
    dplyr::mutate(source = "Genomic\n(reference)"),
  ei_df %>%
    dplyr::select(gene_id, intron_cat, ei_ratio = ei_wt) %>%
    dplyr::mutate(source = "WT\n(effective)"),
  ei_df %>%
    dplyr::select(gene_id, intron_cat, ei_ratio = ei_polr2a) %>%
    dplyr::mutate(source = "POLR2A\n(effective)")
) %>%
  dplyr::filter(ei_ratio > 0) %>%
  dplyr::mutate(
    log10_ei = log10(ei_ratio),
    source   = factor(source, levels = c("Genomic\n(reference)",
                                         "WT\n(effective)",
                                         "POLR2A\n(effective)"))
  )

cat("Plain length-only plot data (WT-anchored):", nrow(length_only_plot_df), "rows,",
    dplyr::n_distinct(length_only_plot_df$gene_id), "genes\n")

length_only_x_range <- range(length_only_plot_df$log10_ei, na.rm = TRUE)

length_only_panels <- lapply(levels(length_only_plot_df$intron_cat), function(cat_name) {
  sub_df <- length_only_plot_df %>% dplyr::filter(intron_cat == cat_name)
  build_ei_panel_generic(sub_df, cat_name, x_range = length_only_x_range)
})

length_only_final <- (wrap_plots(length_only_panels, nrow = 1) + patchwork::plot_layout(guides = "collect")) +
  plot_annotation(
    title    = "Genomic vs nascent exon/intron ratio by intron length",
    subtitle = paste0(
      "Grey dashed = genomic E/I (exon bp / intron bp)  |  ",
      "Solid = nascent exon RPKM / intron RPKM  |  no expression filter  |  ",
      "WT-anchored n = ", nrow(ei_df), " genes"
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(color = "grey40", size = 11))
  ) &
  theme(legend.position = "top")

print(length_only_final)

##--------------------------------------------------------------## 7. Gene-length faceted plot, run twice: High expression, then Low
##   Same panel style as section 6 (density + outline-box strip), but now
##   faceted by intron_cat (Short/Average/Long) instead of expr_bin, with
##   the whole figure pre-filtered to one expression group at a time.

message("\nBuilding gene-length faceted plots, split by expression group...")

length_plot_df <- dplyr::bind_rows(
  ei_df %>%
    dplyr::select(gene_id, expr_bin, intron_cat, ei_ratio = ei_genomic) %>%
    dplyr::mutate(source = "Genomic\n(reference)"),
  ei_df %>%
    dplyr::select(gene_id, expr_bin, intron_cat, ei_ratio = ei_wt) %>%
    dplyr::mutate(source = "WT\n(effective)"),
  ei_df %>%
    dplyr::select(gene_id, expr_bin, intron_cat, ei_ratio = ei_polr2a) %>%
    dplyr::mutate(source = "POLR2A\n(effective)")
) %>%
  dplyr::filter(ei_ratio > 0, !is.na(expr_bin)) %>%
  dplyr::mutate(
    log10_ei = log10(ei_ratio),
    source   = factor(source, levels = source_levels)
  )

build_length_faceted_plot <- function(expr_filter, expr_label) {
  sub_all <- length_plot_df %>% dplyr::filter(expr_bin == expr_filter)
  
  ## Only build panels for intron_cat levels that actually have rows in this
  ## expression subset — an empty panel would break range() on empty data
  cat_levels <- levels(sub_all$intron_cat)
  cat_levels <- cat_levels[cat_levels %in% unique(as.character(sub_all$intron_cat))]
  if (length(cat_levels) < length(levels(sub_all$intron_cat))) {
    missing_cats <- setdiff(levels(sub_all$intron_cat), cat_levels)
    message("  Note: ", expr_label, " expression has no genes in ",
            paste(missing_cats, collapse = ", "), " — panel(s) skipped")
  }
  
  ## Shared x-axis range across all length panels WITHIN this figure (Short/
  ## Average/Long), so the three facets are directly comparable to each
  ## other. The High-expression figure and Low-expression figure are still
  ## separate calls, so each gets its own range unless you want those two
  ## unified as well — see note below the two calls.
  shared_x_range <- range(sub_all$log10_ei, na.rm = TRUE)
  
  panels <- lapply(cat_levels, function(cat_name) {
    sub_df <- sub_all %>% dplyr::filter(intron_cat == cat_name)
    build_ei_panel_generic(sub_df, cat_name, x_range = shared_x_range)
  })
  
  n_genes <- ei_df %>% dplyr::filter(expr_bin == expr_filter) %>% nrow()
  
  (wrap_plots(panels, nrow = 1) + patchwork::plot_layout(guides = "collect")) +
    plot_annotation(
      title    = paste0("Genomic vs nascent exon/intron ratio by intron length \u2014 ",
                        expr_label, " expression"),
      subtitle = paste0(
        "Grey dashed = genomic E/I (exon bp / intron bp)  |  ",
        "Solid = nascent exon RPKM / intron RPKM  |  ",
        expr_label, " WT TPM tertile only  |  ",
        "WT-anchored n = ", n_genes, " genes"
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(color = "grey40", size = 11))
    ) &
    theme(legend.position = "top")
}

message("\n── Gene length faceted plot: High expression only ──────────────────")
length_plot_high <- build_length_faceted_plot("High", "High")
print(length_plot_high)

message("\n── Gene length faceted plot: Low/Unexpressed only ───────────────────")
length_plot_low <- build_length_faceted_plot("Low/Unexpressed", "Low/Unexpressed")
print(length_plot_low)

## Summary statistics for both length-faceted plots
length_summary_df <- length_plot_df %>%
  dplyr::filter(expr_bin %in% c("High","Low/Unexpressed")) %>%
  dplyr::group_by(expr_bin, intron_cat, source) %>%
  dplyr::summarise(
    median_log10_ei = median(log10_ei, na.rm = TRUE),
    mean_log10_ei   = mean(log10_ei,   na.rm = TRUE),
    n               = dplyr::n(),
    .groups         = "drop"
  )

cat("\nMedian log10(E/I) per expression group, intron length category, and source:\n")
print(length_summary_df)

message("Plot complete.")