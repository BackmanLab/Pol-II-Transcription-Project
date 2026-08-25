##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# ATAC-seq exon/intron accessibility, stratified by local gene
# density, with Pol II KO delta tracks, a continuous dose-response
# check, and correlation against EU-seq nascent transcription change.
#
# Three density/context metrics, each producing a Q1 (low) vs Q4 (high)
# split (E/I ratio column added after this header was first written):
#   1. Gene-count density  — # expressed genes within 1Mb window
#   2. Expression density  — mean TPM of expressed genes in that window
#   3. E/I ratio            — genomic exon/intron length ratio (per gene, from GTF)
#
# Per interval type (exon, intron), for each density metric:
#   Top:    metagene signal (DMSO solid / KO dashed, colored by Q1/Q4)
#   Bottom: log2(KO/DMSO) delta track
#
# Plus:
#   - Per-gene log2(KO/DMSO) delta vs. continuous local gene density
#     (dose-response scan instead of just the Q1/Q4 extremes)
#   - Per-gene log2(KO/DMSO) delta vs. EU-seq nascent transcription
#     log2FC (Aux/WT), to see if accessibility loss tracks with how
#     much a gene's own transcription dropped
#
# The intergenic accessibility analysis (original sections 11-13)
# is carried over unchanged and is independent of everything above.
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(GenomeInfoDb)
  library(AnnotationDbi)
  library(scales)
})

##--------------------------------------------------------------## 1. Parameters

ATAC_DIR   <- "/projects/b1042/BackmanLab/Lucas/120525_pol2/data/atac"
DMSO_BW    <- file.path(ATAC_DIR, "pol2_dmso.bw")
KO_BW      <- file.path(ATAC_DIR, "pol2_ko.bw")

GTF_PATH   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
TPM_FILE   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/HCT116_RNAseq_quantification.tsv"

## RPKM table has per-replicate columns, not a precomputed log2FC — the lfc
## is built from it directly in section 14 below.
EU_RPKM_PATH     <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline/counts/euseq_rpkm_analysis.txt"
POLR2A_REP_COLS  <- c("POLR2A_1", "POLR2A_2", "POLR2A_3")
WT_REP_COLS      <- c("WT_1", "WT_2", "WT_3")

DENSITY_WINDOW <- 1e6     ## 1Mb window for local density metrics
N_BINS_EXON    <- 50      ## bins per exon
N_BINS_INTRON  <- 50      ## bins per intron
TPM_THRESHOLD  <- 1       ## minimum TPM for expressed genes
PSEUDOCOUNT    <- 0.01    ## added before log2(KO/DMSO) to avoid log(0); ~2-3
## orders of magnitude below the smallest observed
## exon/intron enrichment values

STD_CHRS <- paste0("chr", c(1:22, "X", "Y"))

## Facet strip labels for the new Q1/Q4 signal + delta plots
## Colors for the Q1/Q4 overlay in the delta panels
QUARTILE_COLORS <- c(Q1 = "#0072B2", Q4 = "#D55E00")
QUARTILE_LABELS <- c(Q1 = "Q1 (low)", Q4 = "Q4 (high)")

## Kept for the intergenic analysis at the end of the script (unchanged)
DENSITY_COLORS <- c("Q1 (Low density)"  = "#0072B2",
                    "Q4 (High density)" = "#D55E00")

base_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.spacing    = unit(0.8, "lines")
  )

metagene_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle    = element_text(color = "grey45", size = 9.5),
    axis.title       = element_text(size = 11),
    axis.title.y     = element_text(margin = margin(r = 8)),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey85", fill = NA, linewidth = 0.4),
    panel.spacing    = unit(1, "lines"),
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    strip.text       = element_text(face = "bold", size = 10),
    plot.margin      = margin(8, 15, 4, 15)
  )

##--------------------------------------------------------------## 2. Helper functions (new analyses)

## Vectorized per-interval tiling
tile_intervals_fast <- function(gr, n_bins) {
  gr <- gr[width(gr) >= n_bins]
  if (length(gr) == 0) return(NULL)
  tiled              <- GenomicRanges::tile(gr, n = n_bins)
  tiled_flat         <- unlist(tiled)
  tiled_flat$gene_id <- rep(gr$gene_id, each = n_bins)
  tiled_flat$bin_idx <- rep(seq_len(n_bins), times = length(gr))
  tiled_flat
}

## Quartile assignment, reused for both density metrics below
assign_density_quartiles <- function(df, value_col) {
  v <- df[[value_col]]
  q <- stats::quantile(v, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  df$density_quartile <- dplyr::case_when(
    v < q[1] ~ "Q1",
    v < q[2] ~ "Q2",
    v < q[3] ~ "Q3",
    TRUE     ~ "Q4"
  )
  df$density_quartile <- factor(df$density_quartile, levels = c("Q1", "Q2", "Q3", "Q4"))
  df
}

## Per-gene, per-bin DMSO/KO signal for one interval type, over ALL supplied
## genes at once (computed once; quartile grouping is applied afterward via
## summarise_profile(), so tile()/binnedAverage() never runs more than once
## per interval type regardless of how many groupings you slice it by).
compute_bin_signal <- function(intervals_gr, gene_ids, n_bins, dmso_rle, ko_rle, label = NULL) {
  
  gr_sub <- intervals_gr[intervals_gr$gene_id %in% gene_ids]
  gr_sub <- gr_sub[as.character(seqnames(gr_sub)) %in% STD_CHRS]
  gr_sub <- GenomeInfoDb::keepSeqlevels(gr_sub, value = STD_CHRS, pruning.mode = "coarse")
  
  if (length(gr_sub) == 0) {
    warning("No intervals found", if (!is.null(label)) paste0(" for: ", label))
    return(NULL)
  }
  
  if (!is.null(label)) {
    cat("  ", label, ": ", length(gr_sub), " intervals from ",
        length(unique(gr_sub$gene_id)), " genes\n", sep = "")
  }
  
  bins_gr <- tile_intervals_fast(gr_sub, n_bins)
  if (is.null(bins_gr) || length(bins_gr) == 0) return(NULL)
  
  common_chrs <- intersect(seqlevels(bins_gr), names(dmso_rle))
  bins_gr     <- GenomeInfoDb::keepSeqlevels(bins_gr, value = common_chrs, pruning.mode = "coarse")
  dmso_sub    <- dmso_rle[common_chrs]
  ko_sub      <- ko_rle[common_chrs]
  
  bins_gr$dmso <- GenomicRanges::binnedAverage(bins = bins_gr, numvar = dmso_sub, varname = "score")$score
  bins_gr$ko   <- GenomicRanges::binnedAverage(bins = bins_gr, numvar = ko_sub,   varname = "score")$score
  
  data.frame(
    gene_id = bins_gr$gene_id,
    bin_idx = bins_gr$bin_idx,
    dmso    = bins_gr$dmso,
    ko      = bins_gr$ko,
    stringsAsFactors = FALSE
  )
}

## Collapse a per-gene-per-bin signal table down to a Q1/Q4 metagene profile,
## given a gene_id -> density_group ("Q1"/"Q4") lookup table.
summarise_profile <- function(bins_df, group_df, n_bins) {
  bins_df %>%
    dplyr::inner_join(group_df, by = "gene_id") %>%
    dplyr::group_by(bin_idx, density_group) %>%
    dplyr::summarise(
      mean_dmso = mean(dmso, na.rm = TRUE),
      mean_ko   = mean(ko,   na.rm = TRUE),
      n_genes   = dplyr::n_distinct(gene_id),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      bin_pct       = (bin_idx - 0.5) / n_bins * 100,
      density_group = factor(density_group, levels = c("Q1", "Q4"))
    ) %>%
    tidyr::pivot_longer(
      cols      = c(mean_dmso, mean_ko),
      names_to  = "condition",
      values_to = "mean_signal"
    ) %>%
    dplyr::mutate(condition = dplyr::recode(condition, mean_dmso = "DMSO", mean_ko = "KO"))
}

## log2(KO/DMSO) delta track from a summarise_profile() output
compute_delta <- function(profile_df) {
  profile_df %>%
    dplyr::select(bin_idx, bin_pct, density_group, condition, mean_signal) %>%
    tidyr::pivot_wider(names_from = condition, values_from = mean_signal) %>%
    dplyr::mutate(log2_delta = log2((KO + PSEUDOCOUNT) / (DMSO + PSEUDOCOUNT)))
}

## Metagene signal panel: Q1 and Q4 as separate facet columns, color = condition
make_signal_plot <- function(profile_df, title) {
  
  n_df <- profile_df %>%
    dplyr::filter(condition == "DMSO") %>%
    dplyr::distinct(density_group, n_genes) %>%
    dplyr::mutate(label = paste0("n=", scales::comma(n_genes)), x_pos = 3, y_pos = Inf)
  
  ggplot(profile_df, aes(x = bin_pct, y = mean_signal,
                         color = condition, linetype = condition, group = condition)) +
    geom_line(linewidth = 0.9) +
    geom_text(data = n_df, aes(x = x_pos, y = y_pos, label = label),
              hjust = 0, vjust = 1.6, size = 3.4, color = "grey40", inherit.aes = FALSE) +
    scale_color_manual(values = c(DMSO = "#0072B2", KO = "#D55E00")) +
    scale_linetype_manual(values = c(DMSO = "solid", KO = "dashed"),
                          labels = c(DMSO = "DMSO (control)", KO = "KO (Pol II degradation)")) +
    scale_x_continuous(labels = function(x) paste0(x, "%"),
                       breaks = c(0, 25, 50, 75, 100), expand = c(0.02, 0)) +
    facet_wrap(~ density_group, nrow = 1, labeller = as_labeller(QUARTILE_LABELS)) +
    labs(x = NULL, y = "Mean enrichment", title = title, color = "Condition", linetype = "Condition") +
    metagene_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(8, 15, 10, 15))
}

## Delta panel: log2(KO/DMSO), Q1 vs Q4 overlaid on one axis (not faceted) so
## the comparison between density groups is direct, with a shaded fill down
## to zero for each line
make_delta_plot <- function(delta_df, x_label) {
  ggplot(delta_df, aes(x = bin_pct, y = log2_delta, color = density_group, fill = density_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_ribbon(aes(ymin = pmin(log2_delta, 0), ymax = pmax(log2_delta, 0)),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = QUARTILE_COLORS, labels = QUARTILE_LABELS) +
    scale_fill_manual(values = QUARTILE_COLORS, guide = "none") +
    scale_x_continuous(labels = function(x) paste0(x, "%"),
                       breaks = c(0, 25, 50, 75, 100), expand = c(0.02, 0)) +
    labs(x = x_label, y = expression(log[2]~(KO/DMSO)), color = "Gene density") +
    metagene_theme +
    theme(plot.margin = margin(12, 15, 4, 15))
}

## Continuous dose-response: per-gene delta vs. log10(local gene density).
## Points are binned means +/- SE (equal-N bins) so the trend is readable;
## the regression line and R^2/p-value are still fit on the full per-gene data.
plot_dose_response <- function(df, region_label, n_bins_dose = 20) {
  fit  <- lm(log2_delta ~ log10_neighbors, data = df)
  smry <- summary(fit)
  r2   <- smry$r.squared
  pval <- smry$coefficients["log10_neighbors", "Pr(>|t|)"]
  label_txt <- paste0("Linear fit: R\u00b2 = ", signif(r2, 3), "   ", scales::pvalue(pval, add_p = TRUE))
  
  binned_df <- df %>%
    dplyr::mutate(bin = dplyr::ntile(log10_neighbors, n_bins_dose)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      x_mean  = mean(log10_neighbors, na.rm = TRUE),
      y_mean  = mean(log2_delta, na.rm = TRUE),
      y_se    = sd(log2_delta, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )
  
  ggplot(binned_df, aes(x = x_mean, y = y_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.9,
                color = "#D55E00", fill = "#D55E00", alpha = 0.15, linewidth = 0.7) +
    geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se), width = 0, color = "grey35") +
    geom_point(color = "grey20", size = 2) +
    annotate("text", x = -Inf, y = Inf, label = label_txt, hjust = -0.05, vjust = 1.5, size = 3.6) +
    labs(x = "log10(neighboring expressed genes within 1Mb + 1)",
         y = expression(log[2]~(KO/DMSO)),
         title = region_label) +
    metagene_theme
}

## Correlation: per-gene accessibility delta vs. EU-seq nascent transcription delta.
## Same binned-mean treatment as plot_dose_response, for the same reason.
plot_eu_correlation <- function(df, region_label, n_bins_dose = 20) {
  fit  <- lm(log2_delta ~ eu_lfc, data = df)
  smry <- summary(fit)
  r2   <- smry$r.squared
  pval <- smry$coefficients["eu_lfc", "Pr(>|t|)"]
  label_txt <- paste0("Linear fit: R\u00b2 = ", signif(r2, 3), "   ", scales::pvalue(pval, add_p = TRUE))
  
  binned_df <- df %>%
    dplyr::mutate(bin = dplyr::ntile(eu_lfc, n_bins_dose)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      x_mean  = mean(eu_lfc, na.rm = TRUE),
      y_mean  = mean(log2_delta, na.rm = TRUE),
      y_se    = sd(log2_delta, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )
  
  ggplot(binned_df, aes(x = x_mean, y = y_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.9,
                color = "#0072B2", fill = "#0072B2", alpha = 0.15, linewidth = 0.7) +
    geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se), width = 0, color = "grey35") +
    geom_point(color = "grey20", size = 2) +
    annotate("text", x = -Inf, y = Inf, label = label_txt, hjust = -0.05, vjust = 1.5, size = 3.6) +
    labs(x = "Nascent transcription log2FC (EU-seq, Aux/WT)",
         y = expression(ATAC~log[2]~(KO/DMSO)),
         title = region_label) +
    metagene_theme
}

## Signal trend vs. a continuous local density metric: DMSO (solid) and KO
## (dashed) loess curves + confidence band, no raw points. Title is just the
## density metric (short, centered, so 3 columns don't collide); region
## (Exon/Intron) goes on the y-axis label instead, since that's rotated
## text that doesn't compete for horizontal space with neighboring panels.
plot_signal_vs_density <- function(df, x_col, metric_title, region_label) {
  df_long <- df %>%
    tidyr::pivot_longer(cols = c(mean_dmso, mean_ko), names_to = "condition", values_to = "signal") %>%
    dplyr::mutate(condition = dplyr::recode(condition, mean_dmso = "DMSO", mean_ko = "KO"))
  
  ggplot(df_long, aes(x = .data[[x_col]], y = signal,
                      color = condition, linetype = condition, fill = condition)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1, alpha = 0.15) +
    scale_color_manual(values = c(DMSO = "#0072B2", KO = "#D55E00")) +
    scale_fill_manual(values = c(DMSO = "#0072B2", KO = "#D55E00"), guide = "none") +
    scale_linetype_manual(values = c(DMSO = "solid", KO = "dashed"),
                          labels = c(DMSO = "DMSO (control)", KO = "KO (Pol II degradation)")) +
    labs(x = NULL, y = paste0(region_label, " \u2014 mean enrichment"), title = metric_title,
         color = "Condition", linetype = "Condition") +
    metagene_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title  = element_text(hjust = 0.5),
          plot.margin = margin(8, 12, 10, 15))
}

## Delta trend vs. the same continuous density metric: loess curve with the
## area shaded down to zero (not a symmetric SE ribbon), styled after the
## reference ΔS(d) panel — warm tan fill, rust line, minimal gridlines
plot_delta_vs_density <- function(df, x_col, x_label, n_grid = 200) {
  fit <- stats::loess(stats::as.formula(paste0("log2_delta ~ ", x_col)), data = df, span = 0.75)
  
  x_grid  <- seq(min(df[[x_col]], na.rm = TRUE), max(df[[x_col]], na.rm = TRUE), length.out = n_grid)
  pred_df <- setNames(data.frame(x_grid), x_col)
  pred_df$y_fit <- as.numeric(predict(fit, newdata = pred_df))
  
  ggplot(pred_df, aes(x = .data[[x_col]], y = y_fit)) +
    geom_ribbon(aes(ymin = pmin(y_fit, 0), ymax = pmax(y_fit, 0)), fill = "#F2C79E", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_line(color = "#C1591A", linewidth = 1.1) +
    labs(x = x_label, y = expression(log[2]~(KO/DMSO))) +
    metagene_theme +
    theme(panel.grid.major.x = element_blank(), panel.border = element_blank(),
          plot.margin = margin(12, 12, 4, 15))
}

## Give a set of plots a common y-axis range (reads back each plot's actual
## rendered y range, including geom_smooth/ribbon extent, then re-applies it
## to all of them via coord_cartesian so they zoom to the same window rather
## than recomputing stats)
apply_shared_ylim <- function(plots) {
  ranges <- lapply(plots, function(p) {
    ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
  })
  shared <- range(unlist(ranges))
  lapply(plots, function(p) p + coord_cartesian(ylim = shared))
}

##--------------------------------------------------------------## 3. Build gene annotations from GTF

message("Building gene annotations from GTF...")

txdb     <- GenomicFeatures::makeTxDbFromGFF(GTF_PATH,
                                             format     = "gtf",
                                             taxonomyId = 9606)
genes_gr <- GenomicFeatures::genes(txdb)
genes_gr$gene_id <- sub("\\..*", "", genes_gr$gene_id)
GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"

cat("Total genes in GTF:", length(genes_gr), "\n")

##--------------------------------------------------------------## 4. Compute exon intervals

message("Computing exon intervals...")

exons_by_gene      <- GenomicRanges::reduce(
  GenomicFeatures::exonsBy(txdb, by = "gene")
)
exons_flat         <- unlist(exons_by_gene)
GenomeInfoDb::seqlevelsStyle(exons_flat) <- "UCSC"
exons_flat$gene_id <- sub("\\..*", "", names(exons_flat))
names(exons_flat)  <- NULL

cat("Exon intervals:", length(exons_flat), "\n")

##--------------------------------------------------------------## 5. Compute intron intervals (fast via TxDb)

message("Computing intron intervals...")

introns_by_tx <- GenomicFeatures::intronsByTranscript(txdb,
                                                      use.names = TRUE)

tx2gene <- AnnotationDbi::select(
  txdb,
  keys    = names(introns_by_tx),
  columns = "GENEID",
  keytype = "TXNAME"
) %>%
  dplyr::filter(!is.na(GENEID))

introns_flat           <- unlist(introns_by_tx)
introns_flat$tx_name   <- names(introns_flat)
introns_flat$gene_id   <- tx2gene$GENEID[
  match(introns_flat$tx_name, tx2gene$TXNAME)]
introns_flat           <- introns_flat[!is.na(introns_flat$gene_id)]
introns_flat$gene_id   <- sub("\\..*", "", introns_flat$gene_id)
names(introns_flat)    <- NULL
GenomeInfoDb::seqlevelsStyle(introns_flat) <- "UCSC"

introns_gr         <- unlist(GenomicRanges::reduce(
  split(introns_flat, introns_flat$gene_id)
))
introns_gr$gene_id <- names(introns_gr)
names(introns_gr)  <- NULL

cat("Intron intervals:", length(introns_gr), "\n")

##--------------------------------------------------------------## 6. Load expression data

message("Loading TPM data...")

tpm_df <- read_tsv(TPM_FILE, show_col_types = FALSE) %>%
  dplyr::mutate(
    gene_id = sub(".*\\|(ENSG[0-9]+)\\.[0-9]+\\|.*", "\\1", target_id)
  ) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(tpm >= TPM_THRESHOLD)

expressed_ids <- tpm_df$gene_id
cat("Expressed genes (TPM >=", TPM_THRESHOLD, "):",
    length(expressed_ids), "\n")

##--------------------------------------------------------------## 7. Local density metrics (gene-count and expression-based)

message("Computing local density metrics (1Mb window)...")

expr_genes_gr <- genes_gr[genes_gr$gene_id %in% expressed_ids]

midpoints_gr <- GRanges(
  seqnames = seqnames(expr_genes_gr),
  ranges   = IRanges(
    start = (start(expr_genes_gr) + end(expr_genes_gr)) %/% 2,
    end   = (start(expr_genes_gr) + end(expr_genes_gr)) %/% 2
  ),
  gene_id  = expr_genes_gr$gene_id
)

window_gr <- GRanges(
  seqnames = seqnames(expr_genes_gr),
  ranges   = IRanges(
    start = pmax(1, (start(expr_genes_gr) + end(expr_genes_gr)) %/% 2 - DENSITY_WINDOW %/% 2),
    end   =          (start(expr_genes_gr) + end(expr_genes_gr)) %/% 2 + DENSITY_WINDOW %/% 2
  ),
  gene_id  = expr_genes_gr$gene_id
)

## Neighbor pairs within 1Mb, self excluded — reused for BOTH density metrics
hits <- findOverlaps(window_gr, midpoints_gr, ignore.strand = TRUE)

neighbor_pairs_df <- data.frame(
  gene_id          = window_gr$gene_id[queryHits(hits)],
  neighbor_gene_id = midpoints_gr$gene_id[subjectHits(hits)],
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(gene_id != neighbor_gene_id)

## Metric 1: gene-count density — number of neighboring expressed genes
density_df <- neighbor_pairs_df %>%
  dplyr::count(gene_id, name = "n_neighbors") %>%
  assign_density_quartiles("n_neighbors")

cat("Gene-count density summary:\n")
print(summary(density_df$n_neighbors))
cat("Quartile counts:\n")
print(table(density_df$density_quartile))

## Metric 2: expression density — mean TPM of neighboring expressed genes
expr_density_df <- neighbor_pairs_df %>%
  dplyr::left_join(tpm_df, by = c("neighbor_gene_id" = "gene_id")) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(mean_neighbor_tpm = mean(tpm, na.rm = TRUE), .groups = "drop") %>%
  assign_density_quartiles("mean_neighbor_tpm")

cat("\nExpression density (mean neighbor TPM) summary:\n")
print(summary(expr_density_df$mean_neighbor_tpm))
cat("Quartile counts:\n")
print(table(expr_density_df$density_quartile))

count_density_groups <- density_df %>%
  dplyr::filter(density_quartile %in% c("Q1", "Q4")) %>%
  dplyr::transmute(gene_id, density_group = as.character(density_quartile))

expr_density_groups <- expr_density_df %>%
  dplyr::filter(density_quartile %in% c("Q1", "Q4")) %>%
  dplyr::transmute(gene_id, density_group = as.character(density_quartile))

##--------------------------------------------------------------## 7b. Genomic E/I ratio (exon length / intron length, from GTF structure)
## Structural property of the gene itself, independent of expression or
## local neighborhood -- exons_by_gene and genes_gr are already built above.

message("Computing genomic E/I ratio...")

exon_lengths <- sum(BiocGenerics::width(exons_by_gene))
names(exon_lengths) <- sub("\\..*", "", names(exon_lengths))   ## strip version to match genes_gr$gene_id

gene_lengths <- setNames(BiocGenerics::width(genes_gr), genes_gr$gene_id)

ei_common  <- intersect(names(exon_lengths), names(gene_lengths))
intron_len <- pmax(gene_lengths[ei_common] - exon_lengths[ei_common], 0)
ei_ratio   <- exon_lengths[ei_common] / (intron_len + 1)

genomic_ei_df <- data.frame(
  gene_id          = ei_common,
  ei_genomic_ratio = as.numeric(ei_ratio),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(ei_genomic_ratio), ei_genomic_ratio > 0)

cat("Genomic E/I ratio computed for", nrow(genomic_ei_df), "genes\n")
print(summary(genomic_ei_df$ei_genomic_ratio))

##--------------------------------------------------------------## 8. Pre-import bigWigs once as RleList

message("Pre-importing bigWigs...")

dmso_rle <- rtracklayer::import.bw(DMSO_BW, as = "RleList")
ko_rle   <- rtracklayer::import.bw(KO_BW,   as = "RleList")

dmso_rle <- dmso_rle[names(dmso_rle) %in% STD_CHRS]
ko_rle   <- ko_rle[names(ko_rle)   %in% STD_CHRS]

message("BigWigs loaded.")

##--------------------------------------------------------------## 9. Per-gene, per-bin ATAC signal (exon & intron, all expressed genes)

message("Computing exon signal (all expressed genes)...")
exon_bins_df <- compute_bin_signal(exons_flat, expressed_ids, N_BINS_EXON,
                                   dmso_rle, ko_rle, label = "Exons")

message("Computing intron signal (all expressed genes)...")
intron_bins_df <- compute_bin_signal(introns_gr, expressed_ids, N_BINS_INTRON,
                                     dmso_rle, ko_rle, label = "Introns")

##--------------------------------------------------------------## 10. Metagene profiles + KO/DMSO delta, both density metrics

exon_count_profile   <- summarise_profile(exon_bins_df,   count_density_groups, N_BINS_EXON)
intron_count_profile <- summarise_profile(intron_bins_df, count_density_groups, N_BINS_INTRON)
exon_count_delta     <- compute_delta(exon_count_profile)
intron_count_delta   <- compute_delta(intron_count_profile)

exon_expr_profile    <- summarise_profile(exon_bins_df,   expr_density_groups, N_BINS_EXON)
intron_expr_profile  <- summarise_profile(intron_bins_df, expr_density_groups, N_BINS_INTRON)
exon_expr_delta      <- compute_delta(exon_expr_profile)
intron_expr_delta    <- compute_delta(intron_expr_profile)

##--------------------------------------------------------------## 11. Plot — Figure A: gene-count density quartiles

p_exon_signal_ct   <- make_signal_plot(exon_count_profile,   "Exon intervals")
p_exon_delta_ct    <- make_delta_plot(exon_count_delta,      "Exon position (%)")
p_intron_signal_ct <- make_signal_plot(intron_count_profile, "Intron intervals")
p_intron_delta_ct  <- make_delta_plot(intron_count_delta,    "Intron position (%)")

fig_count_density <- (p_exon_signal_ct / p_exon_delta_ct / p_intron_signal_ct / p_intron_delta_ct) +
  patchwork::plot_layout(heights = c(2, 1, 2, 1), guides = "collect") &
  theme(legend.position = "top")

fig_count_density <- fig_count_density +
  patchwork::plot_annotation(
    title    = "ATAC-seq accessibility: exon & intron signal with Pol II KO delta",
    subtitle = "Quartiles by local GENE-COUNT density (neighboring expressed genes within 1Mb)",
    theme    = theme(plot.title = element_text(face = "bold", size = 15),
                     plot.subtitle = element_text(color = "grey40", size = 10.5))
  )

print(fig_count_density)

##--------------------------------------------------------------## 12. Plot — Figure B: expression density quartiles

p_exon_signal_ex   <- make_signal_plot(exon_expr_profile,   "Exon intervals")
p_exon_delta_ex    <- make_delta_plot(exon_expr_delta,      "Exon position (%)")
p_intron_signal_ex <- make_signal_plot(intron_expr_profile, "Intron intervals")
p_intron_delta_ex  <- make_delta_plot(intron_expr_delta,    "Intron position (%)")

fig_expr_density <- (p_exon_signal_ex / p_exon_delta_ex / p_intron_signal_ex / p_intron_delta_ex) +
  patchwork::plot_layout(heights = c(2, 1, 2, 1), guides = "collect") &
  theme(legend.position = "top")

fig_expr_density <- fig_expr_density +
  patchwork::plot_annotation(
    title    = "ATAC-seq accessibility: exon & intron signal with Pol II KO delta",
    subtitle = "Quartiles by local EXPRESSION density (mean TPM of neighbors within 1Mb)",
    theme    = theme(plot.title = element_text(face = "bold", size = 15),
                     plot.subtitle = element_text(color = "grey40", size = 10.5))
  )

print(fig_expr_density)

##--------------------------------------------------------------## 13. Per-gene delta vs. continuous local gene density (dose-response)

message("Building per-gene accessibility delta table...")

gene_delta_df <- dplyr::bind_rows(
  exon_bins_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(mean_dmso = mean(dmso, na.rm = TRUE),
                     mean_ko   = mean(ko,   na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(region = "Exon"),
  intron_bins_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(mean_dmso = mean(dmso, na.rm = TRUE),
                     mean_ko   = mean(ko,   na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(region = "Intron")
) %>%
  dplyr::mutate(log2_delta = log2((mean_ko + PSEUDOCOUNT) / (mean_dmso + PSEUDOCOUNT))) %>%
  dplyr::inner_join(density_df %>% dplyr::select(gene_id, n_neighbors), by = "gene_id") %>%
  dplyr::inner_join(expr_density_df %>% dplyr::select(gene_id, mean_neighbor_tpm), by = "gene_id") %>%
  dplyr::inner_join(genomic_ei_df, by = "gene_id") %>%
  dplyr::mutate(
    log10_neighbors      = log10(n_neighbors + 1),
    log10_expr_density   = log10(mean_neighbor_tpm + 1),
    log10_ei             = log10(ei_genomic_ratio)
  )

p_dose_exon   <- plot_dose_response(gene_delta_df %>% dplyr::filter(region == "Exon"),   "Exon")
p_dose_intron <- plot_dose_response(gene_delta_df %>% dplyr::filter(region == "Intron"), "Intron")

fig_dose_response <- (p_dose_exon | p_dose_intron) +
  patchwork::plot_annotation(
    title = "Accessibility change vs. local gene density (continuous, all expressed genes)"
  )

print(fig_dose_response)

##--------------------------------------------------------------## 13b. Signal & delta vs. continuous local density (loess only, no points)
##
## One combined figure: gene-count density, expression density, and genomic
## E/I ratio as three side-by-side columns. Rows are Exon signal / Exon
## delta / Intron signal / Intron delta. Delta panels use the tan-fill
## loess style. No raw points.

message("Building signal + delta vs. continuous density figure...")

exon_df   <- gene_delta_df %>% dplyr::filter(region == "Exon")
intron_df <- gene_delta_df %>% dplyr::filter(region == "Intron")

COUNT_X_COL   <- "log10_neighbors"
COUNT_X_LABEL <- "log10(gene-count density)"
EXPR_X_COL    <- "log10_expr_density"
EXPR_X_LABEL  <- "log10(expression density)"
EI_X_COL      <- "log10_ei"
EI_X_LABEL    <- "log10(E/I ratio)"

p_exon_sig_count   <- plot_signal_vs_density(exon_df,   COUNT_X_COL, "Gene-count density", "Exon")
p_exon_sig_expr    <- plot_signal_vs_density(exon_df,   EXPR_X_COL,  "Expression density", "Exon")
p_exon_sig_ei      <- plot_signal_vs_density(exon_df,   EI_X_COL,    "E/I ratio",          "Exon")
p_exon_delta_count <- plot_delta_vs_density(exon_df,   COUNT_X_COL, COUNT_X_LABEL)
p_exon_delta_expr  <- plot_delta_vs_density(exon_df,   EXPR_X_COL,  EXPR_X_LABEL)
p_exon_delta_ei    <- plot_delta_vs_density(exon_df,   EI_X_COL,    EI_X_LABEL)

p_intron_sig_count   <- plot_signal_vs_density(intron_df, COUNT_X_COL, "Gene-count density", "Intron")
p_intron_sig_expr    <- plot_signal_vs_density(intron_df, EXPR_X_COL,  "Expression density", "Intron")
p_intron_sig_ei      <- plot_signal_vs_density(intron_df, EI_X_COL,    "E/I ratio",          "Intron")
p_intron_delta_count <- plot_delta_vs_density(intron_df, COUNT_X_COL, COUNT_X_LABEL)
p_intron_delta_expr  <- plot_delta_vs_density(intron_df, EXPR_X_COL,  EXPR_X_LABEL)
p_intron_delta_ei    <- plot_delta_vs_density(intron_df, EI_X_COL,    EI_X_LABEL)

## Match y-axis range within each row so the three density metrics are
## directly comparable (each row-group gets its own shared range; signal
## and delta are on different scales so they are not matched to each other)
exon_sig_shared   <- apply_shared_ylim(list(p_exon_sig_count, p_exon_sig_expr, p_exon_sig_ei))
p_exon_sig_count  <- exon_sig_shared[[1]]
p_exon_sig_expr   <- exon_sig_shared[[2]]
p_exon_sig_ei     <- exon_sig_shared[[3]]

exon_delta_shared  <- apply_shared_ylim(list(p_exon_delta_count, p_exon_delta_expr, p_exon_delta_ei))
p_exon_delta_count <- exon_delta_shared[[1]]
p_exon_delta_expr  <- exon_delta_shared[[2]]
p_exon_delta_ei    <- exon_delta_shared[[3]]

intron_sig_shared  <- apply_shared_ylim(list(p_intron_sig_count, p_intron_sig_expr, p_intron_sig_ei))
p_intron_sig_count <- intron_sig_shared[[1]]
p_intron_sig_expr  <- intron_sig_shared[[2]]
p_intron_sig_ei    <- intron_sig_shared[[3]]

intron_delta_shared  <- apply_shared_ylim(list(p_intron_delta_count, p_intron_delta_expr, p_intron_delta_ei))
p_intron_delta_count <- intron_delta_shared[[1]]
p_intron_delta_expr  <- intron_delta_shared[[2]]
p_intron_delta_ei    <- intron_delta_shared[[3]]

fig_trend_combined <- (
  (p_exon_sig_count     | p_exon_sig_expr     | p_exon_sig_ei)     /
    (p_exon_delta_count   | p_exon_delta_expr   | p_exon_delta_ei)   /
    (p_intron_sig_count   | p_intron_sig_expr   | p_intron_sig_ei)   /
    (p_intron_delta_count | p_intron_delta_expr | p_intron_delta_ei)
) +
  patchwork::plot_layout(heights = c(2, 1, 2, 1), guides = "collect") &
  theme(legend.position = "top")

fig_trend_combined <- fig_trend_combined +
  patchwork::plot_annotation(
    title    = "ATAC accessibility vs. local gene density",
    subtitle = "Loess trend + 95% confidence band, all expressed genes (no per-gene points shown)",
    theme    = theme(plot.title = element_text(face = "bold", size = 15),
                     plot.subtitle = element_text(color = "grey40", size = 10.5))
  )

print(fig_trend_combined)

##--------------------------------------------------------------## 14. Correlation with nascent transcription (EU-seq)

message("Loading EU-seq RPKM table and computing nascent transcription log2FC...")

eu_rpkm <- readr::read_tsv(EU_RPKM_PATH, show_col_types = FALSE)

stopifnot(all(c("gene_id", POLR2A_REP_COLS, WT_REP_COLS) %in% names(eu_rpkm)))

## No precomputed log2FC in this file — averaged across replicates here,
## same pseudocount used everywhere else in the script for consistency
eu_df <- eu_rpkm %>%
  dplyr::mutate(
    gene_id     = sub("\\..*", "", gene_id),
    mean_polr2a = rowMeans(dplyr::across(dplyr::all_of(POLR2A_REP_COLS))),
    mean_wt     = rowMeans(dplyr::across(dplyr::all_of(WT_REP_COLS))),
    eu_lfc      = log2((mean_polr2a + PSEUDOCOUNT) / (mean_wt + PSEUDOCOUNT))
  ) %>%
  dplyr::select(gene_id, eu_lfc)

gene_delta_eu_df <- gene_delta_df %>%
  dplyr::inner_join(eu_df, by = "gene_id")

cat("Genes with both ATAC delta and EU-seq lfc:",
    dplyr::n_distinct(gene_delta_eu_df$gene_id), "of",
    dplyr::n_distinct(gene_delta_df$gene_id), "\n")

p_eu_exon   <- plot_eu_correlation(gene_delta_eu_df %>% dplyr::filter(region == "Exon"),   "Exon")
p_eu_intron <- plot_eu_correlation(gene_delta_eu_df %>% dplyr::filter(region == "Intron"), "Intron")

fig_eu_correlation <- (p_eu_exon | p_eu_intron) +
  patchwork::plot_annotation(
    title = "Accessibility change vs. nascent transcription change (EU-seq)"
  )

print(fig_eu_correlation)

message("Exon/intron density + dose-response + EU-seq analysis complete.")


##################################################################
## The section below (intergenic accessibility) is unchanged from
## the original script and does not depend on anything above.
##################################################################

##--------------------------------------------------------------## 15. Intergenic accessibility analysis

message("Computing intergenic regions...")

MIN_INTERGENIC <- 500       ## minimum bp
MAX_INTERGENIC <- 500000    ## maximum bp (500kb cutoff)

## Get expressed gene bodies on standard chromosomes
expr_genes_ucsc <- genes_gr[genes_gr$gene_id %in% expressed_ids]
seqlevelsStyle(expr_genes_ucsc) <- "UCSC"
expr_genes_ucsc <- expr_genes_ucsc[
  as.character(seqnames(expr_genes_ucsc)) %in% STD_CHRS]
expr_genes_ucsc <- GenomeInfoDb::keepSeqlevels(
  expr_genes_ucsc, value = STD_CHRS, pruning.mode = "coarse")

## Sort by chromosome and position
expr_genes_ucsc <- sort(expr_genes_ucsc)

## Build intergenic gaps between adjacent expressed genes
chrom_gaps <- GenomicRanges::gaps(
  GenomicRanges::reduce(expr_genes_ucsc, ignore.strand = TRUE),
  end = seqlengths(expr_genes_ucsc)
)

chrom_gaps <- chrom_gaps[as.character(strand(chrom_gaps)) == "*" |
                           as.character(strand(chrom_gaps)) == "+"]
chrom_gaps <- chrom_gaps[as.character(seqnames(chrom_gaps)) %in% STD_CHRS]

## Size filter
intergenic_gr           <- chrom_gaps[
  width(chrom_gaps) >= MIN_INTERGENIC &
    width(chrom_gaps) <= MAX_INTERGENIC
]
intergenic_gr$region_id <- paste0(seqnames(intergenic_gr), ":",
                                  start(intergenic_gr), "-",
                                  end(intergenic_gr))

cat("Intergenic regions after size filter:",
    length(intergenic_gr), "\n")
cat("Size range: ",
    min(width(intergenic_gr)), "-",
    max(width(intergenic_gr)), " bp\n")

## Verify both flanking genes are expressed
left_check <- GRanges(
  seqnames = seqnames(intergenic_gr),
  ranges   = IRanges(start = start(intergenic_gr) - 1,
                     end   = start(intergenic_gr) - 1)
)
right_check <- GRanges(
  seqnames = seqnames(intergenic_gr),
  ranges   = IRanges(start = end(intergenic_gr) + 1,
                     end   = end(intergenic_gr) + 1)
)

expr_expanded <- expr_genes_ucsc
start(expr_expanded) <- start(expr_expanded) - 1
end(expr_expanded)   <- end(expr_expanded)   + 1

left_hit  <- countOverlaps(left_check,  expr_expanded,
                           ignore.strand = TRUE) > 0
right_hit <- countOverlaps(right_check, expr_expanded,
                           ignore.strand = TRUE) > 0

intergenic_gr <- intergenic_gr[left_hit & right_hit]
intergenic_gr$gene_id <- intergenic_gr$region_id

cat("Intergenic regions flanked by expressed genes:",
    length(intergenic_gr), "\n")

## Assign density quartiles to intergenic regions

message("Assigning density quartiles to intergenic regions...")

intergenic_mid <- GRanges(
  seqnames = seqnames(intergenic_gr),
  ranges   = IRanges(
    start = (start(intergenic_gr) + end(intergenic_gr)) %/% 2,
    end   = (start(intergenic_gr) + end(intergenic_gr)) %/% 2
  ),
  region_id = intergenic_gr$region_id
)

intergenic_window <- GRanges(
  seqnames = seqnames(intergenic_gr),
  ranges   = IRanges(
    start = pmax(1,
                 (start(intergenic_gr) + end(intergenic_gr)) %/% 2 -
                   DENSITY_WINDOW %/% 2),
    end   =
      (start(intergenic_gr) + end(intergenic_gr)) %/% 2 +
      DENSITY_WINDOW %/% 2
  ),
  region_id = intergenic_gr$region_id
)

seqlevelsStyle(midpoints_gr) <- "UCSC"
ig_hits <- findOverlaps(intergenic_window, midpoints_gr,
                        ignore.strand = TRUE)

ig_density_df <- data.frame(
  region_id        = intergenic_window$region_id[queryHits(ig_hits)],
  neighbor_gene_id = midpoints_gr$gene_id[subjectHits(ig_hits)]
) %>%
  dplyr::count(region_id, name = "n_neighbors") %>%
  dplyr::mutate(
    density_quartile = dplyr::case_when(
      n_neighbors < quantile(n_neighbors, 0.25) ~ "Q1 (Low density)",
      n_neighbors < quantile(n_neighbors, 0.50) ~ "Q2",
      n_neighbors < quantile(n_neighbors, 0.75) ~ "Q3",
      TRUE                                       ~ "Q4 (High density)"
    ),
    density_quartile = factor(density_quartile,
                              levels = c("Q1 (Low density)", "Q2",
                                         "Q3", "Q4 (High density)"))
  )

cat("\nIntergenic density quartile counts:\n")
print(table(ig_density_df$density_quartile))

ig_q1_ids <- ig_density_df %>%
  dplyr::filter(density_quartile == "Q1 (Low density)") %>%
  dplyr::pull(region_id)
ig_q4_ids <- ig_density_df %>%
  dplyr::filter(density_quartile == "Q4 (High density)") %>%
  dplyr::pull(region_id)

cat("Intergenic Q1:", length(ig_q1_ids), "\n")
cat("Intergenic Q4:", length(ig_q4_ids), "\n")

## Legacy helpers, kept only for this section so the intergenic plot's
## output stays identical to the original script

build_profile_legacy <- function(intervals_gr, gene_set, n_bins,
                                 density_label, dmso_rle, ko_rle) {
  
  gr_sub <- intervals_gr[intervals_gr$gene_id %in% gene_set]
  gr_sub <- gr_sub[as.character(seqnames(gr_sub)) %in% STD_CHRS]
  gr_sub <- GenomeInfoDb::keepSeqlevels(gr_sub,
                                        value        = STD_CHRS,
                                        pruning.mode = "coarse")
  
  if (length(gr_sub) == 0) {
    warning("No intervals found for: ", density_label)
    return(NULL)
  }
  
  cat("  ", density_label, ": ", length(gr_sub),
      " intervals from ",
      length(unique(gr_sub$gene_id)), " genes\n", sep = "")
  
  bins_gr <- tile_intervals_fast(gr_sub, n_bins)
  if (is.null(bins_gr) || length(bins_gr) == 0) return(NULL)
  
  common_chrs <- intersect(seqlevels(bins_gr), names(dmso_rle))
  bins_gr     <- GenomeInfoDb::keepSeqlevels(bins_gr,
                                             value        = common_chrs,
                                             pruning.mode = "coarse")
  dmso_sub <- dmso_rle[common_chrs]
  ko_sub   <- ko_rle[common_chrs]
  
  bins_gr$dmso <- GenomicRanges::binnedAverage(
    bins = bins_gr, numvar = dmso_sub, varname = "score")$score
  bins_gr$ko   <- GenomicRanges::binnedAverage(
    bins = bins_gr, numvar = ko_sub,   varname = "score")$score
  
  data.frame(
    gene_id = bins_gr$gene_id,
    bin_idx = bins_gr$bin_idx,
    dmso    = bins_gr$dmso,
    ko      = bins_gr$ko,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(bin_idx) %>%
    dplyr::summarise(
      mean_dmso = mean(dmso, na.rm = TRUE),
      mean_ko   = mean(ko,   na.rm = TRUE),
      n_genes   = dplyr::n_distinct(gene_id),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      bin_pct       = (bin_idx - 0.5) / n_bins * 100,
      density_group = density_label
    ) %>%
    tidyr::pivot_longer(
      cols      = c(mean_dmso, mean_ko),
      names_to  = "condition",
      values_to = "mean_signal"
    ) %>%
    dplyr::mutate(
      condition = dplyr::recode(condition,
                                mean_dmso = "DMSO",
                                mean_ko   = "KO")
    )
}

combine_profiles_legacy <- function(q1, q4) {
  dplyr::bind_rows(q1, q4) %>%
    dplyr::mutate(
      density_group = factor(density_group,
                             levels = c("Q1 (Low density)",
                                        "Q4 (High density)"))
    )
}

make_metagene_plot_legacy <- function(df, title, x_label = "Position (%)") {
  
  n_df <- df %>%
    dplyr::filter(condition == "DMSO") %>%
    dplyr::distinct(density_group, n_genes) %>%
    dplyr::group_by(density_group) %>%
    dplyr::slice_max(n_genes, n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      label = paste0("n=", scales::comma(n_genes)),
      x_pos = 3,
      y_pos = Inf
    )
  
  ggplot(df, aes(x        = bin_pct,
                 y        = mean_signal,
                 color    = condition,
                 linetype = condition,
                 group    = condition)) +
    geom_line(linewidth = 0.8) +
    geom_text(data        = n_df,
              aes(x = x_pos, y = y_pos, label = label),
              hjust       = 0, vjust = 1.5,
              size        = 3.5, color = "grey40",
              inherit.aes = FALSE) +
    scale_color_manual(
      values = c("DMSO"  = "#0072B2",
                 "KO" = "#D55E00")
    ) +
    scale_linetype_manual(
      values = c(DMSO = "solid", KO = "dashed"),
      labels = c(DMSO = "DMSO (control)", KO = "KO (Pol II degradation)")
    ) +
    scale_x_continuous(labels = function(x) paste0(x, "%"),
                       breaks = c(0, 25, 50, 75, 100)) +
    facet_wrap(~ density_group, nrow = 1) +
    labs(x        = x_label,
         y        = "Mean enrichment",
         title    = title,
         color    = "Gene density",
         linetype = "Condition") +
    base_theme
}

message("\nBuilding intergenic profiles...")

ig_df <- combine_profiles_legacy(
  build_profile_legacy(intergenic_gr, ig_q1_ids,
                       50, "Q1 (Low density)",  dmso_rle, ko_rle),
  build_profile_legacy(intergenic_gr, ig_q4_ids,
                       50, "Q4 (High density)", dmso_rle, ko_rle)
)

p_intergenic <- make_metagene_plot_legacy(
  ig_df,
  "Intergenic regions",
  "Intergenic position (%)"
)

print(p_intergenic)

message("Intergenic analysis complete.")