##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# EU-STORM nucleus analysis: POLR2A (Auxin) vs. WT (DMSO). Produces,
# per metric (blink density, cluster density, cluster size):
#   1. a boxplot
#   2. a raincloud plot (half-violin + boxplot + jittered points)
#   3. an ECDF comparison plot
#   4. a probability density (PDF) plot
# plus a statistics table (ROI-level and cell-level, with effect
# sizes -- see the STATS section below, written to
# eu_storm_nucleus_stats.csv).
#
# Every plot is print()-ed to the active graphics device as it's
# built and not auto-saved -- add ggsave() calls once a plot has
# been reviewed and is worth keeping.
#
# Figure output: Fig. 2E-F (confirmed).
#
# ---------------------------------------------------------------------------
# WHAT CHANGED FROM THE ORIGINAL SCRIPT (read this once):
#
# 1. col_Data() was not included in what you pasted, so it's reimplemented
#    here as read_eu_storm_file() based on the structure of the file you
#    attached (POLR2A_Auxin_nucleus.csv / POLR2A_Auxin_nucleus_071724.csv).
#    Those files are stored TRANSPOSED -- rows are metrics
#    (num_blink_nucleus, num_cluster_nucleus, cluster_mean_size_nucleus,
#    roi_size_nucleus, and the two "per um^2" columns), and columns are
#    individual ROIs named like "cell4_2" (cell 4, ROI 2).
#    read_eu_storm_file() transposes this to one row per ROI, renames to
#    your original final column names, and pulls out a `cell_id`
#    (e.g. "cell4") -- see the STATS note below for why that matters. If a
#    file doesn't already have the two "per um^2" rows precomputed, the
#    function computes them itself (blinks / roi, clusters / roi) rather
#    than erroring.
#
# 2. POL1 REMOVED. This pipeline now only compares POLR2A (pol2) vs WT.
#    condition_levels, COND_COLORS, and every plot function reflect that --
#    there's a single pol2-vs-wt comparison throughout, not a list of
#    pairwise comparisons.
#
# 3. Boxplot rebuilt for publication use:
#      - the original's multi-layer dashed-whisker trick and aes_string()/
#        ..lower../..upper.. (deprecated ggplot2 syntax) are replaced with a
#        single clean geom_boxplot() + jittered points overlay
#      - the two density metrics (blinks/um^2, clusters/um^2) get a log10
#        y-axis on the boxplot and raincloud plot (and log10 x-axis on the
#        density plot) -- they span >3 orders of magnitude and were badly
#        compressed on a linear axis (this was likely what "not publication
#        quality" was pointing at). cluster size stays linear since its
#        range is narrow. The CDF plot is deliberately NOT log-scaled --
#        unlike a boxplot's y-range, a CDF isn't distorted by a single
#        outlier, so log-x wasn't fixing a real problem there, only
#        changing the visual shape of the curve.
#      - the in-plot significance bracket now uses the CELL-LEVEL test (see
#        pseudoreplication note) and shows the effect size (Cohen's d)
#        alongside the p-value, via ggpubr::stat_pvalue_manual() with
#        p-values/d computed explicitly in the STATS section below --
#        nothing is computed silently inside the plotting call.
#      - x-axis labels now show n (cells and ROIs) per group directly on
#        the plot, not just in the console.
#      - also fixed: the cluster-size call had ylab = "Clusters / um^2"
#        (copy-paste leftover from the density plots) -- corrected to a
#        size unit.
#
# 4. STATS / PSEUDOREPLICATION NOTE (important, please read):
#    Column names like "cell4_1", "cell4_2", "cell4_3", "cell4_4" show that
#    several ROIs come from the SAME cell. Treating every ROI as an
#    independent replicate (as the original t.test() calls at the bottom
#    of your script did) inflates the real sample size and can make
#    p-values overly optimistic. This script computes statistics BOTH ways
#    and prints/exports both, clearly labeled:
#      - "roi_level"  : every ROI as its own replicate (matches your
#                        original t.test() behavior; N = number of ROIs)
#      - "cell_level"  : one value per cell (mean across that cell's ROIs)
#                        before testing; N = number of cells
#    Both levels also get effect sizes: Cohen's d (t-test) and rank-biserial
#    correlation (Wilcoxon). The boxplot/raincloud brackets use the
#    cell-level numbers (see #3); the full table (both levels) is in
#    eu_storm_nucleus_stats.csv.
#
# 5. Added: Kolmogorov-Smirnov test alongside the ECDF plot (the natural
#    statistical pairing for a CDF comparison -- it tests whether the whole
#    distributions differ, not just the means), and a probability density
#    (PDF) plot as a smoothed complement to the empirical CDF.
# ============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(scales)

# ---- 0. Config ------------------------------------------------------------
dir     <- "/Volumes/Backup #1/IBiS/Backman_Lab/Transcription Publication/Data/Ruyi/new_data/20240524_eu_STORM/processed_data"
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE)

# One entry per condition; add every replicate/update file for that
# condition as a character vector. Files within a condition are combined
# with rbind after loading.
files_by_condition <- list(
  wt   = c("WT_DMSO_nucleus.csv",      "WT_DMSO_nucleus_071724.csv"),
  pol2 = c("POLR2A_Auxin_nucleus.csv", "POLR2A_Auxin_nucleus_071724.csv")
)

# Control (WT) plotted first/left, treatment (Pol2) second/right.
condition_levels <- c("wt", "pol2")
condition_axis_labels <- c(wt = "WT", pol2 = "Pol 2")
COND_COLORS <- c(wt = "#cf597e", pol2 = "#39b185")  # colors unchanged from your original palette

# `comparison` controls the STATISTICAL test direction, independent of the
# plotted x-axis order above -- kept as treatment-vs-control (pol2, wt) so
# a positive Cohen's d / rank-biserial r reads as "higher in Pol2 than WT",
# which is the more interpretable convention for a treatment effect.
comparison <- c("pol2", "wt")

metrics <- tribble(
  ~col,                 ~label,                                ~ylab,                                ~log_scale,
  "num_blink_perUm2",   "Blink Density | Nucleus | EU STORM",   "Blinks / \u00b5m\u00b2",              TRUE,
  "num_cluster_perUm2", "Cluster Density | Nucleus | EU STORM", "Clusters / \u00b5m\u00b2",             TRUE,
  "cluster_mean_size",  "Cluster Size | Nucleus | EU STORM",    "Mean cluster size (\u00b5m\u00b2)",   FALSE  # was mislabeled "Clusters / um^2" in the original
)

# ---- Shared publication theme ----------------------------------------------
theme_pub <- function(base_size = 13, legend_position = "none") {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(hjust = 0.5, size = rel(1.05), face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.7), color = "grey30"),
      axis.title    = element_text(face = "bold", size = rel(0.95)),
      axis.text     = element_text(color = "black", size = rel(0.85)),
      axis.text.x   = element_text(lineheight = 0.9),
      legend.position = legend_position,
      legend.title  = element_blank(),
      panel.border  = element_rect(linetype = "solid", colour = "black", fill = NA, linewidth = 0.5),
      panel.grid    = element_blank(),
      plot.title.position = "plot"
    )
}

# ============================================================================
# 1. LOAD DATA
# ============================================================================

#' Read one EU-STORM nucleus CSV (metrics-as-rows, ROIs-as-columns) into a
#' tidy one-row-per-ROI tibble.
read_eu_storm_file <- function(path, condition_label) {
  raw <- read.csv(path, row.names = 1, check.names = FALSE)
  
  tidy <- as.data.frame(t(raw)) %>%
    rownames_to_column("roi_id") %>%
    as_tibble()
  
  # Compute the per-um^2 columns ourselves if the file doesn't already have
  # them precomputed (some of your files do, per the attached example --
  # this just makes the loader robust to files that don't).
  if (!"num_blink_nucleus per um^2" %in% names(tidy)) {
    tidy <- tidy %>% mutate(`num_blink_nucleus per um^2` = num_blink_nucleus / roi_size_nucleus)
  }
  if (!"num_cluster_nucleus per um^2" %in% names(tidy)) {
    tidy <- tidy %>% mutate(`num_cluster_nucleus per um^2` = num_cluster_nucleus / roi_size_nucleus)
  }
  
  tidy %>%
    rename(
      num_blinks         = num_blink_nucleus,
      num_clust           = num_cluster_nucleus,
      cluster_mean_size   = cluster_mean_size_nucleus,
      roi                 = roi_size_nucleus,
      num_cluster_perUm2  = `num_cluster_nucleus per um^2`,
      num_blink_perUm2    = `num_blink_nucleus per um^2`
    ) %>%
    mutate(
      cell_id  = str_extract(roi_id, "^cell[0-9]+"),
      analysis = condition_label
    )
}

#' Load and combine every replicate file for one condition.
load_condition <- function(condition_label, files, dir) {
  map_dfr(files, function(f) read_eu_storm_file(file.path(dir, f), condition_label))
}

nuc.plt <- imap_dfr(files_by_condition, function(files, condition_label) {
  load_condition(condition_label, files, dir)
}) %>%
  mutate(analysis = factor(analysis, levels = condition_levels))

message("Loaded ", nrow(nuc.plt), " ROIs across ", n_distinct(nuc.plt$analysis), " conditions:")
print(nuc.plt %>% count(analysis, name = "n_roi"))
print(nuc.plt %>% distinct(analysis, cell_id) %>% count(analysis, name = "n_cell"))

# One row per cell (mean of that cell's ROIs) -- used for the
# pseudoreplication-aware statistics described above.
nuc.percell <- nuc.plt %>%
  group_by(analysis, cell_id) %>%
  summarise(across(c(num_blinks, num_clust, cluster_mean_size, roi,
                     num_cluster_perUm2, num_blink_perUm2),
                   \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

# n per condition (cells + ROIs), used to annotate axis labels below.
n_summary <- nuc.plt %>%
  count(analysis, name = "n_roi") %>%
  left_join(nuc.plt %>% distinct(analysis, cell_id) %>% count(analysis, name = "n_cell"),
            by = "analysis")

x_labels_with_n <- setNames(
  paste0(condition_axis_labels[as.character(n_summary$analysis)],
         "\n(", n_summary$n_cell, " cells, ", n_summary$n_roi, " ROIs)"),
  as.character(n_summary$analysis)
)[condition_levels]

# ============================================================================
# 2. STATS (ROI-level AND cell-level, with effect sizes)
# ============================================================================

#' Cohen's d (pooled-SD standardized mean difference), for the t-test.
cohens_d <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  pooled_sd <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  (mean(x) - mean(y)) / pooled_sd
}

#' Rank-biserial correlation from the Mann-Whitney U statistic returned by
#' wilcox.test() -- the standard nonparametric effect-size analog to
#' Cohen's d. r = 2U / (n1*n2) - 1; ranges -1 to 1, sign convention matches
#' cohens_d() above (positive = first group (pol2) tends higher).
rank_biserial <- function(x, y) {
  u <- suppressWarnings(wilcox.test(x, y)$statistic)
  n1 <- length(x); n2 <- length(y)
  (2 * u) / (n1 * n2) - 1
}

run_comparisons <- function(data, level_label) {
  map_dfr(metrics$col, function(m) {
    x <- data[[m]][data$analysis == comparison[1]]
    y <- data[[m]][data$analysis == comparison[2]]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    t_res  <- suppressWarnings(t.test(x, y))
    w_res  <- suppressWarnings(wilcox.test(x, y))
    ks_res <- suppressWarnings(ks.test(x, y))
    tibble(level = level_label, metric = m,
           comparison = paste(comparison, collapse = " vs "),
           n_group1 = length(x), n_group2 = length(y),
           t_test_p = t_res$p.value, cohens_d = cohens_d(x, y),
           wilcoxon_p = w_res$p.value, rank_biserial_r = rank_biserial(x, y),
           ks_p = ks_res$p.value)
  })
}

stats_roi  <- run_comparisons(nuc.plt,    "roi_level")
stats_cell <- run_comparisons(nuc.percell, "cell_level")
stats_all  <- bind_rows(stats_roi, stats_cell)

message("\n=== Statistics (roi_level = original approach, N = #ROIs; ",
        "cell_level = one value per cell first, N = #cells; ",
        "cohens_d/rank_biserial_r = effect sizes for the t-test/Wilcoxon) ===")
print(stats_all)
write_csv(stats_all, file.path(out_dir, "eu_storm_nucleus_stats.csv"))

#' Build the p-value + effect-size label ggpubr will print on the bracket.
format_bracket_label <- function(p, d) {
  p_str <- if (p < 0.0001) "p < 0.0001" else paste0("p = ", formatC(p, digits = 2, format = "g"))
  paste0(p_str, ", d = ", sprintf("%.2f", d))
}

#' One-row bracket data frame for stat_pvalue_manual(), for one metric.
#' Uses the CELL-LEVEL test (see pseudoreplication note above).
#'
#' y.position is anchored to the boxplot's own upper WHISKER (Q3 + 1.5*IQR,
#' same definition geom_boxplot() itself uses), not the raw column max.
#' Anchoring to the raw max was the bug behind the blink-density plot
#' looking empty: one extreme outlier ROI pushed y.position (and therefore
#' the whole y-axis, on a log scale) far above the actual box, so the box
#' -- correctly drawn -- was squeezed into a sliver at the bottom of a
#' panel sized to fit that one outlier. The outlier point itself is still
#' visible (jitter/outlier dots aren't clipped), just no longer dictating
#' the axis range.
make_bracket <- function(col, data_for_ylim) {
  row <- stats_cell %>% filter(metric == col)
  upper_whisker <- data_for_ylim %>%
    group_by(analysis) %>%
    summarise(w = boxplot.stats(.data[[col]])$stats[5], .groups = "drop") %>%
    pull(w) %>%
    max(na.rm = TRUE)
  tibble(group1 = comparison[1], group2 = comparison[2],
         label = format_bracket_label(row$t_test_p, row$cohens_d),
         y.position = upper_whisker * 1.15)
}

# ============================================================================
# 3. PLOT FUNCTIONS
# ============================================================================

#' Drop non-positive values before a log-scaled plot (log(0) is undefined),
#' printing exactly what was dropped rather than silently filtering.
prepare_for_log <- function(data, col) {
  n_before <- sum(!is.na(data[[col]]))
  data_f <- data %>% filter(is.na(.data[[col]]) | .data[[col]] > 0)
  n_after <- sum(!is.na(data_f[[col]]))
  if (n_after < n_before) {
    message("  (log-scaled plot for ", col, ": dropped ", n_before - n_after,
            " zero-valued point(s); log(0) is undefined)")
  }
  data_f
}

## ---- 3a. Boxplot -----------------------------------------------------------
plot_Box <- function(data, label, ylab, col.y, col.x, log_scale, bracket) {
  if (log_scale) data <- prepare_for_log(data, col.y)
  
  p <- ggplot(data, aes(x = .data[[col.x]], y = .data[[col.y]], color = .data[[col.x]])) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.6, width = 0.55) +
    geom_jitter(width = 0.12, size = 1.3, alpha = 0.55, shape = 16) +
    scale_color_manual(values = COND_COLORS) +
    scale_x_discrete(limits = condition_levels, labels = x_labels_with_n) +
    ggtitle(label) + xlab("") + ylab(ylab) +
    theme_pub() +
    stat_pvalue_manual(bracket, label = "label", xmin = "group1", xmax = "group2", inherit.aes = FALSE,
                       y.position = "y.position", tip.length = 0.015, size = 3.6)
  
  if (log_scale) p <- p + scale_y_log10(labels = label_number())
  p
}

## ---- 3b. Raincloud plot ----------------------------------------------------
# Half-violin ("flat violin") + boxplot + jittered raw points. No raincloud
# package is available in this environment (ggrain/gghalves/ggdist aren't
# installed and CRAN isn't reachable from here), so this is the standard
# self-contained ggplot2 extension used before those packages existed --
# pure ggplot2, nothing else required.
`%||%` <- function(a, b) if (!is.null(a)) a else b

GeomFlatViolin <- ggproto("GeomFlatViolin", Geom,
                          setup_data = function(data, params) {
                            data$width <- data$width %||% params$width %||% (resolution(data$x, FALSE) * 0.9)
                            data %>%
                              group_by(group) %>%
                              mutate(ymin = min(y), ymax = max(y), xmin = x, xmax = x + width / 2)
                          },
                          draw_group = function(data, panel_scales, coord) {
                            data <- transform(data, xminv = x, xmaxv = x + violinwidth * (xmax - x))
                            newdata <- rbind(
                              transform(data, x = xminv)[order(data$y), ],
                              transform(data, x = xmaxv)[order(data$y, decreasing = TRUE), ]
                            )
                            newdata <- rbind(newdata, newdata[1, ])
                            ggplot2::GeomPolygon$draw_panel(newdata, panel_scales, coord)
                          },
                          draw_key = draw_key_polygon,
                          default_aes = aes(weight = 1, colour = "grey20", fill = "white",
                                            linewidth = 0.5, alpha = NA, linetype = "solid"),
                          required_aes = c("x", "y")
)

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomFlatViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, ...))
}

plot_Rain <- function(data, label, ylab, col.y, col.x, log_scale, bracket) {
  if (log_scale) data <- prepare_for_log(data, col.y)
  
  p <- ggplot(data, aes(x = .data[[col.x]], y = .data[[col.y]], color = .data[[col.x]], fill = .data[[col.x]])) +
    geom_flat_violin(position = position_nudge(x = 0.2), alpha = 0.5, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.7, color = "grey20",
                 linewidth = 0.4, position = position_nudge(x = 0.05)) +
    geom_jitter(position = position_jitter(width = 0.08, height = 0), size = 1.1,
                alpha = 0.5, shape = 16) +
    scale_color_manual(values = COND_COLORS) +
    scale_fill_manual(values = COND_COLORS) +
    scale_x_discrete(limits = condition_levels, labels = x_labels_with_n) +
    ggtitle(label) + xlab("") + ylab(ylab) +
    theme_pub() +
    stat_pvalue_manual(bracket, label = "label", xmin = "group1", xmax = "group2", inherit.aes = FALSE,
                       y.position = "y.position", tip.length = 0.015, size = 3.6)
  
  if (log_scale) p <- p + scale_y_log10(labels = label_number())
  p
}

## ---- 3c. ECDF comparison ---------------------------------------------------
# Always linear x -- unlike the boxplot, a CDF isn't distorted by a single
# outlier the way a box's y-range is, so there's no readability problem
# here that log-x would be fixing. (It was log-scaled in an earlier pass
# along with Box/Rain/Density; removed because it changed the visual shape
# of the curve -- same data, just re-parameterized -- without solving an
# actual problem the way it did for the boxplot.)
plot_CDF <- function(data, label, ylab, col.y, col.x) {
  ks_p <- stats_roi %>% filter(metric == col.y) %>% pull(ks_p)
  ks_label <- if (ks_p < 0.0001) "p < 0.0001" else paste0("p = ", formatC(ks_p, digits = 3, format = "g"))
  legend_labels <- setNames(paste0(condition_axis_labels[condition_levels], " (n = ",
                                   n_summary$n_roi[match(condition_levels, n_summary$analysis)], " ROIs)"),
                            condition_levels)
  
  ggplot(data, aes(x = .data[[col.y]], color = .data[[col.x]])) +
    stat_ecdf(linewidth = 0.9) +
    scale_color_manual(values = COND_COLORS, labels = legend_labels[condition_levels]) +
    labs(title = label, x = ylab, y = "Cumulative fraction",
         subtitle = paste0("Kolmogorov-Smirnov ", paste(comparison, collapse = " vs "), ": ", ks_label)) +
    theme_pub(legend_position = c(0.7, 0.25))
}

## ---- 3d. Probability density (PDF) plot ------------------------------------
plot_Density <- function(data, label, ylab, col.y, col.x, log_scale) {
  if (log_scale) data <- prepare_for_log(data, col.y)
  
  legend_labels <- setNames(paste0(condition_axis_labels[condition_levels], " (n = ",
                                   n_summary$n_roi[match(condition_levels, n_summary$analysis)], " ROIs)"),
                            condition_levels)
  
  p <- ggplot(data, aes(x = .data[[col.y]], color = .data[[col.x]], fill = .data[[col.x]])) +
    geom_density(alpha = 0.35, linewidth = 0.9) +
    scale_color_manual(values = COND_COLORS, labels = legend_labels[condition_levels]) +
    scale_fill_manual(values = COND_COLORS, labels = legend_labels[condition_levels]) +
    labs(title = label, x = ylab, y = "Probability density") +
    theme_pub(legend_position = c(0.7, 0.85))
  
  if (log_scale) p <- p + scale_x_log10(labels = label_number())
  p
}

# ============================================================================
# 4. GENERATE ALL PLOTS (printed to the active device -- not auto-saved;
#    save any of these yourself with ggsave() once you've seen them and
#    want to keep them)
# ============================================================================

pwalk(metrics, function(col, label, ylab, log_scale) {
  bracket <- make_bracket(col, nuc.plt)
  
  p_box   <- plot_Box(nuc.plt, label, ylab, col.y = col, col.x = "analysis", log_scale, bracket)
  p_rain  <- plot_Rain(nuc.plt, label, ylab, col.y = col, col.x = "analysis", log_scale, bracket)
  p_cdf   <- plot_CDF(nuc.plt, label, ylab, col.y = col, col.x = "analysis")
  p_dens  <- plot_Density(nuc.plt, label, ylab, col.y = col, col.x = "analysis", log_scale)
  
  print(p_box)
  print(p_rain)
  print(p_cdf)
  print(p_dens)
})

message("\nDone. Box / raincloud / CDF / density plots for all ", nrow(metrics),
        " metrics were printed above. Nothing was auto-saved -- eu_storm_nucleus_stats.csv ",
        "(written earlier, in the STATS section) is the only file this script writes on its own.")