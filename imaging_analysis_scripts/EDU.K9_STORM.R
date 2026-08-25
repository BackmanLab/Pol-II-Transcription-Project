##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# STORM domain classification scatter: density vs. domain radius.
# Main scatter (density vs. domain radius), with marginal boxplots
# (per-category distribution) on top and right. Quadrant thresholds
# are the wildtype median of each axis; WT points are plain grey
# (not classified into a colored category on the scatter); the
# treatment condition's points are colored by their maturation
# category. No point-size weighting -- every dot is the same size.
# Two versions are produced for each of EDU_STORM/K9_STORM:
# thresholds/classification computed on the raw values, and a
# second version computed on log10-transformed values (see the note
# below on what that does and doesn't change).
#
# Figure output: Fig. 8C-D (probable); Fig. 8F (confirmed);
# SI Fig. 10D-E (confirmed).
#
# ---------------------------------------------------------------------------
# COLUMNS USED:
#   density                -> x-axis, used directly, no transform.
#   seedCluster_size       -> domain_radius_nm = normalization_radius_nm *
#                              sqrt(seedCluster_size). Solved from your
#                              confirmed area formula
#                              (area = seedCluster_size * pi * R^2, R = 80):
#                              r = sqrt(area/pi)
#                                = sqrt(seedCluster_size * pi * R^2 / pi)
#                                = sqrt(seedCluster_size * R^2)
#                                = R * sqrt(seedCluster_size)
#                              (the pi cancels algebraically -- verified
#                              numerically too). domain_area_nm2 is kept as
#                              a back-derived reference-only column (not
#                              used by any plot/threshold/classification).
#   condition              -> used to identify Wildtype vs treatment rows,
#                              and to filter out any other conditions in the
#                              file that aren't part of this comparison.
#                              EDU_STORM and K9_STORM use different values
#                              here (see condition_values in CONFIG) -- EDU
#                              compares WT-DMSO vs POLR2A-Aux (dropping
#                              POLR1A-Aux/POLR1A-DMSO); K9 compares WT vs
#                              Auxin (dropping DMSO).
#
# Quadrant thresholds = median(density) and median(domain_radius_nm), each
# computed from ONLY the Wildtype rows, then applied to classify every
# domain (both conditions) into one of:
#   Large Compact   (density >= cut, radius >= cut)
#   Small Compact   (density >= cut, radius <  cut)
#   Small Uncompact (density <  cut, radius <  cut)
#   Large Uncompact (density <  cut, radius >= cut)
#
# NOTE on the log-transformed second pass: mathematically, computing the
# median on log10(x) and converting back (10^median(log10(x))) gives the
# EXACT SAME threshold value as median(x) directly -- median commutes with
# any monotonic transform, verified numerically (all.equal() confirms it to
# floating-point precision). So the log version's CLASSIFICATION is
# identical to the raw version's, point for point, with one exception: any
# domain with density <= 0 can't be log-transformed and is dropped from the
# log version only (radius should always be > 0 since it's derived from a
# presumably-positive seedCluster_size). The practical value of the log
# version is entirely visual -- log-scaled axes make heavily right-skewed
# STORM count/density data (which this consistently has been, every time
# we've looked at it in this project) much easier to read -- not a
# different classification outcome.
# ============================================================================

library(arrow)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(scales)

# ---- Config ----------------------------------------------------------------
in_dir <- "/Volumes/Backup #1/IBiS/Backman_Lab/Transcription Publication/Pol2 Submission/revised submission/new data/nico/data"
files <- list(
  edu = "EDU_STORM.parquet",
  k9  = "K9_STORM.parquet"
)
normalization_radius_nm <- 80

CONDITION_COL <- "condition"

# EDU_STORM and K9_STORM use completely different condition naming (EDU:
# POLR1A-Aux/POLR1A-DMSO/POLR2A-Aux/WT-DMSO; K9: Auxin/DMSO/WT) -- so these
# are per-file, keyed the same way as `files` above.
condition_values <- list(
  edu = list(wt = "WT-DMSO", treatment = "POLR2A-Aux"),
  k9  = list(wt = "WT",      treatment = "Auxin")
)

CATEGORY_COLORS <- c(
  "Large Compact"   = "#7CAE5B",
  "Small Compact"   = "#8877C4",
  "Small Uncompact" = "#D9A441",
  "Large Uncompact" = "#D33F8C",
  "Wildtype"        = "grey60"
)
category_order <- c("Large Compact", "Small Compact", "Small Uncompact", "Large Uncompact")

# ---- Shared publication theme -----------------------------------------------
theme_pub <- function(base_size = 13, legend_position = "right") {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = rel(1.05), hjust = 0),
      axis.title    = element_text(face = "bold", size = rel(0.95)),
      axis.text     = element_text(color = "black", size = rel(0.85)),
      legend.position = legend_position,
      legend.title  = element_blank(),
      panel.grid    = element_blank()
    )
}

# ============================================================================
# 0. INSPECT RAW DATA (quick sanity check before the real pipeline runs below)
# ============================================================================

iwalk(files, function(f, label) {
  raw <- read_parquet(file.path(in_dir, f))
  cat("\n==== ", label, ": ", f, " ====\n", sep = "")
  cat("dim:", paste(dim(raw), collapse = " x "), "\n")
  print(glimpse(raw))
  cat("\nSummary of numeric columns:\n")
  print(summary(raw))
  if (CONDITION_COL %in% names(raw)) {
    cat("\n'", CONDITION_COL, "' counts:\n", sep = "")
    print(table(raw[[CONDITION_COL]]))
  }
})

# ============================================================================
# 1. LOAD + DERIVE COLUMNS
# ============================================================================

#' Load one STORM domain parquet file, derive domain_radius_nm, identify
#' Wildtype rows, and classify every domain into a maturation category using
#' quadrant thresholds set from the WILDTYPE median of each axis.
#'
#' @param log_transform if TRUE, thresholds are computed on log10(density)/
#'   log10(domain_radius_nm) instead of the raw values (see NOTE at top of
#'   file for why this produces the same classification either way, modulo
#'   dropping non-positive density rows that can't be logged).
#' Returns a harmonized set of columns regardless of which file was loaded
#' (EDU_STORM and K9_STORM otherwise carry different dataset-specific raw
#' columns alongside these, but share the same 21-column schema overall).
load_storm_domains <- function(path, wt_value, treatment_value, log_transform = FALSE) {
  df <- read_parquet(path)
  message("\nLoaded '", basename(path), "': ", nrow(df), " domains.",
          if (log_transform) " [log-transformed pass]" else " [raw pass]")
  message("  Columns present: ", paste(names(df), collapse = ", "))
  
  # -- x-axis: density, used directly --------------------------------------
  if (!"density" %in% names(df)) {
    stop("Expected column 'density' not found in ", basename(path),
         ". Available columns: ", paste(names(df), collapse = ", "))
  }
  
  # -- y-axis: domain radius, solved from your confirmed area formula --------
  if (!"seedCluster_size" %in% names(df)) {
    stop("Expected column 'seedCluster_size' not found in ", basename(path),
         ". Available columns: ", paste(names(df), collapse = ", "))
  }
  df$domain_radius_nm <- normalization_radius_nm * sqrt(df$seedCluster_size)
  # domain_area_nm2: NOT used by anything below -- kept only as a
  # back-derived reference column, the reverse of the situation last version
  # (radius was the reference column then; area is now).
  df$domain_area_nm2 <- df$seedCluster_size * pi * normalization_radius_nm^2
  
  # -- keep only the two conditions this plot compares -----------------------
  if (!CONDITION_COL %in% names(df)) {
    stop("CONDITION_COL = '", CONDITION_COL, "' not found. Available columns: ",
         paste(names(df), collapse = ", "))
  }
  cond_values <- unique(as.character(df[[CONDITION_COL]]))
  message("  '", CONDITION_COL, "' values found: ", paste(cond_values, collapse = ", "))
  
  if (!(wt_value %in% cond_values)) {
    stop("wt_value = '", wt_value, "' not present in this file. Values found: ",
         paste(cond_values, collapse = ", "))
  }
  if (!(treatment_value %in% cond_values)) {
    stop("treatment_value = '", treatment_value, "' not present in this file. Values found: ",
         paste(cond_values, collapse = ", "))
  }
  
  n_before <- nrow(df)
  df <- df %>% dplyr::filter(.data[[CONDITION_COL]] %in% c(wt_value, treatment_value))
  message("  Keeping only '", wt_value, "' and '", treatment_value, "': ", nrow(df),
          " of ", n_before, " domains (dropped ", n_before - nrow(df),
          " from other conditions: ", paste(setdiff(cond_values, c(wt_value, treatment_value)), collapse = ", "), ").")
  
  df$is_wt <- as.character(df[[CONDITION_COL]]) == wt_value
  message("  Wildtype ('", wt_value, "'): ", sum(df$is_wt), " domains. Treatment ('",
          treatment_value, "'): ", sum(!df$is_wt), " domains.")
  
  # Defensive: force plain numeric R vectors. arrow::read_parquet() should
  # already return standard vectors, not lazy/chunked Arrow types, but given
  # the arrow::select() S4-generic collision hit earlier in this task, I'm
  # not fully trusting that every column behaves like a normal R vector for
  # every operation without seeing your real data myself.
  df$density <- as.numeric(df$density)
  df$domain_radius_nm <- as.numeric(df$domain_radius_nm)
  
  if (log_transform) {
    n_pre <- nrow(df)
    df <- df %>% dplyr::filter(density > 0, domain_radius_nm > 0)
    if (nrow(df) < n_pre) {
      message("  [log pass] Dropped ", n_pre - nrow(df),
              " domain(s) with non-positive density or radius (log undefined).")
    }
    x_metric <- log10(df$density)
    y_metric <- log10(df$domain_radius_nm)
  } else {
    x_metric <- df$density
    y_metric <- df$domain_radius_nm
  }
  
  # -- quadrant thresholds: WILDTYPE median of each axis ----------------------
  # Diagnostic FIRST: the actual WT-only distribution feeding the median, so
  # a threshold that looks surprising in the plot can be checked directly
  # against real numbers instead of taken on faith.
  message("  Wildtype-only density summary: ", paste(names(summary(df$density[df$is_wt])),
                                                     round(summary(df$density[df$is_wt]), 4), sep = "=", collapse = ", "))
  message("  Wildtype-only domain_radius_nm summary: ", paste(names(summary(df$domain_radius_nm[df$is_wt])),
                                                              round(summary(df$domain_radius_nm[df$is_wt]), 1), sep = "=", collapse = ", "))
  
  x_cut <- median(x_metric[df$is_wt], na.rm = TRUE)
  y_cut <- median(y_metric[df$is_wt], na.rm = TRUE)
  message("  Quadrant thresholds (Wildtype median", if (log_transform) ", log10-space" else "",
          "): x @ ", round(x_cut, 4), ", y @ ", round(y_cut, 4))
  
  df <- df %>%
    dplyr::mutate(category = case_when(
      x_metric >= x_cut & y_metric >= y_cut ~ "Large Compact",
      x_metric >= x_cut & y_metric <  y_cut ~ "Small Compact",
      x_metric <  x_cut & y_metric <  y_cut ~ "Small Uncompact",
      x_metric <  x_cut & y_metric >= y_cut ~ "Large Uncompact",
      TRUE ~ NA_character_
    ))
  df$category <- factor(df$category, levels = category_order)
  
  # Display grouping for the SCATTER only: Wildtype points are always plain
  # grey regardless of their own quadrant; treatment points are colored by
  # category. (category itself is still computed for every row -- used by
  # the marginal boxplots and available if you want it later.)
  df$display_group <- ifelse(df$is_wt, "Wildtype", as.character(df$category))
  df$display_group <- factor(df$display_group, levels = c("Wildtype", category_order))
  
  message("  Rows with unclassified (NA) category, will be dropped: ", sum(is.na(df$category)), " of ", nrow(df))
  
  # -- direct crosstab, WT vs treatment, per category --------------------------
  message("  Category crosstab (Wildtype vs treatment), before dropping NA:")
  print(table(category = df$category, is_wt = df$is_wt, useNA = "ifany"))
  
  df <- df %>%
    dplyr::filter(!is.na(category)) %>%
    dplyr::select(any_of("cell_id"), condition = all_of(CONDITION_COL), is_wt,
                  density, domain_radius_nm, domain_area_nm2, category, display_group)
  
  # Threshold attributes are stored on the ORIGINAL scale (not log10) even
  # for the log pass, since the plot always displays raw density/radius --
  # log_scale = TRUE on the plot just changes the axis, not the data, so the
  # dashed lines need to sit at the real (unlogged) threshold values.
  attr(df, "x_cut") <- if (log_transform) 10^x_cut else x_cut
  attr(df, "y_cut") <- if (log_transform) 10^y_cut else y_cut
  attr(df, "log_transform") <- log_transform
  df
}

# ============================================================================
# 2. PLOT: scatter + marginal boxplots
# ============================================================================

#' Build the full composite figure (main scatter + top/right marginal
#' boxplots) for one loaded/derived domain table. Reads the classification
#' cutoffs from attr(df, "x_cut")/attr(df, "y_cut") -- set by
#' load_storm_domains() from the Wildtype median -- so the dashed lines
#' always match whatever was actually used to classify. Set log_scale = TRUE
#' to display density/radius on log10 axes (visual only -- see NOTE at top
#' of file; does not change which points are in which category).
plot_domain_classification <- function(df, title, log_scale = FALSE) {
  x_cut <- attr(df, "x_cut")
  y_cut <- attr(df, "y_cut")
  
  p_main <- ggplot(df, aes(x = density, y = domain_radius_nm)) +
    geom_point(data = df %>% dplyr::filter(is_wt), aes(color = display_group),
               size = 1.4, alpha = 0.7) +
    geom_point(data = df %>% dplyr::filter(!is_wt), aes(color = display_group),
               size = 1.4, alpha = 0.7) +
    geom_vline(xintercept = x_cut, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
    geom_hline(yintercept = y_cut, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
    scale_color_manual(values = CATEGORY_COLORS) +
    labs(x = "Density", y = "Domain Radius (nm)", color = NULL) +
    theme_pub() +
    theme(legend.position = c(0.16, 0.85), legend.background = element_blank())
  
  p_top <- ggplot(df, aes(x = density, y = category, color = display_group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.5, show.legend = FALSE,
                 position = position_dodge(width = 0.8)) +
    geom_vline(xintercept = x_cut, linetype = "dashed", color = "firebrick", linewidth = 0.4) +
    scale_color_manual(values = CATEGORY_COLORS) +
    labs(x = NULL, y = NULL, title = title) +
    theme_pub(legend_position = "none") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  p_right <- ggplot(df, aes(x = category, y = domain_radius_nm, color = display_group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.5, show.legend = FALSE,
                 position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = y_cut, linetype = "dashed", color = "firebrick", linewidth = 0.4) +
    scale_color_manual(values = CATEGORY_COLORS) +
    labs(x = NULL, y = NULL) +
    theme_pub(legend_position = "none") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  if (log_scale) {
    p_main  <- p_main + scale_x_log10(labels = label_number()) + scale_y_log10(labels = label_number())
    p_top   <- p_top + scale_x_log10(labels = label_number(), limits = range(df$density[df$density > 0], na.rm = TRUE))
    p_right <- p_right + scale_y_log10(labels = label_number(), limits = range(df$domain_radius_nm[df$domain_radius_nm > 0], na.rm = TRUE))
  } else {
    p_top   <- p_top + scale_x_continuous(limits = range(df$density, na.rm = TRUE))
    p_right <- p_right + scale_y_continuous(limits = range(df$domain_radius_nm, na.rm = TRUE))
  }
  
  p_spacer <- plot_spacer()
  
  (p_top + p_spacer + p_main + p_right) +
    plot_layout(widths = c(4, 1), heights = c(1, 4), design = "AB\nCD")
}

# ============================================================================
# 3. RUN (printed to the active device, not auto-saved)
# ============================================================================

domain_data <- imap(files, function(f, label) {
  load_storm_domains(file.path(in_dir, f),
                     wt_value = condition_values[[label]]$wt,
                     treatment_value = condition_values[[label]]$treatment,
                     log_transform = FALSE)
})

print(plot_domain_classification(domain_data$edu, "EdU-STORM domains"))
print(plot_domain_classification(domain_data$k9,  "H3K9-STORM domains"))

# -- second pass: thresholds/classification computed on log10-transformed
# data (see NOTE at top of file -- classification is mathematically
# identical to the raw pass except for dropped non-positive rows; this is
# primarily a more legible visualization for heavily right-skewed data).
domain_data_log <- imap(files, function(f, label) {
  load_storm_domains(file.path(in_dir, f),
                     wt_value = condition_values[[label]]$wt,
                     treatment_value = condition_values[[label]]$treatment,
                     log_transform = TRUE)
})

print(plot_domain_classification(domain_data_log$edu, "EdU-STORM domains (log10 thresholds + axes)", log_scale = TRUE))
print(plot_domain_classification(domain_data_log$k9,  "H3K9-STORM domains (log10 thresholds + axes)", log_scale = TRUE))

# ============================================================================
# 4. QUADRANT FRACTIONS: treatment/Wildtype ratio per category, with CI
# ============================================================================

#' Bar chart: for each category, the ratio of treatment-condition domain
#' count to Wildtype domain count in that quadrant, with a 95% CI error bar,
#' colored by quadrant (matching the scatter's category colors) and the
#' ratio value labeled above each bar.
#'
#' CI construction: a ratio of two independent counts doesn't have a
#' textbook CI the way a single proportion does, and can legitimately exceed
#' 1 (more treatment domains than WT domains in a quadrant) -- which a naive
#' prop.test(treatment_count, wt_count) would choke on, since that call
#' expects the first count to be a SUBSET of the second (successes out of
#' trials), not two independent totals. Instead: reframe as p = treatment /
#' (treatment + wt), which IS a valid proportion in [0,1], get its CI via
#' prop.test(), then back-transform p and its CI bounds through
#' ratio = p / (1-p) to return to the treatment/WT ratio scale. This stays
#' well-defined even when the ratio is above 1 (verified numerically).
plot_quadrant_fractions <- function(df, title, treatment_label) {
  counts <- df %>%
    dplyr::count(category, is_wt, .drop = FALSE) %>%
    tidyr::pivot_wider(names_from = is_wt, values_from = n, values_fill = 0,
                       names_prefix = "n_")
  # n_wt / n_wt columns from pivot_wider are named n_TRUE/n_FALSE (is_wt values)
  names(counts)[names(counts) == "n_TRUE"]  <- "n_wt"
  names(counts)[names(counts) == "n_FALSE"] <- "n_treatment"
  
  if (any(counts$n_wt == 0)) {
    message("  WARNING: at least one category has zero Wildtype domains -- its ratio/CI is undefined (Inf).")
  }
  
  counts <- counts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      total = n_treatment + n_wt,
      p_hat = n_treatment / total,
      p_ci  = list(if (total > 0) prop.test(n_treatment, total, correct = FALSE)$conf.int else c(NA_real_, NA_real_))
    ) %>%
    dplyr::mutate(
      p_lo = p_ci[1],
      p_hi = p_ci[2]
    ) %>%
    dplyr::mutate(
      ratio       = n_treatment / n_wt,
      ratio_lower = p_lo / (1 - p_lo),
      ratio_upper = p_hi / (1 - p_hi)
    ) %>%
    dplyr::ungroup()
  
  message("  Quadrant fractions (", treatment_label, " / Wildtype), with 95% CI:")
  print(counts %>% dplyr::select(category, n_wt, n_treatment, ratio, ratio_lower, ratio_upper))
  
  ggplot(counts, aes(x = category, y = ratio, fill = category)) +
    geom_col(color = "black", linewidth = 1, show.legend = FALSE) +
    geom_errorbar(aes(ymin = ratio_lower, ymax = ratio_upper), width = 0.2, linewidth = 1) +
    geom_label(aes(label = round(ratio, 3)), fill = "white", vjust = -0.3, label.size = 1) +
    scale_fill_manual(values = CATEGORY_COLORS) +
    labs(title = title, x = NULL, y = paste0("Fraction (", treatment_label, " / Wildtype)")) +
    theme_pub(legend_position = "none")
}

print(plot_quadrant_fractions(domain_data$edu, "EdU-STORM: quadrant fractions", condition_values$edu$treatment))
print(plot_quadrant_fractions(domain_data$k9,  "H3K9-STORM: quadrant fractions", condition_values$k9$treatment))

# ============================================================================
# EXTRA STORM domain plots: ridgelines, CDFs, radius-bin barplot
#
# Depends on objects created by domain_classification_plot.R -- run that
# script first (in the same R session) before this one:
#   domain_data          (list: edu, k9 -- raw-pass harmonized data)
#   condition_values      (list: edu, k9 -- $wt / $treatment condition labels)
#   category_order, CATEGORY_COLORS, theme_pub()
#
# Adds, for both density and domain_radius_nm, for both EDU_STORM and
# K9_STORM: a ridgeline density plot (Wildtype vs treatment) and a CDF plot
# with a Kolmogorov-Smirnov test. Also adds a radius-bin barplot (domain
# counts split into Small/Average/Large radius, Wildtype vs treatment).
#
# Two-color scheme used throughout this file (Wildtype vs one treatment
# condition -- not the same as the 4-category quadrant colors used
# elsewhere): grey60 for Wildtype, "#8877C4" (the same purple used for the
# quadrant fraction bars' treatment fill) for the treatment condition.
# ============================================================================

library(ggridges)
library(dplyr)
library(ggplot2)
library(scales)

TWO_COND_COLORS_BASE <- c(wildtype = "grey60", treatment = "#8877C4")

#' Add a `condition_label` factor column (Wildtype vs the real treatment
#' name) to a harmonized domain_data data frame, for the plots below.
add_condition_label <- function(df, treatment_label, wt_first = TRUE) {
  lvls <- if (wt_first) c("Wildtype", treatment_label) else c(treatment_label, "Wildtype")
  df %>%
    dplyr::mutate(condition_label = factor(ifelse(is_wt, "Wildtype", treatment_label), levels = lvls))
}

# ============================================================================
# Ridgeline density plots
# ============================================================================

#' One ridgeline plot: `metric_col` (density or domain_radius_nm),
#' Wildtype vs treatment, one ridge per condition.
plot_ridgeline <- function(df, metric_col, metric_label, title, treatment_label) {
  df <- add_condition_label(df, treatment_label, wt_first = FALSE)  # treatment ridge on top, WT below -- matches typical ridgeline reading order
  colors <- setNames(c(TWO_COND_COLORS_BASE[["wildtype"]], TWO_COND_COLORS_BASE[["treatment"]]),
                     c("Wildtype", treatment_label))
  
  ggplot(df, aes(x = .data[[metric_col]], y = condition_label, fill = condition_label)) +
    geom_density_ridges(alpha = 0.75, color = "black", linewidth = 0.4, scale = 1.15) +
    scale_fill_manual(values = colors) +
    labs(title = title, x = metric_label, y = NULL, fill = NULL) +
    theme_pub(legend_position = "none")
}

# ============================================================================
# CDF plots (with Kolmogorov-Smirnov test)
# ============================================================================

#' One ECDF comparison plot for `metric_col`, Wildtype vs treatment, with a
#' KS test p-value in the subtitle (same convention as the EU-STORM
#' nucleus-analysis script's CDF plots).
plot_metric_cdf <- function(df, metric_col, metric_label, title, treatment_label) {
  df <- add_condition_label(df, treatment_label, wt_first = TRUE)
  colors <- setNames(c(TWO_COND_COLORS_BASE[["wildtype"]], TWO_COND_COLORS_BASE[["treatment"]]),
                     c("Wildtype", treatment_label))
  
  x <- df[[metric_col]][df$is_wt]
  y <- df[[metric_col]][!df$is_wt]
  ks_p <- suppressWarnings(ks.test(x, y)$p.value)
  ks_label <- if (ks_p < 0.0001) "p < 0.0001" else paste0("p = ", formatC(ks_p, digits = 3, format = "g"))
  
  ggplot(df, aes(x = .data[[metric_col]], color = condition_label)) +
    stat_ecdf(linewidth = 0.9) +
    scale_color_manual(values = colors) +
    labs(title = title, x = metric_label, y = "Cumulative fraction", color = NULL,
         subtitle = paste0("Wildtype vs ", treatment_label, ", Kolmogorov-Smirnov: ", ks_label)) +
    theme_pub(legend_position = c(0.78, 0.25))
}

# ============================================================================
# Run: ridgeline + CDF for density and domain_radius_nm, EDU and K9
# ============================================================================

metric_specs <- list(
  list(col = "density",          label = "Density"),
  list(col = "domain_radius_nm", label = "Domain Radius (nm)")
)
dataset_specs <- list(
  list(key = "edu", label = "EdU-STORM"),
  list(key = "k9",  label = "H3K9-STORM")
)

for (ds in dataset_specs) {
  df <- domain_data[[ds$key]]
  treatment_label <- condition_values[[ds$key]]$treatment
  for (m in metric_specs) {
    print(plot_ridgeline(df, m$col, m$label,
                         paste0(ds$label, ": ", m$label, " (Wildtype vs ", treatment_label, ")"),
                         treatment_label))
    print(plot_metric_cdf(df, m$col, m$label,
                          paste0(ds$label, ": ", m$label, " CDF"),
                          treatment_label))
  }
}

# ============================================================================
# Radius-bin barplot: domains split into Small / Average / Large radius
# ============================================================================
# Fixed nm thresholds (NOT Wildtype-median-derived like the rest of this
# pipeline) taken directly from your own ChromSTEM reference script:
#   Rmax <= 80         -> Small
#   80 < Rmax <= 120    -> Average
#   Rmax >= 120         -> Large
# Worth knowing this is a different, absolute-threshold classification
# scheme than the quadrant categories used elsewhere in this project.

RADIUS_BIN_BREAKS <- c(-Inf, 80, 120, Inf)
RADIUS_BIN_LABELS <- c("Small", "Average", "Large")

#' Grouped bar chart: number of domains per radius bin, Wildtype next to
#' treatment, count labeled above each bar.
plot_radius_bins <- function(df, title, treatment_label) {
  df <- df %>%
    dplyr::mutate(radius_bin = cut(domain_radius_nm, breaks = RADIUS_BIN_BREAKS,
                                   labels = RADIUS_BIN_LABELS, right = TRUE)) %>%
    add_condition_label(treatment_label, wt_first = TRUE)
  
  counts <- df %>% dplyr::count(radius_bin, condition_label, .drop = FALSE)
  colors <- setNames(c(TWO_COND_COLORS_BASE[["wildtype"]], TWO_COND_COLORS_BASE[["treatment"]]),
                     c("Wildtype", treatment_label))
  
  ggplot(counts, aes(x = factor(radius_bin, levels = RADIUS_BIN_LABELS), y = n, fill = condition_label)) +
    geom_col(position = position_dodge(width = 0.75), color = "black", linewidth = 1, width = 0.65) +
    geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.4,
              size = 4.5, fontface = "bold") +
    scale_fill_manual(values = colors) +
    labs(title = title, x = NULL, y = "Number of domains", fill = NULL) +
    theme_pub()
}

print(plot_radius_bins(domain_data$edu, "EdU-STORM: domains by radius bin", condition_values$edu$treatment))
print(plot_radius_bins(domain_data$k9,  "H3K9-STORM: domains by radius bin", condition_values$k9$treatment))

# ============================================================================
# ============================================================================
# WT vs Treatment scatter (Radius vs Density), colored by WT-derived radius
# tercile, plus a matching bar chart -- recreates the reference figure style
# (two side-by-side scatter panels + WT-calibrated break lines), minus the
# regression line / R^2 / slope annotation (dropped per request).
#
# Depends on domain_data, condition_values, theme_pub() from
# domain_classification_plot.R -- run that script first.
#
# NOTE ON AXES: the reference image has Radius on x and Density on y --
# this is the OPPOSITE orientation from the main quadrant script's scatter
# (which has density on x, radius on y). Matched here to the reference
# image exactly; flag it if that was meant to be consistent with the other
# plot instead, and I'll swap it.
#
# NOTE ON THE BREAK LINES: fixed absolute thresholds (<=80nm / 80-120nm /
# >=120nm), the SAME scheme used in domain_extra_plots.R's radius-bin
# barplot (from your ChromSTEM reference script). This was originally
# WT-percentile-derived (33rd/66th percentile of WT's own radius
# distribution) -- switched to fixed thresholds for consistency with the
# other generated plots. Tradeoff: WT no longer shows a perfectly flat
# count across bins by mathematical construction the way the percentile
# version did (that flatness wasn't a bug -- splitting a population by its
# own terciles guarantees near-equal bins for that population -- it just
# wasn't the comparison wanted here).
# ============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(patchwork)
library(scales)

RADIUS_BIN_COLORS <- c(
  "Small radius"   = "#7CAE5B",
  "Average radius" = "#8877C4",
  "Large radius"   = "#D9A441"
)
radius_bin_order <- c("Small radius", "Average radius", "Large radius")

#' Classify every domain (both conditions) into a FIXED radius bin. Uses the
#' same absolute thresholds as domain_extra_plots.R's radius-bin barplot
#' (from your ChromSTEM reference script: <=80nm / 80-120nm / >=120nm) --
#' NOT Wildtype-percentile-derived like the previous version of this
#' function. Switched per request, for consistency with the other generated
#' plots; the tradeoff is WT no longer shows a perfectly flat count across
#' bins by construction (it showed 2094/2094/2094 before purely because the
#' breaks were WT's own terciles -- a mathematical guarantee of that
#' approach, not a bug, but not the comparison wanted here either).
RADIUS_BIN_BREAKS_FIXED <- c(-Inf, 80, 120, Inf)

add_wt_radius_bins <- function(df) {
  message("  Radius bin thresholds (fixed): <=", RADIUS_BIN_BREAKS_FIXED[2], "nm, ",
          RADIUS_BIN_BREAKS_FIXED[2], "-", RADIUS_BIN_BREAKS_FIXED[3], "nm, >",
          RADIUS_BIN_BREAKS_FIXED[3], "nm")
  df <- df %>%
    dplyr::mutate(radius_bin = case_when(
      domain_radius_nm <= RADIUS_BIN_BREAKS_FIXED[2] ~ "Small radius",
      domain_radius_nm >  RADIUS_BIN_BREAKS_FIXED[2] & domain_radius_nm <= RADIUS_BIN_BREAKS_FIXED[3] ~ "Average radius",
      domain_radius_nm >  RADIUS_BIN_BREAKS_FIXED[3] ~ "Large radius",
      TRUE ~ NA_character_
    ))
  df$radius_bin <- factor(df$radius_bin, levels = radius_bin_order)
  attr(df, "break_lo") <- RADIUS_BIN_BREAKS_FIXED[2]
  attr(df, "break_hi") <- RADIUS_BIN_BREAKS_FIXED[3]
  df
}

#' Build the 3-panel composite: WT scatter | treatment scatter | bar chart.
#' WT scatter is plain grey; treatment scatter is colored by radius_bin;
#' both share the same two WT-derived break lines (grey dashed). No
#' regression line, R^2, or slope annotation.
plot_wt_treatment_radius <- function(raw_df, title_prefix, treatment_label) {
  df <- add_wt_radius_bins(raw_df)
  break_lo <- attr(df, "break_lo")
  break_hi <- attr(df, "break_hi")
  
  # SHARED axis limits, computed from the FULL combined data (both
  # conditions together) -- not each panel's own range. This is the actual
  # fix: scatter_base() previously called scale_x_log10()/scale_y_log10()
  # with no explicit limits, so each panel independently auto-scaled to
  # its OWN data's min/max. WT and treatment don't have identical data
  # ranges, so the two panels ended up on different axis scales, making
  # them impossible to compare directly by eye -- a point drawn "high" in
  # one panel wasn't necessarily at the same actual value as an
  # equally-high-looking point in the other. Small pad (via range() on
  # log-transformed data, mapped back) keeps points from sitting exactly
  # on the panel edge.
  x_range_shared <- range(df$domain_radius_nm, na.rm = TRUE)
  y_range_shared <- range(df$density, na.rm = TRUE)
  
  scatter_base <- function(d, title) {
    ggplot(d, aes(x = domain_radius_nm, y = density)) +
      geom_vline(xintercept = c(break_lo, break_hi), linetype = "dashed",
                 color = "grey50", linewidth = 0.5) +
      scale_x_log10(labels = label_number(), limits = x_range_shared) +
      scale_y_log10(labels = label_number(), limits = y_range_shared) +
      labs(title = title, x = "Radius (nm)", y = expression("Density (1/nm"^2*")")) +
      theme_pub(legend_position = "none")
  }
  
  p_wt <- scatter_base(df %>% dplyr::filter(is_wt), "WT") +
    geom_point(color = "grey60", alpha = 0.7, size = 1.4)
  
  # Label each bin's nm range at the bottom of the WT panel, since WT points
  # aren't colored by bin (unlike the treatment panel) -- without this,
  # there's nothing in the WT panel itself indicating what the two dashed
  # lines divide it into. x-positions are the geometric mean of each bin's
  # boundaries (appropriate for a log10 x-axis); y-position sits just above
  # the WT panel's own minimum density value.
  wt_density <- df$density[df$is_wt]
  wt_radius  <- df$domain_radius_nm[df$is_wt]
  label_y    <- min(wt_density, na.rm = TRUE) * 1.05
  label_x <- c(
    sqrt(min(wt_radius, na.rm = TRUE) * break_lo),  # center of <=80nm bin
    sqrt(break_lo * break_hi),                       # center of 80-120nm bin
    sqrt(break_hi * max(wt_radius, na.rm = TRUE))    # center of >=120nm bin
  )
  bin_labels <- c(paste0("<=", break_lo, "nm"),
                  paste0(break_lo, "-", break_hi, "nm"),
                  paste0(">=", break_hi, "nm"))
  
  p_wt <- p_wt +
    annotate("text", x = label_x, y = label_y, label = bin_labels,
             color = "black", size = 3.4, vjust = 1)
  
  p_treat <- scatter_base(df %>% dplyr::filter(!is_wt), treatment_label) +
    geom_point(aes(color = radius_bin), alpha = 0.7, size = 1.4) +
    scale_color_manual(values = RADIUS_BIN_COLORS)
  
  # Bar chart: WT (grey) next to treatment (colored per bin), same break
  # classification as the scatters above. bar_key is an explicit factor
  # with "Wildtype" listed first so it dodges to the LEFT of the treatment
  # bar in every group (position_dodge default-sorts alphabetically
  # otherwise, which put Wildtype on the right since e.g. "Small radius" <
  # "Wildtype" alphabetically -- control-on-the-left is a standing
  # preference here).
  bar_counts <- df %>%
    dplyr::mutate(bar_key = factor(ifelse(is_wt, "Wildtype", as.character(radius_bin)),
                                   levels = c("Wildtype", radius_bin_order))) %>%
    dplyr::count(radius_bin, bar_key, .drop = FALSE) %>%
    # .drop = FALSE (needed so a genuinely-empty WT-or-treatment bin still
    # shows as a labeled 0 bar) also generates nonsensical cross-bin rows
    # from bar_key's other radius-bin levels (e.g. radius_bin = "Small
    # radius" with bar_key = "Average radius") -- those don't correspond to
    # any real data and must be dropped, or they'd add two invisible extra
    # dodge slots per x position and throw off the WT/treatment bar widths.
    dplyr::filter(bar_key == "Wildtype" | as.character(bar_key) == as.character(radius_bin))
  bar_colors <- c(Wildtype = "grey60", RADIUS_BIN_COLORS)
  
  # Both the CI (below) and the significance test (further down) move to
  # CELL-LEVEL statistics: domains are clustered within cells, not
  # independent draws, so pooling all domains together (as the earlier
  # domain-level prop.test/poisson.test versions did) understates the true
  # uncertainty and overstates significance -- the same pseudoreplication
  # issue as the EU-STORM nucleus script's ROI-vs-cell distinction. This is
  # a hard requirement now: cell_id must be present for these two things to
  # be computed correctly.
  if (!"cell_id" %in% names(df)) {
    stop("cell_id is required for cell-level CI/significance testing (see comments above ",
         "make_bracket... this section) and isn't present in this data.")
  }
  
  # One row per (cell_id, is_wt, radius_bin), with the FRACTION of that
  # cell's own domains falling in that bin. .drop = FALSE keeps radius_bin
  # combinations a cell happens to have zero domains in (fraction = 0)
  # rather than silently dropping that cell from the bin's average.
  per_cell_frac <- df %>%
    dplyr::count(cell_id, is_wt, radius_bin, .drop = FALSE) %>%
    dplyr::group_by(cell_id, is_wt) %>%
    dplyr::mutate(cell_total = sum(n), fraction = n / cell_total) %>%
    dplyr::ungroup()
  
  # ---- 1) CI on EACH bar: one-sample question, per condition ----
  # Mean +/- 95% CI of the per-cell fraction, across that condition's own
  # cells (t-distribution, appropriate for the likely-small number of cells
  # per condition). Rescaled to the count scale (x N_total) for display
  # alongside the bar height, same as before.
  bar_counts <- bar_counts %>%
    dplyr::mutate(N_total = ifelse(bar_key == "Wildtype", sum(df$is_wt), sum(!df$is_wt))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cell_fracs = list({
      is_wt_group <- (bar_key == "Wildtype")
      per_cell_frac$fraction[per_cell_frac$is_wt == is_wt_group & per_cell_frac$radius_bin == radius_bin]
    })) %>%
    dplyr::mutate(n_cells_group = length(cell_fracs),
                  ci = list(if (n_cells_group >= 2) t.test(cell_fracs)$conf.int else c(NA_real_, NA_real_))) %>%
    dplyr::mutate(ci_lo = ci[1] * N_total, ci_hi = ci[2] * N_total) %>%
    dplyr::ungroup()
  
  # ---- 2) Per-bin significance test: two-sample question, WT vs treatment ----
  # t-test comparing the per-cell fractions between Wildtype cells and
  # treatment cells, within each bin -- this is the two-group comparison,
  # answering "does the fraction of domains in this bin differ between
  # conditions," properly powered by NUMBER OF CELLS, not number of domains.
  bin_tests <- purrr::map_dfr(radius_bin_order, function(b) {
    x <- per_cell_frac$fraction[per_cell_frac$is_wt & per_cell_frac$radius_bin == b]
    y <- per_cell_frac$fraction[!per_cell_frac$is_wt & per_cell_frac$radius_bin == b]
    test <- suppressWarnings(t.test(x, y))
    tibble::tibble(radius_bin = factor(b, levels = radius_bin_order), p_value = test$p.value,
                   n_cells_wt = length(x), n_cells_treatment = length(y))
  })
  message("  Per-bin significance test (Wildtype vs ", treatment_label,
          "), two-sample t-test on per-cell fractions, no multiple-comparison correction:")
  print(bin_tests)
  
  bin_tests <- bin_tests %>%
    dplyr::left_join(
      bar_counts %>% dplyr::group_by(radius_bin) %>%
        dplyr::summarise(y_top = max(ci_hi, n, na.rm = TRUE), .groups = "drop"),
      by = "radius_bin"
    ) %>%
    dplyr::mutate(
      label = ifelse(p_value < 0.0001, "p < 0.0001", paste0("p = ", formatC(p_value, digits = 3, format = "g")))
    )
  # Fixed ABSOLUTE clearance above each bin's own bar top, sized off the
  # GLOBAL max across all bins -- not each bin's own height. The previous
  # version used y_top * 1.12 (proportional to that bin's OWN height),
  # which gave short bins almost no real clearance (e.g. 12% of a ~1200
  # bar is only ~140 units) while the count label above it needs roughly
  # the same absolute space regardless of bar height -- the two collided
  # and rendered as garbled overlapping text for every bin except the
  # tallest one.
  global_y_max <- max(bar_counts$ci_hi, bar_counts$n, na.rm = TRUE)
  bin_tests <- bin_tests %>%
    dplyr::mutate(y_label = y_top + global_y_max * 0.12)
  
  p_bar <- ggplot(bar_counts, aes(x = radius_bin, y = n, fill = bar_key)) +
    geom_col(position = position_dodge(width = 0.75), color = "black", linewidth = 1, width = 0.65) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(width = 0.75),
                  width = 0.2, linewidth = 0.7, na.rm = TRUE) +
    geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.4,
              size = 4, fontface = "bold") +
    geom_text(data = bin_tests, aes(x = radius_bin, y = y_label, label = label),
              inherit.aes = FALSE, size = 3.6, fontface = "bold") +
    scale_fill_manual(values = bar_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = paste0("Wildtype vs ", treatment_label), x = NULL, y = "Number of domains", fill = NULL) +
    theme_pub(legend_position = "none")
  
  (p_wt | p_treat | p_bar) +
    patchwork::plot_annotation(title = title_prefix)
}

# ============================================================================
# Run (printed, not auto-saved)
# ============================================================================

print(plot_wt_treatment_radius(domain_data$edu, "EdU-STORM", condition_values$edu$treatment))
print(plot_wt_treatment_radius(domain_data$k9,  "H3K9-STORM", condition_values$k9$treatment))

# ============================================================================
# Density KDE curves (faceted by radius bin) + 6-line log-log best-fit plot
#
# Depends on:
#   - domain_data, condition_values, theme_pub() from
#     domain_classification_plot.R
#   - add_wt_radius_bins(), RADIUS_BIN_COLORS, radius_bin_order from
#     domain_wt_vs_treatment_radius.R
# Run both of those first (in that order), then this file.
# ============================================================================

library(dplyr)
library(ggplot2)
library(scales)

# ============================================================================
# 1. KDE / probability density curves of `density`, faceted by radius bin
# ============================================================================

#' Probability density (KDE) curves of the `density` metric, one facet per
#' radius bin (Small/Average/Large), Wildtype and treatment overlaid within
#' each facet. Wildtype stays grey (consistent with every other plot in
#' this project); treatment is colored by that bin's own color.
plot_density_kde_by_bin <- function(raw_df, treatment_label) {
  df <- add_wt_radius_bins(raw_df) %>%
    dplyr::mutate(fill_key = factor(ifelse(is_wt, "Wildtype", as.character(radius_bin)),
                                    levels = c("Wildtype", radius_bin_order)))
  fill_colors <- c(Wildtype = "grey60", RADIUS_BIN_COLORS)
  
  ggplot(df, aes(x = density, color = fill_key, fill = fill_key)) +
    geom_density(alpha = 0.35, linewidth = 0.9) +
    facet_wrap(~radius_bin, scales = "free_y") +
    scale_color_manual(values = fill_colors) +
    scale_fill_manual(values = fill_colors) +
    labs(x = expression("Density (1/nm"^2*")"), y = "Probability density",
         color = NULL, fill = NULL,
         title = paste0("Density distribution by radius bin: Wildtype vs ", treatment_label)) +
    theme_pub(legend_position = "bottom")
}

# ============================================================================
# 2. 6-line log-log best-fit plot: one linear fit per (bin, condition),
#    no raw points
# ============================================================================

#' Density vs radius, log-log linear fit only (no points), one line per
#' radius bin x condition combination (3 bins x 2 conditions = 6 lines),
#' all on one shared panel.
#'
#' Color = radius bin (matches the scatter's category colors); linetype =
#' condition (dashed = Wildtype, solid = treatment). Note this is a
#' deliberate departure from the "Wildtype = grey" convention used
#' elsewhere in this project: with 6 lines sharing one panel, three grey
#' Wildtype lines would be indistinguishable from each other, so condition
#' is encoded via linetype instead and color is reserved for the bin.
#'
#' Confirmed (see chat) that scale_x_log10()/scale_y_log10() cause
#' geom_smooth(method = "lm") to fit the regression IN log10-log10 space
#' (a proper power-law fit), not on the raw scale merely re-displayed on
#' log axes -- verified numerically against a known power-law slope before
#' writing this.
plot_density_radius_fits <- function(raw_df, treatment_label) {
  df <- add_wt_radius_bins(raw_df) %>%
    dplyr::mutate(condition_label = factor(ifelse(is_wt, "Wildtype", treatment_label),
                                           levels = c("Wildtype", treatment_label)))
  linetypes <- setNames(c("dashed", "solid"), c("Wildtype", treatment_label))
  
  ggplot(df, aes(x = domain_radius_nm, y = density, color = radius_bin, linetype = condition_label)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    scale_x_log10(labels = label_number()) +
    scale_y_log10(labels = label_number()) +
    scale_color_manual(values = RADIUS_BIN_COLORS) +
    scale_linetype_manual(values = linetypes) +
    labs(x = "Radius (nm)", y = expression("Density (1/nm"^2*")"),
         color = "Radius bin", linetype = "Condition",
         title = paste0("Density vs Radius: log-log linear fit per bin, Wildtype vs ", treatment_label)) +
    theme_pub(legend_position = "right")
}

# ============================================================================
# Run (printed, not auto-saved)
# ============================================================================

print(plot_density_kde_by_bin(domain_data$edu, condition_values$edu$treatment))
print(plot_density_kde_by_bin(domain_data$k9,  condition_values$k9$treatment))

print(plot_density_radius_fits(domain_data$edu, condition_values$edu$treatment))
print(plot_density_radius_fits(domain_data$k9,  condition_values$k9$treatment))

# ============================================================================
# KDE + matching bar count, arranged as a 2x3 grid (one column per radius
# bin: Small/Average/Large), density curve on top, domain count bar below.
#
# Depends on:
#   - domain_data, condition_values, theme_pub() from
#     domain_classification_plot.R
#   - add_wt_radius_bins(), RADIUS_BIN_COLORS, radius_bin_order from
#     domain_wt_vs_treatment_radius.R
# Run both of those first, then this file.
#
# Median lines are drawn on the TREATMENT density curve only (not
# Wildtype), colored to match that bin's own color -- per request.
# ============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scales)

#' One column's KDE panel: Wildtype (grey) + treatment (bin color) density
#' curves for a single radius bin, with a dashed median line on the
#' treatment curve only, colored to match it. Reports both medians as text
#' in the upper-right corner.
plot_kde_one_bin <- function(df, bin_name, treatment_label, show_y_title = FALSE) {
  d <- df %>% dplyr::filter(radius_bin == bin_name)
  bin_color <- RADIUS_BIN_COLORS[[bin_name]]
  wt_median <- median(d$density[d$is_wt], na.rm = TRUE)
  treatment_median <- median(d$density[!d$is_wt], na.rm = TRUE)
  median_label <- paste0("WT: ", signif(wt_median, 3), "\n", treatment_label, ": ", signif(treatment_median, 3))
  
  ggplot(d, aes(x = density, color = is_wt, fill = is_wt)) +
    geom_density(alpha = 0.35, linewidth = 0.9) +
    geom_vline(xintercept = treatment_median, color = bin_color, linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = Inf, y = Inf, label = median_label, hjust = 1.05, vjust = 1.2,
             size = 3, lineheight = 0.9) +
    scale_color_manual(values = c(`TRUE` = "grey60", `FALSE` = bin_color), guide = "none") +
    scale_fill_manual(values = c(`TRUE` = "grey60", `FALSE` = bin_color), guide = "none") +
    labs(title = bin_name, x = expression("Density (1/nm"^2*")"),
         y = if (show_y_title) "Probability density" else NULL) +
    theme_pub(legend_position = "none") +
    theme(axis.title.x = element_text(size = rel(0.6), face = "plain"))
}

#' One column's bar panel: Wildtype vs treatment domain count, single
#' radius bin, with a cell-level CI and the per-bin significance test
#' already established in domain_wt_vs_treatment_radius.R (same logic,
#' scoped to one bin so it can sit directly under that bin's KDE panel).
plot_bar_one_bin <- function(df, bin_name, treatment_label) {
  bin_color <- RADIUS_BIN_COLORS[[bin_name]]
  per_cell_count <- df %>%
    dplyr::filter(radius_bin == bin_name) %>%
    dplyr::count(cell_id, is_wt, name = "domain_count")
  # cells with genuinely zero domains in this bin need to be represented
  # too (count = 0), or they'd be silently dropped from that condition's
  # mean/CI -- reconstruct the full cell list per condition and left_join.
  all_cells <- df %>% dplyr::distinct(cell_id, is_wt)
  per_cell_count <- all_cells %>%
    dplyr::left_join(per_cell_count, by = c("cell_id", "is_wt")) %>%
    dplyr::mutate(domain_count = tidyr::replace_na(domain_count, 0))
  
  summarise_group <- function(is_wt_val) {
    x <- per_cell_count$domain_count[per_cell_count$is_wt == is_wt_val]
    ci <- if (length(x) >= 2) t.test(x)$conf.int else c(NA_real_, NA_real_)
    tibble::tibble(is_wt = is_wt_val, n = sum(x), n_cells = length(x),
                   ci_lo = ci[1] * length(x), ci_hi = ci[2] * length(x))
  }
  bar_data <- dplyr::bind_rows(summarise_group(TRUE), summarise_group(FALSE)) %>%
    dplyr::mutate(group_label = factor(ifelse(is_wt, "Wildtype", treatment_label),
                                       levels = c("Wildtype", treatment_label)))
  
  x_wt <- per_cell_count$domain_count[per_cell_count$is_wt]
  x_treatment <- per_cell_count$domain_count[!per_cell_count$is_wt]
  p_value <- suppressWarnings(t.test(x_wt, x_treatment)$p.value)
  p_label <- if (p_value < 0.0001) "p < 0.0001" else paste0("p = ", formatC(p_value, digits = 3, format = "g"))
  
  bar_colors <- setNames(c("grey60", bin_color), c("Wildtype", treatment_label))
  
  ggplot(bar_data, aes(x = group_label, y = n, fill = group_label)) +
    geom_col(color = "black", linewidth = 1, width = 0.6, show.legend = FALSE) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.7, na.rm = TRUE) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.8, fontface = "bold") +
    annotate("text", x = 1.5, y = max(bar_data$ci_hi, bar_data$n, na.rm = TRUE) * 1.15,
             label = p_label, size = 3.2, fontface = "bold") +
    scale_fill_manual(values = bar_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
    labs(x = NULL, y = "Number of domains") +
    theme_pub(legend_position = "none")
}

#' Classify every domain (both conditions) into a WT-derived radius
#' tercile: break points are the 33rd/66th percentile of ONLY the
#' Wildtype rows' domain_radius_nm, applied to both WT and treatment.
#' This is the SAME approach used briefly in an earlier version of
#' domain_wt_vs_treatment_radius.R before that file was switched to fixed
#' 80/120nm thresholds for consistency across the project -- brought back
#' here as a K9-specific second version, per request, rather than
#' replacing the fixed-threshold version already in that file.
add_wt_radius_tertiles <- function(df) {
  wt_radius <- df$domain_radius_nm[df$is_wt]
  break_lo <- as.numeric(quantile(wt_radius, 1/3, na.rm = TRUE))
  break_hi <- as.numeric(quantile(wt_radius, 2/3, na.rm = TRUE))
  message("  WT-derived radius tercile breaks (K9): ", round(break_lo, 1), " nm, ", round(break_hi, 1), " nm",
          " (vs the fixed 80nm/120nm breaks used elsewhere in this project)")
  
  df %>%
    dplyr::mutate(radius_bin = dplyr::case_when(
      domain_radius_nm <  break_lo ~ "Small radius",
      domain_radius_nm >= break_lo & domain_radius_nm < break_hi ~ "Average radius",
      domain_radius_nm >= break_hi ~ "Large radius",
      TRUE ~ NA_character_
    )) %>%
    dplyr::mutate(radius_bin = factor(radius_bin, levels = radius_bin_order))
}

#' Classify every domain (both conditions) into a FIXED radius bin at
#' 40nm/80nm -- a third, K9-specific threshold scheme (<=40nm / 40-80nm /
#' >80nm), distinct from both the 80/120nm fixed thresholds and the
#' Wildtype-tercile thresholds above.
add_radius_bins_40_80 <- function(df) {
  df %>%
    dplyr::mutate(radius_bin = dplyr::case_when(
      domain_radius_nm <= 40 ~ "Small radius",
      domain_radius_nm >  40 & domain_radius_nm <= 80 ~ "Average radius",
      domain_radius_nm >  80 ~ "Large radius",
      TRUE ~ NA_character_
    )) %>%
    dplyr::mutate(radius_bin = factor(radius_bin, levels = radius_bin_order))
}

#' Full 2x3 grid: KDE row on top, matching bar row below, one column per
#' radius bin. `binning_fn` selects how radius_bin gets computed --
#' add_wt_radius_bins (fixed 80/120nm thresholds), add_wt_radius_tertiles
#' (K9 Wildtype's own terciles), or add_radius_bins_40_80 (fixed
#' 40/80nm thresholds) -- everything downstream (KDE panels, bar panels,
#' CI, significance tests) is identical either way, since those functions
#' only look at the resulting radius_bin column, not how it was derived.
plot_kde_bar_grid <- function(raw_df, treatment_label, title, binning_fn = add_wt_radius_bins) {
  df <- binning_fn(raw_df)
  
  kde_panels <- purrr::imap(radius_bin_order, function(b, i) {
    plot_kde_one_bin(df, b, treatment_label, show_y_title = (i == 1))
  })
  bar_panels <- purrr::map(radius_bin_order, function(b) {
    plot_bar_one_bin(df, b, treatment_label)
  })
  
  wrap_plots(c(kde_panels, bar_panels), ncol = 3, nrow = 2, byrow = TRUE) +
    patchwork::plot_annotation(title = title)
}

# ============================================================================
# Run (printed, not auto-saved) -- K9 only, per request
# ============================================================================

# Version 1: fixed 80nm/120nm thresholds (same as the rest of this project)
print(plot_kde_bar_grid(domain_data$k9, condition_values$k9$treatment,
                        "H3K9-STORM (fixed 80nm/120nm thresholds)",
                        binning_fn = add_wt_radius_bins))

# Version 2: K9 Wildtype's own tercile breaks
print(plot_kde_bar_grid(domain_data$k9, condition_values$k9$treatment,
                        "H3K9-STORM (Wildtype-tercile thresholds)",
                        binning_fn = add_wt_radius_tertiles))

# Version 3: fixed 40nm/80nm thresholds
print(plot_kde_bar_grid(domain_data$k9, condition_values$k9$treatment,
                        "H3K9-STORM (fixed 40nm/80nm thresholds)",
                        binning_fn = add_radius_bins_40_80))

# ============================================================================
# WT scatter | Treatment scatter | [KDE row / Bar row] x 3 columns, combined
#
# Replaces the simple single bar chart on the right of the WT/Treatment
# scatter composite (domain_wt_vs_treatment_radius.R) with the full faceted
# KDE + count grid (domain_kde_bar_grid.R) instead -- per request.
#
# Uses the FIXED 40nm/80nm thresholds (add_radius_bins_40_80 from
# domain_kde_bar_grid.R) -- not the 80/120nm or Wildtype-tercile schemes
# also defined there -- per request.
#
# Depends on:
#   - domain_data, condition_values, theme_pub() from
#     domain_classification_plot.R
#   - RADIUS_BIN_COLORS, radius_bin_order from
#     domain_wt_vs_treatment_radius.R
#   - add_radius_bins_40_80(), plot_kde_one_bin(), plot_bar_one_bin() from
#     domain_kde_bar_grid.R
# Run all three (in that order), then this file.
#
# K9 only, matching how domain_kde_bar_grid.R was scoped.
# ============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

#' WT and treatment scatter panels (density vs radius, log-log, colored by
#' radius bin for treatment / plain grey for WT). Rebuilt fresh here rather
#' than reused from domain_wt_vs_treatment_radius.R, since that file builds
#' these as a nested/local function inside plot_wt_treatment_radius() and
#' doesn't expose them separately. `binning_fn` and its matching
#' `break_lo`/`break_hi` control which threshold scheme's dashed lines are
#' drawn, so this and the KDE/bar grid stay visually consistent.
build_scatter_panels <- function(raw_df, treatment_label, binning_fn, break_lo, break_hi) {
  df <- binning_fn(raw_df)
  
  scatter_base <- function(d, title) {
    ggplot(d, aes(x = domain_radius_nm, y = density)) +
      geom_vline(xintercept = c(break_lo, break_hi), linetype = "dashed", color = "grey50", linewidth = 0.5) +
      scale_x_log10(labels = label_number()) +
      scale_y_log10(labels = label_number()) +
      labs(title = title, x = "Radius (nm)", y = expression("Density (1/nm"^2*")")) +
      theme_pub(legend_position = "none")
  }
  
  p_wt <- scatter_base(df %>% dplyr::filter(is_wt), "WT") +
    geom_point(color = "grey60", alpha = 0.7, size = 1.4)
  p_treat <- scatter_base(df %>% dplyr::filter(!is_wt), treatment_label) +
    geom_point(aes(color = radius_bin), alpha = 0.7, size = 1.4) +
    scale_color_manual(values = RADIUS_BIN_COLORS)
  
  list(p_wt = p_wt, p_treat = p_treat)
}

#' Full composite: WT scatter | Treatment scatter | (3-column KDE-over-bar
#' grid). binning_fn/break_lo/break_hi flow through to both halves, so the
#' scatter panels' dashed lines and the KDE/bar grid's bins always match.
plot_scatter_kde_bar_combo <- function(raw_df, treatment_label, title,
                                       binning_fn = add_radius_bins_40_80,
                                       break_lo = 40, break_hi = 80) {
  scatters <- build_scatter_panels(raw_df, treatment_label, binning_fn, break_lo, break_hi)
  df <- binning_fn(raw_df)
  
  kde_panels <- purrr::imap(radius_bin_order, function(b, i) {
    plot_kde_one_bin(df, b, treatment_label, show_y_title = (i == 1))
  })
  bar_panels <- purrr::map(radius_bin_order, function(b) {
    plot_bar_one_bin(df, b, treatment_label)
  })
  grid_part <- wrap_plots(c(kde_panels, bar_panels), ncol = 3, nrow = 2, byrow = TRUE)
  
  (scatters$p_wt | scatters$p_treat | grid_part) +
    plot_layout(widths = c(1, 1, 2.6)) +
    patchwork::plot_annotation(title = title)
}

# ============================================================================
# Run (printed, not auto-saved) -- K9 only, fixed 40nm/80nm thresholds
# ============================================================================

print(plot_scatter_kde_bar_combo(domain_data$k9, condition_values$k9$treatment, "H3K9-STORM"))
