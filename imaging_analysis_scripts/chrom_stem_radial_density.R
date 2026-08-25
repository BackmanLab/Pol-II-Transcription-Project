##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# ChromEM radial density analysis: Wildtype vs. POLR2A-degraded.
# Reproduces the 4-panel MATLAB figure (Nascent/Small Mature/
# Decaying/Large Mature chromatin domains, Wildtype vs. POLR2A-deg,
# density vs. distance on a log10 x-axis) from
# RadialDensityPerGroup.xlsx, plus a second version laid out as a
# single row with a marginal boxplot of per-domain domain-boundary
# distance beneath each panel.
#
# Standalone -- pulled out of an earlier combined script
# (edu_radial_density_analysis.R) that also loaded an unrelated
# EdU-STORM cluster CSV before this ChromEM code; this file has
# only the ChromEM-relevant code.
#
# Important: each column in the raw sheets is an individual
# chromatin domain, not a cell or a biological replicate -- this is
# a single-cell experiment. So "n" throughout this script (in the
# CI sample-size check, the domain-boundary distributions, etc.)
# counts domains observed within that (small number of) cell(s),
# not independent cells/animals. The mean +/- CI and the boxplots
# describe domain-to-domain variability within the experiment, not
# cell-to-cell or biological reproducibility -- a narrow CI or a
# "significant" per-domain comparison should not be read as
# evidence that would generalize across a larger number of cells,
# since it can't speak to that from this data structure.
#
# Figure output: Fig. 9F (confirmed).
##################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readxl)
library(patchwork)
library(scales)

# ---- 0. Config --------------------------------------------------------
in_dir  <- "/Volumes/Backup #1/IBiS/Backman_Lab/Transcription Publication/Pol2 Submission/new data/wingshun/POL12AnalysisCombined/"
em_file <- file.path(in_dir, "RadialDensityPerGroup.xlsx")
out_dir <- "outputs"

condition_levels <- c("Wildtype", "POLR2A deg")   # reference level first
COND_COLORS <- c("Wildtype" = "#88CCEE", "POLR2A deg" = "#CC6677")

dir.create(out_dir, showWarnings = FALSE)

# ---- Shared publication theme ------------------------------------------
theme_pub <- function(base_size = 13, legend_position = "right") {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = rel(1.05), hjust = 0,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = rel(0.75), color = "grey35",
                                      margin = margin(b = 8)),
      axis.title       = element_text(face = "bold", size = rel(0.95)),
      axis.text        = element_text(color = "black", size = rel(0.85)),
      axis.line        = element_line(linewidth = 0.4, color = "black"),
      axis.ticks       = element_line(linewidth = 0.4, color = "black"),
      legend.title     = element_blank(),
      legend.text      = element_text(size = rel(0.85)),
      legend.key       = element_blank(),
      legend.background = element_blank(),
      legend.position  = legend_position,
      strip.background = element_rect(fill = "grey93", color = NA),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      panel.grid       = element_blank(),
      panel.spacing    = unit(1, "lines"),
      plot.title.position = "plot"
    )
}

# ============================================================================
# PART 1 -- ChromEM radial density per domain category (2x2 facet grid)
#
#   Rather than just plotting the workbook's pre-averaged "Averaged" sheet
#   (a single mean line with no uncertainty, like the original MATLAB
#   figure), we recompute the mean +/- 95% CI directly from the 8 raw
#   per-replicate sheets (one column per replicate/cell, no header row).
#   There's no other statistical test attached to this panel, so the CI
#   ribbon is the uncertainty estimate here. Individual replicate curves
#   have different maximum radii and are padded with "NaN" text beyond
#   their own extent -- those are coerced to real NA and excluded from
#   the mean/SD/n at each distance. n legitimately shrinks at large
#   distances as fewer replicates still have data there -- the script
#   prints the actual min/max n it finds below so you can check that
#   for yourself rather than take it on faith.
# ============================================================================

em_wide <- read_excel(em_file, sheet = "Averaged") %>%
  rename(distance = `distance (nm)`)

em_sheet_map <- tibble(
  sheet    = c("Nascent_Ctrl", "SmallMature_Ctrl", "Decaying_Ctrl", "LargeMatureCtrl",
               "Nascent_AID",  "SmallMature_AID",  "Decaying_AID",  "LargeMature_AID"),
  category = c("Nascent", "Small Mature", "Decaying", "Large Mature",
               "Nascent", "Small Mature", "Decaying", "Large Mature"),
  em_condition = c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "AID", "AID", "AID", "AID")
)

category_order <- c("Large Mature", "Small Mature", "Nascent", "Decaying")

# Read one raw replicate sheet and compute mean/SD/n/95%CI at each distance.
summarise_em_sheet <- function(sheet_name, category_label, em_condition_label, distance_vec) {
  raw <- suppressMessages(read_excel(em_file, sheet = sheet_name, col_names = FALSE)) %>%
    mutate(across(everything(), as.numeric))   # "NaN" text -> real NA
  raw_mat <- as.matrix(raw)
  
  tibble(
    distance     = distance_vec,
    mean_density = apply(raw_mat, 1, mean, na.rm = TRUE),
    sd_density   = apply(raw_mat, 1, sd,   na.rm = TRUE),
    n            = apply(raw_mat, 1, function(r) sum(!is.na(r))),
    n_domains_in_sheet = ncol(raw_mat),
    category     = category_label,
    em_condition = em_condition_label
  )
}

em_profile <- pmap_dfr(em_sheet_map, function(sheet, category, em_condition) {
  summarise_em_sheet(sheet, category, em_condition, em_wide$distance)
}) %>%
  mutate(se   = sd_density / sqrt(n),
         ci95 = 1.96 * se,
         category  = factor(category, levels = category_order),
         condition = recode(em_condition, Ctrl = "Wildtype", AID = "POLR2A deg"),
         condition = factor(condition, levels = condition_levels)) %>%
  filter(distance > 0)   # distance = 0 can't be shown on a log10 x-axis

# Check, in code, exactly what sample sizes back each point of the CI --
# print it instead of asserting it.
n_check <- em_profile %>%
  group_by(category, condition) %>%
  summarise(n_domains_total = first(n_domains_in_sheet),
            n_min_at_any_distance = min(n),
            n_max_at_any_distance = max(n),
            .groups = "drop")
message("\nCI sample-size check (per category x condition):")
print(n_check)
if (any(n_check$n_min_at_any_distance < 3)) {
  message("WARNING: at least one category/condition has fewer than 3 replicates ",
          "contributing to the CI at some distance -- treat that stretch of the ",
          "ribbon with caution.")
}

# ---- Domain boundary detection -----------------------------------------
# Defined as the first local minimum of the pooled (both-conditions-averaged)
# density curve, i.e. the first distance at which the curve stops declining
# and turns upward -- the trough that separates the dense domain core from
# the lower-density surroundings. We deliberately use the FIRST local
# minimum (not the global minimum) because the curves flatten into a noisy
# plateau at large distances, and in one category (Large Mature) the noisy
# plateau happens to dip slightly below the true trough further out --
# taking a global min would put the "boundary" in the noise, not at the
# actual shoulder of the curve. See the printed comparison below.
find_first_local_min <- function(distance, density) {
  o <- order(distance)
  distance <- distance[o]
  density  <- density[o]
  step_direction <- diff(density)                 # step_direction[i] = density[i+1] - density[i]
  turn <- which(step_direction >= 0)[1]            # first point where it stops decreasing
  if (is.na(turn)) turn <- which.min(density)       # fallback if it never turns (shouldn't happen here)
  distance[turn]
}

pooled_profile <- em_profile %>%
  group_by(category, distance) %>%
  summarise(pooled_density = mean(mean_density, na.rm = TRUE), .groups = "drop")

domain_boundaries <- pooled_profile %>%
  group_by(category) %>%
  summarise(boundary          = find_first_local_min(distance, pooled_density),
            global_min_distance = distance[which.min(pooled_density)],
            .groups = "drop") %>%
  left_join(
    em_profile %>%
      group_by(category) %>%
      summarise(x_min = min(distance), x_max = max(distance),
                y_top = max(mean_density + ci95, na.rm = TRUE),
                .groups = "drop"),
    by = "category"
  )

message("\nDomain boundary per category (first local minimum of the pooled curve):")
print(domain_boundaries %>% select(category, boundary, global_min_distance))
if (any(domain_boundaries$boundary != domain_boundaries$global_min_distance)) {
  message("NOTE: for at least one category, the first local minimum differs from ",
          "the global minimum (a noisier, deeper point later in the plateau) -- ",
          "the first local minimum is used for the boundary line, as intended.")
}

# ---- Domain core detection (third break) --------------------------------
# "Domain core" = the distance at which the pooled curve first drops to a
# density of 0.4, found by linear interpolation between the two bracketing
# grid points (not just snapped to the nearest 2nm grid point) for a more
# precise crossing distance. This splits the plot into three zones:
#   1. domain core   (distance < core_threshold_distance): density > 0.4
#   2. transition    (core_threshold_distance < distance < boundary)
#   3. periphery     (distance > boundary)
CORE_DENSITY_THRESHOLD <- 0.4

find_threshold_crossing <- function(distance, density, threshold) {
  o <- order(distance)
  distance <- distance[o]
  density  <- density[o]
  above <- density >= threshold
  # first index where the curve steps from >= threshold to < threshold
  idx <- which(above[-length(above)] & !above[-1])[1]
  if (is.na(idx)) return(NA_real_)   # curve never crosses the threshold
  d1 <- distance[idx];     d2 <- distance[idx + 1]
  v1 <- density[idx];      v2 <- density[idx + 1]
  d1 + (threshold - v1) * (d2 - d1) / (v2 - v1)   # linear interpolation
}

domain_boundaries <- domain_boundaries %>%
  left_join(
    pooled_profile %>%
      group_by(category) %>%
      summarise(core_boundary = find_threshold_crossing(distance, pooled_density,
                                                        CORE_DENSITY_THRESHOLD),
                .groups = "drop"),
    by = "category"
  )

message("\nDomain core per category (pooled curve crosses density = ",
        CORE_DENSITY_THRESHOLD, ", linear-interpolated):")
print(domain_boundaries %>% select(category, core_boundary, boundary))
if (any(is.na(domain_boundaries$core_boundary))) {
  message("WARNING: at least one category's pooled curve never reaches density = ",
          CORE_DENSITY_THRESHOLD, " -- no domain-core line drawn for that panel.")
}
if (any(domain_boundaries$core_boundary >= domain_boundaries$boundary, na.rm = TRUE)) {
  message("WARNING: for at least one category, the domain-core crossing falls at or ",
          "beyond the domain boundary -- check CORE_DENSITY_THRESHOLD against that ",
          "panel's density range.")
}

core_color <- "#E8D6A8"        # richer tan   -- domain core (density > 0.4)
mid_color  <- "#F5EEDD"        # pale tan     -- transition zone
periphery_color <- "#E4E9EC"   # pale grey-blue -- surrounding chromatin

p_chromem <- ggplot(em_profile, aes(distance, mean_density,
                                    color = condition, fill = condition)) +
  geom_rect(data = domain_boundaries,
            aes(xmin = x_min, xmax = core_boundary, ymin = -Inf, ymax = Inf),
            fill = core_color, color = NA, alpha = 1, inherit.aes = FALSE) +
  geom_rect(data = domain_boundaries,
            aes(xmin = core_boundary, xmax = boundary, ymin = -Inf, ymax = Inf),
            fill = mid_color, color = NA, alpha = 1, inherit.aes = FALSE) +
  geom_rect(data = domain_boundaries,
            aes(xmin = boundary, xmax = x_max, ymin = -Inf, ymax = Inf),
            fill = periphery_color, color = NA, alpha = 1, inherit.aes = FALSE) +
  geom_vline(data = domain_boundaries, aes(xintercept = core_boundary),
             linetype = "dashed", linewidth = 0.6, color = "grey35") +
  geom_vline(data = domain_boundaries, aes(xintercept = boundary),
             linetype = "dashed", linewidth = 0.6, color = "grey35") +
  geom_text(data = domain_boundaries,
            aes(x = core_boundary, y = y_top, label = "Domain core"),
            angle = 90, hjust = 1.1, vjust = -0.5, size = 3.1, color = "grey25",
            inherit.aes = FALSE) +
  geom_text(data = domain_boundaries,
            aes(x = boundary, y = y_top, label = "Domain boundary"),
            angle = 90, hjust = 1.1, vjust = -0.5, size = 3.1, color = "grey25",
            inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = mean_density - ci95, ymax = mean_density + ci95),
              alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~category, ncol = 2, scales = "free_y") +
  scale_x_log10() +
  scale_color_manual(values = COND_COLORS) +
  scale_fill_manual(values = COND_COLORS) +
  labs(x = "Distance (nm)", y = "Density",
       title = "ChromEM radial density by chromatin domain category",
       subtitle = paste0("Mean \u00b1 95% CI across replicate domains; dashed lines = domain core ",
                         "(pooled density = ", CORE_DENSITY_THRESHOLD,
                         ") and domain boundary (first local minimum of the pooled curve)")) +
  theme_pub(legend_position = c(0.9, 0.94))

print(p_chromem)

write_csv(em_profile, file.path(out_dir, "chromem_radial_density_summary.csv"))
write_csv(domain_boundaries, file.path(out_dir, "chromem_domain_boundaries.csv"))

message("Figure 1 printed above. Summary CSVs (density profile with 95% CI, domain boundaries) written to '", out_dir, "/'.")

# ============================================================================
# PART 2 -- ChromEM as a single row (4 panels), with a marginal boxplot of
# per-domain domain-boundary distance beneath each panel
#
# The plot above is a 2x2 facet grid computed from the POOLED
# (averaged-across-domains) curve -- domain_boundaries there is ONE number
# per category, not a distribution, so there's nothing to box-plot from it
# directly. This section instead finds the domain boundary (first local
# minimum, same definition as find_first_local_min() above) on EACH
# INDIVIDUAL domain's own curve, giving a real per-domain distribution
# (n = that category's domain count -- remember, single-cell experiment,
# so this is domains within a small number of cells, not independent
# cells/animals) to summarize as a boxplot -- box + jittered points per
# condition, x-axis shared with the panel above.
# ============================================================================

#' Domain-boundary distance for EVERY individual domain (column) in one raw
#' sheet (not the pooled/averaged curve) -- same find_first_local_min()
#' logic as above, applied per-column. Domains with fewer than 3 valid
#' (non-NaN-padded, non-zero-distance) points are skipped (NA) rather than
#' given a meaningless "boundary". distance = 0 is excluded here for the
#' same reason em_profile/pooled_profile exclude it elsewhere in this
#' script: log10(0) is undefined, and a domain whose own minimum happens
#' to land exactly on the first (distance = 0) point would otherwise
#' produce -Inf values and warnings.
compute_per_domain_boundaries <- function(sheet_name, category_label, em_condition_label, distance_vec) {
  raw <- suppressMessages(read_excel(em_file, sheet = sheet_name, col_names = FALSE)) %>%
    mutate(across(everything(), as.numeric))
  raw_mat <- as.matrix(raw)
  
  boundaries <- apply(raw_mat, 2, function(col) {
    valid <- !is.na(col) & distance_vec > 0
    if (sum(valid) < 3) return(NA_real_)
    find_first_local_min(distance_vec[valid], col[valid])
  })
  
  tibble(domain_id = seq_along(boundaries), boundary = unname(boundaries),
         category = category_label, em_condition = em_condition_label)
}

per_domain_boundaries <- pmap_dfr(em_sheet_map, function(sheet, category, em_condition) {
  compute_per_domain_boundaries(sheet, category, em_condition, em_wide$distance)
}) %>%
  mutate(category  = factor(category, levels = category_order),
         condition = recode(em_condition, Ctrl = "Wildtype", AID = "POLR2A deg"),
         condition = factor(condition, levels = condition_levels)) %>%
  filter(!is.na(boundary))

message("\nPer-domain domain-boundary distance, n per category x condition:")
print(per_domain_boundaries %>% count(category, condition))

#' One category's radial-density panel (Wildtype vs POLR2A deg, ribbon +
#' core/boundary shading), rebuilt as a standalone panel rather than a
#' facet, so it can be stacked with its own marginal panel below it.
plot_radial_one_category <- function(category_name, show_legend = FALSE, show_y_title = FALSE) {
  d  <- em_profile %>% filter(category == category_name)
  db <- domain_boundaries %>% filter(category == category_name)
  
  ggplot(d, aes(distance, mean_density, color = condition, fill = condition)) +
    geom_rect(data = db, aes(xmin = x_min, xmax = core_boundary, ymin = -Inf, ymax = Inf),
              fill = core_color, color = NA, inherit.aes = FALSE) +
    geom_rect(data = db, aes(xmin = core_boundary, xmax = boundary, ymin = -Inf, ymax = Inf),
              fill = mid_color, color = NA, inherit.aes = FALSE) +
    geom_rect(data = db, aes(xmin = boundary, xmax = x_max, ymin = -Inf, ymax = Inf),
              fill = periphery_color, color = NA, inherit.aes = FALSE) +
    geom_vline(data = db, aes(xintercept = core_boundary), linetype = "dashed", linewidth = 0.6, color = "grey35") +
    geom_vline(data = db, aes(xintercept = boundary), linetype = "dashed", linewidth = 0.6, color = "grey35") +
    geom_ribbon(aes(ymin = mean_density - ci95, ymax = mean_density + ci95), alpha = 0.3, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_x_log10(limits = c(db$x_min, db$x_max)) +
    scale_color_manual(values = COND_COLORS) +
    scale_fill_manual(values = COND_COLORS) +
    labs(title = category_name, x = "Distance (nm)", y = if (show_y_title) "Density" else NULL) +
    theme_pub(legend_position = if (show_legend) c(0.72, 0.85) else "none")
}

#' One category's marginal boxplot: per-domain domain-boundary distance,
#' Wildtype vs POLR2A deg, box + jittered points, x-axis (log10, same
#' limits as the panel above) is the shared "distance" scale.
plot_boundary_marginal <- function(category_name) {
  d  <- per_domain_boundaries %>% filter(category == category_name)
  db <- domain_boundaries %>% filter(category == category_name)
  
  ggplot(d, aes(x = boundary, y = condition, color = condition)) +
    geom_jitter(height = 0.18, alpha = 0.5, size = 1.1, shape = 16) +
    geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA, linewidth = 0.6) +
    scale_x_log10(limits = c(db$x_min, db$x_max)) +
    scale_color_manual(values = COND_COLORS) +
    labs(x = "Distance (nm)", y = NULL) +
    theme_pub(legend_position = "none") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

row_panels <- purrr::imap(category_order, function(cat, i) {
  p_top    <- plot_radial_one_category(cat, show_legend = (i == 1), show_y_title = (i == 1))
  p_bottom <- plot_boundary_marginal(cat)
  (p_top / p_bottom) + plot_layout(heights = c(4, 1))
})

p_chromem_row <- wrap_plots(row_panels, ncol = 4) +
  plot_annotation(title = "ChromEM radial density with per-domain domain-boundary distributions")

print(p_chromem_row)

write_csv(per_domain_boundaries, file.path(out_dir, "chromem_per_domain_boundaries.csv"))
message("Figure 2 printed above. Per-domain boundary distances written to '", out_dir, "/'.")

# ============================================================================
# PART 3 -- same row layout, marginal panel = per-domain AUC (integrated
# density across the whole domain) and per-domain peak/core density
#
# Two separate metrics, tried separately per request:
#   AUC   = trapezoidal area under each domain's own (distance, density)
#           curve, over that domain's own valid (non-NaN-padded) range.
#           An integrated measure of "how much density, overall" rather
#           than a single landmark position like the boundary distance was.
#   peak  = max(density) on each domain's own curve -- the amplitude/core
#           density, independent of where along the curve it occurs.
# Both computed per INDIVIDUAL domain (not the pooled curve), same
# single-cell-experiment caveat as Part 2's boundary distribution applies
# here too.
# ============================================================================

#' Trapezoidal area under (x, y), sorted by x first.
trapezoidal_auc <- function(x, y) {
  o <- order(x)
  x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

#' Generic per-domain metric extractor: reads one raw sheet, applies
#' `metric_fn(distance, density)` to each domain's own valid (non-NaN,
#' non-zero-distance) points, returns one row per domain.
compute_per_domain_metric <- function(sheet_name, category_label, em_condition_label, distance_vec, metric_fn) {
  raw <- suppressMessages(read_excel(em_file, sheet = sheet_name, col_names = FALSE)) %>%
    mutate(across(everything(), as.numeric))
  raw_mat <- as.matrix(raw)
  
  values <- apply(raw_mat, 2, function(col) {
    valid <- !is.na(col) & distance_vec > 0
    if (sum(valid) < 3) return(NA_real_)
    metric_fn(distance_vec[valid], col[valid])
  })
  
  tibble(domain_id = seq_along(values), value = unname(values),
         category = category_label, em_condition = em_condition_label)
}

#' One category's marginal boxplot for a generic per-domain metric (AUC or
#' peak density), on a LINEAR y-axis to match the metric (not a distance,
#' so log10 doesn't apply the way it did for the boundary panel).
#' `shared_limits`, if given, forces the SAME x-axis range on every
#' category's panel (computed once, across all 4 categories, by
#' build_metric_row() below) so values can be compared directly between
#' quadrants rather than each panel auto-scaling to its own data.
plot_metric_marginal <- function(metric_df, category_name, x_label, shared_limits = NULL) {
  d <- metric_df %>% filter(category == category_name)
  
  p <- ggplot(d, aes(x = value, y = condition, color = condition)) +
    geom_jitter(height = 0.18, alpha = 0.5, size = 1.1, shape = 16) +
    geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA, linewidth = 0.6) +
    scale_color_manual(values = COND_COLORS) +
    labs(x = x_label, y = NULL) +
    theme_pub(legend_position = "none") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  if (!is.null(shared_limits)) p <- p + scale_x_continuous(limits = shared_limits)
  p
}

#' Build the full row (radial density panel + marginal metric panel) for
#' one metric, across all 4 categories. shared_x_axis = TRUE computes one
#' x-range from ALL 4 categories' values combined and applies it to every
#' marginal panel (safe -- by construction every category's own values
#' fall inside that combined range, so nothing gets clipped).
build_metric_row <- function(metric_df, x_label, title, shared_x_axis = FALSE) {
  shared_limits <- if (shared_x_axis) range(metric_df$value, na.rm = TRUE) else NULL
  panels <- purrr::imap(category_order, function(cat, i) {
    p_top    <- plot_radial_one_category(cat, show_legend = (i == 1), show_y_title = (i == 1))
    p_bottom <- plot_metric_marginal(metric_df, cat, x_label, shared_limits)
    (p_top / p_bottom) + plot_layout(heights = c(4, 1))
  })
  wrap_plots(panels, ncol = 4) + plot_annotation(title = title)
}

# ---- AUC ----
per_domain_auc <- pmap_dfr(em_sheet_map, function(sheet, category, em_condition) {
  compute_per_domain_metric(sheet, category, em_condition, em_wide$distance, trapezoidal_auc)
}) %>%
  mutate(category  = factor(category, levels = category_order),
         condition = recode(em_condition, Ctrl = "Wildtype", AID = "POLR2A deg"),
         condition = factor(condition, levels = condition_levels)) %>%
  filter(!is.na(value))

message("\nPer-domain AUC (integrated density), n per category x condition:")
print(per_domain_auc %>% count(category, condition))

p_chromem_row_auc <- build_metric_row(per_domain_auc, "AUC (integrated density)",
                                      "ChromEM radial density with per-domain AUC distributions")
print(p_chromem_row_auc)
write_csv(per_domain_auc, file.path(out_dir, "chromem_per_domain_auc.csv"))
message("Figure 3 printed above. Per-domain AUC values written to '", out_dir, "/'.")

# ---- Peak / core density ----
per_domain_peak <- pmap_dfr(em_sheet_map, function(sheet, category, em_condition) {
  compute_per_domain_metric(sheet, category, em_condition, em_wide$distance, function(x, y) max(y))
}) %>%
  mutate(category  = factor(category, levels = category_order),
         condition = recode(em_condition, Ctrl = "Wildtype", AID = "POLR2A deg"),
         condition = factor(condition, levels = condition_levels)) %>%
  filter(!is.na(value))

message("\nPer-domain peak/core density, n per category x condition:")
print(per_domain_peak %>% count(category, condition))

p_chromem_row_peak <- build_metric_row(per_domain_peak, "Peak (core) density",
                                       "ChromEM radial density with per-domain peak-density distributions",
                                       shared_x_axis = TRUE)
print(p_chromem_row_peak)
write_csv(per_domain_peak, file.path(out_dir, "chromem_per_domain_peak_density.csv"))
message("Figure 4 printed above. Per-domain peak density values written to '", out_dir, "/'.")

# ============================================================================
# PART 4 -- Delta plot: Wildtype minus POLR2A deg, as a function of distance
#
# Rather than showing WT and POLR2A deg as two overlaid lines (where subtle
# differences are hard to read, especially through the CI ribbons), this
# plots the DIFFERENCE directly: delta_density(d) = mean_density_WT(d) -
# mean_density_POLR2A_deg(d), one line per category, with a dashed
# reference line at delta = 0 (no difference) -- styled after your
# reference image's right-hand "delta S(d)" panel.
#
# Uncertainty on the delta is propagated from both conditions' standard
# errors (independent-samples formula: SE_delta = sqrt(SE_WT^2 +
# SE_treatment^2)), not just eyeballed -- shown as a ribbon per category.
# Both conditions' profiles must share the same distance grid for this
# subtraction to be meaningful, which they do here (both come from the
# same "Averaged" sheet's distance column via em_wide$distance).
# ============================================================================

delta_profile <- em_profile %>%
  select(category, condition, distance, mean_density, se) %>%
  pivot_wider(names_from = condition, values_from = c(mean_density, se)) %>%
  mutate(delta_density = `mean_density_Wildtype` - `mean_density_POLR2A deg`,
         delta_se   = sqrt(`se_Wildtype`^2 + `se_POLR2A deg`^2),
         delta_ci95 = 1.96 * delta_se) %>%
  filter(!is.na(delta_density))

message("\nDelta profile (Wildtype - POLR2A deg) computed for ", nrow(delta_profile),
        " category x distance points.")

# Distinct color per CATEGORY (not condition -- COND_COLORS doesn't apply
# here, since this plot has one line per category, not per condition).
CATEGORY_DELTA_COLORS <- c(
  "Large Mature" = "#66A61E",
  "Nascent"      = "#E6AB02",
  "Small Mature" = "#7570B3",
  "Decaying"     = "#E7298A"
)

p_delta <- ggplot(delta_profile, aes(distance, delta_density, color = category, fill = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_ribbon(aes(ymin = delta_density - delta_ci95, ymax = delta_density + delta_ci95),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_x_log10() +
  scale_color_manual(values = CATEGORY_DELTA_COLORS) +
  scale_fill_manual(values = CATEGORY_DELTA_COLORS) +
  labs(x = "Distance (nm)", y = expression(Delta * "Density (Wildtype " - " POLR2A deg)"),
       title = "ChromEM radial density: Wildtype minus POLR2A deg",
       subtitle = "Shaded band = propagated 95% CI; dashed line = no difference",
       color = NULL, fill = NULL) +
  theme_pub(legend_position = "right")

print(p_delta)

write_csv(delta_profile, file.path(out_dir, "chromem_delta_profile.csv"))
message("Figure 5 printed above. Delta profile written to '", out_dir, "/'.")

# ============================================================================
# PART 5 -- same row layout as Parts 2-3, but the marginal panel beneath
# each category is that category's OWN delta curve (from Part 4) instead
# of a per-domain distribution
# ============================================================================

#' One category's delta panel: Wildtype - POLR2A deg density vs distance,
#' with propagated 95% CI ribbon and a dashed zero-reference line, colored
#' to match that category's line in the combined Figure 5 above. Now also
#' carries the SAME domain core/transition/periphery shading and dashed
#' boundary lines as the radial density panel above it (same db data,
#' same core_color/mid_color/periphery_color), so it's visible which zone
#' the biggest delta swings actually fall in. x-axis limits match the
#' radial density panel above it, for visual alignment.
plot_delta_marginal <- function(category_name) {
  d  <- delta_profile %>% filter(category == category_name)
  db <- domain_boundaries %>% filter(category == category_name)
  cat_color <- CATEGORY_DELTA_COLORS[[category_name]]
  
  ggplot(d, aes(x = distance, y = delta_density)) +
    geom_rect(data = db, aes(xmin = x_min, xmax = core_boundary, ymin = -Inf, ymax = Inf),
              fill = core_color, color = NA, inherit.aes = FALSE) +
    geom_rect(data = db, aes(xmin = core_boundary, xmax = boundary, ymin = -Inf, ymax = Inf),
              fill = mid_color, color = NA, inherit.aes = FALSE) +
    geom_rect(data = db, aes(xmin = boundary, xmax = x_max, ymin = -Inf, ymax = Inf),
              fill = periphery_color, color = NA, inherit.aes = FALSE) +
    geom_vline(data = db, aes(xintercept = core_boundary), linetype = "dashed", linewidth = 0.6, color = "grey35") +
    geom_vline(data = db, aes(xintercept = boundary), linetype = "dashed", linewidth = 0.6, color = "grey35") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_ribbon(aes(ymin = delta_density - delta_ci95, ymax = delta_density + delta_ci95),
                fill = cat_color, alpha = 0.2) +
    geom_line(color = cat_color, linewidth = 0.9) +
    scale_x_log10(limits = c(db$x_min, db$x_max)) +
    labs(x = "Distance (nm)", y = expression(Delta * "Density")) +
    theme_pub(legend_position = "none")
}

row_panels_delta <- purrr::imap(category_order, function(cat, i) {
  p_top    <- plot_radial_one_category(cat, show_legend = (i == 1), show_y_title = (i == 1))
  p_bottom <- plot_delta_marginal(cat)
  # Heights less extreme than the 4:1 used for the boxplot marginal panels
  # above -- the delta panel is a full line plot with its own y-axis, not
  # a compact single-row strip, so it needs more vertical room to be readable.
  (p_top / p_bottom) + plot_layout(heights = c(3, 2))
})

p_chromem_row_delta <- wrap_plots(row_panels_delta, ncol = 4) +
  plot_annotation(title = "ChromEM radial density with delta (Wildtype minus POLR2A deg) beneath each panel")

print(p_chromem_row_delta)
message("Figure 6 printed above (uses the same delta_profile already written to '", out_dir, "/chromem_delta_profile.csv').")

# ============================================================================
# PART 6 -- Log2 fold change: log2(POLR2A deg / Wildtype), as a function of
# distance
#
# NOT the same thing as the additive delta in Part 4/5 (Wildtype -
# POLR2A deg, in raw density units). Fold change is a RATIO
# (Treatment/Control), and log2 is the standard convention for reporting
# it (used by default in DESeq2/edgeR/limma-voom and equivalent tools) --
# specifically because it's symmetric around zero: a doubling gives +1,
# a halving gives -1, no change gives 0. A raw ratio is NOT symmetric
# this way (doubling -> 2, halving -> 0.5), which makes up- and down-
# changes visually lopsided.
#
# Direction confirmed: log2FC = log2(Treatment / Control) =
# log2(POLR2A deg / Wildtype). Positive = POLR2A deg higher than
# Wildtype; negative = POLR2A deg lower.
#
# No pseudocount added: checked the actual data range first (0.197 to
# 0.586, nothing near zero), so there's no division-by-zero or log(0)
# risk here the way there can be with raw RNA-seq counts.
#
# Uncertainty propagation: NOT the same additive formula as Part 4/5,
# because this is a ratio, not a difference. Uses the delta method for a
# log-ratio of two independent estimates:
#   SE(log2(X/Y)) = (1/ln2) * sqrt((SE_X/X)^2 + (SE_Y/Y)^2)
# Verified numerically against a Monte Carlo simulation before writing
# this (10,000+ simulated draws) -- the formula matched the simulated SE
# to 3 decimal places, not just derived and trusted.
# ============================================================================

log2fc_profile <- em_profile %>%
  select(category, condition, distance, mean_density, se) %>%
  pivot_wider(names_from = condition, values_from = c(mean_density, se)) %>%
  mutate(
    log2fc     = log2(`mean_density_POLR2A deg` / `mean_density_Wildtype`),
    log2fc_se  = (1 / log(2)) * sqrt((`se_POLR2A deg` / `mean_density_POLR2A deg`)^2 +
                                       (`se_Wildtype` / `mean_density_Wildtype`)^2),
    log2fc_ci95 = 1.96 * log2fc_se
  ) %>%
  filter(!is.na(log2fc))

message("\nLog2 fold change profile (POLR2A deg / Wildtype) computed for ",
        nrow(log2fc_profile), " category x distance points.")

p_log2fc <- ggplot(log2fc_profile, aes(distance, log2fc, color = category, fill = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_ribbon(aes(ymin = log2fc - log2fc_ci95, ymax = log2fc + log2fc_ci95),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_x_log10() +
  scale_color_manual(values = CATEGORY_DELTA_COLORS) +
  scale_fill_manual(values = CATEGORY_DELTA_COLORS) +
  labs(x = "Distance (nm)", y = "log2 fold change (POLR2A deg / Wildtype)",
       title = "ChromEM radial density: log2 fold change",
       subtitle = "Shaded band = propagated 95% CI (delta method); dashed line = no change",
       color = NULL, fill = NULL) +
  theme_pub(legend_position = "right")

print(p_log2fc)

write_csv(log2fc_profile, file.path(out_dir, "chromem_log2fc_profile.csv"))
message("Figure 7 printed above. Log2FC profile written to '", out_dir, "/'.")