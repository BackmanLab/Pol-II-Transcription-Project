##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# General Hi-C QC (interactive version). Designed to be run
# section-by-section in RStudio. No auto-saving. Plots render in
# the RStudio viewer / graphics device. Modify aesthetics, labels,
# and colour palettes inline before plotting.
#
# Sections:
#   0.  Libraries & paths
#   1.  Sample manifest
#   2.  Parser
#   3.  Parse all samples
#   4.  Colour palettes & theme  <- edit here for aesthetics
#   5a. Plot A -- alignment fate
#   5b. Plot B -- contact classification
#   6.  Summary table
#   7.  Save (manual, when ready)
#
# Figure output: SI Fig. 5A, C-F (probable).
##################################################################


# ── 0.  Libraries & paths ───────────────────────────────────────────────────
library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)

BASE   <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
OUTDIR <- file.path(BASE, "qc_plots")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# ── 1.  Sample manifest ─────────────────────────────────────────────────────
manifest <- tribble(
  ~dataset,      ~condition,    ~replicate, ~inter_txt,
  # ds_wt_pol2
  "ds_wt_pol2",  "WT",          "Rep1", file.path(BASE, "juicer_ds_wt_pol2_genomics/WT_HCT116_CTRL_Rep1/aligned/inter.txt"),
  "ds_wt_pol2",  "WT",          "Rep2", file.path(BASE, "juicer_ds_wt_pol2_genomics/WT_HCT116_CTRL_Rep2/aligned/inter.txt"),
  "ds_wt_pol2",  "Pol2 6hr",    "Rep1", file.path(BASE, "juicer_ds_wt_pol2_genomics/Pol2_6hrs_Aux_Rep1/aligned/inter.txt"),
  "ds_wt_pol2",  "Pol2 6hr",    "Rep2", file.path(BASE, "juicer_ds_wt_pol2_genomics/Pol2_6hrs_Aux_Rep2/aligned/inter.txt"),
  # ds_all
  "ds_all",      "WT",          "Rep1", file.path(BASE, "juicer_ds_all/WT_HCT116_CTRL_Rep1/aligned/inter.txt"),
  "ds_all",      "WT",          "Rep2", file.path(BASE, "juicer_ds_all/WT_HCT116_CTRL_Rep2/aligned/inter.txt"),
  "ds_all",      "Pol2 6hr",    "Rep1", file.path(BASE, "juicer_ds_all/Pol2_6hrs_Aux_Rep1/aligned/inter.txt"),
  "ds_all",      "Pol2 6hr",    "Rep2", file.path(BASE, "juicer_ds_all/Pol2_6hrs_Aux_Rep2/aligned/inter.txt"),
  "ds_all",      "ActD 1hr",    "Rep1", file.path(BASE, "juicer_ds_all/1hr_ActD_Rep1/aligned/inter.txt"),
  "ds_all",      "ActD 1hr",    "Rep2", file.path(BASE, "juicer_ds_all/1hr_ActD_Rep2/aligned/inter.txt"),
  # unsampled
  "unsampled",   "WT",          "Rep1", file.path(BASE, "juicer_unsampled/WT_HCT116_CTRL_Rep1/aligned/inter.txt"),
  "unsampled",   "WT",          "Rep2", file.path(BASE, "juicer_unsampled/WT_HCT116_CTRL_Rep2/aligned/inter.txt"),
  "unsampled",   "Pol2 6hr",    "Rep1", file.path(BASE, "juicer_unsampled/Pol2_6hrs_Aux_Rep1/aligned/inter.txt"),
  "unsampled",   "Pol2 6hr",    "Rep2", file.path(BASE, "juicer_unsampled/Pol2_6hrs_Aux_Rep2/aligned/inter.txt"),
  "unsampled",   "ActD 1hr",    "Rep1", file.path(BASE, "juicer_unsampled/1hr_ActD_Rep1/aligned/inter.txt"),
  "unsampled",   "ActD 1hr",    "Rep2", file.path(BASE, "juicer_unsampled/1hr_ActD_Rep2/aligned/inter.txt")
)

# Condition display order — edit to reorder bars
COND_ORDER <- c("WT", "Pol2 6hr", "ActD 1hr")


# ── 2.  Parser ───────────────────────────────────────────────────────────────
parse_inter_txt <- function(path) {
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  lines <- readLines(path)
  n <- function(pattern) {
    hit <- grep(pattern, lines, value = TRUE, perl = TRUE)
    if (!length(hit)) return(NA_real_)
    # Extract first number after the colon, not the first number on the line
    m <- regmatches(hit[1], regexpr("(?<=:\\s{0,10})[0-9,]+", hit[1], perl = TRUE))
    if (!length(m)) return(NA_real_)
    as.numeric(gsub(",", "", m))
  }
  tibble(
    sequenced      = n("^Sequenced Read Pairs:"),
    unmapped       = n("One or both reads unmapped:"),
    normal_pairs   = n("2 alignments \\(A\\.\\.\\.B\\):"),
    chimeric_pairs = n("2 alignments \\(A1"),
    multi_align    = n("^3 or more alignments:"),
    total_unique   = n("^Total Unique:"),
    duplicates     = n("^Total Duplicates:"),
    intra_frag     = n("^Intra-fragment Reads:"),
    below_mapq     = n("^Below MAPQ Threshold:"),
    inter_chrom    = n("^Inter-chromosomal:"),
    sr_500bp       = n("<500BP:"),
    sr_500_5kb     = n("500BP-5kB:"),
    sr_5_20kb      = n("5kB-20kB:"),
    long_range     = n("^Long Range"),
    lib_complexity = n("^Library Complexity Estimate"),
    ligation_motif = n("^Ligation Motif Present:")
  )
}


# ── 3.  Parse all samples ────────────────────────────────────────────────────
stats_raw <- manifest %>%
  filter(file.exists(inter_txt)) %>%
  mutate(stats = map(inter_txt, parse_inter_txt)) %>%
  unnest(stats) %>%
  mutate(
    condition = factor(condition, levels = COND_ORDER),
    sample    = paste0(condition, "\n", replicate)
  )

# Inspect raw numbers
stats_raw %>%
  select(dataset, condition, replicate, sequenced, total_unique,
         duplicates, long_range, lib_complexity) %>%
  print(n = Inf)


# ── 4.  Palettes & theme  ←  EDIT THIS SECTION ──────────────────────────────

# Plot A: alignment fate — bottom to top stacking order
pal_align <- c(
  "Unmapped"       = "#d9d9d9",
  "Normal pairs"   = "#6baed6",
  "Chimeric pairs" = "#2171b5",
  "3+ alignments"  = "#08306b"
)

# Plot B: contact classification — bottom to top stacking order
pal_contacts <- c(
  "Duplicates"         = "#d9d9d9",
  "Intra-fragment"     = "#fdae6b",
  "Below MAPQ"         = "#fdd0a2",
  "Inter-chrom"        = "#74c476",
  "Short (<500bp)"     = "#fecc5c",
  "Short (500bp-5kb)"  = "#fd8d3c",
  "Short (5kb-20kb)"   = "#e31a1c",
  "Long range (>20kb)" = "#1a1a2e"
)

# Base theme — adjust fonts, grid, etc.
base_theme <- theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    axis.title       = element_text(size = 11),
    legend.position  = "right",
    legend.key.size  = unit(0.45, "cm"),
    legend.text      = element_text(size = 9),
    legend.title     = element_blank(),
    strip.background = element_rect(fill = "#f5f5f5", colour = "grey80"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.spacing    = unit(0.6, "cm"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40")
  )

BAR_WIDTH <- 0.72   # bar width (0–1)


# ── 5a. Plot A — alignment fate ──────────────────────────────────────────────
# Change DATASET_A to switch between "ds_wt_pol2", "ds_all", "unsampled"

DATASET_A <- "ds_wt_pol2"
DATASET_A <- "ds_all"
DATASET_A <- "unsampled"

align_long <- stats_raw %>%
  filter(dataset == DATASET_A) %>%
  transmute(
    condition, replicate,
    "Unmapped"       = unmapped       / sequenced * 100,
    "Normal pairs"   = normal_pairs   / sequenced * 100,
    "Chimeric pairs" = chimeric_pairs / sequenced * 100,
    "3+ alignments"  = multi_align    / sequenced * 100
  ) %>%
  pivot_longer(names(pal_align), names_to = "category", values_to = "pct") %>%
  mutate(category = factor(category, levels = names(pal_align)))

p_align <- ggplot(align_long,
                  aes(x = replicate, y = pct, fill = category)) +
  geom_col(width = BAR_WIDTH, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(pct >= 3, sprintf("%.1f%%", pct), "")),
    position = position_stack(vjust = 0.5),
    colour   = "white",
    fontface = "bold",
    size     = 3.2
  ) +
  scale_fill_manual(values = pal_align) +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    limits = c(0, 105),
    expand = expansion(mult = c(0, 0))
  ) +
  facet_wrap(~ condition, nrow = 1, scales = "free_x") +
  labs(
    title    = paste0("Alignment fate — ", DATASET_A),
    subtitle = "% of sequenced read pairs",
    x        = NULL,
    y        = "% of sequenced read pairs"
  ) +
  base_theme

p_align


# ── 5b. Plot B — contact classification ─────────────────────────────────────
# Change DATASET_B to switch groups

DATASET_B <- "ds_wt_pol2"
DATASET_B <- "ds_all"
DATASET_B <- "unsampled"

contact_long <- stats_raw %>%
  filter(dataset == DATASET_B) %>%
  transmute(
    condition, replicate,
    total = total_unique,
    "Duplicates"         = duplicates  / total * 100,
    "Intra-fragment"     = intra_frag  / total * 100,
    "Below MAPQ"         = below_mapq  / total * 100,
    "Inter-chrom"        = inter_chrom / total * 100,
    "Short (<500bp)"     = sr_500bp    / total * 100,
    "Short (500bp-5kb)"  = sr_500_5kb  / total * 100,
    "Short (5kb-20kb)"   = sr_5_20kb   / total * 100,
    "Long range (>20kb)" = long_range  / total * 100
  ) %>%
  pivot_longer(names(pal_contacts), names_to = "category", values_to = "pct") %>%
  mutate(category = factor(category, levels = names(pal_contacts)))

p_contacts <- ggplot(contact_long,
                     aes(x = replicate, y = pct, fill = category)) +
  geom_col(width = BAR_WIDTH, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(pct >= 3, sprintf("%.1f%%", pct), "")),
    position = position_stack(vjust = 0.5),
    colour   = "white",
    fontface = "bold",
    size     = 3.2
  ) +
  scale_fill_manual(values = pal_contacts) +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    limits = c(0, 105),
    expand = expansion(mult = c(0, 0))
  ) +
  facet_wrap(~ condition, nrow = 1, scales = "free_x") +
  labs(
    title    = paste0("Contact classification — ", DATASET_B),
    subtitle = "% of total unique mapped pairs",
    x        = NULL,
    y        = "% of unique mapped pairs"
  ) +
  base_theme

p_contacts


# ── Optional: combined layout ─────────────────────────────────────────────────
p_combined <- p_align / p_contacts +
  plot_annotation(
    title = paste0("Hi-C QC — ", DATASET_A),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

p_combined


# ── 6.  Summary table ────────────────────────────────────────────────────────
summary_tbl <- stats_raw %>%
  mutate(
    pct_hic       = (inter_chrom + sr_500bp + sr_500_5kb +
                       sr_5_20kb + long_range) / sequenced * 100,
    pct_longrange = long_range / sequenced * 100,
    pct_dups      = duplicates / (total_unique + duplicates) * 100,
    lib_M         = lib_complexity / 1e6
  ) %>%
  select(dataset, condition, replicate, sequenced,
         total_unique, pct_dups, pct_hic, pct_longrange, lib_M) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2)))

print(summary_tbl, n = Inf)

# write_tsv(summary_tbl, file.path(OUTDIR, "hic_qc_summary.tsv"))


# ── 5c. Plot C — QC metrics panel ────────────────────────────────────────────
# Four grouped bar charts: library complexity, % duplicates,
# % long-range contacts, % ligation motif
# Each plot: x = condition, dodge = replicate, facet = dataset
# Change DATASETS_C to subset which groups to show

DATASETS_C <- c("ds_wt_pol2", "ds_all", "unsampled")  # ← subset as needed

# Colour for replicates
pal_rep <- c("Rep1" = "#2171b5", "Rep2" = "#6baed6")

# Shared theme for metric plots — slightly more compact
metric_theme <- base_theme +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 10),
    panel.spacing = unit(0.4, "cm")
  )

# Build tidy metrics dataframe
metrics_long <- stats_raw %>%
  filter(dataset %in% DATASETS_C) %>%
  mutate(
    lib_complexity_M  = lib_complexity / 1e6,
    pct_ligation      = ligation_motif / sequenced * 100,
    pct_dups          = duplicates / (total_unique + duplicates) * 100,
    pct_longrange     = long_range / sequenced * 100
  ) %>%
  select(dataset, condition, replicate,
         lib_complexity_M, pct_dups, pct_longrange, pct_ligation) %>%
  pivot_longer(
    cols      = c(lib_complexity_M, pct_dups, pct_longrange, pct_ligation),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    lib_complexity_M = "Library complexity (M)",
                    pct_dups         = "% Duplicates",
                    pct_longrange    = "% Long-range (>20kb)",
                    pct_ligation     = "% Ligation motif"
    ),
    metric = factor(metric, levels = c(
      "Library complexity (M)",
      "% Duplicates",
      "% Long-range (>20kb)",
      "% Ligation motif"
    )),
    dataset = factor(dataset, levels = DATASETS_C)
  )

make_metric_plot <- function(data, metric_name, y_label, y_limits = NULL) {
  d <- data %>% filter(metric == metric_name)
  p <- ggplot(d, aes(x = condition, y = value, fill = replicate)) +
    geom_col(
      position = position_dodge(width = 0.75),
      width    = 0.65,
      colour   = "white",
      linewidth = 0.2
    ) +
    geom_text(
      aes(label = round(value, 1)),
      position = position_dodge(width = 0.75),
      hjust    = 1.2,     # pushes text just inside the top of the bar
      angle    = 90,
      size     = 2.8,
      colour   = "white"
    ) +
    scale_fill_manual(values = pal_rep, name = NULL) +
    facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
    labs(title = metric_name, x = NULL, y = y_label) +
    metric_theme
  
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(
      limits = y_limits,
      expand = expansion(mult = c(0, 0.1))
    )
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  }
  p
}

p_lib   <- make_metric_plot(metrics_long, "Library complexity (M)",  "Unique molecules (millions)")
p_dups  <- make_metric_plot(metrics_long, "% Duplicates",            "% of mappable pairs")
p_lr    <- make_metric_plot(metrics_long, "% Long-range (>20kb)",    "% of sequenced pairs")
p_motif <- make_metric_plot(metrics_long, "% Ligation motif",        "% of sequenced pairs")

# View individually
p_lib
p_dups
p_lr
p_motif

# Or combined 2×2 grid
p_metrics <- (p_lib | p_dups) / (p_lr | p_motif) +
  plot_annotation(
    title = "Hi-C QC metrics",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

p_metrics

# ── 5d. Plot D — Absolute contact pairs (log10) ───────────────────────────
# Intra and inter-chromosomal contact counts per sample, log10 scale.
# Compares depth across unsampled and downsampled groups.

DATASETS_D <- c("ds_wt_pol2", "ds_all", "unsampled")

# Option 2: Deep red / teal
pal_contact_type <- c("Intra-chromosomal" = "#C0392B",   # deep red
                      "Inter-chromosomal" = "#16A085")   # teal


# Colour for replicates
pal_rep <- c("Rep1" = "#C0392B", "Rep2" = "#16A085")

contacts_abs <- stats_raw %>%
  filter(dataset %in% DATASETS_D) %>%
  mutate(
    intra = sr_500bp + sr_500_5kb + sr_5_20kb + long_range,
    inter = inter_chrom
  ) %>%
  select(dataset, condition, replicate, intra, inter) %>%
  pivot_longer(c(intra, inter),
               names_to  = "contact_type",
               values_to = "count") %>%
  mutate(
    contact_type = recode(contact_type,
                          intra = "Intra-chromosomal",
                          inter = "Inter-chromosomal"),
    contact_type = factor(contact_type,
                          levels = c("Intra-chromosomal", "Inter-chromosomal")),
    dataset = factor(dataset, levels = DATASETS_D)
  )

p_contacts_abs <- ggplot(contacts_abs,
                         aes(x = replicate, y = count, fill = contact_type,
                             alpha = contact_type)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width    = 0.65,
    colour   = "white",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = scales::comma(round(count / 1e6, 1), suffix = "M")),
    position = position_dodge(width = 0.75),
    hjust    = 1.15,
    angle    = 90,
    size     = 2.5,
    colour   = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(values = pal_contact_type, name = NULL) +
  scale_alpha_manual(values = c("Intra-chromosomal" = 1,
                                "Inter-chromosomal" = 0.75),
                     guide = "none") +
  scale_y_log10(
    labels = scales::label_log(digits = 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_grid(condition ~ dataset, scales = "free_x") +
  labs(
    title    = "Hi-C contact pairs per sample",
    subtitle = "Intra- and inter-chromosomal contacts (log₁₀ scale)",
    x        = NULL,
    y        = "Contact pairs (log₁₀)"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

p_contacts_abs

# ── 5e. Plot E — Total contacts as % of WT ───────────────────────────────
# Total Hi-C contacts (intra + inter) per sample expressed as a percentage
# of the mean WT contact count within each dataset group.
# Shows how contact depth in Pol2/ActD compares to WT after downsampling.

DATASETS_E <- c("ds_wt_pol2", "ds_all", "unsampled")

contacts_pct_wt <- stats_raw %>%
  filter(dataset %in% DATASETS_E) %>%
  mutate(
    total_contacts = sr_500bp + sr_500_5kb + sr_5_20kb + long_range + inter_chrom,
    dataset = factor(dataset, levels = DATASETS_E)
  ) %>%
  group_by(replicate) %>%
  mutate(
    wt_unsampled = total_contacts[dataset == "unsampled" & condition == "WT"],
    pct_of_wt    = total_contacts / wt_unsampled * 100
  ) %>%
  ungroup()

p_pct_wt <- ggplot(contacts_pct_wt,
                   aes(x = condition, y = pct_of_wt, fill = replicate)) +
  geom_hline(yintercept = 100, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_col(
    position = position_dodge(width = 0.75),
    width    = 0.65,
    colour   = "white",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_of_wt)),
    position = position_dodge(width = 0.75),
    hjust    = 1.15,
    angle    = 90,
    size     = 2.8,
    colour   = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(values = pal_rep, name = NULL) +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
  labs(
    title    = "Total Hi-C contacts relative to WT",
    subtitle = "% of mean WT contact count within each dataset group",
    x        = NULL,
    y        = "% of WT contacts"
  ) +
  base_theme

p_pct_wt


