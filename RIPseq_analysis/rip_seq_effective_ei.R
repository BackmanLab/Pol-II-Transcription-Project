##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# RIP-seq Pol2 (P2AB) IP enrichment analysis — exon and intron
# RPKM, log2(P2AB / Input) computed two ways:
#   log2(mean P2AB / WT_Input)
#   log2(mean P2AB / P2KO_Input)
#
# P2AB replicates (Rep2, Rep3) are AVERAGED AT THE RPKM LEVEL
# before the ratio is taken (chosen over the alternative of
# ratio-per-replicate-then-average-log2FC).
#
# WT_Input and P2KO_Input are both single, unreplicated samples —
# no replicate handling needed on the denominator side.
#
# EZH2 RPKM is present in the source matrices for reference/QC but
# not used in either ratio here.
##################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(BiocGenerics)
})

##--------------------------------------------------------------## 1. Parameters

RIP_DIR <- "/projects/b1042/BackmanLab/Lucas/022025_RIPseq/ripseq_pipeline"
EXON_RPKM_FILE   <- file.path(RIP_DIR, "exon_counts",   "ripseq_exon_rpkm.txt")
INTRON_RPKM_FILE <- file.path(RIP_DIR, "intron_counts", "ripseq_intron_rpkm.txt")

PSEUDOCOUNT    <- 1     # standard RNA-seq convention: log2(RPKM + 1)
RPKM_THRESHOLD <- 0.5   # minimum mean P2AB RPKM to trust a gene's ratio —
# same threshold value used in the EU-seq E/I script
# (used only in the diagnostics section below — the
# effective E/I section further down uses a
# different, exon-specific filtering rule instead)
GTF_PATH       <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
## NOTE: a fixed EI_PSEUDOCOUNT was previously defined here (1e-3). It's
## been replaced by an empirically-computed INTRON_PSEUDOCOUNT in the
## effective E/I section below (the minimum non-zero intron RPKM actually
## observed in this dataset, pooled across P2AB and WT_Input) — see the
## reasoning there. No fixed exon-side pseudocount is used at all: exon
## RPKM==0 rows are removed rather than pseudocounted.

## Kept for palette continuity with the EU-seq scripts, even though this
## script's plots key off Input identity (WT_Input/P2KO_Input) rather
## than a WT/POLR2A condition pair.
INPUT_COLORS <- c(WT_Input = "#88CCEE", P2KO_Input = "#CC6677")

base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    strip.text       = element_text(face = "bold", size = 11)
  )

##--------------------------------------------------------------## 2. Load RPKM matrices

message("Loading RIP-seq RPKM matrices...")

exon_rpkm <- read_tsv(EXON_RPKM_FILE, show_col_types = FALSE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id))
intron_rpkm <- read_tsv(INTRON_RPKM_FILE, show_col_types = FALSE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id))

cat("Exon RPKM matrix:   ", nrow(exon_rpkm),   "genes x", ncol(exon_rpkm) - 1,   "samples\n")
cat("Intron RPKM matrix: ", nrow(intron_rpkm), "genes x", ncol(intron_rpkm) - 1, "samples\n")
cat("Exon RPKM columns:   ", paste(colnames(exon_rpkm)[-1],   collapse = ", "), "\n")
cat("Intron RPKM columns: ", paste(colnames(intron_rpkm)[-1], collapse = ", "), "\n")

##--------------------------------------------------------------## 3. Average P2AB replicates at the RPKM level, compute log2 ratios
## log2fc = log2((mean_P2AB + PSEUDOCOUNT) / (Input + PSEUDOCOUNT))
## Done independently for exon and intron signal, and independently
## against each of the two Input denominators.

compute_ratios <- function(rpkm_df, interval_label) {
  rpkm_df %>%
    dplyr::mutate(
      mean_P2AB = rowMeans(dplyr::select(., P2AB_Rep2, P2AB_Rep3), na.rm = TRUE),
      log2fc_wt_input   = log2((mean_P2AB + PSEUDOCOUNT) / (WT_Input   + PSEUDOCOUNT)),
      log2fc_p2ko_input = log2((mean_P2AB + PSEUDOCOUNT) / (P2KO_Input + PSEUDOCOUNT))
    ) %>%
    dplyr::select(gene_id, mean_P2AB, WT_Input, P2KO_Input,
                  log2fc_wt_input, log2fc_p2ko_input) %>%
    dplyr::rename_with(~ paste0(interval_label, "_", .x), -gene_id)
}

exon_ratios   <- compute_ratios(exon_rpkm,   "exon")
intron_ratios <- compute_ratios(intron_rpkm, "intron")

message("\nExon mean_P2AB RPKM summary:")
print(summary(exon_ratios$exon_mean_P2AB))
message("Intron mean_P2AB RPKM summary:")
print(summary(intron_ratios$intron_mean_P2AB))

##--------------------------------------------------------------## 4. Merge exon + intron into one table, apply detection filter

ripseq_ratios <- exon_ratios %>%
  dplyr::inner_join(intron_ratios, by = "gene_id")

cat("\nGenes with both exon and intron RPKM computed:", nrow(ripseq_ratios), "\n")

## Filter to genes with detectable P2AB signal. OR logic (exon OR intron
## clears the bar) rather than AND — a gene with strong exon-only or
## intron-only Pol2 binding still passes. Switch to AND for a stricter bar.
ripseq_ratios_filtered <- ripseq_ratios %>%
  dplyr::filter(exon_mean_P2AB >= RPKM_THRESHOLD | intron_mean_P2AB >= RPKM_THRESHOLD)

cat("Genes passing P2AB RPKM >=", RPKM_THRESHOLD, "in exon OR intron:",
    nrow(ripseq_ratios_filtered), "\n")

##--------------------------------------------------------------## 5. Diagnostic summaries

message("\n== log2(P2AB / WT_Input) — exon ==")
print(summary(ripseq_ratios_filtered$exon_log2fc_wt_input))
message("== log2(P2AB / WT_Input) — intron ==")
print(summary(ripseq_ratios_filtered$intron_log2fc_wt_input))
message("== log2(P2AB / P2KO_Input) — exon ==")
print(summary(ripseq_ratios_filtered$exon_log2fc_p2ko_input))
message("== log2(P2AB / P2KO_Input) — intron ==")
print(summary(ripseq_ratios_filtered$intron_log2fc_p2ko_input))

## Correlation between the two Input denominators — if WT_Input and
## P2KO_Input give very different enrichment calls for the same genes,
## that's worth knowing before trusting either one alone.
cor_exon <- cor(ripseq_ratios_filtered$exon_log2fc_wt_input,
                ripseq_ratios_filtered$exon_log2fc_p2ko_input,
                use = "complete.obs")
cor_intron <- cor(ripseq_ratios_filtered$intron_log2fc_wt_input,
                  ripseq_ratios_filtered$intron_log2fc_p2ko_input,
                  use = "complete.obs")
cat("\nCorrelation between WT_Input-based and P2KO_Input-based log2FC:\n")
cat("  Exon:   r =", round(cor_exon,   3), "\n")
cat("  Intron: r =", round(cor_intron, 3), "\n")

##--------------------------------------------------------------## 6. Diagnostic plots

message("\nPlotting log2FC distributions...")

## Long format across all 4 log2FC distributions (exon/intron x WT/P2KO Input)
plot_df <- dplyr::bind_rows(
  ripseq_ratios_filtered %>% dplyr::transmute(gene_id, interval = "Exon",   input = "WT_Input",   log2fc = exon_log2fc_wt_input),
  ripseq_ratios_filtered %>% dplyr::transmute(gene_id, interval = "Exon",   input = "P2KO_Input", log2fc = exon_log2fc_p2ko_input),
  ripseq_ratios_filtered %>% dplyr::transmute(gene_id, interval = "Intron", input = "WT_Input",   log2fc = intron_log2fc_wt_input),
  ripseq_ratios_filtered %>% dplyr::transmute(gene_id, interval = "Intron", input = "P2KO_Input", log2fc = intron_log2fc_p2ko_input)
) %>%
  dplyr::mutate(
    interval = factor(interval, levels = c("Exon","Intron")),
    input    = factor(input,    levels = c("WT_Input","P2KO_Input"))
  )

p_log2fc <- ggplot(plot_df, aes(x = log2fc, fill = input, color = input)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = INPUT_COLORS, name = NULL) +
  scale_color_manual(values = INPUT_COLORS, name = NULL) +
  facet_wrap(~ interval, nrow = 1, scales = "free_y") +
  labs(
    x        = "log2(P2AB / Input)",
    y        = "Density",
    title    = "Pol2 RIP-seq enrichment over Input",
    subtitle = paste0("Mean of P2AB_Rep2/Rep3 RPKM  |  P2AB RPKM \u2265 ", RPKM_THRESHOLD,
                      " (exon or intron)  |  n = ", nrow(ripseq_ratios_filtered), " genes")
  ) +
  base_theme

print(p_log2fc)

## Exon vs intron log2FC scatter (WT_Input-based) — does enrichment
## concentrate in exonic or intronic signal for a given gene?
p_exon_v_intron <- ripseq_ratios_filtered %>%
  dplyr::select(gene_id, exon_log2fc_wt_input, intron_log2fc_wt_input) %>%
  ggplot(aes(x = exon_log2fc_wt_input, y = intron_log2fc_wt_input)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#44AA99") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  labs(
    x        = "log2(P2AB / WT_Input) \u2014 exon",
    y        = "log2(P2AB / WT_Input) \u2014 intron",
    title    = "Exonic vs intronic Pol2 enrichment (vs WT_Input)",
    subtitle = paste0("n = ", nrow(ripseq_ratios_filtered), " genes  |  dashed line = y = x")
  ) +
  base_theme +
  theme(legend.position = "none")

print(p_exon_v_intron)

message("\nDiagnostics complete — continuing to effective E/I analysis below.")

##--------------------------------------------------------------## 7. Genomic E/I reference + intron length categories (from GTF)
## Same computation as the EU-seq effective E/I script, reused here
## so RIP-seq E/I values are compared against the identical genomic
## reference and length categorization. Intronless genes are KEPT
## (not filtered out) — their ei_genomic reduces to roughly their
## exon length in bp, and they're included deliberately per an
## explicit decision, not an oversight.

message("\nComputing genomic E/I reference from GTF...")

txdb          <- GenomicFeatures::makeTxDbFromGFF(GTF_PATH, format = "gtf", taxonomyId = 9606)
exons_by_gene <- GenomicRanges::reduce(GenomicFeatures::exonsBy(txdb, by = "gene"))
exon_lengths  <- sum(BiocGenerics::width(exons_by_gene))
genes_gr      <- GenomicFeatures::genes(txdb)
gene_lengths  <- setNames(BiocGenerics::width(genes_gr), genes_gr$gene_id)

gtf_common <- intersect(names(exon_lengths), names(gene_lengths))
intron_len <- pmax(gene_lengths[gtf_common] - exon_lengths[gtf_common], 0)
ei_genomic <- exon_lengths[gtf_common] / (intron_len + 1)

INTRON_CATS <- c("Short (<1kb)","Average (1-10kb)","Long (>10kb)")

genomic_df <- data.frame(
  gene_id       = sub("\\..*", "", gtf_common),
  intron_len_bp = as.numeric(intron_len),
  ei_genomic    = as.numeric(ei_genomic),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(ei_genomic), ei_genomic > 0, intron_len_bp >= 0) %>%   ## >=0 keeps intronless genes
  dplyr::mutate(
    intron_cat = dplyr::case_when(
      intron_len_bp <  1000  ~ "Short (<1kb)",
      intron_len_bp <= 10000 ~ "Average (1-10kb)",
      TRUE                   ~ "Long (>10kb)"
    ),
    intron_cat = factor(intron_cat, levels = INTRON_CATS)
  )

n_intronless <- sum(genomic_df$intron_len_bp == 0)
cat("Genomic E/I reference computed for:", nrow(genomic_df), "genes\n")
cat("  (of which", n_intronless, "are intronless — kept in, not excluded)\n")
print(table(genomic_df$intron_cat))

##--------------------------------------------------------------## 8. Empirical intron pseudocount
## The minimum non-zero intron RPKM actually observed in this
## dataset, pooled across P2AB and WT_Input, so "zero intron signal"
## sits at the same y-position in both distributions in the plot
## below (a difference in pseudocount alone shouldn't masquerade as
## a biological difference between P2AB and Input).

intron_values_pooled <- c(ripseq_ratios$intron_mean_P2AB, ripseq_ratios$intron_WT_Input)
intron_values_pooled <- intron_values_pooled[!is.na(intron_values_pooled) & intron_values_pooled > 0]
INTRON_PSEUDOCOUNT <- min(intron_values_pooled)

cat("\nIntron pseudocount (empirical minimum non-zero intron RPKM,",
    "pooled P2AB + WT_Input):", INTRON_PSEUDOCOUNT, "\n")
cat("  (", length(intron_values_pooled), "non-zero intron RPKM values contributed to this minimum )\n")

##--------------------------------------------------------------## 9. Effective E/I — anchored on the WT_Input baseline
## ONE gene set, defined by Input (baseline) exon signal:
## exon_WT_Input > 0. That same gene set is used to compute BOTH
## ei_WT_Input and ei_P2AB — a gene qualifies via its Input baseline,
## since Input defines what's actually being transcribed/present,
## which is the reference the IP is being asked about.
##
## If P2AB exon RPKM is exactly 0 for an Input-qualifying gene, the
## gene is KEPT (not dropped) and ei_P2AB is computed as exactly 0 —
## the IP failing to pull down a gene that Input shows is present is
## itself a real result, not something to discard. These rows can't
## appear on the log10 plots below (log10(0) is undefined) and get
## excluded there via the existing ei_ratio>0 filter, but they're
## counted explicitly here first so that exclusion is never silent.
##
## Intron RPKM is NEVER removed for being zero — "no observed intron
## signal" is the finding this whole analysis exists to detect. It
## always gets INTRON_PSEUDOCOUNT before the ratio is taken.

ei_df <- ripseq_ratios %>%
  dplyr::filter(exon_WT_Input > 0) %>%
  dplyr::mutate(
    ei_WT_Input = exon_WT_Input  / (intron_WT_Input  + INTRON_PSEUDOCOUNT),
    ei_P2AB     = exon_mean_P2AB / (intron_mean_P2AB + INTRON_PSEUDOCOUNT)
  ) %>%
  dplyr::inner_join(genomic_df %>% dplyr::select(gene_id, intron_cat), by = "gene_id") %>%
  dplyr::filter(!is.na(intron_cat))

cat("\nWT_Input-anchored gene set (exon WT_Input RPKM > 0):", nrow(ei_df), "genes\n")

## CHECK: how many Input-anchored genes have ZERO P2AB exon signal?
n_p2ab_exon_zero <- sum(ei_df$exon_mean_P2AB == 0)
cat("Of these, genes with ZERO P2AB exon signal (kept in ei_df with",
    "ei_P2AB=0 exactly; excluded from log10 plots below, not from the data):",
    n_p2ab_exon_zero, " (", round(100 * n_p2ab_exon_zero / nrow(ei_df), 1), "% )\n")

##--------------------------------------------------------------## 10. Effective E/I plot — Genomic / Input / P2AB, by intron length
## Built as a reusable function with an include_input SWITCH: TRUE
## produces the original 3-line comparison (Genomic/Input/P2AB),
## FALSE drops the Input line entirely (Genomic/P2AB only). Same
## panel-drawing code (build_ripseq_ei_panel_generic) is reused for
## both — only the data fed into it and the active source levels
## change based on the switch.
##
## Also runs with two different genomic reference definitions:
##   (a) Genomic reference = full GTF annotation, unrestricted
##   (b) Genomic reference = restricted to genes expressed in Input
##       (now exactly ei_df's own gene set — CONFIRMED as the
##       approach to build forward from)

message("\nBuilding RIP-seq effective E/I comparison plots, by intron length...")

ei_source_colors <- c(
  "Genomic\n(reference)" = "grey50",
  "Input\n(WT)"          = "#88CCEE",
  "P2AB\n(Pol2 IP)"      = "#CC6677"
)
ei_source_linetypes <- c(
  "Genomic\n(reference)" = "dashed",
  "Input\n(WT)"          = "solid",
  "P2AB\n(Pol2 IP)"      = "solid"
)
ei_source_levels <- c("Genomic\n(reference)", "Input\n(WT)", "P2AB\n(Pol2 IP)")

## Standalone, reusable density+outline-box panel builder — same pattern
## as the EU-seq script's build_ei_panel_generic. active_levels controls
## which source(s) actually appear (legend, y-axis order); ei_source_colors
## /ei_source_linetypes can safely carry unused entries (scale_*_manual
## matches by name, extras are simply ignored — verified, not an error).
build_ripseq_ei_panel_generic <- function(sub_df, panel_title, x_range, active_levels) {
  p_density <- ggplot(sub_df, aes(x = log10_ei, color = source, linetype = source)) +
    geom_density(linewidth = 0.9, adjust = 1.2) +
    geom_vline(xintercept = 0, color = "grey70", linetype = "dotted", linewidth = 0.4) +
    scale_color_manual(values = ei_source_colors, name = NULL) +
    scale_linetype_manual(values = ei_source_linetypes, name = NULL) +
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
    scale_color_manual(values = ei_source_colors) +
    scale_x_continuous(limits = x_range) +
    scale_y_discrete(limits = rev(active_levels)) +
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

## active_source_levels(): one-line helper so every plot function below
## computes the same thing the same way instead of repeating the
## setdiff() call — the actual switch logic lives in exactly one place.
active_source_levels <- function(include_input) {
  if (include_input) ei_source_levels else setdiff(ei_source_levels, "Input\n(WT)")
}

build_ripseq_ei_comparison <- function(genomic_ref_df, genomic_label, title_suffix, include_input = TRUE) {
  
  active_levels <- active_source_levels(include_input)
  
  plot_df <- dplyr::bind_rows(
    genomic_ref_df %>%
      dplyr::transmute(gene_id, intron_cat, ei_ratio = ei_genomic, source = "Genomic\n(reference)"),
    if (include_input) {
      ei_df %>% dplyr::transmute(gene_id, intron_cat, ei_ratio = ei_WT_Input, source = "Input\n(WT)")
    } else NULL,
    ei_df %>%
      dplyr::transmute(gene_id, intron_cat, ei_ratio = ei_P2AB, source = "P2AB\n(Pol2 IP)")
  ) %>%
    dplyr::filter(ei_ratio > 0) %>%
    dplyr::mutate(
      log10_ei = log10(ei_ratio),
      source   = factor(source, levels = active_levels)
    )
  
  panels <- lapply(levels(plot_df$intron_cat), function(cat_name) {
    sub_df  <- plot_df %>% dplyr::filter(intron_cat == cat_name)
    x_range <- range(sub_df$log10_ei, na.rm = TRUE)
    build_ripseq_ei_panel_generic(sub_df, cat_name, x_range, active_levels)
  })
  
  final_plot <- (patchwork::wrap_plots(panels, nrow = 1) + patchwork::plot_layout(guides = "collect")) +
    patchwork::plot_annotation(
      title    = paste0("RIP-seq effective E/I: Genomic vs ",
                        if (include_input) "Input vs " else "",
                        "Pol2 IP (P2AB), by intron length \u2014 ", title_suffix),
      subtitle = paste0(
        "Grey dashed = genomic E/I (", genomic_label, ")  |  ",
        "Solid = RIP-seq exon RPKM / intron RPKM  |  Input-anchored, P2AB exon=0 kept  |  ",
        "Genomic n = ", dplyr::n_distinct(genomic_ref_df$gene_id),
        "  |  WT_Input-anchored n = ", nrow(ei_df)
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(color = "grey40", size = 11))
    ) &
    theme(legend.position = "top")
  
  print(final_plot)
  invisible(list(plot_df = plot_df, final_plot = final_plot))
}

## (a) Genomic reference = full GTF annotation, unrestricted (3-line)
message("\n\u2500\u2500 Version A: genomic reference = all annotated genes \u2500\u2500")
result_all_genes <- build_ripseq_ei_comparison(
  genomic_ref_df = genomic_df,
  genomic_label  = "all annotated genes",
  title_suffix   = "genomic ref = all genes"
)

## (b) Genomic reference = restricted to genes expressed in Input —
## now exactly ei_df's own gene set, since ei_df is WT_Input-anchored.
## CONFIRMED as the approach to build forward from.
genomic_df_input_expressed <- genomic_df %>%
  dplyr::filter(gene_id %in% ei_df$gene_id)

cat("\nGenomic reference (Input-expressed only):",
    nrow(genomic_df_input_expressed), "genes",
    "(vs", nrow(genomic_df), "in the unrestricted version)\n")

message("\n\u2500\u2500 Version B: genomic reference = genes expressed in Input only (3-line) \u2500\u2500")
result_input_expressed <- build_ripseq_ei_comparison(
  genomic_ref_df = genomic_df_input_expressed,
  genomic_label  = "Input-expressed genes only",
  title_suffix   = "genomic ref = Input-expressed genes"
)

## NEW: same Version B genomic reference, but Input line excluded —
## just Genomic reference + P2AB
message("\n\u2500\u2500 Version B: genomic reference = genes expressed in Input only (P2AB only, no Input line) \u2500\u2500")
result_input_expressed_p2ab_only <- build_ripseq_ei_comparison(
  genomic_ref_df = genomic_df_input_expressed,
  genomic_label  = "Input-expressed genes only",
  title_suffix   = "genomic ref = Input-expressed genes, P2AB only",
  include_input  = FALSE
)

message("\nSection 10 complete. Key objects so far: ei_df, genomic_df,",
        " genomic_df_input_expressed, result_all_genes, result_input_expressed,",
        " result_input_expressed_p2ab_only")

##--------------------------------------------------------------## 11. Pseudocount sensitivity — how much does pseudocount choice
## inflate E/I for high-exon, zero-intron genes?
## Sweeps a range of candidate intron pseudocounts and recomputes E/I
## for the genes most exposed to pseudocount-driven inflation: real
## P2AB exon signal, but zero observed P2AB intron signal. If these
## genes' ratios move a lot across the sweep, the pseudocount choice
## is doing a lot of the work in the final number, not the biology.

message("\nBuilding pseudocount sensitivity analysis...")

## Log-spaced sweep from well below the smallest realistic value up to
## a value comparable to typical intron RPKM, plus the actual value used
## in the main analysis for direct reference.
pseudocount_sweep <- sort(unique(c(10^seq(-4, 1, length.out = 25), INTRON_PSEUDOCOUNT)))

n_zero_intron_p2ab <- sum(ripseq_ratios$exon_mean_P2AB > 0 & ripseq_ratios$intron_mean_P2AB == 0)
cat("\nGenes with P2AB exon signal but ZERO P2AB intron signal:", n_zero_intron_p2ab, "\n")

inflation_candidates <- ripseq_ratios %>%
  dplyr::filter(exon_mean_P2AB > 0, intron_mean_P2AB == 0) %>%
  dplyr::arrange(dplyr::desc(exon_mean_P2AB)) %>%
  dplyr::slice_head(n = 20)

cat("Top 20 selected for sensitivity plot, P2AB exon RPKM range:",
    round(min(inflation_candidates$exon_mean_P2AB), 2), "-",
    round(max(inflation_candidates$exon_mean_P2AB), 2), "\n")

sensitivity_df <- tidyr::expand_grid(
  gene_id     = inflation_candidates$gene_id,
  pseudocount = pseudocount_sweep
) %>%
  dplyr::left_join(
    inflation_candidates %>% dplyr::select(gene_id, exon_mean_P2AB, intron_mean_P2AB),
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    ei_at_pseudocount = exon_mean_P2AB / (intron_mean_P2AB + pseudocount),
    log10_ei          = log10(ei_at_pseudocount)
  )

p_pseudocount_sensitivity <- ggplot(sensitivity_df, aes(x = pseudocount, y = log10_ei, group = gene_id)) +
  geom_line(alpha = 0.4, color = "#CC6677", linewidth = 0.5) +
  geom_vline(xintercept = INTRON_PSEUDOCOUNT, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  annotate("text", x = INTRON_PSEUDOCOUNT,
           y = max(sensitivity_df$log10_ei, na.rm = TRUE),
           label = paste0("Used pseudocount\n(", signif(INTRON_PSEUDOCOUNT, 3), ")"),
           hjust = -0.1, vjust = 1, size = 3, color = "grey30") +
  scale_x_log10(labels = scales::comma) +
  labs(
    x        = "Candidate intron pseudocount (log10 scale)",
    y        = "log10(exon RPKM / (intron RPKM + pseudocount))",
    title    = "Pseudocount sensitivity \u2014 high P2AB exon, zero P2AB intron genes",
    subtitle = paste0("Top ", nrow(inflation_candidates),
                      " genes by P2AB exon RPKM among those with intron RPKM = 0  |  ",
                      "each line = one gene's E/I across the pseudocount sweep")
  ) +
  base_theme +
  theme(legend.position = "none")

print(p_pseudocount_sensitivity)

top_inflation_gene <- inflation_candidates %>% dplyr::slice_head(n = 1)
example_pseudocounts <- c(0.0001, 0.001, 0.01, 0.1, 1)
cat("\nWorst-case example (gene:", top_inflation_gene$gene_id,
    ", exon RPKM =", round(top_inflation_gene$exon_mean_P2AB, 2), ", intron RPKM = 0):\n")
for (pc in example_pseudocounts) {
  ei_val <- top_inflation_gene$exon_mean_P2AB / (top_inflation_gene$intron_mean_P2AB + pc)
  cat(sprintf("  pseudocount=%.4f  ->  E/I=%.1f  (log10=%.2f)\n", pc, ei_val, log10(ei_val)))
}

##--------------------------------------------------------------## 12. Impact of removing zero-intron genes entirely
## Compares the current approach (KEEP zero-intron genes, pseudocount
## them — confirmed as the right call based on the sensitivity plot
## above and the by-length E/I plots: the general trend holds either
## way, and removing intron-less genes was explicitly rejected)
## against the alternative (REMOVE them entirely).

message("\nChecking impact of removing zero-intron genes entirely...")

ei_intron_nonzero <- ei_df %>%
  dplyr::filter(intron_WT_Input > 0, intron_mean_P2AB > 0)

cat("\nOf", nrow(ei_df), "WT_Input-anchored genes:\n")
cat("  ", nrow(ei_df) - sum(ei_df$intron_WT_Input > 0), "have zero Input intron signal\n")
cat("  ", nrow(ei_df) - sum(ei_df$intron_mean_P2AB > 0), "have zero P2AB intron signal\n")
cat("  ", nrow(ei_intron_nonzero), "have non-zero intron signal in BOTH (retained if zero-intron genes are removed)\n")

intron_removal_df <- dplyr::bind_rows(
  ei_df %>%
    dplyr::filter(ei_P2AB > 0) %>%
    dplyr::transmute(gene_id, log10_ei = log10(ei_P2AB), source = "P2AB",
                     intron_filter = "Zero intron KEPT (current)"),
  ei_intron_nonzero %>%
    dplyr::filter(ei_P2AB > 0) %>%
    dplyr::transmute(gene_id, log10_ei = log10(ei_P2AB), source = "P2AB",
                     intron_filter = "Zero intron REMOVED"),
  ei_df %>%
    dplyr::transmute(gene_id, log10_ei = log10(ei_WT_Input), source = "Input",
                     intron_filter = "Zero intron KEPT (current)"),
  ei_intron_nonzero %>%
    dplyr::transmute(gene_id, log10_ei = log10(ei_WT_Input), source = "Input",
                     intron_filter = "Zero intron REMOVED")
) %>%
  dplyr::filter(!is.infinite(log10_ei)) %>%
  dplyr::mutate(
    source        = factor(source, levels = c("Input","P2AB")),
    intron_filter = factor(intron_filter, levels = c("Zero intron KEPT (current)", "Zero intron REMOVED"))
  )

p_intron_removal_impact <- ggplot(intron_removal_df, aes(x = log10_ei, fill = intron_filter, color = intron_filter)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values = c("Zero intron KEPT (current)" = "#CC6677",
                               "Zero intron REMOVED"        = "grey50"), name = NULL) +
  scale_color_manual(values = c("Zero intron KEPT (current)" = "#CC6677",
                                "Zero intron REMOVED"        = "grey50"), name = NULL) +
  facet_wrap(~ source, nrow = 1, scales = "free_y") +
  labs(
    x        = "log10 E/I",
    y        = "Density",
    title    = "Impact of removing zero-intron genes on the E/I distribution",
    subtitle = "Rose = current approach (zero-intron genes kept, pseudocounted)  |  Grey = zero-intron genes removed entirely"
  ) +
  base_theme

print(p_intron_removal_impact)

message("\nSensitivity/impact analysis complete. Key objects: sensitivity_df, intron_removal_df")

##--------------------------------------------------------------## 13. Load WT TPM and compute expression bins
## First use of expression stratification in this script. Same
## tertile convention as the EU-seq E/I script: bottom-tertile-of-
## expressed genes + all TPM==0/NA genes fold into "Low/Unexpressed".

EXPR_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
EXPR_CATS <- c("Low/Unexpressed","Medium","High")

message("\nLoading WT TPM and computing expression tertiles...")

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
    expr_bin = factor(expr_bin, levels = EXPR_CATS)
  )

cat("Expression bins (all genes with TPM data):\n")
print(table(tpm_df$expr_bin))

ei_df <- ei_df %>%
  dplyr::left_join(tpm_df %>% dplyr::select(gene_id, expr_bin), by = "gene_id")

cat("\nWT_Input-anchored genes matched to a TPM expression bin:",
    sum(!is.na(ei_df$expr_bin)), "/", nrow(ei_df), "\n")
cat("Genes per expression bin:\n")
print(table(ei_df$expr_bin, useNA = "no"))

##--------------------------------------------------------------## 14. Effective E/I plot — Genomic / Input / P2AB, by expression level
## Same include_input switch as section 10. Genomic reference =
## Input-expressed genes only (the confirmed approach). Only
## Low/Unexpressed and High are shown (Medium computed/counted above
## but not plotted), matching the EU-seq script's convention.

message("\nBuilding RIP-seq effective E/I comparison, by expression level...")

build_ripseq_ei_expr_comparison <- function(include_input = TRUE) {
  
  active_levels <- active_source_levels(include_input)
  
  plot_df <- dplyr::bind_rows(
    genomic_df_input_expressed %>%
      dplyr::left_join(tpm_df %>% dplyr::select(gene_id, expr_bin), by = "gene_id") %>%
      dplyr::transmute(gene_id, expr_bin, ei_ratio = ei_genomic, source = "Genomic\n(reference)"),
    if (include_input) {
      ei_df %>% dplyr::transmute(gene_id, expr_bin, ei_ratio = ei_WT_Input, source = "Input\n(WT)")
    } else NULL,
    ei_df %>%
      dplyr::transmute(gene_id, expr_bin, ei_ratio = ei_P2AB, source = "P2AB\n(Pol2 IP)")
  ) %>%
    dplyr::filter(ei_ratio > 0, !is.na(expr_bin), expr_bin %in% c("Low/Unexpressed","High")) %>%
    dplyr::mutate(
      log10_ei = log10(ei_ratio),
      source   = factor(source, levels = active_levels),
      expr_bin = factor(expr_bin, levels = c("Low/Unexpressed","High"))
    )
  
  x_range <- range(plot_df$log10_ei, na.rm = TRUE)
  
  panels <- lapply(levels(plot_df$expr_bin), function(bin_name) {
    sub_df <- plot_df %>% dplyr::filter(expr_bin == bin_name)
    build_ripseq_ei_panel_generic(sub_df, bin_name, x_range, active_levels)
  })
  
  final_plot <- (patchwork::wrap_plots(panels, nrow = 1) + patchwork::plot_layout(guides = "collect")) +
    patchwork::plot_annotation(
      title    = paste0("RIP-seq effective E/I: Genomic vs ",
                        if (include_input) "Input vs " else "",
                        "Pol2 IP (P2AB), by expression level"),
      subtitle = paste0(
        "Grey dashed = genomic E/I (Input-expressed genes only)  |  ",
        "Solid = RIP-seq exon RPKM / intron RPKM  |  Input-anchored, P2AB exon=0 kept  |  ",
        "Low/Unexpressed and High WT TPM tertiles only"
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(color = "grey40", size = 11))
    ) &
    theme(legend.position = "top")
  
  print(final_plot)
  invisible(list(plot_df = plot_df, final_plot = final_plot))
}

message("\n\u2500\u2500 Expression-level E/I: Genomic + Input + P2AB \u2500\u2500")
result_expr_with_input <- build_ripseq_ei_expr_comparison(include_input = TRUE)

message("\n\u2500\u2500 Expression-level E/I: Genomic + P2AB only (no Input line) \u2500\u2500")
result_expr_p2ab_only <- build_ripseq_ei_expr_comparison(include_input = FALSE)

##--------------------------------------------------------------## 15. Effective E/I plots — Genomic / Input / P2AB, by intron length,
## split by expression (High-only, then Low-only). Same include_input
## switch as sections 10 and 14. The plain length-only plot from
## section 10 is UNCHANGED and stays in — these are additional.

message("\nBuilding RIP-seq length-faceted plots, split by expression group...")

build_length_by_expr_source_df <- function(include_input) {
  dplyr::bind_rows(
    genomic_df_input_expressed %>%
      dplyr::left_join(tpm_df %>% dplyr::select(gene_id, expr_bin), by = "gene_id") %>%
      dplyr::transmute(gene_id, expr_bin, intron_cat, ei_ratio = ei_genomic, source = "Genomic\n(reference)"),
    if (include_input) {
      ei_df %>% dplyr::transmute(gene_id, expr_bin, intron_cat, ei_ratio = ei_WT_Input, source = "Input\n(WT)")
    } else NULL,
    ei_df %>%
      dplyr::transmute(gene_id, expr_bin, intron_cat, ei_ratio = ei_P2AB, source = "P2AB\n(Pol2 IP)")
  ) %>%
    dplyr::filter(ei_ratio > 0, !is.na(expr_bin)) %>%
    dplyr::mutate(
      log10_ei = log10(ei_ratio),
      source   = factor(source, levels = active_source_levels(include_input))
    )
}

## Precomputed once per switch state, not once per (expr_filter x switch)
## call — 2 computations total instead of 4
ei_length_by_expr_df_full     <- build_length_by_expr_source_df(include_input = TRUE)
ei_length_by_expr_df_p2ab_only <- build_length_by_expr_source_df(include_input = FALSE)

build_ripseq_length_by_expr_plot <- function(expr_filter, expr_label, include_input = TRUE) {
  
  active_levels <- active_source_levels(include_input)
  full_df <- if (include_input) ei_length_by_expr_df_full else ei_length_by_expr_df_p2ab_only
  sub_all <- full_df %>% dplyr::filter(expr_bin == expr_filter)
  
  if (nrow(sub_all) == 0) {
    message("  Skipping ", expr_label, " — no data.")
    return(invisible(NULL))
  }
  
  ## Only build panels for intron_cat levels that actually have rows —
  ## an empty panel would break range() on empty data
  cat_levels <- levels(sub_all$intron_cat)
  cat_levels <- cat_levels[cat_levels %in% unique(as.character(sub_all$intron_cat))]
  if (length(cat_levels) < length(levels(sub_all$intron_cat))) {
    missing_cats <- setdiff(levels(sub_all$intron_cat), cat_levels)
    message("  Note: ", expr_label, " expression has no genes in ",
            paste(missing_cats, collapse = ", "), " — panel(s) skipped")
  }
  
  shared_x_range <- range(sub_all$log10_ei, na.rm = TRUE)
  
  panels <- lapply(cat_levels, function(cat_name) {
    sub_df <- sub_all %>% dplyr::filter(intron_cat == cat_name)
    build_ripseq_ei_panel_generic(sub_df, cat_name, shared_x_range, active_levels)
  })
  
  n_genes <- ei_df %>% dplyr::filter(expr_bin == expr_filter) %>% nrow()
  
  final_plot <- (patchwork::wrap_plots(panels, nrow = 1) + patchwork::plot_layout(guides = "collect")) +
    patchwork::plot_annotation(
      title    = paste0("RIP-seq effective E/I: Genomic vs ",
                        if (include_input) "Input vs " else "",
                        "Pol2 IP (P2AB), by intron length \u2014 ", expr_label, " expression"),
      subtitle = paste0(
        "Grey dashed = genomic E/I (Input-expressed genes only)  |  ",
        "Solid = RIP-seq exon RPKM / intron RPKM  |  Input-anchored, P2AB exon=0 kept  |  ",
        expr_label, " WT TPM tertile only  |  n = ", n_genes, " genes"
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(color = "grey40", size = 11))
    ) &
    theme(legend.position = "top")
  
  print(final_plot)
  invisible(final_plot)
}

# lengthxlowexpr
message("\n\u2500\u2500 RIP-seq length-by-expression: High expression, with Input \u2500\u2500")
result_length_high <- build_ripseq_length_by_expr_plot("High", "High", include_input = TRUE)

message("\n\u2500\u2500 RIP-seq length-by-expression: Low/Unexpressed, with Input \u2500\u2500")
result_length_low <- build_ripseq_length_by_expr_plot("Low/Unexpressed", "Low/Unexpressed", include_input = TRUE)

message("\n\u2500\u2500 RIP-seq length-by-expression: High expression, P2AB only \u2500\u2500")
result_length_high_p2ab_only <- build_ripseq_length_by_expr_plot("High", "High", include_input = FALSE)

message("\n\u2500\u2500 RIP-seq length-by-expression: Low/Unexpressed, P2AB only \u2500\u2500")
result_length_low_p2ab_only <- build_ripseq_length_by_expr_plot("Low/Unexpressed", "Low/Unexpressed", include_input = FALSE)

message("\nAnalysis complete. Key objects: ei_df, result_length_high, result_length_low,",
        " result_length_high_p2ab_only, result_length_low_p2ab_only")