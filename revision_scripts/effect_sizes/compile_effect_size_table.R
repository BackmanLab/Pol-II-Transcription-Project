compile_effect_size_table <- function(
    in_dir       = "/projects/b1042/BackmanLab/Lucas/102726_effect_size",
    out_dir      = "/home/lmc0633/Transcription Paper/101426_revision_scripts/effect sizes/effect_sizes_tables_plots",
    out_file     = file.path(out_dir, "Effect_sizes_SI_table.csv"),
    display_file = file.path(out_dir, "Effect_sizes_SI_table_display.csv"),
    # the CSV written by each effect-size script, relative to in_dir, in the order they appear in the table
    sources = c(PWS                   = "PWS/PWS_effect_sizes.csv",
                ChromSTEM             = "chromSTEM/ChromSTEM_effect_sizes.csv",
                ChromSTEM_domain_loss = "chromSTEM/ChromSTEM_domain_loss.csv",
                SMLM                  = "SMLM/SMLM_effect_sizes.csv",
                EUSTORM               = "SMLM/EUSTORM_effect_sizes.csv",
                EUseq                 = "EUseq/EUseq_effect_sizes.csv",
                TAD                   = "HiC/HiC_TAD_insulation_effect_sizes.csv",
                Compartments          = "HiC/HiC_compartment_effect_sizes.csv",
                Loops                 = "HiC/HiC_loop_effect_sizes.csv",
                Scaling               = "HiC/HiC_scaling_effect_sizes.csv")) {
  
  # ---- input checks: every file must exist and carry the columns written by its script ------------
  paths <- setNames(file.path(in_dir, sources), names(sources))
  if (any(!file.exists(paths))) stop("Input file(s) not found (run the matching effect-size script first):\n", paste(paths[!file.exists(paths)], collapse = "\n"))
  need <- list(
    PWS = c("assay", "comparison", "measure", "n_reps_ctrl", "n_reps_trt", "n_roi_ctrl", "n_roi_trt", "hedges_g", "ci_lo", "ci_hi"),
    ChromSTEM = c("assay", "comparison", "measure", "n_ctrl", "n_trt", "hedges_g", "ci_lo", "ci_hi"),
    ChromSTEM_domain_loss = c("assay", "classification", "comparison", "category", "n_ctrl", "n_trt", "ratio", "ci_lo", "ci_hi", "pct_change"),
    SMLM = c("assay", "comparison", "measure", "size_bins_nm", "n_cells_ctrl", "n_cells_trt", "hedges_g", "ci_lo", "ci_hi"),
    EUSTORM = c("assay", "comparison", "measure", "panel", "level", "n_ctrl", "n_trt", "n_cells_ctrl", "n_cells_trt", "hedges_g", "ci_lo", "ci_hi"),
    EUseq = c("assay", "comparison", "expr_group", "n_genes", "median_lfc", "lfc_lo", "lfc_hi", "median_hedges_g", "g_lo", "g_hi", "frac_down", "frac_down_lo", "frac_down_hi"),
    TAD = c("assay", "comparison", "tad_positions", "ctcf_status", "expr_bin", "n_tads", "hedges_g", "ci_lo", "ci_hi"),
    Compartments = c("assay", "comparison", "expr_bin", "n_chroms", "pct_change_cs", "pct_lo", "pct_hi", "hedges_g_paired", "g_lo", "g_hi"),
    Loops = c("assay", "ds_group", "comparison", "loop_set", "expr_bin", "n_loops", "hedges_g", "g_lo", "g_hi", "frac_lost", "frac_lost_lo", "frac_lost_hi",
              "lost_gained_ratio", "ratio_lo", "ratio_hi"),
    Scaling = c("assay", "method", "comparison", "family", "stratum", "n_genes", "n_bins_ctrl", "n_bins_trt", "delta_S_wt_minus_trt", "dS_lo", "dS_hi",
                "hedges_g", "g_lo", "g_hi", "grad_wt_minus_trt", "grad_lo", "grad_hi"))
  if (!setequal(names(need), names(sources))) stop("sources and the column requirements disagree: ", paste(setdiff(union(names(need), names(sources)), intersect(names(need), names(sources))), collapse = ", "))
  dat <- lapply(names(sources), function(s) {
    d <- read.csv(paths[[s]], stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(d) == 0) stop("No rows in ", paths[[s]])
    if (any(!need[[s]] %in% names(d))) stop("Columns missing in ", paths[[s]], ": ", paste(need[[s]][!need[[s]] %in% names(d)], collapse = ", "),
                                            "\nColumns present: ", paste(names(d), collapse = ", "))
    d
  })
  names(dat) <- names(sources)
  for (o in unique(dirname(c(out_file, display_file)))) {   # the export folder is created if it does not exist yet, and must be writable
    dir.create(o, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(o) || file.access(o, 2) != 0) stop("Cannot write to the export folder: ", o)
  }
  
  # ---- which perturbation a comparison is: every comparison label names the treated condition -------------
  drug_of <- function(cmp, src) {
    out <- ifelse(grepl("actd", cmp, ignore.case = TRUE), "ActD", ifelse(grepl("polr2a|aux", cmp, ignore.case = TRUE), "POLR2A", NA_character_))
    if (anyNA(out)) stop("Comparison label(s) in ", src, " that are neither POLR2A-degraded nor ActD: ", paste(unique(cmp[is.na(out)]), collapse = ", "))
    out
  }
  pm <- function(a, b, what) sprintf("%s %d vs %d", what, as.integer(a), as.integer(b))   # "control vs treated" counts
  
  # ---- one long table: a row per (source, set, measure, stratum, effect type, perturbation) -------------
  long <- list()
  add <- function(src, d, set, measure, stratum, effect_type, effect, lo, hi, unit, digits, note) {
    n <- nrow(d)
    rep1 <- function(x) if (length(x) == 1) rep(x, n) else x
    long[[length(long) + 1]] <<- data.frame(
      source = src, assay = d$assay, set = rep1(set), measure = rep1(measure), stratum = rep1(stratum), effect_type = effect_type,
      drug = drug_of(d$comparison, src), effect = as.numeric(effect), lo = as.numeric(lo), hi = as.numeric(hi),
      unit = rep1(unit), digits = digits, note = note, stringsAsFactors = FALSE)
  }
  G <- "Hedges' g (treated - control)"
  d <- dat$PWS
  add("PWS", d, "", d$measure, "All", G, d$hedges_g, d$ci_lo, d$ci_hi,
      paste0(pm(d$n_roi_ctrl, d$n_roi_trt, "ROIs"), "; replicates ", d$n_reps_ctrl, " vs ", d$n_reps_trt), 2,
      "Unit: ROIs; replicates fixed and ROIs resampled within each replicate")
  d <- dat$ChromSTEM
  add("ChromSTEM", d, "", d$measure, "All", G, d$hedges_g, d$ci_lo, d$ci_hi, pm(d$n_ctrl, d$n_trt, "domains"), 2,
      "Unit: domains from one imaged cell per condition; domains resampled")
  d <- dat$ChromSTEM_domain_loss
  add("ChromSTEM_domain_loss", d, d$classification, "Domain count", d$category, "Ratio of domain counts (treated / control)", d$ratio, d$ci_lo, d$ci_hi,
      pm(d$n_ctrl, d$n_trt, "domains"), 2, "Exact conditional-binomial CI; ratio < 1 = domains lost")
  d <- dat$SMLM
  add("SMLM", d, d$size_bins_nm, d$measure, "All", G, d$hedges_g, d$ci_lo, d$ci_hi, pm(d$n_cells_ctrl, d$n_cells_trt, "cells"), 2,
      "Unit: cells; cells resampled within condition; size bins in nm (small / average / large)")
  d <- dat$EUSTORM
  add("EUSTORM", d, paste0(d$panel, ", ", d$level, " level"), d$measure, "All", G, d$hedges_g, d$ci_lo, d$ci_hi,
      paste0(pm(d$n_ctrl, d$n_trt, ifelse(d$level == "roi", "ROIs", "cells")), "; cells ", d$n_cells_ctrl, " vs ", d$n_cells_trt), 2,
      "Unit: ROIs (level = roi) or per-cell means (level = cell); cells resampled within condition and replicate")
  d <- dat$EUseq
  add("EUseq", d, "", "Median log2 fold change (scaled NoEU-subtracted RPKM)", d$expr_group, "Median log2 fold change (treated - control)",
      d$median_lfc, d$lfc_lo, d$lfc_hi, paste0("genes ", d$n_genes), 2, "chrM-scaled; bootstrap over genes; ActD EU-seq discarded (too few reads)")
  add("EUseq", d, "", "Median per-gene Hedges' g", d$expr_group, G, d$median_hedges_g, d$g_lo, d$g_hi, paste0("genes ", d$n_genes), 2,
      "Per-gene g = LFC / within-group pooled SD (n = 3 vs 3); not comparable in size with the other g values")
  add("EUseq", d, "", "Fraction of genes with LFC < -1", d$expr_group, "Fraction of genes (exact binomial CI)", d$frac_down, d$frac_down_lo, d$frac_down_hi,
      paste0("genes ", d$n_genes), 3, "")
  d <- dat$TAD
  add("TAD", d, paste0("TAD positions: ", d$tad_positions, "; CTCF: ", d$ctcf_status), "Boundary insulation score (z)", d$expr_bin, G, d$hedges_g, d$ci_lo, d$ci_hi,
      paste0("TADs ", d$n_tads), 2, "Unit: TADs, paired across conditions; higher insulation z = weaker insulation")
  d <- dat$Compartments
  add("Compartments", d, "", "Genome-wide compartment strength (CS)", d$expr_bin, "% change in CS (treated vs control)", d$pct_change_cs, d$pct_lo, d$pct_hi,
      paste0("chromosomes ", d$n_chroms), 1, "Bootstrap over chromosomes; negative = weaker compartmentalization")
  add("Compartments", d, "", "Compartment strength per chromosome", d$expr_bin, "Hedges' g, paired (treated - control)", d$hedges_g_paired, d$g_lo, d$g_hi,
      paste0("chromosomes ", d$n_chroms), 2, "g_z = mean / SD of the per-chromosome differences")
  d <- dat$Loops
  lset <- paste0(d$ds_group, "; ", d$loop_set)
  add("Loops", d, lset, "Loop strength, log10(SCALE contacts + 1)", d$expr_bin, G, d$hedges_g, d$g_lo, d$g_hi, paste0("loops ", d$n_loops), 2,
      "Unit: loops, paired; negative = weaker loops")
  add("Loops", d, lset, "Fraction of loops lost", d$expr_bin, "Fraction of loops (exact binomial CI)", d$frac_lost, d$frac_lost_lo, d$frac_lost_hi,
      paste0("loops ", d$n_loops), 3, "")
  add("Loops", d, lset, "Lost : gained loops", d$expr_bin, "Ratio of lost to gained loops (exact CI)", d$lost_gained_ratio, d$ratio_lo, d$ratio_hi,
      paste0("loops ", d$n_loops), 2, "")
  d <- dat$Scaling
  sset <- paste0(d$method, "; ", d$family)
  isg <- is.finite(d$hedges_g); isr <- is.finite(d$grad_wt_minus_trt)
  unit_sc <- paste0("genes ", d$n_genes, "; bins/pools ", d$n_bins_ctrl, " vs ", d$n_bins_trt)
  add("Scaling", d[isg, ], sset[isg], "Contact-scaling slope S", d$stratum[isg], G, d$hedges_g[isg], d$g_lo[isg], d$g_hi[isg], unit_sc[isg], 2,
      "Unit: bins (method 1) or pools (method 2), resampled; positive = S less negative")
  add("Scaling", d[isg, ], sset[isg], "Contact-scaling slope S", d$stratum[isg], "Delta S (Wildtype - treated)", d$delta_S_wt_minus_trt[isg], d$dS_lo[isg], d$dS_hi[isg],
      unit_sc[isg], 3, "Sign of Fig. 4F-G; negative = S rises after treatment")
  add("Scaling", d[isr, ], sset[isr], "Contact-scaling slope S", d$stratum[isr], "Gradient of Delta S across percentile (Wildtype - treated)", d$grad_wt_minus_trt[isr],
      d$grad_lo[isr], d$grad_hi[isr], unit_sc[isr], 3, "Genes resampled within bins; change in Delta S from the 0th to the 100th percentile")
  long <- do.call(rbind, long)
  
  # ---- self-checks on the live data ----------------------------------------------------------------------
  cat("\nRows read and rows entering the table, per source:\n")
  print(data.frame(source = names(dat), rows_in = sapply(dat, nrow), long_rows = as.integer(table(factor(long$source, levels = names(dat)))), row.names = NULL))
  if (any(is.finite(long$lo) & is.finite(long$hi) & long$lo > long$hi)) stop("Rows with lo > hi: ", paste(unique(long$source[is.finite(long$lo) & is.finite(long$hi) & long$lo > long$hi]), collapse = ", "))
  ci_ok <- is.finite(long$lo) & is.finite(long$hi) & is.finite(long$effect)
  out_ci <- ci_ok & (long$effect < long$lo - 1e-9 | long$effect > long$hi + 1e-9)
  cat("Estimates outside their own CI:", sum(out_ci), "of", sum(ci_ok), "rows with a CI\n")
  if (any(out_ci)) { print(long[out_ci, c("source", "set", "measure", "stratum", "effect_type", "drug", "effect", "lo", "hi")], digits = 3); warning(sum(out_ci), " estimate(s) outside their CI; see the list above") }
  key <- paste(long$source, long$assay, long$set, long$measure, long$stratum, long$effect_type, sep = "\r")
  if (any(duplicated(paste(key, long$drug)))) stop("More than one row for the same table cell and perturbation: ", paste(unique(sub("\r", " | ", key[duplicated(paste(key, long$drug))])), collapse = "\n"))
  cat("Self-check passed: every perturbation has at most one row per table cell, every comparison is POLR2A-degraded or ActD\n")
  
  # ---- wide table, POLR2A-degraded and ActD side by side --------------------------------------------------
  fmt <- function(e, l, h, dg) ifelse(is.finite(e), ifelse(is.finite(l) & is.finite(h),
                                                           sprintf("%.*f [%.*f, %.*f]", dg, e, dg, l, dg, h), sprintf("%.*f", dg, e)), "n/a")
  ukey <- unique(key)
  one <- function(drug) {
    i <- match(paste(ukey, drug), paste(key, long$drug))
    data.frame(effect = long$effect[i], lo = long$lo[i], hi = long$hi[i], unit = long$unit[i], stringsAsFactors = FALSE)
  }
  p <- one("POLR2A"); a <- one("ActD")
  first <- match(ukey, key)
  wide <- data.frame(
    source = long$source[first], assay = long$assay[first], set = long$set[first], measure = long$measure[first], stratum = long$stratum[first],
    effect_type = long$effect_type[first],
    POLR2A_effect = p$effect, POLR2A_lo = p$lo, POLR2A_hi = p$hi, POLR2A_unit = p$unit,
    ActD_effect = a$effect, ActD_lo = a$lo, ActD_hi = a$hi, ActD_unit = a$unit,
    note = long$note[first], stringsAsFactors = FALSE)
  dg <- long$digits[first]
  wide$POLR2A_display <- fmt(wide$POLR2A_effect, wide$POLR2A_lo, wide$POLR2A_hi, dg)
  wide$ActD_display   <- fmt(wide$ActD_effect, wide$ActD_lo, wide$ActD_hi, dg)
  stopifnot(nrow(wide) == length(ukey))
  empty <- wide$POLR2A_display == "n/a" & wide$ActD_display == "n/a"   # e.g. a Lost:Gained ratio with no gained loops in either perturbation
  cat("Rows with no estimate for either perturbation (dropped):", sum(empty), "\n")
  wide <- wide[!empty, , drop = FALSE]
  cat("Rows with POLR2A-degraded only:", sum(wide$ActD_display == "n/a"), "| ActD only:", sum(wide$POLR2A_display == "n/a"), "| both:", sum(wide$ActD_display != "n/a" & wide$POLR2A_display != "n/a"), "\n")
  
  disp <- data.frame(
    Assay = wide$assay, Set = wide$set, Measure = wide$measure, Stratum = wide$stratum, `Effect size` = wide$effect_type,
    `POLR2A-degraded` = wide$POLR2A_display, `ActD (1 hr)` = wide$ActD_display,
    `N (POLR2A-degraded)` = ifelse(is.na(wide$POLR2A_unit), "", wide$POLR2A_unit), `N (ActD)` = ifelse(is.na(wide$ActD_unit), "", wide$ActD_unit),
    Note = wide$note, check.names = FALSE, stringsAsFactors = FALSE)
  
  cat("\nColumns of the SI table: effect size [95% CI]; g = treated - control for every Hedges' g; units differ by assay (see N and Note)\n\n")
  print(disp[, c("Assay", "Measure", "Stratum", "Effect size", "POLR2A-degraded", "ActD (1 hr)")], right = FALSE, row.names = FALSE)
  write.csv(wide, out_file, row.names = FALSE)
  write.csv(disp, display_file, row.names = FALSE)
  cat("Wrote:", out_file, "\nWrote:", display_file, "\n")
  invisible(wide)
}

effect_table <- compile_effect_size_table()