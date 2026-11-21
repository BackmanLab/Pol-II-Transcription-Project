##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Fork: pooled loop change across the full range of each variable, at finer
# resolution than thirds. Wildtype loops (Wildtype contacts >= MIN_WT_CONTACTS) whose
# anchor overlaps an expressed gene (TPM > 0 in Wildtype cells; the highest-TPM gene
# that overlaps either anchor) are
# ordered by one variable at a time, and the contacts are pooled:
#   pooled log2FC = log2(sum of POLR2A-degraded contacts / sum of Wildtype contacts)
# with a bootstrap over loops for the 95% interval. Pooling avoids the per-loop
# counting noise of log2FC. Variables: the change in EU-seq and in total RNA-seq after
# POLR2A degradation over the loop anchors and over the whole anchor gene
# (log2((POLR2A-degraded + pc) / (Wildtype + pc)), pc = 1st percentile of the positive
# values, from the hg38 chrM-scaled bigWigs) and the changes in H3K27ac, H3K4me3 and
# H3K9me3. Loops with no Wildtype signal in either assay are removed.
#   1. Deciles of each variable (low = lowest value; for the change variables, low =
#      largest decrease), pooled log2FC with 95% interval.
#   2. Top decile minus bottom decile, bootstrap 95% interval and p.
#   3. A rolling curve: the pooled log2FC in a window of the loops ranked by the
#      variable (window = 2 x ROLL_HALF_WIDTH of the loops, moved in steps of ROLL_STEP),
#      with a bootstrap 95% band.
# Two plots: (A) one panel per variable (value of the variable on x); (B) all variables
# on one plot, each ranked as a percentile (0 = largest decrease, 100 = smallest decrease
# or an increase), coloured by variable.
# A third plot (C, at the end) uses the EU-seq LEVEL instead of the change as the variable:
# the signal over the anchors and the whole anchor gene in POLR2A-degraded cells (what is left
# after degradation) and in Wildtype cells (where it started), same loops and same pooling.
# A fourth plot (D, at the end) repeats plot A for loops called in POLR2A-degraded cells
# (shared and degraded-only; degraded contacts >= MIN_POL2_CONTACTS) with an expressed anchor gene.
# Needs: loop_regression_main.R (sourced below) and the EU-seq and total RNA-seq bigWigs.
##################################################################

# ── Cores available to this session, passed on to loop_regression_main.R ──────
n_cores <- suppressWarnings(as.integer(c(
  Sys.getenv("SLURM_CPUS_PER_TASK", NA),
  Sys.getenv("SLURM_CPUS_ON_NODE", NA),
  sub("\\(.*", "", Sys.getenv("SLURM_JOB_CPUS_PER_NODE", NA)),
  system("nproc", intern = TRUE),
  parallel::detectCores())))
n_cores <- n_cores[!is.na(n_cores) & n_cores >= 1][1]
stopifnot(length(n_cores) == 1)
Sys.setenv(SLURM_CPUS_PER_TASK = as.character(n_cores))
message("Cores available: ", n_cores, "  (SLURM_CPUS_PER_TASK set to ", Sys.getenv("SLURM_CPUS_PER_TASK"), ")")

ASSAY <- "hic"

# ── Plot output ───────────────────────────────────────────────────────────────
# Every plot is also saved here as an SVG (the folder must already exist).
PLOT_DIR <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/chip-seq/plots/regression_pooled_loops"
if (!dir.exists(PLOT_DIR)) stop("PLOT_DIR not found: ", PLOT_DIR)
rm(features) 

# ── Load the shared main script from the directory this script is in ──────────
# (skipped if its objects are already in this session)
find_main_script <- function(main_name = "loop_regression_main.R") {
  cand <- character(0)
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) cand <- c(cand, dirname(f))                       # Rscript
  for (i in rev(seq_len(sys.nframe()))) {                               # source()
    of <- sys.frame(i)$ofile
    if (is.character(of)) { cand <- c(cand, dirname(of)); break }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) cand <- c(cand, dirname(p))                          # RStudio editor
  }
  cand <- c(cand, getwd())                                              # working directory
  hit  <- file.path(normalizePath(cand, mustWork = FALSE), main_name)
  hit  <- hit[file.exists(hit)]
  if (length(hit) == 0) {
    stop("Missing input file: ", main_name, " (looked in: ",
         paste(unique(normalizePath(cand, mustWork = FALSE)), collapse = ", "),
         "). Save it next to this script or setwd() to that directory.")
  }
  hit[1]
}
MAIN_SCRIPT <- find_main_script()
main_stamp <- function(f) paste(file.info(f)$mtime, file.info(f)$size)   # changes whenever main is saved
if (exists("features") && exists("wt_reg") && exists("run_lm_models") &&
    identical(get0("MAIN_STAMP"), main_stamp(MAIN_SCRIPT)) &&
    identical(get0("FEATURES_ASSAY"),
              get0("ASSAY", ifnotfound = Sys.getenv("LOOP_ASSAY", unset = "hic")))) {
  message("Using features already in this session (assay: ", FEATURES_ASSAY, ", loops: ",
          DS_GROUP, ", maps: ", HIC_GROUP, "). Remove 'features' to rebuild.")
} else {
  source(MAIN_SCRIPT)
  MAIN_STAMP <- main_stamp(MAIN_SCRIPT)
}


# ── Transcription change at the anchors and over the anchor gene ───────────────
# Loops whose anchor overlaps a gene (highest-TPM gene, tol_bp = 0) get the change in
# EU-seq and in total RNA-seq after POLR2A degradation, over the anchors and over
# the anchor gene (EU-seq: NoEU-subtracted bigWigs, already x eu_sf; total RNA-seq:
# scaled_mito bigWigs).
EU_BASE <- "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq"
EU_BW   <- c(POLR2A = file.path(EU_BASE, "euseq/coverage/subtracted/POLR2A_8hrs_aux.subtracted.bigWig"),
             WT     = file.path(EU_BASE, "euseq/coverage/subtracted/WT.subtracted.bigWig"))
TOT_BW  <- c(POLR2A = file.path(EU_BASE, "total/coverage/scaled_mito/POLR2A_8hrs_aux.scaled.bigWig"),
             WT     = file.path(EU_BASE, "total/coverage/scaled_mito/WT.scaled.bigWig"))
check_files(c(EU_BW, TOT_BW))

bn_reg <- wt_reg %>%
  dplyr::select(-dplyr::any_of(c("anchor_gene_tpm", "expr_bin", "anchor_gene_id"))) %>%
  link_loops_to_genes(gene_gr, tol_bp = 0L) %>%
  dplyr::filter(!is.na(anchor_gene_id))
n_gene_loops <- nrow(bn_reg)
# Only loops whose anchor gene is expressed in Wildtype cells (TPM > 0, the same
# "expressed" rule as the expression bins in loop_regression_main.R)
bn_reg <- dplyr::filter(bn_reg, anchor_gene_tpm > 0)
message("\nWildtype loops with WT contacts >= ", MIN_WT_CONTACTS, ": ", nrow(wt_reg),
        "   with an anchor overlapping a gene: ", n_gene_loops,
        "   anchor gene expressed (TPM > 0): ", nrow(bn_reg))
stopifnot(nrow(bn_reg) > 0, !anyNA(match(bn_reg$anchor_gene_id, gene_gr$gene_id)))

bn_ug    <- unique(bn_reg$anchor_gene_id)
bn_ug_gr <- gene_gr[match(bn_ug, gene_gr$gene_id)]
get_signals <- function(BW, label) {
  message("  ", label)
  anchors <- sapply(names(BW), function(nm) {
    message("    ", nm, " : ", BW[[nm]])
    a1 <- anchor_signal(BW[[nm]], bn_reg$chr1, bn_reg$start1, bn_reg$end1)
    a2 <- anchor_signal(BW[[nm]], bn_reg$chr2, bn_reg$start2, bn_reg$end2)
    (a1 + a2) / 2
  })
  genes <- sapply(names(BW), function(nm) {
    g <- anchor_signal(BW[[nm]], as.character(seqnames(bn_ug_gr)), start(bn_ug_gr) - 1L, end(bn_ug_gr))
    g[match(bn_reg$anchor_gene_id, bn_ug)]
  })
  stopifnot(nrow(anchors) == nrow(bn_reg), nrow(genes) == nrow(bn_reg),
            all(is.finite(anchors)), all(is.finite(genes)))
  if (any(anchors < 0) || any(genes < 0)) {
    stop("Negative ", label, " signal in the bigWigs (anchors ", sum(anchors < 0), ", genes ",
         sum(genes < 0), "); the hg38 chrM-scaled bigWigs are not negative, so the paths point at the wrong files")
  }
  list(anchor = anchors, gene = genes)
}
message("\n── Signal over the loop anchors and the anchor gene ──────────────────────")
eu_s  <- get_signals(EU_BW,  "EU-seq")
tot_s <- get_signals(TOT_BW, "Total RNA-seq")

# Loops with no Wildtype signal in either assay are removed (same loops for both)
wt_pos   <- function(s) s$anchor[, "WT"] > 0 & s$gene[, "WT"] > 0
keep_eu  <- wt_pos(eu_s)
keep_tot <- wt_pos(tot_s)
keep     <- keep_eu & keep_tot
message("Loops removed for no Wildtype signal: EU-seq ", sum(!keep_eu),
        ", total RNA-seq ", sum(!keep_tot), ", either ", sum(!keep),
        "   kept: ", sum(keep), " of ", length(keep))
stopifnot(sum(keep) > 0)
sub_s  <- function(s) list(anchor = s$anchor[keep, , drop = FALSE], gene = s$gene[keep, , drop = FALSE])
bn_reg <- bn_reg[keep, ]
eu_s   <- sub_s(eu_s)
tot_s  <- sub_s(tot_s)

sig_log2 <- function(m, label) {
  pc <- pc_for(as.vector(m))
  message("  ", label, ": pseudocount = ", signif(pc, 3))
  log2((m[, "POLR2A"] + pc) / (m[, "WT"] + pc))
}
message("Pseudocounts (1st percentile of the positive values):")
bn_reg$eu_anchor  <- sig_log2(eu_s$anchor,  "EU-seq anchors")
bn_reg$eu_gene    <- sig_log2(eu_s$gene,    "EU-seq anchor gene")
bn_reg$tot_anchor <- sig_log2(tot_s$anchor, "Total RNA-seq anchors")
bn_reg$tot_gene   <- sig_log2(tot_s$gene,   "Total RNA-seq anchor gene")

# ── Variables, loop sets and settings ──────────────────────────────────────────
bin_vars <- c(eu_anchor  = "EU-seq change, anchors",
              eu_gene    = "EU-seq change, whole gene",
              tot_anchor = "Total RNA-seq change, anchors",
              tot_gene   = "Total RNA-seq change, whole gene",
              setNames(paste(MARKS, "change"), paste0(MARKS, "_delta")))
stopifnot(all(names(bin_vars) %in% colnames(bn_reg)),
          all(is.finite(as.matrix(bn_reg[, names(bin_vars)]))),
          all(is.finite(bn_reg$wt_contacts)), all(is.finite(bn_reg$pol2_contacts)),
          all(bn_reg$wt_contacts > 0), all(bn_reg$pol2_contacts >= 0))

loop_sets <- list("Loops with an expressed anchor gene" = bn_reg)
set_tag <- c("Loops with an expressed anchor gene" = "pooled-expr")   # short name used in the SVG file names
stopifnot(all(names(loop_sets) %in% names(set_tag)))
message("\nLoops by set: ", paste(names(loop_sets), sapply(loop_sets, nrow), sep = " = ", collapse = ", "))

N_BOOT         <- 1000L
N_BINS         <- 10L      # deciles
MIN_BIN_N      <- 100L     # smallest bin that is reported
ROLL_HALF_WIDTH <- 0.10    # rolling window = 2 x this share of the loops
ROLL_STEP      <- 0.02
roll_centres   <- seq(ROLL_HALF_WIDTH, 1 - ROLL_HALF_WIDTH, by = ROLL_STEP)

pooled_lfc <- function(wt, pol) log2(sum(pol) / sum(wt))
pooled_boot <- function(wt, pol) {
  n <- length(wt)
  b <- replicate(N_BOOT, { i <- sample.int(n, n, replace = TRUE); pooled_lfc(wt[i], pol[i]) })
  list(est = pooled_lfc(wt, pol), boot = b, n = n, wt = sum(wt), pol = sum(pol))
}
# Pooled log2FC in each rolling window of loops ordered by x (window ends by rank)
roll_est <- function(x, wt, pol, with_x = TRUE) {
  n  <- length(x)
  o  <- order(x)
  cw <- c(0, cumsum(wt[o])); cp <- c(0, cumsum(pol[o]))
  i_lo <- pmax(1L, round((roll_centres - ROLL_HALF_WIDTH) * n) + 1L)
  i_hi <- pmin(n,  round((roll_centres + ROLL_HALF_WIDTH) * n))
  list(est = log2((cp[i_hi + 1L] - cp[i_lo]) / (cw[i_hi + 1L] - cw[i_lo])),
       x   = if (with_x) sapply(seq_along(i_lo), function(k) median(x[o][i_lo[k]:i_hi[k]])) else NULL,
       i_lo = i_lo, i_hi = i_hi)
}

dec_rows <- list(); ext_rows <- list(); roll_rows <- list(); ref_rows <- list()
set.seed(CV_SEED)
for (st in names(loop_sets)) {
  d <- loop_sets[[st]]
  if (nrow(d) < N_BINS * MIN_BIN_N) {
    message("Skipped '", st, "': ", nrow(d), " loops, fewer than ", N_BINS * MIN_BIN_N)
    next
  }
  r_all <- pooled_boot(d$wt_contacts, d$pol2_contacts)
  ref_rows[[st]] <- data.frame(set = st, ref = r_all$est)
  message("\n── ", st, ": ", format(nrow(d), big.mark = ","), " loops, pooled log2FC = ", round(r_all$est, 3),
          "  [", round(quantile(r_all$boot, 0.025), 3), ", ", round(quantile(r_all$boot, 0.975), 3), "]")
  for (v in names(bin_vars)) {
    x <- d[[v]]
    g <- dplyr::ntile(x, N_BINS)
    rs <- lapply(seq_len(N_BINS), function(k) pooled_boot(d$wt_contacts[g == k], d$pol2_contacts[g == k]))
    stopifnot(sum(sapply(rs, `[[`, "n")) == nrow(d),
              isTRUE(all.equal(sum(sapply(rs, `[[`, "wt")), sum(d$wt_contacts))))
    for (k in seq_len(N_BINS)) {
      if (rs[[k]]$n < MIN_BIN_N) next
      dec_rows[[length(dec_rows) + 1]] <- data.frame(
        set = st, variable = bin_vars[[v]], decile = k, n = rs[[k]]$n, median_value = median(x[g == k]),
        sum_wt = rs[[k]]$wt, sum_pol2 = rs[[k]]$pol, pooled_log2FC = rs[[k]]$est,
        lo = unname(quantile(rs[[k]]$boot, 0.025)), hi = unname(quantile(rs[[k]]$boot, 0.975)))
    }
    dd <- rs[[N_BINS]]$boot - rs[[1]]$boot
    ext_rows[[length(ext_rows) + 1]] <- data.frame(
      set = st, variable = bin_vars[[v]], top_minus_bottom_decile = rs[[N_BINS]]$est - rs[[1]]$est,
      lo = unname(quantile(dd, 0.025)), hi = unname(quantile(dd, 0.975)),
      boot_p = min(1, 2 * min(sum(dd <= 0) + 1, sum(dd >= 0) + 1) / (N_BOOT + 1)))
    
    # rolling curve, with the loops resampled and re-ranked in every replicate
    re <- roll_est(x, d$wt_contacts, d$pol2_contacts)
    mid <- which.min(abs(roll_centres - 0.5))
    sel <- rank(x, ties.method = "first") > (re$i_lo[mid] - 1L) & rank(x, ties.method = "first") <= re$i_hi[mid]
    stopifnot(isTRUE(all.equal(re$est[mid], pooled_lfc(d$wt_contacts[sel], d$pol2_contacts[sel]))))
    bt <- replicate(N_BOOT, { i <- sample.int(nrow(d), nrow(d), replace = TRUE)
    roll_est(x[i], d$wt_contacts[i], d$pol2_contacts[i], with_x = FALSE)$est })
    roll_rows[[length(roll_rows) + 1]] <- data.frame(
      set = st, variable = bin_vars[[v]], percentile = roll_centres, value = re$x, pooled_log2FC = re$est,
      lo = apply(bt, 1, quantile, 0.025), hi = apply(bt, 1, quantile, 0.975))
  }
}
stopifnot(length(dec_rows) > 0)
dec_tab  <- do.call(rbind, dec_rows)
ext_tab  <- do.call(rbind, ext_rows)
roll_tab <- do.call(rbind, roll_rows)
ref_tab  <- do.call(rbind, ref_rows)

options(max.print = 1e5)
message("\nPooled log2(sum degraded / sum Wildtype contacts) by decile")
print(dec_tab, row.names = FALSE, digits = 3, width = 250)
message("\nTop decile minus bottom decile (pooled log2FC), bootstrap 95% interval and p")
print(ext_tab, row.names = FALSE, digits = 3, width = 250)

# ── Plots ──────────────────────────────────────────────────────────────────────
library(ggplot2)
col_main <- "#3b6b8f"
for (st in unique(dec_tab$set)) {
  dd <- dec_tab[dec_tab$set == st, ]
  rr <- roll_tab[roll_tab$set == st, ]
  dd$variable <- factor(dd$variable, levels = unname(bin_vars))
  rr$variable <- factor(rr$variable, levels = unname(bin_vars))
  n_lab <- paste0("N = ", format(nrow(loop_sets[[st]]), big.mark = ","), " loops; points = deciles, ",
                  "line and band = rolling window (", 100 * 2 * ROLL_HALF_WIDTH, "% of loops); ",
                  "dashed line = all loops in this set pooled")
  p <- ggplot() +
    geom_hline(yintercept = ref_tab$ref[ref_tab$set == st], linetype = "dashed", colour = "grey50") +
    geom_ribbon(data = rr, aes(value, ymin = lo, ymax = hi), fill = col_main, alpha = 0.18) +
    geom_line(data = rr, aes(value, pooled_log2FC), colour = col_main, linewidth = 0.8) +
    geom_linerange(data = dd, aes(median_value, ymin = lo, ymax = hi), colour = "grey25", linewidth = 0.5) +
    geom_point(data = dd, aes(median_value, pooled_log2FC), colour = "grey25", size = 1.8) +
    facet_wrap(~ variable, scales = "free_x", ncol = 4) +
    labs(title = paste0("Pooled loop change across each variable: ", st), subtitle = n_lab,
         x = "Value of the variable (log2)", y = "Pooled log2(sum degraded / sum Wildtype contacts)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
  print(p)
  save_plot_svg(p, set_tag[[st]], "by-var", w = 16, h = 8)
  
  # All variables on one plot: x = percentile of each variable (rolling window centres,
  # so the curves span ROLL_HALF_WIDTH to 1 - ROLL_HALF_WIDTH)
  var_cols <- setNames(c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7"),
                       levels(rr$variable))
  stopifnot(length(var_cols) == 7)
  p2 <- ggplot() +
    geom_hline(yintercept = ref_tab$ref[ref_tab$set == st], linetype = "dashed", colour = "grey50") +
    geom_ribbon(data = rr, aes(100 * percentile, ymin = lo, ymax = hi, fill = variable), alpha = 0.07) +
    geom_line(data = rr, aes(100 * percentile, pooled_log2FC, colour = variable), linewidth = 0.9) +
    scale_colour_manual(values = var_cols, name = NULL) +
    scale_fill_manual(values = var_cols, name = NULL) +
    labs(title = paste0("Pooled loop change by percentile of each variable: ", st),
         subtitle = paste0("N = ", format(nrow(loop_sets[[st]]), big.mark = ","), " loops; lines = rolling window (",
                           100 * 2 * ROLL_HALF_WIDTH, "% of loops), bands = bootstrap 95% interval; ",
                           "dashed line = all loops pooled"),
         x = "Percentile of the variable (0 = largest decrease)",
         y = "Pooled log2(sum degraded / sum Wildtype contacts)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9), legend.position = "bottom")
  print(p2)
  save_plot_svg(p2, set_tag[[st]], "all-vars", w = 11, h = 7)
}

# ══════════════════════════════════════════════════════════════════════════════
# EU-seq LEVEL as the variable (instead of the change)
# Same loops, same pooled log2(sum degraded / sum Wildtype contacts), same deciles and
# rolling curve. The x variable is how much nascent transcription there is, not how much
# it changed:
#   POLR2A-degraded cells: log2(EU-seq signal + pc)   = what is left after degradation
#   Wildtype cells:        log2(EU-seq signal + pc)   = where the gene started
# (same pseudocount as the change variables, so change = degraded level - Wildtype level).
# Shown over the loop anchors and over the whole anchor gene. Separate plot: *_by-var-level.svg
# ══════════════════════════════════════════════════════════════════════════════
message("\n── EU-seq level as the variable ─────────────────────────────────────────")
pc_eu_anchor <- pc_for(as.vector(eu_s$anchor))
pc_eu_gene   <- pc_for(as.vector(eu_s$gene))
bn_reg$eu_anchor_deg <- log2(eu_s$anchor[, "POLR2A"] + pc_eu_anchor)
bn_reg$eu_anchor_wt  <- log2(eu_s$anchor[, "WT"]     + pc_eu_anchor)
bn_reg$eu_gene_deg   <- log2(eu_s$gene[, "POLR2A"]   + pc_eu_gene)
bn_reg$eu_gene_wt    <- log2(eu_s$gene[, "WT"]       + pc_eu_gene)

lvl_vars <- c(eu_anchor_deg = "EU-seq, POLR2A-degraded cells, anchors",
              eu_gene_deg   = "EU-seq, POLR2A-degraded cells, whole gene",
              eu_anchor_wt  = "EU-seq, Wildtype cells, anchors",
              eu_gene_wt    = "EU-seq, Wildtype cells, whole gene")

# Checks on the live data: level - level = the change used above; no missing values
stopifnot(all(names(lvl_vars) %in% colnames(bn_reg)),
          all(is.finite(as.matrix(bn_reg[, names(lvl_vars)]))),
          isTRUE(all.equal(bn_reg$eu_anchor_deg - bn_reg$eu_anchor_wt, bn_reg$eu_anchor)),
          isTRUE(all.equal(bn_reg$eu_gene_deg   - bn_reg$eu_gene_wt,   bn_reg$eu_gene)))
message("Check passed: degraded level - Wildtype level = EU-seq change (anchors and whole gene)")
message("Loops with no EU-seq signal in POLR2A-degraded cells (tied at the lowest value): anchors ",
        sum(eu_s$anchor[, "POLR2A"] == 0), ", whole gene ", sum(eu_s$gene[, "POLR2A"] == 0),
        " of ", nrow(bn_reg))
message("Correlation of the degraded level with the Wildtype level (anchors / whole gene): ",
        round(cor(bn_reg$eu_anchor_deg, bn_reg$eu_anchor_wt), 3), " / ",
        round(cor(bn_reg$eu_gene_deg, bn_reg$eu_gene_wt), 3))
message("Correlation of the degraded level with the change (anchors / whole gene): ",
        round(cor(bn_reg$eu_anchor_deg, bn_reg$eu_anchor), 3), " / ",
        round(cor(bn_reg$eu_gene_deg, bn_reg$eu_gene), 3))

lvl_sets <- list("Loops with an expressed anchor gene" = bn_reg)
stopifnot(all(names(lvl_sets) %in% names(set_tag)))

lvl_dec <- list(); lvl_ext <- list(); lvl_roll <- list()
set.seed(CV_SEED)
for (st in names(lvl_sets)) {
  d <- lvl_sets[[st]]
  stopifnot(nrow(d) >= N_BINS * MIN_BIN_N)
  for (v in names(lvl_vars)) {
    x  <- d[[v]]
    g  <- dplyr::ntile(x, N_BINS)
    rs <- lapply(seq_len(N_BINS), function(k) pooled_boot(d$wt_contacts[g == k], d$pol2_contacts[g == k]))
    stopifnot(sum(sapply(rs, `[[`, "n")) == nrow(d),
              isTRUE(all.equal(sum(sapply(rs, `[[`, "wt")), sum(d$wt_contacts))))
    for (k in seq_len(N_BINS)) {
      if (rs[[k]]$n < MIN_BIN_N) next
      lvl_dec[[length(lvl_dec) + 1]] <- data.frame(
        set = st, variable = lvl_vars[[v]], decile = k, n = rs[[k]]$n, median_value = median(x[g == k]),
        sum_wt = rs[[k]]$wt, sum_pol2 = rs[[k]]$pol, pooled_log2FC = rs[[k]]$est,
        lo = unname(quantile(rs[[k]]$boot, 0.025)), hi = unname(quantile(rs[[k]]$boot, 0.975)))
    }
    dd <- rs[[N_BINS]]$boot - rs[[1]]$boot
    lvl_ext[[length(lvl_ext) + 1]] <- data.frame(
      set = st, variable = lvl_vars[[v]], top_minus_bottom_decile = rs[[N_BINS]]$est - rs[[1]]$est,
      lo = unname(quantile(dd, 0.025)), hi = unname(quantile(dd, 0.975)),
      boot_p = min(1, 2 * min(sum(dd <= 0) + 1, sum(dd >= 0) + 1) / (N_BOOT + 1)))
    
    # rolling curve, loops resampled and re-ranked in every replicate (as above)
    re  <- roll_est(x, d$wt_contacts, d$pol2_contacts)
    mid <- which.min(abs(roll_centres - 0.5))
    sel <- rank(x, ties.method = "first") > (re$i_lo[mid] - 1L) & rank(x, ties.method = "first") <= re$i_hi[mid]
    stopifnot(isTRUE(all.equal(re$est[mid], pooled_lfc(d$wt_contacts[sel], d$pol2_contacts[sel]))))
    bt <- replicate(N_BOOT, { i <- sample.int(nrow(d), nrow(d), replace = TRUE)
    roll_est(x[i], d$wt_contacts[i], d$pol2_contacts[i], with_x = FALSE)$est })
    lvl_roll[[length(lvl_roll) + 1]] <- data.frame(
      set = st, variable = lvl_vars[[v]], percentile = roll_centres, value = re$x, pooled_log2FC = re$est,
      lo = apply(bt, 1, quantile, 0.025), hi = apply(bt, 1, quantile, 0.975))
  }
}
lvl_dec_tab  <- do.call(rbind, lvl_dec)
lvl_ext_tab  <- do.call(rbind, lvl_ext)
lvl_roll_tab <- do.call(rbind, lvl_roll)

message("\nPooled log2(sum degraded / sum Wildtype contacts) by decile of the EU-seq level")
print(lvl_dec_tab, row.names = FALSE, digits = 3, width = 250)
message("\nTop decile minus bottom decile of the EU-seq level (pooled log2FC), bootstrap 95% interval and p")
print(lvl_ext_tab, row.names = FALSE, digits = 3, width = 250)

for (st in unique(lvl_dec_tab$set)) {
  dd <- lvl_dec_tab[lvl_dec_tab$set == st, ]
  rr <- lvl_roll_tab[lvl_roll_tab$set == st, ]
  dd$variable <- factor(dd$variable, levels = unname(lvl_vars))
  rr$variable <- factor(rr$variable, levels = unname(lvl_vars))
  n_lab <- paste0("N = ", format(nrow(lvl_sets[[st]]), big.mark = ","), " loops; points = deciles, ",
                  "line and band = rolling window (", 100 * 2 * ROLL_HALF_WIDTH, "% of loops); ",
                  "dashed line = all loops in this set pooled")
  p3 <- ggplot() +
    geom_hline(yintercept = ref_tab$ref[ref_tab$set == st], linetype = "dashed", colour = "grey50") +
    geom_ribbon(data = rr, aes(value, ymin = lo, ymax = hi), fill = col_main, alpha = 0.18) +
    geom_line(data = rr, aes(value, pooled_log2FC), colour = col_main, linewidth = 0.8) +
    geom_linerange(data = dd, aes(median_value, ymin = lo, ymax = hi), colour = "grey25", linewidth = 0.5) +
    geom_point(data = dd, aes(median_value, pooled_log2FC), colour = "grey25", size = 1.8) +
    facet_wrap(~ variable, scales = "free_x", ncol = 2) +
    labs(title = paste0("Pooled loop change across the EU-seq level: ", st), subtitle = n_lab,
         x = "EU-seq signal (log2)", y = "Pooled log2(sum degraded / sum Wildtype contacts)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
  print(p3)
  save_plot_svg(p3, set_tag[[st]], "by-var-level", w = 11, h = 8)
}


# ══════════════════════════════════════════════════════════════════════════════
# GAINED loops: loops called in POLR2A-degraded cells
# Mirror of the Wildtype-loop set above: loops called in POLR2A-degraded cells (shared and
# degraded-only; degraded contacts >= MIN_POL2_CONTACTS) whose anchor overlaps a gene
# expressed in Wildtype cells (TPM > 0). Same variables (change in EU-seq, total RNA-seq,
# H3K27ac, H3K4me3, H3K9me3; same pseudocounts as the Wildtype-loop set), same pooled
# log2(sum degraded / sum Wildtype contacts), deciles and rolling curve.
# Selecting loops on their degraded contacts pushes the pooled ratio up, as selecting on
# Wildtype contacts pushes it down above, so the slope across the variable is the readout,
# not the level. Separate plot: *_pooled-gained-expr_by-var.svg
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Loops called in POLR2A-degraded cells ────────────────────────────────")
gn_reg <- features %>%
  dplyr::filter(source %in% c("shared", "pol2_only"), pol2_contacts >= MIN_POL2_CONTACTS) %>%
  dplyr::select(-dplyr::any_of(c("anchor_gene_tpm", "expr_bin", "anchor_gene_id"))) %>%
  link_loops_to_genes(gene_gr, tol_bp = 0L)
n_gn_called <- nrow(gn_reg)
gn_reg <- gn_reg %>% dplyr::filter(!is.na(anchor_gene_id))
n_gn_gene <- nrow(gn_reg)
gn_reg <- dplyr::filter(gn_reg, anchor_gene_tpm > 0)
message("Degraded-cell loops with degraded contacts >= ", MIN_POL2_CONTACTS, ": ", n_gn_called,
        " (shared ", sum(features$source == "shared" & features$pol2_contacts >= MIN_POL2_CONTACTS),
        ", degraded-only ", sum(features$source == "pol2_only" & features$pol2_contacts >= MIN_POL2_CONTACTS),
        ")   with an anchor overlapping a gene: ", n_gn_gene, "   anchor gene expressed (TPM > 0): ", nrow(gn_reg))
stopifnot(nrow(gn_reg) > 0, !anyNA(match(gn_reg$anchor_gene_id, gene_gr$gene_id)))

gn_ug    <- unique(gn_reg$anchor_gene_id)
gn_ug_gr <- gene_gr[match(gn_ug, gene_gr$gene_id)]
get_signals_gn <- function(BW, label) {
  message("  ", label)
  anchors <- sapply(names(BW), function(nm) {
    message("    ", nm, " : ", BW[[nm]])
    a1 <- anchor_signal(BW[[nm]], gn_reg$chr1, gn_reg$start1, gn_reg$end1)
    a2 <- anchor_signal(BW[[nm]], gn_reg$chr2, gn_reg$start2, gn_reg$end2)
    (a1 + a2) / 2
  })
  genes <- sapply(names(BW), function(nm) {
    g <- anchor_signal(BW[[nm]], as.character(seqnames(gn_ug_gr)), start(gn_ug_gr) - 1L, end(gn_ug_gr))
    g[match(gn_reg$anchor_gene_id, gn_ug)]
  })
  stopifnot(nrow(anchors) == nrow(gn_reg), nrow(genes) == nrow(gn_reg),
            all(is.finite(anchors)), all(is.finite(genes)))
  if (any(anchors < 0) || any(genes < 0)) {
    stop("Negative ", label, " signal in the bigWigs (anchors ", sum(anchors < 0), ", genes ",
         sum(genes < 0), "); the hg38 chrM-scaled bigWigs are not negative, so the paths point at the wrong files")
  }
  list(anchor = anchors, gene = genes)
}
message("\n── Signal over the loop anchors and the anchor gene (degraded-cell loops) ──")
gn_eu  <- get_signals_gn(EU_BW,  "EU-seq")
gn_tot <- get_signals_gn(TOT_BW, "Total RNA-seq")

# Loops with no Wildtype signal in either assay are removed (same rule as for the Wildtype loops)
gn_keep_eu  <- wt_pos(gn_eu)
gn_keep_tot <- wt_pos(gn_tot)
gn_keep     <- gn_keep_eu & gn_keep_tot
message("Loops removed for no Wildtype signal: EU-seq ", sum(!gn_keep_eu),
        ", total RNA-seq ", sum(!gn_keep_tot), ", either ", sum(!gn_keep),
        "   kept: ", sum(gn_keep), " of ", length(gn_keep))
stopifnot(sum(gn_keep) > 0)
gn_reg <- gn_reg[gn_keep, ]
gn_eu  <- list(anchor = gn_eu$anchor[gn_keep, , drop = FALSE],  gene = gn_eu$gene[gn_keep, , drop = FALSE])
gn_tot <- list(anchor = gn_tot$anchor[gn_keep, , drop = FALSE], gene = gn_tot$gene[gn_keep, , drop = FALSE])

# Same pseudocounts as the Wildtype-loop set (1st percentile of the positive values there),
# so the variables are on the same scale in both plots
pc_ref <- c(eu_anchor  = pc_for(as.vector(eu_s$anchor)),  eu_gene  = pc_for(as.vector(eu_s$gene)),
            tot_anchor = pc_for(as.vector(tot_s$anchor)), tot_gene = pc_for(as.vector(tot_s$gene)))
message("Pseudocounts taken from the Wildtype-loop set: ",
        paste(names(pc_ref), signif(pc_ref, 3), collapse = ", "))
sig_log2_pc <- function(m, pc) log2((m[, "POLR2A"] + pc) / (m[, "WT"] + pc))
gn_reg$eu_anchor  <- sig_log2_pc(gn_eu$anchor,  pc_ref[["eu_anchor"]])
gn_reg$eu_gene    <- sig_log2_pc(gn_eu$gene,    pc_ref[["eu_gene"]])
gn_reg$tot_anchor <- sig_log2_pc(gn_tot$anchor, pc_ref[["tot_anchor"]])
gn_reg$tot_gene   <- sig_log2_pc(gn_tot$gene,   pc_ref[["tot_gene"]])

stopifnot(all(names(bin_vars) %in% colnames(gn_reg)),
          all(is.finite(as.matrix(gn_reg[, names(bin_vars)]))),
          all(is.finite(gn_reg$wt_contacts)), all(is.finite(gn_reg$pol2_contacts)),
          all(gn_reg$wt_contacts >= 0), all(gn_reg$pol2_contacts >= MIN_POL2_CONTACTS),
          sum(gn_reg$wt_contacts) > 0)
message("Loops in the set by source: shared ", sum(gn_reg$source == "shared"),
        ", degraded-only ", sum(gn_reg$source == "pol2_only"),
        "   Gained class (log2FC > ", LOG2FC_THRESHOLD, "): ", sum(gn_reg$loop_class == "Gained"))

# Pooled curves for one loop set (deciles, top minus bottom decile, rolling curve), as above
pooled_curves <- function(d, vars, st) {
  dec <- list(); ext <- list(); rol <- list()
  for (v in names(vars)) {
    x  <- d[[v]]
    g  <- dplyr::ntile(x, N_BINS)
    rs <- lapply(seq_len(N_BINS), function(k) pooled_boot(d$wt_contacts[g == k], d$pol2_contacts[g == k]))
    stopifnot(sum(sapply(rs, `[[`, "n")) == nrow(d),
              isTRUE(all.equal(sum(sapply(rs, `[[`, "wt")), sum(d$wt_contacts))))
    for (k in seq_len(N_BINS)) {
      if (rs[[k]]$n < MIN_BIN_N) next
      dec[[length(dec) + 1]] <- data.frame(
        set = st, variable = vars[[v]], decile = k, n = rs[[k]]$n, median_value = median(x[g == k]),
        sum_wt = rs[[k]]$wt, sum_pol2 = rs[[k]]$pol, pooled_log2FC = rs[[k]]$est,
        lo = unname(quantile(rs[[k]]$boot, 0.025)), hi = unname(quantile(rs[[k]]$boot, 0.975)))
    }
    dd <- rs[[N_BINS]]$boot - rs[[1]]$boot
    ext[[length(ext) + 1]] <- data.frame(
      set = st, variable = vars[[v]], top_minus_bottom_decile = rs[[N_BINS]]$est - rs[[1]]$est,
      lo = unname(quantile(dd, 0.025)), hi = unname(quantile(dd, 0.975)),
      boot_p = min(1, 2 * min(sum(dd <= 0) + 1, sum(dd >= 0) + 1) / (N_BOOT + 1)))
    re  <- roll_est(x, d$wt_contacts, d$pol2_contacts)
    mid <- which.min(abs(roll_centres - 0.5))
    sel <- rank(x, ties.method = "first") > (re$i_lo[mid] - 1L) & rank(x, ties.method = "first") <= re$i_hi[mid]
    stopifnot(isTRUE(all.equal(re$est[mid], pooled_lfc(d$wt_contacts[sel], d$pol2_contacts[sel]))))
    bt <- replicate(N_BOOT, { i <- sample.int(nrow(d), nrow(d), replace = TRUE)
    roll_est(x[i], d$wt_contacts[i], d$pol2_contacts[i], with_x = FALSE)$est })
    rol[[length(rol) + 1]] <- data.frame(
      set = st, variable = vars[[v]], percentile = roll_centres, value = re$x, pooled_log2FC = re$est,
      lo = apply(bt, 1, quantile, 0.025), hi = apply(bt, 1, quantile, 0.975))
  }
  list(dec = do.call(rbind, dec), ext = do.call(rbind, ext), roll = do.call(rbind, rol))
}

gn_name <- "Degraded-cell loops with an expressed anchor gene"
gn_tag  <- c("Degraded-cell loops with an expressed anchor gene" = "pooled-gained-expr")
stopifnot(nrow(gn_reg) >= N_BINS * MIN_BIN_N)
set.seed(CV_SEED)
gn_all <- pooled_boot(gn_reg$wt_contacts, gn_reg$pol2_contacts)
message("\n── ", gn_name, ": ", format(nrow(gn_reg), big.mark = ","), " loops, pooled log2FC = ",
        round(gn_all$est, 3), "  [", round(quantile(gn_all$boot, 0.025), 3), ", ",
        round(quantile(gn_all$boot, 0.975), 3), "]")
gn_res <- pooled_curves(gn_reg, bin_vars, gn_name)

message("\nPooled log2(sum degraded / sum Wildtype contacts) by decile, degraded-cell loops")
print(gn_res$dec, row.names = FALSE, digits = 3, width = 250)
message("\nTop decile minus bottom decile (pooled log2FC), degraded-cell loops, bootstrap 95% interval and p")
print(gn_res$ext, row.names = FALSE, digits = 3, width = 250)

dd <- gn_res$dec
rr <- gn_res$roll
dd$variable <- factor(dd$variable, levels = unname(bin_vars))
rr$variable <- factor(rr$variable, levels = unname(bin_vars))
n_lab <- paste0("N = ", format(nrow(gn_reg), big.mark = ","), " loops; points = deciles, ",
                "line and band = rolling window (", 100 * 2 * ROLL_HALF_WIDTH, "% of loops); ",
                "dashed line = all loops in this set pooled")
p4 <- ggplot() +
  geom_hline(yintercept = gn_all$est, linetype = "dashed", colour = "grey50") +
  geom_ribbon(data = rr, aes(value, ymin = lo, ymax = hi), fill = col_main, alpha = 0.18) +
  geom_line(data = rr, aes(value, pooled_log2FC), colour = col_main, linewidth = 0.8) +
  geom_linerange(data = dd, aes(median_value, ymin = lo, ymax = hi), colour = "grey25", linewidth = 0.5) +
  geom_point(data = dd, aes(median_value, pooled_log2FC), colour = "grey25", size = 1.8) +
  facet_wrap(~ variable, scales = "free_x", ncol = 4) +
  labs(title = paste0("Pooled loop change across each variable: ", gn_name), subtitle = n_lab,
       x = "Value of the variable (log2)", y = "Pooled log2(sum degraded / sum Wildtype contacts)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
print(p4)
save_plot_svg(p4, gn_tag[[gn_name]], "by-var", w = 16, h = 8)


# ══════════════════════════════════════════════════════════════════════════════
# Most loss vs least loss loops (top and bottom of loop log2FC)
# Wildtype loops with an expressed anchor gene (the set of plot A) are ranked by their own log2FC
# (log2((degraded + 1) / (Wildtype + 1) contacts)) and split into n_grp equal groups; only the lowest
# group (most loss) and the highest group (least loss) are kept. For each variable:
#   line plot: pooled log2(sum degraded / sum Wildtype contacts) across deciles (points) and a rolling
#     curve (line and band) of the variable, shown as % change in the contact ratio relative to the
#     group's own pooled ratio, 100 * (2^(pooled log2FC - pooled log2FC of the group) - 1)
# Needs the columns named in `vars` in bn_reg (eu_gene, tot_gene = change in EU-seq and total RNA-seq over
# the whole anchor gene, log2(degraded / Wildtype), negative = decrease; eu_gene_wt = Wildtype EU-seq signal)
# and pooled_curves() above.
# ══════════════════════════════════════════════════════════════════════════════
loss_groups <- function(n_grp, grp_word, vars, tag, type, x_lab) {
  lab <- c(paste0("Most loss (lowest ", grp_word, ")"), paste0("Least loss (highest ", grp_word, ")"))
  col <- setNames(c("#D55E00", "#0072B2"), lab)
  stopifnot(all(c(names(vars), "log2fc") %in% colnames(bn_reg)),
            all(is.finite(as.matrix(bn_reg[, c(names(vars), "log2fc")]))))
  g <- dplyr::ntile(bn_reg$log2fc, n_grp)
  stopifnot(!anyNA(g), all(table(g) >= N_BINS * MIN_BIN_N),
            max(bn_reg$log2fc[g == 1]) <= min(bn_reg$log2fc[g == n_grp]))
  message("\n── ", tag, " (", paste(unname(vars), collapse = " | "), ") ─────────────────────────")
  
  set.seed(CV_SEED)
  dec <- list(); rol <- list(); ref <- setNames(numeric(2), lab)
  for (i in 1:2) {
    d     <- bn_reg[g == c(1, n_grp)[i], ]
    d_all <- pooled_boot(d$wt_contacts, d$pol2_contacts)
    ref[i] <- d_all$est
    message("\n", lab[i], ": ", format(nrow(d), big.mark = ","), " loops, loop log2FC from ",
            round(min(d$log2fc), 2), " to ", round(max(d$log2fc), 2), " (median ", round(median(d$log2fc), 2),
            "), pooled log2FC = ", round(d_all$est, 3), "  [", round(quantile(d_all$boot, 0.025), 3), ", ",
            round(quantile(d_all$boot, 0.975), 3), "]")
    res <- pooled_curves(d, vars, lab[i])
    message("Top decile minus bottom decile of the variable (pooled log2FC), bootstrap 95% interval and p")
    print(res$ext, row.names = FALSE, digits = 3, width = 250)
    dec[[i]] <- res$dec
    rol[[i]] <- res$roll
  }
  dd <- do.call(rbind, dec)
  rr <- do.call(rbind, rol)
  dd$variable <- factor(dd$variable, levels = unname(vars))
  rr$variable <- factor(rr$variable, levels = unname(vars))
  dd$group    <- factor(dd$set, levels = lab)
  rr$group    <- factor(rr$set, levels = lab)
  stopifnot(!anyNA(dd$group), !anyNA(rr$group))
  to_pct <- function(y, grp) 100 * (2^(y - ref[as.character(grp)]) - 1)
  dd$pct <- to_pct(dd$pooled_log2FC, dd$group); dd$pct_lo <- to_pct(dd$lo, dd$group); dd$pct_hi <- to_pct(dd$hi, dd$group)
  rr$pct <- to_pct(rr$pooled_log2FC, rr$group); rr$pct_lo <- to_pct(rr$lo, rr$group); rr$pct_hi <- to_pct(rr$hi, rr$group)
  stopifnot(all(is.finite(c(dd$pct, rr$pct, rr$pct_lo, rr$pct_hi))))
  message("\nRange of the rolling curve (% change relative to the group's pooled ratio):")
  print(aggregate(pct ~ group + variable, data = rr, FUN = function(z) round(range(z), 1)))
  
  th <- theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 9),
          legend.position = "bottom")
  n_lab <- paste0("N = ", format(sum(table(g)[c(1, n_grp)]), big.mark = ","), " loops in the lowest and highest ",
                  grp_word, " of loop log2FC (", format(min(table(g)), big.mark = ","), " to ",
                  format(max(table(g)), big.mark = ","), " per ", grp_word, ")")
  p <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_ribbon(data = rr, aes(value, ymin = pct_lo, ymax = pct_hi, fill = group), alpha = 0.08) +
    geom_line(data = rr, aes(value, pct, colour = group), linewidth = 0.9) +
    geom_point(data = dd, aes(median_value, pct, colour = group), size = 1.6, alpha = 0.8) +
    facet_wrap(~ variable, scales = "free_x", ncol = length(vars)) +
    scale_colour_manual(values = col, name = NULL) +
    scale_fill_manual(values = col, name = NULL) +
    labs(title = paste0("Pooled loop change, most vs least loss (", grp_word, "s)"),
         subtitle = paste0(n_lab, "\npoints = deciles, line = rolling window (", 100 * 2 * ROLL_HALF_WIDTH,
                           "% of the loops in the group), band = bootstrap 95% interval",
                           "\ny = change in contact ratio relative to the group's own pooled ratio"),
         x = x_lab, y = "% contact change") +
    th
  print(p)
  save_plot_svg(p, paste0("pooled-expr-loss-", tag), type, w = if (length(vars) == 1) 6 else 11, h = 6)
  invisible(NULL)
}

chg_vars <- c(eu_gene  = "EU-seq (nascent) change, whole gene",
              tot_gene = "Total RNA-seq change, whole gene")
chg_x    <- "Change after POLR2A degradation (log2 degraded / Wildtype; negative = decrease)"
loss_groups(3, "third",   chg_vars, "thirds-top-bottom",   "by-var-change-gene", chg_x)
loss_groups(4, "quarter", chg_vars, "quarters-top-bottom", "by-var-change-gene", chg_x)
loss_groups(3, "third", c(eu_gene_wt = "EU-seq, Wildtype cells, whole gene"), "thirds-top-bottom",
            "by-var-wt-level-gene", "Wildtype EU-seq signal (log2)")

# ══════════════════════════════════════════════════════════════════════════════
# Pooled loop change, all loops
# All loops = the Wildtype-called loops (Wildtype contacts >= MIN_WT_CONTACTS) and the loops called in
# POLR2A-degraded cells (degraded contacts >= MIN_POL2_CONTACTS), each with an anchor gene expressed in
# Wildtype cells (the sets bn_reg and gn_reg above); a loop called in both is counted once.
# Variables: change in EU-seq (nascent) and total RNA-seq over the whole anchor gene, and the change in
# H3K27ac, H3K4me3 and H3K9me3, each log2(POLR2A-degraded / Wildtype); negative = decrease.
# Two plots:
#   (1) all variables on one plot, x = percentile of the variable among all loops (0 = largest decrease),
#       one coloured line per variable; y = % change in the contact ratio relative to the pooled ratio of
#       all loops, 100 * (2^(pooled log2FC - pooled log2FC of all loops) - 1):   *_pooled-loops-all_all-vars-pct.svg
#   (2) one panel per variable, x = the change, y = pooled log2(sum degraded / sum Wildtype contacts):
#                                                                                *_pooled-loops-all_by-var-log2.svg
# Needs bn_reg, gn_reg and pooled_curves() above.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Pooled loop change, all loops ──────────────────────────────────────────")
all_vars <- c(eu_gene       = "EU-seq (nascent) change, whole gene",
              tot_gene      = "Total RNA-seq change, whole gene",
              H3K27ac_delta = "H3K27ac change",
              H3K4me3_delta = "H3K4me3 change",
              H3K9me3_delta = "H3K9me3 change")
all_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "source", "wt_contacts", "pol2_contacts",
              names(all_vars))
stopifnot(all(all_cols %in% colnames(bn_reg)), all(all_cols %in% colnames(gn_reg)))

loop_key <- function(x) paste(x$chr1, x$start1, x$end1, x$chr2, x$start2, x$end2, sep = ":")
all_raw <- rbind(as.data.frame(bn_reg[, all_cols]), as.data.frame(gn_reg[, all_cols]))
all_key <- loop_key(all_raw)
# a loop called in both is in both tables: its values must be identical in the two
for (v in c("source", "wt_contacts", "pol2_contacts", names(all_vars))) {
  same <- tapply(all_raw[[v]], all_key, function(z) length(unique(round(if (is.numeric(z)) z else as.numeric(factor(z)), 8))) == 1)
  stopifnot(all(same))
}
all_loops <- all_raw[!duplicated(all_key), ]
message("Check passed: loops in both tables have identical values (", sum(duplicated(all_key)), " loops counted once)")
message("All loops: ", format(nrow(all_loops), big.mark = ","), "  (Wildtype-called ", format(nrow(bn_reg), big.mark = ","),
        ", POLR2A-degraded-called ", format(nrow(gn_reg), big.mark = ","), ")")
stopifnot(nrow(all_loops) >= N_BINS * MIN_BIN_N)

set.seed(CV_SEED)
all_pool <- pooled_boot(all_loops$wt_contacts, all_loops$pol2_contacts)
message("Pooled log2FC of all loops = ", round(all_pool$est, 3), "  [", round(quantile(all_pool$boot, 0.025), 3),
        ", ", round(quantile(all_pool$boot, 0.975), 3), "]")
all_res <- pooled_curves(all_loops, all_vars, "All loops")
message("\nTop decile minus bottom decile of the variable (pooled log2FC), bootstrap 95% interval and p")
print(all_res$ext, row.names = FALSE, digits = 3, width = 250)
dd <- all_res$dec
rr <- all_res$roll
dd$variable <- factor(dd$variable, levels = unname(all_vars))
rr$variable <- factor(rr$variable, levels = unname(all_vars))

x_chg   <- "Change after POLR2A degradation (log2 fold change)"
x_pct   <- "Percentile of change (0 = largest decrease)"
var_col <- setNames(c("#E69F00", "#009E73", "#D55E00", "#0072B2", "#CC79A7"), unname(all_vars))
n_all   <- format(nrow(all_loops), big.mark = ",")

# (1) all variables on one plot, y = % change relative to the pooled ratio of all loops
to_pct <- function(y) 100 * (2^(y - all_pool$est) - 1)
rr$pct <- to_pct(rr$pooled_log2FC); rr$pct_lo <- to_pct(rr$lo); rr$pct_hi <- to_pct(rr$hi)
stopifnot(all(is.finite(c(rr$pct, rr$pct_lo, rr$pct_hi))))
message("\nRange of the rolling curve (% change relative to the pooled ratio of all loops):")
print(aggregate(pct ~ variable, data = rr, FUN = function(z) round(range(z), 1)))
p_a <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_ribbon(data = rr, aes(100 * percentile, ymin = pct_lo, ymax = pct_hi, fill = variable), alpha = 0.07) +
  geom_line(data = rr, aes(100 * percentile, pct, colour = variable), linewidth = 0.9) +
  scale_colour_manual(values = var_col, name = NULL) +
  scale_fill_manual(values = var_col, name = NULL) +
  labs(title = "Pooled loop change by percentile of each variable: all loops",
       subtitle = paste0("N = ", n_all, " loops; line = rolling window (", 100 * 2 * ROLL_HALF_WIDTH,
                         "% of loops), band = bootstrap 95% interval",
                         "\ny = change in contact ratio relative to the pooled ratio of all loops"),
       x = x_pct, y = "% contact change") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9), legend.position = "bottom")
print(p_a)
save_plot_svg(p_a, "pooled-loops-all", "all-vars-pct", w = 11, h = 7)

# (2) one panel per variable, bare log2 on y
p_s <- ggplot() +
  geom_hline(yintercept = all_pool$est, linetype = "dashed", colour = "grey50") +
  geom_ribbon(data = rr, aes(value, ymin = lo, ymax = hi), fill = col_main, alpha = 0.18) +
  geom_line(data = rr, aes(value, pooled_log2FC), colour = col_main, linewidth = 0.8) +
  geom_linerange(data = dd, aes(median_value, ymin = lo, ymax = hi), colour = "grey25", linewidth = 0.5) +
  geom_point(data = dd, aes(median_value, pooled_log2FC), colour = "grey25", size = 1.8) +
  facet_wrap(~ variable, scales = "free_x", ncol = 3) +
  labs(title = "Pooled loop change across each variable: all loops",
       subtitle = paste0("N = ", n_all, " loops; points = deciles, line and band = rolling window (",
                         100 * 2 * ROLL_HALF_WIDTH, "% of loops); dashed line = all loops pooled; ",
                         "x: negative = decrease"),
       x = x_chg, y = "Contact change (log2)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
print(p_s)
save_plot_svg(p_s, "pooled-loops-all", "by-var-log2", w = 12, h = 8)


# ══════════════════════════════════════════════════════════════════════════════
# Pooled loop change versus transcription level in Wildtype and POLR2A-degraded cells, all loops
# The Wildtype-level plot above plus the same curves with x = transcription level of the anchor gene in
# POLR2A-degraded cells (whole gene, log2(signal + pseudocount), the same pseudocounts as the change
# variables), each as the percentile among all loops (0 = lowest level). One colour per assay (nascent, total),
# solid line = Wildtype cells, dashed line = POLR2A-degraded cells. y as above.
#   *_pooled-loops-all_wt-deg-level-pct.svg
# Needs the objects of the Wildtype-level block above.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Pooled loop change versus Wildtype and POLR2A-degraded transcription level, all loops ─")
dg_vars <- c(eu_gene_deg  = "Nascent (EU-seq), POLR2A-degraded cells",
             tot_gene_deg = "Total RNA-seq, POLR2A-degraded cells")

deg_level <- function(m, pc) log2(m[, "POLR2A"] + pc)
stopifnot(  # degraded level - Wildtype level = the change columns, in both loop sets
  isTRUE(all.equal(unname(deg_level(eu_s$gene,  pc_ref[["eu_gene"]])  - wt_level(eu_s$gene,  pc_ref[["eu_gene"]])),  bn_reg$eu_gene)),
  isTRUE(all.equal(unname(deg_level(tot_s$gene, pc_ref[["tot_gene"]]) - wt_level(tot_s$gene, pc_ref[["tot_gene"]])), bn_reg$tot_gene)),
  isTRUE(all.equal(unname(deg_level(gn_eu$gene,  pc_ref[["eu_gene"]])  - wt_level(gn_eu$gene,  pc_ref[["eu_gene"]])),  gn_reg$eu_gene)),
  isTRUE(all.equal(unname(deg_level(gn_tot$gene, pc_ref[["tot_gene"]]) - wt_level(gn_tot$gene, pc_ref[["tot_gene"]])), gn_reg$tot_gene)))
message("Check passed: degraded level - Wildtype level = the EU-seq and total RNA-seq change columns (both loop sets)")
dg_raw <- rbind(
  data.frame(key = loop_key(bn_reg), eu_gene_deg  = unname(deg_level(eu_s$gene,  pc_ref[["eu_gene"]])),
             tot_gene_deg = unname(deg_level(tot_s$gene, pc_ref[["tot_gene"]]))),
  data.frame(key = loop_key(gn_reg), eu_gene_deg  = unname(deg_level(gn_eu$gene,  pc_ref[["eu_gene"]])),
             tot_gene_deg = unname(deg_level(gn_tot$gene, pc_ref[["tot_gene"]]))))
for (v in names(dg_vars)) {   # a loop in both tables has the same degraded level in both
  same <- tapply(dg_raw[[v]], dg_raw$key, function(z) length(unique(round(z, 8))) == 1)
  stopifnot(all(same))
}
dg_raw <- dg_raw[!duplicated(dg_raw$key), ]
dg_loops <- all_loops
for (v in names(dg_vars)) dg_loops[[v]] <- dg_raw[[v]][match(loop_key(all_loops), dg_raw$key)]
stopifnot(!anyNA(dg_loops[, names(dg_vars)]), all(is.finite(as.matrix(dg_loops[, names(dg_vars)]))),
          nrow(dg_loops) == nrow(all_loops), nrow(dg_loops) == nrow(dg_raw))
message("All loops: ", format(nrow(dg_loops), big.mark = ","),
        "   with no signal in POLR2A-degraded cells (tied at the lowest level): EU-seq ",
        sum(dg_loops$eu_gene_deg == min(dg_loops$eu_gene_deg)), ", total RNA-seq ",
        sum(dg_loops$tot_gene_deg == min(dg_loops$tot_gene_deg)),
        "   correlation of EU-seq and total RNA-seq level: ",
        round(cor(dg_loops$eu_gene_deg, dg_loops$tot_gene_deg), 3))

set.seed(CV_SEED)
dg_res <- pooled_curves(dg_loops, dg_vars, "All loops")
message("\nTop decile minus bottom decile of the POLR2A-degraded level (pooled log2FC), bootstrap 95% interval and p")
print(dg_res$ext, row.names = FALSE, digits = 3, width = 250)
dg_rr <- dg_res$roll
dg_rr$variable <- factor(dg_rr$variable, levels = unname(dg_vars))
dg_rr$pct <- to_pct(dg_rr$pooled_log2FC); dg_rr$pct_lo <- to_pct(dg_rr$lo); dg_rr$pct_hi <- to_pct(dg_rr$hi)
stopifnot(all(is.finite(c(dg_rr$pct, dg_rr$pct_lo, dg_rr$pct_hi))))

# both conditions together: colour = assay, line type = condition (Wildtype first)
both <- rbind(lv_rr, dg_rr)
both$assay     <- factor(ifelse(grepl("Nascent", both$variable), "Nascent (EU-seq)", "Total RNA-seq"),
                         levels = c("Nascent (EU-seq)", "Total RNA-seq"))
both$condition <- factor(ifelse(grepl("Wildtype", both$variable), "Wildtype cells", "POLR2A-degraded cells"),
                         levels = c("Wildtype cells", "POLR2A-degraded cells"))
both$grp <- paste(both$assay, both$condition)
stopifnot(!anyNA(both$assay), !anyNA(both$condition), length(unique(both$grp)) == 4)
message("\nRange of the rolling curve (% change relative to the pooled ratio of all loops):")
print(aggregate(pct ~ assay + condition, data = both, FUN = function(z) round(range(z), 1)))

assay_col <- c("Nascent (EU-seq)" = "#E69F00", "Total RNA-seq" = "#009E73")
p_both <- ggplot(both, aes(group = grp)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_ribbon(aes(100 * percentile, ymin = pct_lo, ymax = pct_hi, fill = assay), alpha = 0.06) +
  geom_line(aes(100 * percentile, pct, colour = assay, linetype = condition), linewidth = 0.9) +
  scale_colour_manual(values = assay_col, name = NULL) +
  scale_fill_manual(values = assay_col, name = NULL) +
  scale_linetype_manual(values = c("Wildtype cells" = "solid", "POLR2A-degraded cells" = "22"), name = NULL) +
  labs(title = "Pooled loop change versus transcription level, Wildtype and POLR2A-degraded cells: all loops",
       subtitle = paste0("N = ", n_all, " loops; line = rolling window (", 100 * 2 * ROLL_HALF_WIDTH,
                         "% of loops), band = bootstrap 95% interval",
                         "\ny = change in contact ratio relative to the pooled ratio of all loops"),
       x = "Percentile of transcription level (0 = lowest)", y = "% contact change") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9), legend.position = "bottom")
print(p_both)
save_plot_svg(p_both, "pooled-loops-all", "wt-deg-level-pct", w = 10, h = 6)
