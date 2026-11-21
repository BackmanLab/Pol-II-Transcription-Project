eustorm_effect_sizes <- function(
    # the four processed_data CSVs of EUstorm.R (copy them from the processed_data folder of 20240524_eu_STORM to in_dir)
  in_dir   = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/SMLM",
  out_file = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/SMLM/EUSTORM_effect_sizes.csv",
  # files_by_condition of EUstorm.R: the two files of a condition are the two biological replicates (original and 071724)
  files_by_condition = list(wt   = c("WT_DMSO_nucleus.csv",      "WT_DMSO_nucleus_071724.csv"),
                            pol2 = c("POLR2A_Auxin_nucleus.csv", "POLR2A_Auxin_nucleus_071724.csv")),
  control = "wt", treated = "pol2",
  # Fig. 2E (cluster size) and Fig. 2F (blink density); column names as in `metrics` of EUstorm.R
  metrics = data.frame(col   = c("cluster_mean_size", "num_blink_perUm2"),
                       panel = c("Fig. 2E", "Fig. 2F"), stringsAsFactors = FALSE),
  # Fig. 2E-F legend: N (the files hold 28 / 18 ROIs from 15 / 11 cells) and mean clusters per cell (printed for comparison only)
  expected_rois     = c(wt = 28, pol2 = 18),
  legend_clusters   = c(wt = 802.40, pol2 = 549.64),
  n_boot = 10000,
  seed   = 102726) {
  
  library(dplyr)
  library(tibble)
  library(stringr)
  library(purrr)
  
  # ---- input checks ------------------------------------------------------------
  all_files <- unlist(files_by_condition)
  paths <- file.path(in_dir, all_files)
  if (any(!file.exists(paths))) stop("Input file(s) not found (copy the four EU-STORM processed_data CSVs to ", in_dir, "):\n", paste(paths[!file.exists(paths)], collapse = "\n"))
  if (!all(c(control, treated) %in% names(files_by_condition))) stop("control / treated not in files_by_condition")
  if (any(!c(control, treated) %in% names(expected_rois)) || any(!c(control, treated) %in% names(legend_clusters))) stop("expected_rois / legend_clusters need names ", control, " and ", treated)
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  
  # ---- read_eu_storm_file() of EUstorm.R (unchanged apart from the column check) -------------
  read_eu_storm_file <- function(path, condition_label) {
    raw <- read.csv(path, row.names = 1, check.names = FALSE)
    tidy <- as.data.frame(t(raw)) %>%
      tibble::rownames_to_column("roi_id") %>%
      tibble::as_tibble()
    need <- c("num_blink_nucleus", "num_cluster_nucleus", "cluster_mean_size_nucleus", "roi_size_nucleus")
    if (any(!need %in% names(tidy))) stop("Columns missing in ", basename(path), ": ", paste(need[!need %in% names(tidy)], collapse = ", "),
                                          "\nColumns present: ", paste(names(tidy), collapse = ", "))
    if (!"num_blink_nucleus per um^2" %in% names(tidy)) {
      tidy <- tidy %>% dplyr::mutate(`num_blink_nucleus per um^2` = num_blink_nucleus / roi_size_nucleus)
    }
    if (!"num_cluster_nucleus per um^2" %in% names(tidy)) {
      tidy <- tidy %>% dplyr::mutate(`num_cluster_nucleus per um^2` = num_cluster_nucleus / roi_size_nucleus)
    }
    tidy %>%
      dplyr::rename(
        num_blinks         = num_blink_nucleus,
        num_clust          = num_cluster_nucleus,
        cluster_mean_size  = cluster_mean_size_nucleus,
        roi                = roi_size_nucleus,
        num_cluster_perUm2 = `num_cluster_nucleus per um^2`,
        num_blink_perUm2   = `num_blink_nucleus per um^2`
      ) %>%
      dplyr::mutate(
        cell_id  = stringr::str_extract(roi_id, "^cell[0-9]+"),
        analysis = condition_label
      )
  }
  
  # cohens_d() of EUstorm.R (used to check the Hedges' g computed below)
  cohens_d <- function(x, y) {
    n1 <- length(x); n2 <- length(y)
    pooled_sd <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
    (mean(x) - mean(y)) / pooled_sd
  }
  hedges_g <- function(ctrl, trt, strict = FALSE) {   # treated minus control
    n1 <- length(ctrl); n2 <- length(trt)
    sp <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / (n1 + n2 - 2))
    if (!is.finite(sp) || sp == 0) { if (strict) stop("Pooled SD is zero or undefined"); return(NA_real_) }
    ((mean(trt) - mean(ctrl)) / sp) * (1 - 3 / (4 * (n1 + n2) - 9))
  }
  
  # ---- load: one row per ROI, tagged with condition and replicate file ---------------------------
  nuc <- purrr::imap_dfr(files_by_condition, function(files, cond) {
    purrr::map_dfr(files, function(f) read_eu_storm_file(file.path(in_dir, f), cond) %>% dplyr::mutate(file = f))
  })
  if (anyNA(nuc$cell_id)) stop("ROI names without a 'cellN' prefix: ", paste(head(nuc$roi_id[is.na(nuc$cell_id)]), collapse = ", "))
  # the cell number repeats across the replicate files, so a cell is identified by condition + file + cell number
  nuc$cell_key <- paste(nuc$analysis, nuc$file, nuc$cell_id, sep = "|")
  cat("\nROIs and cells per condition and replicate file:\n")
  print(as.data.frame(nuc %>% dplyr::group_by(analysis, file) %>% dplyr::summarise(n_roi = dplyr::n(), n_cell = dplyr::n_distinct(cell_key), .groups = "drop")))
  
  # ---- self-check on the live data: ROI counts must reproduce the N of the Fig. 2E-F legend -------------
  n_rois  <- sapply(c(control, treated), function(cn) sum(nuc$analysis == cn))
  n_cells <- sapply(c(control, treated), function(cn) dplyr::n_distinct(nuc$cell_key[nuc$analysis == cn]))
  n_cells_merged <- sapply(c(control, treated), function(cn) dplyr::n_distinct(paste(nuc$analysis, nuc$cell_id)[nuc$analysis == cn]))
  cat("ROIs:", paste(names(n_rois), n_rois, collapse = "  "), "| legend N:", paste(names(expected_rois), expected_rois, collapse = "  "),
      "\nCells (condition + file + cell number):", paste(names(n_cells), n_cells, collapse = "  "),
      "| if cell numbers were merged across the two files (as cell_id in EUstorm.R):", paste(names(n_cells_merged), n_cells_merged, collapse = "  "), "\n")
  if (any(n_rois != expected_rois[names(n_rois)])) stop("ROI counts differ from the legend N (expected_rois); see the counts printed above")
  cat("Self-check passed: ROI counts match the N of the legend\n")
  
  # ---- one row per cell: mean over that cell's ROIs, as nuc.percell of EUstorm.R -----------------
  percell <- nuc %>%
    dplyr::group_by(analysis, file, cell_key) %>%
    dplyr::summarise(clusters_total = sum(num_clust, na.rm = TRUE), n_roi = dplyr::n(),   # before dplyr::across(), which overwrites num_clust with its mean
                     dplyr::across(c(num_blinks, num_clust, cluster_mean_size, roi, num_cluster_perUm2, num_blink_perUm2), \(x) mean(x, na.rm = TRUE)),
                     .groups = "drop")
  for (m in metrics$col) if (any(!is.finite(nuc[[m]])) || any(!is.finite(percell[[m]]))) stop("Non-finite value of ", m)
  
  # ---- legend value 'mean clusters per cell' (comparison only) -----------------------------------
  cand <- sapply(c(control, treated), function(cn) {
    pc <- percell[percell$analysis == cn, ]; ro <- nuc[nuc$analysis == cn, ]
    c(mean_over_ROIs = mean(ro$num_clust, na.rm = TRUE), mean_of_cell_means = mean(pc$num_clust), mean_of_cell_totals = mean(pc$clusters_total))
  })
  cat("Mean clusters per cell as the legend gives it (", paste(names(legend_clusters), legend_clusters, collapse = "  "), ") against what these files give:\n", sep = "")
  print(round(cand, 2))
  if (!any(apply(abs(sweep(cand, 2, legend_clusters[colnames(cand)])), 1, max) < 0.01)) cat("  -> none of the three definitions reproduces the legend values\n")
  
  # ---- Hedges' g per metric at two levels ----------------------------------------------------------
  # roi  : every ROI is a unit (the N of the figure and of the legend; ROIs of one cell are not independent)
  # cell : one value per cell, the mean of its ROIs (as the cell_level tests of EUstorm.R)
  # bootstrap: cells are resampled with replacement within each condition and replicate file (all ROIs of a resampled cell are taken
  # for the ROI level), so the interval accounts for ROIs sharing a cell
  set.seed(seed)
  cell_boot <- function(df, m, level) {   # one resample of the cells of one condition; returns the values entering g
    lst <- lapply(split(df, df$file), function(d) {
      keys <- unique(d$cell_key)
      pick <- keys[sample.int(length(keys), replace = TRUE)]
      if (level == "cell") d$val[match(pick, d$cell_key)] else unlist(lapply(pick, function(k) d$val[d$cell_key == k]), use.names = FALSE)
    })
    unlist(lst, use.names = FALSE)
  }
  rows <- list()
  for (level in c("roi", "cell")) for (i in seq_len(nrow(metrics))) {
    m <- metrics$col[i]
    src <- if (level == "roi") nuc else percell
    cd  <- src[src$analysis == control, c("file", "cell_key", m)]; names(cd)[3] <- "val"
    td  <- src[src$analysis == treated, c("file", "cell_key", m)]; names(td)[3] <- "val"
    ctrl <- cd$val; trt <- td$val
    g <- hedges_g(ctrl, trt, strict = TRUE)
    # g must equal the repo's Cohen's d times the Hedges correction, and the pooled-t conversion
    j <- 1 - 3 / (4 * (length(ctrl) + length(trt)) - 9)
    chk_t <- unname(t.test(trt, ctrl, var.equal = TRUE)$statistic)
    stopifnot(isTRUE(all.equal(g, cohens_d(trt, ctrl) * j)),
              isTRUE(all.equal(g, chk_t * sqrt(1 / length(trt) + 1 / length(ctrl)) * j)))
    bg <- replicate(n_boot, hedges_g(cell_boot(cd, m, level), cell_boot(td, m, level)))
    ci <- quantile(bg, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    if (!(ci[1] <= g && g <= ci[2])) warning("Observed g outside its bootstrap interval: ", m, " (", level, ")")
    rows[[length(rows) + 1]] <- data.frame(
      assay = "SMLM_EU", comparison = "POLR2A-Aux_vs_WT-DMSO", measure = m, panel = metrics$panel[i], level = level,
      n_ctrl = length(ctrl), n_trt = length(trt),
      n_cells_ctrl = dplyr::n_distinct(cd$cell_key), n_cells_trt = dplyr::n_distinct(td$cell_key),
      mean_ctrl = mean(ctrl), mean_trt = mean(trt), hedges_g = g, ci_lo = ci[1], ci_hi = ci[2],
      welch_p = t.test(trt, ctrl)$p.value, ks_p = suppressWarnings(ks.test(trt, ctrl)$p.value), stringsAsFactors = FALSE)
  }
  cat("Self-check passed: g matches the Cohen's d of EUstorm.R x the Hedges correction and the pooled t-statistic conversion (both levels)\n")
  
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  cat("\nhedges_g = POLR2A-degraded minus Wildtype; level roi = ROIs as units (the N of the figure), cell = one mean per cell; ks_p as in the CDF panels (roi level)\n\n")
  print(res, digits = 3)
  write.csv(res, out_file, row.names = FALSE)
  cat("Wrote:", out_file, "\n")
  invisible(res)
}

eustorm_effect <- eustorm_effect_sizes()