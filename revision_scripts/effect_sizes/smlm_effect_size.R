smlm_effect_sizes <- function(
    in_dir   = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/SMLM",
    out_file = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/SMLM/SMLM_effect_sizes.csv",
    # control / treated labels are the `condition` values used in EDU.K9_STORM.R
    # breaks = radius cutoffs in nm: Small <= b1, Average b1-b2, Large > b2 (same as EDU.K9_STORM.R)
    datasets = list(
      EdU = list(file = "EDU_STORM.parquet", control = "WT-DMSO", treated = "POLR2A-Aux", breaks = c(80, 120),
                 expected_cells   = c(ctrl = 11, trt = 10),
                 expected_domains = c(ctrl = 6282, trt = 2937),
                 expected_bins    = list(ctrl = c(4054, 1197, 1031), trt = c(1915, 589, 433))),
      K9  = list(file = "K9_STORM.parquet", control = "WT", treated = "Auxin", breaks = c(40, 80),
                 expected_cells   = c(ctrl = 29, trt = 24),
                 expected_domains = c(ctrl = 9534, trt = 7525),
                 expected_bins    = list(ctrl = c(2332, 4668, 2534), trt = c(1707, 3599, 2219)))),
    normalization_radius_nm = 80,
    n_boot = 10000,
    seed   = 102726) {
  
  # ---- input checks ------------------------------------------------------------
  paths <- file.path(in_dir, vapply(datasets, function(d) d$file, character(1)))
  if (any(!file.exists(paths))) stop("Input file(s) not found: ", paste(paths[!file.exists(paths)], collapse = ", "))
  if (!requireNamespace("arrow", quietly = TRUE)) stop("R package 'arrow' is not installed")  # ARROWCHECK
  
  bin_labels <- c("Small", "Average", "Large")
  
  # ---- Hedges' g (treated minus control) ----------------------------------------
  hedges_g <- function(ctrl, trt) {
    n1 <- length(ctrl)
    n2 <- length(trt)
    sp <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / (n1 + n2 - 2))
    if (!is.finite(sp) || sp == 0) stop("Pooled SD is zero or undefined")
    d <- (mean(trt) - mean(ctrl)) / sp
    j <- 1 - 3 / (4 * (n1 + n2) - 9)
    d * j
  }
  
  set.seed(seed)
  rows <- list()
  
  for (nm in names(datasets)) {
    ds <- datasets[[nm]]
    f  <- file.path(in_dir, ds$file)
    x  <- as.data.frame(arrow::read_parquet(f))  # READPARQUET
    for (cc in c("condition", "cell_id", "seedCluster_size", "density")) {
      if (!cc %in% colnames(x)) stop("Column '", cc, "' not found in ", ds$file, "\nColumns present: ", paste(colnames(x), collapse = ", "))
    }
    x$condition <- as.character(x$condition)
    cat("\n==== ", nm, " (", ds$file, ") ====\nDomains per condition:\n", sep = ""); print(table(x$condition))
    for (cn in c(ds$control, ds$treated)) if (!cn %in% x$condition) stop("Condition '", cn, "' not found in ", ds$file)
    
    x <- x[x$condition %in% c(ds$control, ds$treated), ]
    if (anyNA(x[, c("cell_id", "seedCluster_size", "density")])) stop("NA values in cell_id, seedCluster_size or density in ", ds$file)
    if (any(x$seedCluster_size <= 0)) stop("Non-positive seedCluster_size in ", ds$file)
    
    # same definitions as EDU.K9_STORM.R: radius = R * sqrt(seedCluster_size), R = 80 nm; bins right-closed
    x$radius_nm <- normalization_radius_nm * sqrt(x$seedCluster_size)
    x$bin <- cut(x$radius_nm, breaks = c(-Inf, ds$breaks, Inf), labels = bin_labels, right = TRUE)
    x$cell_key <- paste(x$condition, x$cell_id, sep = "|")   # cell ids repeat across conditions
    
    # ---- self-checks on the live data: counts must reproduce the manuscript ----------
    n_cells <- c(ctrl = length(unique(x$cell_key[x$condition == ds$control])),
                 trt  = length(unique(x$cell_key[x$condition == ds$treated])))
    n_dom   <- c(ctrl = sum(x$condition == ds$control), trt = sum(x$condition == ds$treated))
    bins_ctrl <- as.numeric(table(x$bin[x$condition == ds$control]))
    bins_trt  <- as.numeric(table(x$bin[x$condition == ds$treated]))
    cat("Cells:", paste(names(n_cells), n_cells, collapse = "  "), "| Domains:", paste(names(n_dom), n_dom, collapse = "  "), "\n")
    cat("Domains by size bin (", paste(bin_labels, collapse = "/"), "; cutoffs ", paste(ds$breaks, collapse = "/"), " nm)\n", sep = "")
    cat("  control:", bins_ctrl, "\n  treated:", bins_trt, "\n")
    if (any(n_cells != ds$expected_cells))     stop(nm, ": cell counts differ from expected_cells")
    if (any(n_dom != ds$expected_domains))     stop(nm, ": domain counts differ from expected_domains")
    if (any(bins_ctrl != ds$expected_bins$ctrl) || any(bins_trt != ds$expected_bins$trt)) stop(nm, ": size-bin counts differ from expected_bins")
    cat("Self-check passed: cell, domain and size-bin counts match the manuscript\n")
    
    # ---- one row per cell (zero counts kept) ----------------------------------------------
    cells <- unique(x[, c("cell_key", "condition")])
    cells <- cells[order(cells$cell_key), ]
    cnt <- table(factor(x$cell_key, levels = cells$cell_key), x$bin)
    pc <- data.frame(cell_key = cells$cell_key, condition = cells$condition,
                     domains_total   = as.numeric(rowSums(cnt)),
                     domains_small   = as.numeric(cnt[, "Small"]),
                     domains_average = as.numeric(cnt[, "Average"]),
                     domains_large   = as.numeric(cnt[, "Large"]),
                     mean_radius_nm  = as.numeric(tapply(x$radius_nm, factor(x$cell_key, levels = cells$cell_key), mean)),
                     mean_density    = as.numeric(tapply(x$density,   factor(x$cell_key, levels = cells$cell_key), mean)))
    stopifnot(sum(pc$domains_total) == sum(n_dom),
              all(pc$domains_total == pc$domains_small + pc$domains_average + pc$domains_large))
    ctrl_pc <- pc[pc$condition == ds$control, ]
    trt_pc  <- pc[pc$condition == ds$treated, ]
    
    # g must equal the pooled-variance t statistic conversion (live data)
    chk_t <- unname(t.test(trt_pc$domains_total, ctrl_pc$domains_total, var.equal = TRUE)$statistic)
    chk_d <- chk_t * sqrt(1 / nrow(trt_pc) + 1 / nrow(ctrl_pc))
    chk_j <- 1 - 3 / (4 * (nrow(trt_pc) + nrow(ctrl_pc)) - 9)
    stopifnot(isTRUE(all.equal(hedges_g(ctrl_pc$domains_total, trt_pc$domains_total), chk_d * chk_j)))
    cat("Self-check passed: g matches the pooled t-statistic conversion\n")
    
    # ---- Hedges' g per measure, bootstrap over cells within each condition ------------------
    measures <- c("domains_total", "domains_small", "domains_average", "domains_large", "mean_radius_nm", "mean_density")
    for (m in measures) {
      ctrl <- ctrl_pc[[m]]
      trt  <- trt_pc[[m]]
      g    <- hedges_g(ctrl, trt)
      bg   <- replicate(n_boot, hedges_g(sample(ctrl, replace = TRUE), sample(trt, replace = TRUE)))
      ci   <- quantile(bg, c(0.025, 0.975), names = FALSE)
      stopifnot(ci[1] <= g, g <= ci[2])
      rows[[length(rows) + 1]] <- data.frame(
        assay      = paste0("SMLM_", nm),
        comparison = paste0(ds$treated, "_vs_", ds$control),
        measure    = m,
        size_bins_nm = paste0("<=", ds$breaks[1], " / ", ds$breaks[1], "-", ds$breaks[2], " / >", ds$breaks[2]),
        n_cells_ctrl = length(ctrl),
        n_cells_trt  = length(trt),
        mean_ctrl    = mean(ctrl),
        mean_trt     = mean(trt),
        hedges_g     = g,
        ci_lo        = ci[1],
        ci_hi        = ci[2],
        stringsAsFactors = FALSE)
    }
  }
  
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  cat("\n")
  print(res, digits = 3)
  write.csv(res, out_file, row.names = FALSE)
  cat("Wrote:", out_file, "\n")
  invisible(res)
}

smlm_effect <- smlm_effect_sizes()