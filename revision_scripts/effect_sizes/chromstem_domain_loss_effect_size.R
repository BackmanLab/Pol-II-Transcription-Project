chromstem_effect_sizes <- function(
    in_dir   = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/chromSTEM",
    out_file = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/chromSTEM/ChromSTEM_effect_sizes.csv",
    actd_file        = "Chrom_STEM_Domains.csv",   # columns Cond, D, Rf, CVC
    polr2a_file      = "alldomains_081725.csv",    # POLR2A domains: columns D, Rmax, Rmin, PackEff, Density, cond (same file as the domain-loss script)
    polr2a_cond      = c(ctrl = "Pol2_DMSO", trt = "Pol2_6hrs_Aux"),   # DMSO control (N=86 in Fig. 9) and POLR2A-degraded (N=61)
    polr2a_xlsx_check = "POL2AID.xlsx",            # optional cross-check of the degraded-cell values (columns D_all, rfit_max_all, Density_all)
    # manuscript measure name = column name in each file
    actd_cols   = c(D = "D",     Rf = "Rf",   CVC = "CVC"),
    polr2a_cols = c(D = "D",     Rf = "Rmax", CVC = "Density"),
    expected_n  = c(ActD_ctrl = 71, ActD_trt = 48, POLR2A_ctrl = 86, POLR2A_trt = 61),
    n_boot = 10000,
    seed   = 102726) {
  
  # ---- input files ------------------------------------------------------------
  paths <- file.path(in_dir, c(actd_file, polr2a_file))
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files) > 0) stop("Input file(s) not found: ", paste(missing_files, collapse = ", "))
  
  read_any <- function(path) {
    if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is needed to read ", path)
      as.data.frame(readxl::read_excel(path))
    } else {
      read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    }
  }
  check_cols <- function(d, cols, label) {
    miss <- setdiff(cols, colnames(d))
    if (length(miss) > 0) stop("Columns not found in ", label, ": ", paste(miss, collapse = ", "),
                               "\nColumns present: ", paste(colnames(d), collapse = ", "))
    for (m in cols) if (!is.numeric(d[[m]])) stop("Column '", m, "' in ", label, " is not numeric")
  }
  
  # ---- ActD vs DMSO (one CSV, Cond column) -----------------------------------
  actd <- read_any(paths[1])
  check_cols(actd, c(actd_cols), actd_file)
  if (!"Cond" %in% colnames(actd)) stop("Column 'Cond' not found in ", actd_file, "\nColumns present: ", paste(colnames(actd), collapse = ", "))
  cat("Cond values in ", actd_file, ":\n", sep = ""); print(table(actd$Cond))
  stopifnot(all(c("DMSO", "ActD") %in% actd$Cond))
  actd_ctrl <- actd[actd$Cond == "DMSO", actd_cols, drop = FALSE]
  actd_trt  <- actd[actd$Cond == "ActD", actd_cols, drop = FALSE]
  
  # ---- POLR2A-degraded vs DMSO (one CSV, cond column) --------------------------------
  p_all <- read_any(paths[2])
  check_cols(p_all, polr2a_cols, polr2a_file)
  if (!"cond" %in% colnames(p_all)) stop("Column 'cond' not found in ", polr2a_file, "\nColumns present: ", paste(colnames(p_all), collapse = ", "))
  cat("cond values in ", polr2a_file, ":\n", sep = ""); print(table(p_all$cond))
  if (any(!polr2a_cond %in% p_all$cond)) stop("cond value(s) not found in ", polr2a_file, ": ", paste(setdiff(polr2a_cond, p_all$cond), collapse = ", "))
  p_ctrl <- p_all[p_all$cond == polr2a_cond[["ctrl"]], polr2a_cols, drop = FALSE]
  p_trt  <- p_all[p_all$cond == polr2a_cond[["trt"]],  polr2a_cols, drop = FALSE]
  
  # all-zero rows are placeholders, not domains: drop and report; any other zero or NA stops the run
  drop_empty <- function(d, label) {
    empty <- apply(d, 1, function(r) all(!is.na(r) & r == 0))
    cat(label, ": ", sum(empty), " all-zero row(s) dropped\n", sep = "")
    d <- d[!empty, , drop = FALSE]
    if (anyNA(d)) stop("NA values in ", label)
    if (any(d == 0)) stop("Zero values remain in ", label, " after dropping all-zero rows")
    d
  }
  p_ctrl <- drop_empty(p_ctrl, polr2a_cond[["ctrl"]])
  p_trt  <- drop_empty(p_trt,  polr2a_cond[["trt"]])
  for (nm in c("actd_ctrl", "actd_trt")) if (anyNA(get(nm))) stop("NA values in ", nm)
  
  # cross-check on the live data: the degraded-cell values of the earlier xlsx export must be the same numbers as in the CSV (same measures, same domains)
  xl_path <- file.path(in_dir, polr2a_xlsx_check)
  if (file.exists(xl_path)) {
    xl <- read_any(xl_path)
    xl_cols <- c(D = "D_all", Rf = "rfit_max_all", CVC = "Density_all")
    check_cols(xl, xl_cols, polr2a_xlsx_check)
    xl <- drop_empty(xl[, xl_cols, drop = FALSE], polr2a_xlsx_check)
    same <- nrow(xl) == nrow(p_trt) && all(sapply(seq_along(xl_cols), function(i) isTRUE(all.equal(sort(xl[[i]]), sort(p_trt[[i]]), tolerance = 1e-6))))
    if (same) cat("Cross-check passed: ", polr2a_cond[["trt"]], " values in ", polr2a_file, " equal those in ", polr2a_xlsx_check, " (D, Rf, CVC)\n", sep = "")
    else warning("The ", polr2a_cond[["trt"]], " values in ", polr2a_file, " differ from ", polr2a_xlsx_check, " (rows: ", nrow(p_trt), " vs ", nrow(xl), "); the effect sizes use ", polr2a_file)
  } else cat("Cross-check skipped: ", xl_path, " not found\n", sep = "")
  
  # group sizes must match the figure legends
  got_n <- c(ActD_ctrl = nrow(actd_ctrl), ActD_trt = nrow(actd_trt),
             POLR2A_ctrl = nrow(p_ctrl), POLR2A_trt = nrow(p_trt))
  cat("\nDomains per group:\n"); print(got_n)
  stopifnot(identical(names(got_n), names(expected_n)))
  if (any(got_n != expected_n)) stop("Domain counts differ from expected_n: ",
                                     paste(names(got_n)[got_n != expected_n], collapse = ", "))
  
  # ---- Hedges' g --------------------------------------------------------------
  hedges_g <- function(ctrl, trt) {
    n1 <- length(ctrl)
    n2 <- length(trt)
    sp <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / (n1 + n2 - 2))
    d  <- (mean(trt) - mean(ctrl)) / sp
    j  <- 1 - 3 / (4 * (n1 + n2) - 9)
    d * j
  }
  
  # Self-check on the live data: d from the pooled-variance t statistic must equal d computed directly.
  chk_ctrl <- actd_ctrl[[actd_cols[1]]]
  chk_trt  <- actd_trt[[actd_cols[1]]]
  chk_t    <- unname(t.test(chk_trt, chk_ctrl, var.equal = TRUE)$statistic)
  chk_d    <- chk_t * sqrt(1 / length(chk_trt) + 1 / length(chk_ctrl))
  chk_j    <- 1 - 3 / (4 * (length(chk_trt) + length(chk_ctrl)) - 9)
  stopifnot(isTRUE(all.equal(hedges_g(chk_ctrl, chk_trt), chk_d * chk_j)))
  cat("\nSelf-check passed: g matches the pooled t-statistic conversion\n")
  
  # ---- bootstrap: domains resampled within each group (no tomogram or cell ids in the files) ----
  groups <- list(
    list(comparison = "ActD_1hr_vs_DMSO",         cols = actd_cols,   ctrl = actd_ctrl, trt = actd_trt),
    list(comparison = "POLR2A_deg_6hr_vs_DMSO",   cols = polr2a_cols, ctrl = p_ctrl,    trt = p_trt))
  
  set.seed(seed)
  rows <- list()
  for (gp in groups) {
    for (m in names(gp$cols)) {
      ctrl <- gp$ctrl[[gp$cols[[m]]]]
      trt  <- gp$trt[[gp$cols[[m]]]]
      g    <- hedges_g(ctrl, trt)
      bg   <- replicate(n_boot, hedges_g(sample(ctrl, replace = TRUE), sample(trt, replace = TRUE)))
      ci   <- quantile(bg, c(0.025, 0.975), names = FALSE)
      stopifnot(ci[1] <= g, g <= ci[2])
      rows[[length(rows) + 1]] <- data.frame(
        assay         = "ChromSTEM",
        comparison    = gp$comparison,
        measure       = m,
        source_column = gp$cols[[m]],
        n_ctrl        = length(ctrl),
        n_trt         = length(trt),
        mean_ctrl     = mean(ctrl),
        mean_trt      = mean(trt),
        hedges_g      = g,
        ci_lo         = ci[1],
        ci_hi         = ci[2]
      )
    }
  }
  
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  print(res, digits = 3)
  write.csv(res, out_file, row.names = FALSE)
  cat("Wrote:", out_file, "\n")
  invisible(res)
}

chromstem_effect <- chromstem_effect_sizes()