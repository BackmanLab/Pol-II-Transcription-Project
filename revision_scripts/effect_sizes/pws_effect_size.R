pws_effect_sizes <- function(
    data_file    = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/PWS/Pols_PWS_replot_042624.csv",
    out_file     = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/PWS/PWS_effect_sizes.csv",
    measure_cols = c("D", "FMM", "Diffusion"),
    control_name = "DMSO",
    expected_replicates = c(DMSO = 3, POLR2A_deg_6hr = 2, ActD_1hr = 3),
    # One row per biological replicate = one Path label in the CSV.
    replicate_map = data.frame(
      condition = c("DMSO", "DMSO", "DMSO",
                    "POLR2A_deg_6hr", "POLR2A_deg_6hr",
                    "ActD_1hr", "ActD_1hr", "ActD_1hr"),
      Path      = c("WT HCT116\\DMSO Only_1", "WT HCT116\\DMSO Only_2", "WT HCT116\\DMSO Only_Rep1",
                    "Auxin_6hr", "Auxin_6hrs",
                    "WT HCT116\\ActD_1hr_4uM_1", "WT HCT116\\Act_1hr_4uM_2", "WT HCT116\\ActD_5ugmL_Rep1"),
      stringsAsFactors = FALSE),
    n_boot = 10000,
    seed   = 102726) {
  
  # ---- read and check columns -----------------------------------------------
  if (!file.exists(data_file)) stop("Input file not found: ", data_file)
  x <- read.csv(data_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  
  needed_cols <- c("Path", measure_cols)
  missing_cols <- setdiff(needed_cols, colnames(x))
  if (length(missing_cols) > 0) {
    stop("Columns not found in CSV: ", paste(missing_cols, collapse = ", "),
         "\nColumns present: ", paste(colnames(x), collapse = ", "))
  }
  for (m in measure_cols) if (!is.numeric(x[[m]])) stop("Column '", m, "' is not numeric")
  
  cat("ROIs per Path label:\n"); print(table(x$Path))
  
  # ---- assign biological replicates (one per Path label) ----------------------
  if (anyDuplicated(replicate_map$Path) > 0) stop("A Path label appears more than once in replicate_map")
  not_found <- setdiff(replicate_map$Path, x$Path)
  if (length(not_found) > 0) stop("Path labels in replicate_map not found in CSV: ", paste(not_found, collapse = ", "))
  
  n_reps <- table(replicate_map$condition)
  stopifnot(identical(sort(names(n_reps)), sort(names(expected_replicates))))
  stopifnot(all(n_reps[names(expected_replicates)] == expected_replicates))
  if (!control_name %in% names(expected_replicates)) stop("control_name is not a condition in replicate_map")
  cat("\nReplicates per condition:\n"); print(n_reps)
  
  excluded <- table(x$Path[!x$Path %in% replicate_map$Path])
  cat("\nPath labels not used (excluded):\n"); print(excluded)
  
  x$rep_id    <- match(x$Path, replicate_map$Path)
  x$condition <- replicate_map$condition[x$rep_id]
  used <- x[!is.na(x$rep_id), ]
  
  # FMM scale check on the rows actually used
  if ("FMM" %in% measure_cols) {
    big <- !is.na(used$FMM) & used$FMM > 1e-15
    if (any(big)) stop(sum(big), " rows used have FMM > 1e-15 (expected ~1e-19): ",
                       paste(unique(used$Path[big]), collapse = "; "))
  }
  
  # structure check: ROIs and means per replicate
  rep_summary <- do.call(rbind, lapply(seq_len(nrow(replicate_map)), function(r) {
    d <- used[used$rep_id == r, ]
    data.frame(condition = replicate_map$condition[r], rep = r, Path = replicate_map$Path[r],
               ROIs = nrow(d), D_mean = mean(d$D, na.rm = TRUE), Diffusion_mean = mean(d$Diffusion, na.rm = TRUE))
  }))
  cat("\nReplicate structure:\n"); print(rep_summary, row.names = FALSE, digits = 3)
  cat("\nROIs per condition:\n")
  print(aggregate(ROIs ~ condition, data = rep_summary, FUN = sum))
  
  treated_names <- setdiff(names(expected_replicates), control_name)
  
  # ---- Hedges' g ------------------------------------------------------------
  hedges_g <- function(ctrl, trt) {
    n1 <- length(ctrl)
    n2 <- length(trt)
    sp <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / (n1 + n2 - 2))
    d  <- (mean(trt) - mean(ctrl)) / sp
    j  <- 1 - 3 / (4 * (n1 + n2) - 9)
    d * j
  }
  
  # Self-check on the live data: d from the pooled-variance t statistic must equal d computed directly.
  m1 <- measure_cols[1]
  chk_ctrl <- used[[m1]][used$condition == control_name & !is.na(used[[m1]])]
  chk_trt  <- used[[m1]][used$condition == treated_names[1] & !is.na(used[[m1]])]
  chk_t    <- unname(t.test(chk_trt, chk_ctrl, var.equal = TRUE)$statistic)
  chk_d    <- chk_t * sqrt(1 / length(chk_trt) + 1 / length(chk_ctrl))
  chk_j    <- 1 - 3 / (4 * (length(chk_trt) + length(chk_ctrl)) - 9)
  stopifnot(isTRUE(all.equal(hedges_g(chk_ctrl, chk_trt), chk_d * chk_j)))
  cat("\nSelf-check passed: g matches the pooled t-statistic conversion\n")
  
  # ---- bootstrap: replicates fixed, ROIs resampled within each replicate ----
  # (2-3 replicates give only 3-10 distinct replicate-level resamples, so replicates are not resampled)
  vals_by_rep <- function(cond, m) {
    lapply(which(replicate_map$condition == cond), function(r) {
      v <- used[[m]][used$rep_id == r]
      v[!is.na(v)]
    })
  }
  draw <- function(rep_list) {
    unlist(lapply(rep_list, function(v) v[sample.int(length(v), replace = TRUE)]), use.names = FALSE)
  }
  flat <- function(rep_list) unlist(rep_list, use.names = FALSE)
  
  set.seed(seed)
  rows <- list()
  for (trt_name in treated_names) {
    for (m in measure_cols) {
      ctrl_l <- vals_by_rep(control_name, m)
      trt_l  <- vals_by_rep(trt_name, m)
      g  <- hedges_g(flat(ctrl_l), flat(trt_l))
      bg <- replicate(n_boot, hedges_g(draw(ctrl_l), draw(trt_l)))
      ci <- quantile(bg, c(0.025, 0.975), names = FALSE)
      stopifnot(ci[1] <= g, g <= ci[2])
      rows[[length(rows) + 1]] <- data.frame(
        assay       = "PWS",
        comparison  = paste0(trt_name, "_vs_", control_name),
        measure     = m,
        n_reps_ctrl = length(ctrl_l),
        n_reps_trt  = length(trt_l),
        n_roi_ctrl  = length(flat(ctrl_l)),
        n_roi_trt   = length(flat(trt_l)),
        mean_ctrl   = mean(flat(ctrl_l)),
        mean_trt    = mean(flat(trt_l)),
        hedges_g    = g,
        ci_lo       = ci[1],
        ci_hi       = ci[2]
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

pws_effect <- pws_effect_sizes()