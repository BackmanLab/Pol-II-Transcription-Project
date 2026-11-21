euseq_effect_sizes <- function(
    # hg38 EU-seq rerun (chrM-scaled), same locations as EU_seq_zscore_analysis_updated.R
  base_dir   = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq",
  anal_matrix = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/euseq_pipeline/counts/euseq_rpkm_analysis.txt",
  sf_file    = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/scale_factors_mito.txt",   # chrM scale factors: columns condition, eu_sf
  # ENCODE HCT116 RSEM gene quantification (ENCFF435PHM), same table as the Hi-C effect-size scripts
  expr_file  = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference/ENCFF435PHM.tsv",
  out_file   = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/EUseq/EUseq_effect_sizes.csv",
  polr2a_samples = c("POLR2A_1", "POLR2A_2", "POLR2A_3"),
  wt_samples     = c("WT_1", "WT_2", "WT_3"),
  sf_wt_label    = "WT",
  sf_polr2a_label = "POLR2A_8hrs_aux",
  rpkm_threshold = 1.0,    # RPKM_THRESHOLD (mean scaled WT RPKM)
  lfc_threshold  = 1,      # LFC_THRESHOLD
  manuscript_up = 69, manuscript_down = 535,   # EU-seq counts quoted in the manuscript text (printed for comparison only)
  n_boot = 10000,
  seed   = 102726) {
  
  library(dplyr)
  library(readr)
  
  # ---- input checks ------------------------------------------------------------
  inputs <- c(anal_matrix, sf_file, expr_file)
  if (any(!file.exists(inputs))) stop("Input file(s) not found:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  
  # ---- chrM scale factors (EU_seq_zscore_analysis_updated.R, section 2c) ---------------------
  sf_tbl <- read.delim(sf_file, stringsAsFactors = FALSE)
  print(sf_tbl)
  if (any(!c("condition", "eu_sf") %in% colnames(sf_tbl))) stop("Columns 'condition' and 'eu_sf' required in ", basename(sf_file), "\nColumns present: ", paste(colnames(sf_tbl), collapse = ", "))
  if (any(!c(sf_wt_label, sf_polr2a_label) %in% sf_tbl$condition)) stop("Conditions '", sf_wt_label, "' and '", sf_polr2a_label, "' required in ", basename(sf_file), "\nConditions present: ", paste(sf_tbl$condition, collapse = ", "))
  if (sf_tbl$eu_sf[sf_tbl$condition == sf_wt_label] != 1) stop("WT scale factor is not 1 in ", basename(sf_file))
  sf_polr2a <- sf_tbl$eu_sf[sf_tbl$condition == sf_polr2a_label]
  sf_vec <- c(setNames(rep(sf_polr2a, length(polr2a_samples)), polr2a_samples), setNames(rep(1, length(wt_samples)), wt_samples))
  cat(sprintf("EU-seq chrM scale factor: POLR2A = %.6f (WT = 1)\n", sf_polr2a))
  
  # ---- analysis matrix (NoEU-subtracted RPKM), scaled ------------------------------------
  anal_mat <- read_tsv(anal_matrix, show_col_types = FALSE)
  need <- c("gene_id", polr2a_samples, wt_samples)
  if (any(!need %in% colnames(anal_mat))) stop("Columns missing in ", basename(anal_matrix), ": ", paste(setdiff(need, colnames(anal_mat)), collapse = ", "), "\nColumns present: ", paste(colnames(anal_mat), collapse = ", "))
  anal_mat$gene_id <- sub("\\..*", "", anal_mat$gene_id)
  if (anyDuplicated(anal_mat$gene_id) > 0) stop("Duplicated gene_id values in ", basename(anal_matrix))
  anal_fpkm <- as.matrix(anal_mat[, c(polr2a_samples, wt_samples)]); mode(anal_fpkm) <- "numeric"
  rownames(anal_fpkm) <- anal_mat$gene_id
  if (any(!is.finite(anal_fpkm)) || any(anal_fpkm < 0)) stop("Non-finite or negative RPKM values in ", basename(anal_matrix))
  cat("Analysis matrix:", nrow(anal_fpkm), "genes x", ncol(anal_fpkm), "samples\n")
  anal_scaled <- sweep(anal_fpkm, 2, sf_vec[colnames(anal_fpkm)], FUN = "*")
  
  # ---- per-gene statistics (section 5 of EU_seq_zscore_analysis.R) --------------------------
  log2_anal <- log2(anal_scaled + 1)
  wt_mean   <- rowMeans(anal_scaled[, wt_samples])
  keep      <- wt_mean >= rpkm_threshold
  log2_filt <- log2_anal[keep, ]
  cat("Genes passing RPKM filter (mean scaled WT RPKM >=", rpkm_threshold, "):", nrow(log2_filt), "\n")
  polr2a_vals <- log2_filt[, polr2a_samples]
  wt_vals     <- log2_filt[, wt_samples]
  n1 <- length(polr2a_samples); n2 <- length(wt_samples)
  lfc <- rowMeans(polr2a_vals) - rowMeans(wt_vals)                          # treated minus control, as in the repo
  z_repo <- lfc / (apply(cbind(polr2a_vals, wt_vals), 1, sd) + 1e-10)       # repo effect size: LFC / SD of all replicates
  # Hedges' g per gene: LFC / within-group pooled SD x small-sample correction (treated minus control)
  sp <- sqrt(((n1 - 1) * apply(polr2a_vals, 1, var) + (n2 - 1) * apply(wt_vals, 1, var)) / (n1 + n2 - 2))
  j  <- 1 - 3 / (4 * (n1 + n2) - 9)
  g  <- ifelse(sp > 0, lfc / sp * j, NA_real_)
  cat("Genes with zero within-group SD (g undefined, excluded from g):", sum(is.na(g)), "\n")
  # self-check on live data: per-gene g must equal the pooled t-statistic conversion
  i_chk <- which(!is.na(g))[1]
  chk_t <- unname(t.test(polr2a_vals[i_chk, ], wt_vals[i_chk, ], var.equal = TRUE)$statistic)
  stopifnot(isTRUE(all.equal(unname(g[i_chk]), chk_t * sqrt(1 / n1 + 1 / n2) * j)))
  cat("Self-check passed: per-gene g matches the pooled t-statistic conversion\n")
  stopifnot(isTRUE(all.equal(unname(lfc), unname(rowMeans(log2(anal_scaled[keep, polr2a_samples] + 1)) - rowMeans(log2(anal_scaled[keep, wt_samples] + 1))))))
  
  # ---- expression groups (section 6 of EU_seq_zscore_analysis.R: TPM < 1 unexpressed, tertiles of TPM >= 1) ----
  expr_raw <- read_tsv(expr_file, show_col_types = FALSE)
  if (any(!c("gene_id", "TPM") %in% colnames(expr_raw))) stop("Columns 'gene_id' and 'TPM' required in ", basename(expr_file), "\nColumns present: ", paste(colnames(expr_raw), collapse = ", "))
  expr <- expr_raw %>% dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>% dplyr::group_by(gene_id) %>% dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")
  tpm_exp <- expr$TPM[expr$TPM >= 1]
  q1 <- quantile(tpm_exp, 1/3, na.rm = TRUE); q2 <- quantile(tpm_exp, 2/3, na.rm = TRUE)
  cat("Expression group thresholds (TPM): Unexpressed <1 | Low <", round(q1, 2), "| Medium <", round(q2, 2), "| High >=", round(q2, 2), "\n")
  expr <- expr %>% dplyr::mutate(expr_group = dplyr::case_when(TPM < 1 ~ "Unexpressed", TPM < q1 ~ "Low", TPM < q2 ~ "Medium", TRUE ~ "High"))
  res_df <- data.frame(gene_id = rownames(log2_filt), lfc = lfc, g = g, z_repo = z_repo, stringsAsFactors = FALSE) %>%
    dplyr::left_join(expr %>% dplyr::select(gene_id, expr_group), by = "gene_id")
  cat("Genes passing the RPKM filter with a TPM value:", sum(!is.na(res_df$expr_group)), "/", nrow(res_df), "\n")
  if (mean(!is.na(res_df$expr_group)) < 0.8) stop("Fewer than 80% of EU-seq genes matched the TPM table; check gene_id formats")
  cat("Genes with |LFC| >", lfc_threshold, "(all genes passing the filter): Up =", sum(res_df$lfc > lfc_threshold), " Down =", sum(res_df$lfc < -lfc_threshold),
      "  | manuscript text: Up =", manuscript_up, " Down =", manuscript_down, "\n")
  
  # ---- summary per expression group, percentile bootstrap over genes --------------------------
  set.seed(seed)
  boot_ci <- function(x, fun = median) {
    n <- length(x)
    b <- replicate(n_boot, fun(x[sample.int(n, replace = TRUE)]))
    quantile(b, c(0.025, 0.975), names = FALSE)
  }
  groups <- list(All = rep(TRUE, nrow(res_df)), Low = res_df$expr_group %in% "Low", Medium = res_df$expr_group %in% "Medium", High = res_df$expr_group %in% "High")
  rows <- list()
  for (gn in names(groups)) {
    d  <- res_df[groups[[gn]], ]
    dg <- d$g[!is.na(d$g)]
    ci_l <- boot_ci(d$lfc); ci_g <- boot_ci(dg)
    n_down <- sum(d$lfc < -lfc_threshold); n_up <- sum(d$lfc > lfc_threshold)
    ci_dn <- binom.test(n_down, nrow(d))$conf.int
    rows[[gn]] <- data.frame(
      assay = "EUseq", comparison = "POLR2A_deg_8hr_vs_WT", expr_group = gn, n_genes = nrow(d),
      median_lfc = median(d$lfc), lfc_lo = ci_l[1], lfc_hi = ci_l[2],
      median_hedges_g = median(dg), g_lo = ci_g[1], g_hi = ci_g[2], median_z_repo = median(d$z_repo),
      n_down = n_down, n_up = n_up, frac_down = n_down / nrow(d), frac_down_lo = ci_dn[1], frac_down_hi = ci_dn[2], stringsAsFactors = FALSE)
    stopifnot(ci_l[1] <= median(d$lfc), median(d$lfc) <= ci_l[2], ci_g[1] <= median(dg), median(dg) <= ci_g[2])
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  stopifnot(res$n_genes[res$expr_group == "All"] == nrow(res_df))
  cat("\nmedian_lfc: median per-gene log2 fold change in scaled NoEU-subtracted RPKM, POLR2A-degraded minus Wildtype (negative = less nascent transcription)\n")
  cat("median_hedges_g: median of per-gene Hedges' g (LFC / within-group pooled SD, n = 3 vs 3); median_z_repo: median of the repo's LFC / SD(all six replicates)\n")
  cat("frac_down: fraction of genes with LFC < -", lfc_threshold, " (exact binomial CI); CIs for the medians are percentile bootstraps over genes\n\n", sep = "")
  print(res, digits = 3)
  write.csv(res, out_file, row.names = FALSE)
  cat("Wrote:", out_file, "\n")
  invisible(list(summary = res, genes = res_df))
}

euseq <- euseq_effect_sizes()