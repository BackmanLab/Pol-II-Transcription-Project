hic_compartment_effect_sizes <- function(
    # ds_all = WT, POLR2A-degraded and ActD at equal depth (same files as saddle_analysis.R / compartment_analysis.R, DS_GROUP = "ds_all")
  base_dir   = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp",
  # comp_bed: eigenvector bedgraphs written by run_compartments.sh (compartments/<label>/<label>_compartments.bedgraph)
  conditions = list(
    "Wildtype"   = list(hic = "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic",
                        comp_bed = "compartments/ds_all_WT/ds_all_WT_compartments.bedgraph"),
    "POLR2A deg" = list(hic = "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic",
                        comp_bed = "compartments/ds_all_Pol2/ds_all_Pol2_compartments.bedgraph"),
    "ActD 1hr"   = list(hic = "mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic",
                        comp_bed = "compartments/ds_all_ActD/ds_all_ActD_compartments.bedgraph")),
  control_label  = "Wildtype",
  treated_labels = c(POLR2A_deg_6hr = "POLR2A deg", ActD_1hr = "ActD 1hr"),
  # ENCODE HCT116 RSEM gene quantification (ENCFF435PHM), same table as the TAD and loop effect-size scripts
  expr_file  = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference/ENCFF435PHM.tsv",
  deseq_file = "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv",
  out_file   = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/HiC/HiC_compartment_effect_sizes.csv",
  # the compartment-strength plots are printed and also exported here as SVGs (the folder must already exist)
  plot_dir   = "/home/lmc0633/Transcription Paper/101426_revision_scripts/effect sizes/effect_sizes_tables_plots",
  panel_w_in = 3.2,        # SVG: width of one facet panel (inches)
  panel_h_in = 3.8,        # SVG: height of the figure (inches)
  resolution  = 500000,    # RESOLUTION
  saddle_bins = 20L,       # SADDLE_BINS
  corner_size = 5,         # corner_size in compartment_strength()
  straw_norm  = "SCALE",   # STRAW_NORM
  min_chroms  = 10,        # strata with fewer chromosomes than this are skipped
  n_boot = 10000,
  seed   = 102726) {
  
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  
  # ---- input checks ------------------------------------------------------------
  for (pk in c("GenomicRanges", "IRanges", "strawr", "parallel")) if (!requireNamespace(pk, quietly = TRUE)) stop("R package '", pk, "' is not installed")
  library(GenomicRanges)
  all_labels <- c(control_label, unname(treated_labels))
  if (any(!all_labels %in% names(conditions))) stop("conditions must be named: ", paste(all_labels, collapse = ", "))
  
  hic_paths  <- sapply(all_labels, function(l) file.path(base_dir, conditions[[l]]$hic))
  comp_paths <- sapply(all_labels, function(l) file.path(base_dir, conditions[[l]]$comp_bed))
  inputs <- c(hic_paths, comp_paths, expr_file, deseq_file)
  if (any(!file.exists(inputs))) stop("Input file(s) not found:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
  if (!dir.exists(plot_dir)) stop("Export directory not found: ", plot_dir)
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  cat("Eigenvector bedgraphs:\n"); print(data.frame(path = comp_paths, modified = format(file.mtime(comp_paths)), row.names = all_labels))
  
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK"); slurm_ntasks <- Sys.getenv("SLURM_NTASKS")
  n_cores <- if (nchar(slurm_cpus) > 0) as.integer(slurm_cpus) else if (nchar(slurm_ntasks) > 0) as.integer(slurm_ntasks) else 2L
  
  # ---- gene expression (WT TPM tertiles of expressed genes, as in saddle_analysis.R / compartment_analysis.R) ----
  BINS <- c("Low", "Medium", "High")
  expr_raw <- read_tsv(expr_file, show_col_types = FALSE)
  if (any(!c("gene_id", "TPM") %in% colnames(expr_raw))) stop("Columns 'gene_id' and 'TPM' required in ", basename(expr_file), "\nColumns present: ", paste(colnames(expr_raw), collapse = ", "))
  expr <- expr_raw %>%
    dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
    dplyr::group_by(gene_id) %>% dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")
  deseq_res <- read_csv(deseq_file, show_col_types = FALSE)
  for (cc in c("gene_id", "gene_seq_start", "gene_seq_end", "seq_name")) {
    if (!cc %in% colnames(deseq_res)) stop("Column '", cc, "' not found in ", basename(deseq_file), "\nColumns present: ", paste(colnames(deseq_res), collapse = ", "))
  }
  deseq_res <- deseq_res %>%
    dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
    left_join(expr, by = "gene_id")
  cat("DESeq genes with a TPM value:", sum(!is.na(deseq_res$TPM)), "/", nrow(deseq_res), "\n")
  if (mean(!is.na(deseq_res$TPM)) < 0.8) stop("Fewer than 80% of genes in the DESeq table matched the TPM table; check gene_id formats")
  expressed  <- deseq_res %>% filter(!is.na(TPM), TPM > 0) %>% pull(TPM)
  tpm_medium <- quantile(expressed, 1/3, na.rm = TRUE)
  tpm_high   <- quantile(expressed, 2/3, na.rm = TRUE)
  deseq_res <- deseq_res %>%
    dplyr::mutate(expr_bin = dplyr::case_when(is.na(TPM) | TPM == 0 ~ "Low", TPM <= tpm_medium ~ "Low", TPM <= tpm_high ~ "Medium", TRUE ~ "High"))
  
  # ---- eigenvectors (load_eigenvectors / check_eigen_sign copied from saddle_analysis.R) ----
  # chromosome names are stored without "chr" here and mapped to the names each .hic uses when reading contacts
  load_eigenvectors <- function(path, label) {
    df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "eigen"))
    if (!is.numeric(df$eigen) || !is.numeric(df$start)) stop("Non-numeric start/eigen column in ", path, "; first rows:\n", paste(utils::capture.output(print(head(df, 3))), collapse = "\n"))
    raw_chr <- unique(as.character(df$chr))
    df <- df %>% dplyr::mutate(chr = sub("^chr", "", as.character(chr))) %>% filter(!is.na(eigen))
    # run_compartments.sh turns every *_eigen.txt in the eigenvectors folder into bedgraph rows, so stale files with the other naming style would duplicate bins
    if (anyDuplicated(sprintf("%s_%.0f", df$chr, df$start)) > 0) stop(basename(path), " has duplicated (chromosome, start) rows after removing 'chr'; chromosome names present: ", paste(raw_chr, collapse = ", "),
                                                                      "\nA stale eigenvector set is probably in the eigenvectors folder; remove it and re-run run_compartments.sh for this condition")
    if (length(raw_chr) != 23) stop(basename(path), " has ", length(raw_chr), " chromosome names (expected 23): ", paste(raw_chr, collapse = ", "))
    if (any(df$start %% resolution != 0)) stop("Bedgraph starts in ", basename(path), " are not multiples of ", resolution, " bp")
    message("  ", label, ": ", nrow(df), " bins across ", length(unique(df$chr)), " chromosomes")
    df
  }
  check_eigen_sign <- function(eigen_df, deseq_res, sample_name) {
    gene_coords <- deseq_res %>%
      filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name), expr_bin == "High") %>%
      distinct(gene_id, .keep_all = TRUE)
    if (nrow(gene_coords) == 0) stop("No High-expression genes with coordinates in the DESeq table; cannot check eigenvector sign")
    gene_gr <- GRanges(seqnames = paste0("chr", sub("^chr", "", gene_coords$seq_name)), ranges = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end))
    bin_gr  <- GRanges(seqnames = paste0("chr", eigen_df$chr), ranges = IRanges(eigen_df$start + 1, eigen_df$end))
    hits      <- findOverlaps(bin_gr, gene_gr)
    high_bins <- unique(queryHits(hits))
    if (length(high_bins) == 0) stop("No eigenvector bin overlaps a High-expression gene for ", sample_name, "; chromosome names do not match between bedgraph and DESeq table")
    mean_high <- mean(eigen_df$eigen[high_bins], na.rm = TRUE)
    pct_pos   <- mean(eigen_df$eigen[high_bins] > 0, na.rm = TRUE) * 100
    message("  ", sample_name, ": mean eigen (High-expr bins) = ", round(mean_high, 4), "  (", round(pct_pos, 1), "% positive)")
    if (mean_high < 0) {
      message("  *** Sign flipped - correcting ***")
      eigen_df <- eigen_df %>% dplyr::mutate(eigen = -eigen)
    }
    eigen_df
  }
  message("\nLoading eigenvectors...")
  eigen_list <- lapply(setNames(all_labels, all_labels), function(l) load_eigenvectors(comp_paths[[l]], l))
  message("\nChecking eigenvector signs...")
  eigen_list <- lapply(setNames(all_labels, all_labels), function(l) check_eigen_sign(eigen_list[[l]], deseq_res, l))
  
  # ---- expression bins of the 500 kb genomic bins (assign_expr_bins_to_genome copied from compartment_analysis.R) ----
  assign_expr_bins_to_genome <- function(eigen_df, deseq_res) {
    gene_coords <- deseq_res %>% filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name)) %>% distinct(gene_id, .keep_all = TRUE)
    gene_gr <- GRanges(seqnames = paste0("chr", sub("^chr", "", gene_coords$seq_name)), ranges = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end), TPM = gene_coords$TPM)
    bin_gr  <- GRanges(seqnames = paste0("chr", eigen_df$chr), ranges = IRanges(eigen_df$start + 1, eigen_df$end))
    hits <- findOverlaps(bin_gr, gene_gr)
    if (length(hits) == 0) stop("No eigenvector bin overlaps a gene; chromosome names do not match between bedgraph and DESeq table")
    bin_tpm <- data.frame(bin_idx = queryHits(hits), TPM = gene_gr$TPM[subjectHits(hits)]) %>%
      dplyr::group_by(bin_idx) %>% dplyr::summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")
    bin_expressed  <- bin_tpm %>% filter(!is.na(mean_TPM), mean_TPM > 0) %>% pull(mean_TPM)
    bin_tpm_medium <- quantile(bin_expressed, 1/3, na.rm = TRUE)
    bin_tpm_high   <- quantile(bin_expressed, 2/3, na.rm = TRUE)
    message("  Bin-level TPM thresholds: Low<=", round(bin_tpm_medium, 2), "  Med<=", round(bin_tpm_high, 2), "  High>", round(bin_tpm_high, 2))
    eigen_df %>%
      dplyr::mutate(bin_idx = row_number()) %>%
      left_join(bin_tpm, by = "bin_idx") %>%
      dplyr::mutate(expr_bin = dplyr::case_when(is.na(mean_TPM) | mean_TPM == 0 ~ "Low", mean_TPM <= bin_tpm_medium ~ "Low", mean_TPM <= bin_tpm_high ~ "Medium", TRUE ~ "High"),
                    expr_bin = factor(expr_bin, levels = BINS)) %>%
      dplyr::select(-bin_idx)
  }
  message("\nAssigning expression bins to 500kb bins...")
  eigen_list <- lapply(setNames(all_labels, all_labels), function(l) { message(l, ":"); assign_expr_bins_to_genome(eigen_list[[l]], deseq_res) })
  for (l in all_labels) { cat(l, "expression-bin counts:\n"); print(table(eigen_list[[l]]$expr_bin)) }
  
  # ---- saddle (compute_saddle_chr copied from saddle_analysis.R; compartment_strength copied, unchanged) ----
  compute_saddle_chr <- function(hic_file, hic_chrom, chrom, eigen_chr, n_bins, res, norm) {
    library(strawr); library(dplyr)
    eigen_valid <- eigen_chr %>% filter(!is.na(eigen)) %>% dplyr::arrange(eigen) %>% dplyr::mutate(quantile_bin = ntile(eigen, n_bins))
    if (nrow(eigen_valid) < n_bins) return(NULL)
    raw <- tryCatch(strawr::straw(norm, hic_file, as.character(hic_chrom), as.character(hic_chrom), "BP", res), error = function(e) NULL)
    if (is.null(raw) || nrow(raw) == 0) return(NULL)
    raw <- raw %>% dplyr::mutate(bin1 = x / res, bin2 = y / res)
    bin_to_quantile <- setNames(eigen_valid$quantile_bin, as.character(eigen_valid$start / res))
    raw <- raw %>%
      dplyr::mutate(q1 = bin_to_quantile[as.character(bin1)], q2 = bin_to_quantile[as.character(bin2)], dist_bins = abs(bin2 - bin1)) %>%
      filter(!is.na(q1), !is.na(q2))
    if (nrow(raw) == 0) return(NULL)
    expected_by_dist <- raw %>% dplyr::group_by(dist_bins) %>% dplyr::summarise(expected = mean(counts, na.rm = TRUE), .groups = "drop")
    raw <- raw %>% left_join(expected_by_dist, by = "dist_bins") %>% dplyr::mutate(oe = counts / expected) %>% filter(!is.na(oe), is.finite(oe))
    if (nrow(raw) == 0) return(NULL)
    saddle_upper <- raw %>% dplyr::group_by(q1, q2) %>% dplyr::summarise(mean_oe = mean(oe, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
    bind_rows(saddle_upper, saddle_upper %>% dplyr::rename(q1 = q2, q2 = q1)) %>%
      distinct(q1, q2, .keep_all = TRUE) %>% dplyr::mutate(chr = as.character(chrom))
  }
  compartment_strength <- function(saddle_df, n_bins = saddle_bins, corner = corner_size) {
    lo <- 1:corner; hi <- (n_bins - corner + 1):n_bins
    get_mean <- function(q1r, q2r) saddle_df %>% filter(q1 %in% q1r, q2 %in% q2r) %>% summarise(m = mean(mean_oe, na.rm = TRUE)) %>% pull(m)
    aa <- get_mean(hi, hi); bb <- get_mean(lo, lo); ab <- get_mean(lo, hi); ba <- get_mean(hi, lo)
    list(AA = aa, BB = bb, AB = ab, BA = ba, CS = (aa + bb) / (ab + ba))
  }
  
  STRATA <- c("All", BINS)   # All = every bin (same saddle function without the expression filter); Low/Medium/High = expression tertiles of the 500 kb bins, as in compartment_analysis.R
  per_chr <- list(); genome_cs <- list(); sad_mat <- list()   # sad_mat: chromosome x (q1,q2) matrix of mean O/E per condition and stratum
  for (l in all_labels) {
    hic_names   <- strawr::readHicChroms(hic_paths[[l]])$name
    hic_name_of <- setNames(hic_names, sub("^chr", "", hic_names))
    e_all <- eigen_list[[l]]
    if (any(!unique(e_all$chr) %in% names(hic_name_of))) stop("Chromosomes in the ", l, " bedgraph not found in the .hic: ", paste(setdiff(unique(e_all$chr), names(hic_name_of)), collapse = ", "),
                                                              "\n.hic chromosomes: ", paste(hic_names, collapse = ", "))
    for (st in STRATA) {
      e_sub <- if (st == "All") e_all else e_all %>% filter(as.character(expr_bin) == st)
      message(l, " | ", st, ": ", nrow(e_sub), " bins")
      if (nrow(e_sub) < saddle_bins * 2) next
      chroms <- unique(e_sub$chr)
      res_list <- parallel::mclapply(setNames(chroms, chroms), function(ch) {
        e_chr <- e_sub %>% filter(chr == ch)
        if (nrow(e_chr) < saddle_bins) return(NULL)
        compute_saddle_chr(hic_paths[[l]], hic_name_of[[ch]], ch, e_chr, saddle_bins, resolution, straw_norm)
      }, mc.cores = n_cores, mc.preschedule = FALSE)
      bad <- vapply(res_list, function(r) inherits(r, "try-error"), logical(1))
      if (any(bad)) stop("Saddle computation failed for ", l, " | ", st, " on chromosomes: ", paste(names(res_list)[bad], collapse = ", "))
      res_list <- Filter(Negate(is.null), res_list)
      if (length(res_list) == 0) next
      key <- paste(l, st)
      per_chr[[key]] <- dplyr::bind_rows(lapply(res_list, function(s) data.frame(sample = l, stratum = st, chr = s$chr[1], CS = compartment_strength(s)$CS, stringsAsFactors = FALSE)))
      m <- matrix(NA_real_, nrow = length(res_list), ncol = saddle_bins^2, dimnames = list(names(res_list), NULL))
      for (ch in names(res_list)) { s <- res_list[[ch]]; m[ch, (s$q1 - 1) * saddle_bins + s$q2] <- s$mean_oe }
      sad_mat[[key]] <- m
      # genome-wide saddle exactly as saddle_analysis.R: mean O/E per (q1, q2) across chromosomes, then CS
      gs <- dplyr::bind_rows(res_list) %>% dplyr::group_by(q1, q2) %>% dplyr::summarise(mean_oe = mean(mean_oe, na.rm = TRUE), .groups = "drop")
      cs <- compartment_strength(gs)
      genome_cs[[key]] <- data.frame(sample = l, stratum = st, n_chroms = length(res_list), AA = cs$AA, BB = cs$BB, AB = cs$AB, BA = cs$BA, CS = cs$CS, stringsAsFactors = FALSE)
    }
  }
  per_chr   <- dplyr::bind_rows(per_chr)
  genome_cs <- dplyr::bind_rows(genome_cs)
  cat("\nGenome-wide compartment strength (CS = (AA+BB)/(AB+BA), saddle_analysis.R definition):\n"); print(genome_cs, digits = 3)
  
  # ---- structure self-checks on live data ----------------------------------------
  n_all <- genome_cs %>% filter(stratum == "All") %>% pull(n_chroms)
  if (length(n_all) != length(all_labels) || any(n_all != 23)) stop("Expected 23 chromosomes (1-22, X) in the All stratum for every condition; got: ", paste(n_all, collapse = ", "))
  if (any(!is.finite(genome_cs$CS)) || any(genome_cs$CS <= 1)) stop("Genome-wide CS is not finite or not > 1 in some condition/stratum; the saddle orientation or quantile assignment is wrong")
  if (any(!is.finite(per_chr$CS[per_chr$stratum == "All"]))) stop("Non-finite per-chromosome CS in the All stratum")
  # CS recomputed from the stored chromosome x (q1,q2) matrices must equal the dplyr genome-wide CS above (the bootstrap below uses the matrices)
  lo <- 1:corner_size; hi <- (saddle_bins - corner_size + 1):saddle_bins
  cell <- function(q1r, q2r) as.vector(outer(q1r - 1, q2r, function(a, b) a * saddle_bins + b))
  idx_hh <- cell(hi, hi); idx_ll <- cell(lo, lo); idx_lh <- cell(lo, hi); idx_hl <- cell(hi, lo)
  cs_from_vec <- function(v) (mean(v[idx_hh], na.rm = TRUE) + mean(v[idx_ll], na.rm = TRUE)) / (mean(v[idx_lh], na.rm = TRUE) + mean(v[idx_hl], na.rm = TRUE))
  genome_vec <- function(m) colMeans(m, na.rm = TRUE)
  for (k in names(sad_mat)) {
    gk <- genome_cs[paste(genome_cs$sample, genome_cs$stratum) == k, ]
    stopifnot(isTRUE(all.equal(cs_from_vec(genome_vec(sad_mat[[k]])), gk$CS)))
  }
  cat("Self-check passed: 23 chromosomes per condition (All stratum), genome-wide CS finite and > 1, matrix CS equals saddle_analysis.R CS\n")
  
  # ---- effect sizes -------------------------------------------------------------------
  # pct_change_cs = CS_treated / CS_control - 1, CS = genome-wide CS (chromosome-averaged saddle, as in saddle_analysis.R).
  # hedges_g_paired = Hedges' g_z on per-chromosome CS (mean / SD of the per-chromosome differences, small-sample corrected, df = n - 1).
  # Both: percentile bootstrap, chromosomes resampled jointly across conditions.
  hedges_gz <- function(diffs) {
    n <- length(diffs); s <- sd(diffs)
    if (!is.finite(s) || s == 0) stop("SD of the per-chromosome differences is zero or undefined")
    (mean(diffs) / s) * (1 - 3 / (4 * (n - 1) - 1))
  }
  set.seed(seed)
  rows <- list(); chk_done <- FALSE
  for (st in STRATA) {
    d <- per_chr %>% filter(stratum == st, is.finite(CS))
    common <- Reduce(intersect, lapply(all_labels, function(l) d$chr[d$sample == l]))   # same chromosomes in every condition
    if (length(common) < min_chroms) { message("Skipping stratum ", st, ": ", length(common), " chromosomes (< ", min_chroms, ")"); next }
    for (tn in names(treated_labels)) {
      tl <- treated_labels[[tn]]
      mc <- sad_mat[[paste(control_label, st)]][common, , drop = FALSE]
      mt <- sad_mat[[paste(tl, st)]][common, , drop = FALSE]
      cc <- setNames(d$CS[d$sample == control_label], d$chr[d$sample == control_label])[common]
      ct <- setNames(d$CS[d$sample == tl], d$chr[d$sample == tl])[common]
      n <- length(common)
      cs_c <- cs_from_vec(genome_vec(mc)); cs_t <- cs_from_vec(genome_vec(mt))
      pct <- cs_t / cs_c - 1
      gz  <- hedges_gz(ct - cc)
      if (!chk_done) {   # self-check on live data: g_z must equal the paired t-statistic conversion
        chk_t <- unname(t.test(ct, cc, paired = TRUE)$statistic)
        stopifnot(isTRUE(all.equal(gz, chk_t / sqrt(n) * (1 - 3 / (4 * (n - 1) - 1)))))
        cat("Self-check passed: g_z matches the paired t-statistic conversion\n")
        chk_done <- TRUE
      }
      bs <- replicate(n_boot, {
        i <- sample.int(n, replace = TRUE)
        c(cs_from_vec(genome_vec(mt[i, , drop = FALSE])) / cs_from_vec(genome_vec(mc[i, , drop = FALSE])) - 1,
          hedges_gz(ct[i] - cc[i]))
      })
      ci_p <- quantile(bs[1, ], c(0.025, 0.975), names = FALSE, na.rm = TRUE)
      ci_g <- quantile(bs[2, ], c(0.025, 0.975), names = FALSE, na.rm = TRUE)
      if (pct < ci_p[1] || pct > ci_p[2] || gz < ci_g[1] || gz > ci_g[2]) warning("Point estimate outside its bootstrap CI for ", tn, " | ", st)
      rows[[length(rows) + 1]] <- data.frame(
        assay = "HiC_compartment_strength", comparison = paste0(tn, "_vs_WT"), expr_bin = st, n_chroms = n,
        cs_ctrl = cs_c, cs_trt = cs_t, pct_change_cs = 100 * pct, pct_lo = 100 * ci_p[1], pct_hi = 100 * ci_p[2],
        hedges_g_paired = gz, g_lo = ci_g[1], g_hi = ci_g[2], stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  cat("\npct_change_cs: % change in genome-wide compartment strength, treated vs Wildtype; negative = weaker compartmentalization\n")
  cat("hedges_g_paired: paired Hedges' g on per-chromosome CS (SD of the differences); chromosomes are the resampling unit for both CIs\n\n")
  print(res, digits = 3)
  write.csv(res, out_file, row.names = FALSE)
  cat("Wrote:", out_file, "\n")
  
  # ---- compartment-strength plots: plot_cs_faceted() and plot_corner_scores() of compartment_analysis.R (section 15) ----
  # Same CS definition as above; bars/points = genome-wide CS and per-corner log2 O/E (genome_cs) for every condition, one facet per expression stratum.
  # Deviations from the repo functions: COND_LEVELS / COND_COLORS are all_labels / cond_colors (the repo's colours), the "All" stratum is added as a facet,
  # and plot_corner_scores takes ncol_wrap (the repo fixes 3 facet columns).
  cond_colors <- setNames(c("#88CCEE", "#CC6677", "#44AA99"), c("Wildtype", "POLR2A deg", "ActD 1hr"))[all_labels]
  if (anyNA(names(cond_colors))) stop("No plot colour defined for: ", paste(setdiff(all_labels, c("Wildtype", "POLR2A deg", "ActD 1hr")), collapse = ", "))
  cs_df <- genome_cs %>%
    dplyr::mutate(category = stratum, log2_AA = log2(AA), log2_BB = log2(BB), log2_AB = log2(AB), log2_BA = log2(BA))
  categories <- STRATA[STRATA %in% cs_df$category]
  stopifnot(isTRUE(all.equal(cs_df$CS, (cs_df$AA + cs_df$BB) / (cs_df$AB + cs_df$BA))), all(is.finite(cs_df$log2_AA + cs_df$log2_BB + cs_df$log2_AB + cs_df$log2_BA)))
  cat("Self-check passed: plotted CS equals (AA+BB)/(AB+BA) of the genome-wide saddle; all corner values finite and > 0\n")
  plot_cs_faceted <- function(cs_df, categories, title = "Compartment strength by category", ncol_wrap = 3) {
    plot_df <- cs_df %>%
      filter(category %in% categories) %>%
      dplyr::mutate(category = factor(category, levels = categories),
                    sample = factor(.data$sample, levels = all_labels))
    ggplot(plot_df, aes(x = sample, y = CS, fill = sample, color = sample)) +
      geom_col(alpha = 0.8, width = 0.6) +
      geom_text(aes(label = round(CS, 2)), vjust = -0.4, size = 3, fontface = "bold") +
      scale_fill_manual(values = cond_colors, guide = "none") +
      scale_color_manual(values = cond_colors, guide = "none") +
      scale_y_continuous(name = "Compartment strength\n(AA+BB)/(AB+BA)", expand = expansion(mult = c(0, 0.2))) +
      facet_wrap(~category, ncol = ncol_wrap, scales = "free_y") +
      labs(title = title, x = NULL) +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            strip.text = element_text(face = "bold", size = 10),
            strip.background = element_blank(),
            axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
            plot.margin = margin(6, 6, 6, 6))
  }
  plot_corner_scores <- function(cs_df, categories, title = "Per-corner log2 O/E scores", ncol_wrap = 3) {
    plot_df <- cs_df %>%
      filter(category %in% categories) %>%
      dplyr::select(sample, category, log2_AA, log2_BB, log2_AB, log2_BA) %>%
      tidyr::pivot_longer(cols = c(log2_AA, log2_BB, log2_AB, log2_BA), names_to = "corner", values_to = "log2_oe") %>%
      dplyr::mutate(corner = gsub("log2_", "", corner),
                    corner = factor(corner, levels = c("AA", "BB", "AB", "BA")),
                    category = factor(category, levels = categories),
                    sample = factor(.data$sample, levels = all_labels))
    ggplot(plot_df, aes(x = corner, y = log2_oe, color = sample, group = sample)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
      geom_line(linewidth = 0.8) + geom_point(size = 3) +
      geom_text(aes(label = round(log2_oe, 2)), vjust = -0.8, size = 2.5, show.legend = FALSE) +
      scale_color_manual(values = cond_colors, name = NULL) +
      scale_y_continuous(name = "log2 O/E") +
      facet_wrap(~category, ncol = ncol_wrap, scales = "free_y") +
      labs(title = title, x = "Compartment interaction") +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            strip.text = element_text(face = "bold", size = 10),
            strip.background = element_blank(),
            legend.position = "bottom",
            axis.text.x = element_text(size = 9, face = "bold"),
            plot.margin = margin(6, 6, 6, 6))
  }
  save_svg <- function(p, name) {   # printed to the graphics device and exported as SVG
    print(p)
    f <- file.path(plot_dir, name)
    w <- panel_w_in * length(categories)
    svg(f, width = w, height = panel_h_in)
    print(p)
    dev.off()
    message(sprintf("Wrote %s  (%.1f x %.1f in)", f, w, panel_h_in))
    invisible(f)
  }
  p_cs     <- plot_cs_faceted(cs_df, categories, title = "Compartment strength by expression bin (All = every 500 kb bin)", ncol_wrap = length(categories))
  p_corner <- plot_corner_scores(cs_df, categories, title = "Per-corner log2 O/E by expression bin (All = every 500 kb bin)", ncol_wrap = length(categories))
  save_svg(p_cs,     "HiC_compartment_strength_by_expression_bin.svg")
  save_svg(p_corner, "HiC_compartment_corner_scores_by_expression_bin.svg")
  
  invisible(list(summary = res, per_chr = per_chr, genome_cs = genome_cs, plots = list(cs = p_cs, corner = p_corner)))
}

hic_compartment <- hic_compartment_effect_sizes()