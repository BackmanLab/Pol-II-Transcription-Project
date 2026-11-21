hic_tad_insulation_effect_sizes <- function(
    # ds_all = WT, POLR2A-degraded and ActD at equal depth (same files as insulation_analysis.R, DS_GROUP = "ds_all")
  base_dir   = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp",
  # insulation scores: loaded from ins_rds if it exists, otherwise computed from hic_files (as in
  # insulation_analysis.R, DS_GROUP = "ds_all") and saved to ins_rds / ins_rds_raw
  ins_rds    = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/insulation/insulation_zscored_ds_all_SCALE.rds",
  ins_rds_raw = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/insulation/insulation_raw_ds_all_SCALE.rds",
  hic_files  = c("Wildtype"   = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic",
                 "POLR2A deg" = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic",
                 "ActD 1hr"   = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic"),
  ins_window   = 10,       # INS_WINDOW (bins)
  straw_norm   = "SCALE",  # STRAW_NORM
  tad_files  = c(wt   = "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30_contact_domains/10000_blocks.bedpe",
                 pol2 = "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30_contact_domains/10000_blocks.bedpe"),
  # ENCODE HCT116 RSEM gene quantification (ENCFF435PHM, GENCODE V29 IDs), same table as the EU-seq hg38 script
  expr_file  = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference/ENCFF435PHM.tsv",
  deseq_file = "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv",
  # SCREEN Registry-V4 cCRE bed with CTCF data (column 10 = cCRE class; CA-CTCF defines CTCF-anchored TADs)
  ctcf_file  = "/projects/b1042/BackmanLab/Lucas/loop_chromatin_regression/data/classification/ENCFF431JDU_ENCFF964OOU_ENCFF787LMI_ENCFF388PVO.bed",
  out_file   = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/HiC/HiC_TAD_insulation_effect_sizes.csv",
  # the profile plots are printed and also exported here as SVGs (one per TAD positions x CTCF status; the folder must already exist)
  plot_dir   = "/home/lmc0633/Transcription Paper/101426_revision_scripts/effect sizes/effect_sizes_tables_plots",
  panel_w_in = 5,          # SVG: width of one profile panel (inches)
  panel_h_in = 3.4,        # SVG: height of one profile panel (inches)
  control_label   = "Wildtype",
  treated_labels   = c(POLR2A_deg_6hr = "POLR2A deg", ActD_1hr = "ActD 1hr"),
  resolution   = 25000,    # RESOLUTION in insulation_analysis.R
  tad_min_kb   = 50,       # TAD_MIN_KB
  tad_max_kb   = 5000,     # TAD_MAX_KB
  gene_tol_bp  = 5000,     # GENE_TOL_BP
  ctcf_tol_bp  = 10000,    # CTCF_TOL_BP
  min_tads     = 10,
  flank_kb     = 500,      # INS_FLANK_KB (profile plots)
  roll_window  = 5,        # ROLL_WINDOW (profile plots)
  n_boot = 10000,
  seed   = 102726) {
  
  library(dplyr)
  library(tidyr)
  library(readr)
  
  # ---- input checks ------------------------------------------------------------
  for (pk in c("GenomicRanges", "IRanges", "ggplot2", "patchwork", "zoo")) if (!requireNamespace(pk, quietly = TRUE)) stop("R package '", pk, "' is not installed")
  library(ggplot2)
  library(GenomicRanges)
  tad_paths <- file.path(base_dir, tad_files)
  inputs <- c(tad_paths, expr_file, deseq_file, ctcf_file, if (!file.exists(ins_rds) && !file.exists(ins_rds_raw)) hic_files)
  if (any(!file.exists(inputs))) stop("Input file(s) not found:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
  if (!dir.exists(plot_dir)) stop("Export directory not found: ", plot_dir)
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  names(tad_paths) <- names(tad_files)
  
  # ---- insulation scores (cached z-scored RDS from insulation_analysis.R) ------------
  all_labels <- c(control_label, unname(treated_labels))
  if (!file.exists(ins_rds)) {
    # ---- compute insulation (function copied from insulation_analysis.R, compute_genome_insulation) ----
    for (pk in c("strawr", "data.table", "parallel")) if (!requireNamespace(pk, quietly = TRUE)) stop("R package '", pk, "' is not installed")
    if (any(!all_labels %in% names(hic_files))) stop("hic_files must be named: ", paste(all_labels, collapse = ", "))
    slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK"); slurm_ntasks <- Sys.getenv("SLURM_NTASKS")
    n_cores <- max(1L, (if (nchar(slurm_cpus) > 0) as.integer(slurm_cpus) else if (nchar(slurm_ntasks) > 0) as.integer(slurm_ntasks) else 2L) - 1L)
    compute_genome_insulation <- function(hic_file, sample_name) {
      .compute_chr <- function(hic_file, chrom, resolution, window, norm) {
        library(strawr); library(data.table)
        raw <- tryCatch(strawr::straw(norm, hic_file, as.character(chrom), as.character(chrom), "BP", resolution), error = function(e) NULL)
        if (is.null(raw) || nrow(raw) == 0) return(NULL)
        dt <- data.table::as.data.table(raw)
        dt[, bin1 := as.integer(x / resolution)]
        dt[, bin2 := as.integer(y / resolution)]
        all_bins <- sort(unique(c(dt$bin1, dt$bin2)))
        if (length(all_bins) < 2 * window + 1) return(NULL)
        dt[, b_lo := pmax(bin1, bin2 - window) + 1L]
        dt[, b_hi := pmin(bin1 + window - 1L, bin2)]
        dt_contrib <- dt[b_lo <= b_hi]
        if (nrow(dt_contrib) == 0) return(NULL)
        dt_exp   <- dt_contrib[, .(b = seq.int(b_lo, b_hi), counts = counts), by = seq_len(nrow(dt_contrib))]
        win_sums <- dt_exp[, .(window_sum = sum(counts, na.rm = TRUE)), by = b]
        result   <- merge(data.table(b = all_bins), win_sums, by = "b", all.x = TRUE)
        data.frame(chr = as.character(chrom), start = result$b * resolution, end = result$b * resolution + resolution,
                   bin_idx = result$b, window_sum = result$window_sum, stringsAsFactors = FALSE)
      }
      hic_chroms <- strawr::readHicChroms(hic_file)
      bare <- sub("^chr", "", hic_chroms$name)   # file may name chromosomes "1" or "chr1"
      chroms <- hic_chroms$name[grepl("^[0-9]+$", bare) | bare == "X"]
      if (length(chroms) == 0) stop("No autosomes/X found in ", hic_file, "\nChromosome names in the file: ", paste(hic_chroms$name, collapse = ", "))
      message("Computing insulation scores: ", sample_name, "  norm=", straw_norm, "  res=", resolution / 1000, "kb  window=", ins_window, " bins  cores=", n_cores)
      res_list <- parallel::mclapply(chroms, function(chrom) .compute_chr(hic_file, chrom, resolution, ins_window, straw_norm),
                                     mc.cores = n_cores, mc.preschedule = FALSE)
      message("  ", sum(!sapply(res_list, is.null)), "/", length(chroms), " chromosomes OK")
      if (all(sapply(res_list, is.null))) stop("strawr returned no data for any chromosome in ", hic_file, " (norm ", straw_norm, ", ", resolution, " bp)")
      if (any(sapply(res_list, is.null))) stop("Insulation failed for ", sum(sapply(res_list, is.null)), " chromosome(s) of ", sample_name)
      dplyr::bind_rows(res_list) %>% dplyr::mutate(ins_score = window_sum, sample = sample_name)
    }
    dir.create(dirname(ins_rds), showWarnings = FALSE, recursive = TRUE)
    if (file.exists(ins_rds_raw)) {
      message("Loading cached raw insulation scores: ", ins_rds_raw)
      ins_all <- readRDS(ins_rds_raw)
    } else {
      ins_all <- dplyr::bind_rows(lapply(all_labels, function(l) compute_genome_insulation(hic_files[[l]], l))) %>%
        dplyr::mutate(sample = factor(sample, levels = all_labels))
      saveRDS(ins_all, ins_rds_raw)
      message("Raw scores saved: ", ins_rds_raw)
    }
    ins_all <- ins_all %>%
      dplyr::group_by(sample, chr) %>%
      dplyr::mutate(ins_zscore = (ins_score - mean(ins_score, na.rm = TRUE)) / sd(ins_score, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sample = factor(sample, levels = all_labels))
    saveRDS(ins_all, ins_rds)
    message("Z-scored insulation saved: ", ins_rds)
  }
  ins <- readRDS(ins_rds)
  need <- c("chr", "start", "sample", "ins_zscore")
  if (any(!need %in% colnames(ins))) stop("Columns missing in insulation RDS: ", paste(setdiff(need, colnames(ins)), collapse = ", "),
                                          "\nColumns present: ", paste(colnames(ins), collapse = ", "))
  ins$sample <- as.character(ins$sample)
  ins$chr    <- sub("^chr", "", as.character(ins$chr))
  if (any(!all_labels %in% ins$sample)) stop("Sample labels not in insulation RDS: ", paste(setdiff(all_labels, ins$sample), collapse = ", "),
                                             "\nSamples present: ", paste(unique(ins$sample), collapse = ", "))
  cat("Insulation bins per sample:\n"); print(table(ins$sample))
  ins_key <- lapply(setNames(all_labels, all_labels), function(s) {
    d <- ins[ins$sample == s, ]
    if (anyDuplicated(paste(d$chr, d$start)) > 0) stop("Duplicated (chr, start) rows for sample ", s)
    setNames(d$ins_zscore, paste(d$chr, d$start))
  })
  
  # ---- expression bins (same definition as insulation_analysis.R: WT TPM tertiles of expressed genes) ----
  BINS <- c("Low", "Medium", "High")
  expr_raw <- read_tsv(expr_file, show_col_types = FALSE)
  if (any(!c("gene_id", "TPM") %in% colnames(expr_raw))) stop("Columns 'gene_id' and 'TPM' required in ", basename(expr_file), "\nColumns present: ", paste(colnames(expr_raw), collapse = ", "))
  expr <- expr_raw %>%
    dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%      # drop Ensembl version / _PAR_Y suffix
    dplyr::group_by(gene_id) %>% dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")
  deseq_res <- read_csv(deseq_file, show_col_types = FALSE)
  for (cc in c("gene_id", "symbol", "gene_seq_start", "gene_seq_end", "seq_name")) {
    if (!cc %in% colnames(deseq_res)) stop("Column '", cc, "' not found in ", basename(deseq_file), "\nColumns present: ", paste(colnames(deseq_res), collapse = ", "))
  }
  deseq_res <- deseq_res %>%
    dplyr::mutate(gene_id = sub("\\..*", "", gene_id),
                  symbol  = dplyr::case_when(!is.na(symbol) & symbol != "" ~ symbol, TRUE ~ gene_id)) %>%
    left_join(expr, by = "gene_id")
  cat("DESeq genes with a TPM value:", sum(!is.na(deseq_res$TPM)), "/", nrow(deseq_res), "\n")
  if (mean(!is.na(deseq_res$TPM)) < 0.8) stop("Fewer than 80% of genes in the DESeq table matched the TPM table; check gene_id formats")
  expressed  <- deseq_res %>% filter(!is.na(TPM), TPM > 0) %>% pull(TPM)
  tpm_medium <- quantile(expressed, 1/3, na.rm = TRUE)
  tpm_high   <- quantile(expressed, 2/3, na.rm = TRUE)
  cat("\nExpression bin thresholds (WT TPM tertiles): Low <=", round(tpm_medium, 2), "| Medium <=", round(tpm_high, 2), "| High >", round(tpm_high, 2), "\n")
  deseq_res <- deseq_res %>%
    dplyr::mutate(expr_bin = dplyr::case_when(is.na(TPM) | TPM == 0 ~ "Low", TPM <= tpm_medium ~ "Low", TPM <= tpm_high ~ "Medium", TRUE ~ "High"),
                  expr_bin = factor(expr_bin, levels = BINS))
  gene_coords <- deseq_res %>%
    filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name)) %>%
    distinct(gene_id, .keep_all = TRUE)
  gene_gr <- GRanges(seqnames = sub("^chr", "", as.character(gene_coords$seq_name)),
                     ranges = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end),
                     symbol = gene_coords$symbol, gene_id = gene_coords$gene_id,
                     expr_bin = as.character(gene_coords$expr_bin), TPM = gene_coords$TPM)
  
  # ---- TADs (functions copied from insulation_analysis.R) -------------------------------
  load_arrowhead <- function(path, min_size_kb = tad_min_kb, max_size_kb = tad_max_kb) {
    all_lines  <- readLines(path)
    data_lines <- all_lines[!grepl("^#", all_lines)]
    df <- read.table(text = paste(data_lines, collapse = "\n"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
    colnames(df)[1:6] <- c("chr1", "x1", "x2", "chr2", "y1", "y2")
    df %>% dplyr::mutate(tad_size_kb = (x2 - x1) / 1000, chr1 = sub("^chr", "", as.character(chr1))) %>%
      filter(!is.na(x1), !is.na(x2), x2 > x1, tad_size_kb >= min_size_kb, tad_size_kb <= max_size_kb)
  }
  link_tads_anchor <- function(tads_df, gene_gr, tol_bp = gene_tol_bp) {
    tads_df <- tads_df %>% dplyr::mutate(tad_idx = row_number())
    gene_gr_exp <- trim(resize(gene_gr, width(gene_gr) + 2 * tol_bp, fix = "center"))
    left_gr  <- GRanges(tads_df$chr1, IRanges(pmax(0, tads_df$x1 - tol_bp), tads_df$x1 + tol_bp))
    right_gr <- GRanges(tads_df$chr1, IRanges(pmax(0, tads_df$x2 - tol_bp), tads_df$x2 + tol_bp))
    hits_left  <- findOverlaps(left_gr,  gene_gr_exp)
    hits_right <- findOverlaps(right_gr, gene_gr_exp)
    all_pairs <- bind_rows(data.frame(tad_idx = queryHits(hits_left),  gene_idx = subjectHits(hits_left)),
                           data.frame(tad_idx = queryHits(hits_right), gene_idx = subjectHits(hits_right))) %>%
      distinct(tad_idx, gene_idx) %>%
      dplyr::mutate(expr_bin = gene_gr$expr_bin[gene_idx], TPM = gene_gr$TPM[gene_idx])
    tad_annotations <- all_pairs %>% dplyr::group_by(tad_idx) %>%
      dplyr::slice_max(TPM, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::select(tad_idx, expr_bin, TPM)
    tads_df %>% left_join(tad_annotations, by = "tad_idx") %>% dplyr::mutate(expr_bin = replace_na(expr_bin, "Low"))
  }
  n_fields <- length(strsplit(readLines(ctcf_file, n = 1), "\t")[[1]])
  if (n_fields < 10) stop("cCRE bed has ", n_fields, " columns; column 10 (cCRE class) is required: ", ctcf_file)
  ctcf_sites <- readr::read_tsv(ctcf_file, col_names = FALSE, show_col_types = FALSE,
                                col_types = readr::cols_only(X1 = readr::col_character(), X2 = readr::col_double(),
                                                             X3 = readr::col_double(), X10 = readr::col_character())) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, ccre_type = X10)
  cat("cCRE classes in the bed:\n"); print(sort(table(ctcf_sites$ccre_type), decreasing = TRUE))
  ctcf_sites <- ctcf_sites %>% filter(ccre_type == "CA-CTCF")
  if (nrow(ctcf_sites) == 0) stop("No 'CA-CTCF' cCREs in ", basename(ctcf_file), "; this bed cannot define CTCF-anchored TADs (classes printed above)")
  ctcf_gr <- GRanges(sub("^chr", "", ctcf_sites$chr), IRanges(ctcf_sites$start, ctcf_sites$end))
  cat("CA-CTCF sites:", length(ctcf_gr), "\n")
  filter_ctcf_tads <- function(tads_df, ctcf_gr, tol_bp = ctcf_tol_bp) {
    left_gr  <- GRanges(tads_df$chr1, IRanges(pmax(0, tads_df$x1 - tol_bp), tads_df$x1 + tol_bp))
    right_gr <- GRanges(tads_df$chr1, IRanges(pmax(0, tads_df$x2 - tol_bp), tads_df$x2 + tol_bp))
    both_ctcf <- intersect(unique(queryHits(findOverlaps(left_gr, ctcf_gr))), unique(queryHits(findOverlaps(right_gr, ctcf_gr))))
    list(ctcf = tads_df[both_ctcf, ], non_ctcf = tads_df[setdiff(seq_len(nrow(tads_df)), both_ctcf), ])
  }
  
  # ---- Hedges' g (treated minus control) -------------------------------------------------
  hedges_g <- function(ctrl, trt) {
    n1 <- length(ctrl)
    n2 <- length(trt)
    sp <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / (n1 + n2 - 2))
    if (!is.finite(sp) || sp == 0) stop("Pooled SD is zero or undefined")
    d <- (mean(trt) - mean(ctrl)) / sp
    j <- 1 - 3 / (4 * (n1 + n2) - 9)
    d * j
  }
  
  # ---- score every TAD subset in every condition at fixed TAD positions ------------------------
  # One value per TAD per condition = mean z-scored insulation at its two boundaries
  # (nearest insulation grid point to each boundary coordinate)
  grid <- function(x) round(x / resolution) * resolution
  tad_score <- function(tads_df, label) {
    v <- ins_key[[label]]
    left  <- v[paste(tads_df$chr1, grid(tads_df$x1))]
    right <- v[paste(tads_df$chr1, grid(tads_df$x2))]
    unname((left + right) / 2)
  }
  
  # Mean z-scored insulation vs distance from the TAD boundaries (same construction as
  # get_boundary_profiles() in insulation_analysis.R: both boundaries of every TAD, offsets
  # relative to the exact boundary, mean +/- SE per offset, rolling mean over roll_window offsets)
  boundary_profiles <- function(tads_df, label) {
    v <- ins_key[[label]]
    bnd_chr <- c(tads_df$chr1, tads_df$chr1)
    bnd_pos <- c(tads_df$x1, tads_df$x2)
    n_off  <- 2 * flank_kb * 1000 / resolution + 1
    first  <- ceiling((bnd_pos - flank_kb * 1000) / resolution)
    starts <- outer(first, 0:(n_off - 1), "+") * resolution
    inside <- as.vector(starts <= bnd_pos + flank_kb * 1000)
    z      <- unname(v[paste(rep(bnd_chr, n_off), as.vector(starts))])
    d <- data.frame(offset_kb = (as.vector(starts) - rep(bnd_pos, n_off)) / 1000, z = z)[inside & !is.na(z), ]
    if (nrow(d) == 0) stop("No insulation scores found around the boundaries (check chromosome naming)")
    d %>% dplyr::group_by(offset_kb) %>%
      dplyr::summarise(mean_score = mean(z), se_score = sd(z) / sqrt(dplyr::n()), n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(offset_kb) %>%
      dplyr::mutate(mean_score = zoo::rollmean(mean_score, k = roll_window, fill = NA, align = "center"),
                    se_score   = zoo::rollmean(se_score,   k = roll_window, fill = NA, align = "center"),
                    sample = label)
  }
  prof_store <- list()
  g_store <- list()
  
  set.seed(seed)
  rows <- list()
  chk_done <- FALSE
  for (pos in names(tad_paths)) {
    tads <- link_tads_anchor(load_arrowhead(tad_paths[[pos]]), gene_gr)
    sp <- filter_ctcf_tads(tads, ctcf_gr)
    cat("\n", pos, " TADs: ", nrow(tads), " (CTCF-anchored ", nrow(sp$ctcf), ", non-CTCF ", nrow(sp$non_ctcf), ")\n", sep = "")
    if (nrow(tads) == 0) stop("No TADs read for ", pos)
    for (ctcf_status in c("CTCF", "non-CTCF")) {
      sub_tads <- if (ctcf_status == "CTCF") sp$ctcf else sp$non_ctcf
      for (bin in c("All", BINS)) {
        d <- if (bin == "All") sub_tads else sub_tads[sub_tads$expr_bin == bin, ]
        if (nrow(d) < min_tads) { cat("  skipped ", pos, "/", ctcf_status, "/", bin, ": ", nrow(d), " TADs\n", sep = ""); next }
        ctrl_all <- tad_score(d, control_label)
        key <- paste(pos, ctcf_status, bin, sep = "|")
        prof_store[[key]] <- list(n = nrow(d), prof = lapply(setNames(all_labels, all_labels), function(s) boundary_profiles(d, s)))
        for (tn in names(treated_labels)) {
          trt_all <- tad_score(d, treated_labels[[tn]])
          ok <- !is.na(ctrl_all) & !is.na(trt_all)
          if (sum(ok) < min_tads) { cat("  skipped ", pos, "/", ctcf_status, "/", bin, "/", tn, ": ", sum(ok), " TADs with scores\n", sep = ""); next }
          if (sum(!ok) > 0.05 * nrow(d)) stop("More than 5% of TADs lack insulation scores for ", pos, "/", ctcf_status, "/", bin, " (", sum(!ok), " of ", nrow(d), ")")
          ctrl <- ctrl_all[ok]
          trt  <- trt_all[ok]
          g <- hedges_g(ctrl, trt)
          if (!chk_done) {   # self-check on live data: g must equal the pooled t-statistic conversion
            chk_t <- unname(t.test(trt, ctrl, var.equal = TRUE)$statistic)
            chk_d <- chk_t * sqrt(1 / length(trt) + 1 / length(ctrl))
            chk_j <- 1 - 3 / (4 * (length(trt) + length(ctrl)) - 9)
            stopifnot(isTRUE(all.equal(g, chk_d * chk_j)))
            cat("Self-check passed: g matches the pooled t-statistic conversion\n")
            chk_done <- TRUE
          }
          n <- length(ctrl)
          bg <- replicate(n_boot, { i <- sample.int(n, replace = TRUE); hedges_g(ctrl[i], trt[i]) })   # TADs resampled jointly (paired)
          ci <- quantile(bg, c(0.025, 0.975), names = FALSE)
          stopifnot(ci[1] <= g, g <= ci[2])
          rows[[length(rows) + 1]] <- data.frame(
            assay = "HiC_TAD_insulation", comparison = paste0(tn, "_vs_WT"),
            tad_positions = pos, ctcf_status = ctcf_status, expr_bin = bin,
            n_tads = n, mean_ctrl = mean(ctrl), mean_trt = mean(trt),
            hedges_g = g, ci_lo = ci[1], ci_hi = ci[2], stringsAsFactors = FALSE)
          g_store[[paste(key, tn, sep = "|")]] <- c(g = g, lo = ci[1], hi = ci[2])
        }
      }
    }
  }
  
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  cat("\nHigher insulation z-score = more contacts crossing the boundary (weaker insulation)\n\n")
  print(res, digits = 3)
  write.csv(res, out_file, row.names = FALSE)
  cat("Wrote:", out_file, "\n")
  
  # ---- profile plots: POLR2A-degraded (left) next to ActD (right), one figure per TAD subset ----
  line_cols <- c("Wildtype" = "#88CCEE", "POLR2A deg" = "#CC6677", "ActD 1hr" = "#4DAF4A")
  line_cols <- line_cols[all_labels]
  plot_names <- c("POLR2A deg" = "POLR2A-degraded", "ActD 1hr" = "ActD (1 hr)", "Wildtype" = "Wildtype")
  figs <- list()
  for (pos in names(tad_paths)) {
    for (ctcf_status in c("CTCF", "non-CTCF")) {
      keys <- paste(pos, ctcf_status, c("All", BINS), sep = "|")
      keys <- keys[keys %in% names(prof_store)]
      if (length(keys) == 0) next
      all_prof <- dplyr::bind_rows(lapply(keys, function(k) dplyr::bind_rows(prof_store[[k]]$prof)))
      ylim <- range(c(all_prof$mean_score - all_prof$se_score, all_prof$mean_score + all_prof$se_score), na.rm = TRUE)
      panels <- list()
      for (k in keys) {
        bin <- sub(".*\\|", "", k)
        for (tn in names(treated_labels)) {
          lab <- treated_labels[[tn]]
          gi  <- g_store[[paste(k, tn, sep = "|")]]
          pd  <- dplyr::bind_rows(prof_store[[k]]$prof[c(control_label, lab)])
          pd$sample <- factor(pd$sample, levels = all_labels)
          panels[[length(panels) + 1]] <- ggplot(pd, aes(x = offset_kb, y = mean_score, color = sample, fill = sample)) +
            geom_ribbon(aes(ymin = mean_score - se_score, ymax = mean_score + se_score), alpha = 0.15, color = NA, na.rm = TRUE) +
            geom_line(linewidth = 0.9, na.rm = TRUE) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
            scale_color_manual(values = line_cols, breaks = all_labels, labels = plot_names[all_labels], name = NULL, drop = FALSE) +
            scale_fill_manual(values = line_cols, breaks = all_labels, labels = plot_names[all_labels], name = NULL, drop = FALSE) +
            scale_x_continuous(name = "Distance from TAD boundary (kb)", breaks = seq(-flank_kb, flank_kb, by = 250)) +
            scale_y_continuous(name = "Mean insulation score (z-score)", limits = ylim) +
            labs(title = paste0(bin, ": ", plot_names[[lab]], " vs Wildtype"),
                 subtitle = if (is.null(gi)) paste0(prof_store[[k]]$n, " TADs") else
                   sprintf("%d TADs | g = %.2f [%.2f, %.2f]", prof_store[[k]]$n, gi[["g"]], gi[["lo"]], gi[["hi"]])) +
            theme_classic(base_size = 11) +
            theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
                  plot.subtitle = element_text(size = 8, color = "gray40", hjust = 0.5), legend.position = "bottom")
        }
      }
      fig <- patchwork::wrap_plots(panels, ncol = length(treated_labels)) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
      fig <- fig + patchwork::plot_annotation(title = paste0("Insulation around TAD boundaries | ", pos, " TAD positions | ", ctcf_status, " TADs | ds_all"))
      print(fig)
      # SVG export: same figure, one panel column per treatment, one panel row per expression bin, plus room for title and legend
      f <- file.path(plot_dir, paste0("HiC_TAD_insulation_profiles_", pos, "_TADs_", gsub("[^A-Za-z0-9._-]+", "-", ctcf_status), ".svg"))
      w <- panel_w_in * length(treated_labels)
      h <- panel_h_in * ceiling(length(panels) / length(treated_labels)) + 0.9
      svg(f, width = w, height = h)
      print(fig)
      dev.off()
      message(sprintf("Wrote %s  (%.1f x %.1f in)", f, w, h))
      figs[[paste(pos, ctcf_status, sep = "|")]] <- fig
    }
  }
  invisible(list(table = res, plots = figs))
}

hic_tad_insulation <- hic_tad_insulation_effect_sizes()