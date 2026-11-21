hic_scaling_effect_sizes <- function(
    # ds_all = Wildtype, POLR2A-degraded and ActD at equal depth
  base_dir   = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp",
  # ---- method 1: scaling_analysis.R (Fig. 4F-G as published) -------------------------------
  # contacts_<x>_genes.rds are the gene contacts cached by scaling_analysis.R (DS_GROUP = "ds_all"); the ActD file
  # is extracted with the same function if it does not exist yet
  rds_dir    = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/ei_scaling/rds/ds_all",   # RDS_DIR
  gtf_file   = "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf",   # GTF_FILE of scaling_analysis.R
  # ENCFF435PHM.tsv is the TPM table of scaling_controls_main.R and of the other Hi-C effect-size scripts
  expr_file  = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference/ENCFF435PHM.tsv",
  # ---- method 2: scaling_controls_main.R (equal-bp E/I and neighborhood-TPM pools) ----------------
  # Method 2 uses the GTF_FILE / EXPR_FILE defined in scaling_controls_main.R (Ensembl 115 chr GTF, ENCFF435PHM.tsv); method 1 uses the GTF above with the same TPM table.
  # rds_base_b caches the contact profiles (RDS_BASE/long/profiles_<wt|pol2|actd>_s1.rds). The cache is keyed on the regions only,
  # not on the .hic file, so it must not be a folder holding profiles from a different set of .hic files (e.g. mega_ds_wt_pol2).
  main_file  = "/home/lmc0633/Transcription Paper/101426_revision_scripts/micro-c/scripts/scaling_controls_main.R",
  rds_base_b = "/projects/b1042/BackmanLab/Lucas/scaling_controls_rds/hic_ds_all",
  conditions = list("Wildtype"   = list(key = "wt",   hic = "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic", rds = "contacts_wt_genes.rds"),
                    "POLR2A deg" = list(key = "pol2", hic = "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic",  rds = "contacts_pol2_genes.rds"),
                    "ActD 1hr"   = list(key = "actd", hic = "mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic",       rds = "contacts_actd_genes.rds")),
  out_file   = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/HiC/HiC_scaling_effect_sizes.csv",
  control_label  = "Wildtype",
  treated_labels = c(POLR2A_deg_6hr = "POLR2A deg", ActD_1hr = "ActD 1hr"),
  resolution    = 5000L,     # RESOLUTION
  straw_norm    = "SCALE",   # STRAW_NORM
  ei_bins       = 100L,      # EI_BINS
  min_gene_size = 100000L,   # MIN_GENE_SIZE
  fit_max_bp    = 10^6,      # slope fit window: 5 kb to 1 Mb
  min_bins      = 10,
  n_boot = 10000,
  seed   = 102726) {
  
  library(dplyr)
  library(tidyr)
  library(readr)
  library(data.table)
  
  # ---- input checks ------------------------------------------------------------
  if (!control_label %in% names(conditions)) stop("control_label not in conditions")
  if (any(!treated_labels %in% names(conditions))) stop("treated_labels not in conditions")
  rds_paths <- setNames(file.path(rds_dir, sapply(conditions, `[[`, "rds")), names(conditions))
  hic_paths <- setNames(file.path(base_dir, sapply(conditions, `[[`, "hic")), names(conditions))
  ei_rds <- file.path(rds_dir, "ei_table.rds")
  need_extract <- names(conditions)[!file.exists(rds_paths)]
  inputs <- c(expr_file, main_file, hic_paths, if (!file.exists(ei_rds)) gtf_file)
  if (any(!file.exists(inputs))) stop("Input file(s) not found (set main_file= to the location of scaling_controls_main.R):\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK"); slurm_ntasks <- Sys.getenv("SLURM_NTASKS")
  n_cores <- max(1L, (if (nchar(slurm_cpus) > 0) as.integer(slurm_cpus) else if (nchar(slurm_ntasks) > 0) as.integer(slurm_ntasks) else 2L) - 1L)   # get_slurm_cores, as in insulation_analysis.R
  message("Using ", n_cores, " cores")
  old_rng <- RNGkind(); RNGkind("L'Ecuyer-CMRG"); on.exit(do.call(RNGkind, as.list(old_rng)), add = TRUE)
  std_chroms <- paste0("chr", c(1:22, "X"))
  trt_names <- names(treated_labels)
  
  # ---- fast slope fits (shared by both methods) ----------------------------------------------
  # Per gene: summed counts at each 5 kb separation below 1 Mb, and floor(max observed separation / resolution).
  # Pooling a bin = summing gene rows, so a resample of genes needs no re-read of the contacts.
  # The constant col_sum of compute_slopes divides P(s) in every bin alike, which moves the intercept and leaves the slope unchanged.
  n_k <- ceiling(fit_max_bp / resolution)           # separations k = 0 .. n_k - 1 (x resolution bp)
  k_idx <- seq_len(n_k) - 1L
  lx <- log10(pmax(k_idx * resolution, 1))
  slope_xy <- function(x, y) { xm <- mean(x); sum((x - xm) * (y - mean(y))) / sum((x - xm)^2) }
  # method 1 = compute_slopes of scaling_analysis.R: possible.contacts assigned by row position
  fit_repo <- function(cnt, tp) {
    if (!is.finite(tp) || tp < 3) return(NA_real_)
    present <- cnt > 0
    pos  <- cumsum(present)
    keep <- present & k_idx >= 1L & pos <= floor(tp - 1 + 1e-10) + 1
    if (sum(keep) < 3) return(NA_real_)
    slope_xy(lx[keep], log10(cnt[keep] / ((tp - pos[keep] + 1) * resolution)))
  }
  # method 2 = pooled_dist / fit_pool of scaling_controls_main.R ("original" mode): possible = total_possible * res - index
  fit_ctrl <- function(cnt, tp) {
    if (!is.finite(tp)) return(NA_real_)
    possible <- (tp - k_idx) * resolution
    keep <- cnt > 0 & k_idx >= 1L & possible >= resolution
    if (sum(keep) < 3) return(NA_real_)
    slope_xy(lx[keep], log10(cnt[keep] / possible[keep]))
  }
  bin_slopes <- function(red, idx_list, fit) {
    vapply(idx_list, function(ix) fit(colSums(red$M[ix, , drop = FALSE]), mean(red$maxk[ix])), numeric(1))
  }
  
  # ---- Hedges' g (treated minus control; the bins / pools are the units, as in the Welch t-tests of Fig. 4F-G) ----
  hedges_g <- function(ctrl, trt, strict = FALSE) {
    n1 <- length(ctrl); n2 <- length(trt)
    sp <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / (n1 + n2 - 2))
    if (!is.finite(sp) || sp == 0) { if (strict) stop("Pooled SD is zero or undefined"); return(NA_real_) }
    ((mean(trt) - mean(ctrl)) / sp) * (1 - 3 / (4 * (n1 + n2) - 9))
  }
  
  # a gene set: rows = indices into the per-condition matrices, bins = rows per bin / pool, xp = percentile of the variable per row
  make_set <- function(rows, binid, nb, strata, family, n_red, gradient = NULL, xp = NULL) {
    b <- split(rows, factor(binid, levels = seq_len(nb)))
    xf <- NULL
    if (!is.null(xp)) { xf <- rep(NA_real_, n_red); xf[rows] <- xp }
    list(rows = rows, bins = b[sapply(b, length) > 0], nb = nb, strata = strata, family = family, gradient = gradient, xp = xf)
  }
  
  # ---- one method: observed effect sizes + gene bootstrap ------------------------------------------
  # delta_S and the gradient follow the sign of Fig. 4F-G (Wildtype minus treated); g is treated minus Wildtype like every other assay.
  # gradient = slope of the per-bin delta-S (Wildtype minus treated) on the median percentile (ecdf) of the variable in the bin,
  #   i.e. the change in delta-S from the 0th to the 100th percentile (single_var_deltas / single_var_effect of the scaling_controls_* scripts)
  chk_g <- FALSE
  run_method <- function(method, red, gene_sets, fit) {
    slope_matrix <- function(gs, resample) {
      idx <- if (resample) lapply(gs$bins, function(ix) ix[sample.int(length(ix), replace = TRUE)]) else gs$bins
      S <- sapply(names(red), function(cn) bin_slopes(red[[cn]], idx, fit))   # the same genes are drawn for every condition
      if (!is.null(gs$xp)) attr(S, "x") <- vapply(idx, function(ix) median(gs$xp[ix]), numeric(1))
      S
    }
    stats_for <- function(S, bins_sel, tl) {
      ctrl <- S[bins_sel, control_label]; trt <- S[bins_sel, treated_labels[[tl]]]
      ctrl <- ctrl[is.finite(ctrl)]; trt <- trt[is.finite(trt)]
      if (length(ctrl) < min_bins || length(trt) < min_bins) return(c(g = NA_real_, dS = NA_real_))
      c(g = hedges_g(ctrl, trt), dS = mean(ctrl) - mean(trt))
    }
    grad_fun <- function(S, tl) {
      d <- S[, control_label] - S[, treated_labels[[tl]]]
      x <- attr(S, "x")
      ok <- is.finite(d) & is.finite(x)
      if (sum(ok) < min_bins) return(NA_real_)
      slope_xy(x[ok], d[ok])
    }
    grad_stats <- function(S, gs) setNames(sapply(trt_names, function(tl) grad_fun(S, tl)), paste(gs$gradient, trt_names, "grad", sep = "|"))
    
    rows <- list()
    for (s in names(gene_sets)) {
      gs <- gene_sets[[s]]
      S <- slope_matrix(gs, FALSE)
      for (st in names(gs$strata)) {
        bs <- intersect(names(gs$bins), gs$strata[[st]])
        for (tl in trt_names) {
          ctrl <- S[bs, control_label]; trt <- S[bs, treated_labels[[tl]]]
          ctrl <- ctrl[is.finite(ctrl)]; trt <- trt[is.finite(trt)]
          if (length(ctrl) < min_bins || length(trt) < min_bins) { cat("  skipped ", st, " / ", tl, ": ", length(ctrl), " and ", length(trt), " bins with slopes\n", sep = ""); next }
          g <- hedges_g(ctrl, trt, strict = TRUE)
          if (!chk_g) {   # self-check on live data: g must equal the pooled t-statistic conversion
            chk_t <- unname(t.test(trt, ctrl, var.equal = TRUE)$statistic)
            chk_d <- chk_t * sqrt(1 / length(trt) + 1 / length(ctrl))
            chk_j <- 1 - 3 / (4 * (length(trt) + length(ctrl)) - 9)
            stopifnot(isTRUE(all.equal(g, chk_d * chk_j)))
            cat("Self-check passed: g matches the pooled t-statistic conversion\n")
            chk_g <<- TRUE
          }
          rows[[length(rows) + 1]] <- data.frame(
            assay = "HiC_contact_scaling", method = method, comparison = paste0(tl, "_vs_WT"), family = gs$family, stratum = st,
            n_genes = sum(sapply(gs$bins[bs], length)), n_bins_ctrl = length(ctrl), n_bins_trt = length(trt),
            mean_slope_ctrl = mean(ctrl), mean_slope_trt = mean(trt), delta_S_wt_minus_trt = mean(ctrl) - mean(trt),
            welch_p = t.test(trt, ctrl)$p.value, hedges_g = g, grad_wt_minus_trt = NA_real_, stringsAsFactors = FALSE)
        }
      }
      if (!is.null(gs$gradient)) for (tl in trt_names) {
        ctrl <- S[, control_label]; trt <- S[, treated_labels[[tl]]]
        rows[[length(rows) + 1]] <- data.frame(
          assay = "HiC_contact_scaling", method = method, comparison = paste0(tl, "_vs_WT"), family = gs$family, stratum = gs$gradient,
          n_genes = length(gs$rows), n_bins_ctrl = sum(is.finite(ctrl)), n_bins_trt = sum(is.finite(trt)),
          mean_slope_ctrl = mean(ctrl, na.rm = TRUE), mean_slope_trt = mean(trt, na.rm = TRUE),
          delta_S_wt_minus_trt = mean(ctrl, na.rm = TRUE) - mean(trt, na.rm = TRUE),
          welch_p = NA_real_, hedges_g = NA_real_, grad_wt_minus_trt = grad_fun(S, tl), stringsAsFactors = FALSE)
      }
    }
    res <- do.call(rbind, rows)
    
    # bootstrap (one seed per resample, seed + r, so the result does not depend on the number of cores)
    #  g and delta_S: the bins / pools are the units of g, so they are resampled (with replacement, within each stratum, the same bins for every
    #    condition). Resampling genes inside bins instead adds the sampling noise of every bin's slope a second time to the between-bin SD that
    #    sits in the denominator of g, which biases the interval low; the CHECK line below quantifies that share.
    #  gradient: genes resampled with replacement within each bin, the same genes in every condition, as bootstrap_single() of the scaling_controls_* scripts
    ci_list <- list()
    for (s in names(gene_sets)) {
      gs <- gene_sets[[s]]
      S0 <- slope_matrix(gs, FALSE)
      part_a <- NULL
      if (length(gs$strata) > 0) {
        message("Bootstrap (", method, " / ", s, "): ", n_boot, " resamples of bins within each stratum")
        part_a <- do.call(rbind, lapply(seq_len(n_boot), function(r) {
          set.seed(seed + r)
          unlist(lapply(names(gs$strata), function(st) {
            bs <- intersect(names(gs$bins), gs$strata[[st]])
            bs <- sample(bs, length(bs), replace = TRUE)
            unlist(lapply(trt_names, function(tl) { v <- stats_for(S0, bs, tl); setNames(v, paste(st, tl, names(v), sep = "|")) }))
          }))
        }))
      }
      part_b <- NULL
      if (!is.null(gs$gradient)) {
        message("Bootstrap (", method, " / ", s, "): ", n_boot, " resamples of genes within ", length(gs$bins), " bins on ", n_cores, " cores")
        chunks <- split(seq_len(n_boot), rep_len(seq_len(max(1L, n_cores)), n_boot))
        out <- parallel::mclapply(chunks, function(ch) {
          try(do.call(rbind, lapply(ch, function(r) { set.seed(seed + r); grad_stats(slope_matrix(gs, TRUE), gs) })), silent = TRUE)
        }, mc.cores = n_cores)
        bad <- sapply(out, function(o) inherits(o, "try-error") || is.null(o))
        if (any(bad)) stop("Bootstrap worker failed in set ", s, ": ", as.character(out[[which(bad)[1]]]))
        part_b <- do.call(rbind, out)
        # CHECK on live data: variance of one bin's slope from resampling its genes (200 resamples) as a share of the between-bin variance of S
        arr <- vapply(seq_len(200), function(r) { set.seed(seed + 50000L + r); as.vector(slope_matrix(gs, TRUE)) }, numeric(length(S0)))
        noise <- colMeans(matrix(apply(arr, 1, var, na.rm = TRUE), nrow(S0)), na.rm = TRUE)
        obs_v <- apply(S0, 2, var, na.rm = TRUE)
        cat(sprintf("  CHECK %s / %s: gene-resampling variance of a bin's S = %s of the between-bin variance (%s)\n", method, s,
                    paste(sprintf("%.0f%%", 100 * noise / obs_v), collapse = " / "), paste(colnames(S0), collapse = " / ")))
      }
      ci_list[[s]] <- do.call(cbind, Filter(Negate(is.null), list(part_a, part_b)))
    }
    boot <- do.call(cbind, ci_list)
    q <- function(key) quantile(boot[, key], c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    res$g_lo <- res$g_hi <- res$dS_lo <- res$dS_hi <- res$grad_lo <- res$grad_hi <- NA_real_
    for (i in seq_len(nrow(res))) {
      tl <- sub("_vs_WT$", "", res$comparison[i])
      if (is.na(res$hedges_g[i])) {   # percentile-gradient row
        cg <- q(paste(res$stratum[i], tl, "grad", sep = "|"))
        res$grad_lo[i] <- cg[1]; res$grad_hi[i] <- cg[2]
        if (!(cg[1] <= res$grad_wt_minus_trt[i] && res$grad_wt_minus_trt[i] <= cg[2])) warning("Observed gradient outside its bootstrap interval: ", res$stratum[i], " / ", tl)
        next
      }
      cg <- q(paste(res$stratum[i], tl, "g", sep = "|")); cd <- q(paste(res$stratum[i], tl, "dS", sep = "|"))
      res$g_lo[i] <- cg[1]; res$g_hi[i] <- cg[2]; res$dS_lo[i] <- cd[1]; res$dS_hi[i] <- cd[2]
      if (!(cg[1] <= res$hedges_g[i] && res$hedges_g[i] <= cg[2])) warning("Observed g outside its bootstrap interval: ", res$stratum[i], " / ", tl)
    }
    res
  }
  
  # =============================================================================================
  # METHOD 1: scaling_analysis.R, Fig. 4F-G as published (100 quantile bins of E/I and of expression)
  # =============================================================================================
  # ---- 1. E/I table (scaling_analysis.R section 1, unchanged) -------------------------
  if (file.exists(ei_rds)) {
    message("Loading cached EI table: ", ei_rds)
    ei_df <- readRDS(ei_rds)
  } else {
    for (pk in c("GenomicFeatures", "GenomicRanges", "GenomeInfoDb")) if (!requireNamespace(pk, quietly = TRUE)) stop("R package '", pk, "' is not installed")
    library(GenomicFeatures); library(GenomicRanges); library(GenomeInfoDb)
    message("Building TxDb from GTF...")
    txdb     <- makeTxDbFromGFF(gtf_file, format = "gtf")
    exons_gr <- exons(txdb) %>% GenomicRanges::reduce()
    genes_gr <- genes(txdb)
    seqlevelsStyle(genes_gr) <- "UCSC"
    seqlevelsStyle(exons_gr) <- "UCSC"
    genes_gr <- keepSeqlevels(genes_gr, intersect(seqlevels(genes_gr), std_chroms), pruning.mode = "coarse")
    exons_gr <- keepSeqlevels(exons_gr, intersect(seqlevels(exons_gr), std_chroms), pruning.mode = "coarse")
    genes_gr$gene_id <- sub("\\..*", "", names(genes_gr))
    hits     <- findOverlaps(genes_gr, exons_gr)
    overlaps <- GenomicRanges::pintersect(genes_gr[queryHits(hits)], exons_gr[subjectHits(hits)])
    exon_bp_per_gene <- tapply(width(overlaps), queryHits(hits), sum)
    exon_bp_vec <- rep(0L, length(genes_gr))
    exon_bp_vec[as.integer(names(exon_bp_per_gene))] <- as.integer(exon_bp_per_gene)
    ei_df <- data.frame(gene_id = genes_gr$gene_id, chr = as.character(seqnames(genes_gr)),
                        start = start(genes_gr), end = end(genes_gr), gene_bp = width(genes_gr),
                        exon_bp = exon_bp_vec, stringsAsFactors = FALSE) %>%
      dplyr::mutate(intron_bp = gene_bp - exon_bp, EI = exon_bp / pmax(intron_bp, 1), size = end - start) %>%
      dplyr::filter(exon_bp > 0, !is.infinite(EI), EI > 0, EI != 1)
    saveRDS(ei_df, ei_rds)
    message("Saved EI table: ", nrow(ei_df), " genes")
  }
  
  # ---- 2. expression bins (scaling_analysis.R section 2, unchanged) ---------------------
  expr <- read_tsv(expr_file, show_col_types = FALSE)
  if (any(!c("gene_id", "TPM") %in% colnames(expr))) stop("Columns 'gene_id' and 'TPM' required in ", basename(expr_file), "\nColumns present: ", paste(colnames(expr), collapse = ", "))
  expr <- expr %>% dplyr::select(gene_id, TPM) %>% dplyr::mutate(gene_id = sub("\\..*", "", gene_id))
  expressed  <- expr %>% dplyr::filter(TPM > 0) %>% dplyr::pull(TPM)
  tpm_medium <- quantile(expressed, 1/3, na.rm = TRUE)
  tpm_high   <- quantile(expressed, 2/3, na.rm = TRUE)
  expr <- expr %>%
    dplyr::mutate(expr_bin = dplyr::case_when(is.na(TPM) | TPM == 0 ~ "Low", TPM <= tpm_medium ~ "Low", TPM <= tpm_high ~ "Medium", TRUE ~ "High"),
                  expr_bin = factor(expr_bin, levels = c("Low", "Medium", "High")))
  ei_df <- ei_df %>%
    left_join(expr, by = "gene_id") %>%
    dplyr::mutate(TPM = replace_na(TPM, 0), expr_bin = replace_na(as.character(expr_bin), "Low"),
                  expr_bin = factor(expr_bin, levels = c("Low", "Medium", "High")))
  
  # ---- 5. E/I bins, single genes >= 100 kb (section 5) and 100 expression bins (section 15a), unchanged ----
  ei_binned <- ei_df %>%
    dplyr::filter(size >= min_gene_size) %>%
    dplyr::mutate(EI_bin = as.integer(ntile(EI, ei_bins)),
                  expr_bin_100 = as.integer(ntile(TPM, 100L)))
  message("EI binned single genes (size >= 100kb): ", nrow(ei_binned), " genes")
  
  # ---- 6. contact extraction (section 6, unchanged; only used for conditions without a cache) ----
  extract_contacts_for_genes <- function(hic_file, genes_df, resolution, norm, n_cores) {
    library(strawr)
    bin_cols <- intersect(names(genes_df), c("EI_bin", "expr_bin", "pol2_bin", "gene_id"))
    # deviation from scaling_analysis.R: the repo pastes "chr" in front of the names stored in the .hic, which gives "chrchr1" when the file
    # already names its chromosomes "chr1"; here the stored name is used for straw() and "chrN" for matching the genes
    hic_names <- strawr::readHicChroms(hic_file) %>% dplyr::filter(name != "ALL", !grepl("_", name)) %>% dplyr::pull(name)
    hmap <- setNames(hic_names, paste0("chr", sub("^chr", "", hic_names)))
    chroms <- intersect(unique(genes_df$chr), names(hmap))
    if (length(chroms) == 0) stop("No chromosome of the genes matches the .hic file.\n  .hic names: ", paste(head(hic_names, 6), collapse = ", "),
                                  "\n  gene chr values: ", paste(head(unique(genes_df$chr), 6), collapse = ", "))
    message("  Extracting contacts for ", length(chroms), " chromosomes...")
    results <- parallel::mclapply(chroms, function(chrom) {
      chrom_bare <- hmap[[chrom]]
      raw <- tryCatch(strawr::straw(norm, hic_file, chrom_bare, chrom_bare, "BP", resolution),
                      error = function(e) { message("  Warning ", chrom, ": ", e$message); NULL })
      if (is.null(raw) || nrow(raw) == 0) return(NULL)
      genes_chr <- genes_df %>% dplyr::filter(chr == chrom)
      if (nrow(genes_chr) == 0) return(NULL)
      raw_dt <- data.table::as.data.table(raw)
      raw_dt <- raw_dt[!is.nan(counts) & !is.na(counts)]
      if (nrow(raw_dt) == 0) return(NULL)
      gene_contacts <- lapply(seq_len(nrow(genes_chr)), function(i) {
        g   <- genes_chr[i, ]
        sub <- raw_dt[x >= g$start & x <= g$end]
        if (nrow(sub) == 0) return(NULL)
        sub[, chr := chrom]
        for (col in bin_cols) sub[, (col) := g[[col]]]
        sub
      })
      data.table::rbindlist(Filter(Negate(is.null), gene_contacts))
    }, mc.cores = n_cores, mc.preschedule = FALSE)
    data.table::rbindlist(Filter(Negate(is.null), results))
  }
  
  # ---- 7. slopes (compute_slopes, section 7, unchanged; used for the self-check) ----
  compute_slopes <- function(contacts_dt, filter_col, resolution) {
    contacts_dt <- data.table::copy(contacts_dt)
    contacts_dt[, index := abs(x - y)]
    contacts_dt <- contacts_dt[!is.nan(counts) & !is.na(counts)]
    bin_list <- split(contacts_dt, contacts_dt[[filter_col]])
    bin_list <- Filter(function(x) nrow(x) > 0, bin_list)
    dist_list <- lapply(names(bin_list), function(b) {
      df <- bin_list[[b]]
      gene_max_dist  <- df[, .(max_dist = max(index)), by = gene_id]
      n_gene         <- floor(gene_max_dist$max_dist / resolution)
      total_possible <- mean(n_gene, na.rm = TRUE)
      dist <- df[, .(counts = sum(counts)), by = index]
      data.table::setorder(dist, index)
      n <- nrow(dist)
      dist[, possible.contacts := seq(total_possible * resolution, resolution, by = -resolution)[seq_len(n)]]
      dist[, contact.prob      := counts / possible.contacts]
      dist[, bin_val           := as.integer(b)]
      dist
    })
    all_dist <- data.table::rbindlist(dist_list)
    col_sum  <- sum(all_dist$contact.prob, na.rm = TRUE)
    if (col_sum == 0) return(NULL)
    dplyr::bind_rows(lapply(dist_list, function(dist) {
      b <- dist$bin_val[1]
      dist[, contact.prob := contact.prob / col_sum]
      dist[index == 0, index := 1L]
      range <- dist[, .(distance = log10(index), prob = log10(contact.prob))]
      range <- range[is.finite(prob) & distance > 0 & distance < 6]
      if (nrow(range) < 3) return(NULL)
      model <- lm(prob ~ distance, data = range, na.action = "na.exclude")
      data.frame(bin_val = as.integer(b), slope = coef(model)[2])
    }))
  }
  
  reduce_contacts <- function(ct) {
    dt <- data.table::as.data.table(ct)[, .(gene_id, x, y, counts)]
    dt <- dt[!is.nan(counts) & !is.na(counts)]
    dt[, index := abs(x - y)]
    maxd <- dt[, .(max_dist = max(index)), by = gene_id]
    sub  <- dt[index < fit_max_bp]
    if (any(sub$index %% resolution != 0)) stop("Contact separations are not multiples of the resolution; check 'resolution'")
    sub[, k := as.integer(index / resolution)]
    agg <- sub[, .(counts = sum(counts)), by = .(gene_id, k)]
    ids <- maxd$gene_id
    M <- matrix(0, length(ids), n_k)
    M[cbind(match(agg$gene_id, ids), agg$k + 1L)] <- agg$counts
    list(ids = ids, M = M, maxk = floor(maxd$max_dist / resolution))
  }
  
  red1 <- list()
  chk_done <- FALSE
  for (cn in names(conditions)) {
    ct <- NULL
    if (file.exists(rds_paths[[cn]])) {
      message(cn, ": loading cached contacts: ", rds_paths[[cn]])
      ct <- readRDS(rds_paths[[cn]])
      if (!all(c("x", "y", "counts", "gene_id", "EI_bin", "expr_bin") %in% names(ct))) {
        message("  cache has no contact columns (", paste(names(ct), collapse = ", "), "; an earlier failed extraction); extracting again"); ct <- NULL }
    }
    if (is.null(ct)) {
      message(cn, ": extracting contacts (", hic_paths[[cn]], ")")
      ct <- extract_contacts_for_genes(hic_paths[[cn]], ei_binned, resolution, straw_norm, n_cores)
      if (nrow(ct) == 0) stop("Extraction returned no contacts for ", cn, " (", hic_paths[[cn]], "); nothing was saved")
      ct[, cond := cn]
      saveRDS(ct, rds_paths[[cn]])
      message("  Saved: ", rds_paths[[cn]])
    }
    ct <- data.table::as.data.table(ct)
    if (any(!c("x", "y", "counts", "gene_id", "EI_bin", "expr_bin") %in% names(ct))) stop("Columns missing in ", rds_paths[[cn]], "\nColumns present: ", paste(names(ct), collapse = ", "))
    # the strata stored with the contacts must be the strata rebuilt here
    gb <- ct[!duplicated(gene_id), .(gene_id, EI_bin, expr_bin = as.character(expr_bin))]
    mm <- match(gb$gene_id, ei_binned$gene_id)
    agree <- mean(!is.na(mm) & gb$EI_bin == ei_binned$EI_bin[mm] & gb$expr_bin == as.character(ei_binned$expr_bin[mm]))
    cat(sprintf("%s: %d genes with contacts; EI/expression bins agree with the rebuilt table for %.1f%%\n", cn, nrow(gb), 100 * agree))
    if (agree < 0.99) stop("Fewer than 99% of genes in ", rds_paths[[cn]], " have the same EI_bin / expr_bin as the table rebuilt from ", basename(expr_file), "; the cache is stale or came from other inputs")
    red1[[cn]] <- reduce_contacts(ct)
    if (!chk_done) {   # self-check on live data: the fast slopes must equal compute_slopes
      ib <- gb$EI_bin[match(red1[[cn]]$ids, gb$gene_id)]
      fast <- bin_slopes(red1[[cn]], split(seq_along(red1[[cn]]$ids), ib), fit_repo)
      rep_s <- suppressMessages(compute_slopes(ct, "EI_bin", resolution))
      d <- fast[as.character(rep_s$bin_val)] - rep_s$slope
      cat(sprintf("Self-check (%s, %d E/I bins): max |fast slope - compute_slopes| = %.2e\n", cn, nrow(rep_s), max(abs(d), na.rm = TRUE)))
      if (any(is.na(d)) || max(abs(d)) > 1e-8) stop("Fast slopes do not reproduce compute_slopes from scaling_analysis.R")
      chk_done <- TRUE
    }
    rm(ct); invisible(gc())
  }
  
  # common genes, same row order in every condition (the design is paired by gene)
  common <- Reduce(intersect, lapply(red1, `[[`, "ids"))
  cat("Genes with contacts in every condition:", length(common), "; per condition:",
      paste(names(red1), sapply(red1, function(r) length(r$ids)), sep = " = ", collapse = ", "), "\n")
  if (length(common) < 0.95 * max(sapply(red1, function(r) length(r$ids)))) stop("More than 5% of genes are missing from at least one condition")
  for (cn in names(red1)) {
    o <- match(common, red1[[cn]]$ids)
    red1[[cn]] <- list(ids = common, M = red1[[cn]]$M[o, , drop = FALSE], maxk = red1[[cn]]$maxk[o])
  }
  gt <- ei_binned[match(common, ei_binned$gene_id), c("gene_id", "EI", "TPM", "EI_bin", "expr_bin", "expr_bin_100")]
  gt$expr_bin <- as.character(gt$expr_bin)
  n1 <- length(common)
  
  # Fig. 4F: E/I tertiles of the 100 bins; Fig. 4G: Low / High expression (100 E/I bins each); Fig. 4G bottom: 100 expression bins
  ei_tert <- setNames(as.integer(ntile(seq_len(ei_bins), 3L)), seq_len(ei_bins))   # tertiles of the unique bins (1-34, 35-67, 68-100)
  all_bins <- as.character(seq_len(ei_bins))
  p_ei  <- ecdf(ei_binned$EI)(gt$EI)     # percentile of the variable among the binned genes (x_p of the scaling_controls_* scripts)
  p_tpm <- ecdf(ei_binned$TPM)(gt$TPM)
  sets1 <- list(
    EI       = make_set(seq_len(n1), gt$EI_bin, ei_bins,
                        list("All genes" = all_bins,
                             "Low E/I (bottom tertile of E/I bins)" = names(ei_tert)[ei_tert == 1L],
                             "High E/I (top tertile of E/I bins)"   = names(ei_tert)[ei_tert == 3L]),
                        "E/I ratio (Fig. 4F)", n1, "E/I percentile (100 bins)", p_ei),
    EXPR100  = make_set(seq_len(n1), gt$expr_bin_100, 100L, list(), "Expression (Fig. 4G, bottom)", n1, "Expression percentile (100 bins)", p_tpm),
    LowExpr  = make_set(which(gt$expr_bin == "Low"),  gt$EI_bin[gt$expr_bin == "Low"],  ei_bins,
                        list("Low expression (all E/I bins)" = all_bins), "Expression (Fig. 4G)", n1),
    HighExpr = make_set(which(gt$expr_bin == "High"), gt$EI_bin[gt$expr_bin == "High"], ei_bins,
                        list("High expression (all E/I bins)" = all_bins), "Expression (Fig. 4G)", n1))
  res1 <- run_method("scaling_analysis.R (Fig. 4F-G, quantile bins)", red1, sets1, fit_repo)
  rm(red1); invisible(gc())
  
  # =============================================================================================
  # METHOD 2: scaling_controls_main.R (equal-bp E/I pools and neighborhood-TPM pools of expressed genes >= 100 kb)
  # =============================================================================================
  message("\nLoading ", main_file)
  source(main_file)   # defines the functions and constants used below, in the global environment
  cond_b <- setNames(lapply(names(conditions), function(cn) list(display = cn, color = NA_character_, hic = hic_paths[[cn]])),
                     sapply(conditions, `[[`, "key"))
  assign("CONDITIONS", cond_b, envir = globalenv())   # read by get_profiles_all / extract_profiles
  assign("ASSAY", "Hi-C", envir = globalenv())
  cfg <- ANALYSIS_CFG$long
  if (cfg$res != resolution || cfg$fit_min != resolution || cfg$fit_max != fit_max_bp || cfg$strata[[1]][1] != min_gene_size)
    stop("ANALYSIS_CFG$long (res / fit window / minimum gene size) does not match this script's resolution / fit_max_bp / min_gene_size")
  if (any(!file.exists(c(GTF_FILE, EXPR_FILE)))) stop("Input file(s) of scaling_controls_main.R not found:\n", paste(c(GTF_FILE, EXPR_FILE)[!file.exists(c(GTF_FILE, EXPR_FILE))], collapse = "\n"))
  dir.create(rds_base_b, showWarnings = FALSE, recursive = TRUE)
  cs <- lapply(cond_b, function(cd) {
    ct <- strawr::readHicChroms(cd$hic)
    mp <- hic_chrom_map(cd$hic)
    data.frame(chr = names(mp), length = ct$length[match(mp, ct$name)], stringsAsFactors = FALSE)
  })
  chrom_sizes <- cs[[1]]
  for (k in cs) stopifnot(identical(k$length, chrom_sizes$length))
  for (cd in cond_b) check_hic_access(cd$hic, cfg, chrom_sizes)
  chrom_len <- setNames(chrom_sizes$length, chrom_sizes$chr)
  
  # annotation and neighborhoods exactly as in run_scaling_analysis() / run_stratum()
  all_genes <- load_gene_table(GTF_FILE, rds_base_b) %>% dplyr::left_join(load_expression(EXPR_FILE), by = "gene_id")
  all_genes <- add_neighborhood(all_genes, NBHD_BP, NONEXPR_TPM_MAX)
  genes_s <- all_genes %>%
    dplyr::filter(exon_bp > 0, !is.infinite(EI), EI > 0, EI != 1) %>%
    dplyr::filter(size >= cfg$strata[[1]][1], size < cfg$strata[[1]][2]) %>%
    dplyr::mutate(n_x_bins = floor(end / cfg$res) - ceiling(start / cfg$res) + 1) %>%
    dplyr::filter(n_x_bins >= 1)
  expressed_all <- genes_s %>% dplyr::filter(!is.na(TPM), TPM > 0)
  pl <- make_expressed_pools(expressed_all)
  cat("Method 2: expressed genes >= 100 kb:", nrow(expressed_all), "| E/I pools:", pl$n_ei, "| neighborhood-TPM pools:", pl$n_tpm,
      "| median bp per pool:", format(round(pl$pool_bp), big.mark = ","), "\n(compare with the 'E/I bins' line of the scaling_controls_hic.R log; run_stratum() can shrink the pools to fit the control supply, this script does not)\n")
  regions <- expressed_all %>%
    dplyr::transmute(region_id = gene_id, chr, start, end, size, n_x_bins, type = "expressed") %>% as.data.frame()
  rds_dir_b <- file.path(rds_base_b, "long")
  dir.create(rds_dir_b, showWarnings = FALSE, recursive = TRUE)
  ex <- get_profiles_all(regions, cfg, chrom_len, 1L, rds_dir_b)
  pass <- Reduce(intersect, lapply(ex, function(e) e$meta[n_x_present / n_x_bins >= MIN_COVERAGE]$region_id))   # coverage in every condition
  cat("Coverage >=", MIN_COVERAGE, "in every condition:", length(pass), "of", nrow(regions), "expressed genes\n")
  
  red2 <- lapply(ex, function(e) {
    pr  <- e$profile[index >= cfg$fit_min & index < cfg$fit_max]
    if (any(pr$index %% cfg$res != 0)) stop("Profile separations are not multiples of the resolution")
    ids <- e$meta$region_id
    ci  <- match(pr$region_id, ids)
    if (anyNA(ci)) stop("Profile rows without metadata")
    M <- matrix(0, length(ids), n_k)
    M[cbind(ci, as.integer(pr$index / cfg$res) + 1L)] <- pr$counts
    list(ids = ids, M = M, maxk = floor(e$meta$max_index / cfg$res))
  })
  common2 <- Reduce(intersect, c(lapply(red2, `[[`, "ids"), list(pass)))
  red2 <- lapply(red2, function(r) { o <- match(common2, r$ids); list(ids = common2, M = r$M[o, , drop = FALSE], maxk = r$maxk[o]) })
  ei_set  <- pl$ei_set  %>% dplyr::filter(gene_id %in% common2)   # pools are assigned before the coverage filter, as in run_stratum()
  tpm_set <- pl$tpm_set %>% dplyr::filter(gene_id %in% common2)
  
  # self-check on live data: the fast pooled slopes must equal compute_pool_slopes of scaling_controls_main.R
  pool_map <- dplyr::bind_rows(ei_set  %>% dplyr::transmute(region_id = gene_id, kind = "EI",  pool),
                               tpm_set %>% dplyr::transmute(region_id = gene_id, kind = "TPM", pool))
  ord_ei  <- split(match(ei_set$gene_id,  common2), ei_set$pool)
  ord_tpm <- split(match(tpm_set$gene_id, common2), tpm_set$pool)
  for (n in names(ex)) {
    ref  <- compute_pool_slopes(ex[[n]]$profile[region_id %in% common2], ex[[n]]$meta[region_id %in% common2], pool_map, cfg)
    mine <- c(setNames(bin_slopes(red2[[n]], ord_ei,  fit_ctrl), paste0("EI.",  names(ord_ei))),
              setNames(bin_slopes(red2[[n]], ord_tpm, fit_ctrl), paste0("TPM.", names(ord_tpm))))
    d <- mine[paste0(ref$kind, ".", ref$pool)] - ref$slope
    cat(sprintf("Self-check (%s, %d pools): max |fast slope - compute_pool_slopes| = %.2e\n", cond_b[[n]]$display, nrow(ref), max(abs(d), na.rm = TRUE)))
    if (any(is.na(d)) || max(abs(d)) > 1e-8) stop("Fast pooled slopes do not reproduce compute_pool_slopes from scaling_controls_main.R")
  }
  
  names(red2) <- sapply(cond_b, `[[`, "display")[names(red2)]   # the profiles are keyed wt / pol2 / actd; run_method uses the display names
  stopifnot(setequal(names(red2), names(conditions)))
  tert_of <- function(n) setNames(as.integer(dplyr::ntile(seq_len(n), 3)), seq_len(n))   # pool_info in run_stratum()
  te <- tert_of(pl$n_ei); tt <- tert_of(pl$n_tpm)
  n2 <- length(common2)
  sets2 <- list(
    EIpool  = make_set(match(ei_set$gene_id, common2), ei_set$pool, pl$n_ei,
                       list("All E/I pools" = as.character(seq_len(pl$n_ei)),
                            "Low E/I (bottom tertile of pools)" = names(te)[te == 1L],
                            "High E/I (top tertile of pools)"   = names(te)[te == 3L]),
                       "E/I ratio pools", n2, "E/I percentile (pools)", ecdf(ei_set$EI)(ei_set$EI)),
    TPMpool = make_set(match(tpm_set$gene_id, common2), tpm_set$pool, pl$n_tpm,
                       list("All neighborhood-TPM pools" = as.character(seq_len(pl$n_tpm)),
                            "Low neighborhood TPM (bottom tertile of pools)" = names(tt)[tt == 1L],
                            "High neighborhood TPM (top tertile of pools)"   = names(tt)[tt == 3L]),
                       "Neighborhood-TPM pools", n2, "Neighborhood-TPM percentile (pools)", ecdf(tpm_set$TPM_nbhd)(tpm_set$TPM_nbhd)))
  res2 <- run_method("scaling_controls_main.R (equal-bp E/I and neighborhood-TPM pools)", red2, sets2, fit_ctrl)
  
  # ---- output -------------------------------------------------------------------------------
  res <- rbind(res1, res2)
  res <- res[, c("assay", "method", "comparison", "family", "stratum", "n_genes", "n_bins_ctrl", "n_bins_trt", "mean_slope_ctrl", "mean_slope_trt",
                 "delta_S_wt_minus_trt", "dS_lo", "dS_hi", "welch_p", "hedges_g", "g_lo", "g_hi", "grad_wt_minus_trt", "grad_lo", "grad_hi")]
  rownames(res) <- NULL
  cat("\nhedges_g = treated minus Wildtype (positive = S less negative, i.e. flatter P(s)); units are bins (method 1) or pools (method 2)",
      "\ndelta_S_wt_minus_trt = Wildtype minus treated mean slope, the sign of Fig. 4F-G (negative = S rises after treatment)",
      "\ngrad_wt_minus_trt = change in delta_S (Wildtype minus treated) from the 0th to the 100th percentile of the variable (percentile rows only)\n\n")
  print(res, digits = 3)
  write.csv(res, out_file, row.names = FALSE)
  cat("Wrote:", out_file, "\n")
  invisible(res)
}

hic_scaling <- hic_scaling_effect_sizes()