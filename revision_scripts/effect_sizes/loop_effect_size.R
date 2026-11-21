hic_loop_effect_sizes <- function(
    base_dir   = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp",
    out_file   = "/projects/b1042/BackmanLab/Lucas/102726_effect_size/HiC/HiC_loop_effect_sizes.csv",
    # repo functions (load_narrowpeak, classify_anchors_narrowpeak, GROUP_MAP), sourced into their own environment
    ep_funcs_file = "~/Transcription Paper/Hic_loop_analysis/final_scripts/EP_analysis_functions.R",
    # Anchor classification = classify_anchors_narrowpeak() as in hic_loop_lfc_analysis.R / final_APA_analysis.R:
    # untreated (DMSO) peaks for control + shared loops; 6 h IAA peaks for treated-only loops of the POLR2A-degraded
    # comparison; untreated peaks for ActD loops (no ActD ChIP/ATAC exists, as in final_APA_analysis.R)
    peaks_dir  = "/projects/b1042/BackmanLab/Lucas/loop_chromatin_regression/data",
    peak_files = list(wt = c(atac = "ATAC_DMSO.narrowPeak", h3k4me3 = "H3K4me3_DMSO.narrowPeak",
                             h3k27ac = "H3K27ac_DMSO.narrowPeak", h3k4me1 = "H3K4me1_DMSO.narrowPeak"),
                      ko = c(atac = "ATAC_IAA.narrowPeak", h3k4me3 = "H3K4me3_IAA.narrowPeak",
                             h3k27ac = "H3K27ac_IAA.narrowPeak", h3k4me1 = "H3K4me1_IAA.narrowPeak")),
    # TxDb is built from the Ensembl 115 GTF
    gtf_file   = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference/Homo_sapiens.GRCh38.115.chr.gtf",
    tss_region = c(-3000, 3000),  # as in the repo's classify_anchors_narrowpeak() calls
    # ENCODE HCT116 RSEM gene quantification (ENCFF435PHM, GENCODE V29 IDs)
    expr_file  = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference/ENCFF435PHM.tsv",
    deseq_file = "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv",
    # Downsampling groups (same switch as hic_loop_lfc_analysis.R): loop calls, .hic files, depth-matched conditions
    groups = list(
      ds_all = list(
        control = list(label = "Wildtype",
                       loops = "loops/ds_all/WT_HCT116_CTRL/mega/merged_loops.bedpe",
                       hic   = "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
        treated = list(POLR2A_deg_6hr = list(label = "POLR2A deg",
                                             loops = "loops/ds_all/Pol2_6hrs_Aux/mega/merged_loops.bedpe",
                                             hic   = "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic", peaks = "ko"),
                       ActD_1hr = list(label = "ActD 1hr",
                                       loops = "loops/ds_all/1hr_ActD/mega/merged_loops.bedpe",
                                       hic   = "mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic", peaks = "wt"))),
      ds_wt_pol2 = list(
        control = list(label = "Wildtype",
                       loops = "loops/ds_wt_pol2/WT_HCT116_CTRL/mega/merged_loops.bedpe",
                       hic   = "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
        treated = list(POLR2A_deg_6hr = list(label = "POLR2A deg",
                                             loops = "loops/ds_wt_pol2/Pol2_6hrs_Aux/mega/merged_loops.bedpe",
                                             hic   = "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic", peaks = "ko")))),
    straw_norm       = "SCALE",   # STRAW_NORM
    straw_resolution = 5000L,     # STRAW_RESOLUTION (resolution used for loop calling)
    min_obs          = 0,         # MIN_OBS_WT / MIN_OBS_POL2 for downsampled groups
    max_fdr          = 0.1,       # filter_loops default
    log2fc_threshold = 0.5,       # LOG2FC_THRESHOLD
    min_contacts     = 5,         # MIN_WT_CONTACTS / MIN_POL2_CONTACTS
    pseudocount      = 1,         # PSEUDOCOUNT
    shared_slop_bp   = 5000L,     # ANCHOR_SLOP_BP (shared-loop matching)
    class_slop_bp    = 2500L,     # classify_anchors_narrowpeak() default slop_bp (as called in the repo)
    gene_tol_bp      = 5000,      # GENE_TOL_BP
    min_loops        = 20,
    n_boot = 10000,
    seed   = 102726) {
  
  library(dplyr)
  library(tidyr)
  library(readr)
  
  # ---- input checks ------------------------------------------------------------
  for (pk in c("GenomicRanges", "IRanges", "strawr", "ChIPseeker", "GenomicFeatures", "GenomeInfoDb")) if (!requireNamespace(pk, quietly = TRUE)) stop("R package '", pk, "' is not installed")
  library(GenomicRanges)
  ep_funcs_file <- path.expand(ep_funcs_file)
  grp_paths <- unlist(lapply(groups, function(g) c(file.path(base_dir, g$control$loops), file.path(base_dir, g$control$hic),
                                                   unlist(lapply(g$treated, function(t) c(file.path(base_dir, t$loops), file.path(base_dir, t$hic)))))))
  peak_paths <- file.path(peaks_dir, unlist(peak_files))
  inputs <- unique(c(ep_funcs_file, peak_paths, gtf_file, expr_file, deseq_file, grp_paths))
  if (any(!file.exists(inputs))) stop("Input file(s) not found:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  
  # ---- repo functions (own environment so their constants do not overwrite ours) ----
  ep <- new.env()
  sys.source(ep_funcs_file, envir = ep)
  for (fn in c("load_narrowpeak", "classify_anchors_narrowpeak")) if (!exists(fn, envir = ep, inherits = FALSE)) stop(fn, " not found in ", ep_funcs_file)
  if (!exists("GROUP_MAP", envir = ep, inherits = FALSE)) stop("GROUP_MAP not found in ", ep_funcs_file)
  
  # ---- expression bins (WT TPM tertiles of expressed genes; same as hic_loop_lfc_analysis.R) ----
  BINS <- c("Low", "Medium", "High")
  expr_raw <- read_tsv(expr_file, show_col_types = FALSE)
  if (any(!c("gene_id", "TPM") %in% colnames(expr_raw))) stop("Columns 'gene_id' and 'TPM' required in ", basename(expr_file), "\nColumns present: ", paste(colnames(expr_raw), collapse = ", "))
  expr <- expr_raw %>% dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
    dplyr::group_by(gene_id) %>% dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")
  deseq_res <- read_csv(deseq_file, show_col_types = FALSE)
  for (cc in c("gene_id", "symbol", "gene_seq_start", "gene_seq_end", "seq_name")) {
    if (!cc %in% colnames(deseq_res)) stop("Column '", cc, "' not found in ", basename(deseq_file), "\nColumns present: ", paste(colnames(deseq_res), collapse = ", "))
  }
  deseq_res <- deseq_res %>% dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>% left_join(expr, by = "gene_id")
  cat("DESeq genes with a TPM value:", sum(!is.na(deseq_res$TPM)), "/", nrow(deseq_res), "\n")
  if (mean(!is.na(deseq_res$TPM)) < 0.8) stop("Fewer than 80% of genes in the DESeq table matched the TPM table; check gene_id formats")
  expressed  <- deseq_res %>% filter(!is.na(TPM), TPM > 0) %>% pull(TPM)
  tpm_medium <- quantile(expressed, 1/3, na.rm = TRUE)
  tpm_high   <- quantile(expressed, 2/3, na.rm = TRUE)
  cat("Expression bin thresholds (WT TPM tertiles): Low <=", round(tpm_medium, 2), "| Medium <=", round(tpm_high, 2), "| High >", round(tpm_high, 2), "\n")
  deseq_res <- deseq_res %>%
    dplyr::mutate(expr_bin = dplyr::case_when(is.na(TPM) | TPM == 0 ~ "Low", TPM <= tpm_medium ~ "Low", TPM <= tpm_high ~ "Medium", TRUE ~ "High"))
  gene_coords <- deseq_res %>% filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name)) %>% distinct(gene_id, .keep_all = TRUE)
  gene_gr <- GRanges(seqnames = as.character(gene_coords$seq_name), ranges = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end),
                     expr_bin = gene_coords$expr_bin, TPM = gene_coords$TPM)
  
  # ---- narrowPeaks + TxDb (as loaded in hic_loop_lfc_analysis.R) ---------------------------
  peaks <- lapply(peak_files, function(v) lapply(v, function(f) ep$load_narrowpeak(file.path(peaks_dir, f))))
  cat("narrowPeak counts:\n"); print(sapply(peaks, function(p) sapply(p, length)))
  cat("Building TxDb from ", basename(gtf_file), "\n", sep = "")
  np_txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
  tx_chr  <- any(grepl("^chr", GenomeInfoDb::seqlevels(np_txdb)))
  chr_fmt <- function(x) if (tx_chr) paste0("chr", sub("^chr", "", x)) else sub("^chr", "", x)   # ChIPseeker needs loop and TxDb chromosome styles to match
  
  # ---- loop functions (copied from hic_loop_lfc_analysis.R) ------------------------------
  load_hiccups_bedpe <- function(path) {
    all_lines   <- readLines(path)
    header_line <- all_lines[grepl("^#chr1", all_lines)][1]
    if (is.na(header_line)) stop("No '#chr1' header line in ", path)
    col_names   <- strsplit(sub("^#", "", header_line), "\t")[[1]]
    data_lines  <- all_lines[!grepl("^#", all_lines)]
    if (length(data_lines) == 0) stop("No data lines in ", path)
    df <- read.table(text = paste(data_lines, collapse = "\n"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
    if (ncol(df) == length(col_names)) colnames(df) <- col_names
    need <- c("observed", "expectedDonut", "fdrDonut", "fdrBL", "fdrH", "fdrV", "centroid1", "centroid2")
    if (any(!need %in% colnames(df))) stop("Columns missing in ", path, ": ", paste(setdiff(need, colnames(df)), collapse = ", "))
    data.frame(chr1 = as.character(df[[1]]), start1 = as.integer(df[[2]]), end1 = as.integer(df[[3]]),
               chr2 = as.character(df[[4]]), start2 = as.integer(df[[5]]), end2 = as.integer(df[[6]]),
               observed = as.numeric(df[["observed"]]), expectedDonut = as.numeric(df[["expectedDonut"]]),
               fdrDonut = as.numeric(df[["fdrDonut"]]), fdrBL = as.numeric(df[["fdrBL"]]),
               fdrH = as.numeric(df[["fdrH"]]), fdrV = as.numeric(df[["fdrV"]]),
               centroid1 = as.integer(df[["centroid1"]]), centroid2 = as.integer(df[["centroid2"]])) %>%
      filter(!is.na(start1), !is.na(start2))
  }
  filter_loops <- function(df) {
    df %>% filter(observed >= min_obs, fdrDonut <= max_fdr, fdrBL <= max_fdr, fdrH <= max_fdr, fdrV <= max_fdr)
  }
  find_shared_loops <- function(loops_a, loops_b, slop = shared_slop_bp) {
    make_gr <- function(df) list(a1 = GRanges(df$chr1, IRanges(df$start1, df$end1)), a2 = GRanges(df$chr2, IRanges(df$start2, df$end2)))
    gr_a <- make_gr(loops_a); gr_b <- make_gr(loops_b)
    expand <- function(gr) resize(gr, width(gr) + 2 * slop, fix = "center")
    hits_a1 <- findOverlaps(expand(gr_a$a1), expand(gr_b$a1))
    hits_a2 <- findOverlaps(expand(gr_a$a2), expand(gr_b$a2))
    shared <- dplyr::inner_join(as.data.frame(hits_a1) %>% dplyr::rename(a_idx = queryHits, b_idx = subjectHits),
                                as.data.frame(hits_a2) %>% dplyr::rename(a_idx = queryHits, b_idx = subjectHits), by = c("a_idx", "b_idx")) %>%
      dplyr::distinct(a_idx, b_idx)
    list(shared = shared, a_only = setdiff(seq_len(nrow(loops_a)), shared$a_idx), b_only = setdiff(seq_len(nrow(loops_b)), shared$b_idx))
  }
  # contacts at each loop's own centroid pair, snapped to the strawr resolution (as extract_loop_contacts in the repo);
  # chromosome names are matched to whatever naming the .hic file uses ("1" or "chr1")
  extract_loop_contacts <- function(loops_df, hic_file) {
    snap <- function(x) as.integer(floor(x / straw_resolution) * straw_resolution)
    results <- rep(NA_real_, nrow(loops_df))
    c1_pos <- snap(loops_df$centroid1); c2_pos <- snap(loops_df$centroid2)
    hic_names <- strawr::readHicChroms(hic_file)$name
    hic_name_of <- setNames(hic_names, sub("^chr", "", hic_names))
    chr1_clean <- sub("^chr", "", loops_df$chr1); chr2_clean <- sub("^chr", "", loops_df$chr2)
    if (any(!c(chr1_clean, chr2_clean) %in% names(hic_name_of))) stop("Loop chromosomes not in ", basename(hic_file), ": ",
                                                                      paste(setdiff(c(chr1_clean, chr2_clean), names(hic_name_of)), collapse = ", "))
    for (cp in unique(paste(chr1_clean, chr2_clean, sep = "|"))) {
      parts <- strsplit(cp, "\\|")[[1]]
      idx   <- which(chr1_clean == parts[1] & chr2_clean == parts[2])
      loc1  <- paste0(hic_name_of[[parts[1]]], ":", min(c1_pos[idx]), ":", max(c1_pos[idx]) + straw_resolution)
      loc2  <- paste0(hic_name_of[[parts[2]]], ":", min(c2_pos[idx]), ":", max(c2_pos[idx]) + straw_resolution)
      mat <- tryCatch(strawr::straw(straw_norm, hic_file, loc1, loc2, "BP", straw_resolution, matrix = "observed"),
                      error = function(e) stop("strawr failed for ", loc1, " / ", loc2, " in ", basename(hic_file), ": ", e$message))
      if (nrow(mat) == 0) next
      lookup <- setNames(mat$counts, paste(pmin(mat$x, mat$y), pmax(mat$x, mat$y), sep = "|"))
      vals   <- lookup[paste(pmin(c1_pos[idx], c2_pos[idx]), pmax(c1_pos[idx], c2_pos[idx]), sep = "|")]
      results[idx] <- unname(vals)
    }
    n_hit <- sum(!is.na(results) & results > 0)
    cat("  ", basename(dirname(dirname(dirname(hic_file)))), ": signal at ", n_hit, " / ", nrow(loops_df), " loop positions (", round(100 * n_hit / nrow(loops_df), 1), "%)\n", sep = "")
    if (n_hit < 0.5 * nrow(loops_df)) stop("Fewer than 50% of loop positions have contacts in ", hic_file, "; check resolution / norm / coordinates")
    results
  }
  link_loops_to_genes <- function(loops_df) {
    fix_chr <- function(gr) {
      has_chr_ref <- any(grepl("^chr", as.character(seqnames(gene_gr)))); has_chr_gr <- any(grepl("^chr", as.character(seqnames(gr))))
      if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
      else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
      suppressWarnings(gr)
    }
    a1 <- fix_chr(GRanges(loops_df$chr1, IRanges(loops_df$start1, loops_df$end1)))
    a2 <- fix_chr(GRanges(loops_df$chr2, IRanges(loops_df$start2, loops_df$end2)))
    gene_exp <- suppressWarnings(trim(resize(gene_gr, width(gene_gr) + 2 * gene_tol_bp, fix = "center")))
    h1 <- findOverlaps(a1, gene_exp); h2 <- findOverlaps(a2, gene_exp)
    pairs <- bind_rows(data.frame(loop_idx = queryHits(h1), gene_idx = subjectHits(h1)), data.frame(loop_idx = queryHits(h2), gene_idx = subjectHits(h2))) %>%
      distinct(loop_idx, gene_idx) %>% dplyr::mutate(expr_bin = gene_gr$expr_bin[gene_idx], TPM = gene_gr$TPM[gene_idx]) %>%
      dplyr::group_by(loop_idx) %>% dplyr::slice_max(TPM, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::select(loop_idx, expr_bin, TPM)
    loops_df %>% left_join(pairs, by = "loop_idx") %>% dplyr::mutate(expr_bin = replace_na(expr_bin, "Unannotated"))
  }
  
  # ---- effect-size functions ---------------------------------------------------------------
  hedges_g <- function(ctrl, trt) {
    n1 <- length(ctrl); n2 <- length(trt)
    sp <- sqrt(((n1 - 1) * var(ctrl) + (n2 - 1) * var(trt)) / (n1 + n2 - 2))
    if (!is.finite(sp) || sp == 0) stop("Pooled SD is zero or undefined")
    ((mean(trt) - mean(ctrl)) / sp) * (1 - 3 / (4 * (n1 + n2) - 9))
  }
  
  set.seed(seed)
  rows <- list()
  chk_done <- FALSE
  for (gname in names(groups)) {
    g <- groups[[gname]]
    cat("\n================ ", gname, " ================\n", sep = "")
    loops_ctrl <- filter_loops(load_hiccups_bedpe(file.path(base_dir, g$control$loops))) %>% dplyr::mutate(loop_idx = row_number())
    cat("Control loops:", nrow(loops_ctrl), "\n")
    for (tn in names(g$treated)) {
      tr <- g$treated[[tn]]
      cat("\n--- ", tr$label, " vs ", g$control$label, " (", gname, ") ---\n", sep = "")
      loops_trt <- filter_loops(load_hiccups_bedpe(file.path(base_dir, tr$loops))) %>% dplyr::mutate(loop_idx = row_number())
      cat("Treated loops:", nrow(loops_trt), "\n")
      ov <- find_shared_loops(loops_ctrl, loops_trt)
      cat("Control loops matching more than one treated loop (counted once here): ", nrow(ov$shared) - length(unique(ov$shared$a_idx)), "\n", sep = "")
      # union set: shared -> control coordinates; control-only -> control; treated-only -> treated (as hic_loop_lfc_analysis.R)
      union_loops <- bind_rows(loops_ctrl %>% dplyr::slice(sort(unique(ov$shared$a_idx))) %>% dplyr::mutate(source = "shared"),   # each control loop once, even if it matches several treated loops
                               loops_ctrl %>% dplyr::slice(ov$a_only) %>% dplyr::mutate(source = "ctrl_only"),
                               loops_trt  %>% dplyr::slice(ov$b_only) %>% dplyr::mutate(source = "trt_only")) %>%
        dplyr::mutate(loop_idx = row_number())
      cat("Union loops:", nrow(union_loops), " (shared ", sum(union_loops$source == "shared"), ", control-only ", sum(union_loops$source == "ctrl_only"),
          ", treated-only ", sum(union_loops$source == "trt_only"), ")\n", sep = "")
      stopifnot(nrow(union_loops) == nrow(loops_ctrl) + length(ov$b_only))
      
      ctrl_c <- replace_na(extract_loop_contacts(union_loops, file.path(base_dir, g$control$hic)), 0)
      trt_c  <- replace_na(extract_loop_contacts(union_loops, file.path(base_dir, tr$hic)), 0)
      union_loops <- union_loops %>%
        dplyr::mutate(ctrl_contacts = ctrl_c, trt_contacts = trt_c,
                      log2fc = log2((trt_contacts + pseudocount) / (ctrl_contacts + pseudocount)),
                      loop_class = dplyr::case_when(ctrl_contacts >= min_contacts & log2fc < -log2fc_threshold ~ "Lost",
                                                    trt_contacts >= min_contacts & log2fc > log2fc_threshold ~ "Gained",
                                                    TRUE ~ "Stable"))
      cat("Loop classes (Gained / Stable / Lost):\n"); print(table(factor(union_loops$loop_class, levels = c("Gained", "Stable", "Lost")), union_loops$source))
      
      # anchor types with the repo's classify_anchors_narrowpeak(): untreated peaks for control + shared loops,
      # the treated comparison's peak set (tr$peaks) for treated-only loops
      classify_with <- function(loops_df, set, label) {
        pk <- peaks[[set]]
        inp <- data.frame(loop_idx = loops_df$loop_idx, chr1 = chr_fmt(loops_df$chr1), start1 = loops_df$start1, end1 = loops_df$end1,
                          chr2 = chr_fmt(loops_df$chr2), start2 = loops_df$start2, end2 = loops_df$end2, stringsAsFactors = FALSE)
        out <- ep$classify_anchors_narrowpeak(loops_df = inp, atac_gr = pk$atac, h3k4me3_gr = pk$h3k4me3, h3k27ac_gr = pk$h3k27ac,
                                              txdb = np_txdb, h3k4me1_gr = pk$h3k4me1, tss_region = tss_region,
                                              ccre_label = label, slop_bp = class_slop_bp)
        stopifnot(nrow(out) == nrow(inp), identical(out$loop_idx, inp$loop_idx),
                  all(c("anchor1_ccre", "anchor2_ccre", "loop_ccre_type", "anchor1_zone", "anchor2_zone") %in% colnames(out)))
        # self-check: a chromosome-style mismatch would put every anchor in INTERGENIC
        if (!any(c(out$anchor1_zone, out$anchor2_zone) == "PROMOTER")) stop("No anchor in a promoter zone for ", label, "; check chromosome naming")
        out[, c("loop_idx", "anchor1_ccre", "anchor2_ccre", "loop_ccre_type")]
      }
      is_trt <- union_loops$source == "trt_only"
      cls <- bind_rows(classify_with(union_loops[!is_trt, ], "wt", "untreated_narrowpeak"),
                       classify_with(union_loops[is_trt, ], tr$peaks, paste0(tr$peaks, "_narrowpeak_treated_only")))
      classified <- union_loops %>% left_join(cls, by = "loop_idx")
      stopifnot(nrow(classified) == nrow(union_loops), !anyNA(classified$loop_ccre_type))
      cat("Anchors by type (both anchors pooled):\n"); print(table(c(classified$anchor1_ccre, classified$anchor2_ccre)))
      classified$group <- unname(ep$GROUP_MAP[classified$loop_ccre_type])
      classified$group <- ifelse(is.na(classified$group), "Other", classified$group)
      cat("Loops by type (P-P / E-P / E-E / Other):\n"); print(table(classified$group))
      classified$promoter_anchored <- classified$anchor1_ccre == "Promoter" | classified$anchor2_ccre == "Promoter"
      classified <- link_loops_to_genes(classified)
      cat("Loops linked to a gene (+/-", gene_tol_bp, " bp): ", sum(classified$expr_bin != "Unannotated"), " / ", nrow(classified), "\n", sep = "")
      
      sets <- list("All loops" = rep(TRUE, nrow(classified)), "Promoter-anchored" = classified$promoter_anchored,
                   "P-P" = classified$group == "P-P", "E-P" = classified$group == "E-P", "E-E" = classified$group == "E-E")
      for (set_name in names(sets)) {
        for (bin in c("All", BINS)) {
          d <- classified[sets[[set_name]] & (bin == "All" | classified$expr_bin == bin), ]
          if (nrow(d) < min_loops) { cat("  skipped ", set_name, "/", bin, ": ", nrow(d), " loops\n", sep = ""); next }
          ctrl <- log10(d$ctrl_contacts + 1)
          trt  <- log10(d$trt_contacts + 1)
          gval <- hedges_g(ctrl, trt)
          if (!chk_done) {   # self-check on live data: g must equal the pooled t-statistic conversion
            chk_t <- unname(t.test(trt, ctrl, var.equal = TRUE)$statistic)
            chk_d <- chk_t * sqrt(1 / length(trt) + 1 / length(ctrl))
            chk_j <- 1 - 3 / (4 * (length(trt) + length(ctrl)) - 9)
            stopifnot(isTRUE(all.equal(gval, chk_d * chk_j)))
            cat("Self-check passed: g matches the pooled t-statistic conversion\n")
            chk_done <- TRUE
          }
          n <- nrow(d)
          bg <- replicate(n_boot, { i <- sample.int(n, replace = TRUE); hedges_g(ctrl[i], trt[i]) })   # loops resampled jointly (paired)
          ci <- quantile(bg, c(0.025, 0.975), names = FALSE)
          stopifnot(ci[1] <= gval, gval <= ci[2])
          n_lost <- sum(d$loop_class == "Lost"); n_gained <- sum(d$loop_class == "Gained"); n_stable <- sum(d$loop_class == "Stable")
          stopifnot(n_lost + n_gained + n_stable == n)
          fl <- binom.test(n_lost, n)$conf.int
          if (n_lost + n_gained > 0) {   # Lost:Gained odds, exact CI from the conditional binomial (Lost among changed loops)
            bt <- binom.test(n_lost, n_lost + n_gained)$conf.int
            ratio <- if (n_gained > 0) n_lost / n_gained else NA_real_
            r_lo <- bt[1] / (1 - bt[1]); r_hi <- bt[2] / (1 - bt[2])
          } else { ratio <- r_lo <- r_hi <- NA_real_ }
          rows[[length(rows) + 1]] <- data.frame(
            assay = "HiC_loops", ds_group = gname, comparison = paste0(tn, "_vs_WT"), loop_set = set_name, expr_bin = bin,
            n_loops = n, mean_log10_ctrl = mean(ctrl), mean_log10_trt = mean(trt), hedges_g = gval, g_lo = ci[1], g_hi = ci[2],
            n_lost = n_lost, n_gained = n_gained, n_stable = n_stable,
            frac_lost = n_lost / n, frac_lost_lo = fl[1], frac_lost_hi = fl[2],
            lost_gained_ratio = ratio, ratio_lo = r_lo, ratio_hi = r_hi, stringsAsFactors = FALSE)
        }
      }
    }
  }
  
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  cat("\nhedges_g: treated minus control on log10(SCALE contacts + 1) at union loop positions; negative = weaker loops\n")
  cat("frac_lost = Lost loops / all loops in the set; lost_gained_ratio = Lost / Gained (exact binomial CI on Lost among changed loops)\n\n")
  print(res, digits = 3)
  write.csv(res, out_file, row.names = FALSE)
  cat("Wrote:", out_file, "\n")
  invisible(res)
}

hic_loop <- hic_loop_effect_sizes()