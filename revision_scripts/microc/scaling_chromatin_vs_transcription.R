##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Reviewer comment (indirect effects): do the structural changes (delta-S)
# concentrate at loci with high WT Pol II occupancy and strong nascent loss,
# or do they follow chromatin state?
#
# Two analyses, each with its own heatmap, effect-size plot and CSV
# (S = slope of log10 P(s) vs log10 distance; delta-S = Wildtype - POLR2A deg;
# fits exactly as in scaling_controls_main.R / scaling_analysis.R):
#
#  1. "change":  is delta-S more likely to follow the CHANGE in transcription or the
#     CHANGE in a histone mark? Variables: nascent transcription change (EU-seq,
#     log2 WT/POLR2A deg) and, for H3K27ac, H3K4me3, H3K9me3, log2 DMSO/IAA of the
#     gene-body signal. Heatmap: one row per variable, columns = tertiles of that
#     variable alone (T1 = least decrease / most increase, T3 = largest decrease).
#
#  2. "state":   is delta-S explained by PRE-EXISTING chromatin state plus transcription
#     change? Variables: WT gene-body signal of Pol II and of each mark (DMSO; no change) and
#     nascent transcription change. Heatmap: for each variable, rows = tertile of nascent
#     transcription change (T3 = largest decrease), columns = tertile of the WT level.
#     Pol II (ENCODE POLR2A ChIP-seq, untreated) is included here only: it has no degraded-cell track.
#
# Every bigWig signal is the sum over the ENTIRE gene body divided by gene length (bases without a value
# count as 0), then divided by that bigWig's own genome-wide mean, so all tracks share one scale; the
# genome-wide mean of each bigWig and the median normalized gene signal are printed as checks.
#
# In both, the effect of every variable is estimated on its own: expressed genes
# are binned by that variable, delta-S is computed per bin, and its slope on the
# variable's percentile gives the change in delta-S from the 0th to the 100th
# percentile (gene bootstrap within bins). Plots per analysis: a heatmap with one
# row per variable, the crossed heatmap (transcription tertile x variable tertile)
# with the same data as a line plot below it, the single-variable effect plot, and
# a gene-length check (plot, table and CSV).
#
# Hi-C and Micro-C are run separately: set ASSAY_RUN below. The extracted
# contacts are the same regions as scaling_controls_*.R and reuse their cache
# (RDS_BASE/<analysis>/profiles_*.rds) when it exists.
##################################################################
.find_script_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(dirname(normalizePath(f)))
  }
  a <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (length(a) > 0) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error=function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
.main_file <- file.path(.find_script_dir(), "scaling_controls_main.R")
if (!file.exists(.main_file))
  stop("scaling_controls_main.R not found next to this script (looked in ",
       dirname(.main_file), ")")
source(.main_file)
message("Loaded shared functions from ", .main_file)

# ══════════════════════════════════════════════════════════════════════════════
# Settings
# ══════════════════════════════════════════════════════════════════════════════
ASSAY_RUN <- "Micro-C"            # "Micro-C" or "Hi-C"

# ── nascent transcription (EU-seq): per-gene RPKM, columns gene_id, POLR2A_1..3, WT_1..3
#    hg38 chrM-scaled EU-seq outputs; tertiles are rank based, so the global eu_sf factor does not change them
EUSEQ_FILE   <- "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/euseq/euseq_pipeline/counts/euseq_rpkm_analysis.txt"
EUSEQ_PSEUDO <- 0.1               # RPKM pseudocount

# ── chromatin variables. Analysis "change": log2(DMSO / IAA) of the gene-body mean signal
#    (positive = decrease in POLR2A-degraded cells, same direction as nascent transcription change);
#    analysis "state": the DMSO (Wildtype) gene-body mean signal itself.
#    files = DMSO (Wildtype) bigWig, files_deg = IAA (POLR2A-degraded) bigWig, from the loop-loss regression data dir.
#    region: "body" = mean signal over the gene body (all variables, the question is about gene geometry);
#    "tss" = mean within +/- TSS_FLANK_BP of the TSS is available but not used
TSS_FLANK_BP <- 2000L
MARK_PSEUDO_FRAC <- 0.01   # mark pseudocount = this fraction of the median positive gene-body signal (both conditions pooled)
DATA_DIR <- "/projects/b1042/BackmanLab/Lucas/loop_chromatin_regression/data"
ALL_VARIABLES <- list(
  POL2    = list(region="body", label="Pol II",   # Wildtype level only (analysis 2): ENCSR000EUU, untreated HCT116
                 files=file.path(DATA_DIR, "POLR2A_ENCODE.bw")),
  H3K27ac = list(region="body", label="H3K27ac",
                 files=file.path(DATA_DIR, "H3K27ac_DMSO.bw"), files_deg=file.path(DATA_DIR, "H3K27ac_IAA.bw")),
  H3K4me3 = list(region="body", label="H3K4me3",
                 files=file.path(DATA_DIR, "H3K4me3_DMSO.bw"), files_deg=file.path(DATA_DIR, "H3K4me3_IAA.bw")),
  H3K9me3 = list(region="body", label="H3K9me3",
                 files=file.path(DATA_DIR, "H3K9me3_DMSO.bw"), files_deg=file.path(DATA_DIR, "H3K9me3_IAA.bw"))
)
VARIABLES <- ALL_VARIABLES   # reduced per analysis in run_chromatin_analysis(): "change" needs a degraded-cell bigWig

COL_LEVELS <- c("L1", "L2", "L3", "NE")
COL_LABELS <- c(L1="Loss T1\n(least)", L2="Loss T2", L3="Loss T3\n(most)", NE="Non-expressed\ncontrols")
ROW_COLORS <- c("#9ECAE1", "#4292C6", "#08306B")

if (ASSAY_RUN == "Micro-C") {
  ASSAY <- "Micro-C"
  MICROC_DIR <- "/projects/b1042/BackmanLab/Lucas/microc_encode"
  CONDITIONS <- list(
    wt   = list(display="Wildtype",   color="#88CCEE", hic=file.path(MICROC_DIR, "ENCFF534TYD.hic")),
    pol2 = list(display="POLR2A deg", color="#CC6677", hic=file.path(MICROC_DIR, "ENCFF302CUI.hic")))
  ANALYSES <- c("long", "short")
  RDS_BASE <- "/projects/b1042/BackmanLab/Lucas/scaling_controls_rds/microc"
} else if (ASSAY_RUN == "Hi-C") {
  ASSAY <- "Hi-C"
  HIC_DIR <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/mega_ds_wt_pol2"
  CONDITIONS <- list(
    wt   = list(display="Wildtype",   color="#88CCEE",
                hic=file.path(HIC_DIR, "WT_HCT116_CTRL/mega/aligned/inter_30.hic")),
    pol2 = list(display="POLR2A deg", color="#CC6677",
                hic=file.path(HIC_DIR, "Pol2_6hrs_Aux/mega/aligned/inter_30.hic")))
  ANALYSES <- "long"
  RDS_BASE <- "/projects/b1042/BackmanLab/Lucas/scaling_controls_rds/hic"
} else stop("ASSAY_RUN must be 'Micro-C' or 'Hi-C'")
COMPARTMENT <- "all"

# ══════════════════════════════════════════════════════════════════════════════
# 1. Gene-level variables
# ══════════════════════════════════════════════════════════════════════════════
# strand of every gene (the cached gene table has none)
load_gene_strand <- function(gtf_file, rds_dir) {
  f <- file.path(rds_dir, "gene_strand.rds")
  if (file.exists(f)) return(readRDS(f))
  message("Reading gene strands from the GTF...")
  g <- readr::read_tsv(gtf_file, comment="#", col_names=FALSE, show_col_types=FALSE,
                       col_types=readr::cols(.default="c"))
  g <- g[g$X3 == "gene", ]
  out <- data.frame(gene_id=sub("\\..*", "", sub('.*gene_id "([^"]+)".*', "\\1", g$X9)),
                    strand=g$X7, stringsAsFactors=FALSE)
  out <- out[!duplicated(out$gene_id) & out$strand %in% c("+", "-"), ]
  saveRDS(out, f)
  out
}

# bigWig signal over a set of regions: the SUM of the signal over every base of the region divided by the
# region length (bases without a bigWig value count as 0), i.e. the mean density over the entire region.
# Also returns the genome-wide mean of the same bigWig (sum over all bases of the chromosomes in `chr` / their
# length), used to put every bigWig on one scale (fold over its own genome-wide mean).
bw_gene_signal <- function(bw_file, chr, start, end) {
  rl <- rtracklayer::import(bw_file, as="RleList")
  strip <- !any(grepl("^chr", names(rl)))
  ch  <- if (strip) sub("^chr", "", chr) else chr
  chs <- unique(ch[ch %in% names(rl)])
  gm  <- sum(vapply(chs, function(k) as.numeric(sum(rl[[k]])), numeric(1))) /
    sum(vapply(chs, function(k) as.numeric(length(rl[[k]])), numeric(1)))
  out <- rep(NA_real_, length(chr))
  for (k in chs) {
    i  <- which(ch == k)
    st <- pmax(1L, as.integer(start[i]))
    en <- pmin(as.integer(end[i]), length(rl[[k]]))
    ok <- !is.na(en) & en >= st
    if (!any(ok)) next
    v <- IRanges::Views(rl[[k]], st[ok], en[ok])
    out[i[ok]] <- IRanges::viewSums(v) / (en[ok] - st[ok] + 1)
  }
  list(mean=out, genome_mean=gm)
}

# signal of every mark for every gene: mode "change" = log2((DMSO + pseudo) / (IAA + pseudo)) of the region
# mean signal; mode "state" = the DMSO region mean. Cached by file paths, region type, pseudocount and mode.
# files[[v]] = c(wt=DMSO bigWig, deg=IAA bigWig)
load_gene_signals <- function(all_genes, strand_df, files, rds_dir, mode) {
  sig <- lapply(names(files), function(v) list(file=unname(files[[v]]), region=VARIABLES[[v]]$region,
                                               flank=TSS_FLANK_BP, pseudo_frac=MARK_PSEUDO_FRAC, mode=mode,
                                               norm="entire-gene sum / length, divided by genome-wide mean"))
  f <- file.path(rds_dir, paste0("gene_signals_", mode, ".rds"))
  if (file.exists(f)) {
    obj <- readRDS(f)
    if (identical(obj$sig, sig)) { message("Loading cached gene signals..."); return(obj$tab) }
  }
  g <- all_genes[, c("gene_id", "chr", "start", "end")]
  g$strand <- strand_df$strand[match(g$gene_id, strand_df$gene_id)]
  n_ns <- sum(is.na(g$strand))
  if (n_ns > 0) message("  genes without strand (use gene body for TSS variables, strand +): ", n_ns)
  g$strand[is.na(g$strand)] <- "+"
  tss <- ifelse(g$strand == "+", g$start, g$end)
  tab <- data.frame(gene_id=g$gene_id, stringsAsFactors=FALSE)
  # every bigWig is read once and cached per file (so the second analysis, and any re-run, reuses it);
  # different files are read in parallel
  jobs <- unique(do.call(rbind, lapply(names(files), function(v)
    data.frame(v=v, bw=unname(files[[v]]), stringsAsFactors=FALSE))))
  read_one <- function(i) {
    v <- jobs$v[i]; bw <- jobs$bw[i]
    cf <- file.path(rds_dir, paste0("bw_signal_", gsub("[^A-Za-z0-9]+", "_", basename(bw)), "_",
                                    VARIABLES[[v]]$region, ".rds"))
    key <- list(bw, as.numeric(file.info(bw)$mtime), VARIABLES[[v]]$region, TSS_FLANK_BP,
                nrow(g), sum(as.numeric(g$start)), sum(as.numeric(g$end)))
    if (file.exists(cf)) { obj <- readRDS(cf); if (identical(obj$key, key)) return(obj$x) }
    x <- if (VARIABLES[[v]]$region == "body") bw_gene_signal(bw, g$chr, g$start, g$end)
    else bw_gene_signal(bw, g$chr, tss - TSS_FLANK_BP, tss + TSS_FLANK_BP)
    saveRDS(list(key=key, x=x), cf)
    x
  }
  t_bw <- Sys.time()
  message("  reading ", nrow(jobs), " bigWig(s) (cached per file) on ", min(N_CORES, nrow(jobs)), " cores")
  raw <- parallel::mclapply(seq_len(nrow(jobs)), function(i) try(read_one(i), silent=TRUE),
                            mc.cores=min(N_CORES, nrow(jobs)))
  bad <- vapply(raw, function(x) inherits(x, "try-error"), logical(1))
  if (any(bad)) stop("Could not read bigWig ", jobs$bw[which(bad)[1]], ": ", as.character(raw[[which(bad)[1]]]))
  names(raw) <- jobs$bw
  message("    bigWig signals ready in ", round(as.numeric(difftime(Sys.time(), t_bw, units="mins")), 1), " min")
  # signal over the region divided by the bigWig's genome-wide mean, so every bigWig is on the same scale
  region_mean <- function(bw, v) {
    x <- raw[[bw]]
    message("    ", basename(bw), ": genome-wide mean ", signif(x$genome_mean, 4),
            "  | median gene signal / genome mean = ", signif(median(x$mean / x$genome_mean, na.rm=TRUE), 3))
    x$mean / x$genome_mean
  }
  for (v in names(files)) {
    wt <- region_mean(files[[v]]["wt"], v)
    if (mode == "state") {
      message("  signal ", v, " (", VARIABLES[[v]]$region, "): ", basename(files[[v]]["wt"]), " (Wildtype level)")
      tab[[v]] <- wt
    } else {
      message("  signal ", v, " (", VARIABLES[[v]]$region, "): ", basename(files[[v]]["wt"]), " vs ", basename(files[[v]]["deg"]))
      deg <- region_mean(files[[v]]["deg"], v)
      pos <- c(wt, deg); pos <- pos[!is.na(pos) & pos > 0]
      pseudo <- MARK_PSEUDO_FRAC * median(pos)
      message("    pseudocount = ", MARK_PSEUDO_FRAC, " x median positive normalized signal (", signif(median(pos), 3),
              ") = ", signif(pseudo, 3))
      tab[[v]] <- log2((wt + pseudo) / (deg + pseudo))
    }
  }
  saveRDS(list(sig=sig, tab=tab), f)
  tab
}

# nascent transcription loss = log2(WT / POLR2A deg) of per-gene EU-seq RPKM (positive = loss)
load_euseq_loss <- function(file) {
  x <- data.table::fread(file)
  need <- c("gene_id", paste0("POLR2A_", 1:3), paste0("WT_", 1:3))
  if (!all(need %in% names(x)))
    stop("EU-seq table lacks columns: ", paste(setdiff(need, names(x)), collapse=", "),
         "  (found: ", paste(names(x), collapse=", "), ")")
  deg <- rowMeans(as.matrix(x[, paste0("POLR2A_", 1:3), with=FALSE]))
  wt  <- rowMeans(as.matrix(x[, paste0("WT_", 1:3), with=FALSE]))
  data.frame(gene_id=sub("\\..*", "", x$gene_id),
             loss=log2((wt + EUSEQ_PSEUDO) / (deg + EUSEQ_PSEUDO)),
             stringsAsFactors=FALSE)[!duplicated(sub("\\..*", "", x$gene_id)), ]
}

# right-closed tertiles from fixed cut points (so controls share the expressed cut points)
tert_cuts <- function(x) {
  q <- unname(quantile(x, c(1/3, 2/3), na.rm=TRUE))
  if (q[1] >= q[2]) stop("Tertile cut points are not distinct (", signif(q[1], 3), ", ", signif(q[2], 3),
                         "): too many tied values in this variable")
  q
}
tert_of <- function(x, q) as.integer(cut(x, c(-Inf, q, Inf), labels=FALSE))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Regions: identical to run_stratum() in scaling_controls_main.R so the
#    cached contacts are reused (same seeds, same functions)
# ══════════════════════════════════════════════════════════════════════════════
build_shared_regions <- function(genes_s, all_genes, chrom_sizes, res) {
  genes_s <- genes_s %>% dplyr::mutate(n_x_bins = floor(end / res) - ceiling(start / res) + 1)
  genes_s <- genes_s %>% dplyr::filter(n_x_bins >= 1)
  expressed_all <- genes_s %>% dplyr::filter(!is.na(TPM), TPM > 0)
  zero_tpm  <- genes_s %>% dplyr::filter(!is.na(TPM), TPM <= NONEXPR_TPM_MAX)
  ctrl_cand <- zero_tpm %>% dplyr::filter(n_expr_ctrl == 0)
  if (nrow(expressed_all) < MIN_BINS) stop("Too few expressed genes in this stratum")
  if (nrow(ctrl_cand) == 0) stop("No non-expressed control genes in unexpressed neighborhoods")
  pool_bp_all <- make_expressed_pools(expressed_all)$pool_bp
  ctrl_genes <- match_by_length(expressed_all$size, ctrl_cand, LEN_MATCH_BINS, SEED)
  elig <- build_eligible_space(all_genes, chrom_sizes, NONGENIC_FLANK_BP, CTRL_NBHD_BP, NULL, "all")
  set.seed(SEED)
  n_win <- ceiling(CTRL_POOLS * pool_bp_all * 1.25 / mean(expressed_all$size))
  n_cap <- floor(0.4 * sum(elig$width) / mean(expressed_all$size))
  n_win <- min(n_win, n_cap)
  win_lens <- sample(expressed_all$size, n_win, replace=TRUE)
  wins <- sample_windows(win_lens, elig, SEED)
  wins$n_x_bins <- floor(wins$end / res) - ceiling(wins$start / res) + 1
  wins <- wins %>% dplyr::filter(n_x_bins >= 1)
  wins <- match_by_length(expressed_all$size, wins, LEN_MATCH_BINS, SEED)
  regions <- dplyr::bind_rows(
    expressed_all %>% dplyr::transmute(region_id=gene_id, chr, start, end, size, n_x_bins, type="expressed"),
    ctrl_genes    %>% dplyr::transmute(region_id=gene_id, chr, start, end, size, n_x_bins, type="nonexpr_gene"),
    wins          %>% dplyr::transmute(region_id, chr, start, end, size, n_x_bins, type="window")
  ) %>% as.data.frame()
  stopifnot(!anyDuplicated(regions$region_id))
  list(regions=regions, expressed_all=expressed_all, ctrl_genes=ctrl_genes)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Cells, pools, delta-S, bootstrap
# ══════════════════════════════════════════════════════════════════════════════
# cell_genes: data.frame(grid, cell, mrow, ccol, gene_id, size)
# -> pool-level delta-S for every pool of every cell
cell_pool_deltas <- function(cg, ex, pool_bp, cfg, seed) {
  cg <- as.data.frame(cg)
  parts <- lapply(split(cg, cg$cell), function(d) {
    x <- assign_pools_by_bp(d, pool_bp, seed)
    if (nrow(x) == 0) return(NULL)
    x
  })
  pm <- dplyr::bind_rows(parts)
  if (nrow(pm) == 0) return(NULL)
  pool_map <- data.frame(region_id=pm$gene_id, kind=pm$cell, pool=pm$pool, stringsAsFactors=FALSE)
  sl <- lapply(names(CONDITIONS), function(n)
    compute_pool_slopes(ex[[n]]$profile, ex[[n]]$meta, pool_map, cfg))
  names(sl) <- names(CONDITIONS)
  d <- merge(sl[[1]][, c("kind", "pool", "slope")], sl[[2]][, c("kind", "pool", "slope")],
             by=c("kind", "pool"), suffixes=c("_1", "_2"))
  d$delta <- d$slope_1 - d$slope_2
  cov <- pm %>% dplyr::group_by(cell, pool) %>%
    dplyr::summarise(grid=grid[1], mrow=mrow[1], ccol=ccol[1],
                     n_genes=dplyr::n(), .groups="drop")
  merge(d[, c("kind", "pool", "slope_1", "slope_2", "delta")], cov,
        by.x=c("kind", "pool"), by.y=c("cell", "pool"))
}

# cell medians from pool-level delta-S
summarise_cells <- function(pl) {
  pl %>% dplyr::group_by(grid, mrow, ccol) %>%
    dplyr::summarise(value=median(delta), n_pools=dplyr::n(), .groups="drop")
}

# ── each variable on its own: expressed genes in equal-bp bins of the variable (as the
#    E/I and TPM bins of scaling_controls_main.R), one delta-S per bin, and the slope of
#    delta-S on the variable's percentile = change in delta-S from the 0th to the 100th percentile
make_single_bins <- function(expr, vars) {
  n <- choose_n_bins(sum(expr$size))
  dplyr::bind_rows(lapply(vars, function(v)
    data.frame(variable=v, bin=dplyr::ntile(expr[[v]], n), gene_id=expr$gene_id, size=expr$size,
               x_p=ecdf(expr[[v]])(expr[[v]]), stringsAsFactors=FALSE)))
}

single_var_deltas <- function(sb, ex, cfg) {
  pool_map <- data.frame(region_id=sb$gene_id, kind=sb$variable, pool=sb$bin, stringsAsFactors=FALSE)
  sl <- lapply(names(CONDITIONS), function(n)
    compute_pool_slopes(ex[[n]]$profile, ex[[n]]$meta, pool_map, cfg))
  d <- merge(sl[[1]][, c("kind", "pool", "slope")], sl[[2]][, c("kind", "pool", "slope")],
             by=c("kind", "pool"), suffixes=c("_1", "_2"))
  d$delta <- d$slope_1 - d$slope_2
  xp <- sb %>% dplyr::group_by(variable, bin) %>% dplyr::summarise(x_p=median(x_p), .groups="drop")
  merge(d, xp, by.x=c("kind", "pool"), by.y=c("variable", "bin"))
}

single_var_effect <- function(dd) {
  do.call(rbind, lapply(split(dd, dd$kind), function(d)
    data.frame(variable=d$kind[1], value=unname(coef(lm(delta ~ x_p, data=d))["x_p"]),
               stringsAsFactors=FALSE)))
}

# median delta-S of the bins in each third of the variable's percentile
tert_medians <- function(dd) {
  dd %>% dplyr::mutate(tert=as.integer(cut(x_p, c(-Inf, 1/3, 2/3, Inf), labels=FALSE))) %>%
    dplyr::group_by(variable=kind, tert) %>%
    dplyr::summarise(value=median(delta), .groups="drop")
}

bootstrap_single <- function(sb, ex, cfg, n_boot=N_BOOT, n_cores=N_CORES) {
  sbt <- data.table::as.data.table(sb)
  one <- function(b) {
    set.seed(SEED + 5000L + b)
    rs <- sbt[, .SD[sample.int(.N, replace=TRUE)], by=.(variable, bin)]
    dd <- single_var_deltas(as.data.frame(rs), ex, cfg)
    list(eff=single_var_effect(dd), tert=tert_medians(dd))
  }
  res <- parallel::mclapply(seq_len(n_boot), function(b) try(one(b), silent=TRUE), mc.cores=n_cores)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (all(bad)) stop("Bootstrap failed in every resample; first error: ", as.character(res[[1]]))
  message("    bootstrap: ", sum(!bad), " of ", n_boot, " resamples completed")
  ok <- res[!bad]
  list(eff=dplyr::bind_rows(lapply(ok, `[[`, "eff")) %>% dplyr::group_by(variable) %>%
         dplyr::summarise(lo=quantile(value, 0.025, na.rm=TRUE), hi=quantile(value, 0.975, na.rm=TRUE), .groups="drop"),
       tert=dplyr::bind_rows(lapply(ok, `[[`, "tert")) %>% dplyr::group_by(variable, tert) %>%
         dplyr::summarise(lo=quantile(value, 0.025, na.rm=TRUE), hi=quantile(value, 0.975, na.rm=TRUE), .groups="drop"))
}

# how noisy is each bin's P(s) curve? a line is fit to log10 P(s) vs log10 distance over the top decade
# of the fit window (10^5 - 10^6 bp for the long genes) and its R2 and residual SD are reported
bin_noise <- function(sb, ex, cfg) {
  pool_map <- data.frame(region_id=sb$gene_id, kind=sb$variable, pool=sb$bin, stringsAsFactors=FALSE)
  lo <- cfg$fit_max / 10
  dplyr::bind_rows(lapply(names(CONDITIONS), function(n) {
    d <- pooled_dist(ex[[n]]$profile, ex[[n]]$meta, pool_map, c("kind", "pool"), cfg)
    d <- d[index >= lo & index < cfg$fit_max]
    n_exp <- length(seq(lo, cfg$fit_max - 1, by=cfg$res))
    res <- lapply(split(d, by=c("kind", "pool"), drop=TRUE), function(x) {
      x <- x[counts > 0]
      if (nrow(x) < 3) return(NULL)
      m <- lm(log10(contact.prob) ~ log10(index), data=x)
      data.frame(variable=x$kind[1], bin=x$pool[1], cond=CONDITIONS[[n]]$display,
                 r2=summary(m)$r.squared, resid_sd=sqrt(mean(residuals(m)^2)),
                 frac_present=nrow(x) / n_exp, stringsAsFactors=FALSE)
    })
    dplyr::bind_rows(res)
  }))
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. One stratum
# ══════════════════════════════════════════════════════════════════════════════
run_chromatin_stratum <- function(genes_s, all_genes, gene_vars, chrom_sizes, cfg,
                                  stratum_i, stratum_label, rds_dir) {
  res <- cfg$res
  chrom_len <- setNames(chrom_sizes$length, chrom_sizes$chr)
  message("\n══ ", ASSAY, "  |  stratum ", stratum_label, "  |  chromatin vs transcription  |  res: ", res, " ══")
  sr <- build_shared_regions(genes_s, all_genes, chrom_sizes, res)
  T_STRAT <- Sys.time()
  tm <- function(what) message("    [time] ", what, ": ", round(as.numeric(difftime(Sys.time(), T_STRAT, units="mins")), 1),
                               " min since stratum start")
  ex <- get_profiles_all(sr$regions, cfg, chrom_len, stratum_i, rds_dir)
  tm("contact profiles loaded")
  
  pass <- Reduce(intersect, lapply(ex, function(e)
    e$meta[n_x_present / n_x_bins >= MIN_COVERAGE]$region_id))
  vnames <- names(VARIABLES)
  attach_vars <- function(df) df %>% dplyr::left_join(gene_vars, by="gene_id")
  expr <- sr$expressed_all %>% dplyr::filter(gene_id %in% pass) %>% attach_vars()
  complete <- stats::complete.cases(expr[, c("loss", vnames)])
  message("  expressed with coverage: ", nrow(expr), "  with all variables (EU-seq + ", length(vnames), " marks): ", sum(complete))
  expr <- expr[complete, ]
  if (nrow(expr) < 9 * MIN_POOLS) stop("Too few expressed genes with all variables (", nrow(expr), ")")
  
  # -- CHECK: how collinear are transcription and chromatin state (expressed genes)?
  cm <- cor(expr[, c("loss", vnames, "TPM")], method="spearman")
  message("  CHECK Spearman correlation among expressed genes (loss = nascent loss, log2 WT/deg):")
  print(round(cm, 2))
  
  # -- tertiles of nascent transcription change (columns of the crossed design)
  loss_q <- tert_cuts(expr$loss)
  expr$ccol <- c("L1", "L2", "L3")[tert_of(expr$loss, loss_q)]
  
  # -- cell membership
  cg <- list()
  for (v in vnames) {
    q <- tert_cuts(expr[[v]])
    er <- tert_of(expr[[v]], q)
    e <- data.frame(grid=v, mrow=er, ccol=expr$ccol, gene_id=expr$gene_id, size=expr$size,
                    stringsAsFactors=FALSE)
    cg[[v]] <- e
    message("  ", v, " tertile cuts: ", paste(signif(q, 3), collapse=" | "))
  }
  cg <- dplyr::bind_rows(cg)
  cg$cell <- paste(cg$grid, cg$mrow, cg$ccol)
  
  # -- one pool size for the whole stratum
  cell_bp <- cg %>% dplyr::group_by(grid, mrow, ccol, cell) %>%
    dplyr::summarise(n_genes=dplyr::n(), bp=sum(size), .groups="drop")
  e_bp <- cell_bp$bp[cell_bp$ccol != "NE"]
  n_bp <- cell_bp$bp[cell_bp$ccol == "NE"]
  pool_bp <- max(POOL_MIN_BP, median(e_bp) / 12)
  if (length(n_bp) > 0) {
    ctrl_bp <- median(n_bp) / (MIN_POOLS * 1.25)
    if (ctrl_bp < pool_bp) {
      pool_bp <- max(ctrl_bp, POOL_MIN_BP / 4)
      message("  pool size lowered to fit the control cells: ", format(round(pool_bp), big.mark=","), " bp")
      if (pool_bp < POOL_MIN_BP)
        warning("Pools of ", format(round(pool_bp), big.mark=","), " bp (< POOL_MIN_BP): noisy", call.=FALSE)
    }
  }
  pool_bp <- round(pool_bp)
  cell_bp$n_pools_possible <- floor(cell_bp$bp / pool_bp)
  drop <- cell_bp$n_pools_possible < MIN_POOLS
  message("  pool size: ", format(pool_bp, big.mark=","), " bp   cells kept: ", sum(!drop), " of ", nrow(cell_bp))
  if (any(drop)) {
    message("  dropped cells (< ", MIN_POOLS, " pools): ",
            paste(cell_bp$cell[drop], "(", cell_bp$n_genes[drop], " genes, ", cell_bp$n_pools_possible[drop], " pools)",
                  sep="", collapse="; "))
    cg <- cg[!cg$cell %in% cell_bp$cell[drop], ]
  }
  if (!any(cg$ccol != "NE")) stop("No expressed cell has enough sequence for ", MIN_POOLS, " pools")
  message("  CHECK genes / pools per cell:")
  print(as.data.frame(cell_bp[, c("cell", "n_genes", "bp", "n_pools_possible")]))
  
  # -- observed pool-level delta-S
  keep_ids <- unique(cg$gene_id)
  exk <- lapply(ex, function(e) list(profile=e$profile[region_id %in% keep_ids],
                                     meta=e$meta[region_id %in% keep_ids]))
  pl <- cell_pool_deltas(cg, exk, pool_bp, cfg, SEED)
  if (is.null(pl)) stop("No pools could be formed")
  cells <- summarise_cells(pl)
  tm("cell pools done")
  message("  CHECK pooled S (median over pools): ", COND_LEVELS[1], "=", round(median(pl$slope_1), 3),
          "  ", COND_LEVELS[2], "=", round(median(pl$slope_2), 3),
          "  | share of negative slopes: ", round(mean(pl$slope_1 < 0 & pl$slope_2 < 0), 2))
  
  # -- each variable on its own (expressed genes): chromatin variables and nascent loss alike
  vars_all <- c(vnames, "loss")
  sb <- make_single_bins(expr, vars_all)
  ex_s <- lapply(ex, function(e) list(profile=e$profile[region_id %in% expr$gene_id],
                                      meta=e$meta[region_id %in% expr$gene_id]))
  obs1 <- single_var_deltas(sb, ex_s, cfg)
  message("  Single-variable analysis: ", length(unique(sb$bin)), " equal-bp bins per variable, median bp per bin: ",
          format(round(median(tapply(sb$size, paste(sb$variable, sb$bin), sum))), big.mark=","))
  message("  CHECK pooled S (median over bins): ", COND_LEVELS[1], "=", round(median(obs1$slope_1), 3),
          "  ", COND_LEVELS[2], "=", round(median(obs1$slope_2), 3),
          "  | share of negative slopes: ", round(mean(obs1$slope_1 < 0 & obs1$slope_2 < 0), 2))
  # -- third analysis: is gene length driving the single-variable results?
  len_bin <- sb %>% dplyr::group_by(variable, bin) %>% dplyr::summarise(median_size=median(size), .groups="drop")
  obs1l <- merge(obs1, len_bin, by.x=c("kind", "pool"), by.y=c("variable", "bin"))
  length_check <- do.call(rbind, lapply(vars_all, function(v) {
    d <- obs1l[obs1l$kind == v, ]
    ct <- suppressWarnings(cor.test(d$delta, d$median_size, method="spearman"))
    data.frame(assay=ASSAY, analysis=ANALYSIS, stratum=stratum_label,
               variable=EFFECT_LABS[[v]],
               rho_variable_vs_log10_gene_length=unname(cor(expr[[v]], log10(expr$size), method="spearman")),
               n_bins=nrow(d),
               rho_bin_deltaS_vs_median_gene_length=unname(ct$estimate), p_value=ct$p.value,
               stringsAsFactors=FALSE)
  }))
  message("  CHECK gene length: Spearman rho of each variable with log10 gene length (expressed genes), and of bin delta-S with the bin's median gene length:")
  print(length_check[, -(1:3)], digits=2, row.names=FALSE)
  
  noise <- bin_noise(sb, ex_s, cfg)
  message("  CHECK curve noise per bin: line fit to log10 P(s) over ", format(cfg$fit_max / 10, scientific=FALSE),
          " - ", format(cfg$fit_max, scientific=FALSE), " bp (R2 and residual SD in log10 units)")
  print(as.data.frame(noise %>% dplyr::group_by(variable, cond) %>%
                        dplyr::summarise(n_bins=dplyr::n(), median_R2=round(median(r2), 3), min_R2=round(min(r2), 3),
                                         median_resid_sd=round(median(resid_sd), 4), max_resid_sd=round(max(resid_sd), 4),
                                         median_frac_present=round(median(frac_present), 2), .groups="drop")))
  message("  Bootstrapping single-variable effects: ", N_BOOT, " gene resamples within bins on ", N_CORES, " cores")
  t1 <- Sys.time()
  ci1 <- bootstrap_single(sb, ex_s, cfg)
  tm("single-variable bootstrap done")
  message("    done in ", round(as.numeric(difftime(Sys.time(), t1, units="mins")), 1), " min")
  eff1 <- dplyr::left_join(single_var_effect(obs1), ci1$eff, by="variable") %>%
    dplyr::transmute(grid=variable, stat="single_expressed", value, lo, hi)
  coefs <- eff1
  message("  CHECK observed effect inside its bootstrap interval: ",
          round(mean(coefs$value >= coefs$lo & coefs$value <= coefs$hi, na.rm=TRUE), 2))
  message("  Effects on delta-S (change from the 0th to the 100th percentile; 95% gene-bootstrap). single_expressed = each variable on its own:")
  print(as.data.frame(coefs), digits=3)
  
  # each variable on its own, in thirds of its percentile: median delta-S of the bins in each third
  tert <- dplyr::left_join(tert_medians(obs1), ci1$tert, by=c("variable", "tert"))
  list(label=stratum_label, cells=cells, tert=tert, coefs=coefs, noise=noise, length_check=length_check, length_bins=obs1l, pool_level=pl, cell_bp=cell_bp,
       pool_bp=pool_bp, cor=cm)
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Plots
# ══════════════════════════════════════════════════════════════════════════════
heat_prep <- function(cells, ccols) {
  full <- expand.grid(grid=names(VARIABLES), mrow=1:3, ccol=ccols, stringsAsFactors=FALSE)
  d <- dplyr::left_join(full, cells[, c("grid", "mrow", "ccol", "value")], by=c("grid", "mrow", "ccol"))
  d$grid <- factor(d$grid, levels=names(VARIABLES), labels=unname(EFFECT_LABS[names(VARIABLES)]))
  d$mrow <- factor(d$mrow, levels=1:3, labels=c("T1", "T2", "T3"))
  d
}

heat_theme <- function() {
  theme_classic(base_size=11)+
    theme(plot.title=element_text(face="bold", size=13), axis.text=element_text(color="black"),
          axis.line=element_blank(), strip.background=element_blank(),
          strip.text=element_text(face="bold"))
}

# analysis "state": one facet per mark; x = tertile of the WT mark level, y = tertile of nascent transcription change
plot_heat_main <- function(r, tag, zlim) {
  d <- heat_prep(r$cells, c("L1", "L2", "L3"))
  d$ccol <- factor(d$ccol, levels=c("L1", "L2", "L3"), labels=c("T1 (least)", "T2", "T3 (largest)"))
  ggplot(d, aes(mrow, ccol, fill=value))+
    geom_tile(color="white", linewidth=0.8)+
    geom_text(data=d[!is.na(d$value), ], aes(label=sprintf("%.3f", value),
                                             color=ifelse(value > mean(zlim), "black", "white")), size=3)+
    scale_color_identity()+
    scale_fill_viridis_c(option="C", limits=zlim, name="ΔS", na.value="gray90")+
    facet_grid(. ~ grid)+
    coord_fixed(ratio=0.55)+
    labs(x=if (MODE == "change") "Tertile of the mark's change (T1 = least decrease or increase, T3 = largest decrease)"
         else "Tertile of the Wildtype level (low to high)",
         y=paste0("Tertile of ", tolower(LOSS_NAME), " decrease"),
         title=paste0(paste0(sub(" \\(.*$", "", MODE_TITLE), if (MODE == "change") " (crossed)" else ""), ": ΔS (", COND_LEVELS[1], " − ",
                      COND_LEVELS[2], ")  |  ", tag))+
    heat_theme()
}

# analysis "change": one row per variable, x = tertile of that variable alone
plot_heat_change <- function(r, tag, zlim) {
  d <- r$tert
  d$variable <- factor(EFFECT_LABS[d$variable], levels=rev(unname(EFFECT_LABS[c("loss", names(VARIABLES))])))
  d$tert <- factor(d$tert, levels=1:3, labels=c("T1", "T2", "T3"))
  ggplot(d, aes(tert, variable, fill=value))+
    geom_tile(color="white", linewidth=0.8)+
    geom_text(aes(label=sprintf("%.3f", value), color=ifelse(value > mean(zlim), "black", "white")), size=3)+
    scale_color_identity()+
    scale_fill_viridis_c(option="C", limits=zlim, name="ΔS", na.value="gray90")+
    coord_fixed(ratio=0.55)+
    labs(x="Tertile of the variable's own change (T1 = least decrease or increase, T3 = largest decrease)", y=NULL,
         title=paste0(sub(" \\(.*$", "", MODE_TITLE), " (each variable alone): ΔS (", COND_LEVELS[1], " − ", COND_LEVELS[2], ")  |  ", tag))+
    heat_theme()
}

# ── effect-size table: one row per stratum x term x variable ───────────────────
EFFECT_TERMS <- c(single_expressed = "Each variable on its own")
LOSS_NAME    <- "Nascent transcription"   # reset by run_chromatin_analysis()
EFFECT_LABS  <- c(H3K27ac="H3K27ac", H3K4me3="H3K4me3", H3K9me3="H3K9me3",
                  loss="Nascent transcription")   # reset by run_chromatin_analysis(mode)

# rows: stratum x term x variable
build_effect_table <- function(results, terms=EFFECT_TERMS, vars=c(names(VARIABLES), "loss")) {
  out <- list()
  for (a in names(results)) for (si in seq_along(results[[a]])) {
    r <- results[[a]][[si]]
    full <- expand.grid(grid=vars, stat=names(terms), stringsAsFactors=FALSE)
    d <- dplyr::left_join(full, r$coefs[, c("grid", "stat", "value", "lo", "hi")], by=c("grid", "stat"))
    out[[length(out) + 1]] <- data.frame(
      assay=ASSAY, analysis=a, stratum=r$label, variable=EFFECT_LABS[d$grid],
      term=unname(terms[d$stat]), effect=d$value, ci_low=d$lo, ci_high=d$hi,
      estimable=!is.na(d$value), stringsAsFactors=FALSE)
  }
  tab <- dplyr::bind_rows(out)
  tab$term <- factor(tab$term, levels=terms)
  tab[order(tab$analysis, tab$term, match(tab$variable, EFFECT_LABS)), ]
}

strat_labels <- function(tab) {
  lab <- paste0(ifelse(tab$analysis == "long", "Long genes ", "Short genes "), tab$stratum)
  factor(lab, levels=unique(lab))
}

# change in delta-S per variable, one panel for expressed genes and one for the controls;
# nascent loss is shown as a variable of its own (red)
plot_effect_table <- function(tab) {
  d <- tab[tab$estimable, ]
  d$strat_lab <- factor(strat_labels(d), levels=levels(strat_labels(tab)))
  d$variable <- factor(d$variable, levels=rev(unname(EFFECT_LABS)))
  d$grp <- ifelse(d$variable == EFFECT_LABS[["loss"]], "Nascent", "Chromatin")
  ggplot(d, aes(effect, variable, color=grp))+
    geom_vline(xintercept=0, linetype="dashed", color="gray40")+
    geom_pointrange(aes(xmin=ci_low, xmax=ci_high), size=0.5)+
    facet_grid(strat_lab ~ term, labeller=labeller(term=label_wrap_gen(28)))+
    scale_color_manual(values=c("Chromatin"="#4575B4", "Nascent"="#D73027"), guide="none")+
    labs(x="Change in ΔS from the 0th to the 100th percentile (95% gene bootstrap)", y=NULL,
         title=paste0(MODE_TITLE, ": effect of each variable on ΔS  |  ", ASSAY))+
    theme_classic(base_size=11)+
    theme(plot.title=element_text(face="bold", size=13), axis.text=element_text(color="black"),
          strip.background=element_blank(), strip.text=element_text(face="bold", size=9))
}

# the crossed design as a line plot: x = tertile of the variable, one line per tertile of nascent
# transcription change; line = cell median, dots = the pools behind each cell
plot_cells_xy <- function(r) {
  pl <- r$pool_level
  pl <- pl[pl$ccol %in% c("L1", "L2", "L3"), ]
  lv <- c("L1", "L2", "L3"); lb <- c("T1 (least)", "T2", "T3 (largest)")
  pl$grid <- factor(pl$grid, levels=names(VARIABLES), labels=unname(EFFECT_LABS[names(VARIABLES)]))
  pl$ccol <- factor(pl$ccol, levels=lv, labels=lb)
  cm <- r$cells[r$cells$ccol %in% lv, ]
  cm$grid <- factor(cm$grid, levels=names(VARIABLES), labels=unname(EFFECT_LABS[names(VARIABLES)]))
  cm$ccol <- factor(cm$ccol, levels=lv, labels=lb)
  ggplot(cm, aes(mrow, value, color=ccol, group=ccol))+
    geom_point(data=pl, aes(mrow, delta), position=position_dodge(width=0.45), size=0.9, alpha=0.35)+
    geom_line(linewidth=0.8)+
    geom_point(size=2.2)+
    scale_x_continuous(breaks=1:3, labels=c("T1", "T2", "T3"))+
    scale_color_manual(values=c("#9ECAE1", "#4292C6", "#08306B"), name=paste0(LOSS_NAME, "\ndecrease tertile"))+
    facet_grid(. ~ grid)+
    labs(x=if (MODE == "change") "Tertile of the mark's change (T1 = least decrease or increase, T3 = largest decrease)"
         else "Tertile of the Wildtype level (low to high)",
         y="ΔS (cell median; dots = pools)")+
    heat_theme()
}

# each variable on its own, as a line plot: x = tertile of the variable's own change, y = delta-S;
# dots = the pools (bins) in that third, line = their median (the values in the one-row heatmap);
# nascent transcription change is one of the panels
plot_tertile_lines <- function(r, tag) {
  lv <- c("loss", names(VARIABLES))
  d <- r$length_bins
  d$tert <- as.integer(cut(d$x_p, c(-Inf, 1/3, 2/3, Inf), labels=FALSE))
  d$variable <- factor(EFFECT_LABS[d$kind], levels=unname(EFFECT_LABS[lv]))
  m <- r$tert
  m$variable <- factor(EFFECT_LABS[m$variable], levels=unname(EFFECT_LABS[lv]))
  d$grp <- ifelse(d$kind == "loss", "Nascent", "Chromatin"); m$grp <- ifelse(m$variable == EFFECT_LABS[["loss"]], "Nascent", "Chromatin")
  ggplot(m, aes(tert, value, color=grp, group=variable))+
    geom_point(data=d, aes(tert, delta), position=position_jitter(width=0.12, height=0, seed=1), size=1.2, alpha=0.35)+
    geom_linerange(aes(ymin=lo, ymax=hi), linewidth=0.7)+
    geom_line(linewidth=0.9)+
    geom_point(size=2.6)+
    scale_x_continuous(breaks=1:3, labels=c("T1", "T2", "T3"))+
    scale_color_manual(values=c("Chromatin"="#4575B4", "Nascent"="#D73027"), guide="none")+
    facet_grid(. ~ variable, labeller=labeller(variable=label_wrap_gen(22)))+
    labs(x="Tertile of the variable's own change (T1 = least decrease or increase, T3 = largest decrease)",
         y="ΔS (median with 95% gene-bootstrap interval; dots = pools)",
         title=paste0(sub(" \\(.*$", "", MODE_TITLE), " (each variable alone)  |  ", tag))+
    heat_theme()
}

# gene-length check as a plot: every point is one bin of the single-variable analysis (x = median gene
# length of the bin, y = the bin's ΔS, colour = the bin's percentile of the variable); the line is a
# linear fit; the text in each panel gives the Spearman correlation between the variable and gene length,
# and between the bins' ΔS and their median gene length
plot_gene_length <- function(results) {
  d <- dplyr::bind_rows(lapply(names(results), function(a) dplyr::bind_rows(lapply(results[[a]], function(r) {
    x <- r$length_bins
    x$strat <- paste0(ifelse(a == "long", "Long genes ", "Short genes "), r$label)
    x$variable <- EFFECT_LABS[x$kind]
    lc <- r$length_check
    x$rho_var <- lc$rho_variable_vs_log10_gene_length[match(x$variable, lc$variable)]
    x$rho_bin <- lc$rho_bin_deltaS_vs_median_gene_length[match(x$variable, lc$variable)]
    x$p_bin   <- lc$p_value[match(x$variable, lc$variable)]
    x
  }))))
  d$strat <- factor(d$strat, levels=unique(d$strat))
  d$variable <- factor(d$variable, levels=unname(EFFECT_LABS))
  ann <- d %>% dplyr::group_by(strat, variable) %>%
    dplyr::summarise(x=min(median_size), rho_var=rho_var[1], rho_bin=rho_bin[1], p_bin=p_bin[1], .groups="drop")
  ann$lab <- sprintf("ρ(variable, length) = %.2f\nρ(bin ΔS, length) = %.2f\np = %.2g", ann$rho_var, ann$rho_bin, ann$p_bin)
  # one row per stratum, each with its own x-axis
  one <- function(st, last) {
    ggplot(d[d$strat == st, ], aes(median_size, delta))+
      geom_point(aes(color=x_p), size=1.8)+
      geom_smooth(method="lm", se=FALSE, color="gray30", linewidth=0.6, formula=y ~ x)+
      scale_x_log10(breaks=scales::breaks_log(n=3), labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
      scale_color_viridis_c(option="D", name="Percentile of\nthe variable", limits=c(0, 1))+
      geom_text(data=ann[ann$strat == st, ], aes(x=x, y=Inf, label=lab), hjust=0, vjust=1.15, size=2.6, lineheight=0.95, inherit.aes=FALSE)+
      scale_y_continuous(expand=expansion(mult=c(0.05, 0.35)))+
      facet_wrap(~variable, nrow=1, scales="free", labeller=labeller(variable=label_wrap_gen(22)))+
      labs(x=if (last) "Median gene length of the bin (bp, log scale)" else NULL, y="ΔS of the bin", title=as.character(st))+
      theme_classic(base_size=11)+
      theme(plot.title=element_text(face="bold", size=10), axis.text=element_text(color="black"),
            strip.background=element_blank(), strip.text=element_text(face="bold", size=9))
  }
  sl <- levels(d$strat)
  cowplot::plot_grid(plotlist=lapply(seq_along(sl), function(i) one(sl[i], i == length(sl))), ncol=1)
}

# Save one plot as SVG in PLOT_DIR (no-op when PLOT_DIR is NULL; PLOT_DIR, PLOT_W, PLOT_H come from
# scaling_controls_main.R). File name: <assay>_<mode><tag>_<analysis>_<stratum>_<plot type>.svg
save_chromatin_svg <- function(p, type, analysis=NULL, label=NULL) {
  if (is.null(PLOT_DIR)) return(invisible(NULL))
  parts <- c(ASSAY, paste0(MODE, PLOT_TAG), analysis, if (!is.null(label)) gsub(">=", "ge", label), type)
  f <- file.path(PLOT_DIR, paste0(gsub("[^A-Za-z0-9._-]+", "-", paste(parts, collapse="_")), ".svg"))
  svg(f, width=PLOT_W, height=PLOT_H)
  print(p)
  dev.off()
  message("  saved ", f)
  invisible(f)
}

print_chromatin_plots <- function(r) {
  tag <- paste0(ASSAY, "  |  ", ANALYSIS, "  |  ", r$label)
  if (MODE == "change") {
    p_heat <- plot_heat_change(r, tag, range(r$tert$value, na.rm=TRUE))
    print(p_heat)
    save_chromatin_svg(p_heat, "heatmap-each-variable", ANALYSIS, r$label)
    p_lines <- plot_tertile_lines(r, tag)
    print(p_lines)
    save_chromatin_svg(p_lines, "lines-each-variable", ANALYSIS, r$label)
  }
  # crossed heatmap with the same data as a line plot below it
  p_cross <- cowplot::plot_grid(plot_heat_main(r, tag, range(r$cells$value, na.rm=TRUE)), plot_cells_xy(r),
                                ncol=1, align="v", axis="lr", rel_heights=c(1, 1.25))
  print(p_cross)
  save_chromatin_svg(p_cross, "heatmap-and-lines-crossed", ANALYSIS, r$label)
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Driver
# ══════════════════════════════════════════════════════════════════════════════
# nascent_fun / tag / title / loss_name let another script swap the transcription variable (see
# scaling_chromatin_vs_neighborhood_transcription.R): nascent_fun(all_genes) must return
# data.frame(gene_id, loss); tag is added to the output file names.
init_mode_globals <- function(mode, loss_name, title=NULL) {
  MODE <<- mode
  LOSS_NAME <<- loss_name
  MODE_TITLE <<- if (!is.null(title)) title
  else if (mode == "change") "Analysis 1 (change in transcription vs change in marks)"
  else "Analysis 2 (Wildtype chromatin state + transcription change)"
  VARIABLES <<- if (mode == "change") ALL_VARIABLES[!sapply(ALL_VARIABLES, function(v) is.null(v$files_deg))]
  else ALL_VARIABLES
  mark_lab <- c(POL2="Pol II", H3K27ac="H3K27ac", H3K4me3="H3K4me3", H3K9me3="H3K9me3")
  EFFECT_LABS <<- c(setNames(trimws(gsub(" +", " ", paste(if (mode == "change") "" else "WT", mark_lab[names(VARIABLES)],
                                                          if (mode == "change") "change" else ""))), names(VARIABLES)),
                    loss=paste(loss_name, "change"))
  COND_LEVELS <<- sapply(CONDITIONS, `[[`, "display")
  COND_COLORS <<- setNames(sapply(CONDITIONS, `[[`, "color"), COND_LEVELS)
}

run_chromatin_analysis <- function(mode=c("change", "state"), nascent_fun=NULL, tag="", title=NULL,
                                   loss_name="Nascent transcription") {
  mode <- match.arg(mode)
  init_mode_globals(mode, loss_name, title)
  PLOT_TAG <<- tag
  message("\n############ ", MODE_TITLE, " ############")
  
  # the DMSO and IAA bigWig of each mark
  vfiles <- lapply(VARIABLES, function(v) if (mode == "change") c(wt=v$files, deg=v$files_deg) else c(wt=v$files))
  # every input must exist before anything runs
  input_files <- c(setNames(sapply(CONDITIONS, `[[`, "hic"), paste0(names(CONDITIONS), "_hic")),
                   gtf=GTF_FILE, expression=EXPR_FILE, euseq=EUSEQ_FILE)
  missing_files <- input_files[!file.exists(input_files)]
  miss_v <- names(vfiles)[!sapply(vfiles, function(f) all(file.exists(f)))]
  if (length(missing_files) > 0 || length(miss_v) > 0) {
    msg <- character(0)
    if (length(missing_files) > 0)
      msg <- c(msg, paste0("  ", names(missing_files), ": ", missing_files))
    for (v in miss_v) {
      dirs <- unique(dirname(vfiles[[v]]))
      lst  <- unlist(lapply(dirs[dir.exists(dirs)], function(d) head(list.files(d, pattern="\\.(bw|bigWig)$"), 40)))
      msg <- c(msg, paste0("  ", v, ": missing ", paste(vfiles[[v]][!file.exists(vfiles[[v]])], collapse=" | "),
                           if (length(lst)) paste0("\n      bigWigs in those directories: ", paste(lst, collapse=", "))))
    }
    stop("Missing input files (edit EUSEQ_FILE / DATA_DIR at the top of the script):\n",
         paste(msg, collapse="\n"))
  }
  if (!is.null(PLOT_DIR) && !dir.exists(PLOT_DIR)) stop("PLOT_DIR not found: ", PLOT_DIR)
  if (!requireNamespace("rtracklayer", quietly=TRUE)) stop("rtracklayer is required to read the bigWigs")
  message("All inputs found. Signals (", if (mode == "change") "DMSO / IAA" else "DMSO", "): ",
          paste(names(vfiles), sapply(vfiles, function(f) paste(basename(f), collapse=" / ")), sep="=", collapse="  "))
  message(ASSAY, "  |  analyses: ", paste(ANALYSES, collapse=", "), "  |  cores: ", N_CORES)
  dir.create(RDS_BASE, showWarnings=FALSE, recursive=TRUE)
  
  cs <- lapply(CONDITIONS, function(cd) {
    ct <- strawr::readHicChroms(cd$hic)
    mp <- hic_chrom_map(cd$hic)
    data.frame(chr=names(mp), length=ct$length[match(mp, ct$name)], stringsAsFactors=FALSE)
  })
  chrom_sizes <- cs[[1]]
  for (k in cs) stopifnot(identical(k$length, chrom_sizes$length))
  
  all_genes <- load_gene_table(GTF_FILE, RDS_BASE) %>%
    dplyr::left_join(load_expression(EXPR_FILE), by="gene_id")
  all_genes <- add_neighborhood(all_genes, NBHD_BP, NONEXPR_TPM_MAX)
  all_genes$n_expr_ctrl <- add_neighborhood(all_genes, CTRL_NBHD_BP, CTRL_TPM_MAX)$n_expr_nbhd
  ei_df <- all_genes %>% dplyr::filter(exon_bp > 0, !is.infinite(EI), EI > 0, EI != 1)
  
  strand_df <- load_gene_strand(GTF_FILE, RDS_BASE)
  sigs <- load_gene_signals(all_genes, strand_df, vfiles, RDS_BASE, mode)
  eu   <- if (is.null(nascent_fun)) load_euseq_loss(EUSEQ_FILE) else nascent_fun(all_genes)
  n_match <- sum(ei_df$gene_id[!is.na(ei_df$TPM) & ei_df$TPM > 0] %in% eu$gene_id)
  n_expr  <- sum(!is.na(ei_df$TPM) & ei_df$TPM > 0)
  message("  EU-seq genes: ", nrow(eu), "  matched to expressed analysis genes: ", n_match, " of ", n_expr)
  if (n_match < 0.5 * n_expr) stop("Fewer than half of the expressed genes are in the EU-seq table: gene_id format mismatch?")
  gene_vars <- dplyr::left_join(sigs, eu, by="gene_id")
  message("  CHECK signal summary (all genes; zero and NA fractions):")
  print(do.call(rbind, lapply(names(VARIABLES), function(v)
    data.frame(var=v, median=signif(median(gene_vars[[v]], na.rm=TRUE), 3),
               frac_zero=round(mean(gene_vars[[v]] == 0, na.rm=TRUE), 2),
               frac_NA=round(mean(is.na(gene_vars[[v]])), 2)))))
  
  results <- list()
  for (a in ANALYSES) {
    ANALYSIS <<- a
    cfg <- ANALYSIS_CFG[[a]]
    message("\n════════ ", ASSAY, " ANALYSIS: ", a, " ════════")
    rds_dir <- file.path(RDS_BASE, a)
    dir.create(rds_dir, showWarnings=FALSE, recursive=TRUE)
    for (n in names(CONDITIONS)) check_hic_access(CONDITIONS[[n]]$hic, cfg, chrom_sizes)
    results[[a]] <- list()
    for (si in seq_along(cfg$strata)) {
      lo <- cfg$strata[[si]][1]; hi <- cfg$strata[[si]][2]
      label <- if (is.infinite(hi)) paste0(">=", fmt_bp(lo)) else paste0(fmt_bp(lo), "-", fmt_bp(hi))
      genes_s <- ei_df %>% dplyr::filter(size >= lo, size < hi)
      results[[a]][[si]] <- run_chromatin_stratum(genes_s, all_genes, gene_vars, chrom_sizes, cfg,
                                                  si, label, rds_dir)
      print_chromatin_plots(results[[a]][[si]])
    }
  }
  saveRDS(results, file.path(RDS_BASE, paste0("chromatin_vs_transcription_", mode, tag, "_", ASSAY, ".rds")))
  
  # -- effect-size table: CSV and plot
  tab <- build_effect_table(results)
  csv_file <- file.path(RDS_BASE, paste0("effect_sizes_", mode, tag, "_", ASSAY, ".csv"))
  readr::write_csv(tab, csv_file)
  message("\nEffect-size table written to ", csv_file)
  print(as.data.frame(tab), digits=3)
  p_eff <- plot_effect_table(tab)
  print(p_eff)
  save_chromatin_svg(p_eff, "effect-sizes", "all-strata")
  
  # -- gene-length check: plot, table and CSV
  tab3 <- dplyr::bind_rows(lapply(results, function(ra) dplyr::bind_rows(lapply(ra, `[[`, "length_check"))))
  csv3 <- file.path(RDS_BASE, paste0("gene_length_check_", mode, tag, "_", ASSAY, ".csv"))
  readr::write_csv(tab3, csv3)
  message("\nGene-length check written to ", csv3)
  print(as.data.frame(tab3), digits=2)
  p_len <- plot_gene_length(results)
  print(p_len)
  save_chromatin_svg(p_len, "gene-length-check", "all-strata")
  invisible(results)
}

# (another script can set options(chromatin_script_no_run=TRUE) before source()-ing this file to reuse the functions)
if (!isTRUE(getOption("chromatin_script_no_run"))) {
  # every plot is saved here as an SVG (folder must already exist); size in inches, same for every plot
  PLOT_DIR <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/micro-c/plots/scaling_v_transcription_plots"
  PLOT_W   <- 14
  PLOT_H   <- 10
  results_change <- run_chromatin_analysis("change")   # analysis 1: change in transcription vs change in marks
  results_state  <- run_chromatin_analysis("state")    # analysis 2: Wildtype chromatin state + transcription change
}