##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Shared functions for the contact-scaling control analysis (reviewer
# comment: "show Delta-S for length-matched, non-expressed control
# regions"). Adapted from scaling_analysis.R: pooled P(s) per bin,
# S = slope of log10(P(s)) vs. log10(distance), the same delta-S panels
# and the same loess S(d) panel, with controls added.
#
# This file defines functions and shared constants only. Run one of:
#   scaling_controls_hic.R         Hi-C,     genes >=100kb
#   scaling_controls_microc.R      Micro-C,  genes >=100kb and 1kb-100kb
#   scaling_controls_microc_AB.R   Micro-C,  A and B compartments
# They load this file from their own directory.
#
# Groups
#   Expressed: genes with TPM > 0, pooled into bins of exon/intron ratio
#     and bins of neighborhood TPM (mean TPM of the gene and every gene
#     within NBHD_BP of it); Q1 / Q4 of neighborhood TPM for S(d).
#   Controls (genome-wide and shared by the all-compartment, A and B runs;
#     length-matched to the expressed genes, pooled to the same bp per pool,
#     and drawn ONLY from unexpressed neighborhoods, i.e. no gene with
#     TPM > CTRL_TPM_MAX within CTRL_NBHD_BP):
#       non-expressed genes (TPM <= NONEXPR_TPM_MAX)
#       non-genic windows (no annotated gene +/- NONGENIC_FLANK_BP)
##################################################################
suppressPackageStartupMessages({
  library(strawr)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(parallel)
})

# ══════════════════════════════════════════════════════════════════════════════
# Shared constants
# ══════════════════════════════════════════════════════════════════════════════
STRAW_NORM <- "SCALE"
STD_CHROMS <- paste0("chr", c(1:22, "X"))
get_slurm_cores <- function(default=2L) {
  for (v in c("SLURM_CPUS_PER_TASK", "SLURM_NTASKS", "SLURM_CPUS_ON_NODE")) {
    x <- suppressWarnings(as.integer(Sys.getenv(v)))
    if (!is.na(x) && x > 0L) { message(v, ": ", x, " cores"); return(x) }
  }
  message("No SLURM allocation variables found - using ", default, " cores")
  default
}
N_CORES    <- get_slurm_cores()

HG38_REF  <- "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference"
GTF_FILE  <- file.path(HG38_REF, "Homo_sapiens.GRCh38.115.chr.gtf")  # Ensembl 115, chr-prefixed
EXPR_FILE <- file.path(HG38_REF, "ENCFF435PHM.tsv")                  # ENCODE HCT116 RSEM (gene_id, TPM, ...)

POOL_MIN_BP       <- 2000000L  # min total bp per expressed bin (sets the number of bins)
MIN_BINS          <- 10L
MAX_BINS          <- 100L
CTRL_POOLS        <- 30L       # control pools per control type
MIN_POOLS         <- 5L        # minimum pools required per group
LEN_MATCH_BINS    <- 20L       # length-matching bins (quantiles of expressed lengths)
NONGENIC_FLANK_BP <- 10000L    # flank excluded around every annotated gene
NBHD_BP           <- 100000L   # neighborhood for the expressed-gene neighborhood TPM (gene +/- this distance)
CTRL_NBHD_BP      <- 100000L   # controls: no gene with TPM > CTRL_TPM_MAX within this distance (same scale as NBHD_BP)
CTRL_TPM_MAX      <- 0         # controls: tolerated TPM of genes in the control neighborhood
NBHD_BANDS_BP     <- c(0, 1e5, 2.5e5, 5e5, 1e6, 2e6)  # flanking distance bands for the neighborhood-activity plot
MIN_COVERAGE      <- 0.8       # min fraction of a region's bins with contacts, both conditions
NONEXPR_TPM_MAX   <- 0         # unexpressed: TPM <= this
SEED              <- 42L
N_BOOT            <- 200L      # region resamples for the delta-S(d) bands
PLOT_DIR          <- NULL      # set in the wrapper script to save every plot as SVG (folder must exist)
PLOT_W            <- 12        # saved SVG size (inches), same for every plot
PLOT_H            <- 9
PLOT_H_LOESS      <- 5         # height of the side-by-side S(d) / delta-S(d) plot (its panels have a fixed aspect ratio)

ANALYSIS_CFG <- list(
  long = list(
    res=5000L, strata=list(c(1e5, Inf)), extract="chrom", band_only=FALSE,
    fit_min=5000, fit_max=1e6, profile_max=5e6, possible_mode="original",
    loess_lo=5000, loess_hi=4e6, delta_ylim=c(-0.5, 0.5),
    axis_kb=c(100, 250, 500, 1000, 2000, 3000, 4000),
    axis_lab=c("100kb","250kb","500kb","1Mb","2Mb","3Mb","4Mb")),
  short = list(
    res=500L, strata=list(c(1000, 1e5)), extract="tiles", tile_bp=1e6, band_only=TRUE,
    fit_min=1000, fit_max=1e5, profile_max=1e5, possible_mode="constant",
    loess_lo=1000, loess_hi=1e5, delta_ylim=NULL,
    axis_kb=c(1, 2, 5, 10, 20, 50, 100),
    axis_lab=c("1kb","2kb","5kb","10kb","20kb","50kb","100kb"))
)

GROUP_LEVELS <- c("Expressed, Low E/I", "Expressed, High E/I",
                  "Expressed, Low nbhd TPM", "Expressed, High nbhd TPM",
                  "Non-expressed genes", "Non-genic windows")
SD_LEVELS    <- c("Expressed, Q1 (Low)", "Expressed, Q4 (High)",
                  "Non-expressed genes", "Non-genic windows")
SD_COLORS    <- setNames(c("#4575B4", "#D73027", "#7B3294", "#636363"), SD_LEVELS)
# groups still analysed but not drawn: the intergenic windows are the control shown in the figures
PLOT_DROP_GROUPS <- "Non-expressed genes"
PLOT_LEVELS  <- setdiff(SD_LEVELS, PLOT_DROP_GROUPS)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Annotation, expression, compartments
# ══════════════════════════════════════════════════════════════════════════════
load_gene_table <- function(gtf_file, rds_dir) {
  gene_rds <- file.path(rds_dir, "gene_table.rds")
  if (file.exists(gene_rds)) {
    message("Loading cached gene table...")
    return(readRDS(gene_rds))
  }
  message("Building TxDb from GTF...")
  txdb     <- makeTxDbFromGFF(gtf_file, format="gtf")
  exons_gr <- exons(txdb) %>% GenomicRanges::reduce()
  genes_gr <- genes(txdb)
  seqlevelsStyle(genes_gr) <- "UCSC"
  seqlevelsStyle(exons_gr) <- "UCSC"
  genes_gr <- keepSeqlevels(genes_gr,
                            intersect(seqlevels(genes_gr), STD_CHROMS),
                            pruning.mode="coarse")
  exons_gr <- keepSeqlevels(exons_gr,
                            intersect(seqlevels(exons_gr), STD_CHROMS),
                            pruning.mode="coarse")
  genes_gr$gene_id <- sub("\\..*", "", names(genes_gr))
  
  hits     <- findOverlaps(genes_gr, exons_gr)
  overlaps <- GenomicRanges::pintersect(genes_gr[queryHits(hits)],
                                        exons_gr[subjectHits(hits)])
  exon_bp_per_gene <- tapply(width(overlaps), queryHits(hits), sum)
  exon_bp_vec <- rep(0L, length(genes_gr))
  exon_bp_vec[as.integer(names(exon_bp_per_gene))] <- as.integer(exon_bp_per_gene)
  
  gene_df <- data.frame(
    gene_id = genes_gr$gene_id,
    chr     = as.character(seqnames(genes_gr)),
    start   = start(genes_gr),
    end     = end(genes_gr),
    gene_bp = width(genes_gr),
    exon_bp = exon_bp_vec,
    stringsAsFactors=FALSE
  ) %>%
    dplyr::mutate(
      intron_bp = gene_bp - exon_bp,
      EI        = exon_bp / pmax(intron_bp, 1),
      size      = end - start
    )
  saveRDS(gene_df, gene_rds)
  message("Saved gene table: ", nrow(gene_df), " genes")
  gene_df
}

load_expression <- function(expr_file) {
  read_tsv(expr_file, show_col_types=FALSE) %>%
    dplyr::select(gene_id, TPM) %>%
    dplyr::mutate(gene_id=sub("\\..*", "", gene_id)) %>%
    dplyr::distinct(gene_id, .keep_all=TRUE)
}

# Wildtype Hi-C eigenvector (500kb bedgraph); sign set so that bins
# overlapping highly expressed genes are positive (A), as in saddle_analysis.R
load_eigen_wt <- function(comp_bed, high_genes) {
  stopifnot(file.exists(comp_bed))
  df <- read.table(comp_bed, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names=c("chr", "start", "end", "eigen")) %>%
    dplyr::mutate(chr=paste0("chr", sub("^chr", "", as.character(chr)))) %>%
    dplyr::filter(!is.na(eigen))
  bin_gr  <- GRanges(df$chr, IRanges(df$start + 1L, df$end))
  gene_gr <- GRanges(high_genes$chr, IRanges(high_genes$start, high_genes$end))
  hits      <- findOverlaps(bin_gr, gene_gr)
  high_bins <- unique(queryHits(hits))
  mean_high <- mean(df$eigen[high_bins], na.rm=TRUE)
  message("  Eigenvector: ", nrow(df), " bins; mean eigen (High-expr bins) = ",
          round(mean_high, 4), "  (",
          round(mean(df$eigen[high_bins] > 0, na.rm=TRUE) * 100, 1), "% positive)")
  if (mean_high < 0) {
    message("  Sign flipped")
    df$eigen <- -df$eigen
  }
  df
}

# "A" / "B" when every overlapped eigenvector bin has that sign, else NA
assign_compartment <- function(regions, eigen_df) {
  bin_gr <- GRanges(eigen_df$chr, IRanges(eigen_df$start + 1L, eigen_df$end))
  reg_gr <- GRanges(regions$chr, IRanges(regions$start, regions$end))
  hits   <- findOverlaps(reg_gr, bin_gr)
  sgn    <- sign(eigen_df$eigen[subjectHits(hits)])
  comp   <- rep(NA_character_, nrow(regions))
  all_pos <- tapply(sgn > 0, queryHits(hits), all)
  all_neg <- tapply(sgn < 0, queryHits(hits), all)
  comp[as.integer(names(all_pos))[all_pos]] <- "A"
  comp[as.integer(names(all_neg))[all_neg]] <- "B"
  comp
}

# Neighborhood expression: for every gene, the mean TPM of the gene and of all
# genes within bp of it (gap <= bp, overlapping genes included; genes without
# a TPM count as 0), and the number of genes with TPM > tpm_max in that set
add_neighborhood <- function(genes, bp, tpm_max) {
  tpm0 <- ifelse(is.na(genes$TPM), 0, genes$TPM)
  gr   <- GRanges(genes$chr, IRanges(genes$start, genes$end))
  hits <- findOverlaps(gr, gr, maxgap=bp)
  agg  <- data.table::data.table(q=queryHits(hits), t=tpm0[subjectHits(hits)])[
    , .(s=sum(t), n=.N, n_expr=sum(t > tpm_max)), by=q]
  genes$TPM_nbhd      <- NA_real_
  genes$n_expr_nbhd   <- NA_integer_
  genes$TPM_nbhd[agg$q]    <- agg$s / agg$n
  genes$n_expr_nbhd[agg$q] <- agg$n_expr
  stopifnot(!anyNA(genes$TPM_nbhd), !anyNA(genes$n_expr_nbhd))
  genes
}

# Transcriptional activity flanking each region: for every distance band on both
# sides of the region (measured from its edges), the summed TPM of genes whose
# midpoint falls in the band, per Mb of band (genes without a TPM count as 0).
# reg: region_id, chr, start, end, group
neighborhood_activity <- function(reg, all_genes, chrom_len, edges=NBHD_BANDS_BP) {
  tpm0 <- ifelse(is.na(all_genes$TPM), 0, all_genes$TPM)
  mid  <- GRanges(all_genes$chr,
                  IRanges(as.integer(floor((all_genes$start + all_genes$end) / 2)), width=1L))
  cl   <- chrom_len[reg$chr]
  out  <- vector("list", length(edges) - 1L)
  for (b in seq_len(length(edges) - 1L)) {
    lo <- edges[b]
    hi <- edges[b + 1L]
    sides <- list(
      list(s=pmax(reg$start - hi, 1),   e=reg$start - lo - 1),
      list(s=reg$end + lo + 1,          e=pmin(reg$end + hi, cl)))
    tpm_sum <- numeric(nrow(reg))
    bp      <- numeric(nrow(reg))
    for (sd in sides) {
      ok  <- which(sd$e >= sd$s)
      gr  <- GRanges(reg$chr[ok], IRanges(as.integer(sd$s[ok]), as.integer(sd$e[ok])))
      h   <- suppressWarnings(findOverlaps(gr, mid))
      add <- tapply(tpm0[subjectHits(h)], queryHits(h), sum)
      ix  <- ok[as.integer(names(add))]
      tpm_sum[ix] <- tpm_sum[ix] + as.numeric(add)
      bp[ok]      <- bp[ok] + width(gr)
    }
    out[[b]] <- data.frame(region_id=reg$region_id, group=reg$group,
                           band=paste0(ifelse(lo == 0, "0", fmt_bp(lo)), "-", fmt_bp(hi)),
                           tpm_per_mb=ifelse(bp > 0, tpm_sum / bp * 1e6, NA_real_),
                           stringsAsFactors=FALSE)
  }
  dplyr::bind_rows(out)
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Control regions
# ══════════════════════════════════════════════════════════════════════════════
choose_n_bins <- function(total_bp) {
  as.integer(min(MAX_BINS, max(MIN_BINS, floor(total_bp / POOL_MIN_BP))))
}

# Take, per length bin (quantiles of the expressed lengths), as many control
# genes as there are expressed genes (or all available if fewer)
match_by_length <- function(expr_sizes, ctrl_df, n_bins, seed) {
  set.seed(seed)
  br <- unique(quantile(expr_sizes, probs=seq(0, 1, length.out=n_bins + 1), names=FALSE))
  br[1] <- -Inf
  br[length(br)] <- Inf
  e_bin <- cut(expr_sizes, br, labels=FALSE)
  c_bin <- cut(ctrl_df$size, br, labels=FALSE)
  keep <- integer(0)
  for (b in seq_len(length(br) - 1)) {
    n_e <- sum(e_bin == b)
    idx <- which(c_bin == b)
    n_take <- min(n_e, length(idx))
    if (n_take > 0) keep <- c(keep, idx[sample.int(length(idx), n_take)])
  }
  ctrl_df[sort(keep), , drop=FALSE]
}

# Non-genic space in unexpressed neighborhoods: chromosome minus every annotated
# gene +/- flank, minus every gene with TPM > CTRL_TPM_MAX +/- nbhd_bp (so no
# window comes within nbhd_bp of such a gene); optionally intersected with one
# compartment's bins
build_eligible_space <- function(all_genes, chrom_sizes, flank_bp, nbhd_bp, eigen_df, compartment) {
  chr_gr <- GRanges(chrom_sizes$chr, IRanges(1L, as.integer(chrom_sizes$length)))
  genic  <- GenomicRanges::reduce(GRanges(
    all_genes$chr,
    IRanges(pmax(as.integer(all_genes$start - flank_bp), 1L),
            as.integer(all_genes$end + flank_bp))))
  ex  <- all_genes[!is.na(all_genes$TPM) & all_genes$TPM > CTRL_TPM_MAX, , drop=FALSE]
  act <- GenomicRanges::reduce(GRanges(
    ex$chr,
    IRanges(pmax(as.integer(ex$start - nbhd_bp - 1L), 1L),
            as.integer(ex$end + nbhd_bp + 1L))))
  elig <- suppressWarnings(GenomicRanges::setdiff(GenomicRanges::setdiff(chr_gr, genic), act))
  if (compartment != "all") {
    keep <- if (compartment == "A") eigen_df$eigen > 0 else eigen_df$eigen < 0
    comp_gr <- GenomicRanges::reduce(GRanges(
      eigen_df$chr[keep], IRanges(eigen_df$start[keep] + 1L, eigen_df$end[keep])))
    elig <- suppressWarnings(GenomicRanges::intersect(elig, comp_gr))
  }
  data.frame(chr=as.character(seqnames(elig)), start=start(elig),
             end=end(elig), width=width(elig), stringsAsFactors=FALSE)
}

# Place one window per requested length at a random position inside the
# eligible space; windows may not overlap each other
sample_windows <- function(lengths, eligible, seed, max_tries=50L) {
  set.seed(seed)
  lengths <- sort(lengths, decreasing=TRUE)
  chrs <- unique(eligible$chr)
  acc_s <- setNames(rep(list(numeric(0)), length(chrs)), chrs)
  acc_e <- setNames(rep(list(numeric(0)), length(chrs)), chrs)
  out <- vector("list", length(lengths))
  n_ok <- 0L
  for (L in lengths) {
    cand <- which(eligible$width >= L + 1)
    if (length(cand) == 0L) next
    w <- eligible$width[cand] - L
    for (tr in seq_len(max_tries)) {
      i  <- cand[sample.int(length(cand), 1L, prob=w)]
      s  <- eligible$start[i] + sample.int(eligible$width[i] - L, 1L) - 1L
      e  <- s + L
      ch <- eligible$chr[i]
      if (!any(acc_s[[ch]] <= e & acc_e[[ch]] >= s)) {
        acc_s[[ch]] <- c(acc_s[[ch]], s)
        acc_e[[ch]] <- c(acc_e[[ch]], e)
        n_ok <- n_ok + 1L
        out[[n_ok]] <- data.frame(chr=ch, start=s, end=e, size=L, stringsAsFactors=FALSE)
        break
      }
    }
  }
  message("  Windows placed: ", n_ok, " of ", length(lengths))
  if (n_ok == 0L) stop("No windows could be placed in the eligible space")
  w_df <- dplyr::bind_rows(out[seq_len(n_ok)])
  w_df$region_id <- sprintf("win_%05d", seq_len(nrow(w_df)))
  w_df
}

# Random pools of >= pool_bp total bp; the last incomplete pool is dropped
assign_pools_by_bp <- function(regions, pool_bp, seed) {
  set.seed(seed)
  regions <- regions[sample.int(nrow(regions)), , drop=FALSE]
  pool <- integer(nrow(regions))
  cur  <- 1L
  acc  <- 0
  for (i in seq_len(nrow(regions))) {
    pool[i] <- cur
    acc <- acc + regions$size[i]
    if (acc >= pool_bp) { cur <- cur + 1L; acc <- 0 }
  }
  if (acc > 0) pool[pool == cur] <- NA_integer_
  regions$pool <- pool
  regions[!is.na(regions$pool), , drop=FALSE]
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Contact extraction (aggregated by distance per region)
#    profile: region_id, index (|x - y|), counts   [index <= profile_max]
#    meta:    region_id, max_index, n_x_bins, n_x_present
#    Extraction units are (condition x chromosome) for "chrom" and
#    (condition x chromosome x tile) for "tiles"; all units of all conditions
#    run in one pool of N_CORES workers.
# ══════════════════════════════════════════════════════════════════════════════
# canonical "chrN" name -> chromosome name used inside the .hic file
hic_chrom_map <- function(hic_file) {
  nm  <- strawr::readHicChroms(hic_file)$name
  out <- ifelse(STD_CHROMS %in% nm, STD_CHROMS, sub("^chr", "", STD_CHROMS))
  if (!all(out %in% nm))
    stop("Chromosomes missing from ", basename(hic_file), ": ",
         paste(STD_CHROMS[!out %in% nm], collapse=", "))
  setNames(out, STD_CHROMS)
}

# per-region aggregation of a pixel table (x <= y, one row per pixel)
slice_regions <- function(raw_dt, regs, cfg) {
  raw_dt <- raw_dt[!is.nan(counts) & !is.na(counts)]
  raw_dt[, index := abs(x - y)]
  if (cfg$band_only) raw_dt <- raw_dt[index <= cfg$profile_max]
  data.table::setkey(raw_dt, x)
  xs <- raw_dt$x
  res_list <- lapply(seq_len(nrow(regs)), function(i) {
    g  <- regs[i, ]
    i1 <- findInterval(g$start - 1e-6, xs) + 1L
    i2 <- findInterval(g$end, xs)
    if (i1 > i2) return(NULL)
    sub <- raw_dt[i1:i2]
    list(profile=sub[index <= cfg$profile_max, .(counts=sum(counts)), by=index][, region_id := g$region_id],
         meta=data.table::data.table(region_id=g$region_id, max_index=max(sub$index),
                                     n_x_present=data.table::uniqueN(sub$x)))
  })
  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0) return(NULL)
  list(profile=data.table::rbindlist(lapply(res_list, `[[`, "profile")),
       meta=data.table::rbindlist(lapply(res_list, `[[`, "meta")))
}

# one straw() call for a whole chromosome (every block decoded once)
extract_chrom <- function(hic_file, hic_chrom, regs, cfg, norm) {
  t0  <- Sys.time()
  raw <- strawr::straw(norm, hic_file, hic_chrom, hic_chrom, "BP", cfg$res)
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  n_raw <- nrow(raw)
  out   <- slice_regions(data.table::setDT(raw), regs, cfg)
  message("    ", hic_chrom, ": ", format(n_raw, big.mark=","), " pixels read, ",
          nrow(regs), " regions, ",
          round(as.numeric(difftime(Sys.time(), t0, units="secs"))), " s")
  out
}

# one straw() call for a tile: x in [first start, last end], y up to profile_max beyond
extract_tile <- function(hic_file, hic_chrom, regs, cfg, norm, chrom_len) {
  qs  <- min(regs$start)
  qe  <- max(regs$end)
  far <- min(qe + cfg$profile_max, chrom_len)
  raw <- strawr::straw(norm, hic_file,
                       sprintf("%s:%.0f:%.0f", hic_chrom, qs, qe),
                       sprintf("%s:%.0f:%.0f", hic_chrom, qs, far),
                       "BP", cfg$res)
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  slice_regions(data.table::setDT(raw), regs, cfg)
}

# regions of one chromosome -> tiles of consecutive regions spanning <= tile_bp
split_into_tiles <- function(regs, tile_bp) {
  regs <- regs[order(regs$start), , drop=FALSE]
  tile <- integer(nrow(regs))
  cur  <- 1L
  t0   <- regs$start[1]
  for (i in seq_len(nrow(regs))) {
    if (regs$start[i] - t0 > tile_bp) { cur <- cur + 1L; t0 <- regs$start[i] }
    tile[i] <- cur
  }
  split(regs, tile)
}

extract_profiles <- function(cond_keys, regions, cfg, chrom_len, norm=STRAW_NORM,
                             n_cores=N_CORES) {
  chroms <- intersect(unique(regions$chr), names(chrom_len))
  chroms <- chroms[order(-chrom_len[chroms])]   # largest first for load balance
  hmap   <- setNames(lapply(cond_keys, function(n) hic_chrom_map(CONDITIONS[[n]]$hic)), cond_keys)
  
  units <- list()
  for (chrom in chroms) {
    regs  <- regions[regions$chr == chrom, , drop=FALSE]
    parts <- if (cfg$extract == "tiles") split_into_tiles(regs, cfg$tile_bp) else list(regs)
    for (p in parts) for (n in cond_keys)
      units[[length(units) + 1L]] <- list(cond=n, chrom=chrom, regs=p)
  }
  n_u <- length(units)
  message("  Extracting ", nrow(regions), " regions: ", n_u, " ", cfg$extract, " units (",
          length(cond_keys), " conditions, ", length(chroms), " chromosomes), res=", cfg$res,
          ", cores=", n_cores)
  
  run_unit <- function(u) {
    hic <- CONDITIONS[[u$cond]]$hic
    hc  <- hmap[[u$cond]][[u$chrom]]
    tryCatch(
      if (cfg$extract == "tiles") extract_tile(hic, hc, u$regs, cfg, norm, chrom_len[[u$chrom]])
      else extract_chrom(hic, hc, u$regs, cfg, norm),
      error=function(e) list(error=paste0(u$cond, " ", u$chrom, ": ", conditionMessage(e))))
  }
  
  out <- vector("list", n_u)
  t0  <- Sys.time()
  done <- 0L
  if (cfg$extract == "tiles") {
    # time the first units in this process and project the total
    k  <- min(3L, n_u)
    tp <- system.time(out[seq_len(k)] <- lapply(units[seq_len(k)], run_unit))[["elapsed"]]
    proj <- tp / k * n_u / n_cores / 60
    message(sprintf("    probe: %d tiles in %.1f s  ->  projected %.1f min on %d cores%s",
                    k, tp, proj, n_cores,
                    if (proj > 30) "  (> 30 min: request more cores or raise tile_bp)" else ""))
    done <- k
  }
  chunk <- if (cfg$extract == "tiles") max(n_cores * 20L, 40L) else n_u
  while (done < n_u) {
    idx <- (done + 1L):min(done + chunk, n_u)
    out[idx] <- parallel::mclapply(units[idx], run_unit,
                                   mc.cores=if (cfg$extract == "tiles") n_cores else min(n_cores, 8L),
                                   mc.preschedule=(cfg$extract == "tiles"))
    done <- max(idx)
    el <- as.numeric(difftime(Sys.time(), t0, units="mins"))
    message(sprintf("    %d / %d units  |  %.1f min elapsed  |  ~%.1f min left",
                    done, n_u, el, el * (n_u - done) / done))
  }
  
  errs <- unlist(lapply(out, function(r)
    if (inherits(r, "try-error")) as.character(r) else if (is.list(r) && !is.null(r$error)) r$error))
  if (length(errs) > 0)
    stop(length(errs), " extraction units failed; first: ", errs[1])
  
  res <- lapply(setNames(cond_keys, cond_keys), function(n) {
    keep <- Filter(Negate(is.null), out[vapply(units, function(u) u$cond == n, logical(1))])
    list(profile=data.table::rbindlist(lapply(keep, `[[`, "profile")),
         meta=data.table::rbindlist(lapply(keep, `[[`, "meta")) %>%
           dplyr::left_join(regions %>% dplyr::select(region_id, n_x_bins), by="region_id") %>%
           data.table::as.data.table())
  })
  message("    done in ", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), " min")
  res
}

# cached per condition; a cache is used only if regions and extraction settings match
get_profiles_all <- function(regions, cfg, chrom_len, stratum_i, rds_dir) {
  sig <- list(regions=regions[, c("region_id", "chr", "start", "end")],
              res=cfg$res, profile_max=cfg$profile_max, band_only=cfg$band_only)
  cache_of <- function(n) file.path(rds_dir, sprintf("profiles_%s_s%d.rds", n, stratum_i))
  out  <- list()
  todo <- character(0)
  for (n in names(CONDITIONS)) {
    if (file.exists(cache_of(n))) {
      obj <- readRDS(cache_of(n))
      if (identical(obj$sig, sig)) {
        message("  ", CONDITIONS[[n]]$display, ": loaded cached profiles")
        out[[n]] <- obj[c("profile", "meta")]
        next
      }
      message("  ", CONDITIONS[[n]]$display, ": cache does not match regions, re-extracting")
    }
    todo <- c(todo, n)
  }
  if (length(todo) > 0) {
    ex <- extract_profiles(todo, regions, cfg, chrom_len)
    for (n in todo) {
      out[[n]] <- ex[[n]]
      saveRDS(c(list(sig=sig), ex[[n]]), cache_of(n))
    }
  }
  out[names(CONDITIONS)]
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Pooled P(s) and S
#    Pools are summed over regions per distance. "original" mode keeps the
#    possible-contacts term of scaling_analysis.R; "constant" mode uses the
#    pooled number of x bins. Constant factors (col_sum) shift log10 P(s)
#    only and do not change the slope.
# ══════════════════════════════════════════════════════════════════════════════
pooled_dist <- function(prof, meta, map, keys, cfg) {
  res <- cfg$res
  map <- data.table::as.data.table(map)
  d <- merge(prof, map, by="region_id", allow.cartesian=TRUE)[
    , .(counts=sum(counts)), by=c(keys, "index")]
  m <- merge(meta, map, by="region_id", allow.cartesian=TRUE)[
    , .(total_possible=mean(floor(max_index / res)),
        n_bins=sum(n_x_bins), n_regions=.N), by=keys]
  d <- merge(d, m, by=keys)
  if (cfg$possible_mode == "original") {
    d[, possible := total_possible * res - index]
    d <- d[possible >= res]
  } else {
    d[, possible := as.numeric(n_bins)]
  }
  d[, contact.prob := counts / possible]
  d
}

fit_pool <- function(dd, cfg) {
  res <- cfg$res
  win <- dd[index >= cfg$fit_min & index < cfg$fit_max]
  n_expected <- length(seq(cfg$fit_min, cfg$fit_max - 1, by=res))
  rng <- win[counts > 0, .(distance=log10(index), prob=log10(contact.prob))]
  rng <- rng[is.finite(prob)]
  if (nrow(rng) < 3) return(NULL)
  model <- lm(prob ~ distance, data=rng)
  data.frame(kind=dd$kind[1], pool=dd$pool[1],
             slope=unname(coef(model)[2]), intercept=unname(coef(model)[1]),
             r2=summary(model)$r.squared, n_points=nrow(rng),
             frac_present=nrow(rng) / n_expected,
             med_counts=median(win$counts), stringsAsFactors=FALSE)
}

compute_pool_slopes <- function(prof, meta, pool_map, cfg) {
  d <- pooled_dist(prof, meta, pool_map, c("kind", "pool"), cfg)
  lst <- split(d, by=c("kind", "pool"), drop=TRUE)
  dplyr::bind_rows(lapply(lst, fit_pool, cfg=cfg))
}

compute_ps_groups <- function(prof, meta, group_map, cfg) {
  d <- pooled_dist(prof, meta, group_map, "group", cfg)
  d <- d[counts > 0 & index > 0 & index <= cfg$profile_max]
  as.data.frame(d[, .(group, index, distance=log10(index), prob=log10(contact.prob))])
}

compute_loess_deriv <- function(ps_df, lo_bp, hi_bp, span=0.3) {
  ps_df <- ps_df %>%
    dplyr::filter(distance >= log10(lo_bp), distance <= log10(hi_bp)) %>%
    dplyr::arrange(distance)
  if (nrow(ps_df) < 10) return(NULL)
  fit <- tryCatch(loess(prob ~ distance, data=ps_df, span=span),
                  error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  x_seq  <- seq(min(ps_df$distance), max(ps_df$distance), length.out=200)
  y_pred <- predict(fit, newdata=data.frame(distance=x_seq))
  dy <- diff(y_pred) / diff(x_seq)
  dx <- (x_seq[-1] + x_seq[-length(x_seq)]) / 2
  data.frame(distance=dx, dist_bp=10^dx,
             prob_smooth=predict(fit, newdata=data.frame(distance=dx)),
             slope=dy)
}

# Region-level bootstrap of delta-S(d): within every group, regions are resampled
# with replacement (the same resample for both conditions), pooled P(s), the loess
# derivative and delta-S(d) are recomputed, and every resampled curve is interpolated
# onto the observed distance grid of its group. Returns 2.5% / 97.5% quantiles.
bootstrap_delta_loess <- function(ex, group_map, obs_delta, cfg, n_boot=N_BOOT, n_cores=N_CORES) {
  keep <- unique(group_map$region_id)
  ex_k <- lapply(ex, function(e) list(profile=e$profile[region_id %in% keep],
                                      meta=e$meta[region_id %in% keep]))
  gm <- data.table::as.data.table(group_map)
  obs_grp <- as.character(unique(obs_delta$group))
  one <- function(b) {
    set.seed(SEED + b)
    bm <- gm[, .SD[sample.int(.N, replace=TRUE)], by=group]
    sl <- lapply(ex_k, function(e) {
      compute_ps_groups(e$profile, e$meta, bm, cfg) %>%
        dplyr::group_by(group) %>%
        dplyr::group_modify(~compute_loess_deriv(.x, cfg$loess_lo, cfg$loess_hi)) %>%
        dplyr::ungroup()
    })
    dplyr::bind_rows(lapply(obs_grp, function(g) {
      og <- obs_delta$dist_bp[as.character(obs_delta$group) == g]
      s1 <- sl[[1]][sl[[1]]$group == g, ]
      s2 <- sl[[2]][sl[[2]]$group == g, ]
      if (nrow(s1) < 2 || nrow(s2) < 2) return(NULL)
      data.frame(group=g, dist_bp=og,
                 delta=approx(log10(s1$dist_bp), s1$slope, log10(og))$y -
                   approx(log10(s2$dist_bp), s2$slope, log10(og))$y,
                 stringsAsFactors=FALSE)
    }))
  }
  res <- parallel::mclapply(seq_len(n_boot), function(b) try(one(b), silent=TRUE), mc.cores=n_cores)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (all(bad)) stop("Bootstrap failed in every resample; first error: ", as.character(res[[1]]))
  message("    bootstrap: ", sum(!bad), " of ", n_boot, " resamples completed")
  dplyr::bind_rows(res[!bad]) %>%
    dplyr::group_by(group, dist_bp) %>%
    dplyr::summarise(lo=quantile(delta, 0.025, na.rm=TRUE), hi=quantile(delta, 0.975, na.rm=TRUE),
                     .groups="drop") %>%
    dplyr::mutate(group=factor(group, levels=SD_LEVELS))
}

label_groups <- function(df) {
  df %>% dplyr::mutate(
    group = dplyr::case_when(
      kind == "EI"  & tert == 1 ~ "Expressed, Low E/I",
      kind == "EI"  & tert == 3 ~ "Expressed, High E/I",
      kind == "TPM" & tert == 1 ~ "Expressed, Low nbhd TPM",
      kind == "TPM" & tert == 3 ~ "Expressed, High nbhd TPM",
      kind == "NONEXPR"  ~ "Non-expressed genes",
      kind == "NONGENIC" ~ "Non-genic windows",
      TRUE ~ NA_character_))
}

compare_delta <- function(delta_df) {
  ctrls <- c("Non-expressed genes", "Non-genic windows")
  refs  <- list("Expressed (all E/I bins)" = delta_df$delta[delta_df$kind == "EI"],
                "Expressed (all TPM bins)" = delta_df$delta[delta_df$kind == "TPM"])
  dplyr::bind_rows(lapply(ctrls, function(ctl) {
    b <- delta_df$delta[delta_df$group %in% ctl]
    dplyr::bind_rows(lapply(names(refs), function(rn) {
      a <- refs[[rn]]
      data.frame(control=ctl, reference=rn, n_ref=length(a), n_ctl=length(b),
                 median_deltaS_ref=median(a), median_deltaS_ctl=median(b),
                 diff_of_medians=median(a) - median(b),
                 p_wilcox=tryCatch(wilcox.test(a, b)$p.value, error=function(e) NA_real_),
                 stringsAsFactors=FALSE)
    }))
  }))
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. One length stratum: regions -> contacts -> pools -> S, delta-S, S(d)
# ══════════════════════════════════════════════════════════════════════════════
# E/I bins and neighborhood-TPM bins of a set of expressed genes
make_expressed_pools <- function(expressed) {
  ei_set  <- expressed %>% dplyr::filter(intron_bp > 0)
  n_ei    <- choose_n_bins(sum(ei_set$size))
  n_tpm   <- choose_n_bins(sum(expressed$size))
  ei_set$pool  <- dplyr::ntile(ei_set$EI, n_ei)
  tpm_set <- expressed
  tpm_set$pool <- dplyr::ntile(tpm_set$TPM_nbhd, n_tpm)
  bin_bp  <- c(tapply(ei_set$size, ei_set$pool, sum), tapply(tpm_set$size, tpm_set$pool, sum))
  list(ei_set=ei_set, tpm_set=tpm_set, n_ei=n_ei, n_tpm=n_tpm, pool_bp=median(bin_bp))
}

run_stratum <- function(genes_s, all_genes, chrom_sizes, eigen_df,
                        cfg, stratum_i, stratum_label, rds_dir) {
  res <- cfg$res
  chrom_len <- setNames(chrom_sizes$length, chrom_sizes$chr)
  message("\n══ ", ASSAY, "  |  stratum ", stratum_label, "  |  compartment: ", COMPARTMENT,
          "  |  res: ", res, " ══")
  
  # -- genes in stratum
  genes_s <- genes_s %>%
    dplyr::mutate(n_x_bins = floor(end / res) - ceiling(start / res) + 1)
  message("  genes in stratum: ", nrow(genes_s),
          "  (dropped, no bin start inside gene: ", sum(genes_s$n_x_bins < 1), ")")
  genes_s <- genes_s %>% dplyr::filter(n_x_bins >= 1)
  if (COMPARTMENT != "all") {
    genes_s$comp <- assign_compartment(genes_s, eigen_df)
    message("  compartment calls: A=", sum(genes_s$comp == "A", na.rm=TRUE),
            " B=", sum(genes_s$comp == "B", na.rm=TRUE),
            " mixed/NA=", sum(is.na(genes_s$comp)))
  }
  
  # -- shared control set: genome-wide (all compartments), identical in every run.
  #    Candidates: TPM <= NONEXPR_TPM_MAX with no gene of TPM > CTRL_TPM_MAX within CTRL_NBHD_BP.
  expressed_all <- genes_s %>% dplyr::filter(!is.na(TPM), TPM > 0)
  zero_tpm  <- genes_s %>% dplyr::filter(!is.na(TPM), TPM <= NONEXPR_TPM_MAX)
  ctrl_cand <- zero_tpm %>% dplyr::filter(n_expr_ctrl == 0)
  message("  expressed (all compartments): ", nrow(expressed_all), "  TPM<=", NONEXPR_TPM_MAX, ": ",
          nrow(zero_tpm), "  of which with no gene of TPM>", CTRL_TPM_MAX, " within ",
          CTRL_NBHD_BP / 1000, " kb: ", nrow(ctrl_cand))
  if (nrow(expressed_all) < MIN_BINS) stop("Too few expressed genes in this stratum")
  if (nrow(ctrl_cand) == 0) stop("No non-expressed control genes in unexpressed neighborhoods")
  ctrl_gr <- with(all_genes[!is.na(all_genes$TPM) & all_genes$TPM > CTRL_TPM_MAX, ],
                  GRanges(chr, IRanges(start, end)))
  
  pool_bp_all <- make_expressed_pools(expressed_all)$pool_bp
  ctrl_genes <- match_by_length(expressed_all$size, ctrl_cand, LEN_MATCH_BINS, SEED)
  message("  control genes length-matched: ", nrow(ctrl_genes), " (expressed: ", nrow(expressed_all), ")")
  elig <- build_eligible_space(all_genes, chrom_sizes, NONGENIC_FLANK_BP, CTRL_NBHD_BP, NULL, "all")
  message("  eligible non-genic space (>", CTRL_NBHD_BP / 1000, " kb from any gene of TPM>", CTRL_TPM_MAX,
          "): ", format(sum(elig$width), big.mark=","), " bp")
  set.seed(SEED)
  # request enough windows for CTRL_POOLS pools, but never more than ~40% of the eligible space
  # (a nearly full space places long windows first and leaves out short ones)
  n_win <- ceiling(CTRL_POOLS * pool_bp_all * 1.25 / mean(expressed_all$size))
  n_cap <- floor(0.4 * sum(elig$width) / mean(expressed_all$size))
  if (n_win > n_cap) message("  windows requested capped by eligible space: ", n_win, " -> ", n_cap)
  n_win <- min(n_win, n_cap)
  win_lens <- sample(expressed_all$size, n_win, replace=TRUE)
  wins <- sample_windows(win_lens, elig, SEED)
  wins$n_x_bins <- floor(wins$end / res) - ceiling(wins$start / res) + 1
  wins <- wins %>% dplyr::filter(n_x_bins >= 1)
  n_placed <- nrow(wins)
  wins <- match_by_length(expressed_all$size, wins, LEN_MATCH_BINS, SEED)
  message("  windows after length matching: ", nrow(wins), " (placed: ", n_placed, ")")
  
  # -- CHECK: controls really are far from every gene of TPM > CTRL_TPM_MAX
  n_bad_g <- sum(suppressWarnings(overlapsAny(GRanges(ctrl_genes$chr, IRanges(ctrl_genes$start, ctrl_genes$end)),
                                              ctrl_gr, maxgap=CTRL_NBHD_BP)))
  n_bad_w <- sum(suppressWarnings(overlapsAny(GRanges(wins$chr, IRanges(wins$start, wins$end)),
                                              ctrl_gr, maxgap=CTRL_NBHD_BP)))
  message("  CHECK controls with a gene of TPM>", CTRL_TPM_MAX, " within ", CTRL_NBHD_BP / 1000,
          " kb: genes=", n_bad_g, " windows=", n_bad_w, "  (must be 0)")
  stopifnot(n_bad_g == 0, n_bad_w == 0)
  
  # -- expressed genes of this run (one compartment, or all) and their pools
  expressed <- if (COMPARTMENT == "all") expressed_all else
    expressed_all %>% dplyr::filter(comp == COMPARTMENT)
  message("  expressed in this run (compartment ", COMPARTMENT, "): ", nrow(expressed))
  if (nrow(expressed) < MIN_BINS) stop("Too few expressed genes in this stratum/compartment")
  pl      <- make_expressed_pools(expressed)
  ei_set  <- pl$ei_set
  tpm_set <- pl$tpm_set
  n_ei    <- pl$n_ei
  n_tpm   <- pl$n_tpm
  pool_bp <- pl$pool_bp
  message("  E/I bins: ", n_ei, "  neighborhood-TPM bins: ", n_tpm, "  median bp per bin: ",
          format(pool_bp, big.mark=","))
  
  # -- control supply vs pool size. All pools hold the same bp; when the controls cannot
  #    fill MIN_POOLS pools of the expressed pool size (with 25% margin for the coverage
  #    filter), every group uses smaller pools.
  supply <- min(sum(ctrl_genes$size), sum(wins$size))
  message("  control material (bp): genes=", format(sum(ctrl_genes$size), big.mark=","),
          "  windows=", format(sum(wins$size), big.mark=","),
          "  | pool size (bp): ", format(round(pool_bp), big.mark=","),
          "  | needed for ", MIN_POOLS, " pools + 25% margin: ",
          format(round(MIN_POOLS * 1.25 * pool_bp), big.mark=","))
  if (supply < MIN_POOLS * 1.25 * pool_bp) {
    new_bp <- floor(supply / (MIN_POOLS * 1.25))
    if (new_bp < POOL_MIN_BP)
      warning("Control sequence (", format(supply, big.mark=","), " bp in the smaller control type) ",
              "only supports pools of ", format(new_bp, big.mark=","), " bp, below POOL_MIN_BP = ",
              format(POOL_MIN_BP, big.mark=","), ": control and expressed pools are small and noisy",
              call.=FALSE)
    n_ei  <- as.integer(floor(sum(ei_set$size)    / new_bp))
    n_tpm <- as.integer(floor(sum(expressed$size) / new_bp))
    ei_set$pool  <- dplyr::ntile(ei_set$EI, n_ei)
    tpm_set$pool <- dplyr::ntile(tpm_set$TPM_nbhd, n_tpm)
    bin_bp  <- c(tapply(ei_set$size, ei_set$pool, sum), tapply(tpm_set$size, tpm_set$pool, sum))
    message("  pool size reduced to fit the controls: ", format(round(pool_bp), big.mark=","), " -> ",
            format(round(median(bin_bp)), big.mark=","), " bp  (E/I bins: ", n_ei,
            ", neighborhood-TPM bins: ", n_tpm, ")")
    pool_bp <- median(bin_bp)
  }
  
  # -- regions to extract: all expressed genes + the shared controls (so every run of this
  #    stratum, A, B or all, reads the same cached contacts)
  regions <- dplyr::bind_rows(
    expressed_all %>% dplyr::transmute(region_id=gene_id, chr, start, end, size, n_x_bins, type="expressed"),
    ctrl_genes    %>% dplyr::transmute(region_id=gene_id, chr, start, end, size, n_x_bins, type="nonexpr_gene"),
    wins          %>% dplyr::transmute(region_id, chr, start, end, size, n_x_bins, type="window")
  ) %>% as.data.frame()
  stopifnot(!anyDuplicated(regions$region_id))
  
  # -- contacts for both conditions
  ex <- get_profiles_all(regions, cfg, chrom_len, stratum_i, rds_dir)
  
  # -- coverage filter (region must pass in both conditions)
  pass <- Reduce(intersect, lapply(ex, function(e)
    e$meta[n_x_present / n_x_bins >= MIN_COVERAGE]$region_id))
  used <- regions[regions$type != "expressed" | regions$region_id %in% expressed$gene_id, ]
  message("  coverage >= ", MIN_COVERAGE, " in both conditions: ",
          paste(names(table(used$type[used$region_id %in% pass])),
                table(used$type[used$region_id %in% pass]), sep="=", collapse="  "),
          "  (of ", paste(names(table(used$type)), table(used$type), sep="=", collapse="  "), ")")
  
  ctrl_genes <- ctrl_genes %>% dplyr::filter(gene_id %in% pass) %>%
    dplyr::mutate(region_id=gene_id) %>% assign_pools_by_bp(pool_bp, SEED)
  wins <- wins %>% dplyr::filter(region_id %in% pass) %>% assign_pools_by_bp(pool_bp, SEED) %>%
    dplyr::filter(pool <= CTRL_POOLS)
  ei_set  <- ei_set  %>% dplyr::filter(gene_id %in% pass)
  tpm_set <- tpm_set %>% dplyr::filter(gene_id %in% pass)
  
  pool_map <- dplyr::bind_rows(
    ei_set     %>% dplyr::transmute(region_id=gene_id, kind="EI",  pool),
    tpm_set    %>% dplyr::transmute(region_id=gene_id, kind="TPM", pool),
    ctrl_genes %>% dplyr::transmute(region_id, kind="NONEXPR",  pool),
    wins       %>% dplyr::transmute(region_id, kind="NONGENIC", pool))
  pool_x <- dplyr::bind_rows(
    ei_set  %>% dplyr::group_by(pool) %>% dplyr::summarise(x=median(EI),       .groups="drop") %>%
      dplyr::mutate(kind="EI"),
    tpm_set %>% dplyr::group_by(pool) %>% dplyr::summarise(x=median(TPM_nbhd), .groups="drop") %>%
      dplyr::mutate(kind="TPM"))
  pool_n <- pool_map %>% dplyr::distinct(kind, pool) %>% dplyr::count(kind)
  message("  pools per kind: ", paste(pool_n$kind, pool_n$n, sep="=", collapse="  "))
  if (any(pool_n$n < MIN_POOLS) || nrow(pool_n) < 4)
    stop("Fewer than ", MIN_POOLS, " pools in at least one group: ",
         paste(pool_n$kind, pool_n$n, sep="=", collapse="  "))
  
  # -- CHECK: length matching (median size and KS test against expressed)
  size_tab <- dplyr::bind_rows(
    expressed  %>% dplyr::filter(gene_id %in% pass) %>% dplyr::transmute(grp="Expressed", size),
    ctrl_genes %>% dplyr::transmute(grp="Non-expressed genes", size),
    wins       %>% dplyr::transmute(grp="Non-genic windows", size))
  ks_p <- sapply(c("Non-expressed genes", "Non-genic windows"), function(g)
    suppressWarnings(ks.test(size_tab$size[size_tab$grp == g],
                             size_tab$size[size_tab$grp == "Expressed"])$p.value))
  message("  CHECK length matching (median size, bp): ",
          paste(names(tapply(size_tab$size, size_tab$grp, median)),
                round(tapply(size_tab$size, size_tab$grp, median)), sep="=", collapse="  "),
          "  | KS p vs Expressed: ", paste(names(ks_p), signif(ks_p, 2), sep="=", collapse="  "))
  
  # -- pooled S per condition
  slopes <- dplyr::bind_rows(lapply(names(CONDITIONS), function(n) {
    s <- compute_pool_slopes(ex[[n]]$profile, ex[[n]]$meta, pool_map, cfg)
    s$cond <- CONDITIONS[[n]]$display
    s
  })) %>% dplyr::mutate(cond=factor(cond, levels=COND_LEVELS))
  
  pool_info <- dplyr::bind_rows(
    data.frame(kind="EI",  pool=seq_len(n_ei),  tert=dplyr::ntile(seq_len(n_ei), 3)),
    data.frame(kind="TPM", pool=seq_len(n_tpm), tert=dplyr::ntile(seq_len(n_tpm), 3)),
    data.frame(kind="NONEXPR",  pool=seq_len(max(ctrl_genes$pool)), tert=NA_integer_),
    data.frame(kind="NONGENIC", pool=seq_len(max(wins$pool)),       tert=NA_integer_))
  slopes <- slopes %>% dplyr::left_join(pool_info, by=c("kind", "pool")) %>% label_groups()
  
  # -- CHECK: fit quality on live data
  chk <- slopes %>% dplyr::group_by(kind, cond) %>%
    dplyr::summarise(n_pools=dplyr::n(), median_S=round(median(slope), 3),
                     frac_negative=round(mean(slope < 0), 2),
                     median_frac_present=round(median(frac_present), 3),
                     median_counts_per_point=round(median(med_counts)), .groups="drop")
  message("  CHECK pooled S:"); print(as.data.frame(chk))
  if (any(chk$median_frac_present < 0.95))
    warning("Median fraction of non-zero distance points < 0.95 for at least one group: ",
            "pools are sparse at this resolution")
  if (any(chk$frac_negative < 1))
    warning("Some pooled slopes are not negative")
  if (cfg$possible_mode == "original") {
    pd <- pooled_dist(ex$wt$profile, ex$wt$meta, pool_map, c("kind", "pool"), cfg)
    pr <- pd[index %in% c(cfg$fit_min, cfg$fit_max - cfg$res),
             .(ratio=possible[index == cfg$fit_min] / possible[index == cfg$fit_max - cfg$res]),
             by=.(kind, pool)]
    message("  CHECK possible-contacts term, value at fit_min / value at fit_max (should be ~1): ",
            round(median(pr$ratio, na.rm=TRUE), 3))
  }
  
  # -- delta-S per pool (same pools in both conditions)
  delta_df <- slopes %>%
    dplyr::select(kind, pool, group, cond, slope) %>%
    tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
    dplyr::mutate(delta = .data[[COND_LEVELS[1]]] - .data[[COND_LEVELS[2]]]) %>%
    dplyr::filter(!is.na(delta))
  stopifnot(nrow(delta_df) == nrow(slopes) / length(CONDITIONS))
  stat_delta <- compare_delta(delta_df)
  message("  deltaS (", COND_LEVELS[1], " - ", COND_LEVELS[2], "), controls vs expressed pools:")
  print(as.data.frame(stat_delta), digits=3)
  
  # -- S(d): group-level pooled P(s) and loess derivative; Q1 / Q4 of neighborhood TPM
  exp_k <- expressed %>% dplyr::filter(gene_id %in% pass)
  exp_k$q4 <- dplyr::ntile(exp_k$TPM_nbhd, 4L)
  message("  CHECK neighborhood TPM (gene + genes within ", NBHD_BP / 1000, " kb), median: Q1=",
          round(median(exp_k$TPM_nbhd[exp_k$q4 == 1L]), 3), "  Q4=",
          round(median(exp_k$TPM_nbhd[exp_k$q4 == 4L]), 3), "  | controls: 0 by construction")
  group_map <- dplyr::bind_rows(
    exp_k %>% dplyr::filter(q4 %in% c(1L, 4L)) %>%
      dplyr::transmute(region_id=gene_id,
                       group=ifelse(q4 == 1L, "Expressed, Q1 (Low)", "Expressed, Q4 (High)")),
    ctrl_genes %>% dplyr::transmute(region_id, group="Non-expressed genes"),
    wins       %>% dplyr::transmute(region_id, group="Non-genic windows"))
  loess_df <- dplyr::bind_rows(lapply(names(CONDITIONS), function(n) {
    ps <- compute_ps_groups(ex[[n]]$profile, ex[[n]]$meta, group_map, cfg)
    ps %>% dplyr::group_by(group) %>%
      dplyr::group_modify(~compute_loess_deriv(.x, cfg$loess_lo, cfg$loess_hi)) %>%
      dplyr::ungroup() %>% dplyr::mutate(cond=CONDITIONS[[n]]$display)
  })) %>% dplyr::mutate(cond=factor(cond, levels=COND_LEVELS),
                        group=factor(group, levels=SD_LEVELS))
  delta_loess <- loess_df %>%
    dplyr::select(dist_bp, slope, cond, group) %>%
    tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
    dplyr::mutate(delta = .data[[COND_LEVELS[1]]] - .data[[COND_LEVELS[2]]]) %>%
    dplyr::filter(!is.na(delta))
  
  message("  Bootstrapping delta-S(d): ", N_BOOT, " region resamples on ", N_CORES, " cores")
  t_boot   <- Sys.time()
  delta_ci <- bootstrap_delta_loess(ex, group_map, delta_loess, cfg)
  message("    done in ", round(as.numeric(difftime(Sys.time(), t_boot, units="mins")), 1), " min")
  chk_ci <- dplyr::inner_join(delta_loess, delta_ci, by=c("group", "dist_bp")) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(frac_observed_inside_band=round(mean(delta >= lo & delta <= hi, na.rm=TRUE), 2),
                     median_band_width=round(median(hi - lo, na.rm=TRUE), 3), .groups="drop")
  message("  CHECK bootstrap bands (observed delta-S(d) should lie inside its own band):")
  print(as.data.frame(chk_ci))
  
  # -- CHECK: transcriptional activity flanking each group
  reg_g <- regions %>% dplyr::inner_join(group_map, by="region_id") %>%
    dplyr::select(region_id, chr, start, end, group)
  nbhd_act <- neighborhood_activity(reg_g, all_genes, chrom_len)
  nbhd_act$band <- factor(nbhd_act$band, levels=unique(nbhd_act$band))
  nbhd_tab <- nbhd_act %>% dplyr::group_by(group, band) %>%
    dplyr::summarise(n=dplyr::n(), median_TPM_per_Mb=round(median(tpm_per_mb, na.rm=TRUE), 1),
                     mean_TPM_per_Mb=round(mean(tpm_per_mb, na.rm=TRUE), 1), .groups="drop")
  message("  CHECK activity flanking each group (summed TPM of genes per Mb of flank):")
  print(as.data.frame(tidyr::pivot_wider(nbhd_tab[, c("group", "band", "median_TPM_per_Mb")],
                                         names_from=band, values_from=median_TPM_per_Mb)))
  
  list(label=stratum_label, slopes=slopes, delta=delta_df, stats=stat_delta,
       loess=loess_df, delta_loess=delta_loess, pool_map=pool_map, pool_x=pool_x,
       delta_ci=delta_ci, nbhd_act=nbhd_act)
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Plot helpers
# ══════════════════════════════════════════════════════════════════════════════
# delta-S panel in the style of make_delta_panel() in scaling_analysis.R:
# orange line + area across the expressed bins, 5 evenly spaced x labels
# (median E/I or TPM of the bin). Controls: median (dashed) and IQR (band)
# of delta-S across the control pools.
make_delta_pools <- function(delta_df, pool_x, kind_sel, x_title, ctrl_df, ylim) {
  dd <- delta_df %>% dplyr::filter(kind == kind_sel) %>%
    dplyr::left_join(pool_x, by=c("kind", "pool")) %>% dplyr::arrange(pool)
  idx <- unique(round(seq(1, nrow(dd), length.out=5)))
  ctrl_cols <- c("Non-expressed genes"="#7B3294", "Non-genic windows"="#636363")
  ctrl_cols <- ctrl_cols[!names(ctrl_cols) %in% PLOT_DROP_GROUPS]
  ggplot(dd, aes(x=pool, y=delta))+
    geom_rect(data=ctrl_df, aes(xmin=-Inf, xmax=Inf, ymin=q25, ymax=q75, fill=group),
              inherit.aes=FALSE, alpha=0.18)+
    geom_hline(data=ctrl_df, aes(yintercept=med, color=group),
               linetype="dashed", linewidth=0.6)+
    geom_hline(yintercept=0, linetype="dashed", color="gray40", linewidth=0.5)+
    geom_line(color="#D55E00", linewidth=0.8)+
    geom_area(alpha=0.15, fill="#D55E00")+
    scale_fill_manual(values=ctrl_cols, name=NULL)+
    scale_color_manual(values=ctrl_cols, name=NULL)+
    scale_x_continuous(breaks=dd$pool[idx], labels=sprintf("%.2g", dd$x[idx]))+
    coord_cartesian(ylim=ylim)+
    labs(x=x_title, y="ΔS", title=paste0("ΔS (", COND_LEVELS[1], " − ", COND_LEVELS[2], ")"))+
    theme_classic(base_size=11)+
    theme(plot.title=element_text(face="bold", size=11, hjust=0.5),
          axis.text=element_text(size=9, color="black"),
          axis.title=element_text(size=10, color="black"),
          panel.grid.major.y=element_line(color="grey92", linewidth=0.3),
          panel.grid.major.x=element_blank(),
          legend.position="bottom",
          plot.margin=margin(4, 16, 8, 8))
}

# S(d) and delta-S(d) panels, sized like the paper panels (plot area height / width = SD_ASPECT)
SD_ASPECT <- 0.46

make_sd_main <- function(deriv_df, cfg) {
  deriv_df <- deriv_df %>% dplyr::filter(!group %in% PLOT_DROP_GROUPS) %>% dplyr::mutate(group=droplevels(factor(group)))
  ggplot(deriv_df, aes(x=dist_bp/1000, y=slope, color=group, linetype=cond))+
    geom_hline(yintercept=0, linetype="dotted", color="gray50", linewidth=0.4)+
    geom_line(linewidth=0.8, alpha=0.85)+
    scale_color_manual(values=SD_COLORS, name=NULL)+
    scale_linetype_manual(values=c("solid", "dashed"), name=NULL)+
    scale_x_continuous(breaks=cfg$axis_kb, labels=cfg$axis_lab)+
    labs(x="Genomic distance", y="S(d)")+
    theme_classic(base_size=13)+
    theme(legend.position="top", aspect.ratio=SD_ASPECT,
          axis.text.x=element_text(angle=45, hjust=1),
          panel.grid.major.y=element_line(color="grey92", linewidth=0.3))
}

make_sd_delta <- function(delta_df, cfg, ci_df) {
  delta_df <- delta_df %>% dplyr::filter(!group %in% PLOT_DROP_GROUPS)
  ci_df    <- ci_df    %>% dplyr::filter(!group %in% PLOT_DROP_GROUPS)
  ggplot(delta_df, aes(x=dist_bp/1000, y=delta, color=group))+
    geom_hline(yintercept=0, linetype="dashed", color="gray40", linewidth=0.5)+
    geom_ribbon(data=ci_df, aes(x=dist_bp/1000, ymin=lo, ymax=hi, fill=group),
                inherit.aes=FALSE, alpha=0.15)+
    scale_fill_manual(values=SD_COLORS, guide="none")+
    geom_line(linewidth=0.8)+
    scale_color_manual(values=SD_COLORS, name=NULL)+
    scale_x_continuous(breaks=cfg$axis_kb, labels=cfg$axis_lab)+
    coord_cartesian(ylim=cfg$delta_ylim)+
    labs(x="Genomic distance", y="ΔS(d)")+
    theme_classic(base_size=13)+
    theme(legend.position="none", aspect.ratio=SD_ASPECT,
          axis.text.x=element_text(angle=45, hjust=1),
          panel.grid.major.y=element_line(color="grey92", linewidth=0.3))
}

make_nbhd_plot <- function(act_df, title) {
  lv <- levels(factor(act_df$band, levels=unique(act_df$band)))
  s <- act_df %>%
    dplyr::filter(!is.na(tpm_per_mb), !group %in% PLOT_DROP_GROUPS) %>%
    dplyr::group_by(group, band) %>%
    dplyr::summarise(med=median(tpm_per_mb), q25=quantile(tpm_per_mb, 0.25),
                     q75=quantile(tpm_per_mb, 0.75), .groups="drop") %>%
    dplyr::mutate(band=factor(band, levels=lv), group=factor(group, levels=PLOT_LEVELS))
  dodge <- position_dodge(width=0.55)
  ggplot(s, aes(x=band, y=med, color=group, group=group))+
    geom_line(position=dodge, linewidth=0.6, alpha=0.8)+
    geom_pointrange(aes(ymin=q25, ymax=q75), position=dodge, size=0.45, linewidth=0.6)+
    scale_color_manual(values=SD_COLORS, name=NULL, drop=FALSE)+
    scale_y_continuous(trans=scales::pseudo_log_trans(base=10),
                       breaks=c(0, 1, 10, 100, 1000, 10000),
                       labels=c("0", "1", "10", "100", "1000", "10000"))+
    labs(x="Distance from region edge (both sides)", y="Flanking activity (TPM per Mb)",
         title=title, subtitle="median and interquartile range across regions; genes assigned by midpoint")+
    theme_classic(base_size=13)+
    theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
          plot.subtitle=element_text(color="grey40", size=10, hjust=0.5),
          legend.position="right",
          panel.grid.major.y=element_line(color="grey92", linewidth=0.3))
}

# Save one plot as SVG in PLOT_DIR (no-op when PLOT_DIR is NULL). File name:
# <assay>_<analysis>_<stratum>_compartment-<A|B|all>_<plot type>.svg
save_plot_svg <- function(p, type, label, height=PLOT_H) {
  if (is.null(PLOT_DIR)) return(invisible(NULL))
  nm <- paste(ASSAY, ANALYSIS, gsub(">=", "ge", label), paste0("compartment-", COMPARTMENT), type, sep="_")
  f <- file.path(PLOT_DIR, paste0(gsub("[^A-Za-z0-9._-]+", "-", nm), ".svg"))
  svg(f, width=PLOT_W, height=height)
  print(p)
  dev.off()
  message("  saved ", f)
  invisible(f)
}

print_stratum_plots <- function(res_obj, cfg) {
  tag <- paste0(ASSAY, "  |  ", ANALYSIS, "  |  ", res_obj$label, "  |  compartment: ", COMPARTMENT)
  ctrl_df <- res_obj$delta %>%
    dplyr::filter(group %in% setdiff(c("Non-expressed genes", "Non-genic windows"), PLOT_DROP_GROUPS)) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(q25=quantile(delta, 0.25), med=median(delta), q75=quantile(delta, 0.75),
                     .groups="drop")
  ylim_d <- range(c(res_obj$delta$delta, ctrl_df$q25, ctrl_df$q75), na.rm=TRUE)
  ylim_d <- ylim_d + c(-1, 1) * 0.05 * diff(ylim_d)
  p_d_ei  <- make_delta_pools(res_obj$delta, res_obj$pool_x, "EI",  "E/I ratio (bin median)", ctrl_df, ylim_d)+
    theme(legend.position="none")
  p_d_tpm <- make_delta_pools(res_obj$delta, res_obj$pool_x, "TPM", "Neighborhood TPM (bin median)", ctrl_df, ylim_d)
  title_grob <- cowplot::ggdraw()+
    cowplot::draw_label(paste0("ΔS vs Controls  |  ", tag), fontface="bold", size=13, x=0.01, hjust=0)
  p_delta <- cowplot::plot_grid(title_grob,
                                cowplot::plot_grid(p_d_ei, p_d_tpm, ncol=1, rel_heights=c(1, 1.15)),
                                ncol=1, rel_heights=c(0.06, 1))
  print(p_delta)
  save_plot_svg(p_delta, "deltaS-vs-controls", res_obj$label)
  
  p_sd  <- make_sd_main(res_obj$loess, cfg)
  p_sdd <- make_sd_delta(res_obj$delta_loess, cfg, res_obj$delta_ci)
  leg   <- cowplot::get_legend(p_sd + theme(legend.box="horizontal"))
  sd_title <- cowplot::ggdraw()+
    cowplot::draw_label(paste0("S(d) and ΔS(d)  |  ", tag, "  |  loess span=0.3  |  bands: 95% region bootstrap"),
                        size=10, x=0.01, hjust=0)
  p_loess <- cowplot::plot_grid(sd_title, leg,
                                cowplot::plot_grid(p_sd + theme(legend.position="none"), p_sdd, ncol=2, align="h", axis="tb"),
                                ncol=1, rel_heights=c(0.07, 0.12, 1))
  print(p_loess)
  save_plot_svg(p_loess, "loess-Sd-and-deltaSd", res_obj$label, height=PLOT_H_LOESS)
  p_nbhd <- make_nbhd_plot(res_obj$nbhd_act, paste0("Transcriptional activity flanking each group  |  ", tag))
  print(p_nbhd)
  save_plot_svg(p_nbhd, "flanking-activity", res_obj$label)
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. Checks on the .hic files
# ══════════════════════════════════════════════════════════════════════════════
check_hic_access <- function(hic_file, cfg, chrom_sizes, norm=STRAW_NORM) {
  stopifnot(file.exists(hic_file))
  stopifnot(cfg$res %in% strawr::readHicBpResolutions(hic_file))
  stopifnot(norm %in% strawr::readHicNormTypes(hic_file))
  hc   <- hic_chrom_map(hic_file)[[chrom_sizes$chr[1]]]
  p0   <- floor(chrom_sizes$length[1] * 0.4 / 1e5) * 1e5   # inside the first chromosome
  test <- strawr::straw(norm, hic_file,
                        sprintf("%s:%.0f:%.0f", hc, p0, p0 + 1e5),
                        sprintf("%s:%.0f:%.0f", hc, p0, p0 + 2e5),
                        "BP", cfg$res)
  stopifnot(nrow(test) > 0, any(is.finite(test$counts)))
  message("  CHECK ", basename(hic_file), ": res ", cfg$res, " / ", norm,
          " readable (", nrow(test), " test pixels)")
}

# Raw (NONE) cis contacts at 25kb: total and fraction at distance >= 25kb
summarize_depth <- function(hic_file, chrom_sizes, n_cores=N_CORES) {
  hmap <- hic_chrom_map(hic_file)
  per_chr <- parallel::mclapply(chrom_sizes$chr, function(ch) {
    raw <- strawr::straw("NONE", hic_file, hmap[[ch]], hmap[[ch]], "BP", 25000L)
    dt  <- data.table::as.data.table(raw)
    c(total=sum(dt$counts), long=sum(dt$counts[abs(dt$x - dt$y) >= 25000]))
  }, mc.cores=n_cores, mc.preschedule=FALSE)
  m <- do.call(rbind, per_chr)
  data.frame(cis_contacts=sum(m[, "total"]),
             frac_ge_25kb=sum(m[, "long"]) / sum(m[, "total"]))
}

fmt_bp <- function(x) ifelse(x >= 1e6, paste0(x/1e6, "Mb"),
                             ifelse(x >= 1e3, paste0(x/1e3, "kb"), paste0(x, "bp")))

# ══════════════════════════════════════════════════════════════════════════════
# 8. Driver. The fork scripts set CONDITIONS, ASSAY, RDS_BASE, ANALYSES,
#    COMPARTMENT (and COMP_BED when COMPARTMENT is "A" or "B"), then call this.
# ══════════════════════════════════════════════════════════════════════════════
run_scaling_analysis <- function() {
  need <- c("CONDITIONS", "ASSAY", "RDS_BASE", "ANALYSES", "COMPARTMENT")
  miss <- need[!vapply(need, exists, logical(1))]
  if (length(miss) > 0) stop("Set before run_scaling_analysis(): ", paste(miss, collapse=", "))
  stopifnot(COMPARTMENT %in% c("all", "A", "B"), all(ANALYSES %in% names(ANALYSIS_CFG)))
  COND_LEVELS <<- sapply(CONDITIONS, `[[`, "display")
  COND_COLORS <<- setNames(sapply(CONDITIONS, `[[`, "color"), COND_LEVELS)
  
  # every input file must exist before anything runs
  input_files <- c(setNames(sapply(CONDITIONS, `[[`, "hic"), paste0(names(CONDITIONS), "_hic")),
                   gtf=GTF_FILE, expression=EXPR_FILE)
  if (COMPARTMENT != "all") input_files <- c(input_files, compartments=COMP_BED)
  missing_files <- input_files[!file.exists(input_files)]
  if (length(missing_files) > 0)
    stop("Missing input files:\n", paste0("  ", names(missing_files), ": ", missing_files, collapse="\n"))
  if (!is.null(PLOT_DIR) && !dir.exists(PLOT_DIR)) stop("PLOT_DIR not found: ", PLOT_DIR)
  message("All ", length(input_files), " input files found")
  message(ASSAY, "  |  analyses: ", paste(ANALYSES, collapse=", "),
          "  |  compartment: ", COMPARTMENT, "  |  cores: ", N_CORES)
  dir.create(RDS_BASE, showWarnings=FALSE, recursive=TRUE)
  
  # chromosome sizes from the .hic files (canonical chrN names)
  cs <- lapply(CONDITIONS, function(cd) {
    ct <- strawr::readHicChroms(cd$hic)
    mp <- hic_chrom_map(cd$hic)
    data.frame(chr=names(mp), length=ct$length[match(mp, ct$name)], stringsAsFactors=FALSE)
  })
  chrom_sizes <- cs[[1]]
  for (k in cs) stopifnot(identical(k$length, chrom_sizes$length))
  
  depth_rds <- file.path(RDS_BASE, paste0("depth_summary_", ASSAY, ".rds"))
  if (file.exists(depth_rds)) {
    depth <- readRDS(depth_rds)
  } else {
    depth <- dplyr::bind_rows(lapply(names(CONDITIONS), function(n) {
      d <- summarize_depth(CONDITIONS[[n]]$hic, chrom_sizes)
      d$cond <- CONDITIONS[[n]]$display
      d
    }))
    saveRDS(depth, depth_rds)
  }
  message("\n  CHECK depth (raw cis contacts, 25kb):"); print(depth)
  message("  cis ratio (", COND_LEVELS[2], "/", COND_LEVELS[1], "): ",
          round(depth$cis_contacts[2] / depth$cis_contacts[1], 3))
  
  message("\n── Annotation, expression, neighborhoods ──")
  all_genes <- load_gene_table(GTF_FILE, RDS_BASE) %>%
    dplyr::left_join(load_expression(EXPR_FILE), by="gene_id")
  message("  genes: ", nrow(all_genes), "  TPM missing=", sum(is.na(all_genes$TPM)),
          " (counted as unexpressed in neighborhoods)  zero=", sum(all_genes$TPM == 0, na.rm=TRUE),
          "  >0=", sum(all_genes$TPM > 0, na.rm=TRUE))
  all_genes <- add_neighborhood(all_genes, NBHD_BP, NONEXPR_TPM_MAX)
  all_genes$n_expr_ctrl <- add_neighborhood(all_genes, CTRL_NBHD_BP, CTRL_TPM_MAX)$n_expr_nbhd
  message("  genes with no gene of TPM>", CTRL_TPM_MAX, " within ", CTRL_NBHD_BP / 1000, " kb (incl. self): ",
          sum(all_genes$n_expr_ctrl == 0), "  | median neighborhood TPM of expressed genes: ",
          round(median(all_genes$TPM_nbhd[!is.na(all_genes$TPM) & all_genes$TPM > 0]), 3))
  ei_df <- all_genes %>%
    dplyr::filter(exon_bp > 0, !is.infinite(EI), EI > 0, EI != 1)
  
  eigen_df <- NULL
  if (COMPARTMENT != "all") {
    message("\n── Compartments ──")
    tpm_high <- quantile(ei_df$TPM[ei_df$TPM > 0], 2/3, na.rm=TRUE)
    eigen_df <- load_eigen_wt(COMP_BED, ei_df %>% dplyr::filter(!is.na(TPM), TPM > tpm_high))
  }
  
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
      lo <- cfg$strata[[si]][1]
      hi <- cfg$strata[[si]][2]
      label <- if (is.infinite(hi)) paste0(">=", fmt_bp(lo)) else paste0(fmt_bp(lo), "-", fmt_bp(hi))
      genes_s <- ei_df %>% dplyr::filter(size >= lo, size < hi)
      results[[a]][[si]] <- run_stratum(genes_s, all_genes, chrom_sizes, eigen_df,
                                        cfg, si, label, rds_dir)
      print_stratum_plots(results[[a]][[si]], cfg)
    }
  }
  invisible(results)
}