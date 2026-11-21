# Analysis 3: same as Analysis 1 (change in S vs change in nascent transcription and change in
# H3K27ac / H3K4me3 / H3K9me3), but the nascent transcription change of each gene is replaced by
# the change of its gene neighborhood. Everything else stays gene-based (mark changes, length
# matching, pooling, tertiles, bootstrap).
#
# Neighborhood = the gene plus all genes within NBHD_BP (100 kb, same as add_neighborhood() in
# scaling_controls_main.R). Neighborhood EU-seq = mean RPKM over those genes, per condition
# (Wildtype: WT_1..3, POLR2A deg: POLR2A_1..3), then
#   change = log2((WT_nbhd + EUSEQ_PSEUDO) / (deg_nbhd + EUSEQ_PSEUDO));  positive = decrease.
#
# Settings (ASSAY_RUN, EUSEQ_FILE, DATA_DIR, N_BOOT, ...) live in scaling_chromatin_vs_transcription.R.
# Outputs carry the tag "_neighborhood", so Analysis 1 files are not overwritten.

.find_script_dir2 <- function() {
  f <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
  if (is.null(f)) {
    a <- grep("^--file=", commandArgs(FALSE), value=TRUE)
    if (length(a)) f <- sub("^--file=", "", a[1])
  }
  if (is.null(f) && requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable())
    f <- tryCatch(rstudioapi::getSourceEditorContext()$path, error=function(e) NULL)
  if (!is.null(f) && nzchar(f)) dirname(normalizePath(f)) else getwd()
}
SCRIPT_DIR2 <- .find_script_dir2()
BASE_SCRIPT <- file.path(SCRIPT_DIR2, "scaling_chromatin_vs_transcription.R")
if (!file.exists(BASE_SCRIPT)) stop("Missing: ", BASE_SCRIPT)
options(chromatin_script_no_run=TRUE)
source(BASE_SCRIPT)

# ------------------------------------------------------------------------------------------
# Neighborhood nascent transcription change. Returns data.frame(gene_id, loss) like load_euseq_loss().
# ------------------------------------------------------------------------------------------
load_euseq_neighborhood_loss <- function(file, all_genes, bp=NBHD_BP) {
  if (!file.exists(file)) stop("Missing EU-seq file: ", file)
  d <- read.delim(file, check.names=FALSE, stringsAsFactors=FALSE)
  need <- c("gene_id", paste0("POLR2A_", 1:3), paste0("WT_", 1:3))
  if (!all(need %in% names(d))) stop("EU-seq file lacks columns: ", paste(setdiff(need, names(d)), collapse=", "))
  d$gene_id <- sub("\\.\\d+$", "", d$gene_id)
  d <- d[!duplicated(d$gene_id), ]
  wt  <- rowMeans(as.matrix(d[, paste0("WT_", 1:3)]))
  deg <- rowMeans(as.matrix(d[, paste0("POLR2A_", 1:3)]))
  names(wt) <- names(deg) <- d$gene_id
  
  g <- all_genes
  g$gene_id <- sub("\\.\\d+$", "", g$gene_id)
  gr <- GenomicRanges::GRanges(g$chr, IRanges::IRanges(g$start, g$end))
  ov <- GenomicRanges::findOverlaps(gr, gr, maxgap=bp)       # includes the gene itself
  q <- S4Vectors::queryHits(ov); s <- S4Vectors::subjectHits(ov)
  wt_s  <- unname(wt[g$gene_id[s]]); deg_s <- unname(deg[g$gene_id[s]])
  ok <- !is.na(wt_s) & !is.na(deg_s)
  n_genes <- tabulate(q[ok], nbins=nrow(g))
  wt_n  <- rowsum(wt_s[ok],  q[ok]); deg_n <- rowsum(deg_s[ok], q[ok])
  wt_m  <- rep(NA_real_, nrow(g)); deg_m <- rep(NA_real_, nrow(g))
  idx <- as.integer(rownames(wt_n))
  wt_m[idx] <- wt_n[, 1] / n_genes[idx]; deg_m[idx] <- deg_n[, 1] / n_genes[idx]
  
  message(sprintf("Neighborhood EU-seq (+/-%s bp): median %d genes per neighborhood (genes with EU-seq values); %d of %d genes without any",
                  format(bp, big.mark=","), as.integer(median(n_genes[n_genes > 0])), sum(n_genes == 0), nrow(g)))
  data.frame(gene_id=g$gene_id,
             loss=log2((wt_m + EUSEQ_PSEUDO) / (deg_m + EUSEQ_PSEUDO)),
             stringsAsFactors=FALSE)
}

# every plot is saved here as an SVG (folder must already exist); size in inches, same for every plot
PLOT_DIR <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/micro-c/plots/scaling_v_neighbor_transcription_plots"
PLOT_W   <- 14
PLOT_H   <- 10

results_neighborhood <- run_chromatin_analysis(
  "change",
  nascent_fun=function(all_genes) load_euseq_neighborhood_loss(EUSEQ_FILE, all_genes),
  tag="_neighborhood",
  title="Analysis 3 (change in neighborhood transcription vs change in marks)",
  loss_name="Neighborhood nascent transcription")