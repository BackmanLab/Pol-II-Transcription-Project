##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Metagene analysis of log2(POLR2A/WT) signal for EU-seq and
# total RNA-seq. Stratified by E/I ratio, expression group,
# and DEG status across four plots:
#
# hg38 version: signal from the chrM-scaled hg38 bigWigs
#   EU-seq:        euseq/coverage/log2ratio_bamcompare/POLR2A_vs_WT.log2ratio.bamcompare.bigWig
#                  (generate_euseq_log2ratio_bamcompare_hg38.sh; bamCompare on the
#                  merged EU-seq BAMs, chrM-scaled, 20 bp bins, 60 bp smoothing)
#   total RNA-seq: total/coverage/log2ratio_mito/POLR2A_vs_WT.log2ratio.bigWig
#                  (generate_total_bigwigs_mito_hg38.sh)
# Genes are classified with the previous total RNA-seq DESeq2 results
# (Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv) and the ENCODE HCT116
# expression file (ENCFF435PHM).
#
#   Plot 1: facet by E/I group only (all genes)
#   Plot 2: facet_grid(de_status ~ ei_group) — Up/Down DEGs
#   Plot 3: facet_grid(ei_group ~ expr_group) — all genes
#   Plot 4: facet by expr_group, Up and Down DEGs separate
#
# Two separate profile summaries computed from the same raw signal:
#   profile_ei_df  — grouped by ei_group + expr_group + de_status
#   profile_deg_df — grouped by expr_group + de_status only
#
# Summarisation uses hardcoded group_by columns (not programmatic
# column selection) to avoid dplyr scoping issues with lapply.
#
# All factor stratification columns converted to character before
# rep.int() in bin construction to prevent silent integer coercion.
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(rtracklayer)
  library(plyranges)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(BiocParallel)
  library(patchwork)
  library(scales)
})

reduce <- GenomicRanges::reduce
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename
n      <- dplyr::n

##--------------------------------------------------------------## 1. Parameters

BASE       <- "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq"
GTF_PATH   <- file.path(BASE, "reference/Homo_sapiens.GRCh38.115.chr.gtf")
EXPR_FILE  <- file.path(BASE, "reference/ENCFF435PHM.tsv")
DESEQ_FILE <- file.path(BASE, "reference/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv")

EU_DIR  <- file.path(BASE, "euseq/coverage/log2ratio_bamcompare")
RNA_DIR <- file.path(BASE, "total/coverage/log2ratio_mito")

BW_FILES <- c(
  EU_log2  = file.path(EU_DIR,  "POLR2A_vs_WT.log2ratio.bamcompare.bigWig"),
  RNA_log2 = file.path(RNA_DIR, "POLR2A_vs_WT.log2ratio.bigWig")
)

COND_LEVELS <- c("EU_log2", "RNA_log2")
COLORS <- c(
  EU_log2  = "#E69F00",
  RNA_log2 = "#0072B2"
)

LABELS <- c(
  EU_log2  = "EU-seq log2(POLR2A/WT)",
  RNA_log2 = "Total RNA-seq log2(POLR2A/WT)"
)

LFC_THRESHOLD  <- 1
PADJ_THRESHOLD <- 0.05
TPM_THRESHOLD  <- 1
STD_CHROMS     <- paste0("chr", c(1:22, "X", "Y"))
N_WORKERS      <- 7

N_BINS_UPSTREAM   <- 10
N_BINS_BODY       <- 50
N_BINS_DOWNSTREAM <- 10
N_BINS_INTRON     <- 100
FLANK             <- 2000

MIN_WIDTH_BODY   <- N_BINS_BODY   + 1L
MIN_WIDTH_INTRON <- N_BINS_INTRON + 1L

## Every plot is printed to the live device and then saved here as an SVG
## (the folder must already exist; set to NULL to save nothing)
PLOT_DIR <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/transcriptomics/plots/metagene analysis plots"

##--------------------------------------------------------------## 2. Validate paths

message("Checking file paths...")
missing_bw <- BW_FILES[!file.exists(BW_FILES)]
if (length(missing_bw) > 0) stop("Missing bigWigs:\n",
                                 paste(" -", missing_bw, collapse = "\n"))
if (!file.exists(GTF_PATH))   stop("GTF not found: ",      GTF_PATH)
if (!file.exists(EXPR_FILE))  stop("TPM file not found: ",  EXPR_FILE)
if (!file.exists(DESEQ_FILE)) stop("DESeq2 file not found: ", DESEQ_FILE)
if (!is.null(PLOT_DIR) && !dir.exists(PLOT_DIR)) stop("PLOT_DIR not found: ", PLOT_DIR)
message("  All files found.")

##--------------------------------------------------------------## 3. Build TxDb

message("Building TxDb from GTF...")
txdb <- makeTxDbFromGFF(GTF_PATH, format = "gtf")

genes_gr         <- genes(txdb)
genes_gr$gene_id <- sub("\\..*", "", names(genes_gr))
seqlevelsStyle(genes_gr) <- "UCSC"
genes_gr <- keepSeqlevels(genes_gr,
                          intersect(seqlevels(genes_gr), STD_CHROMS),
                          pruning.mode = "coarse")

exons_by_gene_gr <- exonsBy(txdb, by = "gene")
names(exons_by_gene_gr) <- sub("\\..*", "", names(exons_by_gene_gr))
seqlevelsStyle(exons_by_gene_gr) <- "UCSC"
exons_by_gene_gr <- keepSeqlevels(exons_by_gene_gr,
                                  intersect(seqlevels(exons_by_gene_gr),
                                            STD_CHROMS),
                                  pruning.mode = "coarse")

tx2gene <- AnnotationDbi::select(txdb,
                                 keys    = keys(txdb, "TXNAME"),
                                 columns = c("TXNAME", "GENEID"),
                                 keytype = "TXNAME")
introns_by_tx_gr   <- intronsByTranscript(txdb, use.names = TRUE)
introns_unlisted   <- unlist(introns_by_tx_gr, use.names = TRUE)
intron_gene_ids    <- tx2gene$GENEID[match(names(introns_unlisted),
                                           tx2gene$TXNAME)]
introns_by_gene_gr <- splitAsList(introns_unlisted,
                                  sub("\\..*", "", intron_gene_ids))
seqlevelsStyle(introns_by_gene_gr) <- "UCSC"
introns_by_gene_gr <- keepSeqlevels(introns_by_gene_gr,
                                    intersect(seqlevels(introns_by_gene_gr),
                                              STD_CHROMS),
                                    pruning.mode = "coarse")
rm(introns_by_tx_gr, introns_unlisted, intron_gene_ids); gc()

message("  TxDb built: ", length(genes_gr), " genes")

##--------------------------------------------------------------## 4. E/I ratio (all genes)

message("Computing E/I ratio...")
exon_bp <- sum(width(reduce(exons_by_gene_gr)))

gene_lengths <- dplyr::tibble(
  gene_id     = sub("\\..*", "", names(genes_gr)),
  gene_length = as.numeric(width(genes_gr))
)

exon_bp_df <- dplyr::tibble(
  gene_id = sub("\\..*", "", names(exon_bp)),
  exon_bp = as.numeric(exon_bp)
) %>%
  dplyr::left_join(gene_lengths, by = "gene_id") %>%
  dplyr::mutate(
    intron_bp = gene_length - exon_bp,
    ei_ratio  = dplyr::if_else(intron_bp > 0,
                               exon_bp / intron_bp,
                               NA_real_),
    ei_group  = dplyr::case_when(
      is.na(ei_ratio)                                    ~ "High E/I",
      ei_ratio < quantile(ei_ratio, 0.33, na.rm = TRUE) ~ "Low E/I",
      ei_ratio < quantile(ei_ratio, 0.67, na.rm = TRUE) ~ "Medium E/I",
      TRUE                                               ~ "High E/I"
    ),
    ei_group = factor(ei_group,
                      levels = c("Low E/I", "Medium E/I", "High E/I"))
  )

message("E/I ratio quantile breaks:")
message(paste(
  c("0%", "33%", "67%", "100%"),
  round(quantile(exon_bp_df$ei_ratio, c(0, 0.33, 0.67, 1),
                 na.rm = TRUE), 3),
  sep = ": ", collapse = " | "
))

cat("E/I group counts:\n")
print(table(exon_bp_df$ei_group, useNA = "always"))

##--------------------------------------------------------------## 5. Expression annotation (TPM)

message("Loading TPM data...")
expr_raw <- read_tsv(EXPR_FILE, show_col_types = FALSE)
cat("Expression file columns:", paste(colnames(expr_raw), collapse = ", "), "\n")
stopifnot(all(c("gene_id", "TPM") %in% colnames(expr_raw)))
expr <- expr_raw %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
  dplyr::select(gene_id, TPM) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

genes_meta <- dplyr::tibble(gene_id = genes_gr$gene_id) %>%
  dplyr::left_join(expr, by = "gene_id") %>%
  dplyr::mutate(TPM = tidyr::replace_na(TPM, 0))

tpm_expressed <- genes_meta$TPM[genes_meta$TPM >= TPM_THRESHOLD]
genes_meta <- genes_meta %>%
  dplyr::mutate(
    expr_group = dplyr::case_when(
      TPM < TPM_THRESHOLD                                     ~ "Unexpressed",
      TPM < quantile(tpm_expressed, 1/3, na.rm = TRUE)       ~ "Low",
      TPM < quantile(tpm_expressed, 2/3, na.rm = TRUE)       ~ "Medium",
      TRUE                                                    ~ "High"
    ),
    expr_group = factor(expr_group,
                        levels = c("Unexpressed", "Low", "Medium", "High"))
  )

##--------------------------------------------------------------## 6. DESeq2 annotation

message("Loading DESeq2 results...")
deseq_raw <- read_csv(DESEQ_FILE, show_col_types = FALSE)
cat("DESeq2 file columns:", paste(colnames(deseq_raw), collapse = ", "), "\n")
stopifnot(all(c("gene_id", "padj", "log2FoldChange") %in% colnames(deseq_raw)))
deseq_res <- deseq_raw %>%
  dplyr::mutate(
    gene_id   = sub("\\..*", "", gene_id),
    de_status = dplyr::case_when(
      is.na(padj)                                               ~ "Filtered",
      padj >= PADJ_THRESHOLD                                    ~ "Not DE",
      padj <  PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD ~ "Up",
      padj <  PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD ~ "Down",
      TRUE                                                      ~ "Not DE"
    ),
    de_status = factor(de_status,
                       levels = c("Up", "Not DE", "Down", "Filtered"))
  )

cat("DE status counts:\n")
print(table(deseq_res$de_status))

##--------------------------------------------------------------## 7. Merge all metadata onto genes_gr

genes_meta <- genes_meta %>%
  dplyr::left_join(
    exon_bp_df %>% dplyr::select(gene_id, ei_ratio, ei_group),
    by = "gene_id"
  ) %>%
  dplyr::left_join(
    deseq_res %>% dplyr::select(gene_id, de_status,
                                log2FoldChange, padj),
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    de_status = tidyr::replace_na(as.character(de_status), "Not DE"),
    de_status = factor(de_status,
                       levels = c("Up", "Not DE", "Down", "Filtered"))
  )

genes_gr$TPM            <- genes_meta$TPM
genes_gr$expr_group     <- genes_meta$expr_group
genes_gr$ei_group       <- genes_meta$ei_group
genes_gr$ei_ratio       <- genes_meta$ei_ratio
genes_gr$de_status      <- genes_meta$de_status
genes_gr$log2FoldChange <- genes_meta$log2FoldChange
genes_gr$padj           <- genes_meta$padj

cat("\nE/I x DEG cross-tab:\n")
print(table(genes_gr$ei_group, genes_gr$de_status, useNA = "always"))
cat("\nE/I x Expression cross-tab:\n")
print(table(genes_gr$ei_group, genes_gr$expr_group, useNA = "always"))
cat("\nExpression x DEG cross-tab (DEG calls from ", basename(DESEQ_FILE), "):\n", sep = "")
print(table(genes_gr$expr_group, genes_gr$de_status, useNA = "always"))

##--------------------------------------------------------------## 8. safe_tile() helper

safe_tile <- function(gr, n_bins, label = "features") {
  tiled       <- tile(gr, n = n_bins)
  tile_counts <- lengths(tiled)
  bad         <- tile_counts != n_bins
  if (any(bad)) {
    message("  Dropping ", sum(bad), " ", label,
            " where tile() returned != ", n_bins, " bins")
    gr    <- gr[!bad]
    tiled <- tiled[!bad]
  }
  list(gr = gr, tiled = tiled)
}

##--------------------------------------------------------------## 9. Build bin GRanges
## All factor columns converted to character before rep.int()
## to prevent silent coercion of factor levels to integers.

make_feature_bins <- function(gr, feature_label, n_bins, min_width) {
  if (feature_label == "Upstream") {
    tss_gr  <- GenomicRanges::resize(gr, width = 1, fix = "start")
    feat_gr <- GenomicRanges::flank(tss_gr, width = FLANK,
                                    start = TRUE, both = FALSE)
    feat_gr <- feat_gr[width(feat_gr) > n_bins]
  } else if (feature_label == "Body") {
    feat_gr <- gr[width(gr) > min_width]
  } else if (feature_label == "Downstream") {
    tes_gr  <- GenomicRanges::resize(gr, width = 1, fix = "end")
    feat_gr <- GenomicRanges::flank(tes_gr, width = FLANK,
                                    start = FALSE, both = FALSE)
    feat_gr <- feat_gr[width(feat_gr) > n_bins]
  }
  
  feat_gr <- keepSeqlevels(feat_gr,
                           intersect(seqlevels(feat_gr), STD_CHROMS),
                           pruning.mode = "coarse")
  
  res     <- safe_tile(feat_gr, n_bins, paste(feature_label, "bins"))
  feat_gr <- res$gr
  tiled   <- res$tiled
  
  flat        <- unlist(tiled, use.names = FALSE)
  n_features  <- length(feat_gr)
  bin_idx     <- rep.int(seq_len(n_bins), n_features)
  feature_idx <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus    <- as.logical(strand(feat_gr) == "-")
  minus_mask  <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = feature_label,
    gene_id     = rep.int(feat_gr$gene_id,
                          rep.int(n_bins, n_features)),
    ei_group    = rep.int(as.character(feat_gr$ei_group),
                          rep.int(n_bins, n_features)),
    expr_group  = rep.int(as.character(feat_gr$expr_group),
                          rep.int(n_bins, n_features)),
    de_status   = rep.int(as.character(feat_gr$de_status),
                          rep.int(n_bins, n_features))
  )
  flat
}

make_intron_bins <- function(introns_by_gene, genes_gr,
                             n_bins    = N_BINS_INTRON,
                             min_width = MIN_WIDTH_INTRON) {
  introns_red      <- GenomicRanges::reduce(introns_by_gene)
  flat_introns     <- unlist(introns_red, use.names = TRUE)
  flat_introns$gene_id <- sub("\\..*", "", names(flat_introns))
  flat_introns     <- flat_introns[width(flat_introns) > min_width]
  seqlevelsStyle(flat_introns) <- "UCSC"
  flat_introns <- keepSeqlevels(flat_introns,
                                intersect(seqlevels(flat_introns), STD_CHROMS),
                                pruning.mode = "coarse")
  
  ei_map   <- setNames(as.character(genes_gr$ei_group),   genes_gr$gene_id)
  expr_map <- setNames(as.character(genes_gr$expr_group), genes_gr$gene_id)
  de_map   <- setNames(as.character(genes_gr$de_status),  genes_gr$gene_id)
  
  flat_introns$ei_group   <- ei_map[flat_introns$gene_id]
  flat_introns$expr_group <- expr_map[flat_introns$gene_id]
  flat_introns$de_status  <- de_map[flat_introns$gene_id]
  
  res          <- safe_tile(flat_introns, n_bins, "intron bins")
  flat_introns <- res$gr
  tiled        <- res$tiled
  
  flat        <- unlist(tiled, use.names = FALSE)
  n_features  <- length(flat_introns)
  bin_idx     <- rep.int(seq_len(n_bins), n_features)
  feature_idx <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus    <- as.logical(strand(flat_introns) == "-")
  minus_mask  <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "Intron",
    gene_id     = rep.int(flat_introns$gene_id,
                          rep.int(n_bins, n_features)),
    ei_group    = rep.int(flat_introns$ei_group,
                          rep.int(n_bins, n_features)),
    expr_group  = rep.int(flat_introns$expr_group,
                          rep.int(n_bins, n_features)),
    de_status   = rep.int(flat_introns$de_status,
                          rep.int(n_bins, n_features))
  )
  flat
}

cat("Building upstream bins...\n")
upstream_bins   <- make_feature_bins(genes_gr, "Upstream",
                                     N_BINS_UPSTREAM, N_BINS_UPSTREAM + 1L)
cat("Building body bins...\n")
body_bins       <- make_feature_bins(genes_gr, "Body",
                                     N_BINS_BODY, MIN_WIDTH_BODY)
cat("Building downstream bins...\n")
downstream_bins <- make_feature_bins(genes_gr, "Downstream",
                                     N_BINS_DOWNSTREAM, N_BINS_DOWNSTREAM + 1L)
cat("Building intron bins...\n")
intron_bins     <- make_intron_bins(introns_by_gene_gr, genes_gr)

metagene_bins <- c(upstream_bins, body_bins, downstream_bins)
rm(upstream_bins, body_bins, downstream_bins); gc()

cat("de_status in metagene_bins:\n")
print(table(mcols(metagene_bins)$de_status, useNA = "always"))

##--------------------------------------------------------------## 10. Signal extraction

## mean_signal = length-weighted mean of the bigWig over the positions of the bin that
## carry a bigWig interval (positions without an interval are skipped, not counted as 0;
## bins with no interval are dropped). covered = share of the bin width inside a bigWig
## interval. Overlap widths come from the explicit intersection of bin and interval.
extract_profile_signal_one_chrom <- function(chrom, bw_gr, bins) {
  bw_chrom   <- bw_gr[as.character(seqnames(bw_gr))  == chrom]
  bins_chrom <- bins[as.character(seqnames(bins))     == chrom]
  if (length(bw_chrom) == 0 || length(bins_chrom) == 0) return(NULL)
  bins_chrom$bin_id    <- seq_along(bins_chrom)
  bins_chrom$bin_width <- BiocGenerics::width(bins_chrom)
  tryCatch({
    covered_bins <- plyranges::join_overlap_intersect(bins_chrom, bw_chrom) %>%
      dplyr::as_tibble() %>%
      dplyr::select(bin_id, overlap_width = width, score) %>%
      dplyr::group_by(bin_id) %>%
      dplyr::summarise(mean_signal = sum(score * overlap_width) / sum(overlap_width),
                       covered_bp  = sum(overlap_width),
                       .groups     = "drop")
    dplyr::as_tibble(as.data.frame(mcols(bins_chrom))) %>%
      dplyr::inner_join(covered_bins, by = "bin_id") %>%
      dplyr::mutate(covered = covered_bp / bin_width) %>%
      dplyr::select(gene_id, feature, feature_idx, bin_idx,
                    ei_group, expr_group, de_status, mean_signal, covered)
  }, error = function(e) { message("  ", chrom, ": ", conditionMessage(e)); NULL })
}

extract_profile_signal <- function(bw_gr, bins, n_workers = N_WORKERS) {
  chroms_to_run <- intersect(
    as.character(unique(seqnames(bins))),
    as.character(unique(seqnames(bw_gr)))
  )
  BiocParallel::register(BiocParallel::MulticoreParam(workers = n_workers))
  chrom_results <- BiocParallel::bplapply(chroms_to_run, function(chrom) {
    extract_profile_signal_one_chrom(chrom, bw_gr, bins)
  })
  res <- dplyr::bind_rows(chrom_results)
  if (nrow(res) == 0) stop("No bigWig signal extracted for any chromosome")
  res
}

##--------------------------------------------------------------## 11. Profile summarisation
## Two hardcoded functions with explicit group_by columns.
## Avoids dplyr scoping issues with programmatic column selection
## inside lapply.

add_x_pos <- function(df, sample_name,
                      n_bins_up   = N_BINS_UPSTREAM,
                      n_bins_body = N_BINS_BODY) {
  df %>%
    dplyr::mutate(
      sample = sample_name,
      x_pos  = dplyr::case_when(
        feature == "Upstream"   ~ as.numeric(bin_idx),
        feature == "Body"       ~ as.numeric(bin_idx) + n_bins_up,
        feature == "Downstream" ~ as.numeric(bin_idx) + n_bins_up + n_bins_body,
        feature == "Intron"     ~ as.numeric(bin_idx)
      )
    )
}

## For Plots 1-3: full stratification including ei_group
summarise_profile_ei <- function(df) {
  df %>%
    dplyr::group_by(sample, feature, ei_group,
                    expr_group, de_status, x_pos) %>%
    dplyr::summarise(
      mean    = mean(mean_signal, na.rm = TRUE),
      se      = sd(mean_signal,   na.rm = TRUE) /
        sqrt(sum(!is.na(mean_signal))),
      n_genes = sum(!is.na(mean_signal)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)
}

## For Plot 4: ei_group collapsed, full DEG n per expr bin
summarise_profile_deg <- function(df) {
  df %>%
    dplyr::group_by(sample, feature, expr_group, de_status, x_pos) %>%
    dplyr::summarise(
      mean    = mean(mean_signal, na.rm = TRUE),
      se      = sd(mean_signal,   na.rm = TRUE) /
        sqrt(sum(!is.na(mean_signal))),
      n_genes = sum(!is.na(mean_signal)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)
}

##--------------------------------------------------------------## 12. Run pipeline

raw_signal_list <- lapply(names(BW_FILES), function(sample_name) {
  cat("Extracting signal:", sample_name, "\n")
  
  bw_gr <- plyranges::read_bigwig(BW_FILES[sample_name])
  seqlevelsStyle(bw_gr) <- "UCSC"
  bw_gr <- keepSeqlevels(bw_gr,
                         intersect(seqlevels(bw_gr), STD_CHROMS),
                         pruning.mode = "coarse")
  
  cat("  Metagene bins...\n")
  meta_signal   <- extract_profile_signal(bw_gr, metagene_bins)
  cat("  Intron bins...\n")
  intron_signal <- extract_profile_signal(bw_gr, intron_bins)
  rm(bw_gr); gc()
  
  dplyr::bind_rows(meta_signal, intron_signal) %>%
    add_x_pos(sample_name)
})

## Sanity check before summarisation
cat("de_status in raw_signal_list[[1]]:\n")
print(table(raw_signal_list[[1]]$de_status, useNA = "always"))
cat("Share of each bin inside a bigWig interval (summary over bins with signal):\n")
print(round(do.call(rbind, lapply(setNames(raw_signal_list, names(BW_FILES)),
                                  function(d) summary(d$covered))), 3))

## Summarise for Plots 1-3
cat("Summarising EI profiles...\n")
profile_ei_df <- dplyr::bind_rows(
  lapply(raw_signal_list, summarise_profile_ei)
) %>%
  dplyr::rename(condition = sample) %>%
  dplyr::mutate(condition = factor(condition, levels = COND_LEVELS))

## Summarise for Plot 4
cat("Summarising DEG profiles...\n")
profile_deg_df <- dplyr::bind_rows(
  lapply(raw_signal_list, summarise_profile_deg)
) %>%
  dplyr::rename(condition = sample) %>%
  dplyr::mutate(condition = factor(condition, levels = COND_LEVELS))

rm(raw_signal_list); gc()

## Verify
cat("de_status in profile_ei_df:\n")
print(table(profile_ei_df$de_status, useNA = "always"))
cat("de_status in profile_deg_df:\n")
print(table(profile_deg_df$de_status, useNA = "always"))

##--------------------------------------------------------------## 13. Plot setup

EI_LEVELS   <- c("Low E/I", "Medium E/I", "High E/I")
EXPR_LEVELS <- c("Unexpressed", "Low", "Medium", "High")
DE_LEVELS   <- c("Up", "Down")

N_TOTAL  <- N_BINS_UPSTREAM + N_BINS_BODY + N_BINS_DOWNSTREAM

x_breaks <- c(1,  N_BINS_UPSTREAM + 1,
              N_BINS_UPSTREAM + N_BINS_BODY,
              N_TOTAL)
x_labels <- c(paste0("-", FLANK/1000, "kb"), "TSS",
              "TES", paste0("+", FLANK/1000, "kb"))

boundary_lines <- data.frame(
  xintercept = c(N_BINS_UPSTREAM + 0.5,
                 N_BINS_UPSTREAM + N_BINS_BODY + 0.5)
)

ei_theme <- theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.title     = element_blank(),
    legend.key.size  = unit(0.8, "lines"),
    axis.line        = element_line(color = "grey40"),
    axis.ticks       = element_line(color = "grey40"),
    axis.text        = element_text(color = "grey20"),
    axis.title       = element_text(color = "grey20"),
    panel.spacing    = unit(1.0, "lines"),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

SUBTITLE_BASE <- "log2(POLR2A/WT)  |  EU-seq and total RNA-seq  |  HCT116"

genes_meta_df <- dplyr::as_tibble(mcols(genes_gr)) %>%
  dplyr::mutate(
    ei_group   = factor(ei_group,   levels = EI_LEVELS),
    expr_group = factor(expr_group, levels = EXPR_LEVELS),
    de_status  = factor(de_status,  levels = c("Up", "Not DE",
                                               "Down", "Filtered"))
  )

make_labels <- function(df, group_vars) {
  df %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(label = paste0("n = ", n_genes)) %>%
    dplyr::select(all_of(c(group_vars, "label")))
}

##--------------------------------------------------------------## 14. Core plot helpers

## Save one plot as SVG in PLOT_DIR (no-op when PLOT_DIR is NULL).
## File name: metagene_<plot type>.svg
save_plot_svg <- function(p, type, width, height) {
  if (is.null(PLOT_DIR)) return(invisible(NULL))
  if (!dir.exists(PLOT_DIR)) stop("PLOT_DIR not found: ", PLOT_DIR)
  nm <- gsub("^-+|-+$", "", gsub("[^A-Za-z0-9._-]+", "-", paste("metagene", type, sep = "_")))
  f  <- file.path(PLOT_DIR, paste0(nm, ".svg"))
  svg(f, width = width, height = height)
  print(p)
  dev.off()
  message("  saved ", f)
  invisible(f)
}

make_body_plot <- function(df, label_df, facet_formula,
                           title_str, subtitle_str) {
  ggplot(df, aes(x = x_pos, y = mean,
                 color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.4) +
    geom_vline(data        = boundary_lines,
               aes(xintercept = xintercept),
               linetype    = "dashed", color = "grey50",
               linewidth   = 0.4, inherit.aes = FALSE) +
    geom_text(data        = label_df,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.1, vjust = 1.5,
              size = 3.0, color = "grey40", fontface = "italic") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = COLORS, labels = LABELS) +
    scale_fill_manual(values  = COLORS, labels = LABELS) +
    facet_grid(facet_formula, scales = "free_y") +
    labs(x = NULL, y = "log2(POLR2A/WT)",
         title = title_str, subtitle = subtitle_str) +
    ei_theme
}

make_intron_plot <- function(df, label_df, facet_formula, title_str) {
  ggplot(df, aes(x = x_pos, y = mean,
                 color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.4) +
    geom_text(data        = label_df,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.1, vjust = 1.5,
              size = 3.0, color = "grey40", fontface = "italic") +
    scale_x_continuous(breaks = c(1, N_BINS_INTRON/2, N_BINS_INTRON),
                       labels = c("5'", "50%", "3'")) +
    scale_color_manual(values = COLORS, labels = LABELS) +
    scale_fill_manual(values  = COLORS, labels = LABELS) +
    facet_grid(facet_formula, scales = "free_y") +
    labs(x = "Position along intron (5' to 3')",
         y = "log2(POLR2A/WT)",
         title = title_str,
         subtitle = "Scaled to % of intron length") +
    ei_theme
}

make_body_plot_wrap <- function(df, label_df, title_str, subtitle_str) {
  ggplot(df, aes(x = x_pos, y = mean,
                 color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.4) +
    geom_vline(data        = boundary_lines,
               aes(xintercept = xintercept),
               linetype    = "dashed", color = "grey50",
               linewidth   = 0.4, inherit.aes = FALSE) +
    geom_text(data        = label_df,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.1, vjust = 1.5,
              size = 3.0, color = "grey40", fontface = "italic") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = COLORS, labels = LABELS) +
    scale_fill_manual(values  = COLORS, labels = LABELS) +
    facet_wrap(~ expr_group, nrow = 1, scales = "free_y") +
    labs(x = NULL, y = "log2(POLR2A/WT)",
         title = title_str, subtitle = subtitle_str) +
    ei_theme
}


make_intron_plot_wrap <- function(df, label_df, title_str) {
  ggplot(df, aes(x = x_pos, y = mean,
                 color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.4) +
    geom_text(data        = label_df,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.1, vjust = 1.5,
              size = 3.0, color = "grey40", fontface = "italic") +
    scale_x_continuous(breaks = c(1, N_BINS_INTRON/2, N_BINS_INTRON),
                       labels = c("5'", "50%", "3'")) +
    scale_color_manual(values = COLORS, labels = LABELS) +
    scale_fill_manual(values  = COLORS, labels = LABELS) +
    facet_wrap(~ expr_group, nrow = 1, scales = "free_y") +
    labs(x = "Position along intron (5' to 3')",
         y = "log2(POLR2A/WT)",
         title = title_str,
         subtitle = "Scaled to % of intron length") +
    ei_theme
}

assemble <- function(p_body, p_intron, title_str, subtitle_str) {
  p_body / p_intron +
    plot_layout(heights = c(2, 1.5), guides = "collect") +
    plot_annotation(
      title    = title_str,
      subtitle = subtitle_str,
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey40", size = 11)
      )
    ) & theme(legend.position = "top")
}

##--------------------------------------------------------------## 15. Plot 1 — E/I group only (all genes)

p1_body_df <- profile_ei_df %>%
  dplyr::filter(feature != "Intron") %>%
  dplyr::mutate(ei_group = factor(ei_group, levels = EI_LEVELS)) %>%
  dplyr::group_by(condition, feature, ei_group, x_pos) %>%
  dplyr::summarise(mean = mean(mean, na.rm = TRUE),
                   se   = mean(se,   na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)

p1_intron_df <- profile_ei_df %>%
  dplyr::filter(feature == "Intron") %>%
  dplyr::mutate(ei_group = factor(ei_group, levels = EI_LEVELS)) %>%
  dplyr::group_by(condition, feature, ei_group, x_pos) %>%
  dplyr::summarise(mean = mean(mean, na.rm = TRUE),
                   se   = mean(se,   na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)

labels_p1 <- make_labels(genes_meta_df, "ei_group") %>%
  dplyr::mutate(ei_group = factor(ei_group, levels = EI_LEVELS))

p1 <- assemble(
  make_body_plot(p1_body_df, labels_p1, . ~ ei_group,
                 "Gene body — E/I group (all genes)",
                 paste0("-", FLANK/1000, "kb  |  all genes")),
  make_intron_plot(p1_intron_df, labels_p1, . ~ ei_group,
                   "Intron signal — E/I group (all genes)"),
  "Plot 1 — E/I group", SUBTITLE_BASE
)
print(p1)
save_plot_svg(p1, "plot1_ei-group", width = 12, height = 9)

##--------------------------------------------------------------## 16. Plot 2 — E/I x DEG status

p2_body_df <- profile_ei_df %>%
  dplyr::filter(feature != "Intron",
                de_status %in% c("Up", "Down")) %>%
  dplyr::mutate(ei_group  = factor(ei_group,  levels = EI_LEVELS),
                de_status = factor(de_status, levels = DE_LEVELS)) %>%
  dplyr::group_by(condition, feature, ei_group, de_status, x_pos) %>%
  dplyr::summarise(mean = mean(mean, na.rm = TRUE),
                   se   = mean(se,   na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)

p2_intron_df <- profile_ei_df %>%
  dplyr::filter(feature == "Intron",
                de_status %in% c("Up", "Down")) %>%
  dplyr::mutate(ei_group  = factor(ei_group,  levels = EI_LEVELS),
                de_status = factor(de_status, levels = DE_LEVELS)) %>%
  dplyr::group_by(condition, feature, ei_group, de_status, x_pos) %>%
  dplyr::summarise(mean = mean(mean, na.rm = TRUE),
                   se   = mean(se,   na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)

labels_p2 <- make_labels(
  genes_meta_df %>% dplyr::filter(de_status %in% c("Up", "Down")),
  c("de_status", "ei_group")
) %>%
  dplyr::mutate(ei_group  = factor(ei_group,  levels = EI_LEVELS),
                de_status = factor(de_status, levels = DE_LEVELS))

p2 <- assemble(
  make_body_plot(p2_body_df, labels_p2, de_status ~ ei_group,
                 "Gene body — DEG status x E/I group",
                 paste0("-", FLANK/1000, "kb  |  |LFC| > ",
                        LFC_THRESHOLD, ", padj < ", PADJ_THRESHOLD)),
  make_intron_plot(p2_intron_df, labels_p2, de_status ~ ei_group,
                   "Intron signal — DEG status x E/I group"),
  "Plot 2 — E/I group by DEG status", SUBTITLE_BASE
)
print(p2)
save_plot_svg(p2, "plot2_ei-group-by-deg-status", width = 12, height = 11)

##--------------------------------------------------------------## 17. Plot 3 — E/I x expression group

p3_body_df <- profile_ei_df %>%
  dplyr::filter(feature != "Intron") %>%
  dplyr::mutate(ei_group   = factor(ei_group,   levels = EI_LEVELS),
                expr_group = factor(expr_group, levels = EXPR_LEVELS)) %>%
  dplyr::group_by(condition, feature, ei_group, expr_group, x_pos) %>%
  dplyr::summarise(mean = mean(mean, na.rm = TRUE),
                   se   = mean(se,   na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)

p3_intron_df <- profile_ei_df %>%
  dplyr::filter(feature == "Intron") %>%
  dplyr::mutate(ei_group   = factor(ei_group,   levels = EI_LEVELS),
                expr_group = factor(expr_group, levels = EXPR_LEVELS)) %>%
  dplyr::group_by(condition, feature, ei_group, expr_group, x_pos) %>%
  dplyr::summarise(mean = mean(mean, na.rm = TRUE),
                   se   = mean(se,   na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)

labels_p3 <- make_labels(genes_meta_df, c("ei_group", "expr_group")) %>%
  dplyr::mutate(ei_group   = factor(ei_group,   levels = EI_LEVELS),
                expr_group = factor(expr_group, levels = EXPR_LEVELS))

p3 <- assemble(
  make_body_plot(p3_body_df, labels_p3, ei_group ~ expr_group,
                 "Gene body — E/I group x Expression",
                 paste0("-", FLANK/1000, "kb  |  all genes")),
  make_intron_plot(p3_intron_df, labels_p3, ei_group ~ expr_group,
                   "Intron signal — E/I group x Expression"),
  "Plot 3 — E/I group by expression level", SUBTITLE_BASE
)
print(p3)
save_plot_svg(p3, "plot3_ei-group-by-expression", width = 14, height = 12)

##--------------------------------------------------------------## 18. Plot 4 — DEG status x expression group

labels_p4_up <- make_labels(
  genes_meta_df %>% dplyr::filter(de_status == "Up"), "expr_group"
) %>%
  dplyr::mutate(expr_group = factor(expr_group, levels = EXPR_LEVELS))

labels_p4_down <- make_labels(
  genes_meta_df %>% dplyr::filter(de_status == "Down"), "expr_group"
) %>%
  dplyr::mutate(expr_group = factor(expr_group, levels = EXPR_LEVELS))

up_body_df <- profile_deg_df %>%
  dplyr::filter(feature != "Intron", de_status == "Up") %>%
  dplyr::mutate(expr_group = factor(expr_group, levels = EXPR_LEVELS))

up_intron_df <- profile_deg_df %>%
  dplyr::filter(feature == "Intron", de_status == "Up") %>%
  dplyr::mutate(expr_group = factor(expr_group, levels = EXPR_LEVELS))

down_body_df <- profile_deg_df %>%
  dplyr::filter(feature != "Intron", de_status == "Down") %>%
  dplyr::mutate(expr_group = factor(expr_group, levels = EXPR_LEVELS))

down_intron_df <- profile_deg_df %>%
  dplyr::filter(feature == "Intron", de_status == "Down") %>%
  dplyr::mutate(expr_group = factor(expr_group, levels = EXPR_LEVELS))

## One flat stack of four rows (nesting the Up and Down pairs gives the four
## heights to only two items and squashes the lower half); all rows equal height.
p4 <- make_body_plot_wrap(up_body_df, labels_p4_up,
                          "Upregulated genes — gene body",
                          paste0("-", FLANK/1000, "kb  |  LFC > ",
                                 LFC_THRESHOLD, ", padj < ",
                                 PADJ_THRESHOLD, "  |  HCT116")) /
  make_intron_plot_wrap(up_intron_df, labels_p4_up,
                        "Upregulated genes — intron signal") /
  make_body_plot_wrap(down_body_df, labels_p4_down,
                      "Downregulated genes — gene body",
                      paste0("-", FLANK/1000, "kb  |  LFC < -",
                             LFC_THRESHOLD, ", padj < ",
                             PADJ_THRESHOLD, "  |  HCT116")) /
  make_intron_plot_wrap(down_intron_df, labels_p4_down,
                        "Downregulated genes — intron signal") +
  plot_layout(heights = c(1, 1, 1, 1), guides = "collect") +
  plot_annotation(
    title    = "Plot 4 — DEG status by expression level",
    subtitle = SUBTITLE_BASE,
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey40", size = 11)
    )
  ) & theme(legend.position = "top")

print(p4)
save_plot_svg(p4, "plot4_deg-status-by-expression", width = 14, height = 16)

p4_up <- (
  make_body_plot_wrap(up_body_df, labels_p4_up,
                      "Upregulated genes — gene body",
                      paste0("-", FLANK/1000, "kb  |  LFC > ",
                             LFC_THRESHOLD, ", padj < ",
                             PADJ_THRESHOLD, "  |  HCT116")) /
    make_intron_plot_wrap(up_intron_df, labels_p4_up,
                          "Upregulated genes — intron signal")
) +
  plot_layout(heights = c(2, 1.5), guides = "collect") +
  plot_annotation(
    title    = "Plot 4a — Upregulated genes by expression level",
    subtitle = SUBTITLE_BASE,
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey40", size = 11)
    )
  ) & theme(legend.position = "top")

p4_down <- (
  make_body_plot_wrap(down_body_df, labels_p4_down,
                      "Downregulated genes — gene body",
                      paste0("-", FLANK/1000, "kb  |  LFC < -",
                             LFC_THRESHOLD, ", padj < ",
                             PADJ_THRESHOLD, "  |  HCT116")) /
    make_intron_plot_wrap(down_intron_df, labels_p4_down,
                          "Downregulated genes — intron signal")
) +
  plot_layout(heights = c(2, 1.5), guides = "collect") +
  plot_annotation(
    title    = "Plot 4b — Downregulated genes by expression level",
    subtitle = SUBTITLE_BASE,
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey40", size = 11)
    )
  ) & theme(legend.position = "top")

print(p4_up)
save_plot_svg(p4_up, "plot4a_upregulated-by-expression", width = 14, height = 8)
print(p4_down)
save_plot_svg(p4_down, "plot4b_downregulated-by-expression", width = 14, height = 8)

##--------------------------------------------------------------## 19. Unfaceted DEG summary

profile_nofacet_df <- profile_deg_df %>%
  dplyr::group_by(condition, feature, de_status, x_pos) %>%
  dplyr::summarise(
    mean    = mean(mean, na.rm = TRUE),
    se      = mean(se,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)

label_nofacet_up   <- dplyr::tibble(
  label = paste0("n = ", sum(genes_gr$de_status == "Up")))
label_nofacet_down <- dplyr::tibble(
  label = paste0("n = ", sum(genes_gr$de_status == "Down")))

make_body_plot_nofacet <- function(df, label_df, title_str, subtitle_str) {
  ggplot(df, aes(x = x_pos, y = mean, color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.4) +
    geom_vline(data = boundary_lines, aes(xintercept = xintercept),
               linetype = "dashed", color = "grey50",
               linewidth = 0.4, inherit.aes = FALSE) +
    geom_text(data = label_df, aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.1, vjust = 1.5,
              size = 3.5, color = "grey40", fontface = "italic") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = COLORS, labels = LABELS) +
    scale_fill_manual(values  = COLORS, labels = LABELS) +
    labs(x = NULL, y = "log2(POLR2A/WT)",
         title = title_str, subtitle = subtitle_str) +
    ei_theme
}

make_intron_plot_nofacet <- function(df, label_df, title_str) {
  ggplot(df, aes(x = x_pos, y = mean, color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.4) +
    geom_text(data = label_df, aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.1, vjust = 1.5,
              size = 3.5, color = "grey40", fontface = "italic") +
    scale_x_continuous(breaks = c(1, N_BINS_INTRON/2, N_BINS_INTRON),
                       labels = c("5'", "50%", "3'")) +
    scale_color_manual(values = COLORS, labels = LABELS) +
    scale_fill_manual(values  = COLORS, labels = LABELS) +
    labs(x = "Position along intron (5' to 3')", y = "log2(POLR2A/WT)",
         title = title_str, subtitle = "Scaled to % of intron length") +
    ei_theme
}

p_body_up_nf <- make_body_plot_nofacet(
  profile_nofacet_df %>%
    dplyr::filter(feature != "Intron", de_status == "Up"),
  label_nofacet_up,
  "Upregulated genes — gene body",
  paste0("-", FLANK/1000, "kb  |  LFC > ", LFC_THRESHOLD,
         ", padj < ", PADJ_THRESHOLD, "  |  HCT116")
)

p_intron_up_nf <- make_intron_plot_nofacet(
  profile_nofacet_df %>%
    dplyr::filter(feature == "Intron", de_status == "Up"),
  label_nofacet_up, "Upregulated genes — intron signal"
)

p_body_down_nf <- make_body_plot_nofacet(
  profile_nofacet_df %>%
    dplyr::filter(feature != "Intron", de_status == "Down"),
  label_nofacet_down,
  "Downregulated genes — gene body",
  paste0("-", FLANK/1000, "kb  |  LFC < -", LFC_THRESHOLD,
         ", padj < ", PADJ_THRESHOLD, "  |  HCT116")
)

p_intron_down_nf <- make_intron_plot_nofacet(
  profile_nofacet_df %>%
    dplyr::filter(feature == "Intron", de_status == "Down"),
  label_nofacet_down, "Downregulated genes — intron signal"
)

p_combined_nf <- (p_body_up_nf | p_intron_up_nf) /
  (p_body_down_nf | p_intron_down_nf) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "POLR2A degradation — DEG log2 ratio profiles",
    subtitle = SUBTITLE_BASE,
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey40", size = 11)
    )
  ) & theme(legend.position = "top")

print(p_combined_nf)
save_plot_svg(p_combined_nf, "deg-profiles-unfaceted", width = 12, height = 9)

# ── Unfaceted expression group summary ───────────────────────────────────────
# Summarise collapsing de_status and ei_group
profile_nofacet_expr_df <- profile_ei_df %>%
  dplyr::group_by(condition, feature, expr_group, x_pos) %>%
  dplyr::summarise(
    mean    = mean(mean, na.rm = TRUE),
    se      = mean(se,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    upper      = mean + 1.96*se,
    lower      = mean - 1.96*se,
    expr_group = factor(expr_group, levels = EXPR_LEVELS)
  )

labels_expr_nf <- make_labels(genes_meta_df, "expr_group") %>%
  dplyr::mutate(expr_group = factor(expr_group, levels = EXPR_LEVELS))

p_expr_body <- make_body_plot_wrap(
  profile_nofacet_expr_df %>% dplyr::filter(feature != "Intron"),
  labels_expr_nf,
  "All genes — gene body by expression group",
  paste0("-", FLANK/1000, "kb  |  all genes  |  HCT116")
)

p_expr_intron <- make_intron_plot_wrap(
  profile_nofacet_expr_df %>% dplyr::filter(feature == "Intron"),
  labels_expr_nf,
  "All genes — intron signal by expression group"
)

p_expr_nf <- assemble(
  p_expr_body, p_expr_intron,
  "Expression group — gene body and intron signal",
  SUBTITLE_BASE
)

print(p_expr_nf)
save_plot_svg(p_expr_nf, "expression-group-unfaceted", width = 14, height = 8)

##--------------------------------------------------------------## 20. Figure 3C panel
## Upregulated genes with Unexpressed or Low expression (the two groups pooled, this
## plot only) and Downregulated genes with High expression: gene-body and intron
## profiles, EU-seq vs total RNA-seq. The pooled mean and standard error are computed
## exactly from the per-group mean, SE and n of profile_deg_df (pooled variance =
## within-group + between-group).

UP_GROUPS   <- c("Unexpressed", "Low")
DOWN_GROUPS <- "High"

pool_expr_groups <- function(df) {
  df %>%
    dplyr::filter(n_genes > 0) %>%
    dplyr::mutate(var_i = tidyr::replace_na(se, 0)^2 * n_genes) %>%
    dplyr::group_by(condition, feature, x_pos) %>%
    dplyr::summarise(
      N       = sum(n_genes),
      M       = sum(mean * n_genes) / sum(n_genes),
      V       = (sum((n_genes - 1) * var_i) + sum(n_genes * (mean - M)^2)) / (N - 1),
      .groups = "drop"
    ) %>%
    dplyr::mutate(mean = M, se = sqrt(V / N),
                  upper = mean + 1.96*se, lower = mean - 1.96*se)
}

fig3c_up_df <- profile_deg_df %>%
  dplyr::filter(de_status == "Up", expr_group %in% UP_GROUPS) %>%
  pool_expr_groups()
fig3c_down_df <- profile_deg_df %>%
  dplyr::filter(de_status == "Down", expr_group %in% DOWN_GROUPS)

fig3c_n_up   <- sum(genes_meta_df$de_status == "Up"   & genes_meta_df$expr_group %in% UP_GROUPS,   na.rm = TRUE)
fig3c_n_down <- sum(genes_meta_df$de_status == "Down" & genes_meta_df$expr_group %in% DOWN_GROUPS, na.rm = TRUE)
cat("Figure 3C genes: Up (", paste(UP_GROUPS, collapse = " + "), ") n = ", fig3c_n_up,
    " | Down (", paste(DOWN_GROUPS, collapse = " + "), ") n = ", fig3c_n_down, "\n", sep = "")

make_fig3c_panel <- function(df, feature_type, panel_title, n_label, n_at_top) {
  d <- if (feature_type == "body") dplyr::filter(df, feature != "Intron") else dplyr::filter(df, feature == "Intron")
  p <- ggplot(d, aes(x = x_pos, y = mean, color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 1.1) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.5)
  if (feature_type == "body") {
    p <- p +
      geom_vline(data = boundary_lines, aes(xintercept = xintercept),
                 linetype = "dashed", color = "grey50", linewidth = 0.5) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels)
  } else {
    p <- p + scale_x_continuous(breaks = c(1, N_BINS_INTRON/2, N_BINS_INTRON),
                                labels = c("5'", "50%", "3'"))
  }
  p +
    annotate("text", x = Inf, y = if (n_at_top) Inf else -Inf, label = n_label,
             hjust = 1.1, vjust = if (n_at_top) 1.5 else -0.6,
             size = 4, color = "grey40", fontface = "italic") +
    scale_color_manual(values = COLORS, labels = LABELS) +
    scale_fill_manual(values  = COLORS, labels = LABELS) +
    labs(x = NULL, y = NULL, title = panel_title) +
    theme_classic(base_size = 14) +
    theme(plot.title      = element_text(hjust = 0.5, size = 12),
          legend.title    = element_blank(),
          legend.position = "bottom",
          axis.line       = element_line(color = "grey20"),
          axis.ticks      = element_line(color = "grey20"))
}

fig3c_text <- function(label, angle = 0, face = "plain", size = 5) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = label, angle = angle, fontface = face, size = size) +
    theme_void()
}

fig3c_design <- c(
  patchwork::area(1, 2, 1, 2), patchwork::area(1, 3, 1, 3),   # column titles
  patchwork::area(2, 2, 2, 2), patchwork::area(2, 3, 2, 3),   # gene body
  patchwork::area(3, 2, 3, 2), patchwork::area(3, 3, 3, 3),   # intron
  patchwork::area(2, 1, 3, 1),                                # y label
  patchwork::area(4, 2, 4, 3)                                 # x label
)

p_fig3c <- patchwork::wrap_plots(
  list(
    fig3c_text("Upregulated-Low Expression",   face = "bold", size = 4.5),
    fig3c_text("Downregulated-High Expression", face = "bold", size = 4.5),
    make_fig3c_panel(fig3c_up_df,   "body",   "Gene Body Signal", paste0("n = ", fig3c_n_up),   FALSE),
    make_fig3c_panel(fig3c_down_df, "body",   "Gene Body Signal", paste0("n = ", fig3c_n_down), TRUE),
    make_fig3c_panel(fig3c_up_df,   "intron", "Intronic Signal",  paste0("n = ", fig3c_n_up),   FALSE),
    make_fig3c_panel(fig3c_down_df, "intron", "Intronic Signal",  paste0("n = ", fig3c_n_down), TRUE),
    fig3c_text("log2(POLR2A/WT)", angle = 90, size = 5.5),
    fig3c_text("Position along intron (5' to 3')", size = 5.5)
  ),
  design = fig3c_design, widths = c(0.07, 1, 1), heights = c(0.12, 1, 1, 0.08),
  guides = "collect"
) & theme(legend.position = "bottom")

print(p_fig3c)
save_plot_svg(p_fig3c, "fig3c_up-unexpressed-low_down-high", width = 9, height = 8)

length(genes_gr)
table(genes_meta_df$expr_group, useNA = "always")
# Option: restrict to protein-coding only
genes_gr <- genes_gr[genes_gr$gene_biotype == "protein_coding"]