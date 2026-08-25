##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Gene-body coverage profiles over DEGs (ChIP-seq). Computes
# normalized exon/intron gene-body coverage profiles (TSS -> TES,
# 20 bins) over differentially expressed genes, from ChIP-seq
# bigWig signal (H3K27ac, H3K9me3, POL2 +/- POLR2A-KO; public
# ENCODE data, see Table 5).
#
# Split from the original MAIN_chip_bw_compute_coverage.R, which
# toggled between ChIP-seq and CUT&Tag data via which BW_FILES
# block was left active -- now a fixed, single-purpose script. See
# CT_compute_coverage.R for the CUT&Tag counterpart (same pipeline
# structure, deliberately duplicated rather than factored into
# shared code).
#
# Status: probable exploratory code; no confirmed figure match yet.
##################################################################

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
library(cowplot)
library(scales)

# ── Parameters ────────────────────────────────────────────────────────────────
GTF              <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
EXPR_FILE        <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/HCT116_RNAseq_quantification.tsv"
DESEQ_FILE       <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"

PATH.BW          <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/chip/"
BW_FILES         <- c(H3K27ac    = file.path(PATH.BW, "h3k27ac.bw"),
                      H3K9me3    = file.path(PATH.BW, "h3k9me3.bw"),
                      POL2       = file.path(PATH.BW, "pol2.bw"),
                      H3K27ac_ko = file.path(PATH.BW, "h3k27ac_ko.bw"),
                      H3K9me3_ko = file.path(PATH.BW, "h3k9me3_ko.bw"))
DATA_SOURCE      <- "ChIP-seq"


LFC_THRESHOLD    <- 1
PADJ_THRESHOLD   <- 0.05
TPM_THRESHOLD    <- 1
STD_CHROMS       <- paste0("chr", c(1:22, "X", "Y", "M"))
N_WORKERS        <- 7
N_BINS_TSS       <- 40
N_BINS_TES       <- 40
N_BINS_EXON      <- 50
N_BINS_INTRON    <- 100
FLANK            <- 2000
MIN_WIDTH_EXON   <- 50
MIN_WIDTH_INTRON <- 200

# ── 1. Build TxDb ─────────────────────────────────────────────────────────────
txdb <- makeTxDbFromGFF(GTF, format = "gtf")

# ── 2. Build genes GRanges ────────────────────────────────────────────────────
genes_gr         <- genes(txdb)
genes_gr$gene_id <- sub("\\..*", "", names(genes_gr))
seqlevelsStyle(genes_gr) <- "UCSC"
genes_gr <- keepSeqlevels(genes_gr,
                          intersect(seqlevels(genes_gr), STD_CHROMS),
                          pruning.mode = "coarse")

# ── 3. Build per-gene exon GRangesList ────────────────────────────────────────
exons_by_gene_gr <- exonsBy(txdb, by = "gene")
names(exons_by_gene_gr) <- sub("\\..*", "", names(exons_by_gene_gr))
seqlevelsStyle(exons_by_gene_gr) <- "UCSC"
exons_by_gene_gr <- keepSeqlevels(exons_by_gene_gr,
                                  intersect(seqlevels(exons_by_gene_gr),
                                            STD_CHROMS),
                                  pruning.mode = "coarse")

# ── 4. Build per-gene intron GRangesList ──────────────────────────────────────
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

# ── 5. Compute E/I ratio per gene ─────────────────────────────────────────────
exon_bp      <- sum(width(reduce(exons_by_gene_gr)))
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

# ── 5. Load expression data ───────────────────────────────────────────────────
expr <- read_tsv(EXPR_FILE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "",
                              sapply(strsplit(target_id, "\\|"), `[`, 2))) %>%
  dplyr::select(gene_id, tpm) %>%
  dplyr::rename(TPM = tpm) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

# ── 6. Load DESeq2 results ────────────────────────────────────────────────────
deseq_res <- read_csv(DESEQ_FILE) %>%
  dplyr::mutate(
    gene_id   = sub("\\..*", "", gene_id),
    de_status = dplyr::case_when(
      is.na(padj)                                              ~ "Filtered",
      padj >= PADJ_THRESHOLD                                   ~ "Not DE",
      padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD ~ "Down",
      TRUE                                                     ~ "Not DE"
    ),
    de_status = factor(de_status,
                       levels = c("Up", "Not DE", "Down", "Filtered"))
  )

cat("DE status counts:\n")
print(table(deseq_res$de_status))

# ── 7. Add E/I and DE groups to genes_gr ──────────────────────────────────────
de_ei_joined <- dplyr::tibble(gene_id = genes_gr$gene_id) %>%
  dplyr::left_join(
    deseq_res %>% dplyr::select(gene_id, log2FoldChange, padj, de_status),
    by = "gene_id"
  ) %>%
  dplyr::left_join(
    exon_bp_df %>% dplyr::select(gene_id, ei_ratio, ei_group),
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    de_status = dplyr::if_else(
      is.na(de_status),
      factor("Not DE", levels = levels(deseq_res$de_status)),
      de_status
    )
  )

genes_gr$de_status <- de_ei_joined$de_status
genes_gr$ei_group  <- de_ei_joined$ei_group

cat("DE status counts:\n")
print(table(genes_gr$de_status, useNA = "always"))
cat("E/I group counts:\n")
print(table(genes_gr$ei_group, useNA = "always"))

# ── 8. Build bins ─────────────────────────────────────────────────────────────
make_tss_bins <- function(genes_gr, flank = FLANK, n_bins = N_BINS_TSS) {
  tss_gr      <- GenomicRanges::resize(genes_gr, width = 1, fix = "start")
  tss_flanked <- GenomicRanges::flank(tss_gr, width = flank, both = TRUE)
  tss_flanked <- keepSeqlevels(tss_flanked,
                               intersect(seqlevels(tss_flanked), STD_CHROMS),
                               pruning.mode = "coarse")
  tiled       <- tile(tss_flanked, n = n_bins)
  flat        <- unlist(tiled, use.names = FALSE)
  n_features  <- length(genes_gr)
  bin_idx     <- rep.int(seq_len(n_bins), n_features)
  feature_idx <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus    <- as.logical(strand(genes_gr) == "-")
  minus_mask  <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "TSS",
    gene_id     = rep.int(genes_gr$gene_id, rep.int(n_bins, n_features))
  )
  flat
}

make_tes_bins <- function(genes_gr, flank = FLANK, n_bins = N_BINS_TES) {
  tes_gr      <- GenomicRanges::resize(genes_gr, width = 1, fix = "end")
  tes_flanked <- GenomicRanges::flank(tes_gr, width = flank, both = TRUE)
  tes_flanked <- keepSeqlevels(tes_flanked,
                               intersect(seqlevels(tes_flanked), STD_CHROMS),
                               pruning.mode = "coarse")
  tiled       <- tile(tes_flanked, n = n_bins)
  flat        <- unlist(tiled, use.names = FALSE)
  n_features  <- length(genes_gr)
  bin_idx     <- rep.int(seq_len(n_bins), n_features)
  feature_idx <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus    <- as.logical(strand(genes_gr) == "-")
  minus_mask  <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "TES",
    gene_id     = rep.int(genes_gr$gene_id, rep.int(n_bins, n_features))
  )
  flat
}

make_exon_bins <- function(exons_by_gene, n_bins = N_BINS_EXON,
                           min_width = MIN_WIDTH_EXON) {
  exons_reduced      <- reduce(exons_by_gene)
  flat_exons         <- unlist(exons_reduced, use.names = TRUE)
  flat_exons$gene_id <- sub("\\..*", "", names(flat_exons))
  flat_exons         <- flat_exons[width(flat_exons) >= min_width]
  seqlevelsStyle(flat_exons) <- "UCSC"
  flat_exons <- keepSeqlevels(flat_exons,
                              intersect(seqlevels(flat_exons), STD_CHROMS),
                              pruning.mode = "coarse")
  too_narrow <- width(flat_exons) < n_bins
  if (any(too_narrow)) {
    message(sum(too_narrow), " exons removed: narrower than ", n_bins, " bins")
    flat_exons <- flat_exons[!too_narrow]
  }
  tiled         <- tile(flat_exons, n = n_bins)
  flat          <- unlist(tiled, use.names = FALSE)
  n_features    <- length(flat_exons)
  bin_idx       <- rep.int(seq_len(n_bins), n_features)
  feature_idx   <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus      <- as.logical(strand(flat_exons) == "-")
  minus_mask    <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "Exon",
    gene_id     = rep.int(flat_exons$gene_id, rep.int(n_bins, n_features))
  )
  flat
}

make_intron_bins <- function(introns_by_gene, n_bins = N_BINS_INTRON,
                             min_width = MIN_WIDTH_INTRON) {
  introns_reduced      <- reduce(introns_by_gene)
  flat_introns         <- unlist(introns_reduced, use.names = TRUE)
  flat_introns$gene_id <- sub("\\..*", "", names(flat_introns))
  flat_introns         <- flat_introns[width(flat_introns) >= min_width]
  seqlevelsStyle(flat_introns) <- "UCSC"
  flat_introns <- keepSeqlevels(flat_introns,
                                intersect(seqlevels(flat_introns), STD_CHROMS),
                                pruning.mode = "coarse")
  too_narrow <- width(flat_introns) < n_bins
  if (any(too_narrow)) {
    message(sum(too_narrow), " introns removed: narrower than ",
            n_bins, " bins")
    flat_introns <- flat_introns[!too_narrow]
  }
  tiled         <- tile(flat_introns, n = n_bins)
  flat          <- unlist(tiled, use.names = FALSE)
  n_features    <- length(flat_introns)
  bin_idx       <- rep.int(seq_len(n_bins), n_features)
  feature_idx   <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus      <- as.logical(strand(flat_introns) == "-")
  minus_mask    <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "Intron",
    gene_id     = rep.int(flat_introns$gene_id, rep.int(n_bins, n_features))
  )
  flat
}

cat("Building TSS bins...\n")
tss_bins    <- make_tss_bins(genes_gr)
cat("Building TES bins...\n")
tes_bins    <- make_tes_bins(genes_gr)
cat("Building exon bins...\n")
exon_bins   <- make_exon_bins(exons_by_gene_gr)
cat("Building intron bins...\n")
intron_bins <- make_intron_bins(introns_by_gene_gr)

body_bins <- c(tss_bins, exon_bins, tes_bins)
rm(tss_bins, tes_bins, exon_bins); gc()

cat("Body bins size:",   format(object.size(body_bins),   units = "MB"), "\n")
cat("Intron bins size:", format(object.size(intron_bins), units = "MB"), "\n")

# ── 9. Extract signal per chromosome ─────────────────────────────────────────
extract_profile_signal_one_chrom <- function(chrom, bw_gr, bins) {
  bw_chrom   <- bw_gr[as.character(seqnames(bw_gr))  == chrom]
  bins_chrom <- bins[as.character(seqnames(bins))     == chrom]
  if (length(bw_chrom) == 0 || length(bins_chrom) == 0) return(NULL)
  tryCatch({
    plyranges::join_overlap_inner(bins_chrom, bw_chrom) %>%
      dplyr::as_tibble() %>%
      dplyr::rename(overlap_width = width) %>%
      dplyr::group_by(gene_id, feature, feature_idx, bin_idx) %>%
      dplyr::summarise(
        mean_signal = sum(score * overlap_width) / sum(overlap_width),
        .groups = "drop"
      )
  }, error = function(e) NULL)
}

extract_profile_signal <- function(bw_gr, bins, n_workers = N_WORKERS) {
  chroms_to_run <- intersect(
    as.character(unique(seqnames(bins))),
    as.character(unique(seqnames(bw_gr)))
  )
  BiocParallel::register(BiocParallel::MulticoreParam(workers = n_workers))
  chrom_results <- BiocParallel::bplapply(chroms_to_run, function(chrom) {
    extract_profile_signal_one_chrom(chrom = chrom, bw_gr = bw_gr, bins = bins)
  })
  dplyr::bind_rows(chrom_results)
}

# ── 10. Aggregate ─────────────────────────────────────────────────────────────
compute_profile <- function(signal_df, genes_gr, sample_name,
                            n_bins_tss  = N_BINS_TSS,
                            n_bins_body = N_BINS_EXON,
                            n_bins_tes  = N_BINS_TES) {
  
  genes_meta <- dplyr::as_tibble(mcols(genes_gr)) %>%
    dplyr::select(gene_id, ei_group, de_status)
  
  signal_df %>%
    dplyr::mutate(sample = sample_name) %>%
    dplyr::left_join(genes_meta, by = "gene_id") %>%
    dplyr::filter(!is.na(ei_group),
                  !is.na(de_status),
                  de_status %in% c("Up", "Down")) %>%
    dplyr::mutate(
      x_pos = dplyr::case_when(
        .data$feature == "TSS"    ~ as.numeric(.data$bin_idx),
        .data$feature == "Exon"   ~ as.numeric(.data$bin_idx) + n_bins_tss,
        .data$feature == "TES"    ~ as.numeric(.data$bin_idx) + n_bins_tss +
          n_bins_body,
        .data$feature == "Intron" ~ as.numeric(.data$bin_idx)
      )
    ) %>%
    dplyr::group_by(sample, .data$feature, ei_group, de_status, x_pos) %>%
    dplyr::summarise(
      mean  = mean(mean_signal, na.rm = TRUE),
      se    = sd(mean_signal,   na.rm = TRUE) /
        sqrt(sum(!is.na(mean_signal))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      upper = mean + 1.96 * se,
      lower = mean - 1.96 * se
    )
}

# ── 11. Run pipeline ──────────────────────────────────────────────────────────
profile_results <- lapply(names(BW_FILES), function(sample_name,
                                                    body_bins,
                                                    intron_bins,
                                                    genes_gr,
                                                    n_bins_tss,
                                                    n_bins_body,
                                                    n_bins_tes,
                                                    n_workers) {
  cat("Processing:", sample_name, "\n")
  
  bw_gr <- plyranges::read_bigwig(BW_FILES[sample_name])
  seqlevelsStyle(bw_gr) <- "UCSC"
  bw_gr <- keepSeqlevels(bw_gr,
                         intersect(seqlevels(bw_gr), STD_CHROMS),
                         pruning.mode = "coarse")
  
  cat("  Extracting body signal...\n")
  body_signal   <- extract_profile_signal(bw_gr, body_bins,
                                          n_workers = n_workers)
  cat("  Extracting intron signal...\n")
  intron_signal <- extract_profile_signal(bw_gr, intron_bins,
                                          n_workers = n_workers)
  rm(bw_gr); gc()
  
  all_signal <- dplyr::bind_rows(body_signal, intron_signal)
  rm(body_signal, intron_signal); gc()
  
  cat("  Computing profile...\n")
  compute_profile(all_signal, genes_gr, sample_name,
                  n_bins_tss  = n_bins_tss,
                  n_bins_body = n_bins_body,
                  n_bins_tes  = n_bins_tes)
  
}, body_bins   = body_bins,
intron_bins = intron_bins,
genes_gr    = genes_gr,
n_bins_tss  = N_BINS_TSS,
n_bins_body = N_BINS_EXON,
n_bins_tes  = N_BINS_TES,
n_workers   = N_WORKERS)

profile_df <- dplyr::bind_rows(profile_results)

# ── 12. Shared plot elements ──────────────────────────────────────────────────

# Tol's "Muted" palette ## LETS USE THIS ONE
sample_colors <- c(
  H3K27ac    = "#CC6677",  ## rose
  H3K4me3 = "#CC6677",  ## rose
  H3K27ac_ko = "#DDCC77",  ## sand
  H3K9me3    = "#332288",  ## indigo
  H3K9me3_ko = "#88CCEE",  ## cyan
  POL2       = "#44AA99"   ## teal
)


x_breaks <- c(1,
              N_BINS_TSS / 2,
              N_BINS_TSS + 1,
              N_BINS_TSS + N_BINS_EXON,
              N_BINS_TSS + N_BINS_EXON + N_BINS_TES / 2,
              N_BINS_TSS + N_BINS_EXON + N_BINS_TES)
x_labels <- c(paste0("-", FLANK / 1000, "kb"),
              "TSS",
              "5'", "3'",
              "TES",
              paste0("+", FLANK / 1000, "kb"))

boundary_lines <- data.frame(
  xintercept = c(N_BINS_TSS / 2,                              # TSS midpoint
                 N_BINS_TSS + N_BINS_EXON + N_BINS_TES / 2), # TES midpoint
  label      = c("TSS", "TES")
)

body_data   <- profile_df %>% dplyr::filter(.data$feature != "Intron")
intron_data <- profile_df %>% dplyr::filter(.data$feature == "Intron")

cat("Body data rows:",   nrow(body_data),   "\n")
cat("Intron data rows:", nrow(intron_data), "\n")

# ── 13. Label data ────────────────────────────────────────────────────────────
label_df <- dplyr::as_tibble(mcols(genes_gr)) %>%
  dplyr::filter(de_status %in% c("Up", "Down"),
                !is.na(ei_group)) %>%
  dplyr::group_by(ei_group, de_status) %>%
  dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n_genes))

# ── 14. Plot gene body ────────────────────────────────────────────────────────


p_body <- ggplot(body_data,
                 aes(x = x_pos, y = mean, color = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(data      = boundary_lines,
             aes(xintercept = xintercept),
             linetype  = "dashed",
             color     = "grey50",
             linewidth = 0.4) +
  geom_text(data        = label_df,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3, color = "grey30") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values  = sample_colors) +
  facet_grid(ei_group ~ de_status) +
  labs(x        = NULL,
       y        = paste0("Mean ", DATA_SOURCE, " signal"),
       color    = "Sample", fill = "Sample",
       title    = paste0("Gene body ", DATA_SOURCE, " profile"),
       subtitle = paste0("TSS ±", FLANK / 1000, "kb — scaled exon body — TES ±",
                         FLANK / 1000, "kb — stratified by E/I ratio × DE status")) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "top")

# Set the same y range on both plots
#Y_MIN <- 0
#Y_MAX <- 0.4  # adjust to your signal range
#p_body <- p_body + coord_cartesian(ylim = c(Y_MIN, Y_MAX))

print(p_body)

# ── 15. Plot introns ──────────────────────────────────────────────────────────

p_intron <- ggplot(intron_data,
                   aes(x = x_pos, y = mean, color = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_text(data        = label_df,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3, color = "grey30") +
  scale_x_continuous(breaks = c(1, N_BINS_INTRON / 2, N_BINS_INTRON),
                     labels = c("5'", "50%", "3'")) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values  = sample_colors) +
  facet_grid(ei_group ~ de_status) +
  labs(x        = "Position (% of intron, 5'→3')",
       y        = paste0("Mean ", DATA_SOURCE, " signal"),
       color    = "Sample", fill = "Sample",
       title    = paste0("Intron ", DATA_SOURCE, " profile"),
       subtitle = "Scaled to % of intron length — stratified by E/I ratio × DE status") +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "top")

print(p_intron)

###### -- Just print DE status ──────────────────────────────────────────────────────────

# Reaggregate from existing profile_df — no need to rerun extraction
profile_df_de <- profile_df %>%
  dplyr::group_by(sample, .data$feature, de_status, x_pos) %>%
  dplyr::summarise(
    mean  = mean(mean, na.rm = TRUE),
    se    = mean(se,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    upper = mean + 1.96 * se,
    lower = mean - 1.96 * se
  )

# Label data
label_de <- dplyr::as_tibble(mcols(genes_gr)) %>%
  dplyr::filter(de_status %in% c("Up", "Down")) %>%
  dplyr::group_by(de_status) %>%
  dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n_genes))

body_data_de   <- profile_df_de %>% dplyr::filter(.data$feature != "Intron")
intron_data_de <- profile_df_de %>% dplyr::filter(.data$feature == "Intron")

unique(body_data_de$sample)
#body_data_de<- body_data_de %>% filter(body_data_de$sample == "POL2" | body_data_de$sample == "POL2_ko")


p_body_de <- ggplot(body_data_de,
                    aes(x = x_pos, y = mean, color = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(data      = boundary_lines,
             aes(xintercept = xintercept),
             linetype  = "dashed", color = "grey50", linewidth = 0.4) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3, color = "grey30") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values  = sample_colors) +
  facet_wrap(~ de_status, nrow = 1) +
  labs(x = NULL, y = paste0("Mean ", DATA_SOURCE, " signal"),
       color = "Sample", fill = "Sample",
       title = paste0("Gene body ", DATA_SOURCE, " profile")) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "top")

print(p_body_de)

unique(intron_data_de$sample)
#intron_data_de<- intron_data_de %>% filter(intron_data_de$sample == "H3K9me3" | intron_data_de$sample == "H3K9me3_ko")

p_intron_de <- ggplot(intron_data_de,
                      aes(x = x_pos, y = mean, color = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3, color = "grey30") +
  scale_x_continuous(breaks = c(1, N_BINS_INTRON / 2, N_BINS_INTRON),
                     labels = c("5'", "50%", "3'")) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values  = sample_colors) +
  facet_wrap(~ de_status, nrow = 1) +
  labs(x = "Position (% of intron, 5'→3')",
       y = paste0("Mean ", DATA_SOURCE, " signal"),
       color = "Sample", fill = "Sample",
       title = paste0("Intron ", DATA_SOURCE, " profile")) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "top")

print(p_intron_de)

#########################################################################################
# PRETTY PLOT
#########################################################################################


library(RColorBrewer)

# ── Color palette — divergent, colorblind-friendly ────────────────────────────
# Using a manually curated palette that works well for ChIP marks
sample_colors <- c(
  H3K27ac    = "#E69F00",  ## orange
  H3K9me3    = "#0072B2",  ## blue
  POL2       = "#009E73",  ## bluish-green
  H3K27ac_ko = "#CC79A7",  ## pink/mauve
  H3K9me3_ko = "#D55E00"   ## vermillion
)

# Okabe-Ito 
sample_colors <- c(
  H3K27ac    = "#E69F00",  ## amber
  H3K27ac_ko = "#F0E442",  ## yellow (lighter partner)
  H3K9me3    = "#0072B2",  ## deep blue
  H3K9me3_ko = "#56B4E9",  ## sky blue (lighter partner)
  POL2       = "#009E73"   ## teal
)

# Paul Tol's "Bright" palette
sample_colors <- c(
  H3K27ac    = "#4477AA",  ## blue
  H3K27ac_ko = "#66CCEE",  ## cyan
  H3K9me3    = "#AA3377",  ## purple
  H3K9me3_ko = "#EE6677",  ## rose
  POL2       = "#228833"   ## green
)

# Tol's "Muted" palette ## LETS USE THIS ONE
sample_colors <- c(
  H3K27ac    = "#CC6677",  ## rose
  H3K4me3 = "#CC6677",  ## rose
  H3K27ac_ko = "#DDCC77",  ## sand
  H3K9me3    = "#332288",  ## indigo
  H3K9me3_ko = "#88CCEE",  ## cyan
  POL2       = "#44AA99"   ## teal
)

# Paired high-contrast
sample_colors <- c(
  H3K27ac    = "#D95F02",  ## dark orange
  H3K27ac_ko = "#FDB863",  ## light orange
  H3K9me3    = "#1F78B4",  ## dark blue
  H3K9me3_ko = "#A6CEE3",  ## light blue
  POL2       = "#33A02C"   ## dark green
)

# ── Fix DE status labels ───────────────────────────────────────────────────────

body_data_de   <- profile_df_de %>% dplyr::filter(.data$feature != "Intron")
intron_data_de <- profile_df_de %>% dplyr::filter(.data$feature == "Intron")

body_data_de <- body_data_de %>%
  dplyr::mutate(
    de_status = dplyr::recode(de_status,
                              "Up"   = "Upregulated",
                              "Down" = "Downregulated"),
    de_status = factor(de_status, levels = c("Upregulated", "Downregulated"))
  )

intron_data_de <- intron_data_de %>%
  dplyr::mutate(
    de_status = dplyr::recode(de_status,
                              "Up"   = "Upregulated",
                              "Down" = "Downregulated"),
    de_status = factor(de_status, levels = c("Upregulated", "Downregulated"))
  )

label_de <- label_de %>%
  dplyr::mutate(
    de_status = dplyr::recode(de_status,
                              "Up"   = "Upregulated",
                              "Down" = "Downregulated"),
    de_status = factor(de_status, levels = c("Upregulated", "Downregulated"))
  )

# ── Shared theme ──────────────────────────────────────────────────────────────
chip_theme <- theme_classic(base_size = 13) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 13),
    legend.position   = "top",
    legend.title      = element_blank(),
    legend.key.size   = unit(0.8, "lines"),
    axis.line         = element_line(color = "grey40"),
    axis.ticks        = element_line(color = "grey40"),
    axis.text         = element_text(color = "grey20"),
    axis.title        = element_text(color = "grey20"),
    panel.spacing     = unit(1.2, "lines"),
    plot.title        = element_text(face = "bold", size = 14),
    plot.subtitle     = element_text(color = "grey40", size = 11)
  )

# ── Gene body plot ────────────────────────────────────────────────────────────

unique(body_data_de$sample)
#body_data_de<- body_data_de %>% filter(body_data_de$sample %in% c("H3K27ac", "POL2"))
#body_data_de<- body_data_de %>% filter(body_data_de$sample %in% c("H3K9me3", "H3K9me3_ko"))

p_body_de <- ggplot(body_data_de,
                    aes(x = x_pos, y = mean, color = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(data      = boundary_lines,
             aes(xintercept = xintercept),
             linetype  = "dashed",
             color     = "grey60",
             linewidth = 0.4) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40",
            fontface = "italic") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values  = sample_colors) +
  facet_wrap(~ de_status, nrow = 1) +
  facet_grid(ei_group ~ de_status) +
  labs(x        = NULL,
       y        = paste0("Mean ", DATA_SOURCE, " signal"),
       title    = paste0("Gene body ", DATA_SOURCE, " profile"),
       subtitle = paste0("TSS ±", FLANK / 1000, "kb  |  scaled exon body  |  TES ±",
                         FLANK / 1000, "kb")) +
  chip_theme

print(p_body_de)

# ── Intron plot ───────────────────────────────────────────────────────────────

unique(intron_data_de$sample)
intron_data_de<- intron_data_de %>% filter(intron_data_de$sample %in% c("H3K27ac", "POL2"))
intron_data_de<- intron_data_de %>% filter(intron_data_de$sample %in% c("H3K9me3", "H3K9me3_ko"))


p_intron_de <- ggplot(intron_data_de,
                      aes(x = x_pos, y = mean, color = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40",
            fontface = "italic") +
  scale_x_continuous(breaks = c(1, N_BINS_INTRON / 2, N_BINS_INTRON),
                     labels = c("5'", "50%", "3'")) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values  = sample_colors) +
  facet_wrap(~ de_status, nrow = 1) +
  labs(x        = "Position along intron (5'→3')",
       y        = paste0("Mean ", DATA_SOURCE, " signal"),
       title    = paste0("Intron ", DATA_SOURCE, " profile"),
       subtitle = "Signal scaled to percentage of intron length") +
  chip_theme


print(p_intron_de)

#########################################################################################
# COMPUTE PEAK SUM FOR DEGs VERSUS TPM
#########################################################################################

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
library(cowplot)
library(scales)

# ── 0. Parameters ──────────────────────────────────────────────────────────────
GTF        <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
EXPR_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/HCT116_RNAseq_quantification.tsv"

# ── FOR CHIP ──────────────────────────────────────────────────────────────

PATH.BW          <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/chip/"
BW_FILES         <- c(H3K27ac    = file.path(PATH.BW, "h3k27ac.bw"),
                      H3K9me3    = file.path(PATH.BW, "h3k9me3.bw"),
                      POL2       = file.path(PATH.BW, "pol2.bw"),
                      H3K27ac_ko = file.path(PATH.BW, "h3k27ac_ko.bw"),
                      H3K9me3_ko = file.path(PATH.BW, "h3k9me3_ko.bw"))
DATA_SOURCE      <- "ChIP-seq"

# ── FOR CUT ──────────────────────────────────────────────────────────────

PATH.BW <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/cut/Binned_Coverage"
BW_FILES   <- c(H3K4me3= file.path(PATH.BW ,"K4AB_WT.bw"), 
                #H3K9me3= file.path(PATH.BW ,"WT_K9AB.bw"), 
                POL2 = file.path(PATH.BW ,"WT_P2AB.bw"))
DATA_SOURCE <- "CUT&Tag"

DESEQ_FILE       <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
LFC_THRESHOLD    <- 1
PADJ_THRESHOLD   <- 0.05
TPM_THRESHOLD    <- 1
STD_CHROMS       <- paste0("chr", c(1:22, "X", "Y", "M"))
N_WORKERS        <- 11
N_BINS_TSS       <- 40
N_BINS_TES       <- 40
N_BINS_EXON <- 50   # 1 bin per 2bp for a 100bp exon — reasonable resolution
N_BINS_BODY <- 100  # keep for TSS/TES flanks and introns
N_BINS_INTRON    <- 100
FLANK            <- 2000
MIN_WIDTH_EXON   <- 50
MIN_WIDTH_INTRON <- 200

# ── 1. Build TxDb ─────────────────────────────────────────────────────────────
txdb <- makeTxDbFromGFF(GTF, format = "gtf")

# ── 2. Build genes GRanges ────────────────────────────────────────────────────
genes_gr         <- genes(txdb)
genes_gr$gene_id <- sub("\\..*", "", names(genes_gr))
seqlevelsStyle(genes_gr) <- "UCSC"
genes_gr <- keepSeqlevels(genes_gr,
                          intersect(seqlevels(genes_gr), STD_CHROMS),
                          pruning.mode = "coarse")

# ── 3. Build per-gene exon GRangesList ────────────────────────────────────────
exons_by_gene_gr <- exonsBy(txdb, by = "gene")
names(exons_by_gene_gr) <- sub("\\..*", "", names(exons_by_gene_gr))
seqlevelsStyle(exons_by_gene_gr) <- "UCSC"
exons_by_gene_gr <- keepSeqlevels(exons_by_gene_gr,
                                  intersect(seqlevels(exons_by_gene_gr),
                                            STD_CHROMS),
                                  pruning.mode = "coarse")

# ── 4. Build per-gene intron GRangesList ──────────────────────────────────────
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

# ── 5. Load expression data ───────────────────────────────────────────────────
expr <- read_tsv(EXPR_FILE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "",
                              sapply(strsplit(target_id, "\\|"), `[`, 2))) %>%
  dplyr::select(gene_id, tpm) %>%
  dplyr::rename(TPM = tpm) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

# ── 6. Load DESeq2 results ────────────────────────────────────────────────────
deseq_res <- read_csv(DESEQ_FILE) %>%
  dplyr::mutate(
    gene_id   = sub("\\..*", "", gene_id),
    de_status = dplyr::case_when(
      is.na(padj)                                              ~ "Filtered",
      padj >= PADJ_THRESHOLD                                   ~ "Not DE",
      padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD ~ "Down",
      TRUE                                                     ~ "Not DE"
    ),
    de_status = factor(de_status,
                       levels = c("Up", "Not DE", "Down", "Filtered"))
  )

cat("DE status counts:\n")
print(table(deseq_res$de_status))

# ── 7. Add expression and DE groups to genes_gr ───────────────────────────────
de_expr_joined <- dplyr::tibble(gene_id = genes_gr$gene_id) %>%
  dplyr::left_join(
    deseq_res %>% dplyr::select(gene_id, log2FoldChange, padj, de_status),
    by = "gene_id"
  ) %>%
  dplyr::left_join(expr %>% dplyr::select(gene_id, TPM), by = "gene_id") %>%
  dplyr::mutate(
    TPM        = tidyr::replace_na(TPM, 0),
    expr_group = dplyr::case_when(
      TPM < TPM_THRESHOLD                                             ~ "Unexpressed",
      TPM < quantile(TPM[TPM >= TPM_THRESHOLD], 0.33, na.rm = TRUE) ~ "Low",
      TPM < quantile(TPM[TPM >= TPM_THRESHOLD], 0.67, na.rm = TRUE) ~ "Medium",
      TRUE                                                            ~ "High"
    ),
    expr_group = factor(expr_group,
                        levels = c("Unexpressed", "Low", "Medium", "High")),
    de_status  = dplyr::if_else(
      is.na(de_status),
      factor("Not DE", levels = levels(deseq_res$de_status)),
      de_status
    )
  )

genes_gr$expr_group <- de_expr_joined$expr_group
genes_gr$de_status  <- de_expr_joined$de_status

cat("Expression group counts:\n")
print(table(genes_gr$expr_group, useNA = "always"))
cat("DE status counts on genes_gr:\n")
print(table(genes_gr$de_status, useNA = "always"))

# ── 8. Build bins ─────────────────────────────────────────────────────────────
make_tss_bins <- function(genes_gr, flank = FLANK, n_bins = N_BINS_TSS) {
  
  tss_gr      <- GenomicRanges::resize(genes_gr, width = 1, fix = "start")
  tss_flanked <- GenomicRanges::flank(tss_gr, width = flank, both = TRUE)
  tss_flanked <- keepSeqlevels(tss_flanked,
                               intersect(seqlevels(tss_flanked), STD_CHROMS),
                               pruning.mode = "coarse")
  tiled       <- tile(tss_flanked, n = n_bins)
  flat        <- unlist(tiled, use.names = FALSE)
  n_features  <- length(genes_gr)
  bin_idx     <- rep.int(seq_len(n_bins), n_features)
  feature_idx <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus    <- as.logical(strand(genes_gr) == "-")
  minus_mask  <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "TSS",
    gene_id     = rep.int(genes_gr$gene_id, rep.int(n_bins, n_features))
  )
  flat
}

make_tes_bins <- function(genes_gr, flank = FLANK, n_bins = N_BINS_TES) {
  
  tes_gr      <- GenomicRanges::resize(genes_gr, width = 1, fix = "end")
  tes_flanked <- GenomicRanges::flank(tes_gr, width = flank, both = TRUE)
  tes_flanked <- keepSeqlevels(tes_flanked,
                               intersect(seqlevels(tes_flanked), STD_CHROMS),
                               pruning.mode = "coarse")
  tiled       <- tile(tes_flanked, n = n_bins)
  flat        <- unlist(tiled, use.names = FALSE)
  n_features  <- length(genes_gr)
  bin_idx     <- rep.int(seq_len(n_bins), n_features)
  feature_idx <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus    <- as.logical(strand(genes_gr) == "-")
  minus_mask  <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "TES",
    gene_id     = rep.int(genes_gr$gene_id, rep.int(n_bins, n_features))
  )
  flat
}

make_exon_bins <- function(exons_by_gene, n_bins = N_BINS_EXON,
                           min_width = MIN_WIDTH_EXON) {
  
  exons_reduced      <- reduce(exons_by_gene)
  flat_exons         <- unlist(exons_reduced, use.names = TRUE)
  flat_exons$gene_id <- sub("\\..*", "", names(flat_exons))
  flat_exons         <- flat_exons[width(flat_exons) >= min_width]
  seqlevelsStyle(flat_exons) <- "UCSC"
  flat_exons <- keepSeqlevels(flat_exons,
                              intersect(seqlevels(flat_exons), STD_CHROMS),
                              pruning.mode = "coarse")
  too_narrow <- width(flat_exons) < n_bins
  if (any(too_narrow)) {
    message(sum(too_narrow), " exons removed: narrower than ", n_bins, " bins")
    flat_exons <- flat_exons[!too_narrow]
  }
  tiled         <- tile(flat_exons, n = n_bins)
  flat          <- unlist(tiled, use.names = FALSE)
  n_features    <- length(flat_exons)
  bin_idx       <- rep.int(seq_len(n_bins), n_features)
  feature_idx   <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus      <- as.logical(strand(flat_exons) == "-")
  minus_mask    <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "Exon",
    gene_id     = rep.int(flat_exons$gene_id, rep.int(n_bins, n_features))
  )
  flat
}

make_intron_bins <- function(introns_by_gene, n_bins = N_BINS_INTRON,
                             min_width = MIN_WIDTH_INTRON) {
  
  introns_reduced      <- reduce(introns_by_gene)
  flat_introns         <- unlist(introns_reduced, use.names = TRUE)
  flat_introns$gene_id <- sub("\\..*", "", names(flat_introns))
  flat_introns         <- flat_introns[width(flat_introns) >= min_width]
  seqlevelsStyle(flat_introns) <- "UCSC"
  flat_introns <- keepSeqlevels(flat_introns,
                                intersect(seqlevels(flat_introns), STD_CHROMS),
                                pruning.mode = "coarse")
  too_narrow <- width(flat_introns) < n_bins
  if (any(too_narrow)) {
    message(sum(too_narrow), " introns removed: narrower than ",
            n_bins, " bins")
    flat_introns <- flat_introns[!too_narrow]
  }
  tiled         <- tile(flat_introns, n = n_bins)
  flat          <- unlist(tiled, use.names = FALSE)
  n_features    <- length(flat_introns)
  bin_idx       <- rep.int(seq_len(n_bins), n_features)
  feature_idx   <- rep.int(seq_len(n_features), rep.int(n_bins, n_features))
  is_minus      <- as.logical(strand(flat_introns) == "-")
  minus_mask    <- rep.int(is_minus, rep.int(n_bins, n_features))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = "Intron",
    gene_id     = rep.int(flat_introns$gene_id, rep.int(n_bins, n_features))
  )
  flat
}

cat("Building TSS bins...\n")
tss_bins    <- make_tss_bins(genes_gr)
cat("Building TES bins...\n")
tes_bins    <- make_tes_bins(genes_gr)
cat("Building exon bins...\n")
exon_bins   <- make_exon_bins(exons_by_gene_gr)
cat("Building intron bins...\n")
intron_bins <- make_intron_bins(introns_by_gene_gr)

body_bins <- c(tss_bins, exon_bins, tes_bins)
rm(tss_bins, tes_bins, exon_bins); gc()

cat("Body bins size:",   format(object.size(body_bins),   units = "MB"), "\n")
cat("Intron bins size:", format(object.size(intron_bins), units = "MB"), "\n")

# ── 9. Extract signal per chromosome ─────────────────────────────────────────
extract_profile_signal_one_chrom <- function(chrom, bw_gr, bins) {
  
  bw_chrom   <- bw_gr[as.character(seqnames(bw_gr))   == chrom]
  bins_chrom <- bins[as.character(seqnames(bins))      == chrom]
  
  if (length(bw_chrom) == 0 || length(bins_chrom) == 0) return(NULL)
  
  tryCatch({
    plyranges::join_overlap_inner(bins_chrom, bw_chrom) %>%
      dplyr::as_tibble() %>%
      dplyr::rename(overlap_width = width) %>%
      dplyr::group_by(gene_id, feature, feature_idx, bin_idx) %>%
      dplyr::summarise(
        mean_signal = sum(score * overlap_width) / sum(overlap_width),
        .groups = "drop"
      )
  }, error = function(e) NULL)
}

extract_profile_signal <- function(bw_gr, bins, n_workers = N_WORKERS) {
  
  chroms_to_run <- intersect(
    as.character(unique(seqnames(bins))),
    as.character(unique(seqnames(bw_gr)))
  )
  
  BiocParallel::register(BiocParallel::MulticoreParam(workers = n_workers))
  
  chrom_results <- BiocParallel::bplapply(chroms_to_run, function(chrom) {
    extract_profile_signal_one_chrom(
      chrom = chrom,
      bw_gr = bw_gr,
      bins  = bins
    )
  })
  
  dplyr::bind_rows(chrom_results)
}

# ── 10. Aggregate and normalize ───────────────────────────────────────────────
compute_profile <- function(signal_df, genes_gr, sample_name,
                            n_bins_tss  = N_BINS_TSS,
                            n_bins_body = N_BINS_EXON,
                            n_bins_tes  = N_BINS_TES) {
  
  genes_meta <- dplyr::as_tibble(mcols(genes_gr)) %>%
    dplyr::select(gene_id, expr_group, de_status)
  
  signal_df %>%
    dplyr::mutate(sample = sample_name) %>%
    dplyr::left_join(genes_meta, by = "gene_id") %>%
    dplyr::filter(!is.na(expr_group),
                  !is.na(de_status),
                  de_status %in% c("Up", "Down")) %>%
    dplyr::mutate(
      x_pos = dplyr::case_when(
        .data$feature == "TSS"    ~ as.numeric(.data$bin_idx),
        .data$feature == "Exon"   ~ as.numeric(.data$bin_idx) + n_bins_tss,
        .data$feature == "TES"    ~ as.numeric(.data$bin_idx) + n_bins_tss +
          n_bins_body,
        .data$feature == "Intron" ~ as.numeric(.data$bin_idx)
      )
    ) %>%
    dplyr::group_by(sample, .data$feature, expr_group, de_status, x_pos) %>%
    dplyr::summarise(
      mean  = mean(mean_signal, na.rm = TRUE),
      se    = sd(mean_signal,   na.rm = TRUE) /
        sqrt(sum(!is.na(mean_signal))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      upper = mean + 1.96 * se,
      lower = mean - 1.96 * se
    )
}

# ── 11. Run pipeline ──────────────────────────────────────────────────────────
profile_results <- lapply(names(BW_FILES), function(sample_name,
                                                    body_bins,
                                                    intron_bins,
                                                    genes_gr,
                                                    n_bins_tss,
                                                    n_bins_body,
                                                    n_bins_tes,
                                                    n_workers) {
  cat("Processing:", sample_name, "\n")
  
  bw_gr <- plyranges::read_bigwig(BW_FILES[sample_name])
  seqlevelsStyle(bw_gr) <- "UCSC"
  bw_gr <- keepSeqlevels(bw_gr,
                         intersect(seqlevels(bw_gr), STD_CHROMS),
                         pruning.mode = "coarse")
  
  cat("  Extracting body signal...\n")
  body_signal   <- extract_profile_signal(bw_gr, body_bins,
                                          n_workers = n_workers)
  
  cat("  Extracting intron signal...\n")
  intron_signal <- extract_profile_signal(bw_gr, intron_bins,
                                          n_workers = n_workers)
  
  rm(bw_gr); gc()
  
  all_signal <- dplyr::bind_rows(body_signal, intron_signal)
  rm(body_signal, intron_signal); gc()
  
  cat("  Computing profile...\n")
  compute_profile(all_signal, genes_gr, sample_name,
                  n_bins_tss  = n_bins_tss,
                  n_bins_body = n_bins_body,
                  n_bins_tes  = n_bins_tes)
  
}, body_bins   = body_bins,
intron_bins = intron_bins,
genes_gr    = genes_gr,
n_bins_tss  = N_BINS_TSS,
n_bins_body = N_BINS_EXON,
n_bins_tes  = N_BINS_TES,
n_workers   = N_WORKERS)

profile_df <- dplyr::bind_rows(profile_results)

# ── 12. Plot gene body profile ────────────────────────────────────────────────

# Tol's "Muted" palette ## LETS USE THIS ONE
sample_colors <- c(
  H3K27ac    = "#CC6677",  ## rose
  H3K4me3    = "#CC6677",  ## rose
  H3K27ac_ko = "#DDCC77",  ## sand
  H3K9me3    = "#332288",  ## indigo
  H3K9me3_ko = "#88CCEE",  ## cyan
  POL2       = "#44AA99"   ## teal
)


x_breaks <- c(1,
              N_BINS_TSS / 2,
              N_BINS_TSS + 1,
              N_BINS_TSS + N_BINS_EXON,
              N_BINS_TSS + N_BINS_EXON + N_BINS_TES / 2,
              N_BINS_TSS + N_BINS_EXON + N_BINS_TES)
x_labels <- c(paste0("-", FLANK / 1000, "kb"),
              "TSS",
              "5'", "3'",
              "TES",
              paste0("+", FLANK / 1000, "kb"))

boundary_lines <- data.frame(
  xintercept = c(N_BINS_TSS / 2,                              # TSS midpoint
                 N_BINS_TSS + N_BINS_EXON + N_BINS_TES / 2), # TES midpoint
  label      = c("TSS", "TES")
)

# ── Build label data from genes_gr ────────────────────────────────────────────
label_df <- dplyr::as_tibble(mcols(genes_gr)) %>%
  dplyr::filter(de_status %in% c("Up", "Down")) %>%
  dplyr::group_by(expr_group, de_status) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("DEGs (n) = ", n_genes))

label_df


body_data   <- profile_df %>% dplyr::filter(.data$feature != "Intron")
intron_data <- profile_df %>% dplyr::filter(.data$feature == "Intron")

cat("Body data rows:",   nrow(body_data), "\n")
cat("Intron data rows:", nrow(intron_data), "\n")

# ── Gene body plot with labels ────────────────────────────────────────────────
body_data   <- profile_df %>% dplyr::filter(.data$feature != "Intron")
unique(body_data$sample)

#body_data<- body_data %>% filter(body_data$sample %in% c("H3K27ac", "H3K9me3", "POL2"))
#body_data<- body_data %>% filter(body_data$sample %in% c("H3K27ac","H3K27ac_ko", "POL2"))
#body_data<- body_data %>% filter(body_data$sample %in% c("POL2"))
#body_data<- body_data %>% filter(body_data$sample %in% c( "H3K9me3", "H3K9me3_ko"))

p_body <- ggplot(body_data,
                 aes(x = x_pos, y = mean, color = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(data        = boundary_lines,
             aes(xintercept = xintercept),
             linetype    = "dashed",
             color       = "grey50",
             linewidth   = 0.4,
             inherit.aes = FALSE) +
  geom_text(data        = label_df,
            aes(x       = Inf,
                y       = Inf,
                label   = label),
            inherit.aes = FALSE,
            hjust       = 1.1,
            vjust       = 1.5,
            size        = 3,
            color       = "grey30") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values  = sample_colors) +
  facet_grid(expr_group ~ de_status) +
  labs(x        = NULL,
       y        = paste0("Mean ", DATA_SOURCE, " signal"),
       color    = "Sample",
       fill     = "Sample",
       title    = paste0("Gene body ", DATA_SOURCE, " profile"),
       subtitle = paste0("TSS ±", FLANK / 1000, "kb — scaled exon body — TES ±",
                         FLANK / 1000, "kb")) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "top")

print(p_body)



intron_data <- profile_df %>% dplyr::filter(.data$feature == "Intron")
unique(intron_data$sample)

# pol2.k27.deg.tpm.meta
#intron_data<- intron_data %>% filter(intron_data$sample %in% c("H3K27ac", "H3K9me3", "POL2"))
#intron_data<- intron_data %>% filter(intron_data$sample %in% c("H3K27ac","H3K27ac_ko", "POL2"))
#intron_data<- intron_data %>% filter(intron_data$sample %in% c("POL2"))
#intron_data<- intron_data %>% filter(intron_data$sample %in% c( "H3K9me3", "H3K9me3_ko"))

# ── Intron plot with labels ───────────────────────────────────────────────────
p_intron <- ggplot(intron_data,
                   aes(x = x_pos, y = mean, color = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_text(data        = label_df,
            aes(x       = Inf,
                y       = Inf,
                label   = label),
            inherit.aes = FALSE,
            hjust       = 1.1,
            vjust       = 1.5,
            size        = 3,
            color       = "grey30") +
  scale_x_continuous(breaks = c(1, N_BINS_INTRON / 2, N_BINS_INTRON),
                     labels = c("5'", "50%", "3'")) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values  = sample_colors) +
  facet_grid(expr_group ~ de_status) +
  labs(x        = "Position (% of intron, 5'→3')",
       y        = paste0("Mean ", DATA_SOURCE, " signal"),
       color    = "Sample",
       fill     = "Sample",
       title    = paste0("Intron ", DATA_SOURCE, " profile"),
       subtitle = "Scaled to % of intron length") +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "top")


print(p_intron)

