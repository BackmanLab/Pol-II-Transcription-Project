##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Intron length vs. ChIP-seq peak signal. Relates narrowPeak signal
# (H3K27ac, H3K9me3, POL2 +/- POLR2A-KO) to intron length,
# stratified by DEG status and expression.
#
# Figure output: Fig. 7f (confirmed -- "Cumulative fraction of
# total signal" / "ChIP mark" exact text match).
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
library(scales)

# ── Parameters ────────────────────────────────────────────────────────────────
GTF              <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
EXPR_FILE        <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/HCT116_RNAseq_quantification.tsv"
DESEQ_FILE       <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
PATH.NP          <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks/"
PEAK_FILES       <- c(H3K27ac    = file.path(PATH.NP, "h3k27ac.narrowPeak"),
                      H3K9me3    = file.path(PATH.NP, "h3k9me3.narrowPeak"),
                      POL2       = file.path(PATH.NP, "pol2.narrowPeak"),
                      H3K27ac_ko = file.path(PATH.NP, "h3k27ac_ko.narrowPeak"),
                      H3K9me3_ko = file.path(PATH.NP, "h3k9me3_ko.narrowPeak"))
LFC_THRESHOLD    <- 1
PADJ_THRESHOLD   <- 0.05
TPM_THRESHOLD    <- 1
STD_CHROMS       <- paste0("chr", c(1:22, "X", "Y", "M"))
N_WORKERS        <- 11
MIN_WIDTH_INTRON <- 200

# ── Color palette ─────────────────────────────────────────────────────────────
# Tol's "Muted" palette ## LETS USE THIS ONE
sample_colors <- c(
  H3K27ac    = "#CC6677",  ## rose
  H3K4me3    = "#CC6677",  ## rose
  H3K27ac_ko = "#DDCC77",  ## sand
  H3K9me3    = "#332288",  ## indigo
  H3K9me3_ko = "#88CCEE",  ## cyan
  POL2       = "#44AA99"   ## teal
)

# ── 1. Build TxDb ─────────────────────────────────────────────────────────────
txdb <- makeTxDbFromGFF(GTF, format = "gtf")

# ── 2. Build genes GRanges ────────────────────────────────────────────────────
genes_gr         <- genes(txdb)
genes_gr$gene_id <- sub("\\..*", "", names(genes_gr))
seqlevelsStyle(genes_gr) <- "UCSC"
genes_gr <- keepSeqlevels(genes_gr,
                          intersect(seqlevels(genes_gr), STD_CHROMS),
                          pruning.mode = "coarse")

# ── 3. Build per-gene intron GRangesList ──────────────────────────────────────
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

# ── 4. Build flat intron GRanges carrying gene_id ────────────────────────────
introns_reduced      <- reduce(introns_by_gene_gr)
introns_flat         <- unlist(introns_reduced, use.names = TRUE)
introns_flat$gene_id <- sub("\\..*", "", names(introns_flat))
introns_flat         <- introns_flat[width(introns_flat) >= MIN_WIDTH_INTRON]
seqlevelsStyle(introns_flat) <- "UCSC"
introns_flat <- keepSeqlevels(introns_flat,
                              intersect(seqlevels(introns_flat), STD_CHROMS),
                              pruning.mode = "coarse")
introns_flat$intron_idx   <- seq_along(introns_flat)
introns_flat$intron_width <- width(introns_flat)

cat("Introns for analysis:", length(introns_flat), "\n")

# ── 5. Load expression data ───────────────────────────────────────────────────
expr <- read_tsv(EXPR_FILE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "",
                              sapply(strsplit(target_id, "\\|"), `[`, 2))) %>%
  dplyr::select(gene_id, tpm) %>%
  dplyr::rename(TPM = tpm) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

cat("Overlapping gene IDs:",
    sum(unique(expr$gene_id) %in% unique(genes_gr$gene_id)), "\n")

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
gene_meta <- dplyr::tibble(gene_id = genes_gr$gene_id) %>%
  dplyr::left_join(expr, by = "gene_id") %>%
  dplyr::left_join(
    deseq_res %>% dplyr::select(gene_id, log2FoldChange, padj, de_status),
    by = "gene_id"
  ) %>%
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

genes_gr$expr_group <- gene_meta$expr_group
genes_gr$de_status  <- gene_meta$de_status

cat("Expression group counts:\n")
print(table(genes_gr$expr_group, useNA = "always"))
cat("DE status counts on genes_gr:\n")
print(table(genes_gr$de_status, useNA = "always"))

# ── 8. Build intron metadata map ──────────────────────────────────────────────
intron_meta_map <- dplyr::tibble(
  intron_idx = introns_flat$intron_idx,
  gene_id    = introns_flat$gene_id
) %>%
  dplyr::left_join(
    dplyr::as_tibble(mcols(genes_gr)) %>%
      dplyr::select(gene_id, expr_group, de_status),
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    de_status = dplyr::if_else(
      is.na(de_status),
      factor("Not DE", levels = levels(deseq_res$de_status)),
      de_status
    )
  )

# ── 9. Load narrowPeak files ──────────────────────────────────────────────────
load_narrowpeak <- function(path) {
  read.table(path, header = FALSE, sep = "\t",
             colClasses = c("character", "integer",  "integer",
                            "character", "integer",  "character",
                            "numeric",   "numeric",  "numeric",
                            "integer")) %>%
    dplyr::rename(
      seqnames    = V1, start = V2, end = V3,
      name        = V4, score = V5, strand = V6,
      signalValue = V7, pValue = V8, qValue = V9,
      peak        = V10
    ) %>%
    dplyr::mutate(
      start  = start + 1L,
      strand = dplyr::if_else(strand == ".", "*", strand)
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    {
      seqlevelsStyle(.) <- "UCSC"
      keepSeqlevels(., intersect(seqlevels(.), STD_CHROMS),
                    pruning.mode = "coarse")
    }
}

# ── 10. Overlap peaks with introns per chromosome ─────────────────────────────
get_intron_signal_one_chrom <- function(chrom, peaks_gr, introns_gr) {
  peaks_chrom   <- peaks_gr[as.character(seqnames(peaks_gr))    == chrom]
  introns_chrom <- introns_gr[as.character(seqnames(introns_gr)) == chrom]
  if (length(peaks_chrom) == 0 || length(introns_chrom) == 0) return(NULL)
  tryCatch({
    plyranges::join_overlap_inner(introns_chrom, peaks_chrom) %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(intron_idx, intron_width, gene_id) %>%
      dplyr::summarise(
        mean_signal = mean(signalValue, na.rm = TRUE),
        n_peaks     = dplyr::n(),
        .groups     = "drop"
      )
  }, error = function(e) NULL)
}

# ── 11. Run across all peak files ─────────────────────────────────────────────
scatter_results <- lapply(names(PEAK_FILES), function(sample_name) {
  cat("Processing:", sample_name, "\n")
  
  peaks_gr <- load_narrowpeak(PEAK_FILES[sample_name])
  
  chroms_to_run <- intersect(
    as.character(unique(seqnames(introns_flat))),
    as.character(unique(seqnames(peaks_gr)))
  )
  
  BiocParallel::register(BiocParallel::MulticoreParam(workers = N_WORKERS))
  
  chrom_results <- BiocParallel::bplapply(chroms_to_run, function(chrom) {
    get_intron_signal_one_chrom(
      chrom      = chrom,
      peaks_gr   = peaks_gr,
      introns_gr = introns_flat
    )
  })
  
  rm(peaks_gr); gc()
  
  dplyr::bind_rows(chrom_results) %>%
    dplyr::mutate(sample = sample_name)
})

# ── 12. Join metadata and build scatter_df ────────────────────────────────────
scatter_df <- dplyr::bind_rows(scatter_results) %>%
  dplyr::left_join(
    intron_meta_map %>%
      dplyr::select(intron_idx, expr_group, de_status),
    by = "intron_idx"
  ) %>%
  dplyr::filter(!is.na(expr_group))

cat("Scatter data rows:", nrow(scatter_df), "\n")

unique(scatter_df$sample)
scatter_df <-scatter_df %>% filter(scatter_df$sample %in% c("H3K27ac", "H3K9me3", "POL2" ))

# ── 13. Label data — expression only ─────────────────────────────────────────
label_expr <- scatter_df %>%
  dplyr::distinct(gene_id, expr_group) %>%
  dplyr::group_by(expr_group) %>%
  dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n_genes))

# ── 14. Scatter plot — faceted by expression ──────────────────────────────────
p_scatter <- ggplot(scatter_df,
                    aes(x = intron_width, y = mean_signal, color = sample)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  geom_text(data        = label_expr,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = sample_colors) +
  facet_grid(sample ~ expr_group, scales = "free_y") +
  labs(x        = "Intron length (bp, log scale)",
       y        = "Mean peak signalValue",
       title    = "Intron length vs ChIP-seq peak signal",
       subtitle = "Stratified by expression level") +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

# ── 15. Cumulative signal plot — faceted by expression ────────────────────────
scatter_df_cum <- scatter_df %>%
  dplyr::group_by(sample, expr_group) %>%
  dplyr::arrange(intron_width) %>%
  dplyr::mutate(
    cum_signal = cumsum(mean_signal) / sum(mean_signal)
  ) %>%
  dplyr::ungroup()

p_cumulative <- ggplot(scatter_df_cum,
                       aes(x = intron_width, y = cum_signal,
                           color = sample)) +
  geom_line(linewidth = 0.9) +
  geom_text(data        = label_expr,
            aes(x = Inf, y = 0, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = -0.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_x_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_color_manual(values = sample_colors) +
  facet_wrap(~ expr_group, nrow = 1) +
  labs(x        = "Intron length (bp, log10)",
       y        = "Cumulative fraction of total signal",
       colour   = "ChIP mark",
       title    = "Cumulative ChIP-seq signal over intron lengths",
       subtitle = "Stratified by expression level") +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "right",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

# ── 16. Density plot — faceted by expression ──────────────────────────────────
p_density <- ggplot(scatter_df,
                    aes(x = intron_width, colour = sample)) +
  geom_density(linewidth = 0.8, adjust = 1.5) +
  geom_text(data        = label_expr,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_x_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_colour_manual(values = sample_colors) +
  facet_wrap(~ expr_group, nrow = 1) +
  labs(
    x        = "Intron length (bp, log10)",
    y        = "Density",
    colour   = "ChIP mark",
    title    = "Intron length distribution by ChIP mark",
    subtitle = "Stratified by expression level"
  ) +
  #theme_bw() +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "right",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

# ── 17. Print all plots ───────────────────────────────────────────────────────
print(p_scatter)
print(p_cumulative)
print(p_density)
