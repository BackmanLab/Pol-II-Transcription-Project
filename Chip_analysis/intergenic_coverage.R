##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# ChIP-seq signal over intergenic regions vs. expression. Builds
# intergenic intervals from gene-to-gene gaps, computes ChIP-seq
# (H3K27ac, H3K9me3, POL2 +/- POLR2A-KO) mean signal per interval,
# and relates intergenic signal/distance to upstream gene
# expression tertile.
#
# Analyses:
#   1. p_cumulative -- cumulative fraction of signal vs. intergenic
#      distance, by expression tertile
#   2. p_density -- density of intergenic distance, by expression
#      tertile
#
# Figure output: Fig. 7g (p_density, confirmed). p_cumulative
# structurally matches Fig. 7f but that panel's axis is intronic
# (not intergenic) distance -- likely the same template applied to
# intron gap widths elsewhere; exact source unconfirmed.
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
PATH.NP          <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks/"
PEAK_FILES       <- c(H3K27ac    = file.path(PATH.NP, "h3k27ac.narrowPeak"),
                      H3K9me3    = file.path(PATH.NP, "h3k9me3.narrowPeak"),
                      POL2       = file.path(PATH.NP, "pol2.narrowPeak"),
                      H3K27ac_ko = file.path(PATH.NP, "h3k27ac_ko.narrowPeak"),
                      H3K9me3_ko = file.path(PATH.NP, "h3k9me3_ko.narrowPeak"))
TPM_THRESHOLD    <- 1
STD_CHROMS       <- paste0("chr", c(1:22, "X", "Y", "M"))
N_WORKERS        <- 11
MIN_GAP          <- 250    # minimum intergenic distance to include

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

# ── 3. Load expression data ───────────────────────────────────────────────────
expr <- read_tsv(EXPR_FILE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "",
                              sapply(strsplit(target_id, "\\|"), `[`, 2))) %>%
  dplyr::select(gene_id, tpm) %>%
  dplyr::rename(TPM = tpm) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

# ── 4. Add expression groups to genes_gr ─────────────────────────────────────
expr_joined <- dplyr::tibble(gene_id = genes_gr$gene_id) %>%
  dplyr::left_join(expr, by = "gene_id") %>%
  dplyr::mutate(
    TPM        = tidyr::replace_na(TPM, 0),
    expr_group = dplyr::case_when(
      TPM < TPM_THRESHOLD                                             ~ "Unexpressed",
      TPM < quantile(TPM[TPM >= TPM_THRESHOLD], 0.33, na.rm = TRUE) ~ "Low",
      TPM < quantile(TPM[TPM >= TPM_THRESHOLD], 0.67, na.rm = TRUE) ~ "Medium",
      TRUE                                                            ~ "High"
    ),
    expr_group = factor(expr_group,
                        levels = c("Unexpressed", "Low", "Medium", "High"))
  )

genes_gr$expr_group <- expr_joined$expr_group
genes_gr$TPM        <- expr_joined$TPM

cat("Expression group counts:\n")
print(table(genes_gr$expr_group, useNA = "always"))

# ── 5. Compute intergenic regions ─────────────────────────────────────────────
# Sort genes by chromosome and position
genes_df <- as.data.frame(genes_gr) %>%
  dplyr::as_tibble() %>%
  dplyr::rename(chrom = seqnames) %>%
  dplyr::arrange(chrom, start) %>%
  dplyr::group_by(chrom) %>%
  dplyr::mutate(
    # Gene B attributes (next gene)
    next_gene_id  = dplyr::lead(gene_id),
    next_tss      = dplyr::lead(start),     # TSS of next gene
    next_strand   = dplyr::lead(strand),
    # Intergenic interval
    gap_start     = end + 1L,               # TES of current gene + 1
    gap_end       = dplyr::lead(start) - 1L # TSS of next gene - 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    !is.na(next_gene_id),
    gap_end > gap_start,                    # remove overlapping genes
    (gap_end - gap_start + 1) >= MIN_GAP    # remove very short gaps
  ) %>%
  dplyr::mutate(
    gap_width   = gap_end - gap_start + 1L,
    interval_id = dplyr::row_number()
  )


cat("Intergenic intervals found:", nrow(genes_df), "\n")
summary(genes_df$gap_width)

# Custom human readable breaks
gap_breaks_plot <- c(100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000)
gap_labels_plot <- c("100bp", "500bp", "1kb", "5kb", "10kb", "50kb", "100kb", "500kb", "1Mb")

hist_plot <- ggplot(genes_df, aes(x = gap_width)) +
  geom_histogram(bins = 100,
                 fill = "#0072B2", color = "white",
                 linewidth = 0.15, alpha = 0.85) +
  geom_vline(xintercept = median(genes_df$gap_width, na.rm = TRUE),
             linetype = "dashed", color = "#D55E00",
             linewidth = 0.7) +
  annotate("text",
           x     = median(genes_df$gap_width, na.rm = TRUE) * 1.15,
           y     = Inf,
           label = paste0("median = ",
                          scales::comma(round(
                            median(genes_df$gap_width, na.rm = TRUE)))),
           hjust = 0, vjust = 1.5,
           size = 3.5, color = "#D55E00", fontface = "italic") +
  scale_x_log10(
    breaks = c(100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000),
    labels = c("100bp", "500bp", "1kb", "5kb", "10kb",
               "50kb", "100kb", "500kb", "1Mb")
  ) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x        = "Intergenic distance (log scale)",
       y        = "Number of intervals",
       title    = "Distribution of intergenic distances",
       subtitle = paste0("n = ", scales::comma(nrow(genes_df)),
                         " intergenic regions")) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y     = element_text(size = 11),
    axis.title      = element_text(size = 13),
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(color = "grey40", size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    axis.line       = element_line(color = "grey30")
  )

print(hist_plot)

# ── 6. Add intergenic length groups ───────────────────────────────────────────
gap_breaks <- quantile(genes_df$gap_width,
                       probs = c(0, 0.33, 0.67, 1),
                       na.rm = TRUE)

message("Intergenic gap quantile breaks (bp):")
message(paste(names(gap_breaks), round(gap_breaks),
              sep = ": ", collapse = " | "))

genes_df <- genes_df %>%
  dplyr::mutate(
    gap_group = cut(gap_width,
                    breaks         = gap_breaks,
                    labels         = c("Short", "Medium", "Long"),
                    include.lowest = TRUE),
    gap_group = factor(gap_group,
                       levels = c("Short", "Medium", "Long"))
  )

# ── 7. Build intergenic GRanges ───────────────────────────────────────────────
intergenic_gr <- GenomicRanges::GRanges(
  seqnames = genes_df$chrom,
  ranges   = IRanges::IRanges(start = genes_df$gap_start,
                              end   = genes_df$gap_end),
  strand   = "*"
)
intergenic_gr$interval_id  <- genes_df$interval_id
intergenic_gr$gap_width     <- genes_df$gap_width
intergenic_gr$gap_group     <- genes_df$gap_group
intergenic_gr$upstream_gene <- genes_df$gene_id
intergenic_gr$expr_group    <- genes_df$expr_group
intergenic_gr$TPM           <- genes_df$TPM

seqlevelsStyle(intergenic_gr) <- "UCSC"
intergenic_gr <- keepSeqlevels(intergenic_gr,
                               intersect(seqlevels(intergenic_gr), STD_CHROMS),
                               pruning.mode = "coarse")

# ── 7.b filter overlaps ──────────────────────────────────────────────────
# Convert to GRanges and check for overlaps with genes
test_overlap <- findOverlaps(intergenic_gr, genes_gr)
cat("Intervals overlapping gene bodies:", 
    length(unique(queryHits(test_overlap))), "\n")
cat("Expected: 0\n")

# Keep ranges in gr1 that DO NOT overlap with gr2
intergenic_gr <- intergenic_gr %>% filter_by_non_overlaps(genes_gr)

# Convert to GRanges and check for overlaps with genes
test_overlap <- findOverlaps(intergenic_gr, genes_gr)
cat("Intervals overlapping gene bodies:", 
    length(unique(queryHits(test_overlap))), "\n")
cat("Expected: 0\n")

cat("Intervals remaining after filtering:", length(intergenic_gr), "\n")

# Also recheck the gap width distribution since some intervals may have changed
summary(intergenic_gr$gap_width)
hist(log10(intergenic_gr$gap_width), breaks = 100,
     main = "Intergenic distance distribution (post-filter)",
     xlab = "log10(gap width in bp)")

cat("Intergenic GRanges built:", length(intergenic_gr), "\n")

# Recompute gap groups on filtered intervals
gap_breaks_new <- quantile(intergenic_gr$gap_width,
                           probs = c(0, 0.33, 0.67, 1),
                           na.rm = TRUE)

message("Updated intergenic gap quantile breaks (bp):")
message(paste(names(gap_breaks_new), round(gap_breaks_new),
              sep = ": ", collapse = " | "))

intergenic_gr$gap_group <- cut(intergenic_gr$gap_width,
                               breaks         = gap_breaks_new,
                               labels         = c("Short", "Medium", "Long"),
                               include.lowest = TRUE)

cat("Gap group counts:\n")
print(table(intergenic_gr$gap_group))

# ── 7.c SANITY CHECK ───────────────────────────────────────────────
data.frame(intergenic_gr)

# Export top 20 intervals as BED for loading into genome browser ## Sanity check
data.frame(intergenic_gr)  %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::select(seqnames, start, end, interval_id,
                gap_width, upstream_gene) %>%
  dplyr::mutate(start = start - 1L) %>%  # convert back to 0-based BED
  readr::write_tsv("intergenic_check.bed", col_names = FALSE)

# Convert to GRanges and check for overlaps with genes
test_overlap <- findOverlaps(intergenic_gr, genes_gr)
cat("Intervals overlapping gene bodies:", 
    length(unique(queryHits(test_overlap))), "\n")
cat("Expected: 0\n")

# Check that gap_start is always less than gap_end
cat("Intervals with negative width:", 
    sum(genes_df$gap_width <= 0), "\n")

# Check that upstream TES and downstream TSS bracket the interval
genes_df %>%
  dplyr::mutate(
    tes_matches = (gap_start - 1L) == end,   # gap_start should be TES + 1
    tss_matches = (gap_end + 1L)  == next_tss # gap_end should be TSS - 1
  ) %>%
  dplyr::summarise(
    tes_correct = sum(tes_matches),
    tss_correct = sum(tss_matches),
    total       = dplyr::n()
  )


cat("Actual intervals before gap filter:", 
    genes_df %>% dplyr::filter(gap_end > gap_start) %>% nrow(), "\n")

# ── 8. Load narrowPeak files ──────────────────────────────────────────────────
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

# ── 9. Overlap peaks with intergenic regions per chromosome ───────────────────
get_intergenic_signal_one_chrom <- function(chrom, peaks_gr, intergenic_gr) {
  peaks_chrom      <- peaks_gr[as.character(seqnames(peaks_gr))      == chrom]
  intergenic_chrom <- intergenic_gr[as.character(seqnames(intergenic_gr)) == chrom]
  
  if (length(peaks_chrom) == 0 || length(intergenic_chrom) == 0) return(NULL)
  
  tryCatch({
    plyranges::join_overlap_inner(intergenic_chrom, peaks_chrom) %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(interval_id, gap_width, gap_group,
                      upstream_gene, expr_group, TPM) %>%
      dplyr::summarise(
        mean_signal = mean(signalValue, na.rm = TRUE),
        n_peaks     = dplyr::n(),
        .groups     = "drop"
      )
  }, error = function(e) NULL)
}

# ── 10. Run across all peak files ─────────────────────────────────────────────
intergenic_results <- lapply(names(PEAK_FILES), function(sample_name) {
  cat("Processing:", sample_name, "\n")
  
  peaks_gr <- load_narrowpeak(PEAK_FILES[sample_name])
  
  chroms_to_run <- intersect(
    as.character(unique(seqnames(intergenic_gr))),
    as.character(unique(seqnames(peaks_gr)))
  )
  
  BiocParallel::register(BiocParallel::MulticoreParam(workers = N_WORKERS))
  
  chrom_results <- BiocParallel::bplapply(chroms_to_run, function(chrom) {
    get_intergenic_signal_one_chrom(
      chrom         = chrom,
      peaks_gr      = peaks_gr,
      intergenic_gr = intergenic_gr
    )
  })
  
  rm(peaks_gr); gc()
  
  dplyr::bind_rows(chrom_results) %>%
    dplyr::mutate(sample = sample_name)
})

intergenic_df <- dplyr::bind_rows(intergenic_results)

cat("Intergenic data rows:", nrow(intergenic_df), "\n")

# ── 11. Label data ────────────────────────────────────────────────────────────
label_expr <- intergenic_df %>%
  dplyr::distinct(upstream_gene, expr_group) %>%
  dplyr::group_by(expr_group) %>%
  dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n_genes))

label_gap <- intergenic_df %>%
  dplyr::distinct(interval_id, gap_group) %>%
  dplyr::group_by(gap_group) %>%
  dplyr::summarise(n_intervals = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n_intervals))

# ── 12. Scatter plot ──────────────────────────────────────────────────────────

unique(intergenic_df$sample)
intergenic_df <- intergenic_df %>% filter(sample %in% c("H3K27ac", "H3K9me3", "POL2"))

p_scatter <- ggplot(intergenic_df,
                    aes(x = gap_width, y = mean_signal, color = sample)) +
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
  labs(x        = "Intergenic distance (bp, log scale)",
       y        = "Mean peak signalValue",
       title    = "Intergenic distance vs ChIP-seq peak signal",
       subtitle = "Stratified by upstream gene expression level") +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "none",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

# ── 13. Cumulative signal plot ────────────────────────────────────────────────
intergenic_cum <- intergenic_df %>%
  dplyr::group_by(sample, expr_group) %>%
  dplyr::arrange(gap_width) %>%
  dplyr::mutate(
    cum_signal = cumsum(mean_signal) / sum(mean_signal)
  ) %>%
  dplyr::ungroup()

p_cumulative <- ggplot(intergenic_cum,
                       aes(x = gap_width, y = cum_signal, color = sample)) +
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
  labs(x        = "Intergenic distance (bp, log10)",
       y        = "Cumulative fraction of total signal",
       colour   = "ChIP mark",
       title    = "Cumulative ChIP-seq signal over intergenic distances",
       subtitle = "Stratified by upstream gene expression level") +
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

# ── 14. Density plot ──────────────────────────────────────────────────────────
p_density <- ggplot(intergenic_df,
                    aes(x = gap_width, colour = sample)) +
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
    x        = "Intergenic distance (bp, log10)",
    y        = "Density",
    colour   = "ChIP mark",
    title    = "Intergenic distance distribution by ChIP mark",
    subtitle = "Stratified by upstream gene expression level"
  ) +
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

# ── 15. Gap group plot — faceted by gap length ────────────────────────────────
p_gap_groups <- ggplot(intergenic_df,
                       aes(x = expr_group, y = mean_signal, fill = sample)) +
  geom_violin(alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               fill = "white", color = "grey30") +
  geom_text(data        = label_gap,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_fill_manual(values = sample_colors) +
  facet_grid(sample ~ gap_group, scales = "free_y") +
  labs(x        = "Upstream gene expression",
       y        = "Mean peak signalValue",
       title    = "ChIP-seq signal in intergenic regions",
       subtitle = "Stratified by intergenic distance and upstream gene expression") +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "none",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

# ── 16. Print all plots ───────────────────────────────────────────────────────
print(p_scatter)
print(p_cumulative)
print(p_density) ## this one is amazing
print(p_gap_groups)

## check mean signal
intergenic_df %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    mean_signal = mean(mean_signal, na.rm = TRUE),
    median_signal = median(mean_signal, na.rm = TRUE),
    .groups = "drop"
  )

# ── 17. updated dens plot with quaretiles ───────────────────────────────────────────────────────

library(purrr)
unique(intergenic_df$sample)
#intergenic_df_quart <- intergenic_df %>% filter(intergenic_df$sample == "H3K27ac" | intergenic_df$sample == "H3K9me3")
intergenic_df_quart <- intergenic_df %>% filter(intergenic_df$sample %in% c( "H3K27ac", "H3K9me3", "H3K9me3_ko") )
sample_colors <- c(
  H3K27ac    = "#4575B4",
  H3K9me3    = "#D73027",
  #POL2       = "#1A9850",
  #H3K27ac_ko = "#F46D43",
  H3K9me3_ko = "black"
)

library(purrr)

# ── Compute density and quartile fills per sample per DE status ───────────────
density_df <- intergenic_df_quart  %>%
  dplyr::group_by(sample, expr_group) %>%
  dplyr::group_modify(~ {
    d   <- density(log10(.x$gap_width), adjust = 1.5, n = 512)
    q25 <- log10(quantile(.x$gap_width, 0.25, na.rm = TRUE))
    q50 <- log10(quantile(.x$gap_width, 0.50, na.rm = TRUE))
    q75 <- log10(quantile(.x$gap_width, 0.75, na.rm = TRUE))
    dplyr::tibble(
      x    = d$x,
      y    = d$y,
      q25  = q25,
      q50  = q50,
      q75  = q75,
      zone = dplyr::case_when(
        x < q25             ~ "Q1",
        x >= q25 & x < q50 ~ "Q2",
        x >= q50 & x < q75 ~ "Q3",
        x >= q75            ~ "Q4"
      )
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(zone_key = paste(sample, zone, sep = "_"))

# ── Quartile color palettes ───────────────────────────────────────────────────
quartile_colors <- list(
  #H3K27ac    = c(Q1 = "#FDDBC7", Q2 = "#F4A582", Q3 = "#D6604D", Q4 = "#B2182B"),
  H3K27ac    = c(Q1 = "white", Q2 = "white", Q3 = "white", Q4 = "white"),
  H3K9me3    = c(Q1 = "#FEE8C8", Q2 = "#FDBB84", Q3 = "#E34A33", Q4 = "#B30000"),
  #POL2       = c(Q1 = "#D9F0D3", Q2 = "#A6DBA0", Q3 = "#5AAE61", Q4 = "#1B7837"),
  #H3K27ac_ko = c(Q1 = "#FEE8C8", Q2 = "#FDBB84", Q3 = "#E34A33", Q4 = "#B30000"),
  #H3K9me3_ko = c(Q1 = "#EDF8FB", Q2 = "#B2E2E2", Q3 = "#66C2A4", Q4 = "#238B45")
  H3K9me3_ko = c(Q1 = "white", Q2 = "white", Q3 = "white", Q4 = "white")
)

color_lookup <- purrr::imap_dfr(quartile_colors, function(cols, samp) {
  dplyr::tibble(
    sample   = samp,
    zone     = names(cols),
    fill     = unname(cols),
    zone_key = paste(samp, names(cols), sep = "_")
  )
})

fill_values <- setNames(color_lookup$fill, color_lookup$zone_key)

# ── Median line data ──────────────────────────────────────────────────────────
median_df <- density_df %>%
  dplyr::distinct(sample, expr_group, q50)

# ── Line data — one line per sample per DE status ─────────────────────────────
line_df <- density_df %>%
  dplyr::group_by(sample, expr_group, x) %>%
  dplyr::summarise(y = mean(y), .groups = "drop")

# ── Plot ──────────────────────────────────────────────────────────────────────
p_density_quartile <- ggplot(density_df,
                             aes(x = x, y = y,
                                 group = zone_key,
                                 fill  = zone_key)) +
  geom_area(alpha = 0.8, color = NA) +
  geom_line(data      = line_df,
            aes(x = x, y = y, color = sample, group = sample),
            linewidth = 0.8,
            inherit.aes = FALSE) +
  geom_vline(data      = median_df,
             aes(xintercept = q50, color = sample),
             linetype  = "dashed",
             linewidth = 0.6,
             inherit.aes = FALSE) +
  geom_text(data        = label_expr ,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_color_manual(values = sample_colors) +
  scale_x_continuous(
    breaks = log10(c(100, 1000, 10000,
                     100000, 1000000)),
    labels = c("100bp", "1kb", 
               "10kb", "100kb", "1Mb")
  ) +
  facet_grid(~ de_status, scales = "free_y") +
  labs(
    x        = "Intergenic distance",
    y        = "Density",
    colour   = "ChIP mark",
    title    = "Intergenic distance distribution by ChIP mark",
    subtitle = "Shading shows quartile zones — dashed line shows median"
  ) +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.title       = element_text(size = 13),
    legend.position  = "right",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 11),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

p_density_quartile <- ggplot(density_df,
                             aes(x = x, y = y,
                                 group = zone_key,
                                 fill  = zone_key)) +
  geom_area(alpha = 0.75, color = NA, position = "identity") +
  geom_line(data        = density_df %>%
              dplyr::group_by(sample, expr_group, x) %>%
              dplyr::summarise(y = mean(y), .groups = "drop"),
            aes(x = x, y = y, color = sample, group = sample),
            linewidth   = 0.8,
            inherit.aes = FALSE) +
  geom_text(data        = label_expr ,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_color_manual(values = sample_colors) +
  scale_x_continuous(
    breaks = log10(c(100, 1000, 10000,
                     100000, 1000000)),
    labels = c("100bp", "1kb", 
               "10kb", "100kb", "1Mb")
  ) +
  facet_wrap(~ expr_group, nrow = 1) +
  labs(
    x        = "Intergenic distance",
    y        = "Density",
    colour   = "ChIP mark",
    title    = "Intergenic distance distribution",
    subtitle = "H3K27ac vs H3K9me3 — shading shows quartile zones"
  ) +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.title       = element_text(size = 13),
    legend.position  = "right",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 11),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

print(p_density_quartile)


#########################################################################################
# INTERGENIC COVERAGE X DEG STATUS
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
library(scales)

# ── Parameters ────────────────────────────────────────────────────────────────
GTF              <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
DESEQ_FILE       <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
PATH.NP          <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks/"
PEAK_FILES       <- c(H3K27ac    = file.path(PATH.NP, "h3k27ac.narrowPeak"),
                      H3K9me3    = file.path(PATH.NP, "h3k9me3.narrowPeak"),
                      POL2       = file.path(PATH.NP, "pol2.narrowPeak"),
                      H3K27ac_ko = file.path(PATH.NP, "h3k27ac_ko.narrowPeak"),
                      H3K9me3_ko = file.path(PATH.NP, "h3k9me3_ko.narrowPeak"))
LFC_THRESHOLD    <- 1
PADJ_THRESHOLD   <- 0.05
STD_CHROMS       <- paste0("chr", c(1:22, "X", "Y", "M"))
N_WORKERS        <- 7
MIN_GAP          <- 250    # minimum intergenic distance to include

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

# ── 3. Load DESeq2 results ────────────────────────────────────────────────────
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

# ── 4. Add DE status to genes_gr ──────────────────────────────────────────────
de_joined <- dplyr::tibble(gene_id = genes_gr$gene_id) %>%
  dplyr::left_join(
    deseq_res %>% dplyr::select(gene_id, log2FoldChange, padj, de_status),
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    de_status = dplyr::if_else(
      is.na(de_status),
      factor("Not DE", levels = levels(deseq_res$de_status)),
      de_status
    )
  )

genes_gr$de_status <- de_joined$de_status

cat("DE status counts on genes_gr:\n")
print(table(genes_gr$de_status, useNA = "always"))

# ── 5. Compute intergenic regions ─────────────────────────────────────────────
genes_df <- as.data.frame(genes_gr) %>%
  dplyr::as_tibble() %>%
  dplyr::rename(chrom = seqnames) %>%
  dplyr::arrange(chrom, start) %>%
  dplyr::group_by(chrom) %>%
  dplyr::mutate(
    next_gene_id  = dplyr::lead(gene_id),
    next_tss      = dplyr::lead(start),
    next_strand   = dplyr::lead(strand),
    gap_start     = end + 1L,
    gap_end       = dplyr::lead(start) - 1L
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    !is.na(next_gene_id),
    gap_end > gap_start,
    (gap_end - gap_start + 1) >= MIN_GAP
  ) %>%
  dplyr::mutate(
    gap_width   = gap_end - gap_start + 1L,
    interval_id = dplyr::row_number()
  )

cat("Intergenic intervals before overlap filter:", nrow(genes_df), "\n")

# ── 6. Add intergenic length groups ───────────────────────────────────────────
gap_breaks <- quantile(genes_df$gap_width,
                       probs = c(0, 0.33, 0.67, 1),
                       na.rm = TRUE)

message("Intergenic gap quantile breaks (bp):")
message(paste(names(gap_breaks), round(gap_breaks),
              sep = ": ", collapse = " | "))

genes_df <- genes_df %>%
  dplyr::mutate(
    gap_group = cut(gap_width,
                    breaks         = gap_breaks,
                    labels         = c("Short", "Medium", "Long"),
                    include.lowest = TRUE),
    gap_group = factor(gap_group,
                       levels = c("Short", "Medium", "Long"))
  )

# ── 7. Build intergenic GRanges ───────────────────────────────────────────────
intergenic_gr <- GenomicRanges::GRanges(
  seqnames = genes_df$chrom,
  ranges   = IRanges::IRanges(start = genes_df$gap_start,
                              end   = genes_df$gap_end),
  strand   = "*"
)
intergenic_gr$interval_id  <- genes_df$interval_id
intergenic_gr$gap_width     <- genes_df$gap_width
intergenic_gr$gap_group     <- genes_df$gap_group
intergenic_gr$upstream_gene <- genes_df$gene_id
intergenic_gr$de_status     <- genes_df$de_status

seqlevelsStyle(intergenic_gr) <- "UCSC"
intergenic_gr <- keepSeqlevels(intergenic_gr,
                               intersect(seqlevels(intergenic_gr), STD_CHROMS),
                               pruning.mode = "coarse")

# Remove intervals overlapping gene bodies
intergenic_gr <- intergenic_gr %>%
  plyranges::filter_by_non_overlaps(genes_gr)

cat("Intergenic intervals after filtering:", length(intergenic_gr), "\n")

# Recompute gap groups on filtered intervals
gap_breaks_filtered <- quantile(intergenic_gr$gap_width,
                                probs = c(0, 0.33, 0.67, 1),
                                na.rm = TRUE)

message("Updated gap quantile breaks (bp):")
message(paste(names(gap_breaks_filtered), round(gap_breaks_filtered),
              sep = ": ", collapse = " | "))

intergenic_gr$gap_group <- cut(intergenic_gr$gap_width,
                               breaks         = gap_breaks_filtered,
                               labels         = c("Short", "Medium", "Long"),
                               include.lowest = TRUE)

cat("Gap group counts:\n")
print(table(intergenic_gr$gap_group))

# ── 8. Load narrowPeak files ──────────────────────────────────────────────────
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

# ── 9. Overlap peaks with intergenic regions per chromosome ───────────────────
get_intergenic_signal_one_chrom <- function(chrom, peaks_gr, intergenic_gr) {
  peaks_chrom      <- peaks_gr[as.character(seqnames(peaks_gr))         == chrom]
  intergenic_chrom <- intergenic_gr[as.character(seqnames(intergenic_gr)) == chrom]
  if (length(peaks_chrom) == 0 || length(intergenic_chrom) == 0) return(NULL)
  tryCatch({
    plyranges::join_overlap_inner(intergenic_chrom, peaks_chrom) %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(interval_id, gap_width, gap_group,
                      upstream_gene, de_status) %>%
      dplyr::summarise(
        mean_signal = mean(signalValue, na.rm = TRUE),
        n_peaks     = dplyr::n(),
        .groups     = "drop"
      )
  }, error = function(e) NULL)
}

# ── 10. Run across all peak files ─────────────────────────────────────────────
intergenic_results <- lapply(names(PEAK_FILES), function(sample_name) {
  cat("Processing:", sample_name, "\n")
  
  peaks_gr <- load_narrowpeak(PEAK_FILES[sample_name])
  
  chroms_to_run <- intersect(
    as.character(unique(seqnames(intergenic_gr))),
    as.character(unique(seqnames(peaks_gr)))
  )
  
  BiocParallel::register(BiocParallel::MulticoreParam(workers = N_WORKERS))
  
  chrom_results <- BiocParallel::bplapply(chroms_to_run, function(chrom) {
    get_intergenic_signal_one_chrom(
      chrom         = chrom,
      peaks_gr      = peaks_gr,
      intergenic_gr = intergenic_gr
    )
  })
  
  rm(peaks_gr); gc()
  
  dplyr::bind_rows(chrom_results) %>%
    dplyr::mutate(sample = sample_name)
})

intergenic_df <- dplyr::bind_rows(intergenic_results) %>%
  dplyr::filter(de_status %in% c("Up", "Down", "Not DE")) %>%
  dplyr::mutate(
    de_status = dplyr::recode(de_status,
                              "Up"     = "Upregulated",
                              "Down"   = "Downregulated",
                              "Not DE" = "No Change"),
    de_status = factor(de_status,
                       levels = c("Upregulated", "No Change", "Downregulated"))
  )

cat("Intergenic data rows:", nrow(intergenic_df), "\n")

unique(intergenic_df$sample)
intergenic_df <- intergenic_df %>% filter(sample %in% c("H3K27ac", "H3K9me3", "POL2"))
intergenic_df <- intergenic_df %>% filter(sample %in% c("H3K27ac", "H3K9me3", "H3K9me3_ko"))


# ── 11. Label data ────────────────────────────────────────────────────────────
label_de <- intergenic_df %>%
  dplyr::distinct(upstream_gene, de_status) %>%
  dplyr::group_by(de_status) %>%
  dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n_genes))

label_gap <- intergenic_df %>%
  dplyr::distinct(interval_id, gap_group) %>%
  dplyr::group_by(gap_group) %>%
  dplyr::summarise(n_intervals = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n_intervals))

# ── 12. Scatter plot ──────────────────────────────────────────────────────────
p_scatter <- ggplot(intergenic_df,
                    aes(x = gap_width, y = mean_signal, color = sample)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = sample_colors) +
  facet_grid(sample ~ de_status, scales = "free_y") +
  labs(x        = "Intergenic distance (bp, log scale)",
       y        = "Mean peak signalValue",
       title    = "Intergenic distance vs ChIP-seq peak signal",
       subtitle = "Stratified by upstream gene DE status") +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "none",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

# ── 13. Cumulative signal plot ────────────────────────────────────────────────
intergenic_cum <- intergenic_df %>%
  dplyr::group_by(sample, de_status) %>%
  dplyr::arrange(gap_width) %>%
  dplyr::mutate(
    cum_signal = cumsum(mean_signal) / sum(mean_signal)
  ) %>%
  dplyr::ungroup()

p_cumulative <- ggplot(intergenic_cum,
                       aes(x = gap_width, y = cum_signal, color = sample)) +
  geom_line(linewidth = 0.9) +
  geom_text(data        = label_de,
            aes(x = Inf, y = 0, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = -0.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_x_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_color_manual(values = sample_colors) +
  facet_wrap(~ de_status, nrow = 1) +
  labs(x        = "Intergenic distance (bp, log10)",
       y        = "Cumulative fraction of total signal",
       colour   = "ChIP mark",
       title    = "Cumulative ChIP-seq signal over intergenic distances",
       subtitle = "Stratified by upstream gene DE status") +
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

# ── 14. Density plot ──────────────────────────────────────────────────────────
p_density <- ggplot(intergenic_df,
                    aes(x = gap_width, colour = sample)) +
  geom_density(linewidth = 0.8, adjust = 1.5) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_x_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_colour_manual(values = sample_colors) +
  facet_wrap(~ de_status, nrow = 1) +
  labs(
    x        = "Intergenic distance (bp, log10)",
    y        = "Density",
    colour   = "ChIP mark",
    title    = "Intergenic distance distribution by ChIP mark",
    subtitle = "Stratified by upstream gene DE status"
  ) +
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

# ── 15. Gap group plot ────────────────────────────────────────────────────────
p_gap_groups <- ggplot(intergenic_df,
                       aes(x = de_status, y = mean_signal, fill = sample)) +
  geom_violin(alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               fill = "white", color = "grey30") +
  geom_text(data        = label_gap,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_fill_manual(values = sample_colors) +
  facet_grid(sample ~ gap_group, scales = "free_y") +
  labs(x        = "Upstream gene DE status",
       y        = "Mean peak signalValue",
       title    = "ChIP-seq signal in intergenic regions",
       subtitle = "Stratified by intergenic distance and upstream gene DE status") +
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

# ── 16. Print all plots ───────────────────────────────────────────────────────
print(p_scatter)
print(p_cumulative)
print(p_density)
print(p_gap_groups)

# ── 17. updated dens plot ───────────────────────────────────────────────────────
# Compute median per sample per DE status
median_df <- intergenic_df %>%
  dplyr::group_by(sample, de_status) %>%
  dplyr::summarise(
    median_width = median(gap_width, na.rm = TRUE),
    .groups = "drop"
  )

p_density <- ggplot(intergenic_df,
                    aes(x = gap_width, colour = sample)) +
  geom_density(linewidth = 0.8, adjust = 1.5) +
  geom_vline(data      = median_df,
             aes(xintercept = median_width, colour = sample),
             linetype  = "dashed",
             linewidth = 0.6) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_x_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_colour_manual(values = sample_colors) +
  facet_wrap(~ de_status, nrow = 1) +
  labs(
    x        = "Intergenic distance (bp, log10)",
    y        = "Density",
    colour   = "ChIP mark",
    title    = "Intergenic distance distribution by ChIP mark",
    subtitle = "Dashed lines show median per mark — stratified by upstream gene DE status"
  ) +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

print(p_density)

# ── 17. updated dens plot with quaretiles ───────────────────────────────────────────────────────

library(purrr)
unique(intergenic_df$sample)
intergenic_df_quart <- intergenic_df %>% filter(intergenic_df$sample == "H3K27ac" | intergenic_df$sample == "H3K9me3")

sample_colors <- c(
  H3K27ac    = "#D73027",
  H3K9me3    = "#4575B4"
  #POL2       = "#1A9850",
  #H3K27ac_ko = "#F46D43",
  #H3K9me3_ko = "#74ADD1"
)

library(purrr)

# ── Compute density and quartile fills per sample per DE status ───────────────
density_df <- intergenic_df_quart  %>%
  dplyr::group_by(sample, de_status) %>%
  dplyr::group_modify(~ {
    d   <- density(log10(.x$gap_width), adjust = 1.5, n = 512)
    q25 <- log10(quantile(.x$gap_width, 0.25, na.rm = TRUE))
    q50 <- log10(quantile(.x$gap_width, 0.50, na.rm = TRUE))
    q75 <- log10(quantile(.x$gap_width, 0.75, na.rm = TRUE))
    dplyr::tibble(
      x    = d$x,
      y    = d$y,
      q25  = q25,
      q50  = q50,
      q75  = q75,
      zone = dplyr::case_when(
        x < q25             ~ "Q1",
        x >= q25 & x < q50 ~ "Q2",
        x >= q50 & x < q75 ~ "Q3",
        x >= q75            ~ "Q4"
      )
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(zone_key = paste(sample, zone, sep = "_"))

# ── Quartile color palettes ───────────────────────────────────────────────────
quartile_colors <- list(
  #H3K27ac    = c(Q1 = "#FDDBC7", Q2 = "#F4A582", Q3 = "#D6604D", Q4 = "#B2182B"),
  H3K27ac    = c(Q1 = "white", Q2 = "white", Q3 = "white", Q4 = "white"),
  H3K9me3    = c(Q1 = "#D1E5F0", Q2 = "#92C5DE", Q3 = "#4393C3", Q4 = "#2166AC")
  #POL2       = c(Q1 = "#D9F0D3", Q2 = "#A6DBA0", Q3 = "#5AAE61", Q4 = "#1B7837"),
  #H3K27ac_ko = c(Q1 = "#FEE8C8", Q2 = "#FDBB84", Q3 = "#E34A33", Q4 = "#B30000"),
  #H3K9me3_ko = c(Q1 = "#EDF8FB", Q2 = "#B2E2E2", Q3 = "#66C2A4", Q4 = "#238B45")
)

color_lookup <- purrr::imap_dfr(quartile_colors, function(cols, samp) {
  dplyr::tibble(
    sample   = samp,
    zone     = names(cols),
    fill     = unname(cols),
    zone_key = paste(samp, names(cols), sep = "_")
  )
})

fill_values <- setNames(color_lookup$fill, color_lookup$zone_key)

# ── Median line data ──────────────────────────────────────────────────────────
median_df <- density_df %>%
  dplyr::distinct(sample, de_status, q50)

# ── Line data — one line per sample per DE status ─────────────────────────────
line_df <- density_df %>%
  dplyr::group_by(sample, de_status, x) %>%
  dplyr::summarise(y = mean(y), .groups = "drop")

# ── Plot ──────────────────────────────────────────────────────────────────────
p_density_quartile <- ggplot(density_df,
                             aes(x = x, y = y,
                                 group = zone_key,
                                 fill  = zone_key)) +
  geom_area(alpha = 0.8, color = NA) +
  geom_line(data      = line_df,
            aes(x = x, y = y, color = sample, group = sample),
            linewidth = 0.8,
            inherit.aes = FALSE) +
  geom_vline(data      = median_df,
             aes(xintercept = q50, color = sample),
             linetype  = "dashed",
             linewidth = 0.6,
             inherit.aes = FALSE) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_color_manual(values = sample_colors) +
  scale_x_continuous(
    breaks = log10(c(100, 1000, 10000,
                      100000, 1000000)),
    labels = c("100bp", "1kb", 
               "10kb", "100kb", "1Mb")
  ) +
  facet_grid(~ de_status, scales = "free_y") +
  labs(
    x        = "Intergenic distance",
    y        = "Density",
    colour   = "ChIP mark",
    title    = "Intergenic distance distribution by ChIP mark",
    subtitle = "Shading shows quartile zones — dashed line shows median"
  ) +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.title       = element_text(size = 13),
    legend.position  = "right",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 11),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

p_density_quartile <- ggplot(density_df,
                             aes(x = x, y = y,
                                 group = zone_key,
                                 fill  = zone_key)) +
  geom_area(alpha = 0.75, color = NA, position = "identity") +
  geom_line(data        = density_df %>%
              dplyr::group_by(sample, de_status, x) %>%
              dplyr::summarise(y = mean(y), .groups = "drop"),
            aes(x = x, y = y, color = sample, group = sample),
            linewidth   = 0.8,
            inherit.aes = FALSE) +
  geom_text(data        = label_de,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3.5, color = "grey40", fontface = "italic") +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_color_manual(values = sample_colors) +
  scale_x_continuous(
    breaks = log10(c(100, 1000, 10000,
                      100000, 1000000)),
    labels = c("100bp", "1kb", 
               "10kb", "100kb", "1Mb")
  ) +
  facet_wrap(~ de_status, nrow = 1) +
  labs(
    x        = "Intergenic distance",
    y        = "Density",
    colour   = "ChIP mark",
    title    = "Intergenic distance distribution",
    subtitle = "H3K27ac vs H3K9me3 — shading shows quartile zones"
  ) +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.title       = element_text(size = 13),
    legend.position  = "right",
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 11),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

print(p_density_quartile)

