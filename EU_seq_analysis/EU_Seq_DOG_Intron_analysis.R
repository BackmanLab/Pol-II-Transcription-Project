##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# EU-seq intronic signal and Downstream of Gene (DoG) analysis.
# All objects defined before first use.
# RPKM_THRESHOLD = 0.5
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(ggdist)
  library(patchwork)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(BiocGenerics)
})

##--------------------------------------------------------------## 1. Parameters

CNT_DIR    <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline"
INTRON_MAT <- file.path(CNT_DIR, "intron_dog/intron_counts/euseq_intron_rpkm.txt")
DOG_MAT    <- file.path(CNT_DIR, "intron_dog/dog_counts/euseq_dog_rpkm.txt")
GENE_MAT   <- file.path(CNT_DIR, "counts/euseq_rpkm_analysis.txt")
GTF_PATH   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"

LFC_THRESHOLD  <- 1
RPKM_THRESHOLD <- 0.5
DOG_THRESHOLD  <- NULL   ## set after inspecting DoG distribution
DOG_WINDOW     <- 5000
N_LABELS       <- 20

COL_UP   <- "#D55E00"
COL_DOWN <- "#0072B2"
COL_NS   <- "grey70"

EI_COLORS <- c(
  "Low E/I"    = "#0072B2",
  "Medium E/I" = "#E69F00",
  "High E/I"   = "#CC79A7"
)

POLR2A_SAMPLES <- c("POLR2A_1", "POLR2A_2", "POLR2A_3")
WT_SAMPLES     <- c("WT_1", "WT_2", "WT_3")
ALL_SAMPLES    <- c(POLR2A_SAMPLES, WT_SAMPLES)

CONDITION_MAP <- c(
  POLR2A_1 = "POLR2A", POLR2A_2 = "POLR2A", POLR2A_3 = "POLR2A",
  WT_1 = "WT", WT_2 = "WT", WT_3 = "WT"
)
CONDITION_COLORS <- c(POLR2A = "#7B2D8B", WT = "#1B7837")

TOTAL_RNASEQ_SF <- c(
  POLR2A_1 = 0.747437, POLR2A_2 = 0.747437, POLR2A_3 = 0.747437,
  WT_1 = 1.000000, WT_2 = 1.000000, WT_3 = 1.000000
)

##--------------------------------------------------------------## 2. Helper functions

to_matrix <- function(df) {
  mat <- as.matrix(df[, -1])
  rownames(mat) <- sub("\\..*", "", df[[1]])
  mode(mat) <- "numeric"
  mat
}

apply_sf <- function(mat, sf) {
  sweep(mat, 2, sf[colnames(mat)], FUN = "*")
}

zscore_analysis <- function(mat_scaled, threshold, label) {
  mat_scaled <- mat_scaled[rowSums(mat_scaled) > 0, ]
  wt_mean    <- rowMeans(mat_scaled[, WT_SAMPLES])
  keep       <- wt_mean >= threshold
  cat(sprintf("%s: %d genes passing threshold\n", label, sum(keep)))
  if (sum(keep) == 0) {
    warning("No genes passed threshold for ", label)
    return(data.frame(gene_id = character(0), mean_rpkm = numeric(0),
                      mean_polr2a = numeric(0), mean_wt = numeric(0),
                      lfc = numeric(0), pooled_sd = numeric(0),
                      z_score = numeric(0),
                      de_status = factor(character(0),
                                         levels = c("Up","Not DE","Down"))))
  }
  mat_filt    <- log2(mat_scaled[keep, ] + 1)
  polr2a_vals <- mat_filt[, POLR2A_SAMPLES]
  wt_vals     <- mat_filt[, WT_SAMPLES]
  mean_p      <- rowMeans(polr2a_vals)
  mean_w      <- rowMeans(wt_vals)
  lfc         <- mean_p - mean_w
  sd_pool     <- apply(cbind(polr2a_vals, wt_vals), 1, sd)
  data.frame(
    gene_id     = rownames(mat_filt),
    mean_rpkm   = rowMeans(mat_scaled[keep, ]),
    mean_polr2a = mean_p,
    mean_wt     = mean_w,
    lfc         = lfc,
    pooled_sd   = sd_pool,
    z_score     = lfc / (sd_pool + 1e-10),
    de_status   = factor(dplyr::case_when(
      lfc >  LFC_THRESHOLD ~ "Up",
      lfc < -LFC_THRESHOLD ~ "Down",
      TRUE                 ~ "Not DE"),
      levels = c("Up","Not DE","Down")),
    stringsAsFactors = FALSE
  )
}

base_theme <- theme_classic(base_size = 13) +
  theme(legend.position = "top", legend.title = element_blank(),
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey40", size = 11),
        axis.line     = element_line(color = "grey40"),
        axis.ticks    = element_line(color = "grey40"))

de_colors <- c(Up = COL_UP, Down = COL_DOWN, "Not DE" = COL_NS)

##--------------------------------------------------------------## 3. Load and scale matrices

message("Loading matrices...")
intron_mat <- to_matrix(read_tsv(INTRON_MAT, show_col_types = FALSE))
dog_mat    <- to_matrix(read_tsv(DOG_MAT,    show_col_types = FALSE))
gene_mat   <- to_matrix(read_tsv(GENE_MAT,   show_col_types = FALSE))

intron_scaled <- apply_sf(intron_mat, TOTAL_RNASEQ_SF)
dog_scaled    <- apply_sf(dog_mat,    TOTAL_RNASEQ_SF)
gene_scaled   <- apply_sf(gene_mat,   TOTAL_RNASEQ_SF)

## Exonic RPKM = gene body - intron, floor at 0
common_gi   <- intersect(rownames(gene_scaled), rownames(intron_scaled))
exon_scaled <- gene_scaled[common_gi, ] - intron_scaled[common_gi, ]
exon_scaled[exon_scaled < 0] <- 0

cat("Intron matrix:    ", nrow(intron_scaled), "x", ncol(intron_scaled), "\n")
cat("DoG matrix:       ", nrow(dog_scaled),    "x", ncol(dog_scaled),    "\n")
cat("Gene body matrix: ", nrow(gene_scaled),   "x", ncol(gene_scaled),   "\n")
cat("Exon matrix:      ", nrow(exon_scaled),   "x", ncol(exon_scaled),   "\n")

##--------------------------------------------------------------## 4. Compute genomic E/I groups from GTF
## exon length / intron length per gene — gene architecture
## Completely independent of nascent transcription signal

message("Computing genomic E/I groups from GTF...")
txdb           <- GenomicFeatures::makeTxDbFromGFF(GTF_PATH, format = "gtf",
                                                   taxonomyId = 9606)
exons_by_gene  <- GenomicRanges::reduce(
  GenomicFeatures::exonsBy(txdb, by = "gene"))
exon_lengths   <- sum(BiocGenerics::width(exons_by_gene))
genes_gr       <- GenomicFeatures::genes(txdb)
gene_lengths_g <- setNames(BiocGenerics::width(genes_gr), genes_gr$gene_id)
gtf_common     <- intersect(names(exon_lengths), names(gene_lengths_g))
intron_len_g   <- pmax(gene_lengths_g[gtf_common] - exon_lengths[gtf_common], 0)
ei_genomic     <- exon_lengths[gtf_common] / (intron_len_g + 1)

genomic_ei_df <- data.frame(
  gene_id          = sub("\\..*", "", gtf_common),
  ei_genomic_ratio = as.numeric(ei_genomic),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(ei_genomic_ratio), ei_genomic_ratio > 0) %>%
  dplyr::mutate(
    genomic_ei_group = dplyr::case_when(
      ei_genomic_ratio < quantile(ei_genomic_ratio, 1/3, na.rm = TRUE) ~ "Low E/I",
      ei_genomic_ratio < quantile(ei_genomic_ratio, 2/3, na.rm = TRUE) ~ "Medium E/I",
      TRUE                                                               ~ "High E/I"
    ),
    genomic_ei_group = factor(genomic_ei_group,
                              levels = c("Low E/I","Medium E/I","High E/I"))
  )

cat("Genomic E/I groups: ", nrow(genomic_ei_df), "genes\n")
print(table(genomic_ei_df$genomic_ei_group))

##--------------------------------------------------------------## 5. Nascent E/I ratio
## exon RPKM / intron RPKM — reflects co-transcriptional splicing
## Filter by gene body RPKM threshold rather than requiring both
## exon and intron signal — avoids filtering out High E/I genes
## which have little intronic signal by definition

ei_nascent_common <- intersect(rownames(intron_scaled), rownames(exon_scaled))
ei_nascent_common <- intersect(ei_nascent_common, rownames(gene_scaled))

nascent_ei_df <- data.frame(
  gene_id       = ei_nascent_common,
  intron_wt     = rowMeans(intron_scaled[ei_nascent_common, WT_SAMPLES]),
  exon_wt       = rowMeans(exon_scaled[ei_nascent_common,   WT_SAMPLES]),
  intron_polr2a = rowMeans(intron_scaled[ei_nascent_common, POLR2A_SAMPLES]),
  exon_polr2a   = rowMeans(exon_scaled[ei_nascent_common,   POLR2A_SAMPLES]),
  gene_wt       = rowMeans(gene_scaled[ei_nascent_common,   WT_SAMPLES]),
  stringsAsFactors = FALSE
) %>%
  ## Filter by gene body RPKM — not by exon/intron individually
  ## Using gene body ensures High E/I genes aren't excluded
  dplyr::filter(gene_wt >= RPKM_THRESHOLD) %>%
  dplyr::mutate(
    ## Small pseudocount so ratio is defined even when one
    ## compartment has near-zero signal
    ei_wt     = (exon_wt     + 1e-3) / (intron_wt     + 1e-3),
    ei_polr2a = (exon_polr2a + 1e-3) / (intron_polr2a + 1e-3),
    ei_lfc    = log2(ei_polr2a) - log2(ei_wt)
  ) %>%
  dplyr::filter(!is.na(ei_lfc), !is.infinite(ei_lfc))

cat("Nascent E/I computed for:", nrow(nascent_ei_df), "genes\n")
cat("Summary of ei_lfc:\n")
print(summary(nascent_ei_df$ei_lfc))

##==============================================================##
##                    PART 1 — INTRON ANALYSIS                  ##
##==============================================================##

message("\n========== PART 1: INTRON ANALYSIS ==========\n")

##--------------------------------------------------------------## 6. Intron QC

log2_intron      <- log2(intron_scaled + 1)
log2_intron_sig  <- log2_intron[rowSums(log2_intron) > 0, ]
cat("Genes with any intronic signal:", nrow(log2_intron_sig), "\n")

## QC 1 — Boxplot per sample
box_long <- as.data.frame(log2_intron_sig) %>%
  tibble::rownames_to_column("gene_id") %>%
  tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "log2rpkm") %>%
  dplyr::mutate(condition = factor(CONDITION_MAP[sample],
                                   levels = c("POLR2A","WT")),
                sample    = factor(sample, levels = ALL_SAMPLES))

p_intron_box <- ggplot(box_long, aes(x = sample, y = log2rpkm, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, linewidth = 0.5) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "log2(intronic RPKM + 1)",
       title = "Intronic RPKM per sample",
       subtitle = paste0(nrow(log2_intron_sig), " genes with signal")) +
  base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_intron_box)

## QC 2 — Density
p_intron_dens <- ggplot(box_long %>% dplyr::filter(log2rpkm > 0),
                        aes(x = log2rpkm, color = sample,
                            linetype = condition)) +
  geom_density(linewidth = 0.7, adjust = 1.5) +
  scale_color_manual(values = setNames(
    c(rep("#7B2D8B",3), rep("#1B7837",3)), ALL_SAMPLES)) +
  scale_linetype_manual(values = c(POLR2A = "solid", WT = "solid")) +
  labs(x = "log2(intronic RPKM + 1)", y = "Density",
       title = "Intronic RPKM density per sample") +
  base_theme
print(p_intron_dens)

## QC 3 — Pearson correlation heatmap
keep_qc   <- rowSums(log2_intron >= log2(2)) >= 3
log2_expr <- log2_intron[keep_qc, ]
cor_mat   <- cor(log2_expr, method = "pearson")
cat("Genes for correlation QC:", nrow(log2_expr), "\n")
print(round(cor_mat, 4))

ann_col  <- data.frame(Condition = CONDITION_MAP[colnames(cor_mat)],
                       row.names = colnames(cor_mat))
cor_min  <- floor(min(cor_mat) * 100) / 100
pheatmap(cor_mat,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdYlBu")))(100),
         breaks = seq(cor_min, 1, length.out = 101),
         annotation_col = ann_col, annotation_row = ann_col,
         annotation_colors = list(Condition = CONDITION_COLORS),
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.3f",
         number_color = "grey20", fontsize_number = 8,
         fontsize = 11, cellwidth = 50, cellheight = 50,
         border_color = NA,
         main = paste0("Intronic Pearson correlation\n",
                       nrow(log2_expr), " genes  |  6 samples"))

## QC 4 — PCA
pca_res <- prcomp(t(log2_expr), center = TRUE, scale. = TRUE)
pca_df  <- as.data.frame(pca_res$x[,1:3]) %>%
  dplyr::mutate(sample    = rownames(.),
                condition = factor(CONDITION_MAP[sample],
                                   levels = c("POLR2A","WT")))
pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = paste0("PC1 (",pct_var[1],"%)"), y = paste0("PC2 (",pct_var[2],"%)"),
       title = "Intronic EU-seq PCA") +
  base_theme

p_pca3 <- ggplot(pca_df, aes(x = PC1, y = PC3, color = condition, label = sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = paste0("PC1 (",pct_var[1],"%)"), y = paste0("PC3 (",pct_var[3],"%)"),
       title = "Intronic EU-seq PCA") +
  base_theme

print(p_pca | p_pca3)

## QC 5 — RPKM distribution with threshold line
wt_mean_intron <- rowMeans(intron_scaled[, WT_SAMPLES])
thresh_df <- data.frame(gene_id = rownames(intron_scaled),
                        wt_mean = wt_mean_intron) %>%
  dplyr::filter(wt_mean > 0) %>%
  dplyr::mutate(log2_rpkm = log2(wt_mean + 1))

p_thresh <- ggplot(thresh_df, aes(x = log2_rpkm)) +
  geom_histogram(bins = 100, fill = "#1B7837", color = NA, alpha = 0.7) +
  geom_vline(xintercept = log2(RPKM_THRESHOLD + 1),
             linetype = "dashed", color = "#D55E00", linewidth = 0.7) +
  annotate("text", x = log2(RPKM_THRESHOLD + 1) + 0.05, y = Inf,
           vjust = 1.5, hjust = 0,
           label = paste0("Threshold = ", RPKM_THRESHOLD),
           color = "#D55E00", size = 3.5) +
  labs(x = "log2(mean WT intronic RPKM + 1)", y = "Number of genes",
       title    = "Intronic RPKM distribution",
       subtitle = paste0("Threshold = ", RPKM_THRESHOLD,
                         "  |  Genes passing: ",
                         sum(wt_mean_intron >= RPKM_THRESHOLD))) +
  base_theme
print(p_thresh)

##--------------------------------------------------------------## 7. Fold-change analysis

message("Running fold-change analysis...")
intron_results <- zscore_analysis(intron_scaled, RPKM_THRESHOLD, "Intronic")
gene_results   <- zscore_analysis(gene_scaled,   RPKM_THRESHOLD, "Gene body")
exon_results   <- zscore_analysis(exon_scaled,   RPKM_THRESHOLD, "Exonic")

cat("\nIntronic DE counts:\n"); print(table(intron_results$de_status))
cat("Gene body DE counts:\n"); print(table(gene_results$de_status))

##--------------------------------------------------------------## 8. Intron plots

## --- Plot I-1: Intronic LFC distribution ---
p_intron_lfc <- ggplot(intron_results, aes(x = lfc, fill = de_status)) +
  geom_histogram(bins = 100, color = NA, alpha = 0.8) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  scale_fill_manual(values = de_colors) +
  labs(x = "log2 FC (POLR2A / WT)", y = "Number of genes",
       title    = "Intronic LFC distribution",
       subtitle = paste0(nrow(intron_results), " genes  |  |LFC| > ",
                         LFC_THRESHOLD)) +
  base_theme + theme(legend.title = element_blank())
print(p_intron_lfc)

## --- Plot I-2: Intronic vs gene body RPKM scatter (WT) ---
ig_common <- intersect(rownames(intron_scaled), rownames(gene_scaled))
ig_df <- data.frame(
  gene_id   = ig_common,
  intron_wt = rowMeans(intron_scaled[ig_common, WT_SAMPLES]),
  gene_wt   = rowMeans(gene_scaled[ig_common,   WT_SAMPLES]),
  intron_p  = rowMeans(intron_scaled[ig_common, POLR2A_SAMPLES]),
  gene_p    = rowMeans(gene_scaled[ig_common,   POLR2A_SAMPLES]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(gene_wt >= RPKM_THRESHOLD, gene_p >= RPKM_THRESHOLD) %>%
  dplyr::mutate(
    intron_frac_wt     = intron_wt / gene_wt,
    intron_frac_polr2a = intron_p  / gene_p
  )

p_ig_scatter <- ggplot(ig_df, aes(x = log2(gene_wt+1), y = log2(intron_wt+1))) +
  geom_point(size = 0.6, alpha = 0.3, color = "#1B7837") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  labs(x = "log2(gene body RPKM + 1) — WT",
       y = "log2(intronic RPKM + 1) — WT",
       title    = "Intronic vs gene body RPKM — WT",
       subtitle = paste0("Diagonal = all signal is intronic  |  n=", nrow(ig_df))) +
  base_theme
print(p_ig_scatter)

## --- Plot I-3a: Intronic fraction by genomic E/I group (paired Wilcoxon) ---
intron_frac_base <- ig_df %>%
  dplyr::select(gene_id, intron_frac_wt, intron_frac_polr2a) %>%
  dplyr::filter(intron_frac_wt > 0, intron_frac_wt <= 1,
                intron_frac_polr2a >= 0, intron_frac_polr2a <= 1)

intron_frac_ei <- intron_frac_base %>%
  dplyr::left_join(genomic_ei_df %>% dplyr::select(gene_id, genomic_ei_group),
                   by = "gene_id") %>%
  dplyr::filter(!is.na(genomic_ei_group)) %>%
  tidyr::pivot_longer(cols = c(intron_frac_wt, intron_frac_polr2a),
                      names_to = "condition", values_to = "intron_frac") %>%
  dplyr::mutate(
    condition        = dplyr::recode(condition,
                                     intron_frac_wt     = "WT",
                                     intron_frac_polr2a = "POLR2A"),
    condition        = factor(condition, levels = c("WT","POLR2A")),
    genomic_ei_group = factor(genomic_ei_group,
                              levels = c("Low E/I","Medium E/I","High E/I"))
  )

ei_stats_paired <- intron_frac_base %>%
  dplyr::left_join(genomic_ei_df %>% dplyr::select(gene_id, genomic_ei_group),
                   by = "gene_id") %>%
  dplyr::filter(!is.na(genomic_ei_group)) %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(
    p_value   = wilcox.test(intron_frac_wt, intron_frac_polr2a,
                            paired = TRUE)$p.value,
    median_wt = median(intron_frac_wt,     na.rm = TRUE),
    median_p  = median(intron_frac_polr2a, na.rm = TRUE),
    n         = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = dplyr::case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                             p_adj < 0.05  ~ "*",   TRUE         ~ "ns"),
    y_pos = max(intron_frac_ei$intron_frac, na.rm = TRUE) * 1.05,
    genomic_ei_group = factor(genomic_ei_group,
                              levels = c("Low E/I","Medium E/I","High E/I"))
  )

cat("\nPaired Wilcoxon — intronic fraction by genomic E/I group:\n")
print(ei_stats_paired %>% dplyr::select(genomic_ei_group, n, median_wt,
                                        median_p, p_adj, label))

p_intron_frac_ei <- ggplot(intron_frac_ei,
                           aes(x = condition, y = intron_frac, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2,
               linewidth = 0.5, width = 0.55) +
  geom_text(data = ei_stats_paired,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "grey20") +
  geom_segment(data = ei_stats_paired,
               aes(x = 1, xend = 2, y = y_pos * 0.98, yend = y_pos * 0.98),
               inherit.aes = FALSE, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = CONDITION_COLORS) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, max(intron_frac_ei$intron_frac,
                                       na.rm = TRUE) * 1.12)) +
  facet_wrap(~ genomic_ei_group, nrow = 1) +
  labs(x = NULL, y = "Intronic fraction (intron / gene body RPKM)",
       title    = "Intronic signal fraction by genomic E/I group",
       subtitle = "Paired Wilcoxon  |  * p<0.05  ** p<0.01  *** p<0.001") +
  base_theme + theme(legend.position = "none")
print(p_intron_frac_ei)

## --- Plot I-3b: Intronic fraction by gene body RPKM quartile (paired Wilcoxon) ---
intron_frac_rpkm <- intron_frac_base %>%
  dplyr::left_join(ig_df %>% dplyr::select(gene_id, gene_wt), by = "gene_id") %>%
  dplyr::filter(!is.na(gene_wt)) %>%
  dplyr::mutate(
    rpkm_group = dplyr::case_when(
      gene_wt < quantile(gene_wt, 0.25, na.rm = TRUE) ~ "Q1 (Low)",
      gene_wt < quantile(gene_wt, 0.50, na.rm = TRUE) ~ "Q2",
      gene_wt < quantile(gene_wt, 0.75, na.rm = TRUE) ~ "Q3",
      TRUE                                              ~ "Q4 (High)"
    ),
    rpkm_group = factor(rpkm_group, levels = c("Q1 (Low)","Q2","Q3","Q4 (High)"))
  )

rpkm_stats_paired <- intron_frac_rpkm %>%
  dplyr::group_by(rpkm_group) %>%
  dplyr::summarise(
    p_value   = wilcox.test(intron_frac_wt, intron_frac_polr2a,
                            paired = TRUE)$p.value,
    median_wt = median(intron_frac_wt,     na.rm = TRUE),
    median_p  = median(intron_frac_polr2a, na.rm = TRUE),
    n         = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = dplyr::case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                             p_adj < 0.05  ~ "*",   TRUE         ~ "ns")
  )

intron_frac_rpkm_long <- intron_frac_rpkm %>%
  tidyr::pivot_longer(cols = c(intron_frac_wt, intron_frac_polr2a),
                      names_to = "condition", values_to = "intron_frac") %>%
  dplyr::mutate(
    condition = dplyr::recode(condition,
                              intron_frac_wt     = "WT",
                              intron_frac_polr2a = "POLR2A"),
    condition = factor(condition, levels = c("WT","POLR2A"))
  )

rpkm_stats_paired <- rpkm_stats_paired %>%
  dplyr::mutate(y_pos = max(intron_frac_rpkm_long$intron_frac, na.rm=TRUE) * 1.05)

cat("\nPaired Wilcoxon — intronic fraction by RPKM quartile:\n")
print(rpkm_stats_paired %>% dplyr::select(rpkm_group, n, median_wt,
                                          median_p, p_adj, label))

p_intron_frac_rpkm <- ggplot(intron_frac_rpkm_long,
                             aes(x = condition, y = intron_frac, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2,
               linewidth = 0.5, width = 0.55) +
  geom_text(data = rpkm_stats_paired,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "grey20") +
  geom_segment(data = rpkm_stats_paired,
               aes(x = 1, xend = 2, y = y_pos * 0.98, yend = y_pos * 0.98),
               inherit.aes = FALSE, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = CONDITION_COLORS) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, max(intron_frac_rpkm_long$intron_frac,
                                       na.rm=TRUE) * 1.12)) +
  facet_wrap(~ rpkm_group, nrow = 1) +
  labs(x = NULL, y = "Intronic fraction (intron / gene body RPKM)",
       title    = "Intronic signal fraction by gene body RPKM quartile",
       subtitle = "Paired Wilcoxon  |  * p<0.05  ** p<0.01  *** p<0.001") +
  base_theme + theme(legend.position = "none")
print(p_intron_frac_rpkm)

## --- Plot I-4: Intronic LFC vs gene body LFC scatter ---
intron_gene_df <- intron_results %>%
  dplyr::select(gene_id, lfc_intron = lfc) %>%
  dplyr::inner_join(gene_results %>% dplyr::select(gene_id, lfc_gene = lfc),
                    by = "gene_id") %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      lfc_intron >  LFC_THRESHOLD & lfc_gene >  LFC_THRESHOLD ~ "Both Up",
      lfc_intron < -LFC_THRESHOLD & lfc_gene < -LFC_THRESHOLD ~ "Both Down",
      lfc_intron >  LFC_THRESHOLD & lfc_gene < -LFC_THRESHOLD ~ "Intron Up\nGene Down",
      lfc_intron < -LFC_THRESHOLD & lfc_gene >  LFC_THRESHOLD ~ "Intron Down\nGene Up",
      TRUE                                                      ~ "Not DE"
    ),
    direction = factor(direction,
                       levels = c("Both Up","Both Down","Intron Up\nGene Down",
                                  "Intron Down\nGene Up","Not DE"))
  )

r_ig <- cor(intron_gene_df$lfc_intron, intron_gene_df$lfc_gene,
            use = "complete.obs")
cat("Pearson r (intronic vs gene body LFC):", round(r_ig, 3), "\n")

dir_colors <- c("Both Up" = COL_UP, "Both Down" = COL_DOWN,
                "Intron Up\nGene Down" = "#CC79A7",
                "Intron Down\nGene Up" = "#009E73", "Not DE" = COL_NS)

p_intron_gene <- ggplot(intron_gene_df,
                        aes(x = lfc_gene, y = lfc_intron, color = direction)) +
  geom_point(size = 0.7, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = dir_colors) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(r_ig, 3)),
           hjust = 1.1, vjust = -0.5, size = 4,
           color = "grey30", fontface = "italic") +
  labs(x = "Gene body log2 FC (POLR2A / WT)",
       y = "Intronic log2 FC (POLR2A / WT)",
       title    = "Intronic vs gene body LFC",
       subtitle = paste0("Diagonal = proportional change  |  r = ",
                         round(r_ig, 3))) +
  base_theme
print(p_intron_gene)

## --- Plot I-4b/c: Intronic fraction and E/I ratio as a function of gene length ---
## Build per-gene data frame joining nascent signal with genomic gene length

gene_length_df <- data.frame(
  gene_id    = sub("\\..*", "", names(gene_lengths_g)),
  gene_length = as.numeric(gene_lengths_g),
  stringsAsFactors = FALSE
)

length_signal_df <- ig_df %>%
  dplyr::left_join(gene_length_df, by = "gene_id") %>%
  dplyr::filter(!is.na(gene_length), gene_length > 0) %>%
  dplyr::left_join(
    nascent_ei_df %>% dplyr::select(gene_id, ei_wt, ei_polr2a),
    by = "gene_id"
  ) %>%
  dplyr::filter(!is.na(ei_wt)) %>%
  dplyr::mutate(
    log10_length = log10(gene_length),
    length_rank  = rank(gene_length, ties.method = "average") / n()
  )

cat("Genes for length vs signal plots:", nrow(length_signal_df), "\n")
cat("Gene length range:", range(length_signal_df$gene_length), "\n")

## Pivot to long format for both conditions
length_long <- length_signal_df %>%
  dplyr::select(gene_id, gene_length, log10_length, length_rank,
                intron_frac_wt, intron_frac_polr2a,
                ei_wt, ei_polr2a) %>%
  tidyr::pivot_longer(
    cols = c(intron_frac_wt, intron_frac_polr2a, ei_wt, ei_polr2a),
    names_to  = "metric_condition",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    condition = dplyr::case_when(
      grepl("wt|_wt",      metric_condition) ~ "WT",
      grepl("polr2a|_p$",  metric_condition) ~ "POLR2A"
    ),
    metric = dplyr::case_when(
      grepl("intron_frac", metric_condition) ~ "Intronic fraction\n(intron / gene body RPKM)",
      grepl("ei",          metric_condition) ~ "E/I ratio\n(exon / intron RPKM)"
    ),
    condition = factor(condition, levels = c("WT", "POLR2A"))
  ) %>%
  dplyr::filter(!is.na(condition), !is.na(metric))

## Bin genes into 50 equally-sized bins for mean line
## Two versions: by gene length rank and by log10 gene length

bin_means <- function(df, x_var, n_bins = 50) {
  df %>%
    dplyr::mutate(bin = cut(!!rlang::sym(x_var),
                            breaks = n_bins, labels = FALSE)) %>%
    dplyr::group_by(bin, condition, metric) %>%
    dplyr::summarise(
      mean_value = mean(value, na.rm = TRUE),
      mean_x     = mean(!!rlang::sym(x_var), na.rm = TRUE),
      n          = dplyr::n(),
      .groups    = "drop"
    )
}

bins_rank   <- bin_means(length_long, "length_rank",  n_bins = 50)
bins_log10  <- bin_means(length_long, "log10_length", n_bins = 50)

## --- Plot I-4b: X = gene length rank (0 to 1), faceted by metric ---
p_length_rank <- ggplot(bins_rank,
                        aes(x = mean_x, y = mean_value,
                            color = condition)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = CONDITION_COLORS) +
  facet_wrap(~ metric, nrow = 2, scales = "free_y") +
  labs(
    x        = "Gene length rank (shortest → longest)",
    y        = "Mean value per bin",
    title    = "Nascent signal vs gene length rank",
    subtitle = paste0("50 equal-size bins  |  n=", nrow(length_signal_df),
                      " genes  |  mean across replicates")
  ) +
  base_theme

print(p_length_rank)

## --- Plot I-4c: X = log10(gene length bp), faceted by metric ---
p_length_log10 <- ggplot(bins_log10,
                         aes(x = mean_x, y = mean_value,
                             color = condition)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = CONDITION_COLORS) +
  scale_x_continuous(
    breaks = log10(c(1e3, 1e4, 1e5, 1e6)),
    labels = c("1kb", "10kb", "100kb", "1Mb")
  ) +
  facet_wrap(~ metric, nrow = 2, scales = "free_y") +
  labs(
    x        = "Gene length — log10 scale",
    y        = "Mean value per bin",
    title    = "Nascent signal vs gene length (log10 scale)",
    subtitle = paste0("50 equal-size bins  |  n=", nrow(length_signal_df),
                      " genes  |  mean across replicates")
  ) +
  base_theme

print(p_length_log10)


## ggdist half-violin + boxplot, faceted by genomic E/I group

ei_dist_df <- nascent_ei_df %>%
  dplyr::select(gene_id, ei_lfc) %>%
  dplyr::filter(!is.na(ei_lfc), !is.infinite(ei_lfc)) %>%
  dplyr::left_join(genomic_ei_df %>% dplyr::select(gene_id, genomic_ei_group),
                   by = "gene_id") %>%
  dplyr::filter(!is.na(genomic_ei_group)) %>%
  dplyr::mutate(
    genomic_ei_group = factor(genomic_ei_group,
                              levels = c("Low E/I","Medium E/I","High E/I"))
  )

cat("Genes in ei_dist_df:", nrow(ei_dist_df), "\n")
print(table(ei_dist_df$genomic_ei_group))

## One-sample Wilcoxon vs zero per group
ei_lfc_stats <- ei_dist_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(
    p_value    = wilcox.test(ei_lfc, mu = 0)$p.value,
    median_lfc = median(ei_lfc, na.rm = TRUE),
    n          = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = dplyr::case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                             p_adj < 0.05  ~ "*",   TRUE         ~ "ns"),
    y_pos = 2.7
  )

cat("\nOne-sample Wilcoxon vs zero — nascent E/I LFC:\n")
print(ei_lfc_stats %>% dplyr::select(genomic_ei_group, n, median_lfc,
                                     p_adj, label))

p_ei_dist <- ggplot(ei_dist_df,
                    aes(x = genomic_ei_group, y = ei_lfc,
                        fill = genomic_ei_group, color = genomic_ei_group)) +
  ## Half violin on left
  ggdist::stat_halfeye(adjust = 1.2, width = 0.5, .width = 0,
                       point_color = NA, alpha = 0.7, side = "left") +
  ## Boxplot on right
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "grey20",
               alpha = 0.3, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_text(data = ei_lfc_stats,
            aes(x = genomic_ei_group, y = y_pos,
                label = paste0(label, "\nn=", n)),
            inherit.aes = FALSE, size = 3.5, color = "grey20", vjust = 0) +
  scale_fill_manual(values  = EI_COLORS) +
  scale_color_manual(values = EI_COLORS) +
  scale_x_discrete(labels = c("Low E/I"    = "Low E/I\n(intron-rich)",
                              "Medium E/I" = "Medium E/I",
                              "High E/I"   = "High E/I\n(exon-rich)")) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "Genomic E/I group (exon length / intron length)",
       y = "log2 FC of nascent\nexon RPKM / intron RPKM\n(POLR2A / WT)",
       title    = "Change in nascent exon/intron RPKM ratio after POLR2A degradation",
       subtitle = paste0("Positive = more exonic relative to intronic in POLR2A  |  ",
                         "One-sample Wilcoxon vs zero  |  y clipped ±3")) +
  base_theme + theme(legend.position = "none")
print(p_ei_dist)

## --- Plot I-6a: ei_lfc distribution by genomic E/I group ---
## Top panel: overlapping density lines (normal orientation)
## Bottom panel: vertical boxplots aligned below
## Median marked with color-matched dashed line

EI_LINE_COLORS <- c(
  "Low E/I"    = "#0072B2",
  "Medium E/I" = "#E69F00",
  "High E/I"   = "#CC79A7"
)

ei_plot_df <- nascent_ei_df %>%
  dplyr::select(gene_id, ei_lfc) %>%
  dplyr::filter(!is.na(ei_lfc), !is.infinite(ei_lfc),
                ei_lfc >= -3, ei_lfc <= 3) %>%
  dplyr::left_join(
    genomic_ei_df %>% dplyr::select(gene_id, genomic_ei_group),
    by = "gene_id"
  ) %>%
  dplyr::filter(!is.na(genomic_ei_group)) %>%
  dplyr::mutate(
    genomic_ei_group = factor(genomic_ei_group,
                              levels = c("Low E/I","Medium E/I","High E/I"))
  )

ei_medians <- ei_plot_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(median_lfc = median(ei_lfc, na.rm = TRUE),
                   n          = dplyr::n(),
                   .groups    = "drop")

## One-sample Wilcoxon for significance labels on boxplot
ei_wilcox <- ei_plot_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(
    p_value = wilcox.test(ei_lfc, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    )
  ) %>%
  dplyr::left_join(ei_medians, by = "genomic_ei_group") %>%
  dplyr::mutate(y_pos = 2.8)

## Top panel — density lines, normal orientation
p_ei_dens <- ggplot(ei_plot_df,
                    aes(x = ei_lfc, color = genomic_ei_group)) +
  geom_density(linewidth = 0.8, adjust = 1.2) +
  geom_vline(data = ei_medians,
             aes(xintercept = median_lfc, color = genomic_ei_group),
             linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  scale_color_manual(values = EI_LINE_COLORS) +
  xlim(-3, 3) +
  labs(x = NULL, y = "Density",
       title = "Change in nascent exon/intron RPKM ratio") +
  base_theme +
  theme(legend.position  = "top",
        plot.title       = element_text(face = "bold", size = 13))

## Bottom panel — vertical boxplots
p_ei_box <- ggplot(ei_plot_df,
                   aes(x = genomic_ei_group, y = ei_lfc,
                       color = genomic_ei_group)) +
  geom_boxplot(fill = NA, outlier.shape = NA,
               linewidth = 0.6, width = 0.5) +
  geom_hline(yintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  geom_text(data = ei_wilcox,
            aes(x = genomic_ei_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, color = "grey20") +
  scale_color_manual(values = EI_LINE_COLORS) +
  scale_x_discrete(labels = c("Low E/I"    = "Low E/I\n(intron-rich)",
                              "Medium E/I" = "Medium E/I",
                              "High E/I"   = "High E/I\n(exon-rich)")) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "Genomic E/I group",
       y = "log2 FC of nascent\nexon/intron RPKM\n(POLR2A / WT)") +
  base_theme +
  theme(legend.position = "none")

print(p_ei_dens / p_ei_box +
        patchwork::plot_layout(heights = c(2, 1)))

## --- Plot I-6a-intron: Intronic RPKM LFC by genomic E/I group ---
## Same layout as I-6a but Y = log2 FC of intronic RPKM (POLR2A / WT)
## Shows intron rate of change independent of exon influence

intron_lfc_ei_df <- intron_results %>%
  dplyr::select(gene_id, intron_lfc = lfc) %>%
  dplyr::filter(!is.na(intron_lfc), !is.infinite(intron_lfc),
                intron_lfc >= -3, intron_lfc <= 3) %>%
  dplyr::left_join(
    genomic_ei_df %>% dplyr::select(gene_id, genomic_ei_group),
    by = "gene_id"
  ) %>%
  dplyr::filter(!is.na(genomic_ei_group)) %>%
  dplyr::mutate(
    genomic_ei_group = factor(genomic_ei_group,
                              levels = c("Low E/I","Medium E/I","High E/I"))
  )

intron_lfc_medians <- intron_lfc_ei_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(median_lfc = median(intron_lfc, na.rm = TRUE),
                   .groups    = "drop")

intron_lfc_wilcox <- intron_lfc_ei_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(
    p_value = wilcox.test(intron_lfc, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = 2.8
  )

cat("\nIntronic LFC by genomic E/I — medians and significance:\n")
print(dplyr::left_join(intron_lfc_medians, intron_lfc_wilcox,
                       by = "genomic_ei_group") %>%
        dplyr::select(genomic_ei_group, median_lfc, p_adj, label))

p_intron_lfc_dens <- ggplot(intron_lfc_ei_df,
                            aes(x = intron_lfc, color = genomic_ei_group)) +
  geom_density(linewidth = 0.8, adjust = 1.2) +
  geom_vline(data = intron_lfc_medians,
             aes(xintercept = median_lfc, color = genomic_ei_group),
             linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  scale_color_manual(values = EI_LINE_COLORS) +
  xlim(-3, 3) +
  labs(x = NULL, y = "Density",
       title = "Change in nascent intronic RPKM after POLR2A degradation") +
  base_theme +
  theme(legend.position = "top",
        plot.title      = element_text(face = "bold", size = 13))

p_intron_lfc_box <- ggplot(intron_lfc_ei_df,
                           aes(x = genomic_ei_group, y = intron_lfc,
                               color = genomic_ei_group)) +
  geom_boxplot(fill = NA, outlier.shape = NA,
               linewidth = 0.6, width = 0.5) +
  geom_hline(yintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  geom_text(data = intron_lfc_wilcox,
            aes(x = genomic_ei_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, color = "grey20") +
  scale_color_manual(values = EI_LINE_COLORS) +
  scale_x_discrete(labels = c("Low E/I"    = "Low E/I\n(intron-rich)",
                              "Medium E/I" = "Medium E/I",
                              "High E/I"   = "High E/I\n(exon-rich)")) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "Genomic E/I group",
       y = "log2 FC of nascent\nintronic RPKM\n(POLR2A / WT)") +
  base_theme +
  theme(legend.position = "none")

print(p_intron_lfc_dens / p_intron_lfc_box +
        patchwork::plot_layout(heights = c(2, 1)))


GL_COLORS <- c(
  "Q1 (Short)" = "#0072B2",
  "Q2"         = "#56B4E9",
  "Q3"         = "#E69F00",
  "Q4 (Long)"  = "#D55E00"
)

gl_quartile_df <- nascent_ei_df %>%
  dplyr::select(gene_id, ei_lfc) %>%
  dplyr::filter(!is.na(ei_lfc), !is.infinite(ei_lfc),
                ei_lfc >= -3, ei_lfc <= 3) %>%
  dplyr::left_join(gene_length_df, by = "gene_id") %>%
  dplyr::filter(!is.na(gene_length)) %>%
  dplyr::mutate(
    length_group = dplyr::case_when(
      gene_length < quantile(gene_length, 0.25, na.rm = TRUE) ~ "Q1 (Short)",
      gene_length < quantile(gene_length, 0.50, na.rm = TRUE) ~ "Q2",
      gene_length < quantile(gene_length, 0.75, na.rm = TRUE) ~ "Q3",
      TRUE                                                      ~ "Q4 (Long)"
    ),
    length_group = factor(length_group,
                          levels = c("Q1 (Short)","Q2","Q3","Q4 (Long)"))
  )

gl_medians <- gl_quartile_df %>%
  dplyr::group_by(length_group) %>%
  dplyr::summarise(median_lfc = median(ei_lfc, na.rm = TRUE),
                   n          = dplyr::n(), .groups = "drop")

gl_wilcox <- gl_quartile_df %>%
  dplyr::group_by(length_group) %>%
  dplyr::summarise(
    p_value = wilcox.test(ei_lfc, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = 2.8
  )

p_gl_dens <- ggplot(gl_quartile_df,
                    aes(x = ei_lfc, color = length_group)) +
  geom_density(linewidth = 0.8, adjust = 1.2) +
  geom_vline(data = gl_medians,
             aes(xintercept = median_lfc, color = length_group),
             linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  scale_color_manual(values = GL_COLORS) +
  xlim(-3, 3) +
  labs(x = NULL, y = "Density",
       title = "Change in nascent exon/intron RPKM ratio by gene length") +
  base_theme +
  theme(legend.position  = "top",
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        plot.title       = element_text(face = "bold", size = 13))

p_gl_box <- ggplot(gl_quartile_df,
                   aes(x = length_group, y = ei_lfc,
                       color = length_group)) +
  geom_boxplot(fill = NA, outlier.shape = NA,
               linewidth = 0.6, width = 0.5) +
  geom_hline(yintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  geom_text(data = gl_wilcox,
            aes(x = length_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, color = "grey20") +
  scale_color_manual(values = GL_COLORS) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "Gene length quartile",
       y = "log2 FC of nascent\nexon/intron RPKM\n(POLR2A / WT)") +
  base_theme +
  theme(legend.position = "none")

print(p_gl_dens / p_gl_box +
        patchwork::plot_layout(heights = c(2, 1)))

## --- Plot I-6c: Intronic RPKM vs log10 gene length ---
## Two lines: WT and POLR2A
## X = log10(intronic RPKM + 1), Y = log10(gene length)
## Binned by gene length (every 50 genes)

intron_length_df <- data.frame(
  gene_id      = ig_common,
  intron_wt    = rowMeans(intron_scaled[ig_common, WT_SAMPLES]),
  intron_p     = rowMeans(intron_scaled[ig_common, POLR2A_SAMPLES]),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(gene_length_df, by = "gene_id") %>%
  dplyr::filter(!is.na(gene_length), gene_length > 0,
                intron_wt > 0 | intron_p > 0) %>%
  dplyr::mutate(log10_length = log10(gene_length)) %>%
  dplyr::arrange(gene_length)

n_bins_c <- floor(nrow(intron_length_df) / 50)
cat("Plot I-6c bins:", n_bins_c, "\n")

intron_bins <- intron_length_df %>%
  dplyr::mutate(bin = cut(log10_length, breaks = n_bins_c,
                          labels = FALSE)) %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(
    mean_log10_length  = mean(log10_length, na.rm = TRUE),
    mean_intron_wt     = mean(log2(intron_wt + 1), na.rm = TRUE),
    mean_intron_polr2a = mean(log2(intron_p  + 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(mean_intron_wt, mean_intron_polr2a),
                      names_to = "condition", values_to = "mean_log2_rpkm") %>%
  dplyr::mutate(
    condition = dplyr::recode(condition,
                              mean_intron_wt     = "WT",
                              mean_intron_polr2a = "POLR2A"),
    condition = factor(condition, levels = c("WT","POLR2A"))
  )

p_intron_length <- ggplot(intron_bins,
                          aes(x = mean_log10_length,
                              y = mean_log2_rpkm,
                              color = condition)) +
  geom_line(linewidth = 0.4, alpha = 0.3) +
  geom_smooth(method = "loess", span = 0.2,
              se = FALSE, linewidth = 1.0) +
  scale_color_manual(values = CONDITION_COLORS) +
  scale_x_continuous(
    breaks = log10(c(1e3, 1e4, 1e5, 1e6)),
    labels = c("1kb","10kb","100kb","1Mb")
  ) +
  labs(
    x        = "Gene length — log10 scale",
    y        = "log2(intronic RPKM + 1)",
    title    = "Intronic RPKM vs gene length",
    subtitle = paste0("One bin per 50 genes  |  n=",
                      nrow(intron_length_df),
                      " genes  |  loess smoothed (span=0.2)")
  ) +
  base_theme

print(p_intron_length)




##==============================================================##

message("\n========== PART 2: DoG ANALYSIS ==========\n")

##--------------------------------------------------------------## 9. DoG QC

wt_mean_dog <- rowMeans(dog_scaled[, WT_SAMPLES])
dog_dist_df <- data.frame(gene_id = rownames(dog_scaled),
                          wt_mean = wt_mean_dog) %>%
  dplyr::filter(wt_mean > 0) %>%
  dplyr::mutate(log2_rpkm = log2(wt_mean + 1))

cat("Genes with DoG windows:", nrow(dog_scaled), "\n")
cat("Genes with DoG signal > 0:", nrow(dog_dist_df), "\n")

for (thresh in c(0.01, 0.05, 0.1, 0.5)) {
  cat(sprintf("  >= %.2f: %d genes\n", thresh, sum(wt_mean_dog >= thresh)))
}

## QC 1 — Distribution
p_dog_dist <- ggplot(dog_dist_df, aes(x = log2_rpkm)) +
  geom_histogram(bins = 100, fill = "#E69F00", color = NA, alpha = 0.7) +
  geom_vline(xintercept = log2(c(0.01, 0.05, 0.1, 0.5) + 1),
             linetype = "dashed",
             color = c("grey60","#CC79A7","#E69F00","#D55E00"),
             linewidth = 0.7) +
  annotate("text", x = log2(c(0.01, 0.05, 0.1, 0.5) + 1) + 0.02,
           y = Inf, vjust = 1.5, hjust = 0,
           label = c("0.01","0.05","0.1","0.5"),
           color = c("grey60","#CC79A7","#E69F00","#D55E00"), size = 3.5) +
  labs(x = "log2(mean WT DoG RPKM + 1)", y = "Number of genes",
       title    = paste0("DoG RPKM distribution (",
                         DOG_WINDOW/1000, "kb window)"),
       subtitle = paste0("Genes with signal > 0: ", nrow(dog_dist_df),
                         "  |  Set DOG_THRESHOLD then re-run")) +
  base_theme
print(p_dog_dist)

## QC 2 — Density by condition
dog_long <- as.data.frame(dog_scaled) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::filter(rowSums(dplyr::across(all_of(ALL_SAMPLES))) > 0) %>%
  tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "rpkm") %>%
  dplyr::mutate(condition = factor(CONDITION_MAP[sample], levels = c("POLR2A","WT")),
                log2rpkm  = log2(rpkm + 1)) %>%
  dplyr::filter(log2rpkm > 0)

p_dog_dens <- ggplot(dog_long, aes(x = log2rpkm, color = condition, fill = condition)) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.5) +
  scale_color_manual(values = CONDITION_COLORS) +
  scale_fill_manual(values  = CONDITION_COLORS) +
  labs(x = "log2(DoG RPKM + 1)", y = "Density",
       title = "DoG RPKM density by condition") +
  base_theme
print(p_dog_dens)

## QC 3 — Boxplot per sample
dog_box_df <- as.data.frame(log2(dog_scaled + 1)) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::filter(rowSums(dplyr::across(all_of(ALL_SAMPLES))) > 0) %>%
  tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "log2rpkm") %>%
  dplyr::mutate(condition = factor(CONDITION_MAP[sample], levels = c("POLR2A","WT")),
                sample    = factor(sample, levels = ALL_SAMPLES))

p_dog_box <- ggplot(dog_box_df, aes(x = sample, y = log2rpkm, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, linewidth = 0.5) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "log2(DoG RPKM + 1)", title = "DoG RPKM per sample") +
  base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_dog_box)

DOG_THRESHOLD=0.05

if (is.null(DOG_THRESHOLD)) {
  message("\n*** Set DOG_THRESHOLD in Section 1, then re-run from Section 10 ***")
  stop("DOG_THRESHOLD not set.", call. = FALSE)
}

##--------------------------------------------------------------## 10. DoG analysis plots

message("Running DoG analysis...")
dog_results <- zscore_analysis(dog_scaled, DOG_THRESHOLD, "DoG")
cat("\nDoG DE counts:\n"); print(table(dog_results$de_status))

## --- Plot D-1: Detection rate ---
dog_detect <- data.frame(
  gene_id     = rownames(dog_scaled),
  wt_mean     = rowMeans(dog_scaled[, WT_SAMPLES]),
  polr2a_mean = rowMeans(dog_scaled[, POLR2A_SAMPLES]),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    det_wt  = wt_mean     >= DOG_THRESHOLD,
    det_p   = polr2a_mean >= DOG_THRESHOLD,
    category = dplyr::case_when(
      det_wt &  det_p  ~ "Both",
      det_wt & !det_p  ~ "WT only",
      !det_wt & det_p  ~ "POLR2A only",
      TRUE             ~ "Neither"
    )
  )

detect_counts <- dog_detect %>%
  dplyr::count(category) %>%
  dplyr::mutate(pct = n / sum(n) * 100,
                category = factor(category,
                                  levels = c("Both","WT only",
                                             "POLR2A only","Neither")))
cat("\nDoG detection:\n"); print(detect_counts)

p_dog_detect <- ggplot(detect_counts %>% dplyr::filter(category != "Neither"),
                       aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(n, "\n(", round(pct,1), "%)")),
            vjust = -0.3, size = 3.5, color = "grey20") +
  scale_fill_manual(values = c(Both = "#009E73", "WT only" = "#1B7837",
                               "POLR2A only" = "#7B2D8B")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of genes",
       title    = "DoG signal detection by condition",
       subtitle = paste0("DoG RPKM >= ", DOG_THRESHOLD,
                         "  |  ", DOG_WINDOW/1000, "kb window")) +
  base_theme + theme(legend.position = "none")
print(p_dog_detect)

## --- Plot D-2: DoG MA plot ---
p_dog_ma <- ggplot(dog_results, aes(x = log10(mean_rpkm+1), y = lfc,
                                    color = de_status)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = dog_results %>% dplyr::filter(de_status %in% c("Up","Down")) %>%
      dplyr::arrange(desc(abs(lfc))) %>% dplyr::slice_head(n = N_LABELS),
    aes(label = gene_id), size = 2.8, max.overlaps = 20,
    segment.size = 0.3, segment.color = "grey50", show.legend = FALSE) +
  scale_color_manual(values = de_colors) +
  labs(x = "log10(mean DoG RPKM + 1)", y = "log2 FC (POLR2A / WT)",
       title    = paste0("DoG MA plot — ", DOG_WINDOW/1000, "kb window"),
       subtitle = paste0("|LFC| > ", LFC_THRESHOLD, "  |  Top ",
                         N_LABELS, " genes labeled")) +
  base_theme
print(p_dog_ma)

## --- Plot D-3: DoG/gene body ratio per condition ---
dog_gene_com <- intersect(rownames(dog_scaled), rownames(gene_scaled))
dog_gene_com <- dog_gene_com[
  rowMeans(dog_scaled[dog_gene_com, WT_SAMPLES]) > 0 |
    rowMeans(dog_scaled[dog_gene_com, POLR2A_SAMPLES]) > 0]

dog_ratio_df <- data.frame(
  gene_id     = dog_gene_com,
  dog_wt      = rowMeans(dog_scaled[dog_gene_com, WT_SAMPLES]),
  gene_wt2    = rowMeans(gene_scaled[dog_gene_com, WT_SAMPLES]),
  dog_p       = rowMeans(dog_scaled[dog_gene_com, POLR2A_SAMPLES]),
  gene_p2     = rowMeans(gene_scaled[dog_gene_com, POLR2A_SAMPLES]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(gene_wt2 > 0, gene_p2 > 0) %>%
  dplyr::mutate(ratio_wt = dog_wt / (gene_wt2 + 1e-6),
                ratio_p  = dog_p  / (gene_p2  + 1e-6)) %>%
  tidyr::pivot_longer(cols = c(ratio_wt, ratio_p),
                      names_to = "condition", values_to = "dog_ratio") %>%
  dplyr::mutate(condition = dplyr::recode(condition,
                                          ratio_wt = "WT", ratio_p = "POLR2A"),
                condition = factor(condition, levels = c("WT","POLR2A"))) %>%
  dplyr::filter(dog_ratio > 0)

p_dog_ratio <- ggplot(dog_ratio_df,
                      aes(x = condition, y = log2(dog_ratio + 1e-6),
                          fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2,
               linewidth = 0.5, width = 0.5) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "log2(DoG / gene body RPKM)",
       title    = "DoG signal as fraction of gene body signal",
       subtitle = paste0(DOG_WINDOW/1000, "kb window")) +
  base_theme + theme(legend.position = "none")
print(p_dog_ratio)

## --- Plot D-4: DoG LFC vs gene body LFC ---
dog_gene_lfc <- dog_results %>%
  dplyr::select(gene_id, lfc_dog = lfc, de_status_dog = de_status) %>%
  dplyr::inner_join(gene_results %>% dplyr::select(gene_id, lfc_gene = lfc),
                    by = "gene_id")

r_dg <- cor(dog_gene_lfc$lfc_dog, dog_gene_lfc$lfc_gene, use = "complete.obs")
cat("Pearson r (DoG vs gene body LFC):", round(r_dg, 3), "\n")

p_dog_gene <- ggplot(dog_gene_lfc,
                     aes(x = lfc_gene, y = lfc_dog, color = de_status_dog)) +
  geom_point(size = 0.7, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = de_colors) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(r_dg, 3)),
           hjust = 1.1, vjust = -0.5, size = 4,
           color = "grey30", fontface = "italic") +
  labs(x = "Gene body log2 FC (POLR2A / WT)",
       y = "DoG log2 FC (POLR2A / WT)",
       title    = "DoG vs gene body LFC",
       subtitle = paste0("r = ", round(r_dg, 3))) +
  base_theme
print(p_dog_gene)

## --- Plot D-5: Gene expression of DoG+ vs DoG- ---
dog_char_df <- data.frame(
  gene_id   = rownames(gene_scaled),
  gene_rpkm = rowMeans(gene_scaled[, WT_SAMPLES]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(gene_rpkm > 0) %>%
  dplyr::mutate(
    dog_pos   = gene_id %in% (dog_results %>%
                                dplyr::filter(mean_wt >= DOG_THRESHOLD) %>%
                                dplyr::pull(gene_id)),
    dog_group = factor(ifelse(dog_pos, "DoG positive", "DoG negative"),
                       levels = c("DoG positive","DoG negative"))
  )

cat("\nDoG positive:", sum(dog_char_df$dog_pos), "\n")
cat("DoG negative:", sum(!dog_char_df$dog_pos), "\n")

p_dog_char <- ggplot(dog_char_df,
                     aes(x = dog_group, y = log2(gene_rpkm+1), fill = dog_group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2,
               linewidth = 0.5, width = 0.5) +
  scale_fill_manual(values = c("DoG positive" = "#7B2D8B", "DoG negative" = COL_NS)) +
  labs(x = NULL, y = "log2(gene body RPKM + 1) — WT",
       title    = "Gene body expression: DoG positive vs negative",
       subtitle = paste0("DoG positive = mean WT DoG RPKM >= ", DOG_THRESHOLD)) +
  base_theme + theme(legend.position = "none")
print(p_dog_char)

##--------------------------------------------------------------## 11. DoG extended plots
## Same plots as intronic analysis but for DoG signal

## Build DoG LFC data frame joined with genomic E/I and gene length
dog_lfc_base <- dog_results %>%
  dplyr::select(gene_id, dog_lfc = lfc) %>%
  dplyr::filter(!is.na(dog_lfc), !is.infinite(dog_lfc),
                dog_lfc >= -3, dog_lfc <= 3)

## --- Plot D-6a: DoG LFC by genomic E/I group ---
dog_lfc_ei_df <- dog_lfc_base %>%
  dplyr::left_join(
    genomic_ei_df %>% dplyr::select(gene_id, genomic_ei_group),
    by = "gene_id"
  ) %>%
  dplyr::filter(!is.na(genomic_ei_group)) %>%
  dplyr::mutate(
    genomic_ei_group = factor(genomic_ei_group,
                              levels = c("Low E/I","Medium E/I","High E/I"))
  )

dog_ei_medians <- dog_lfc_ei_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(median_lfc = median(dog_lfc, na.rm = TRUE),
                   .groups    = "drop")

dog_ei_wilcox <- dog_lfc_ei_df %>%
  dplyr::group_by(genomic_ei_group) %>%
  dplyr::summarise(
    p_value = wilcox.test(dog_lfc, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = 2.8
  )

cat("\nDoG LFC by genomic E/I — medians:\n")
print(dplyr::left_join(dog_ei_medians, dog_ei_wilcox,
                       by = "genomic_ei_group") %>%
        dplyr::select(genomic_ei_group, median_lfc, p_adj, label))

p_dog_ei_dens <- ggplot(dog_lfc_ei_df,
                        aes(x = dog_lfc, color = genomic_ei_group)) +
  geom_density(linewidth = 0.8, adjust = 1.2) +
  geom_vline(data = dog_ei_medians,
             aes(xintercept = median_lfc, color = genomic_ei_group),
             linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  scale_color_manual(values = EI_LINE_COLORS) +
  xlim(-3, 3) +
  labs(x = "log2 FC of DoG RPKM (POLR2A / WT)",
       y = "Density",
       title = "Change in DoG RPKM after POLR2A degradation — by genomic E/I group") +
  base_theme +
  theme(legend.position = "top",
        plot.title      = element_text(face = "bold", size = 13))

p_dog_ei_box <- ggplot(dog_lfc_ei_df,
                       aes(x = genomic_ei_group, y = dog_lfc,
                           color = genomic_ei_group)) +
  geom_boxplot(fill = NA, outlier.shape = NA,
               linewidth = 0.6, width = 0.5) +
  geom_hline(yintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  geom_text(data = dog_ei_wilcox,
            aes(x = genomic_ei_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, color = "grey20") +
  scale_color_manual(values = EI_LINE_COLORS) +
  scale_x_discrete(labels = c("Low E/I"    = "Low E/I\n(intron-rich)",
                              "Medium E/I" = "Medium E/I",
                              "High E/I"   = "High E/I\n(exon-rich)")) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "Genomic E/I group",
       y = "log2 FC of DoG RPKM\n(POLR2A / WT)") +
  base_theme + theme(legend.position = "none")

print(p_dog_ei_dens / p_dog_ei_box +
        patchwork::plot_layout(heights = c(2, 1)))

## --- Plot D-6b: DoG LFC by gene length quartile ---
dog_lfc_gl_df <- dog_lfc_base %>%
  dplyr::left_join(gene_length_df, by = "gene_id") %>%
  dplyr::filter(!is.na(gene_length)) %>%
  dplyr::mutate(
    length_group = dplyr::case_when(
      gene_length < quantile(gene_length, 0.25, na.rm = TRUE) ~ "Q1 (Short)",
      gene_length < quantile(gene_length, 0.50, na.rm = TRUE) ~ "Q2",
      gene_length < quantile(gene_length, 0.75, na.rm = TRUE) ~ "Q3",
      TRUE                                                      ~ "Q4 (Long)"
    ),
    length_group = factor(length_group,
                          levels = c("Q1 (Short)","Q2","Q3","Q4 (Long)"))
  )

dog_gl_medians <- dog_lfc_gl_df %>%
  dplyr::group_by(length_group) %>%
  dplyr::summarise(median_lfc = median(dog_lfc, na.rm = TRUE),
                   n          = dplyr::n(), .groups = "drop")

dog_gl_wilcox <- dog_lfc_gl_df %>%
  dplyr::group_by(length_group) %>%
  dplyr::summarise(
    p_value = wilcox.test(dog_lfc, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos = 2.8
  )

cat("\nDoG LFC by gene length quartile — medians:\n")
print(dplyr::left_join(dog_gl_medians, dog_gl_wilcox,
                       by = "length_group") %>%
        dplyr::select(length_group, n, median_lfc, p_adj, label))

p_dog_gl_dens <- ggplot(dog_lfc_gl_df,
                        aes(x = dog_lfc, color = length_group)) +
  geom_density(linewidth = 0.8, adjust = 1.2) +
  geom_vline(data = dog_gl_medians,
             aes(xintercept = median_lfc, color = length_group),
             linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  scale_color_manual(values = GL_COLORS) +
  xlim(-3, 3) +
  labs(x = "log2 FC of DoG RPKM (POLR2A / WT)",
       y = "Density",
       title = "Change in DoG RPKM after POLR2A degradation — by gene length") +
  base_theme +
  theme(legend.position = "top",
        plot.title      = element_text(face = "bold", size = 13))

p_dog_gl_box <- ggplot(dog_lfc_gl_df,
                       aes(x = length_group, y = dog_lfc,
                           color = length_group)) +
  geom_boxplot(fill = NA, outlier.shape = NA,
               linewidth = 0.6, width = 0.5) +
  geom_hline(yintercept = 0, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  geom_text(data = dog_gl_wilcox,
            aes(x = length_group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, color = "grey20") +
  scale_color_manual(values = GL_COLORS) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "Gene length quartile",
       y = "log2 FC of DoG RPKM\n(POLR2A / WT)") +
  base_theme + theme(legend.position = "none")

print(p_dog_gl_dens / p_dog_gl_box +
        patchwork::plot_layout(heights = c(2, 1)))

## --- Plot D-6c: DoG RPKM vs log10 gene length ---
dog_length_df <- data.frame(
  gene_id = rownames(dog_scaled),
  dog_wt  = rowMeans(dog_scaled[, WT_SAMPLES]),
  dog_p   = rowMeans(dog_scaled[, POLR2A_SAMPLES]),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(gene_length_df, by = "gene_id") %>%
  dplyr::filter(!is.na(gene_length), gene_length > 0,
                dog_wt > 0 | dog_p > 0) %>%
  dplyr::mutate(log10_length = log10(gene_length)) %>%
  dplyr::arrange(gene_length)

n_bins_dog <- floor(nrow(dog_length_df) / 50)
cat("Plot D-6c bins:", n_bins_dog, "\n")

dog_length_bins <- dog_length_df %>%
  dplyr::mutate(bin = cut(log10_length, breaks = n_bins_dog,
                          labels = FALSE)) %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(
    mean_log10_length = mean(log10_length, na.rm = TRUE),
    mean_dog_wt       = mean(log2(dog_wt + 1), na.rm = TRUE),
    mean_dog_polr2a   = mean(log2(dog_p  + 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(mean_dog_wt, mean_dog_polr2a),
                      names_to = "condition", values_to = "mean_log2_rpkm") %>%
  dplyr::mutate(
    condition = dplyr::recode(condition,
                              mean_dog_wt     = "WT",
                              mean_dog_polr2a = "POLR2A"),
    condition = factor(condition, levels = c("WT","POLR2A"))
  )

p_dog_length <- ggplot(dog_length_bins,
                       aes(x = mean_log10_length,
                           y = mean_log2_rpkm,
                           color = condition)) +
  geom_line(linewidth = 0.4, alpha = 0.3) +
  geom_smooth(method = "loess", span = 0.2,
              se = FALSE, linewidth = 1.0) +
  scale_color_manual(values = CONDITION_COLORS) +
  scale_x_continuous(
    breaks = log10(c(1e3, 1e4, 1e5, 1e6)),
    labels = c("1kb","10kb","100kb","1Mb")
  ) +
  labs(
    x        = "Gene length — log10 scale",
    y        = "log2(DoG RPKM + 1)",
    title    = "DoG RPKM vs gene length",
    subtitle = paste0("One bin per 50 genes  |  n=",
                      nrow(dog_length_df),
                      " genes  |  loess smoothed (span=0.2)")
  ) +
  base_theme

print(p_dog_length)

## --- Plot D-6d: DoG fraction of gene body RPKM by RPKM quartile ---
## Same as intronic fraction plot but for DoG signal

dog_frac_common <- intersect(rownames(dog_scaled), rownames(gene_scaled))
dog_frac_df <- data.frame(
  gene_id     = dog_frac_common,
  dog_wt      = rowMeans(dog_scaled[dog_frac_common,  WT_SAMPLES]),
  dog_p       = rowMeans(dog_scaled[dog_frac_common,  POLR2A_SAMPLES]),
  gene_wt2    = rowMeans(gene_scaled[dog_frac_common, WT_SAMPLES]),
  gene_p2     = rowMeans(gene_scaled[dog_frac_common, POLR2A_SAMPLES]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(gene_wt2 >= RPKM_THRESHOLD, gene_p2 >= RPKM_THRESHOLD) %>%
  dplyr::mutate(
    dog_frac_wt     = dog_wt / (gene_wt2 + 1e-6),
    dog_frac_polr2a = dog_p  / (gene_p2  + 1e-6)
  ) %>%
  dplyr::filter(dog_frac_wt     >= 0, dog_frac_wt     <= 1,
                dog_frac_polr2a >= 0, dog_frac_polr2a <= 1)

dog_frac_rpkm <- dog_frac_df %>%
  dplyr::mutate(
    rpkm_group = dplyr::case_when(
      gene_wt2 < quantile(gene_wt2, 0.25, na.rm = TRUE) ~ "Q1 (Low)",
      gene_wt2 < quantile(gene_wt2, 0.50, na.rm = TRUE) ~ "Q2",
      gene_wt2 < quantile(gene_wt2, 0.75, na.rm = TRUE) ~ "Q3",
      TRUE                                                ~ "Q4 (High)"
    ),
    rpkm_group = factor(rpkm_group,
                        levels = c("Q1 (Low)","Q2","Q3","Q4 (High)"))
  )

dog_frac_stats <- dog_frac_rpkm %>%
  dplyr::group_by(rpkm_group) %>%
  dplyr::summarise(
    p_value   = wilcox.test(dog_frac_wt, dog_frac_polr2a,
                            paired = TRUE)$p.value,
    median_wt = median(dog_frac_wt,     na.rm = TRUE),
    median_p  = median(dog_frac_polr2a, na.rm = TRUE),
    n         = dplyr::n(), .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    )
  )

dog_frac_long <- dog_frac_rpkm %>%
  tidyr::pivot_longer(cols      = c(dog_frac_wt, dog_frac_polr2a),
                      names_to  = "condition",
                      values_to = "dog_frac") %>%
  dplyr::mutate(
    condition = dplyr::recode(condition,
                              dog_frac_wt     = "WT",
                              dog_frac_polr2a = "POLR2A"),
    condition = factor(condition, levels = c("WT","POLR2A"))
  )

dog_frac_stats <- dog_frac_stats %>%
  dplyr::mutate(
    y_pos = max(dog_frac_long$dog_frac, na.rm = TRUE) * 1.05
  )

cat("\nDoG fraction by RPKM quartile (paired Wilcoxon):\n")
print(dog_frac_stats %>% dplyr::select(rpkm_group, n, median_wt,
                                       median_p, p_adj, label))

p_dog_frac_rpkm <- ggplot(dog_frac_long,
                          aes(x = condition, y = dog_frac,
                              fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2,
               linewidth = 0.5, width = 0.55) +
  geom_text(data = dog_frac_stats,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "grey20") +
  geom_segment(data = dog_frac_stats,
               aes(x = 1, xend = 2,
                   y = y_pos * 0.98, yend = y_pos * 0.98),
               inherit.aes = FALSE, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = CONDITION_COLORS) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, max(dog_frac_long$dog_frac,
                                       na.rm = TRUE) * 1.12)) +
  facet_wrap(~ rpkm_group, nrow = 1) +
  labs(x = NULL,
       y = "DoG fraction (DoG / gene body RPKM)",
       title    = "DoG signal fraction by gene body RPKM quartile",
       subtitle = "Paired Wilcoxon WT vs POLR2A  |  * p<0.05  ** p<0.01  *** p<0.001") +
  base_theme + theme(legend.position = "none")

print(p_dog_frac_rpkm)

## --- Plot I-6a-intron-gl: Intronic RPKM LFC by gene length quartile ---

## DoG RPKM by gene body RPKM quartile — raw RPKM not fraction
dog_rpkm_long <- dog_frac_rpkm %>%
  dplyr::select(gene_id, rpkm_group, dog_wt, dog_p) %>%
  tidyr::pivot_longer(cols      = c(dog_wt, dog_p),
                      names_to  = "condition",
                      values_to = "dog_rpkm") %>%
  dplyr::mutate(
    condition = dplyr::recode(condition,
                              dog_wt = "WT",
                              dog_p  = "POLR2A"),
    condition = factor(condition, levels = c("WT","POLR2A")),
    log2_rpkm = log2(dog_rpkm + 1)
  )

dog_rpkm_stats <- dog_rpkm_long %>%
  dplyr::group_by(rpkm_group) %>%
  dplyr::summarise(
    p_value = wilcox.test(log2_rpkm[condition == "WT"],
                          log2_rpkm[condition == "POLR2A"],
                          paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj  = p.adjust(p_value, method = "BH"),
    label  = dplyr::case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    ),
    y_pos  = max(dog_rpkm_long$log2_rpkm, na.rm = TRUE) * 1.05
  )

p_dog_rpkm_quartile <- ggplot(dog_rpkm_long,
                              aes(x = condition, y = log2_rpkm,
                                  fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2,
               linewidth = 0.5, width = 0.55) +
  geom_text(data = dog_rpkm_stats,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "grey20") +
  geom_segment(data = dog_rpkm_stats,
               aes(x = 1, xend = 2,
                   y = y_pos * 0.98, yend = y_pos * 0.98),
               inherit.aes = FALSE, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = CONDITION_COLORS) +
  scale_y_continuous(limits = c(0, max(dog_rpkm_long$log2_rpkm,
                                       na.rm = TRUE) * 1.12)) +
  facet_wrap(~ rpkm_group, nrow = 1) +
  labs(x = NULL,
       y = "log2(DoG RPKM + 1)",
       title    = "DoG RPKM by gene body RPKM quartile",
       subtitle = "Paired Wilcoxon WT vs POLR2A  |  * p<0.05  ** p<0.01  *** p<0.001") +
  base_theme + theme(legend.position = "none")

print(p_dog_rpkm_quartile)

message("\nAll intron and DoG analysis plots complete.")