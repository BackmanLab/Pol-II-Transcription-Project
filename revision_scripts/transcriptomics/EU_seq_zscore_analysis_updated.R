##################################################################
# Author: Lucas Carter
# Email: lucascarter2025@u.northwestern.edu
# PI: Vadim Backman
# Description:
# EU-seq nascent transcription analysis using FPKM values from
# NoEU-subtracted bedgraphs.
#
# Workflow:
#   QC:      Pearson correlation heatmap + PCA (all 7 samples)
#   Filter:  Inspect FPKM distribution to set expression threshold
#   Z-score: Z-test on log2(FPKM+1) across POLR2A vs WT replicates
#            with BH FDR correction
#   Plots:   MA, volcano, scatter vs total RNA-seq,
#            LFC boxplot by expression group,
#            DEG overlap bar chart,
#            cumulative LFC distribution,
#            RRHO-style rank correlation heatmap
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(BiocGenerics)
})

##--------------------------------------------------------------## 1. Parameters

## hg38 EU-seq rerun (euseq_rpkm_hg38.sh output)
BASE         <- "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq"
CNT_DIR      <- file.path(BASE, "euseq/euseq_pipeline/counts")
QC_MATRIX    <- file.path(CNT_DIR, "euseq_rpkm_qc.txt")
ANAL_MATRIX  <- file.path(CNT_DIR, "euseq_rpkm_analysis.txt")
GTF_PATH     <- file.path(BASE, "reference/Homo_sapiens.GRCh38.115.chr.gtf")
SF_FILE      <- file.path(BASE, "scale_factors_mito.txt")   ## chrM scale factors (compute_scale_factors_mito_hg38.sh)
TOTAL_DESEQ  <- file.path(BASE, "reference/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv")
## HCT116 total RNA-seq gene quantification: ENCODE ENCFF435PHM (experiment ENCSR698RPL,
## RSEM gene quantifications, GRCh38, GENCODE V29)
EXPR_FILE    <- file.path(BASE, "reference/ENCFF435PHM.tsv")

## Analysis parameters — set after inspecting distribution plots
## These are initial defaults; adjust after running Section 4
FPKM_THRESHOLD <- NULL   ## set to NULL to trigger distribution plot first
LFC_THRESHOLD  <- 1
FDR_THRESHOLD  <- 0.05
N_LABELS       <- 25

## Colors — DE status (Okabe-Ito red/blue)
COL_UP     <- "#D55E00"   ## vermillion — upregulated
COL_DOWN   <- "#0072B2"   ## blue — downregulated
COL_NS     <- "grey70"
COL_FILTER <- "grey85"

## Sample metadata
POLR2A_SAMPLES <- c("POLR2A_1", "POLR2A_2", "POLR2A_3")
WT_SAMPLES     <- c("WT_1", "WT_2", "WT_3")
NOEU_SAMPLE    <- "NoEU"

CONDITION_MAP <- c(
  POLR2A_1 = "POLR2A", POLR2A_2 = "POLR2A", POLR2A_3 = "POLR2A",
  WT_1     = "WT",     WT_2     = "WT",     WT_3     = "WT",
  NoEU     = "NoEU"
)

## Colors — conditions (color-blind friendly purple/green, Paul Tol palette)
CONDITION_COLORS <- c(
  POLR2A = "#7B2D8B",   ## purple — POLR2A degraded
  WT     = "#1B7837",   ## green  — wildtype
  NoEU   = "#999999"    ## grey   — NoEU control
)

## chrM (mitochondrial) WT-anchored scale factors
## eu_sf = mean chrM RPM (WT) / mean chrM RPM (POLR2A-degraded), WT = reference (scale = 1)
## Read from scale_factors_mito.txt (compute_scale_factors_mito_hg38.sh) and applied
## to the NoEU-subtracted EU-seq RPKM, consistent with the bigWig normalization —
## EU-seq cannot normalize itself due to variable EU incorporation and pulldown
## efficiency; chrM transcripts act as an endogenous spike-in
sf_tbl <- read.delim(SF_FILE, stringsAsFactors = FALSE)
print(sf_tbl)
stopifnot(all(c("WT", "POLR2A_8hrs_aux") %in% sf_tbl$condition),
          sf_tbl$eu_sf[sf_tbl$condition == "WT"] == 1)
SF_POLR2A <- sf_tbl$eu_sf[sf_tbl$condition == "POLR2A_8hrs_aux"]
TOTAL_RNASEQ_SF <- c(
  POLR2A_1 = SF_POLR2A, POLR2A_2 = SF_POLR2A, POLR2A_3 = SF_POLR2A,
  WT_1     = 1,         WT_2     = 1,         WT_3     = 1
)
message(sprintf("EU-seq chrM scale factor: POLR2A = %.6f (WT = 1)", SF_POLR2A))

##--------------------------------------------------------------## 2. Load matrices

message("Loading FPKM matrices...")
qc_mat   <- read_tsv(QC_MATRIX,   show_col_types = FALSE)
anal_mat <- read_tsv(ANAL_MATRIX, show_col_types = FALSE)

## Strip version suffixes from gene IDs
qc_mat$gene_id   <- sub("\\..*", "", qc_mat$gene_id)
anal_mat$gene_id <- sub("\\..*", "", anal_mat$gene_id)

## Convert to numeric matrices
to_matrix <- function(df) {
  mat <- as.matrix(df[, -1])
  rownames(mat) <- df$gene_id
  mode(mat) <- "numeric"
  mat
}

qc_fpkm   <- to_matrix(qc_mat)
anal_fpkm <- to_matrix(anal_mat)

cat("QC matrix:       ", nrow(qc_fpkm),   "genes x", ncol(qc_fpkm),   "samples\n")
cat("Analysis matrix: ", nrow(anal_fpkm), "genes x", ncol(anal_fpkm), "samples\n")

##--------------------------------------------------------------## 2b. QC plots — raw RPKM (all 7 samples)
## All QC performed on raw RPKM before scaling to assess
## technical quality of the experiment independently of
## biological normalization.

message("Running distributional QC on raw RPKM...")

## log2(RPKM + 1) for QC
log2_qc <- log2(qc_fpkm + 1)

## Filter genes with any signal for QC
log2_qc_filt <- log2_qc[rowSums(log2_qc) > 0, ]
cat("Genes with any signal for QC:", nrow(log2_qc_filt), "\n")

## Apply CPM-style filter — keep genes with log2(RPKM+1) >= 1
## in at least 3 samples
keep_qc <- rowSums(log2_qc >= log2(1 + 1)) >= 3
log2_qc_expressed <- log2_qc[keep_qc, ]
cat("Genes passing expression filter (log2 RPKM+1 >= 1 in >= 3 samples):",
    nrow(log2_qc_expressed), "\n")

##-- QC Plot 1: Boxplot of log2(RPKM+1) per sample
## Quick check that RPM normalization produced comparable
## signal levels across samples

box_df <- as.data.frame(log2_qc_filt) %>%
  tibble::rownames_to_column("gene_id") %>%
  tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "log2rpkm") %>%
  dplyr::mutate(
    condition = CONDITION_MAP[sample],
    condition = factor(condition, levels = c("POLR2A", "WT", "NoEU")),
    sample    = factor(sample, levels = colnames(qc_fpkm))
  )

p_box_qc <- ggplot(box_df, aes(x = sample, y = log2rpkm,
                               fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2,
               linewidth = 0.5) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(
    x        = NULL,
    y        = "log2(RPKM + 1)",
    title    = "EU-seq RPKM distribution per sample",
    subtitle = paste0("Raw RPKM  |  ", nrow(log2_qc_filt),
                      " genes with signal  |  7 samples")
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x   = element_text(angle = 30, hjust = 1),
    legend.title  = element_blank(),
    plot.title    = element_text(face = "bold")
  )

print(p_box_qc)

##-- QC Plot 2: Histogram of log2(RPKM+1) per sample
## Overlapping density curves — check for bimodal distribution
## and consistent shape across replicates

dens_df <- box_df %>%
  dplyr::filter(log2rpkm > 0)

p_dens_qc <- ggplot(dens_df, aes(x = log2rpkm,
                                 color = sample,
                                 linetype = condition)) +
  geom_density(linewidth = 0.7, adjust = 1.5) +
  scale_color_manual(
    values = setNames(
      c(rep("#7B2D8B", 3), rep("#1B7837", 3), "#999999"),
      colnames(qc_fpkm)
    )
  ) +
  scale_linetype_manual(values = c(POLR2A = "solid",
                                   WT     = "solid",
                                   NoEU   = "dashed")) +
  labs(
    x        = "log2(RPKM + 1)",
    y        = "Density",
    title    = "EU-seq RPKM density per sample",
    subtitle = "Raw RPKM  |  genes with RPKM > 0  |  7 samples"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank(),
        plot.title   = element_text(face = "bold"))

print(p_dens_qc)

##-- QC Plot 3: Pearson correlation heatmap

cor_mat <- cor(log2_qc_expressed, method = "pearson")

cat("\nPearson correlation matrix:\n")
print(round(cor_mat, 4))

## Within-condition mean r
for (cond in unique(CONDITION_MAP)) {
  samps <- names(CONDITION_MAP)[CONDITION_MAP == cond]
  samps <- samps[samps %in% colnames(cor_mat)]
  if (length(samps) < 2) next
  vals <- cor_mat[samps, samps]
  cat(sprintf("  Mean within-%s Pearson r: %.4f\n",
              cond, mean(vals[upper.tri(vals)])))
}

annotation_col <- data.frame(
  Condition = CONDITION_MAP[colnames(cor_mat)],
  row.names = colnames(cor_mat)
)
annotation_colors <- list(Condition = CONDITION_COLORS)

cor_min <- floor(min(cor_mat) * 100) / 100
breaks  <- seq(cor_min, 1, length.out = 101)
hm_colors <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(9, "RdYlBu"))
)(100)

pheatmap(
  cor_mat,
  color             = hm_colors,
  breaks            = breaks,
  annotation_col    = annotation_col,
  annotation_row    = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  display_numbers   = TRUE,
  number_format     = "%.3f",
  number_color      = "grey20",
  fontsize_number   = 8,
  fontsize          = 11,
  cellwidth         = 50,
  cellheight        = 50,
  border_color      = NA,
  main              = paste0(
    "EU-seq — Pearson correlation (log2 RPKM+1, raw)\n",
    nrow(log2_qc_expressed), " expressed genes  |  7 samples"
  )
)

##-- QC Plot 4: PCA

pca_res <- prcomp(t(log2_qc_expressed), center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x[, 1:3]) %>%
  dplyr::mutate(
    sample    = rownames(.),
    condition = factor(CONDITION_MAP[sample],
                       levels = c("POLR2A", "WT", "NoEU"))
  )

pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_theme <- theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey40", size = 11),
        legend.title  = element_blank())

p_pca12 <- ggplot(pca_df,
                  aes(x = PC1, y = PC2,
                      color = condition, label = sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(
    x        = paste0("PC1 (", pct_var[1], "% variance)"),
    y        = paste0("PC2 (", pct_var[2], "% variance)"),
    title    = "EU-seq PCA — PC1 vs PC2",
    subtitle = "Raw log2(RPKM+1)  |  7 samples"
  ) + pca_theme

p_pca13 <- ggplot(pca_df,
                  aes(x = PC1, y = PC3,
                      color = condition, label = sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(
    x        = paste0("PC1 (", pct_var[1], "% variance)"),
    y        = paste0("PC3 (", pct_var[3], "% variance)"),
    title    = "EU-seq PCA — PC1 vs PC3",
    subtitle = "Raw log2(RPKM+1)  |  7 samples"
  ) + pca_theme

print(p_pca12 | p_pca13)

##-- QC Plot 5: RPKM distribution with candidate thresholds
## Plot on WT samples only — threshold should be set based
## on WT signal since that is the reference condition

wt_mean_raw <- rowMeans(qc_fpkm[, WT_SAMPLES])

thresh_df <- data.frame(
  gene_id   = rownames(qc_fpkm),
  wt_mean   = wt_mean_raw,
  log2_rpkm = log2(wt_mean_raw + 1)
) %>%
  dplyr::filter(wt_mean > 0)

p_thresh <- ggplot(thresh_df, aes(x = log2_rpkm)) +
  geom_histogram(bins = 100, fill = "#1B7837",
                 color = NA, alpha = 0.7) +
  geom_vline(xintercept = log2(c(0.1, 0.5, 1.0) + 1),
             linetype = "dashed",
             color    = c("#CC79A7", "#E69F00", "#D55E00"),
             linewidth = 0.7) +
  annotate("text",
           x     = log2(c(0.1, 0.5, 1.0) + 1) + 0.05,
           y     = Inf, vjust = 1.5, hjust = 0,
           label = c("RPKM=0.1", "RPKM=0.5", "RPKM=1.0"),
           color = c("#CC79A7", "#E69F00", "#D55E00"),
           size  = 3.5) +
  labs(
    x        = "log2(mean WT RPKM + 1)",
    y        = "Number of genes",
    title    = "EU-seq WT RPKM distribution",
    subtitle = paste0("Raw RPKM  |  genes with WT RPKM > 0: ",
                      nrow(thresh_df),
                      "  |  Use to set RPKM_THRESHOLD")
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_thresh)

for (thresh in c(0.1, 0.5, 1.0)) {
  cat(sprintf("  Genes with mean WT RPKM >= %.1f: %d\n",
              thresh, sum(wt_mean_raw >= thresh)))
}

message("\n*** Inspect QC plots above, then set RPKM_THRESHOLD",
        " in Section 1 before continuing. ***\n")

RPKM_THRESHOLD=1.0

if (is.null(RPKM_THRESHOLD)) {
  stop("RPKM_THRESHOLD not set — inspect distribution plot first.",
       call. = FALSE)
}

##--------------------------------------------------------------## 2c. Apply chrM scale factors to EU-seq FPKM
## Scale factors anchored to WT chrM RPM; applied after the per-nucleotide
## NoEU subtraction (floored at 0).
## Applied per-sample by multiplying each FPKM column by its
## corresponding scale factor. WT samples scale = 1 (no change).
## QC matrix is left unscaled — QC checks raw technical consistency.

message("Applying chrM scale factors to analysis matrix...")

sf_vec <- TOTAL_RNASEQ_SF[colnames(anal_fpkm)]

if (any(is.na(sf_vec))) {
  stop("Scale factor missing for samples: ",
       paste(colnames(anal_fpkm)[is.na(sf_vec)], collapse = ", "))
}

anal_fpkm_scaled <- sweep(anal_fpkm, 2, sf_vec, FUN = "*")

cat("Scale factors applied:\n")
for (s in names(sf_vec)) cat(sprintf("  %-12s %.6f\n", s, sf_vec[s]))
cat("FPKM range before scaling: [",
    round(min(anal_fpkm), 3), ",", round(max(anal_fpkm), 3), "]\n")
cat("FPKM range after scaling:  [",
    round(min(anal_fpkm_scaled), 3), ",", round(max(anal_fpkm_scaled), 3), "]\n")


##--------------------------------------------------------------## 5. Fold-change analysis using Z-score as effect size
## log2(RPKM + 1) transformation
## Z-score = (mean_POLR2A - mean_WT) / pooled_SD
##         = effect size metric (like Cohen's d)
## No p-value thresholding — genes classified by LFC threshold only.
## This is appropriate given high inter-replicate variability
## in EU-seq with n=3, where formal statistical tests are underpowered.

message("Running fold-change analysis with RPKM threshold: ", RPKM_THRESHOLD)

log2_anal <- log2(anal_fpkm_scaled + 1)

## Filter to expressed genes (mean WT RPKM >= threshold)
wt_mean   <- rowMeans(anal_fpkm_scaled[, WT_SAMPLES])
keep      <- wt_mean >= RPKM_THRESHOLD
log2_filt <- log2_anal[keep, ]
cat("Genes passing RPKM filter:", nrow(log2_filt), "\n")

## Per-gene statistics
polr2a_vals <- log2_filt[, POLR2A_SAMPLES]
wt_vals     <- log2_filt[, WT_SAMPLES]

mean_polr2a <- rowMeans(polr2a_vals)
mean_wt     <- rowMeans(wt_vals)
lfc         <- mean_polr2a - mean_wt

## Pooled SD across all 6 replicates — used as effect size denominator
all_vals  <- cbind(polr2a_vals, wt_vals)
pooled_sd <- apply(all_vals, 1, sd)

## Z-score as effect size (signal-to-noise ratio per gene)
z_score <- lfc / (pooled_sd + 1e-10)

## Mean RPKM across all 6 replicates (for MA plot x-axis)
mean_rpkm_all <- rowMeans(anal_fpkm_scaled[keep, ])

## Results data frame — classify by LFC threshold only, no p-value
eu_results <- data.frame(
  gene_id     = rownames(log2_filt),
  mean_rpkm   = mean_rpkm_all,
  mean_polr2a = mean_polr2a,
  mean_wt     = mean_wt,
  lfc         = lfc,
  pooled_sd   = pooled_sd,
  z_score     = z_score,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    de_status = dplyr::case_when(
      lfc >  LFC_THRESHOLD ~ "Up",
      lfc < -LFC_THRESHOLD ~ "Down",
      TRUE                 ~ "Not DE"
    ),
    de_status = factor(de_status, levels = c("Up", "Not DE", "Down"))
  )

cat("\nFold-change analysis complete.\n")
cat("Gene counts by LFC threshold (|LFC| >", LFC_THRESHOLD, "):\n")
print(table(eu_results$de_status))

## Plot Z-score distribution
p_zdist <- ggplot(eu_results, aes(x = z_score)) +
  geom_histogram(bins = 100, fill = "#0072B2", color = NA, alpha = 0.7) +
  geom_vline(xintercept = c(-2, 2),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(
    x        = "Z-score (effect size)",
    y        = "Number of genes",
    title    = "EU-seq Z-score distribution",
    subtitle = paste0(nrow(eu_results), " expressed genes  |  ",
                      "Z-score = effect size, not significance")
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

## Plot LFC distribution colored by threshold
p_lfcdist <- ggplot(eu_results, aes(x = lfc, fill = de_status)) +
  geom_histogram(bins = 100, color = NA, alpha = 0.8) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = c(Up = COL_UP, "Not DE" = COL_NS,
                               Down = COL_DOWN)) +
  labs(
    x        = "log2 fold change (POLR2A / WT)",
    y        = "Number of genes",
    title    = "EU-seq LFC distribution",
    subtitle = paste0(nrow(eu_results), " expressed genes  |  ",
                      "colored by |LFC| > ", LFC_THRESHOLD)
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title   = element_text(face = "bold"),
        legend.title = element_blank())

print(p_zdist | p_lfcdist)

##--------------------------------------------------------------## 6. Load total RNA-seq DESeq2 results + TPM

message("Loading total RNA-seq results and TPM...")

total_df <- read_csv(TOTAL_DESEQ, show_col_types = FALSE) %>%
  dplyr::select(-1) %>%
  dplyr::mutate(
    gene_id         = sub("\\..*", "", gene_id),
    de_status_total = dplyr::case_when(
      is.na(padj)                                               ~ "Filtered",
      padj >= FDR_THRESHOLD                                     ~ "Not DE",
      padj <  FDR_THRESHOLD & log2FoldChange >  LFC_THRESHOLD  ~ "Up",
      padj <  FDR_THRESHOLD & log2FoldChange < -LFC_THRESHOLD  ~ "Down",
      TRUE                                                      ~ "Not DE"
    ),
    de_status_total = factor(de_status_total,
                             levels = c("Up", "Not DE", "Down", "Filtered"))
  ) %>%
  dplyr::select(gene_id,
                lfc_total    = log2FoldChange,
                padj_total   = padj,
                de_status_total)

## ENCODE RSEM gene quantification: one row per gene, columns gene_id, TPM, ...
expr_raw <- read_tsv(EXPR_FILE, show_col_types = FALSE)
cat("Expression file columns:", paste(colnames(expr_raw), collapse = ", "), "\n")
stopifnot(all(c("gene_id", "TPM") %in% colnames(expr_raw)))

expr <- expr_raw %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
  dplyr::select(gene_id, TPM) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")
cat("Genes with TPM:", nrow(expr), "| TPM >= 1:", sum(expr$TPM >= 1), "\n")

tpm_exp <- expr$TPM[expr$TPM >= 1]
expr <- expr %>%
  dplyr::mutate(
    expr_group = dplyr::case_when(
      TPM < 1                                           ~ "Unexpressed",
      TPM < quantile(tpm_exp, 1/3, na.rm = TRUE)       ~ "Low",
      TPM < quantile(tpm_exp, 2/3, na.rm = TRUE)       ~ "Medium",
      TRUE                                              ~ "High"
    ),
    expr_group = factor(expr_group,
                        levels = c("Unexpressed", "Low", "Medium", "High"))
  )

## Merge all results
combined_df <- eu_results %>%
  dplyr::left_join(total_df, by = "gene_id") %>%
  dplyr::left_join(expr,     by = "gene_id")

cat("Combined df dimensions:", nrow(combined_df), "genes\n")
cat("Genes with both EU and total results:",
    sum(!is.na(combined_df$lfc_total)), "\n")

##--------------------------------------------------------------## 6b. Compute E/I ratio from GTF
## Exon/intron ratio per gene — computed from GTF exon features.
## E/I = total exon length / total intron length
## High E/I = exon-dense genes (short introns or many exons)
## Low E/I  = intron-rich genes (long introns dominate gene body)

message("Computing E/I ratios from GTF...")

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
})

txdb <- GenomicFeatures::makeTxDbFromGFF(
  GTF_PATH,
  format = "gtf",
  taxonomyId = 9606
)

## Exon ranges per gene
exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
exons_by_gene <- GenomicRanges::reduce(exons_by_gene)
exon_lengths  <- sum(BiocGenerics::width(exons_by_gene))

## Gene ranges (TSS to TES)
genes_gr    <- GenomicFeatures::genes(txdb)
gene_lengths <- BiocGenerics::width(genes_gr)
names(gene_lengths) <- genes_gr$gene_id

## Intron length = gene length - exon length
common_genes  <- intersect(names(exon_lengths), names(gene_lengths))
intron_lengths <- gene_lengths[common_genes] - exon_lengths[common_genes]
intron_lengths[intron_lengths < 0] <- 0

## E/I ratio — add pseudocount to avoid division by zero
ei_ratio <- exon_lengths[common_genes] /
  (intron_lengths + 1)

## Strip version suffixes
ei_df <- data.frame(
  gene_id  = sub("\\..*", "", common_genes),
  ei_ratio = as.numeric(ei_ratio),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(ei_ratio), ei_ratio > 0) %>%
  dplyr::mutate(
    ei_group = dplyr::case_when(
      ei_ratio < quantile(ei_ratio, 1/3, na.rm = TRUE) ~ "Low E/I",
      ei_ratio < quantile(ei_ratio, 2/3, na.rm = TRUE) ~ "Medium E/I",
      TRUE                                              ~ "High E/I"
    ),
    ei_group = factor(ei_group,
                      levels = c("Low E/I", "Medium E/I", "High E/I"))
  )

cat("E/I ratios computed for", nrow(ei_df), "genes\n")
cat("E/I group distribution:\n")
print(table(ei_df$ei_group))

## Join E/I to combined_df
combined_df <- combined_df %>%
  dplyr::left_join(ei_df %>% dplyr::select(gene_id, ei_ratio, ei_group),
                   by = "gene_id")

cat("Genes with E/I annotation:",
    sum(!is.na(combined_df$ei_group)), "\n")

##--------------------------------------------------------------## 7. Plot theme

base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11),
    axis.line        = element_line(color = "grey40"),
    axis.ticks       = element_line(color = "grey40")
  )

COL_DOWN="#e88471"
COL_UP="#39b185"

de_colors <- c(Up = COL_UP, Down = COL_DOWN, "Not DE" = COL_NS)

##--------------------------------------------------------------## 8. Plot 1 — MA plot (EU-seq) ALL GENES

message("Plot 1: MA plot...")
## Extract gene_id -> gene_symbol from GTF
gtf_genes <- rtracklayer::import(GTF_PATH, feature.type = "gene") %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_symbol = gene_name) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

combined_df <- combined_df %>%
  dplyr::left_join(gtf_genes, by = "gene_id") %>%
  dplyr::mutate(gene_symbol = dplyr::coalesce(gene_symbol, gene_id))

## Count string for subtitle
de_counts   <- combined_df %>%
  dplyr::count(de_status) %>%
  dplyr::filter(de_status %in% c("Up", "Down")) %>%
  dplyr::arrange(desc(de_status))

count_label <- paste(paste0(de_counts$de_status, ": ", de_counts$n), collapse = "  |  ")
total_label <- paste0("n = ", nrow(combined_df), "  |  ", count_label)

## Top genes to label
label_ma <- combined_df %>%
  dplyr::filter(de_status %in% c("Up", "Down")) %>%
  dplyr::arrange(desc(abs(lfc))) %>%
  dplyr::slice_head(n = N_LABELS)

## MA plot
p_ma <- ggplot(
  combined_df,
  aes(x = log10(mean_rpkm + 1), y = lfc, color = de_status)
) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = label_ma, aes(label = gene_symbol),
    size = 2.8, max.overlaps = 20,
    segment.size = 0.3, segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(values = de_colors) +
  labs(
    x        = "log10(mean FPKM + 1)",
    y        = "log2 fold change (POLR2A / WT)",
    title    = "EU-seq MA plot — POLR2A degradation (nascent RNA)",
    subtitle = paste0(total_label,
                      "\n|LFC| > ", LFC_THRESHOLD,
                      "  |  Top ", N_LABELS, " genes labeled by |LFC|")
  ) +
  base_theme

print(p_ma)

##--------------------------------------------------------------## 8b. Plot 1 — MA plot (EU-seq) DEGs in total RNA seq


## Classify concordant DE overlap between EU-seq and total RNA-seq
combined_df <- combined_df %>%
  dplyr::mutate(
    overlap_status = dplyr::case_when(
      de_status == "Up"   & de_status_total == "Up"   ~ "Up (both)",
      de_status == "Down" & de_status_total == "Down" ~ "Down (both)",
      TRUE                                             ~ "Other"
    ),
    overlap_status = factor(overlap_status,
                            levels = c("Up (both)", "Down (both)", "Other"))
  )

## Top concordant genes to label
label_ma_overlap <- combined_df %>%
  dplyr::filter(overlap_status %in% c("Up (both)", "Down (both)")) %>%
  dplyr::arrange(desc(abs(lfc))) %>%
  dplyr::slice_head(n = N_LABELS)

## Count string
overlap_counts   <- combined_df %>%
  dplyr::count(overlap_status) %>%
  dplyr::filter(overlap_status %in% c("Up (both)", "Down (both)")) %>%
  dplyr::arrange(desc(overlap_status))

overlap_label <- paste(paste0(overlap_counts$overlap_status, ": ",
                              overlap_counts$n), collapse = "  |  ")
overlap_total <- paste0("n = ", nrow(combined_df), "  |  ", overlap_label)

## Colors
overlap_colors <- c(
  "Up (both)"   = COL_UP,
  "Down (both)" = COL_DOWN,
  "Other"       = COL_NS
)

## Symmetric y-axis: same distance above and below 0
MA_YLIM <- ceiling(max(abs(combined_df$lfc), na.rm = TRUE))
cat("MA plot y-axis range: +/-", MA_YLIM, "\n")

## MA plot — concordant DEGs only
p_ma_overlap <- ggplot(
  combined_df,
  aes(x = log10(mean_rpkm + 1), y = lfc, color = overlap_status)
) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = label_ma_overlap, aes(label = gene_symbol),
    size = 2.8, max.overlaps = 20,
    segment.size = 0.3, segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(values = overlap_colors) +
  scale_y_continuous(breaks = seq(-floor(MA_YLIM / 2.5) * 2.5, floor(MA_YLIM / 2.5) * 2.5, by = 2.5)) +
  coord_cartesian(ylim = c(-MA_YLIM, MA_YLIM)) +
  labs(
    x        = "log10(mean FPKM + 1)",
    y        = "log2 fold change (POLR2A / WT)",
    title    = "EU-seq MA plot — concordant with total RNA-seq DEGs",
    subtitle = paste0(overlap_total,
                      "\n|LFC| > ", LFC_THRESHOLD,
                      "  |  concordant direction in both EU-seq and total RNA-seq")
  ) +
  base_theme

print(p_ma_overlap)

##--------------------------------------------------------------## 8b. Plot 1b — E/I distribution plots (LFC and RPKM)
## Two density plots side by side:
##   Left:  LFC distribution by E/I group
##   Right: log2(RPKM+1) distribution by E/I group
## Shows whether E/I ratio predicts nascent transcription level
## and/or fold-change response to POLR2A degradation

message("Plot 1b: E/I distribution plots...")

ei_colors <- c(
  "Low E/I"    = "#0072B2",
  "Medium E/I" = "#E69F00",
  "High E/I"   = "#CC79A7"
)

ei_plot_df <- combined_df %>%
  dplyr::filter(!is.na(ei_group), !is.na(lfc), !is.na(mean_rpkm))

## Left: LFC distribution by E/I group
p_ei_lfc <- ggplot(ei_plot_df,
                   aes(x = lfc, color = ei_group, fill = ei_group)) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.2) +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dotted", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = ei_colors) +
  scale_fill_manual(values  = ei_colors) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(
    x        = "log2 fold change (POLR2A / WT)",
    y        = "Density",
    title    = "LFC distribution by E/I group",
    subtitle = "EU-seq  |  POLR2A degradation"
  ) +
  base_theme +
  theme(legend.position = "top",
        legend.title    = element_blank())

## Right: log2(RPKM+1) distribution by E/I group
p_ei_rpkm <- ggplot(ei_plot_df,
                    aes(x     = log2(mean_rpkm + 1),
                        color = ei_group,
                        fill  = ei_group)) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.2) +
  scale_color_manual(values = ei_colors) +
  scale_fill_manual(values  = ei_colors) +
  labs(
    x        = "log2(mean RPKM + 1)",
    y        = "Density",
    title    = "RPKM distribution by E/I group",
    subtitle = "EU-seq  |  mean across all 6 replicates"
  ) +
  base_theme +
  theme(legend.position = "top",
        legend.title    = element_blank())

print(p_ei_lfc | p_ei_rpkm)

##--------------------------------------------------------------## 8c. Plot 1b — E/I distribution plots (LFC and RPKM) with boxplot on margin

## Left: LFC density (same as before)
p_ei_lfc <- ggplot(ei_plot_df,
                   aes(x = lfc, color = ei_group, fill = ei_group)) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.2) +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dotted", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = ei_colors) +
  scale_fill_manual(values  = ei_colors) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(x = NULL, y = "Density",
       title    = "LFC distribution by E/I group",
       subtitle = "EU-seq  |  POLR2A degradation") +
  base_theme +
  theme(legend.position  = "top",
        legend.title     = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

## Left: LFC boxplot strip
p_ei_lfc_box <- ggplot(ei_plot_df,
                       aes(x = lfc, y = ei_group,
                           color = ei_group, fill = ei_group)) +
  geom_boxplot(alpha = 0.15, linewidth = 0.6, outlier.size = 0.3,
               outlier.alpha = 0.2, width = 0.5) +
  scale_color_manual(values = ei_colors) +
  scale_fill_manual(values  = ei_colors) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(x = "log2 fold change (POLR2A / WT)", y = NULL) +
  base_theme +
  theme(legend.position = "none",
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        axis.line.y     = element_blank())

## Right: RPKM density
p_ei_rpkm <- ggplot(ei_plot_df,
                    aes(x = log2(mean_rpkm + 1),
                        color = ei_group, fill = ei_group)) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.2) +
  scale_color_manual(values = ei_colors) +
  scale_fill_manual(values  = ei_colors) +
  labs(x = NULL, y = "Density",
       title    = "RPKM distribution by E/I group",
       subtitle = "EU-seq  |  mean across all 6 replicates") +
  base_theme +
  theme(legend.position  = "top",
        legend.title     = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

## Right: RPKM boxplot strip
p_ei_rpkm_box <- ggplot(ei_plot_df,
                        aes(x = log2(mean_rpkm + 1), y = ei_group,
                            color = ei_group, fill = ei_group)) +
  geom_boxplot(alpha = 0.15, linewidth = 0.6, outlier.size = 0.3,
               outlier.alpha = 0.2, width = 0.5) +
  scale_color_manual(values = ei_colors) +
  scale_fill_manual(values  = ei_colors) +
  labs(x = "log2(mean RPKM + 1)", y = NULL) +
  base_theme +
  theme(legend.position = "none",
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        axis.line.y     = element_blank())

## Stack each density + boxplot, then place side by side
p_left  <- p_ei_lfc  / p_ei_lfc_box  + plot_layout(heights = c(5, 1))
p_right <- p_ei_rpkm / p_ei_rpkm_box + plot_layout(heights = c(5, 1))

print(p_left | p_right)

##--------------------------------------------------------------## 8d. — expression distribution plots (LFC and RPKM)


message("Plot 1b: E/I distribution plots...")

expr_colors <- c(
  "Unexpressed" = "#999999",
  "Low"         = "#0072B2",
  "Medium"      = "#E69F00",
  "High"        = "#CC79A7"
)

expr_plot_df <- combined_df %>%
  dplyr::filter(!is.na(expr_group), !is.na(lfc), !is.na(mean_rpkm),
                expr_group != "Unexpressed") %>%
  dplyr::mutate(expr_group = droplevels(expr_group))

## Left: LFC density
p_expr_lfc <- ggplot(expr_plot_df,
                     aes(x = lfc, color = expr_group, fill = expr_group)) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.2) +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dotted", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = expr_colors) +
  scale_fill_manual(values  = expr_colors) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(x = NULL, y = "Density",
       title    = "LFC distribution by expression group",
       subtitle = "EU-seq  |  POLR2A degradation") +
  base_theme +
  theme(legend.position  = "top",
        legend.title     = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

## Left: LFC boxplot strip
p_expr_lfc_box <- ggplot(expr_plot_df,
                         aes(x = lfc, y = expr_group,
                             color = expr_group, fill = expr_group)) +
  geom_boxplot(alpha = 0.15, linewidth = 0.6, outlier.size = 0.3,
               outlier.alpha = 0.2, width = 0.5) +
  scale_color_manual(values = expr_colors) +
  scale_fill_manual(values  = expr_colors) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(x = "log2 fold change (POLR2A / WT)", y = NULL) +
  base_theme +
  theme(legend.position = "none",
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        axis.line.y     = element_blank())

## Right: RPKM density
p_expr_rpkm <- ggplot(expr_plot_df,
                      aes(x = log2(mean_rpkm + 1),
                          color = expr_group, fill = expr_group)) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.2) +
  scale_color_manual(values = expr_colors) +
  scale_fill_manual(values  = expr_colors) +
  labs(x = NULL, y = "Density",
       title    = "RPKM distribution by expression group",
       subtitle = "EU-seq  |  mean across all 6 replicates") +
  base_theme +
  theme(legend.position  = "top",
        legend.title     = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

## Right: RPKM boxplot strip
p_expr_rpkm_box <- ggplot(expr_plot_df,
                          aes(x = log2(mean_rpkm + 1), y = expr_group,
                              color = expr_group, fill = expr_group)) +
  geom_boxplot(alpha = 0.15, linewidth = 0.6, outlier.size = 0.3,
               outlier.alpha = 0.2, width = 0.5) +
  scale_color_manual(values = expr_colors) +
  scale_fill_manual(values  = expr_colors) +
  labs(x = "log2(mean RPKM + 1)", y = NULL) +
  base_theme +
  theme(legend.position = "none",
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        axis.line.y     = element_blank())

## Combine
p_left  <- p_expr_lfc  / p_expr_lfc_box  + plot_layout(heights = c(5, 1))
p_right <- p_expr_rpkm / p_expr_rpkm_box + plot_layout(heights = c(5, 1))

print(p_left | p_right)

##--------------------------------------------------------------## 8e. — Volcano plot

message("Plot 2: Volcano plot...")

label_vol <- combined_df %>%
  dplyr::filter(de_status %in% c("Up", "Down")) %>%
  dplyr::arrange(desc(abs(lfc))) %>%
  dplyr::slice_head(n = N_LABELS)

p_vol <- ggplot(
  combined_df,
  aes(x = lfc, y = z_score, color = de_status)
) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = c(-2, 2),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = label_vol, aes(label = gene_id),
    size = 2.8, max.overlaps = 20,
    segment.size = 0.3, segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(values = de_colors) +
  labs(
    x        = "log2 fold change (POLR2A / WT)",
    y        = "Z-score (effect size)",
    title    = "EU-seq LFC vs effect size — POLR2A degradation (nascent RNA)",
    subtitle = paste0("|LFC| > ", LFC_THRESHOLD,
                      "  |  Top ", N_LABELS, " genes labeled by |LFC|")
  ) +
  base_theme

print(p_vol)

##--------------------------------------------------------------## 11. Plot 4 — LFC scatter: EU-seq vs total RNA-seq

message("Plot 4: LFC scatter...")
scatter_df <- combined_df %>%
  dplyr::filter(!is.na(lfc), !is.na(lfc_total)) %>%
  dplyr::mutate(
    concordance = dplyr::case_when(
      de_status == "Up"   & de_status_total == "Up"   ~ "Up in both",
      de_status == "Down" & de_status_total == "Down" ~ "Down in both",
      de_status %in% c("Up","Down") &
        de_status_total %in% c("Up","Down")           ~ "Discordant",
      de_status_total %in% c("Up", "Down")            ~ "Steady-state only",
      de_status %in% c("Up", "Down")                  ~ "Nascent only",
      TRUE                                             ~ "Not DE"
    ),
    concordance = factor(concordance,
                         levels = c("Up in both", "Down in both",
                                    "Discordant", "Steady-state only",
                                    "Nascent only", "Not DE"))
  )

## Add these two lines here
n_up_both   <- sum(scatter_df$concordance == "Up in both",   na.rm = TRUE)
n_down_both <- sum(scatter_df$concordance == "Down in both", na.rm = TRUE)

cat("Up in both:  ", n_up_both,   "\n")
cat("Down in both:", n_down_both, "\n")

r_all <- cor(scatter_df$lfc, scatter_df$lfc_total,
             use = "complete.obs", method = "pearson")
cat("Pearson r (EU-seq vs total RNA-seq LFC):", round(r_all, 3), "\n")

conc_colors <- c(
  "Up in both"      = COL_UP,
  "Down in both"    = COL_DOWN,
  "Discordant"      = "#CC6677",
  "Steady-state only" = "#E69F00",
  "Nascent only"    = "#88CCEE",
  "Not DE"          = COL_NS
)

label_scat <- scatter_df %>%
  dplyr::filter(concordance %in% c("Up in both", "Down in both", "Discordant")) %>%
  dplyr::arrange(desc(abs(lfc))) %>%
  dplyr::slice_head(n = N_LABELS)

library(ggside)

margin_df <- scatter_df %>%
  dplyr::filter(concordance %in% c("Up in both", "Down in both",
                                   "Steady-state only", "Nascent only"))

## Equal axes: one symmetric range and one tick set shared by x and y
AX_LIM    <- ceiling(max(abs(c(scatter_df$lfc_total, scatter_df$lfc)), na.rm = TRUE))
AX_TICK   <- floor(AX_LIM / 2) * 2
AX_BREAKS <- seq(-AX_TICK, AX_TICK, by = 2)
cat("Scatter axis range: +/-", AX_LIM, " | ticks:", AX_BREAKS, "\n")

p_scatter <- ggplot(scatter_df,
                    aes(x = lfc_total, y = lfc,
                        color = concordance,
                        fill  = concordance)) +
  geom_point(size = 0.7, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  ggside::geom_xsideboxplot(data = margin_df,
                            aes(y = concordance),
                            orientation = "y",
                            outlier.size = 0.3,
                            outlier.alpha = 0.3,
                            linewidth = 0.4) +
  ggside::geom_ysideboxplot(data = margin_df,
                            aes(x = concordance),
                            orientation = "x",
                            outlier.size = 0.3,
                            outlier.alpha = 0.3,
                            linewidth = 0.4) +
  ggside::scale_xsidey_discrete() +
  ggside::scale_ysidex_discrete(guide = guide_axis(angle = 90)) +
  ggside::theme_ggside_void() +
  scale_color_manual(values = conc_colors) +
  scale_fill_manual(values  = conc_colors) +
  scale_x_continuous(breaks = AX_BREAKS) +
  scale_y_continuous(breaks = AX_BREAKS) +
  coord_cartesian(xlim = c(-AX_LIM, AX_LIM), ylim = c(-AX_LIM, AX_LIM)) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Up in both: n=", n_up_both),
           hjust = 1.1, vjust = 2.0,
           size = 3.5, color = COL_UP, fontface = "bold") +
  annotate("text",
           x = -Inf, y = -Inf,
           label = paste0("Down in both: n=", n_down_both),
           hjust = -0.1, vjust = -1.5,
           size = 3.5, color = COL_DOWN, fontface = "bold") +
  annotate("text",
           x = Inf, y = -Inf,
           label = paste0("r = ", round(r_all, 3)),
           hjust = 1.1, vjust = -0.5,
           size = 4, color = "grey30", fontface = "italic") +
  labs(
    x        = "Total RNA-seq log2 FC (POLR2A / WT)",
    y        = "EU-seq log2 FC (POLR2A / WT)",
    title    = "Nascent vs steady-state fold change",
    subtitle = paste0("Diagonal = perfect concordance  |  r = ",
                      round(r_all, 3))
  ) +
  base_theme +
  theme(ggside.panel.scale = 0.2,
        aspect.ratio = 1)

print(p_scatter)

##--------------------------------------------------------------## 12. Plot 5 — Boxplot LFC by DEG category x E/I group

message("Plot 5: Boxplot...")

box_df <- combined_df %>%
  dplyr::filter(!is.na(lfc), !is.na(ei_group),
                de_status_total %in% c("Up", "Not DE", "Down")) %>%
  dplyr::mutate(
    de_status_total = factor(de_status_total,
                             levels = c("Up", "Not DE", "Down")),
    ei_group = factor(ei_group,
                      levels = c("Low E/I", "Medium E/I", "High E/I"))
  )

## Compute n per category
box_counts <- box_df %>%
  dplyr::group_by(ei_group, de_status_total) %>%
  dplyr::summarise(n = dplyr::n(),
                   label = paste0(dplyr::n()),
                   .groups = "drop")

## Place labels at the median LFC per group so they sit inside the box
box_medians <- box_df %>%
  dplyr::group_by(ei_group, de_status_total) %>%
  dplyr::summarise(median_lfc = median(lfc, na.rm = TRUE),
                   .groups = "drop")

box_counts <- box_counts %>%
  dplyr::left_join(box_medians, by = c("ei_group", "de_status_total"))

## Publication styling (Figure 3E): colors, y-axis marks, strip labels
BOX_COLS  <- c(Up = "#97c9e8", "Not DE" = "#b3b3b3", Down = "#be6b78")
BOX_YLIM  <- c(-5.25, 5.25)          ## clipped view; boxes/outliers beyond are cut, not dropped
BOX_YTICK <- c(-5, 0, 5)
BOX_LW    <- 0.7                     ## box, whisker, axis and strip line width
ei_strip_labels <- c("Low E/I"    = "Intron-rich",
                     "Medium E/I" = "Average E/I",
                     "High E/I"   = "Exon-rich")

p_box <- ggplot(box_df,
                aes(x = de_status_total, y = lfc,
                    fill = de_status_total)) +
  geom_boxplot(outlier.size = 0.7, outlier.alpha = 0.5, outlier.colour = "grey50",
               color = "grey15", linewidth = BOX_LW, width = 0.8, fatten = 1) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_text(data = box_counts,
            aes(x     = de_status_total,
                y     = median_lfc,
                label = label),
            color    = "white",
            size     = 3.3,
            fontface = "bold",
            inherit.aes = FALSE) +
  scale_fill_manual(values = BOX_COLS) +
  scale_y_continuous(breaks = BOX_YTICK) +
  coord_cartesian(ylim = BOX_YLIM) +
  facet_wrap(~ ei_group, nrow = 1, labeller = labeller(ei_group = ei_strip_labels)) +
  labs(
    x = "Total RNA-seq DEG status",
    y = "EU-seq log2 FC (POLR2A / WT)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "none",
    axis.line        = element_line(color = "black", linewidth = BOX_LW),
    axis.ticks       = element_line(color = "black", linewidth = BOX_LW),
    axis.text        = element_text(color = "black", size = 13),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    axis.title       = element_text(color = "black", size = 15),
    strip.background = element_rect(color = "black", fill = "white", linewidth = BOX_LW),
    strip.text       = element_text(color = "black", size = 14),
    panel.spacing.x  = unit(9, "pt")
  )

print(p_box)   ## publication panel is 7.3 x 5.9 in (730 x 594 px at 100 dpi)

##--------------------------------------------------------------## 13. Plot 6 — DEG overlap bar chart

message("Plot 6: DEG overlap...")

overlap_df <- combined_df %>%
  dplyr::filter(!is.na(lfc), !is.na(lfc_total)) %>%
  dplyr::mutate(
    category = dplyr::case_when(
      de_status %in% c("Up","Down") &
        de_status_total %in% c("Up","Down")  ~ "Both nascent\n& steady-state",
      de_status %in% c("Up","Down") &
        !de_status_total %in% c("Up","Down") ~ "Nascent only",
      !de_status %in% c("Up","Down") &
        de_status_total %in% c("Up","Down")  ~ "Steady-state\nonly",
      TRUE                                   ~ "Not DE in either"
    )
  )

overlap_counts <- overlap_df %>%
  dplyr::count(category) %>%
  dplyr::mutate(
    pct      = n / sum(n) * 100,
    category = factor(category,
                      levels = c("Both nascent\n& steady-state",
                                 "Nascent only", "Steady-state\nonly",
                                 "Not DE in either"))
  )

cat("\nDEG overlap counts:\n")
print(overlap_counts)


overlap_colors <- c(
  "Both nascent\n& steady-state" = "#009E73",
  "Nascent only"                 = "#E69F00",
  "Steady-state\nonly"           = "#0072B2",
  "Not DE in either"             = COL_NS
)

p_overlap <- ggplot(
  overlap_counts %>% dplyr::filter(category != "Not DE in either"),
  aes(x = category, y = n, fill = category)
) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(n, "\n(", round(pct, 1), "%)")),
            vjust = -0.3, size = 3.5, color = "grey20") +
  scale_fill_manual(values = overlap_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x        = NULL,
    y        = "Number of genes",
    title    = "DEG overlap: nascent vs steady-state",
    subtitle = paste0("|LFC| > ", LFC_THRESHOLD)
  ) +
  base_theme +
  theme(legend.position = "none")

print(p_overlap)

##--------------------------------------------------------------## 14. Plot 7 — Cumulative LFC distribution

message("Plot 7: Cumulative LFC distribution...")

cdf_df <- combined_df %>%
  dplyr::filter(!is.na(lfc), !is.na(lfc_total),
                de_status_total %in% c("Up", "Not DE", "Down")) %>%
  dplyr::select(gene_id, lfc_eu = lfc, lfc_total, de_status_total) %>%
  tidyr::pivot_longer(cols = c(lfc_eu, lfc_total),
                      names_to = "assay", values_to = "lfc") %>%
  dplyr::mutate(
    assay = dplyr::recode(assay,
                          lfc_eu    = "EU-seq (nascent)",
                          lfc_total = "Total RNA-seq (steady-state)"),
    de_status_total = factor(de_status_total,
                             levels = c("Up", "Not DE", "Down"))
  )

p_cdf <- ggplot(cdf_df, aes(x = lfc, color = assay, linetype = assay)) +
  stat_ecdf(linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dotted",
             color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = c(
    "EU-seq (nascent)"             = "#E69F00",
    "Total RNA-seq (steady-state)" = "#0072B2"
  )) +
  facet_wrap(~ de_status_total, nrow = 1) +
  labs(
    x        = "log2 fold change (POLR2A / WT)",
    y        = "Cumulative fraction",
    title    = "Cumulative LFC distribution: nascent vs steady-state",
    subtitle = "Faceted by steady-state DEG status"
  ) +
  base_theme

print(p_cdf)

##--------------------------------------------------------------## 15. Plot 8 — RRHO-style rank correlation heatmap

message("Plot 8: Rank correlation heatmap...")

rrho_df <- combined_df %>%
  dplyr::filter(!is.na(lfc), !is.na(lfc_total)) %>%
  dplyr::mutate(
    rank_eu    = rank(lfc,       ties.method = "average") / n(),
    rank_total = rank(lfc_total, ties.method = "average") / n()
  )

N_BINS_RRHO <- 20
rrho_grid <- rrho_df %>%
  dplyr::mutate(
    bin_eu    = cut(rank_eu,    breaks = N_BINS_RRHO, labels = FALSE),
    bin_total = cut(rank_total, breaks = N_BINS_RRHO, labels = FALSE)
  ) %>%
  dplyr::group_by(bin_eu, bin_total) %>%
  dplyr::summarise(n_genes = dplyr::n(), .groups = "drop")

p_rrho <- ggplot(rrho_grid,
                 aes(x = bin_total, y = bin_eu, fill = n_genes)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(9, "RdYlBu")),
    name   = "Genes\nper bin"
  ) +
  scale_x_continuous(breaks = c(1, N_BINS_RRHO/2, N_BINS_RRHO),
                     labels = c("Down", "No change", "Up")) +
  scale_y_continuous(breaks = c(1, N_BINS_RRHO/2, N_BINS_RRHO),
                     labels = c("Down", "No change", "Up")) +
  labs(
    x        = "Total RNA-seq rank (POLR2A / WT)",
    y        = "EU-seq rank (POLR2A / WT)",
    title    = "Rank-rank concordance: nascent vs steady-state",
    subtitle = paste0("Genome-wide LFC concordance  |  ",
                      nrow(rrho_df), " genes  |  ",
                      N_BINS_RRHO, "×", N_BINS_RRHO, " grid")
  ) +
  base_theme +
  theme(legend.position = "right")

print(p_rrho)

message("All plots complete.")