##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# QC Pearson correlation heatmap and PCA for EU-seq replicates
# using raw (unsubtracted) counts from featureCounts.
#
# This is a pre-subtraction QC check to verify replicate
# consistency at the raw data level, independent of the
# NoEU background subtraction step.
#
# 7 samples: POLR2A (BM23-25), WT (BM29-31), NoEU (S49)
# Paired-end, reverse stranded (strandSpecific = 2)
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(Rsubread)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
})

##--------------------------------------------------------------## 1. Parameters

BAM_DIR  <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/BAM"
GTF_PATH <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
N_THREADS <- 7
MIN_CPM   <- 1   ## genes with CPM < MIN_CPM in all samples are dropped

BAM_FILES <- c(
  POLR2A_1 = file.path(BAM_DIR, "BM23_S40_filt.sort.bam"),
  POLR2A_2 = file.path(BAM_DIR, "BM24_S41_filt.sort.bam"),
  POLR2A_3 = file.path(BAM_DIR, "BM25_S42_filt.sort.bam"),
  WT_1     = file.path(BAM_DIR, "BM29_S46_filt.sort.bam"),
  WT_2     = file.path(BAM_DIR, "BM30_S47_filt.sort.bam"),
  WT_3     = file.path(BAM_DIR, "BM31_S48_filt.sort.bam"),
  NoEU     = file.path(BAM_DIR, "No_EU_ctrl_S49_filt.sort.bam")
)

CONDITION_MAP <- c(
  POLR2A_1 = "POLR2A", POLR2A_2 = "POLR2A", POLR2A_3 = "POLR2A",
  WT_1     = "WT",     WT_2     = "WT",     WT_3     = "WT",
  NoEU     = "NoEU"
)

CONDITION_COLORS <- c(
  POLR2A = "#7B2D8B",
  WT     = "#1B7837",
  NoEU   = "#999999"
)

##--------------------------------------------------------------## 2. Validate paths

message("Checking paths...")
missing <- BAM_FILES[!file.exists(BAM_FILES)]
if (length(missing) > 0) stop("Missing BAMs:\n",
                              paste(" -", missing, collapse = "\n"))
if (!file.exists(GTF_PATH)) stop("GTF not found: ", GTF_PATH)
message("  All files found.")

##--------------------------------------------------------------## 3. featureCounts
## Paired-end, reverse stranded, count over exons
## EU-seq counts over full gene body (exons + introns) would be
## ideal but featureCounts operates on GTF features — using exons
## here is standard for QC purposes since we just want to check
## replicate consistency, not compute FPKM

message("Running featureCounts (gene-level, full gene body including introns)...")
fc <- featureCounts(
  files               = unname(BAM_FILES),
  annot.ext           = GTF_PATH,
  isGTFAnnotationFile = TRUE,
  GTF.featureType     = "gene",      ## full gene body TSS to TES
  GTF.attrType        = "gene_id",
  useMetaFeatures     = TRUE,
  isPairedEnd         = TRUE,
  strandSpecific      = 2,           ## reverse stranded
  nthreads            = N_THREADS,
  countMultiMappingReads = FALSE, ## was false
  verbose             = FALSE
)

counts <- fc$counts
colnames(counts) <- names(BAM_FILES)
counts <- counts[rowSums(counts) > 0, ]

cat("featureCounts complete.\n")
cat("  Genes with counts:", nrow(counts), "\n")
cat("  Samples:          ", ncol(counts), "\n\n")
cat("Assignment summary:\n")
print(fc$stat)

##--------------------------------------------------------------## 4. Normalize to log2(CPM + 1)

message("Normalizing to log2(CPM + 1)...")
lib_sizes <- colSums(counts)
cat("\nLibrary sizes:\n")
print(lib_sizes)

cpm      <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
keep     <- rowSums(cpm >= MIN_CPM) >= 0
cpm      <- cpm[keep, ]
log2cpm  <- log2(cpm + 1)

cat("\nGenes passing CPM filter (>= ", MIN_CPM, " in >= 1 sample):",
    nrow(log2cpm), "\n")

##--------------------------------------------------------------## 5. Pearson correlation heatmap

message("Plotting correlation heatmap...")
cor_mat <- cor(log2cpm, method = "pearson")

cat("\nPearson correlation matrix:\n")
print(round(cor_mat, 4))

## Within-condition mean r
within_cors <- sapply(unique(CONDITION_MAP), function(cond) {
  samps <- names(CONDITION_MAP)[CONDITION_MAP == cond]
  if (length(samps) < 2) return(NA)
  vals  <- cor_mat[samps, samps]
  mean(vals[upper.tri(vals)])
})
cat("\nMean within-condition Pearson r:\n")
print(round(within_cors, 4))

annotation_col <- data.frame(
  Condition = CONDITION_MAP[colnames(cor_mat)],
  row.names = colnames(cor_mat)
)
annotation_colors <- list(Condition = CONDITION_COLORS)

cor_min <- floor(min(cor_mat) * 100) / 100
breaks  <- seq(cor_min, 1, length.out = 101)
colors  <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(9, "RdYlBu"))
)(100)

pheatmap(
  cor_mat,
  color             = colors,
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
  cellwidth         = 35,
  cellheight        = 35,
  border_color      = NA,
  main              = paste0(
    "EU-seq — Pearson correlation (log2 CPM+1, unsubtracted)\n",
    "Gene-level counts (TSS to TES, introns included)  |  ",
    nrow(log2cpm), " genes  |  7 samples"
  )
)

##--------------------------------------------------------------## 6. PCA

message("Plotting PCA...")
pca_res <- prcomp(t(log2cpm), center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x[, 1:3]) %>%
  dplyr::mutate(
    sample    = rownames(.),
    condition = CONDITION_MAP[sample]
  )

pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

base_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11),
    legend.title  = element_blank()
  )

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
    subtitle = "log2(CPM+1)  |  gene-level (TSS-TES)  |  unsubtracted  |  7 samples"
  ) +
  base_theme

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
    subtitle = "log2(CPM+1)  |  gene-level (TSS-TES)  |  unsubtracted  |  7 samples"
  ) +
  base_theme

print(p_pca12 | p_pca13)

message("EU-seq QC plots complete.")
