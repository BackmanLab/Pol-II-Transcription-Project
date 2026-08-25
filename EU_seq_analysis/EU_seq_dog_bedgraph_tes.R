##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# DoG metagene analysis — profiles EU-seq signal in a window
# around TES (-2kb to +10kb) to show transcriptional readthrough.
#
# Uses subtracted RPM-normalized bedgraphs directly (no bigWig).
# Signal extracted per gene, averaged across replicates per
# condition, then plotted as metagene profiles stratified by
# gene body RPKM expression quartile.
#
# Strandedness handled explicitly:
#   + strand genes: use plus bedgraph
#   - strand genes: use minus bedgraph, flip window coordinates
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(BiocGenerics)
  library(BiocParallel)
})

##--------------------------------------------------------------## 1. Parameters

PIPELINE_DIR <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline"
SUB_DIR      <- file.path(PIPELINE_DIR, "subtracted")
GENE_BED     <- file.path(PIPELINE_DIR, "hg38_genes.bed")
GENE_MAT     <- file.path(PIPELINE_DIR, "counts/euseq_rpkm_analysis.txt")
GTF_PATH     <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"

UPSTREAM    <- 2000    ## bp upstream of TES
DOWNSTREAM  <- 10000  ## bp downstream of TES
BIN_SIZE    <- 100    ## bp per bin
RPKM_THRESHOLD <- 0.5
N_CORES     <- 4

POLR2A_SAMPLES <- c("POLR2A_1","POLR2A_2","POLR2A_3")
WT_SAMPLES     <- c("WT_1","WT_2","WT_3")
ALL_SAMPLES    <- c(POLR2A_SAMPLES, WT_SAMPLES)

CONDITION_COLORS <- c(POLR2A = "#7B2D8B", WT = "#1B7837")

TOTAL_RNASEQ_SF <- c(
  POLR2A_1 = 0.747437, POLR2A_2 = 0.747437, POLR2A_3 = 0.747437,
  WT_1     = 1.000000, WT_2     = 1.000000, WT_3     = 1.000000
)

base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey60"),
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

##--------------------------------------------------------------## 2. Load gene BED and expression data

message("Loading gene annotations...")

genes_df <- read_tsv(GENE_BED,
                     col_names = c("chr","start","end",
                                   "gene_id","score","strand"),
                     show_col_types = FALSE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id))

cat("Genes in BED:", nrow(genes_df), "\n")

## Load gene body RPKM for expression stratification
to_matrix <- function(df) {
  mat <- as.matrix(df[, -1])
  rownames(mat) <- sub("\\..*", "", df[[1]])
  mode(mat) <- "numeric"
  mat
}

gene_mat    <- to_matrix(read_tsv(GENE_MAT, show_col_types = FALSE))
gene_scaled <- sweep(gene_mat, 2, TOTAL_RNASEQ_SF[colnames(gene_mat)],
                     FUN = "*")
wt_mean     <- rowMeans(gene_scaled[, WT_SAMPLES])

## Assign expression quartiles
expressed   <- wt_mean[wt_mean >= RPKM_THRESHOLD]
expr_df <- data.frame(
  gene_id    = names(wt_mean),
  wt_rpkm    = as.numeric(wt_mean),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(wt_rpkm >= RPKM_THRESHOLD) %>%
  dplyr::mutate(
    expr_group = dplyr::case_when(
      wt_rpkm < quantile(expressed, 0.25, na.rm = TRUE) ~ "Q1 (Low)",
      wt_rpkm < quantile(expressed, 0.50, na.rm = TRUE) ~ "Q2",
      wt_rpkm < quantile(expressed, 0.75, na.rm = TRUE) ~ "Q3",
      TRUE                                               ~ "Q4 (High)"
    ),
    expr_group = factor(expr_group,
                        levels = c("Q1 (Low)","Q2","Q3","Q4 (High)"))
  )

cat("Genes passing RPKM threshold:", nrow(expr_df), "\n")
print(table(expr_df$expr_group))

## Join expression groups to gene BED
genes_df <- genes_df %>%
  dplyr::inner_join(expr_df %>% dplyr::select(gene_id, expr_group),
                    by = "gene_id")

cat("Genes with expression annotation:", nrow(genes_df), "\n")

##--------------------------------------------------------------## 2b. Compute gene lengths from GTF

message("Computing gene lengths from GTF...")

txdb         <- GenomicFeatures::makeTxDbFromGFF(GTF_PATH,
                                                 format = "gtf",
                                                 taxonomyId = 9606)
genes_gtf    <- GenomicFeatures::genes(txdb)
gene_len_df  <- data.frame(
  gene_id    = sub("\\..*", "", genes_gtf$gene_id),
  gene_length = BiocGenerics::width(genes_gtf),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    length_cat = dplyr::case_when(
      gene_length <  1000  ~ "Short (<1kb)",
      gene_length <= 10000 ~ "Average (1-10kb)",
      TRUE                 ~ "Long (>10kb)"
    ),
    length_cat = factor(length_cat,
                        levels = c("Short (<1kb)",
                                   "Average (1-10kb)",
                                   "Long (>10kb)"))
  )

cat("\nGene length categories:\n")
print(table(gene_len_df$length_cat))

## Join length categories to genes_df
genes_df <- genes_df %>%
  dplyr::left_join(
    gene_len_df %>% dplyr::select(gene_id, gene_length, length_cat),
    by = "gene_id"
  )

##--------------------------------------------------------------## 3. Define TES windows
## For + strand genes: TES = end coordinate
##   window = [TES - UPSTREAM, TES + DOWNSTREAM]
## For - strand genes: TES = start coordinate
##   window = [TES - DOWNSTREAM, TES + UPSTREAM]
## Both in genomic coordinates — strand handled during extraction

n_bins  <- (UPSTREAM + DOWNSTREAM) / BIN_SIZE
bin_pos <- seq(-UPSTREAM, DOWNSTREAM - BIN_SIZE, by = BIN_SIZE) +
  BIN_SIZE / 2   ## bin centers relative to TES

genes_gr <- genes_df %>%
  dplyr::mutate(
    tes = ifelse(strand == "+", end, start),
    win_start = ifelse(strand == "+",
                       tes - UPSTREAM,
                       tes - DOWNSTREAM),
    win_end   = ifelse(strand == "+",
                       tes + DOWNSTREAM,
                       tes + UPSTREAM),
    win_start = pmax(win_start, 0)
  ) %>%
  dplyr::filter(win_end > win_start)

cat("Genes with valid TES windows:", nrow(genes_gr), "\n")

##--------------------------------------------------------------## 4. Import bedgraphs
## Import once per strand per sample — cache in memory
## Each bedgraph covers the whole genome so import once and
## query multiple times

message("Importing bedgraphs (this may take a few minutes)...")

import_bedgraph <- function(path) {
  message("  Loading: ", basename(path))
  gr <- rtracklayer::import.bedGraph(path)
  ## Ensure score column exists
  if (!"score" %in% names(mcols(gr))) {
    mcols(gr)$score <- mcols(gr)[[1]]
  }
  gr
}

## Import all bedgraphs — store in named list
bg_list <- list()
for (samp in ALL_SAMPLES) {
  bg_list[[paste0(samp, "_plus")]]  <-
    import_bedgraph(file.path(SUB_DIR, paste0(samp, "_plus.sub.bgr")))
  bg_list[[paste0(samp, "_minus")]] <-
    import_bedgraph(file.path(SUB_DIR, paste0(samp, "_minus.sub.bgr")))
}

message("All bedgraphs loaded.")

##--------------------------------------------------------------## 5. Extract signal in TES windows
## For each gene, extract coverage from the appropriate strand
## bedgraph in BIN_SIZE bins across the TES window

extract_gene_signal <- function(gene_row, bg_plus, bg_minus) {
  strand  <- gene_row$strand
  chr     <- gene_row$chr
  ws      <- gene_row$win_start
  we      <- gene_row$win_end
  bg      <- if (strand == "+") bg_plus else bg_minus
  
  ## Create bin ranges
  bins <- GRanges(
    seqnames = chr,
    ranges   = IRanges(
      start = seq(ws, we - BIN_SIZE, by = BIN_SIZE) + 1,
      width = BIN_SIZE
    )
  )
  
  ## Find overlapping bedgraph intervals
  hits  <- findOverlaps(bins, bg)
  score <- rep(0, length(bins))
  
  if (length(hits) > 0) {
    for (i in seq_along(bins)) {
      ovlp_idx <- subjectHits(hits[queryHits(hits) == i])
      if (length(ovlp_idx) == 0) next
      ovlp_gr  <- bg[ovlp_idx]
      ## Weighted mean by overlap length
      ovlp_len <- width(pintersect(bins[i], ovlp_gr,
                                   ignore.strand = TRUE))
      score[i] <- sum(mcols(ovlp_gr)$score * ovlp_len) /
        BIN_SIZE
    }
  }
  
  ## For minus strand genes flip so position 0 is TES
  if (strand == "-") score <- rev(score)
  score
}

## Extract signal for all genes x all samples
message("Extracting signal in TES windows (parallel, ", N_CORES, " cores)...")

BiocParallel::register(BiocParallel::MulticoreParam(workers = N_CORES))

extract_sample_signal <- function(samp) {
  message("  Processing: ", samp)
  bg_plus  <- bg_list[[paste0(samp, "_plus")]]
  bg_minus <- bg_list[[paste0(samp, "_minus")]]
  
  signal_list <- BiocParallel::bplapply(
    seq_len(nrow(genes_gr)),
    function(i) {
      tryCatch(
        extract_gene_signal(genes_gr[i, ], bg_plus, bg_minus),
        error = function(e) rep(0, n_bins)
      )
    }
  )
  
  mat <- do.call(rbind, signal_list)
  rownames(mat) <- genes_gr$gene_id
  cat("  Done:", samp, "— mean signal:",
      round(mean(mat, na.rm = TRUE), 4), "\n")
  mat
}

results <- lapply(ALL_SAMPLES, function(samp) {
  extract_sample_signal(samp)
})
names(results) <- ALL_SAMPLES

message("Signal extraction complete.")

##--------------------------------------------------------------## 6. Average across replicates per condition

message("Averaging across replicates...")

## Mean across POLR2A replicates
polr2a_mean <- (results[["POLR2A_1"]] +
                  results[["POLR2A_2"]] +
                  results[["POLR2A_3"]]) / 3

## Mean across WT replicates
wt_mean_mat <- (results[["WT_1"]] +
                  results[["WT_2"]] +
                  results[["WT_3"]]) / 3

##--------------------------------------------------------------## 7. Build metagene data frame

message("Building metagene profiles...")

build_metagene <- function(signal_mat, condition, genes_df) {
  ## signal_mat: genes x bins
  ## Join with expression groups
  gene_ids <- rownames(signal_mat)
  
  expr_map <- genes_df %>%
    dplyr::select(gene_id, expr_group) %>%
    dplyr::filter(gene_id %in% gene_ids)
  
  result <- lapply(levels(expr_map$expr_group), function(grp) {
    grp_genes <- expr_map$gene_id[expr_map$expr_group == grp]
    grp_genes <- intersect(grp_genes, gene_ids)
    if (length(grp_genes) == 0) return(NULL)
    
    grp_mat  <- signal_mat[grp_genes, , drop = FALSE]
    grp_mean <- colMeans(grp_mat, na.rm = TRUE)
    
    data.frame(
      bin_pos    = bin_pos,
      mean_signal = grp_mean,
      expr_group = grp,
      condition  = condition,
      n_genes    = length(grp_genes),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(result)
}

meta_polr2a <- build_metagene(polr2a_mean, "POLR2A", genes_gr)
meta_wt     <- build_metagene(wt_mean_mat, "WT",     genes_gr)

meta_df <- dplyr::bind_rows(meta_polr2a, meta_wt) %>%
  dplyr::mutate(
    condition  = factor(condition, levels = c("WT","POLR2A")),
    expr_group = factor(expr_group,
                        levels = c("Q1 (Low)","Q2","Q3","Q4 (High)"))
  )

cat("Metagene data frame:", nrow(meta_df), "rows\n")

## Gene counts per group
meta_df %>%
  dplyr::distinct(expr_group, condition, n_genes) %>%
  print()

##--------------------------------------------------------------## 8. Plot DoG metagene

message("Plotting DoG metagene...")

## Expression group colors
EXPR_COLORS <- c(
  "Q1 (Low)"  = "#0072B2",
  "Q2"        = "#56B4E9",
  "Q3"        = "#E69F00",
  "Q4 (High)" = "#D55E00"
)

## --- Plot 1: Faceted by expression group, two lines per facet ---
p_meta_expr <- ggplot(meta_df,
                      aes(x = bin_pos / 1000,
                          y = mean_signal,
                          color = condition,
                          linetype = condition)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "grey40",
             linetype = "dashed", linewidth = 0.5) +
  annotate("rect",
           xmin = -UPSTREAM/1000, xmax = 0,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.3) +
  scale_color_manual(values = CONDITION_COLORS) +
  scale_linetype_manual(values = c(WT = "solid", POLR2A = "solid")) +
  scale_x_continuous(
    breaks = c(-2, 0, 2, 6, 10),
    labels = c("-2kb","TES","+2kb",
               "+6kb","+10kb")
  ) +
  facet_wrap(~ expr_group, nrow = 1, scales = "free_y") +
  labs(
    x        = "Position relative to TES",
    y        = "Mean RPM signal",
    title    = "DoG metagene — EU-seq signal around TES",
    subtitle = paste0("Window: -", UPSTREAM/1000, "kb to +",
                      DOWNSTREAM/1000, "kb  |  Bin size: ",
                      BIN_SIZE, "bp  |  Grey = gene body region")
  ) +
  base_theme

print(p_meta_expr)

## --- Plot 2: Faceted by condition, four lines per facet ---
p_meta_cond <- ggplot(meta_df,
                      aes(x = bin_pos / 1000,
                          y = mean_signal,
                          color = expr_group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "grey40",
             linetype = "dashed", linewidth = 0.5) +
  annotate("rect",
           xmin = -UPSTREAM/1000, xmax = 0,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.3) +
  scale_color_manual(values = EXPR_COLORS) +
  scale_x_continuous(
    breaks = c(-2, 0, 2, 6, 10),
    labels = c("-2kb","TES","+2kb",
               "+6kb","+10kb")
  ) +
  facet_wrap(~ condition, nrow = 1) +
  labs(
    x        = "Position relative to TES",
    y        = "Mean RPM signal",
    title    = "DoG metagene — EU-seq signal around TES",
    subtitle = paste0("Window: -", UPSTREAM/1000, "kb to +",
                      DOWNSTREAM/1000, "kb  |  ",
                      "Colored by expression quartile")
  ) +
  base_theme

print(p_meta_cond)

##--------------------------------------------------------------## 8b. Gene length metagene plots

message("Building gene length metagene profiles...")

## Build metagene stratified by gene length category
build_metagene_length <- function(signal_mat, condition, genes_df) {
  gene_ids <- rownames(signal_mat)
  
  len_map <- genes_df %>%
    dplyr::select(gene_id, length_cat) %>%
    dplyr::filter(!is.na(length_cat), gene_id %in% gene_ids)
  
  result <- lapply(levels(len_map$length_cat), function(cat) {
    cat_genes <- len_map$gene_id[len_map$length_cat == cat]
    cat_genes <- intersect(cat_genes, gene_ids)
    if (length(cat_genes) == 0) return(NULL)
    
    cat_mat  <- signal_mat[cat_genes, , drop = FALSE]
    cat_mean <- colMeans(cat_mat, na.rm = TRUE)
    
    data.frame(
      bin_pos     = bin_pos,
      mean_signal = cat_mean,
      length_cat  = cat,
      condition   = condition,
      n_genes     = length(cat_genes),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(result)
}

meta_len_polr2a <- build_metagene_length(polr2a_mean, "POLR2A", genes_df)
meta_len_wt     <- build_metagene_length(wt_mean_mat, "WT",     genes_df)

meta_len_df <- dplyr::bind_rows(meta_len_polr2a, meta_len_wt) %>%
  dplyr::mutate(
    condition  = factor(condition, levels = c("WT","POLR2A")),
    length_cat = factor(length_cat,
                        levels = c("Short (<1kb)",
                                   "Average (1-10kb)",
                                   "Long (>10kb)"))
  )

cat("\nGene length metagene — gene counts:\n")
meta_len_df %>%
  dplyr::distinct(length_cat, condition, n_genes) %>%
  print()

## Length category colors
LENGTH_COLORS <- c(
  "Short (<1kb)"     = "#0072B2",
  "Average (1-10kb)" = "#E69F00",
  "Long (>10kb)"     = "#D55E00"
)

## --- Plot: Faceted by gene length, WT vs POLR2A lines ---
p_meta_length <- ggplot(meta_len_df,
                        aes(x = bin_pos / 1000,
                            y = mean_signal,
                            color = condition)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "grey40",
             linetype = "dashed", linewidth = 0.5) +
  annotate("rect",
           xmin = -UPSTREAM/1000, xmax = 0,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.3) +
  scale_color_manual(values = CONDITION_COLORS) +
  scale_x_continuous(
    breaks = c(-2,  0, 2, 6,  10),
    labels = c("-2kb","TES","+2kb",
               "+6kb","+10kb")
  ) +
  facet_wrap(~ length_cat, nrow = 1, scales = "free_y") +
  labs(
    x        = "Position relative to TES",
    y        = "Mean RPM signal",
    title    = "DoG metagene by gene length — EU-seq signal around TES",
    subtitle = paste0("Window: -", UPSTREAM/1000, "kb to +",
                      DOWNSTREAM/1000, "kb  |  Grey = gene body region")
  ) +
  base_theme

print(p_meta_length)

## --- Plot: Faceted by condition, three length category lines ---
p_meta_length_cond <- ggplot(meta_len_df,
                             aes(x = bin_pos / 1000,
                                 y = mean_signal,
                                 color = length_cat)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "grey40",
             linetype = "dashed", linewidth = 0.5) +
  annotate("rect",
           xmin = -UPSTREAM/1000, xmax = 0,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.3) +
  scale_color_manual(values = LENGTH_COLORS) +
  scale_x_continuous(
    breaks = c(-2, -1, 0, 2, 4, 6, 8, 10),
    labels = c("-2kb","-1kb","TES","+2kb","+4kb",
               "+6kb","+8kb","+10kb")
  ) +
  facet_wrap(~ condition, nrow = 1) +
  labs(
    x        = "Position relative to TES",
    y        = "Mean RPM signal",
    title    = "DoG metagene by gene length — EU-seq signal around TES",
    subtitle = "Colored by gene length category"
  ) +
  base_theme

print(p_meta_length_cond)


## Identify genes where DoG LFC (POLR2A / WT) >= 1
## Recompute metagene for this subset only

message("Computing high-DoG gene metagene...")

DOG_MAT <- file.path(PIPELINE_DIR,
                     "intron_dog/dog_counts/euseq_dog_rpkm.txt")
DOG_THRESHOLD <- 0.1   ## minimum mean WT DoG RPKM to include gene

## Load and scale DoG matrix
dog_mat    <- to_matrix(read_tsv(DOG_MAT, show_col_types = FALSE))
dog_scaled <- sweep(dog_mat, 2,
                    TOTAL_RNASEQ_SF[colnames(dog_mat)], FUN = "*")
dog_scaled <- dog_scaled[rowSums(dog_scaled) > 0, ]

## Compute DoG LFC per gene
wt_mean_dog     <- rowMeans(log2(dog_scaled[, WT_SAMPLES]     + 1))
polr2a_mean_dog <- rowMeans(log2(dog_scaled[, POLR2A_SAMPLES] + 1))
dog_lfc         <- polr2a_mean_dog - wt_mean_dog

## Filter to genes with DoG LFC >= 1 and detectable WT signal
high_dog_genes <- names(dog_lfc)[
  dog_lfc >= 1 &
    rowMeans(dog_scaled[names(dog_lfc), WT_SAMPLES]) >= DOG_THRESHOLD
]

cat("\nGenes with DoG LFC >= 1:", length(high_dog_genes), "\n")

if (length(high_dog_genes) < 10) {
  warning("Too few high-DoG genes — try lowering DOG_THRESHOLD")
} else {
  ## Subset signal matrices to high-DoG genes
  high_dog_genes_present <- intersect(high_dog_genes,
                                      rownames(polr2a_mean))
  
  cat("High-DoG genes with TES signal:", length(high_dog_genes_present), "\n")
  
  ## Average signal across high-DoG genes per condition
  high_dog_meta <- dplyr::bind_rows(
    data.frame(
      bin_pos     = bin_pos,
      mean_signal = colMeans(polr2a_mean[high_dog_genes_present, ,
                                         drop = FALSE], na.rm = TRUE),
      condition   = "POLR2A",
      stringsAsFactors = FALSE
    ),
    data.frame(
      bin_pos     = bin_pos,
      mean_signal = colMeans(wt_mean_mat[high_dog_genes_present, ,
                                         drop = FALSE], na.rm = TRUE),
      condition   = "WT",
      stringsAsFactors = FALSE
    )
  ) %>%
    dplyr::mutate(condition = factor(condition, levels = c("WT","POLR2A")))
  
  p_meta_high_dog <- ggplot(high_dog_meta,
                            aes(x     = bin_pos / 1000,
                                y     = mean_signal,
                                color = condition)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = 0, color = "grey40",
               linetype = "dashed", linewidth = 0.5) +
    annotate("rect",
             xmin = -UPSTREAM/1000, xmax = 0,
             ymin = -Inf, ymax = Inf,
             fill = "grey90", alpha = 0.3) +
    scale_color_manual(values = CONDITION_COLORS) +
    scale_x_continuous(
      breaks = c(-2, -1, 0, 2, 4, 6, 8, 10),
      labels = c("-2kb","-1kb","TES","+2kb","+4kb",
                 "+6kb","+8kb","+10kb")
    ) +
    labs(
      x        = "Position relative to TES",
      y        = "Mean RPM signal",
      title    = "DoG metagene — genes with DoG LFC >= 1",
      subtitle = paste0("n = ", length(high_dog_genes_present),
                        " genes  |  DoG RPKM threshold = ",
                        DOG_THRESHOLD,
                        "  |  Grey = gene body region")
    ) +
    base_theme
  
  print(p_meta_high_dog)
}

message("DoG metagene analysis complete.")
