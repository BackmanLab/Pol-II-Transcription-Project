##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Mito scaling QC (SI): is chrM a stable enough internal reference
# to scale (A) EU-seq (nascent) and (B) total RNA-seq (steady state)?
#
# Input (written by count_chrM_genes_hg38.sh and the two
# compute_scale_factors_mito*_hg38.sh scripts):
#   BASE/analysis/mito_qc/chrM_gene_counts.tsv
#   BASE/analysis/mito_qc/chrM_bin_counts.tsv
#   BASE/analysis/mito_qc/chrM_genes.bed
#   BASE/scale_factors_mito_samples.txt          (EU-seq)
#   BASE/scale_factors_mito.txt
#   BASE/scale_factors_mito_total_samples.txt    (total RNA-seq)
#   BASE/scale_factors_mito_total.txt
#
# Figure 1 (main SI figure), three separate figures
#   1a  A  chrM reads per gene class (rRNA / mRNA / tRNA / unannotated)
#       B  same, as a fraction of chrM reads
#   1b  C  PCA of the chrM genes, EU-seq        (log2 RPM, genes z-scored)
#       D  PCA of the chrM genes, total RNA-seq
#   1c  E  chrM reads as % of the library, per replicate
#       F  mito scale factor: 9 pairwise replicate ratios + 6 leave-one-out
#   1d  K  which chrM genes the reads come from (% of each sample's chrM reads)
# Figure 2 (supplement), three separate figures
#   2a  G    scale factor computed from rRNA reads / mRNA reads / all chrM
#   2b  H    chrM coverage (RPM per 100 bp) with the gene track
#   2c  I/J  PCA of the within-chrM composition (library depth removed)
#
# NoEU (EU-seq background control) is shown in A/B only: it is not
# part of any factor and has no library size in the factor files.
#
# Plots are print()ed. The RStudio plot pane is small and squeezes the
# layout: set SAVE_FIGS <- TRUE to write each figure at its real size
# (FIG_FMT "svg" or "pdf") into BASE/analysis/mito_qc and view those files.
# Sizes (in): 1a 11 x 5, 1b 11 x 4.5, 1c 11 x 5.5, 1d 11 x 9, 2a 8 x 5.5, 2b 11 x 7, 2c 11 x 5.
##################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})
select <- dplyr::select   ## DESeq2 / GenomicFeatures / AnnotationDbi loaded earlier in the session mask dplyr::select

BASE     <- Sys.getenv("MITO_QC_BASE", unset = "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq")
QC_DIR   <- file.path(BASE, "analysis", "mito_qc")
SAVE_FIGS     <- FALSE
FIG_FMT       <- "svg"    # "svg" or "pdf" (used when SAVE_FIGS is TRUE)
MIN_GENE_READS <- 20      # a chrM gene needs >= this many reads over the 6 samples to enter the PCA
BIN           <- 100      # must match BIN in count_chrM_genes_hg38.sh
CHRM_LEN      <- 16569

ASSAYS    <- c("EU-seq", "Total RNA-seq")
COND_LEV  <- c("WT", "POLR2A_8hrs_aux")
COND_LAB  <- c(WT = "Wildtype cells", POLR2A_8hrs_aux = "POLR2A-degraded cells")
COND_COL  <- c(WT = "#3B6FB6", POLR2A_8hrs_aux = "#D1495B", NoEU = "#999999")
CLASS_LEV <- c("rRNA", "mRNA", "tRNA", "Other", "Unannotated")
CLASS_LAB <- c("rRNA (12S, 16S)", "mRNA (protein-coding)", "tRNA", "Other", "Unannotated")
CLASS_COL <- c("rRNA (12S, 16S)" = "#E69F00", "mRNA (protein-coding)" = "#0072B2",
               "tRNA" = "#009E73", "Other" = "#CC79A7", "Unannotated" = "#999999")
SAMPLE_LEV <- c(paste0("WT_", 1:3), paste0("POLR2A_", 1:3), "NoEU")

theme_qc <- function() {
  theme_classic(base_size = 10) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          legend.key.size = unit(0.35, "cm"),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

##--------------------------------------------------------------## 1. Read + check
rd <- function(f) {
  if (!file.exists(f)) stop("missing input: ", f, "\n  (run count_chrM_genes_hg38.sh / the compute_scale_factors_mito*_hg38.sh scripts first)")
  read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
}

gc <- rd(file.path(QC_DIR, "chrM_gene_counts.tsv"))
bn <- rd(file.path(QC_DIR, "chrM_bin_counts.tsv"))
gb <- read.delim(file.path(QC_DIR, "chrM_genes.bed"), header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("chr", "start", "end", "name", "score", "strand"))

samples <- bind_rows(
  rd(file.path(BASE, "scale_factors_mito_samples.txt")) %>%
    transmute(assay = "EU-seq", sample, condition, library = read1_library, chrM_reads_sf = chrM_reads),
  rd(file.path(BASE, "scale_factors_mito_total_samples.txt")) %>%
    transmute(assay = "Total RNA-seq", sample, condition, library, chrM_reads_sf = chrM_reads)
) %>% mutate(chrM_rpm = chrM_reads_sf / library * 1e6)

sf_file <- c("EU-seq" = rd(file.path(BASE, "scale_factors_mito.txt")) %>%
               filter(condition == "POLR2A_8hrs_aux") %>% pull(eu_sf),
             "Total RNA-seq" = rd(file.path(BASE, "scale_factors_mito_total.txt")) %>%
               filter(condition == "POLR2A_8hrs_aux") %>% pull(mito_sf))

gc <- gc %>%
  mutate(class_f = factor(class, levels = CLASS_LEV, labels = CLASS_LAB),
         assay = factor(assay, levels = ASSAYS))

## the gene-class counts must add up to the chrM reads the factors were computed from
chk <- gc %>% group_by(assay, sample) %>% summarise(n = sum(reads), .groups = "drop") %>%
  inner_join(samples %>% mutate(assay = factor(assay, levels = ASSAYS)), by = c("assay", "sample"))
stopifnot("gene-class counts do not add up to chrM_reads in the scale-factor files" = all(chk$n == chk$chrM_reads_sf))
stopifnot("expected 12 samples with factors (6 EU-seq + 6 total)" = nrow(samples) == 12)
message("Inputs OK: ", nrow(chk), " samples; gene-class totals == chrM_reads of the scale-factor files")

lab_sample <- function(x) gsub("_", " ", x)
cond_f <- function(x) factor(x, levels = COND_LEV, labels = COND_LAB[COND_LEV])
samples <- samples %>% mutate(assay = factor(assay, levels = ASSAYS))

##--------------------------------------------------------------## 2. Scale-factor statistics
pair_ratios <- function(df) {  # df: sample, condition, rpm
  w <- df %>% filter(condition == "WT"); p <- df %>% filter(condition == "POLR2A_8hrs_aux")
  expand.grid(wt = w$sample, pol = p$sample, stringsAsFactors = FALSE) %>%
    mutate(ratio = w$rpm[match(wt, w$sample)] / p$rpm[match(pol, p$sample)])
}
factor_of <- function(df) mean(df$rpm[df$condition == "WT"]) / mean(df$rpm[df$condition == "POLR2A_8hrs_aux"])

stats <- list(); pairs_all <- list(); loo_all <- list()
for (a in ASSAYS) {
  d <- samples %>% filter(assay == a) %>% transmute(sample, condition, rpm = chrM_rpm)
  sf <- factor_of(d)
  if (abs(sf - sf_file[[a]]) > 1e-4 * sf) stop("recomputed factor (", sf, ") differs from the factor file (", sf_file[[a]], ") for ", a)
  pr <- pair_ratios(d) %>% mutate(assay = a)
  lo <- bind_rows(lapply(d$sample, function(s) data.frame(dropped = s, ratio = factor_of(d %>% filter(sample != s))))) %>% mutate(assay = a)
  w <- d$rpm[d$condition == "WT"]; p <- d$rpm[d$condition == "POLR2A_8hrs_aux"]
  tt <- t.test(w, p)
  stats[[a]] <- data.frame(
    assay = a, mito_sf = sf,
    wt_mean_pct = mean(w) / 1e4, wt_cv_pct = 100 * sd(w) / mean(w),
    pol_mean_pct = mean(p) / 1e4, pol_cv_pct = 100 * sd(p) / mean(p),
    fold_wt_over_pol = mean(w) / mean(p), welch_p = tt$p.value,
    pair_min = min(pr$ratio), pair_max = max(pr$ratio), pair_sd = sd(pr$ratio),
    loo_min = min(lo$ratio), loo_max = max(lo$ratio))
  pairs_all[[a]] <- pr; loo_all[[a]] <- lo
}
stats <- bind_rows(stats); pairs_all <- bind_rows(pairs_all); loo_all <- bind_rows(loo_all)
stats$assay <- factor(stats$assay, levels = ASSAYS); pairs_all$assay <- factor(pairs_all$assay, levels = ASSAYS); loo_all$assay <- factor(loo_all$assay, levels = ASSAYS)

cat("\n==== chrM as a scaling reference: summary ====\n")
print(stats %>% mutate(across(where(is.numeric), ~ signif(., 3))), row.names = FALSE)

cat("\n==== chrM reads per class (all samples) ====\n")
print(gc %>% group_by(assay, sample, class_f) %>% summarise(reads = sum(reads), .groups = "drop") %>%
        pivot_wider(names_from = class_f, values_from = reads, values_fill = 0) %>% as.data.frame(), row.names = FALSE)

##--------------------------------------------------------------## 3. Panels A-B: reads per class
cls <- gc %>% group_by(assay, sample, condition, class_f) %>% summarise(reads = sum(reads), .groups = "drop") %>%
  mutate(sample_f = factor(sample, levels = SAMPLE_LEV), cond_col = condition)
cls <- cls %>% filter(class_f %in% unique(class_f[reads > 0]))

pA <- ggplot(cls, aes(sample_f, reads, fill = class_f)) +
  geom_col(width = 0.75) +
  facet_wrap(~ assay, scales = "free", nrow = 1) +
  scale_x_discrete(labels = lab_sample) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = CLASS_COL) +
  labs(title = "A  chrM reads per gene class", x = NULL, y = "chrM reads") +
  theme_qc()

pB <- ggplot(cls, aes(sample_f, reads, fill = class_f)) +
  geom_col(width = 0.75, position = "fill") +
  facet_wrap(~ assay, scales = "free_x", nrow = 1) +
  scale_x_discrete(labels = lab_sample) +
  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = CLASS_COL) +
  labs(title = "B  Composition of chrM reads", x = NULL, y = "Fraction of chrM reads") +
  theme_qc()

##--------------------------------------------------------------## 3b. Panel K: chrM reads per gene
gene_order <- vapply(strsplit(gb$name[order(gb$start)], "|", fixed = TRUE), `[`, "", 2)
gk <- gc %>%
  group_by(assay, sample) %>% mutate(pct = 100 * reads / sum(reads)) %>% ungroup() %>%
  mutate(gene_lab = ifelse(class == "Unannotated", "Unannotated", sub("^MT-", "", gene_name)),
         gene_f = factor(gene_lab, levels = rev(unique(c(sub("^MT-", "", gene_order), "Unannotated")))),
         sample_f = factor(sample, levels = SAMPLE_LEV))
stopifnot("a chrM gene is missing from the gene order" = !any(is.na(gk$gene_f)))

pK <- ggplot(gk, aes(sample_f, gene_f, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(data = gk %>% filter(pct >= 0.5), aes(label = sprintf("%.1f", pct), colour = pct > 12), size = 2.3, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey20")) +
  facet_grid(class_f ~ assay, scales = "free", space = "free") +
  scale_x_discrete(labels = lab_sample) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B", trans = "sqrt", name = "% of chrM reads") +
  labs(title = "K  Which chrM genes the reads come from", x = NULL, y = NULL,
       subtitle = "% of each sample's chrM reads; rows in chrM order within each gene class; blank = under 0.5%") +
  theme_qc() + theme(legend.title = element_text(size = 8), strip.text.y = element_text(angle = 0, hjust = 0),
                     panel.spacing = unit(0.3, "lines"))

##--------------------------------------------------------------## 4. PCA of the chrM genes
rpm_genes <- gc %>%
  filter(class != "Unannotated", sample != "NoEU") %>%
  inner_join(samples %>% mutate(assay = factor(assay, levels = ASSAYS)) %>% select(assay, sample, library), by = c("assay", "sample")) %>%
  mutate(rpm = reads / library * 1e6)

run_pca <- function(a, value = c("rpm", "composition")) {
  value <- match.arg(value)
  d <- rpm_genes %>% filter(assay == a)
  if (value == "composition") {
    d <- d %>% group_by(sample) %>% mutate(v = log2(reads / sum(reads) * 1e6 + 1)) %>% ungroup()
  } else {
    d <- d %>% mutate(v = log2(rpm + 1))
  }
  keep <- d %>% group_by(gene_name) %>% summarise(n = sum(reads), .groups = "drop") %>% filter(n >= MIN_GENE_READS) %>% pull(gene_name)
  m <- d %>% filter(gene_name %in% keep) %>% select(sample, gene_name, v) %>%
    pivot_wider(names_from = gene_name, values_from = v) %>% as.data.frame()
  rownames(m) <- m$sample; m$sample <- NULL; m <- as.matrix(m)
  m <- m[, apply(m, 2, sd) > 0, drop = FALSE]
  pc <- prcomp(m, center = TRUE, scale. = TRUE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  sc <- as.data.frame(pc$x[, 1:2, drop = FALSE]); sc$sample <- rownames(sc)
  sc$condition <- samples$condition[match(paste(a, sc$sample), paste(samples$assay, samples$sample))]
  list(assay = a, value = value, scores = sc, ve = ve, rot = pc$rotation, n_genes = ncol(m))
}

plot_pca <- function(pc, tag, title) {
  ggplot(pc$scores, aes(PC1, PC2, colour = condition, label = lab_sample(sample))) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(size = 3) +
    geom_text_repel(size = 2.8, min.segment.length = Inf, seed = 1, show.legend = FALSE) +
    scale_colour_manual(values = COND_COL, breaks = COND_LEV, labels = COND_LAB[COND_LEV]) +
    scale_x_continuous(expand = expansion(mult = 0.2)) + scale_y_continuous(expand = expansion(mult = 0.2)) +
    labs(title = paste0(tag, "  ", title), subtitle = paste0(pc$n_genes, " chrM genes; n = 3 per condition"),
         x = sprintf("PC1 (%.0f%%)", pc$ve[1]), y = sprintf("PC2 (%.0f%%)", pc$ve[2])) +
    theme_qc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

pca <- list(eu = run_pca("EU-seq", "rpm"), tot = run_pca("Total RNA-seq", "rpm"),
            eu_c = run_pca("EU-seq", "composition"), tot_c = run_pca("Total RNA-seq", "composition"))

cat("\n==== PCA of chrM genes (log2 RPM, genes z-scored): scores and variance explained ====\n")
for (k in c("eu", "tot")) {
  cat("\n", pca[[k]]$assay, " -- PC1 ", sprintf("%.1f", pca[[k]]$ve[1]), "%, PC2 ", sprintf("%.1f", pca[[k]]$ve[2]), "%, ", pca[[k]]$n_genes, " genes\n", sep = "")
  print(pca[[k]]$scores %>% mutate(PC1 = round(PC1, 2), PC2 = round(PC2, 2)) %>% select(sample, condition, PC1, PC2), row.names = FALSE)
  top <- sort(abs(pca[[k]]$rot[, 1]), decreasing = TRUE)[1:min(5, nrow(pca[[k]]$rot))]
  cat("  top PC1 loadings: ", paste(names(top), collapse = ", "), "\n", sep = "")
}

pC <- plot_pca(pca$eu,  "C", "PCA of chrM genes: EU-seq")
pD <- plot_pca(pca$tot, "D", "PCA of chrM genes: total RNA-seq")
pI <- plot_pca(pca$eu_c,  "I", "chrM composition: EU-seq")
pJ <- plot_pca(pca$tot_c, "J", "chrM composition: total RNA-seq")

##--------------------------------------------------------------## 5. Panel E: chrM % of library
sm <- samples %>% mutate(pct = chrM_rpm / 1e4, cond_f = cond_f(condition))
cond_stats <- sm %>% group_by(assay, cond_f) %>% summarise(m = mean(pct), s = sd(pct), .groups = "drop")
ann <- stats %>% transmute(assay, lab = sprintf("factor %.2f\nWelch p = %.2g\nCV %.0f%% / %.0f%%", mito_sf, welch_p, wt_cv_pct, pol_cv_pct))

pE <- ggplot() +
  geom_col(data = cond_stats, aes(cond_f, m, fill = cond_f), width = 0.6, alpha = 0.35, show.legend = FALSE) +
  geom_errorbar(data = cond_stats, aes(cond_f, ymin = m - s, ymax = m + s), width = 0.15) +
  geom_point(data = sm, aes(cond_f, pct, colour = cond_f), size = 2.6, position = position_jitter(width = 0.08, height = 0, seed = 1), show.legend = FALSE) +
  geom_text(data = ann, aes(x = 1.5, y = Inf, label = lab), vjust = 1.15, size = 2.7, lineheight = 0.95) +
  facet_wrap(~ assay, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = setNames(COND_COL[COND_LEV], COND_LAB[COND_LEV])) +
  scale_colour_manual(values = setNames(COND_COL[COND_LEV], COND_LAB[COND_LEV])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "E  chrM share of the library", x = NULL, y = "chrM reads (% of library)",
       subtitle = "bars = mean +/- SD; CV given as WT / POLR2A-degraded") +
  theme_qc()

##--------------------------------------------------------------## 6. Panel F: factor spread
fdat <- bind_rows(pairs_all %>% transmute(assay, type = "Pairwise replicate ratios (9)", ratio),
                  loo_all %>% transmute(assay, type = "Leave-one-sample-out (6)", ratio))
fdat$type <- factor(fdat$type, levels = c("Pairwise replicate ratios (9)", "Leave-one-sample-out (6)"))

pF <- ggplot(fdat, aes(assay, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(aes(shape = type, colour = type), size = 2.3, alpha = 0.85,
             position = position_jitter(width = 0.12, height = 0, seed = 2)) +
  geom_errorbar(data = stats, aes(assay, y = mito_sf, ymin = mito_sf, ymax = mito_sf), width = 0.5, colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
  scale_y_continuous(trans = "log2", breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2)) +
  scale_shape_manual(values = c(16, 17)) +
  scale_colour_manual(values = c("grey35", "#7B4EA3")) +
  labs(title = "F  Stability of the mito scale factor", x = NULL, y = "POLR2A-degraded scale factor\n(WT / POLR2A chrM RPM)",
       subtitle = "black bar = factor used (ratio of means); dashed = 1 (no scaling)") +
  theme_qc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.direction = "vertical")

##--------------------------------------------------------------## 7. Panel G: class-specific factors
cls_rpm <- gc %>% filter(class %in% c("rRNA", "mRNA"), sample != "NoEU") %>%
  group_by(assay, sample, condition, class) %>% summarise(reads = sum(reads), .groups = "drop") %>%
  inner_join(samples %>% mutate(assay = factor(assay, levels = ASSAYS)) %>% select(assay, sample, library), by = c("assay", "sample")) %>%
  mutate(rpm = reads / library * 1e6, grp = recode(class, rRNA = "rRNA reads only", mRNA = "mRNA reads only"))
all_rpm <- samples %>% mutate(assay = factor(assay, levels = ASSAYS)) %>% transmute(assay, sample, condition, grp = "All chrM reads", rpm = chrM_rpm)
g_dat <- bind_rows(cls_rpm %>% select(assay, sample, condition, grp, rpm), all_rpm)

g_pairs <- g_dat %>% group_by(assay, grp) %>% group_modify(~ pair_ratios(.x)) %>% ungroup()
g_sf <- g_dat %>% group_by(assay, grp) %>% summarise(sf = factor_of(data.frame(condition = condition, rpm = rpm)), .groups = "drop")
g_pairs$grp <- factor(g_pairs$grp, levels = c("All chrM reads", "rRNA reads only", "mRNA reads only"))
g_sf$grp    <- factor(g_sf$grp,    levels = levels(g_pairs$grp))

cat("\n==== Mito factor computed from different chrM read classes ====\n")
print(g_sf %>% mutate(sf = round(sf, 3)) %>% arrange(assay, grp) %>% as.data.frame(), row.names = FALSE)

pG <- ggplot(g_pairs, aes(grp, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(colour = "grey35", alpha = 0.8, size = 2, position = position_jitter(width = 0.1, height = 0, seed = 3)) +
  geom_errorbar(data = g_sf, aes(grp, y = sf, ymin = sf, ymax = sf), width = 0.5, linewidth = 0.8, inherit.aes = FALSE) +
  facet_wrap(~ assay, nrow = 1) +
  scale_y_continuous(trans = "log2", breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 4)) +
  labs(title = "G  Does the factor depend on which chrM reads are used?", x = NULL, y = "WT / POLR2A scale factor",
       subtitle = "points = 9 pairwise replicate ratios; bar = ratio of means") +
  theme_qc()

##--------------------------------------------------------------## 8. Panel H: chrM coverage
bins <- seq(0, CHRM_LEN - 1, by = BIN)
cov <- bn %>% filter(sample != "NoEU") %>%
  mutate(assay = factor(assay, levels = ASSAYS)) %>%
  inner_join(samples %>% mutate(assay = factor(assay, levels = ASSAYS)) %>% select(assay, sample, library), by = c("assay", "sample")) %>%
  mutate(rpm = reads / library * 1e6) %>%
  select(assay, sample, condition, bin_start, rpm)
cov <- cov %>% group_by(assay, sample, condition) %>% complete(bin_start = bins, fill = list(rpm = 0)) %>% ungroup() %>%
  mutate(x = bin_start + BIN / 2, line_id = paste(assay, sample))

cov_plot <- function(a, title) {
  ggplot(cov %>% filter(assay == a), aes(x, rpm, colour = condition, group = line_id)) +
    geom_line(linewidth = 0.35, alpha = 0.8) +
    scale_colour_manual(values = COND_COL, breaks = COND_LEV, labels = COND_LAB[COND_LEV]) +
    scale_x_continuous(limits = c(0, CHRM_LEN), expand = c(0, 0)) +
    scale_y_continuous(trans = "sqrt", labels = label_comma(), expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, x = NULL, y = paste0("RPM per ", BIN, " bp")) +
    theme_qc() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
}

gb <- gb %>% mutate(parts = strsplit(name, "|", fixed = TRUE),
                    gene_name = vapply(parts, `[`, "", 2), biotype = vapply(parts, `[`, "", 3),
                    class_f = factor(case_when(biotype == "Mt_rRNA" ~ "rRNA (12S, 16S)", biotype == "protein_coding" ~ "mRNA (protein-coding)",
                                               biotype == "Mt_tRNA" ~ "tRNA", TRUE ~ "Other"), levels = CLASS_LAB),
                    ypos = ifelse(strand == "+", 1, -1))
pTrack <- ggplot(gb) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_rect(aes(xmin = start, xmax = end, ymin = ifelse(strand == "+", 0.05, -1), ymax = ifelse(strand == "+", 1, -0.05), fill = class_f)) +
  geom_text(data = gb %>% filter(biotype != "Mt_tRNA"), aes(x = (start + end) / 2, y = ypos * 0.5, label = sub("^MT-", "", gene_name)),
            size = 1.9, check_overlap = TRUE, colour = "white") +
  scale_fill_manual(values = CLASS_COL) +
  scale_x_continuous(limits = c(0, CHRM_LEN), expand = c(0, 0), labels = label_comma()) +
  scale_y_continuous(breaks = c(-0.5, 0.5), labels = c("- strand", "+ strand")) +
  labs(x = "chrM position (bp)", y = NULL) +
  theme_qc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

pH <- (cov_plot("EU-seq", "H  chrM coverage: EU-seq") / cov_plot("Total RNA-seq", "chrM coverage: total RNA-seq") / pTrack) +
  plot_layout(heights = c(4, 4, 1.2), guides = "collect")

##--------------------------------------------------------------## 9. Assemble + print
fig1a <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "chrM as an internal reference for EU-seq and total RNA-seq scaling",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))
fig1b <- pC | pD
fig1c <- (pE | pF) + plot_layout(guides = "keep")
fig1d <- pK
fig2a <- pG
fig2b <- pH
fig2c <- pI | pJ

write.table(stats, file.path(QC_DIR, "mito_qc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWritten: ", file.path(QC_DIR, "mito_qc_summary.tsv"), "\n", sep = "")

print(fig1a)
print(fig1b)
print(fig1c)
print(fig1d)
print(fig2a)
print(fig2b)
print(fig2c)

if (SAVE_FIGS) {
  if (FIG_FMT == "svg" && !isTRUE(unname(capabilities("cairo")))) {
    message("This R has no cairo: svg() is unavailable, writing PDFs instead.")
    FIG_FMT <- "pdf"
  }
  save_fig <- function(p, name, w, h) {
    f <- file.path(QC_DIR, paste0(name, ".", FIG_FMT))
    if (FIG_FMT == "svg") svg(f, width = w, height = h) else pdf(f, width = w, height = h)
    print(p)
    dev.off()
    cat("Written: ", f, "\n", sep = "")
  }
  save_fig(fig1a, "mito_scaling_qc_fig1a", 11, 5)
  save_fig(fig1b, "mito_scaling_qc_fig1b", 11, 4.5)
  save_fig(fig1c, "mito_scaling_qc_fig1c", 11, 5.5)
  save_fig(fig1d, "mito_scaling_qc_fig1d", 11, 9)
  save_fig(fig2a, "mito_scaling_qc_fig2a", 8, 5.5)
  save_fig(fig2b, "mito_scaling_qc_fig2b", 11, 7)
  save_fig(fig2c, "mito_scaling_qc_fig2c", 11, 5)
}