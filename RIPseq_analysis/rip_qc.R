##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# RIP-seq QC (R companion to ripseq_qc.sh): replicate correlation,
# IP enrichment over Input, and genomic feature distribution of
# signal. Standalone script — reloads bedgraphs and feature BEDs
# fresh from disk rather than assuming a shared session with
# ripseq_feature_beds.R.
#
# Peak-calling-based FRiP is deliberately excluded — see
# ripseq_qc.sh's docstring for why (already tried via Piranha/
# MACS2, found uninformative for this specific IP-vs-multi-
# polymerase-Input comparison). The genomic feature distribution
# below is the substitute used instead.
##################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(GenomicRanges)
  library(plyranges)
  library(rtracklayer)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

##--------------------------------------------------------------## 1. Parameters

BGR_DIR          <- "/projects/b1042/BackmanLab/Lucas/022025_RIPseq/ripseq_pipeline/bedgraphs"
CHROM_SIZES_PATH <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/chroms/hg38.chrom.sizes"
FEATURE_BED_DIR  <- "/projects/b1042/BackmanLab/Lucas/022025_RIPseq/feature_beds"
EUSEQ_BED_DIR    <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline/intron_dog/beds"
GTF_PATH         <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
OUT_DIR          <- "/projects/b1042/BackmanLab/Lucas/022025_RIPseq/ripseq_pipeline/qc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES  <- c("P2AB_Rep2", "P2AB_Rep3", "WT_Input", "EZH2_Rep1", "EZH2_Rep2")
BIN_SIZE <- 1000L   ## genome-wide bin size for correlation/enrichment

## NOT yet validated against this dataset's actual coverage scale —
## if the enrichment distribution looks compressed/uninformative once
## plotted, this is the first parameter to reconsider (bin-level
## bp-summed coverage may sit at a different numeric scale than the
## RPKM-based pseudocount=1 convention used elsewhere in this pipeline).
PSEUDOCOUNT_QC <- 1

base_theme <- theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey40", size = 11))

##--------------------------------------------------------------## 2. Load bedgraphs, build genome-wide bins

read_bedgraph <- function(path) {
  df <- read_tsv(path, col_names = c("chr","start","end","score"),
                 col_types = "ciid", show_col_types = FALSE)
  GRanges(seqnames = df$chr, ranges = IRanges(df$start + 1L, df$end), score = df$score)
}

message("Loading bedgraphs...")
bedgraphs <- setNames(lapply(SAMPLES, function(s) {
  list(plus  = read_bedgraph(file.path(BGR_DIR, paste0(s, "_plus.norm.bgr"))),
       minus = read_bedgraph(file.path(BGR_DIR, paste0(s, "_minus.norm.bgr"))))
}), SAMPLES)
bedgraphs_unstranded <- setNames(lapply(bedgraphs, function(b) c(b$plus, b$minus)), SAMPLES)

message("Loading chrom sizes and building genome-wide bins (", BIN_SIZE, "bp)...")
chrom_sizes <- read_tsv(CHROM_SIZES_PATH, col_names = c("chr","size"), col_types = "ci", show_col_types = FALSE)

## Match chrom_sizes naming to the bedgraph's naming style (recurring
## fix throughout this whole project)
has_chr_bgr <- any(grepl("^chr", as.character(seqnames(bedgraphs_unstranded[[1]]))))
has_chr_cs  <- any(grepl("^chr", chrom_sizes$chr))
if (has_chr_bgr && !has_chr_cs) chrom_sizes$chr <- paste0("chr", chrom_sizes$chr)
if (!has_chr_bgr && has_chr_cs) chrom_sizes$chr <- sub("^chr", "", chrom_sizes$chr)

seqinfo_obj <- Seqinfo(seqnames = chrom_sizes$chr, seqlengths = chrom_sizes$size)
genome_bins <- GenomicRanges::tileGenome(seqinfo_obj, tilewidth = BIN_SIZE, cut.last.tile.in.chrom = TRUE)
mcols(genome_bins)$bin_id <- paste0("bin_", seq_along(genome_bins))
cat("Genome bins:", length(genome_bins), "\n")

##--------------------------------------------------------------## 3. Per-bin coverage per sample (bp-weighted, matches established convention)

compute_bin_coverage <- function(bedgraph_gr) {
  overlaps <- plyranges::join_overlap_intersect(genome_bins, bedgraph_gr)
  if (length(overlaps) == 0) {
    return(tibble(bin_id = mcols(genome_bins)$bin_id, cov = 0))
  }
  cov_tbl <- tibble(bin_id = mcols(overlaps)$bin_id,
                    cov_bp = mcols(overlaps)$score * width(overlaps)) %>%
    dplyr::group_by(bin_id) %>%
    dplyr::summarise(cov_bp = sum(cov_bp), .groups = "drop")
  tibble(bin_id = mcols(genome_bins)$bin_id) %>%
    dplyr::left_join(cov_tbl, by = "bin_id") %>%
    dplyr::mutate(cov = tidyr::replace_na(cov_bp, 0)) %>%
    dplyr::select(bin_id, cov)
}

message("Computing per-bin coverage for each sample...")
bin_cov <- setNames(lapply(SAMPLES, function(s) compute_bin_coverage(bedgraphs_unstranded[[s]])), SAMPLES)

bin_matrix <- bin_cov$P2AB_Rep2 %>% dplyr::rename(P2AB_Rep2 = cov) %>%
  dplyr::left_join(bin_cov$P2AB_Rep3 %>% dplyr::rename(P2AB_Rep3 = cov), by = "bin_id") %>%
  dplyr::left_join(bin_cov$WT_Input  %>% dplyr::rename(WT_Input  = cov), by = "bin_id")

##--------------------------------------------------------------## 4. Replicate correlation (P2AB_Rep2 vs Rep3)

message("Computing replicate correlation...")
## log1p to de-emphasize the huge mass of near-zero bins — standard
## for genome-wide coverage correlation, where most bins have little
## to no signal regardless of library quality.
bin_matrix <- bin_matrix %>%
  dplyr::mutate(log_rep2 = log1p(P2AB_Rep2), log_rep3 = log1p(P2AB_Rep3))

pearson_r  <- cor(bin_matrix$log_rep2, bin_matrix$log_rep3, method = "pearson",  use = "complete.obs")
spearman_r <- cor(bin_matrix$log_rep2, bin_matrix$log_rep3, method = "spearman", use = "complete.obs")
cat("Replicate correlation (P2AB_Rep2 vs Rep3):\n")
cat("  Pearson r  (log1p coverage):", round(pearson_r, 3), "\n")
cat("  Spearman r (log1p coverage):", round(spearman_r, 3), "\n")

## Subsampled for the scatter — plotting all ~3M genome bins directly
## would be a massive overplotted mess and slow to render; correlation
## coefficients above are computed on the FULL data, not the subsample.
p_replicate_corr <- bin_matrix %>%
  dplyr::filter(P2AB_Rep2 > 0 | P2AB_Rep3 > 0) %>%
  dplyr::slice_sample(n = min(200000, nrow(.))) %>%
  ggplot(aes(x = log_rep2, y = log_rep3)) +
  geom_point(size = 0.3, alpha = 0.15, color = "#CC6677") +
  geom_smooth(method = "lm", se = FALSE, color = "grey20", linewidth = 0.6) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = paste0("Pearson r = ", round(pearson_r, 3), "\nSpearman r = ", round(spearman_r, 3)),
           size = 4) +
  labs(x = "log1p(coverage) \u2014 P2AB Rep2", y = "log1p(coverage) \u2014 P2AB Rep3",
       title = "P2AB replicate correlation",
       subtitle = paste0(BIN_SIZE, "bp genome-wide bins  |  ", scales::comma(nrow(bin_matrix)),
                         " bins total, plotted subsample of nonzero bins")) +
  base_theme

##--------------------------------------------------------------## 5. IP enrichment over Input

message("Computing IP enrichment over Input...")
bin_matrix <- bin_matrix %>%
  dplyr::mutate(
    P2AB_mean = rowMeans(cbind(P2AB_Rep2, P2AB_Rep3), na.rm = TRUE),
    log2fc = log2((P2AB_mean + PSEUDOCOUNT_QC) / (WT_Input + PSEUDOCOUNT_QC))
  )

cat("\nIP enrichment over Input (log2 P2AB/Input), all genome bins:\n")
print(summary(bin_matrix$log2fc))

p_enrichment <- ggplot(bin_matrix, aes(x = log2fc)) +
  geom_histogram(bins = 100, fill = "#882255", color = NA, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_log10(labels = scales::comma) +
  labs(x = "log2(P2AB / WT_Input)", y = "Number of genome bins (log10 scale)",
       title = "IP enrichment over Input, genome-wide",
       subtitle = paste0(BIN_SIZE, "bp bins  |  y-axis log-scaled: most bins have ~zero signal",
                         " in both samples (pseudocount collapses these to log2FC~0),",
                         " which would otherwise hide the real enrichment tail")) +
  base_theme

##--------------------------------------------------------------## 6. Genomic feature distribution of signal
## Fraction of total P2AB / Input signal falling in each feature
## category vs everything else — the RIP-seq-appropriate substitute
## for peak-calling-based FRiP here (see ripseq_qc.sh docstring).
##
## NOTE: this reloads exon/intron BEDs fresh from EUSEQ_BED_DIR —
## these are the ORIGINAL, un-subtracted versions, NOT the enhancer-
## subtracted intron intervals used elsewhere in ripseq_feature_beds.R.
## Fine for a QC overview, but the Intron category here won't reflect
## that correction.
##
## These five categories can OVERLAP each other (e.g. an enhancer
## sitting inside an intron) — bars below won't necessarily sum to
## exactly 100%, shown as-is rather than forced to.

message("Loading feature BEDs for genomic distribution...")
read_bed_simple <- function(path, col_names) {
  read_tsv(path, col_names = col_names, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
}

promoter_bed <- read_bed_simple(file.path(FEATURE_BED_DIR, "promoter_regions.bed"), c("chr","start","end","gene_id"))
enhancer_bed <- read_bed_simple(file.path(FEATURE_BED_DIR, "enhancer_regions.bed"), c("chr","start","end","gene_id"))
tes_bed      <- read_bed_simple(file.path(FEATURE_BED_DIR, "tes_windows.bed"),      c("chr","start","end","gene_id"))
exon_plus_bed    <- read_bed_simple(file.path(EUSEQ_BED_DIR, "exons_plus.bed"),         c("chr","start","end","gene_id","score","strand"))
exon_minus_bed   <- read_bed_simple(file.path(EUSEQ_BED_DIR, "exons_minus.bed"),        c("chr","start","end","gene_id","score","strand"))
intron_plus_bed  <- read_bed_simple(file.path(EUSEQ_BED_DIR, "hg38_introns_plus.bed"),  c("chr","start","end","gene_id","score","strand"))
intron_minus_bed <- read_bed_simple(file.path(EUSEQ_BED_DIR, "hg38_introns_minus.bed"), c("chr","start","end","gene_id","score","strand"))

to_gr <- function(df) GRanges(seqnames = df$chr, ranges = IRanges(as.integer(df$start) + 1L, as.integer(df$end)))

feature_grs <- list(
  Promoter = to_gr(promoter_bed),
  Enhancer = to_gr(enhancer_bed),
  Exon     = c(to_gr(exon_plus_bed), to_gr(exon_minus_bed)),
  Intron   = c(to_gr(intron_plus_bed), to_gr(intron_minus_bed)),
  TES      = to_gr(tes_bed)
)

fix_seqlevel_style_qc <- function(gr, ref_gr) {
  has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
  has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
  if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
  else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
  suppressWarnings(gr)
}
feature_grs <- lapply(feature_grs, fix_seqlevel_style_qc, ref_gr = bedgraphs_unstranded[[1]])

compute_feature_signal_fraction <- function(bedgraph_gr, feature_grs) {
  total_signal <- sum(as.numeric(mcols(bedgraph_gr)$score) * width(bedgraph_gr))
  feature_signal <- sapply(names(feature_grs), function(fn) {
    ov <- plyranges::join_overlap_intersect(feature_grs[[fn]], bedgraph_gr)
    if (length(ov) == 0) return(0)
    sum(mcols(ov)$score * width(ov))
  })
  tibble(feature = names(feature_signal), signal = as.numeric(feature_signal)) %>%
    dplyr::mutate(fraction = signal / total_signal)
}

message("Computing genomic distribution of signal per sample...")
dist_p2ab_rep2 <- compute_feature_signal_fraction(bedgraphs_unstranded$P2AB_Rep2, feature_grs) %>% dplyr::mutate(sample = "P2AB_Rep2")
dist_p2ab_rep3 <- compute_feature_signal_fraction(bedgraphs_unstranded$P2AB_Rep3, feature_grs) %>% dplyr::mutate(sample = "P2AB_Rep3")
dist_input     <- compute_feature_signal_fraction(bedgraphs_unstranded$WT_Input,  feature_grs) %>% dplyr::mutate(sample = "WT_Input")

dist_all <- dplyr::bind_rows(dist_p2ab_rep2, dist_p2ab_rep3, dist_input) %>%
  dplyr::mutate(feature = factor(feature, levels = c("Promoter","Enhancer","Exon","Intron","TES")))

cat("\nFraction of total signal by feature type, per sample:\n")
print(dist_all %>% dplyr::select(sample, feature, fraction) %>%
        tidyr::pivot_wider(names_from = feature, values_from = fraction))

FEATURE_QC_COLORS <- c(Promoter = "#882255", Enhancer = "#44AA99", Exon = "#88CCEE",
                       Intron = "#CC6677", TES = "#DDCC77")

p_distribution <- ggplot(dist_all, aes(x = sample, y = fraction, fill = feature)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = FEATURE_QC_COLORS, name = "Feature") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Fraction of total signal",
       title = "Genomic distribution of signal by feature type",
       subtitle = "Categories can overlap (e.g. enhancer-in-intron) \u2014 bars won't necessarily sum to 100%") +
  base_theme

##--------------------------------------------------------------## 6b. EZH2 positive control — literature-established PRC2 target lncRNAs
## Checks whether the whole RIP -> analysis pipeline correctly detects
## a KNOWN, well-documented RNA-protein interaction, as a sanity check
## independent of any genome-wide statistic above.
##
## XIST deliberately excluded despite being the most textbook-famous
## EZH2/PRC2 target: HCT116 is male-derived (established earlier in
## this project), so XIST would not be meaningfully expressed here
## regardless of IP quality — including it risked reading as a failed
## positive control when it's really just a cell-line mismatch.

message("\nChecking EZH2 positive control (literature-established PRC2 target lncRNAs)...")

EZH2_POSITIVE_CONTROLS <- c("HOTAIR", "CDKN2B-AS1", "MEG3")

genes_gr_qc <- rtracklayer::import(GTF_PATH, feature.type = "gene")
positive_control_info <- as.data.frame(genes_gr_qc) %>%
  dplyr::filter(gene_name %in% EZH2_POSITIVE_CONTROLS) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
  dplyr::select(gene_name, gene_id, seqnames, start, end, strand)

cat("Positive control genes found in GTF:\n")
print(positive_control_info)

missing_controls <- setdiff(EZH2_POSITIVE_CONTROLS, positive_control_info$gene_name)
if (length(missing_controls) > 0) {
  cat("WARNING: not found in GTF by this exact gene_name:", paste(missing_controls, collapse = ", "), "\n")
}

if (nrow(positive_control_info) > 0) {
  positive_control_gr <- GRanges(
    seqnames  = positive_control_info$seqnames,
    ranges    = IRanges(positive_control_info$start, positive_control_info$end),
    strand    = positive_control_info$strand,
    gene_name = positive_control_info$gene_name
  )
  positive_control_gr <- fix_seqlevel_style_qc(positive_control_gr, bedgraphs_unstranded[[1]])
  
  ## Strand-specific (matching each gene's own strand) rather than the
  ## unstranded bedgraphs used for the genome-wide checks above — more
  ## precise for a small, specific gene list, avoids picking up
  ## antisense noise.
  compute_gene_coverage_stranded <- function(gr, bedgraph_plus, bedgraph_minus) {
    result <- setNames(rep(0, length(gr)), mcols(gr)$gene_name)
    gr_plus  <- gr[strand(gr) == "+"]
    gr_minus <- gr[strand(gr) == "-"]
    for (pair in list(list(gr_plus, bedgraph_plus), list(gr_minus, bedgraph_minus))) {
      sub_gr <- pair[[1]]; bgr <- pair[[2]]
      if (length(sub_gr) == 0) next
      ov <- plyranges::join_overlap_intersect(sub_gr, bgr)
      if (length(ov) == 0) next
      cov_tbl <- tibble(gene_name = mcols(ov)$gene_name, cov_bp = mcols(ov)$score * width(ov)) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::summarise(cov_bp = sum(cov_bp), .groups = "drop")
      result[cov_tbl$gene_name] <- cov_tbl$cov_bp
    }
    result
  }
  
  gene_widths <- setNames(width(positive_control_gr), mcols(positive_control_gr)$gene_name)
  
  ezh2_rep1_cov <- compute_gene_coverage_stranded(positive_control_gr, bedgraphs$EZH2_Rep1$plus, bedgraphs$EZH2_Rep1$minus)
  ezh2_rep2_cov <- compute_gene_coverage_stranded(positive_control_gr, bedgraphs$EZH2_Rep2$plus, bedgraphs$EZH2_Rep2$minus)
  input_cov     <- compute_gene_coverage_stranded(positive_control_gr, bedgraphs$WT_Input$plus,  bedgraphs$WT_Input$minus)
  ezh2_mean_cov <- (ezh2_rep1_cov + ezh2_rep2_cov) / 2
  
  positive_control_result <- tibble(
    gene_name  = names(gene_widths),
    gene_kb    = gene_widths / 1000,
    ezh2_rpkm  = ezh2_mean_cov[names(gene_widths)] / gene_kb,
    input_rpkm = input_cov[names(gene_widths)] / gene_kb
  ) %>%
    dplyr::mutate(log2fc = log2((ezh2_rpkm + 1) / (input_rpkm + 1)))
  
  cat("\nEZH2 positive control \u2014 enrichment at PRC2 target lncRNAs:\n")
  print(positive_control_result)
  
  p_ezh2_positive_control <- ggplot(positive_control_result, aes(x = gene_name, y = log2fc)) +
    geom_col(fill = "#7B2D8B", width = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = NULL, y = "log2(EZH2 / WT_Input)",
         title = "EZH2 positive control \u2014 literature-established PRC2 target lncRNAs",
         subtitle = "Positive = expected enrichment. Near-zero/negative may mean low expression\nin this cell line, not necessarily a failed IP \u2014 check input_rpkm before concluding either way.") +
    base_theme
  
  ggplot2::ggsave(file.path(OUT_DIR, "ezh2_positive_control.png"), p_ezh2_positive_control,
                  width = 6, height = 6, dpi = 300)
  cat("Saved:", file.path(OUT_DIR, "ezh2_positive_control.png"), "\n")
} else {
  message("  No positive-control genes found in GTF by gene_name \u2014 skipping this check.")
}

##--------------------------------------------------------------## 7. Combine and save

message("Saving QC plots...")
qc_panel <- (p_replicate_corr | p_enrichment) / p_distribution +
  patchwork::plot_layout(heights = c(1, 1))
qc_panel
p_replicate_corr
ggplot2::ggsave(file.path(OUT_DIR, "ripseq_qc_panel.png"), qc_panel,
                width = 14, height = 10, dpi = 300)

cat("\nQC analysis complete. Saved:", file.path(OUT_DIR, "ripseq_qc_panel.png"), "\n")

##--------------------------------------------------------------## 8. Alignment statistics (separate plot)
## Alignment RATE comes from HISAT2 log files (pre-filtering) rather
## than qc_summary.tsv's mapping_pct — that column was computed from
## already-filtered BAMs (every filename in this project has "_filt"
## in it), where mapped reads and total reads are the same set by
## construction, trivially producing ~100% for every sample. HISAT2's
## own "overall alignment rate" line reflects true aligner performance
## on the original reads, before any filtering happened.
##
## Sequencing depth (below) still comes from qc_summary.tsv — total
## read COUNT isn't affected by the filtering issue the way a RATIO
## is, but note "mapped_reads" there has the same post-filter caveat
## as the old mapping_pct did; not rebuilt here since only the rate
## plot was asked for.

message("\nBuilding alignment statistics plot...")

HISAT2_LOG_DIR <- "/projects/b1042/BackmanLab/Lucas/022025_RIPseq/SAM"

## Filenames verified against the actual SAM directory listing.
HISAT2_LOG_FILES <- c(
  P2AB_Rep2  = "WT_P2AB_Rep2_S19_L001_hisat2.log",
  P2AB_Rep3  = "WT_P2AB_Rep3_S20_L001_hisat2.log",
  EZH2_Rep1  = "WT_EZH2_Rep1_S9_L001_hisat2.log",
  EZH2_Rep2  = "WT_EZH2_Rep2_S10_L001_hisat2.log",
  WT_Input   = "WT_Input_S13_L001_hisat2.log",
  P2KO_input = "P2KO_input_S12_L001_hisat2.log",
  WT_IgG     = "WT_IgG_S14_L001_hisat2.log"
)

parse_hisat2_overall_rate <- function(log_path) {
  if (!file.exists(log_path)) return(NA_real_)
  lines <- readLines(log_path)
  rate_line <- lines[grepl("overall alignment rate", lines)]
  if (length(rate_line) == 0) return(NA_real_)
  as.numeric(sub("%.*", "", trimws(rate_line[1])))
}

hisat2_rates <- tibble(
  sample = names(HISAT2_LOG_FILES),
  overall_alignment_rate = sapply(HISAT2_LOG_FILES, function(f) {
    parse_hisat2_overall_rate(file.path(HISAT2_LOG_DIR, f))
  })
) %>%
  dplyr::mutate(sample = factor(sample, levels = sample))

missing_logs <- hisat2_rates %>% dplyr::filter(is.na(overall_alignment_rate))
if (nrow(missing_logs) > 0) {
  cat("WARNING: could not find/parse HISAT2 log for:",
      paste(missing_logs$sample, collapse = ", "), "\n")
}

## Single- vs multi-mapper breakdown, parsed directly from the same
## HISAT2 logs — "aligned concordantly exactly 1 time" (single mapper)
## and "aligned concordantly >1 times" (multi-mapper). Stack height is
## SLIGHTLY LOWER than hisat2_rates$overall_alignment_rate above: HISAT2
## also recovers a small number of reads via discordant/single-mate
## alignment beyond these two concordant categories, which aren't
## counted here — this is the concordant-alignment breakdown
## specifically, not identical to the overall rate.
parse_hisat2_mapper_breakdown <- function(log_path) {
  if (!file.exists(log_path)) {
    return(tibble(category = c("Single mapper", "Multi-mapper"), pct = NA_real_))
  }
  lines <- trimws(readLines(log_path))
  
  extract_pct <- function(pattern) {
    hit <- lines[grepl(pattern, lines)]
    if (length(hit) == 0) return(NA_real_)
    m <- regmatches(hit[1], regexec("\\(([0-9.]+)%\\)", hit[1]))[[1]]
    if (length(m) < 2) return(NA_real_)
    as.numeric(m[2])
  }
  
  tibble(
    category = c("Single mapper", "Multi-mapper"),
    pct      = c(extract_pct("aligned concordantly exactly 1 time"),
                 extract_pct("aligned concordantly >1 times"))
  )
}

mapper_breakdown <- dplyr::bind_rows(lapply(names(HISAT2_LOG_FILES), function(s) {
  parse_hisat2_mapper_breakdown(file.path(HISAT2_LOG_DIR, HISAT2_LOG_FILES[s])) %>%
    dplyr::mutate(sample = s)
})) %>%
  dplyr::mutate(
    sample   = factor(sample, levels = names(HISAT2_LOG_FILES)),
    category = factor(category, levels = c("Multi-mapper", "Single mapper"))
  )

p_mapping_pct <- ggplot(mapper_breakdown, aes(x = sample, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = c("Single mapper" = "#44AA99", "Multi-mapper" = "#CC6677"), name = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "% of reads",
       title = "HISAT2 concordant alignment breakdown per library",
       subtitle = paste0("Single- vs multi-mapping reads, from pre-filtering HISAT2 logs \u2014 stack total is",
                         " slightly below overall alignment rate (excludes discordant/single-mate recovery)")) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

## Sequencing depth computed entirely from the HISAT2 logs too — same
## fix as the rate panel above. qc_summary.tsv's mapped_reads has the
## identical post-filter problem (and that file has its own separate
## corruption issue, the stray "0" row, still unresolved as of last
## check) — using the HISAT2 logs for both total and mapped reads
## avoids depending on that file for this plot at all.
parse_hisat2_total_reads <- function(log_path) {
  if (!file.exists(log_path)) return(NA_real_)
  lines <- readLines(log_path)
  total_line <- lines[grepl("^[0-9]+ reads; of these:", trimws(lines))]
  if (length(total_line) == 0) return(NA_real_)
  as.numeric(sub("^([0-9]+).*", "\\1", trimws(total_line[1])))
}

hisat2_depth <- hisat2_rates %>%
  dplyr::mutate(
    total_reads  = sapply(HISAT2_LOG_FILES, function(f) parse_hisat2_total_reads(file.path(HISAT2_LOG_DIR, f))),
    mapped_reads = round(total_reads * overall_alignment_rate / 100)
  )

p_read_depth <- hisat2_depth %>%
  dplyr::select(sample, total_reads, mapped_reads) %>%
  tidyr::pivot_longer(cols = c(total_reads, mapped_reads), names_to = "metric", values_to = "n_reads") %>%
  dplyr::mutate(metric = factor(metric, levels = c("total_reads", "mapped_reads"),
                                labels = c("Total reads", "Mapped reads"))) %>%
  ggplot(aes(x = sample, y = n_reads, fill = metric)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("Total reads" = "grey70", "Mapped reads" = "#882255"), name = NULL) +
  labs(x = NULL, y = "Read count", title = "Sequencing depth per library",
       subtitle = "From HISAT2 logs \u2014 total reads and true mapped-read count, not post-filter BAM counts") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

alignment_stats_panel <- p_mapping_pct / p_read_depth +
  patchwork::plot_layout(heights = c(1, 1))

alignment_stats_panel

ggplot2::ggsave(file.path(OUT_DIR, "alignment_stats.png"), alignment_stats_panel,
                width = 8, height = 9, dpi = 300)
cat("Saved:", file.path(OUT_DIR, "alignment_stats.png"), "\n")