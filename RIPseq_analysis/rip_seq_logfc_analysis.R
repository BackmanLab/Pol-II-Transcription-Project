##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Computes log2(P2AB/Input) by feature type, faceted by expression
# and length -> Fig. 10d. Shares a "Step A" BED-building preamble
# with rip_seq_intron_dog_heatmap.R (both scripts build the same
# five feature-type BED interval sets independently before
# diverging) -- if editing that shared logic, check both files.
#
# Step A. Builds
# five feature-type BED interval sets for downstream RPKM
# quantification and log2(P2AB/Input) analysis:
#   Promoter, Enhancer, Exon, Intron, TES
#
# Promoter and Enhancer regions are classified from WT ATAC-seq +
# ChIP-seq narrowPeak calls using the SAME geography-first
# (ChIPseeker) + mark-based logic as the APA pipeline's
# classify_anchors_narrowpeak() (EP_analysis_functions.R) —
# adapted here from paired loop anchors to a flat set of individual
# genomic regions (the union of all ATAC/ChIP peaks), rather than
# rebuilt from scratch.
#
# Exon and intron intervals are NOT rebuilt here — this script only
# produces Promoter/Enhancer/TES BEDs. The existing GTF-derived
# exon/intron BEDs already used by ripseq_exon_intron_rpkm.sh
# (090124_euseq/euseq_pipeline/intron_dog/beds/) are reused as-is
# in the next step.
#
# EXPRESSION-BIN LINKAGE (Promoter/Enhancer): taken directly from
# ChIPseeker::annotatePeak()'s own nearest-gene output (geneId,
# distanceToTSS) — this is already computed as part of the
# classification step, so no separate proximity-join is needed.
#
# STRAND HANDLING: Promoter/Enhancer/TES regions are UNSTRANDED —
# the downstream RPKM step (not in this script) sums BOTH strand
# bedgraphs for these, unlike exon/intron which stay strand-
# specific. A regulatory element doesn't have an inherent RNA
# strand of its own (enhancers especially can produce bidirectional
# eRNA), so no strand assumption is made here. This is a deliberate
# choice, not a default left unexamined.
#
# TES WINDOW: DEFAULTED to centered +/-1kb around the TES, NOT yet
# explicitly confirmed. See TES_WINDOW_BP below for how to switch
# to a downstream-only DoG-style window instead.
##################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(ChIPseeker)
  library(plyranges)
  library(ggplot2)
  library(scales)
})

##--------------------------------------------------------------## 1. Parameters

TXDB_PATH <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/GRCh38.115.sqlite"
CHIP_DIR  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"

ATAC_FILE    <- file.path(CHIP_DIR, "atac_wt.narrowPeak")
H3K4ME3_FILE <- file.path(CHIP_DIR, "H3K4me3.narrowPeak")
H3K27AC_FILE <- file.path(CHIP_DIR, "h3k27ac.narrowPeak")
H3K4ME1_FILE <- file.path(CHIP_DIR, "H3K4me1.narrowPeak")

OUT_DIR <- "/projects/b1042/BackmanLab/Lucas/022025_RIPseq/feature_beds"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TSS_REGION    <- c(-3000, 3000)   # matches the APA script's own convention
## NOTE: a PROXIMAL_ENHANCER_MAX_BP distance cap was briefly added here to
## reclassify GENIC-zone "Proximal-Enhancer" calls further than 10kb from
## any TSS as Distal-Enhancer instead. Reverted — GENIC-zone membership
## (anywhere in a gene's exon/intron/UTR/downstream) is sufficient on its
## own for "Proximal," regardless of raw distance to TSS.
TES_WINDOW_BP <- 1000L            # DEFAULT: centered +/-1kb around TES.
# NOT yet confirmed — for a downstream-only
# DoG-style window instead (matching
# euseq_dog_intron_count.sh's DOG_WINDOW=5000,
# TES to TES+5kb, excluding neighboring-gene
# overlap), this whole section 6 would need
# the strand-aware one-sided window logic
# from that script instead of the simple
# centered window built here.

##--------------------------------------------------------------## 2. Load narrowPeak files
## Reused verbatim from EP_analysis_functions.R (APA pipeline)

load_narrowpeak <- function(path) {
  if (!file.exists(path)) stop("narrowPeak not found: ", path)
  df <- read.table(path, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE,
                   col.names = c("chr","start","end","name","score",
                                 "strand","signalValue","pValue","qValue","peak"))
  GRanges(seqnames = df$chr, ranges = IRanges(df$start + 1L, df$end),
          name = df$name, score = df$score)
}

message("Loading narrowPeak files...")
atac_gr    <- load_narrowpeak(ATAC_FILE)
h3k4me3_gr <- load_narrowpeak(H3K4ME3_FILE)
h3k27ac_gr <- load_narrowpeak(H3K27AC_FILE)
h3k4me1_gr <- load_narrowpeak(H3K4ME1_FILE)

cat("ATAC peaks:    ", length(atac_gr),    "\n")
cat("H3K4me3 peaks: ", length(h3k4me3_gr), "\n")
cat("H3K27ac peaks: ", length(h3k27ac_gr), "\n")
cat("H3K4me1 peaks: ", length(h3k4me1_gr), "\n")

##--------------------------------------------------------------## 3. Candidate region set — union of all marks, merged
## Every place ANY mark shows up becomes a candidate to classify.
## Overlapping peaks across marks are merged into one region via
## reduce() so a locus with, say, both ATAC and H3K27ac isn't
## classified twice as two overlapping regions.

candidate_regions <- GenomicRanges::reduce(c(atac_gr, h3k4me3_gr, h3k27ac_gr, h3k4me1_gr))
cat("\nCandidate regions (union of all marks, merged):", length(candidate_regions), "\n")

##--------------------------------------------------------------## 4. Adapted single-region classifier
## Same geography-first (ChIPseeker) + mark-based logic as the APA
## pipeline's classify_anchors_narrowpeak() — extracted from its
## per-anchor internal logic (which processes chr1/start1/end1 and
## chr2/start2/end2 as a PAIR) and applied here to a single flat
## GRanges of individual regions instead.

classify_regions_narrowpeak <- function(regions_gr, atac_gr, h3k4me3_gr, h3k27ac_gr, txdb,
                                        h3k4me1_gr = NULL, tss_region = c(-3000, 3000)) {
  
  fix_gr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    suppressWarnings(gr)
  }
  
  chipseeker_to_zone <- function(annotation) {
    dplyr::case_when(
      grepl("^Promoter", annotation)                 ~ "PROMOTER",
      grepl("UTR|Exon|Intron|Downstream", annotation) ~ "GENIC",
      TRUE                                            ~ "INTERGENIC"
    )
  }
  
  ## Harmonize seqlevel STYLE between regions_gr and the TxDb, then
  ## RESTRICT regions_gr to only chromosomes the TxDb actually has
  ## annotation for. Genome-wide ATAC/ChIP narrowPeak calls routinely
  ## include alt contigs, unplaced scaffolds, or chrM that a GTF-derived
  ## TxDb has no gene annotation for at all — ChIPseeker::annotatePeak()
  ## fails with "invalid subscript" trying to look up a nearest gene on
  ## a chromosome the TxDb doesn't know about, rather than erroring
  ## cleanly. This didn't surface in the APA script's loop-anchor version
  ## because Hi-C-derived loop coordinates are already restricted to
  ## primary chromosomes by construction — narrowPeak calls are not.
  txdb_seqlevels <- GenomeInfoDb::seqlevels(txdb)
  has_chr_txdb    <- any(grepl("^chr", txdb_seqlevels))
  has_chr_regions <- any(grepl("^chr", as.character(seqnames(regions_gr))))
  if (has_chr_txdb && !has_chr_regions) {
    seqlevels(regions_gr) <- paste0("chr", seqlevels(regions_gr))
  } else if (!has_chr_txdb && has_chr_regions) {
    seqlevels(regions_gr) <- gsub("^chr", "", seqlevels(regions_gr))
  }
  regions_gr <- suppressWarnings(regions_gr)
  
  common_seqlevels <- intersect(seqlevels(regions_gr), txdb_seqlevels)
  n_before <- length(regions_gr)
  regions_gr <- GenomeInfoDb::keepSeqlevels(regions_gr, common_seqlevels, pruning.mode = "coarse")
  n_after <- length(regions_gr)
  if (n_after < n_before) {
    message("  Dropped ", n_before - n_after,
            " regions on chromosomes not present in the TxDb",
            " (alt/random contigs, chrM, etc.) — ", n_after, " regions remain")
  }
  
  atac_f    <- fix_gr(atac_gr,    regions_gr)
  h3k4me3_f <- fix_gr(h3k4me3_gr, regions_gr)
  h3k27ac_f <- fix_gr(h3k27ac_gr, regions_gr)
  h3k4me1_f <- if (!is.null(h3k4me1_gr)) fix_gr(h3k4me1_gr, regions_gr) else NULL
  
  message("classify_regions_narrowpeak: ", length(regions_gr), " candidate regions")
  message("  Running ChIPseeker::annotatePeak...")
  suppressWarnings({
    anno <- ChIPseeker::annotatePeak(
      regions_gr,
      tssRegion                 = tss_region,
      TxDb                      = txdb,
      level                     = "transcript",
      assignGenomicAnnotation   = TRUE,
      genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
      verbose                   = FALSE
    )
  })
  anno_df <- as.data.frame(anno)
  zone    <- chipseeker_to_zone(anno_df$annotation)
  message("  Zone breakdown:"); print(table(zone))
  
  midpoints <- (start(regions_gr) + end(regions_gr)) / 2
  
  has_atac    <- overlapsAny(regions_gr, atac_f)
  has_h3k4me3 <- overlapsAny(regions_gr, h3k4me3_f)
  
  ## H3K27ac — disambiguate multi-hit regions by closest peak center
  ## to region midpoint (same convention as classify_anchors_narrowpeak)
  hits_ac <- findOverlaps(regions_gr, h3k27ac_f)
  has_h3k27ac <- rep(FALSE, length(regions_gr))
  if (length(hits_ac) > 0) {
    ac_df <- data.frame(
      region_idx  = queryHits(hits_ac),
      peak_center = (start(h3k27ac_f)[subjectHits(hits_ac)] + end(h3k27ac_f)[subjectHits(hits_ac)]) / 2
    ) %>%
      dplyr::mutate(dist = abs(peak_center - midpoints[region_idx])) %>%
      dplyr::group_by(region_idx) %>%
      dplyr::slice_min(dist, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    has_h3k27ac[ac_df$region_idx] <- TRUE
  }
  
  has_h3k4me1 <- if (!is.null(h3k4me1_f)) overlapsAny(regions_gr, h3k4me1_f) else rep(FALSE, length(regions_gr))
  
  ## Loosened definition: H3K27ac OR H3K4me1 alone qualifies as an
  ## enhancer mark, no ATAC co-occurrence required. The "not overlapping
  ## a TSS/promoter" exclusion is already handled below via the zone
  ## check (PROMOTER zone is excluded from Proximal/Distal-Enhancer
  ## classification, same as before) — nothing needed to change there.
  ## This is more inclusive than the previous ATAC-gated definition: a
  ## region with just one mark and no confirmed accessibility now
  ## counts equally alongside one with all three signals present.
  has_enhancer_marks <- has_h3k27ac | has_h3k4me1
  
  region_type <- dplyr::case_when(
    zone == "PROMOTER"   & has_h3k4me3        ~ "Promoter",
    zone == "PROMOTER"   & has_enhancer_marks ~ "Proximal-Enhancer",
    zone == "PROMOTER"   & has_atac           ~ "Open-Chromatin",
    zone == "GENIC"      & has_enhancer_marks ~ "Proximal-Enhancer",
    zone == "GENIC"      & has_atac           ~ "Open-Chromatin",
    zone == "INTERGENIC" & has_enhancer_marks ~ "Distal-Enhancer",
    zone == "INTERGENIC" & has_atac           ~ "Open-Chromatin",
    TRUE ~ "Unannotated"
  )
  ## GENIC-zone membership (anywhere in a gene's exon/intron/UTR/
  ## downstream) qualifies as "Proximal" on its own, regardless of raw
  ## distance to TSS — no distance-based correction applied here.
  message("  Region type breakdown:"); print(table(region_type))
  
  data.frame(
    chr             = as.character(seqnames(regions_gr)),
    start           = start(regions_gr) - 1L,   # back to 0-based BED
    end             = end(regions_gr),
    region_type     = region_type,
    zone            = zone,
    nearest_gene_id = sub("\\..*", "", anno_df$geneId),
    dist_to_tss     = anno_df$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

message("\nLoading TxDb...")
txdb <- AnnotationDbi::loadDb(TXDB_PATH)

message("Classifying candidate regions...")
classified_regions <- classify_regions_narrowpeak(
  regions_gr = candidate_regions,
  atac_gr    = atac_gr,
  h3k4me3_gr = h3k4me3_gr,
  h3k27ac_gr = h3k27ac_gr,
  txdb       = txdb,
  h3k4me1_gr = h3k4me1_gr,
  tss_region = TSS_REGION
)

##--------------------------------------------------------------## 5. Split into Promoter and Enhancer BEDs
## 4-column BED (chr/start/end/nearest_gene_id) — no strand column,
## since these regions are treated as unstranded downstream (see
## docstring). "name" here means nearest gene, kept for the
## expression-bin join in the analysis step, NOT structural gene
## membership the way exon/intron BEDs' gene_id column means it.

promoter_regions <- classified_regions %>% dplyr::filter(region_type == "Promoter")
enhancer_regions <- classified_regions %>% dplyr::filter(region_type %in% c("Proximal-Enhancer","Distal-Enhancer"))

cat("\nPromoter regions:", nrow(promoter_regions), "\n")
cat("Enhancer regions (proximal + distal):", nrow(enhancer_regions), "\n")
cat("  Proximal:", sum(classified_regions$region_type == "Proximal-Enhancer"), "\n")
cat("  Distal:  ", sum(classified_regions$region_type == "Distal-Enhancer"), "\n")

write.table(
  promoter_regions %>% dplyr::select(chr, start, end, nearest_gene_id),
  file.path(OUT_DIR, "promoter_regions.bed"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  enhancer_regions %>% dplyr::select(chr, start, end, nearest_gene_id),
  file.path(OUT_DIR, "enhancer_regions.bed"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

## Full classified table (all region types, all annotation columns) —
## kept as a QC/reference file, not used directly downstream
write_tsv(classified_regions, file.path(OUT_DIR, "all_classified_regions.tsv"))

cat("\nWritten:\n")
cat(" ", file.path(OUT_DIR, "promoter_regions.bed"), "\n")
cat(" ", file.path(OUT_DIR, "enhancer_regions.bed"), "\n")
cat(" ", file.path(OUT_DIR, "all_classified_regions.tsv"), "(full classification, QC reference)\n")

##--------------------------------------------------------------## 6. TES windows — DEFAULT: centered +/- TES_WINDOW_BP
## NOT yet explicitly confirmed — see TES_WINDOW_BP parameter note.
## gene_id here IS structural (this window belongs to this gene's
## TES), unlike Promoter/Enhancer's nearest_gene_id column above.

message("\nBuilding TES windows (+/- ", TES_WINDOW_BP, "bp, centered)...")

genes_gr <- GenomicFeatures::genes(txdb)
tes_pos  <- ifelse(as.character(strand(genes_gr)) == "+", end(genes_gr), start(genes_gr))

tes_df <- data.frame(
  chr     = as.character(seqnames(genes_gr)),
  start   = pmax(0L, tes_pos - TES_WINDOW_BP - 1L),   # -1 for BED 0-based start
  end     = tes_pos + TES_WINDOW_BP,
  gene_id = sub("\\..*", "", genes_gr$gene_id),
  stringsAsFactors = FALSE
)

cat("TES windows built for:", nrow(tes_df), "genes\n")

write.table(
  tes_df, file.path(OUT_DIR, "tes_windows.bed"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
cat("Written:", file.path(OUT_DIR, "tes_windows.bed"), "\n")

message("\nFeature BED construction complete.")
message("Next step: bedtools map these against the P2AB/WT_Input bedgraphs,")
message("mirroring ripseq_exon_intron_rpkm.sh's compute_interval_rpkm() —")
message("but summing BOTH strand bedgraphs for these unstranded regions,")
message("rather than the strand-specific split used for exon/intron.")

##--------------------------------------------------------------## 7. Diagnostic plots
## Nothing above this point produced a visual — everything was
## console-printed counts. These three plots are the actual
## sanity checks for whether the classification looks reasonable
## before trusting it downstream.

message("\nBuilding diagnostic plots...")

base_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11)
  )

## 7a. Region-type breakdown — how many candidate regions landed in
## each category. If Promoter or Enhancer looks suspiciously tiny or
## huge relative to the others, that's the first thing to notice.
region_type_levels <- c("Promoter","Proximal-Enhancer","Distal-Enhancer",
                        "Open-Chromatin","Unannotated")
region_type_colors <- c(
  "Promoter"           = "#882255",
  "Proximal-Enhancer"  = "#44AA99",
  "Distal-Enhancer"    = "#4393c3",
  "Open-Chromatin"     = "grey60",
  "Unannotated"        = "grey85"
)

p_region_type <- classified_regions %>%
  dplyr::mutate(region_type = factor(region_type, levels = region_type_levels)) %>%
  dplyr::count(region_type) %>%
  ggplot(aes(x = region_type, y = n, fill = region_type)) +
  geom_col() +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = region_type_colors, guide = "none") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL, y = "Number of regions",
    title    = "Candidate region classification breakdown",
    subtitle = paste0("n = ", scales::comma(nrow(classified_regions)),
                      " candidate regions (union of ATAC/H3K4me3/H3K27ac/H3K4me1 peaks)")
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_region_type)

## 7b. Nearest-gene match rate for Promoter and Enhancer specifically
## — these are the two categories whose expression-bin linkage
## downstream depends entirely on nearest_gene_id being non-NA.
match_rate_df <- classified_regions %>%
  dplyr::filter(region_type %in% c("Promoter","Proximal-Enhancer","Distal-Enhancer")) %>%
  dplyr::mutate(
    region_type = factor(region_type, levels = c("Promoter","Proximal-Enhancer","Distal-Enhancer")),
    has_gene    = ifelse(is.na(nearest_gene_id), "No nearest gene (NA)", "Nearest gene assigned")
  ) %>%
  dplyr::count(region_type, has_gene)

p_gene_match <- ggplot(match_rate_df, aes(x = region_type, y = n, fill = has_gene)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("Nearest gene assigned" = "#44AA99",
                               "No nearest gene (NA)"  = "#CC6677"), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = NULL, y = "Fraction of regions",
    title    = "Nearest-gene assignment rate, by region type",
    subtitle = "From ChIPseeker::annotatePeak()'s own geneId output \u2014 this is what the expression-bin join relies on"
  ) +
  base_theme

print(p_gene_match)

## 7c. Distance-to-TSS distribution for Proximal vs Distal enhancers —
## checks whether the zone-based Proximal/Distal split actually
## separates TSS-adjacent from genuinely distal regions the way the
## names imply, rather than overlapping heavily.
p_dist_tss <- classified_regions %>%
  dplyr::filter(region_type %in% c("Proximal-Enhancer","Distal-Enhancer"),
                !is.na(dist_to_tss)) %>%
  dplyr::mutate(region_type = factor(region_type, levels = c("Proximal-Enhancer","Distal-Enhancer"))) %>%
  ggplot(aes(x = abs(dist_to_tss), color = region_type, fill = region_type)) +
  geom_density(alpha = 0.3, linewidth = 0.9) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = c("Proximal-Enhancer" = "#44AA99",
                                "Distal-Enhancer"   = "#4393c3"), name = NULL) +
  scale_fill_manual(values = c("Proximal-Enhancer" = "#44AA99",
                               "Distal-Enhancer"   = "#4393c3"), name = NULL) +
  labs(
    x = "Distance to nearest TSS, bp (log10 scale)",
    y = "Density",
    title    = "Distance-to-TSS: Proximal vs Distal enhancer regions",
    subtitle = "Overlap here is expected \u2014 Proximal includes any GENIC-zone region (in a gene body) regardless of raw TSS distance, not just TSS-adjacent ones"
  ) +
  base_theme

print(p_dist_tss)

message("Diagnostic plots complete.")

##--------------------------------------------------------------## 8. Load RIP-seq bedgraphs
## Reuses the RPM-normalized, strand-specific bedgraphs already on
## disk from ripseq_exon_intron_rpkm.sh — no BAM reprocessing.

RIP_PIPELINE_DIR <- "/projects/b1042/BackmanLab/Lucas/022025_RIPseq/ripseq_pipeline"
BGR_DIR <- file.path(RIP_PIPELINE_DIR, "bedgraphs")

BEDGRAPH_SAMPLES <- c("P2AB_Rep2", "P2AB_Rep3", "WT_Input")

read_bedgraph <- function(path) {
  df <- read_tsv(path, col_names = c("chr","start","end","score"),
                 col_types = "ciid", show_col_types = FALSE)
  GRanges(seqnames = df$chr, ranges = IRanges(df$start + 1L, df$end), score = df$score)
}

message("\nLoading RIP-seq bedgraphs...")
bedgraphs <- setNames(lapply(BEDGRAPH_SAMPLES, function(s) {
  list(
    plus  = read_bedgraph(file.path(BGR_DIR, paste0(s, "_plus.norm.bgr"))),
    minus = read_bedgraph(file.path(BGR_DIR, paste0(s, "_minus.norm.bgr")))
  )
}), BEDGRAPH_SAMPLES)

## Combined (unstranded) bedgraph per sample — for Promoter/Enhancer/TES,
## which are treated as unstranded (see docstring at the top of this script)
bedgraphs_unstranded <- setNames(
  lapply(bedgraphs, function(b) c(b$plus, b$minus)),
  BEDGRAPH_SAMPLES
)

cat("Bedgraph intervals loaded (plus / minus) per sample:\n")
for (s in BEDGRAPH_SAMPLES) {
  cat("  ", s, ": ", length(bedgraphs[[s]]$plus), " / ", length(bedgraphs[[s]]$minus), "\n", sep = "")
}

##--------------------------------------------------------------## 9. Load exon/intron BEDs (reused, strand-specific — NOT rebuilt)

EUSEQ_BED_DIR <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/euseq_pipeline/intron_dog/beds"

read_feature_bed <- function(path) {
  df <- read_tsv(path, col_names = c("chr","start","end","gene_id","score","strand"),
                 col_types = "ciicdc", show_col_types = FALSE)
  GRanges(seqnames = df$chr, ranges = IRanges(df$start + 1L, df$end),
          feature_id = sub("\\..*", "", df$gene_id))
}

message("Loading exon/intron BEDs (reused from EU-seq pipeline)...")
exon_plus_gr    <- read_feature_bed(file.path(EUSEQ_BED_DIR, "exons_plus.bed"))
exon_minus_gr   <- read_feature_bed(file.path(EUSEQ_BED_DIR, "exons_minus.bed"))
intron_plus_gr  <- read_feature_bed(file.path(EUSEQ_BED_DIR, "hg38_introns_plus.bed"))
intron_minus_gr <- read_feature_bed(file.path(EUSEQ_BED_DIR, "hg38_introns_minus.bed"))

## Stable per-original-intron ID, assigned BEFORE section 10c's enhancer
## subtraction can split some introns into multiple fragments. This ID
## must exist prior to the split so it gets replicated onto every
## resulting fragment (subtraction already copies all existing metadata
## onto fragments) — letting RPKM computation later re-merge fragments
## from the same original intron by grouping on this ID, instead of
## treating each fragment as its own unrelated feature.
mcols(intron_plus_gr)$interval_id  <- paste0("intron_plus_",  seq_along(intron_plus_gr))
mcols(intron_minus_gr)$interval_id <- paste0("intron_minus_", seq_along(intron_minus_gr))

##--------------------------------------------------------------## 10. Prepare Promoter/Enhancer/TES GRanges
## Promoter/Enhancer: feature_id is a unique REGION id (each region is
## its own row in the final analysis, not aggregated per gene) —
## nearest_gene_id is carried as separate metadata for the expression/
## length-bin join later, NOT used as the RPKM aggregation key.
## TES: feature_id IS gene_id directly (one window per gene, structural).

promoter_gr <- GRanges(
  seqnames        = promoter_regions$chr,
  ranges          = IRanges(promoter_regions$start + 1L, promoter_regions$end),
  feature_id      = paste0("promoter_", seq_len(nrow(promoter_regions))),
  nearest_gene_id = promoter_regions$nearest_gene_id
)
enhancer_gr <- GRanges(
  seqnames        = enhancer_regions$chr,
  ranges          = IRanges(enhancer_regions$start + 1L, enhancer_regions$end),
  feature_id      = paste0("enhancer_", seq_len(nrow(enhancer_regions))),
  nearest_gene_id = enhancer_regions$nearest_gene_id
)
tes_gr <- GRanges(
  seqnames   = tes_df$chr,
  ranges     = IRanges(tes_df$start + 1L, tes_df$end),
  feature_id = tes_df$gene_id
)

## Lookup tables for joining nearest_gene_id back after RPKM computation
## strips extra metadata (compute_rpkm() below only tracks feature_id)
promoter_gene_lookup <- tibble(
  feature_id      = mcols(promoter_gr)$feature_id,
  nearest_gene_id = mcols(promoter_gr)$nearest_gene_id
)
enhancer_gene_lookup <- tibble(
  feature_id      = mcols(enhancer_gr)$feature_id,
  nearest_gene_id = mcols(enhancer_gr)$nearest_gene_id
)

##--------------------------------------------------------------## 10b. Harmonize chromosome-naming style against the BEDGRAPHS
## promoter_gr/enhancer_gr were harmonized to match the TxDb's OWN
## seqlevel style back in section 4 (the ChIPseeker "invalid subscript"
## fix) — never checked against the bedgraphs. tes_gr inherits the
## TxDb's native style directly from genes(txdb) with no harmonization
## at all. If the TxDb's style differs from the bedgraphs' (e.g. "1" vs
## "chr1"), join_overlap_intersect() below would silently find ZERO
## overlaps for all three, every feature would get WT_Input_rpkm=0, and
## the baseline-anchor filter would wipe them out entirely — while
## Exon/Intron, which are independently hardcoded chr-prefixed by the
## original EU-seq bash script, would be unaffected. This is exactly
## that failure mode; fixed here by harmonizing against the bedgraph
## directly rather than the TxDb.

fix_seqlevel_style <- function(gr, ref_gr) {
  has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
  has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
  if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
  else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
  suppressWarnings(gr)
}

reference_bedgraph <- bedgraphs_unstranded$WT_Input

cat("\nSeqlevel style check (first few seqnames of each, BEFORE fix):\n")
cat("  Bedgraph (reference): ", paste(head(as.character(seqnames(reference_bedgraph)), 3), collapse=", "), "\n")
cat("  promoter_gr:          ", paste(head(as.character(seqnames(promoter_gr)), 3), collapse=", "), "\n")
cat("  enhancer_gr:          ", paste(head(as.character(seqnames(enhancer_gr)), 3), collapse=", "), "\n")
cat("  tes_gr:               ", paste(head(as.character(seqnames(tes_gr)), 3), collapse=", "), "\n")
cat("  exon_plus_gr (for comparison, should already match bedgraph):",
    paste(head(as.character(seqnames(exon_plus_gr)), 3), collapse=", "), "\n")

promoter_gr <- fix_seqlevel_style(promoter_gr, reference_bedgraph)
enhancer_gr <- fix_seqlevel_style(enhancer_gr, reference_bedgraph)
tes_gr      <- fix_seqlevel_style(tes_gr,      reference_bedgraph)

cat("\nAFTER fix:\n")
cat("  promoter_gr:", paste(head(as.character(seqnames(promoter_gr)), 3), collapse=", "), "\n")
cat("  enhancer_gr:", paste(head(as.character(seqnames(enhancer_gr)), 3), collapse=", "), "\n")
cat("  tes_gr:     ", paste(head(as.character(seqnames(tes_gr)), 3), collapse=", "), "\n")

##--------------------------------------------------------------## 10c. Subtract Enhancer regions from Intron territory
## An enhancer sitting inside an intron would otherwise have its RIP-
## seq coverage counted TWICE — once toward Enhancer RPKM, once toward
## Intron RPKM, over the exact same base pairs. Plain GenomicRanges::
## setdiff() can't fix this directly: applied to whole GRanges objects
## it collapses everything into a flat disjoint set and DROPS metadata
## entirely, so a single intron split by an enhancer in its middle
## would lose track of which gene either resulting piece belongs to.
##
## This processes each AFFECTED intron individually (only introns that
## actually overlap an enhancer — the rest pass through untouched) and
## re-attaches that intron's own feature_id/gene_id to every piece that
## survives the subtraction, including cases where an intron gets split
## into two pieces by an enhancer in its middle, or fully consumed by
## an enhancer spanning its entire length (that intron then contributes
## nothing, correctly, rather than double-counting the enhancer's
## territory as intron signal too).

subtract_enhancer_from_introns <- function(intron_gr, enhancer_gr) {
  t0 <- Sys.time()
  enhancer_gr_matched <- fix_seqlevel_style(enhancer_gr, intron_gr)
  cat("  [timing] seqlevel fix:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
  
  t1 <- Sys.time()
  hits <- GenomicRanges::findOverlaps(intron_gr, enhancer_gr_matched, ignore.strand = TRUE)
  overlapping_idx <- unique(S4Vectors::queryHits(hits))
  cat("  [timing] findOverlaps:", round(as.numeric(Sys.time() - t1, units = "secs"), 1), "s\n")
  
  cat("  Introns overlapping >=1 enhancer region:", length(overlapping_idx),
      "/", length(intron_gr), "\n")
  
  if (length(overlapping_idx) == 0) {
    return(intron_gr)
  }
  
  unaffected_gr <- intron_gr[-overlapping_idx]
  affected_gr   <- intron_gr[overlapping_idx]
  bp_before     <- sum(width(affected_gr))
  
  ## Build enhancer_groups via ONE bulk subset + one vectorized split(),
  ## NOT a per-intron lapply()+GRangesList() — that combination is a
  ## well-known Bioconductor slow path: building a GRangesList by
  ## looping and combining thousands of individually-built pieces one
  ## at a time, instead of one bulk operation followed by one
  ## C-optimized regroup. This was the actual remaining bottleneck —
  ## psetdiff() itself was already vectorized, but its INPUT was still
  ## being assembled through exactly the kind of loop it was meant to
  ## replace.
  t2 <- Sys.time()
  grouping_factor        <- factor(S4Vectors::queryHits(hits), levels = overlapping_idx)
  relevant_enhancer_hits <- enhancer_gr_matched[S4Vectors::subjectHits(hits)]
  enhancer_groups        <- S4Vectors::split(relevant_enhancer_hits, grouping_factor)
  cat("  [timing] build enhancer_groups:", round(as.numeric(Sys.time() - t2, units = "secs"), 1), "s\n")
  
  ## VERIFICATION NOTE: psetdiff() is a less commonly used part of the
  ## GenomicRanges API and this couldn't be run/tested directly before
  ## handing it over — compare the two diagnostic numbers below against
  ## a previous verified run on the same data; they should match exactly.
  t3 <- Sys.time()
  remaining_grl <- GenomicRanges::psetdiff(affected_gr, enhancer_groups, ignore.strand = TRUE)
  cat("  [timing] psetdiff:", round(as.numeric(Sys.time() - t3, units = "secs"), 1), "s\n")
  
  ## Reattach each original intron's own metadata to every resulting
  ## piece — an intron split into 2 pieces gets its metadata repeated
  ## twice; one fully consumed by an enhancer contributes 0 pieces.
  t4 <- Sys.time()
  n_pieces       <- S4Vectors::elementNROWS(remaining_grl)
  remaining_flat <- unlist(remaining_grl, use.names = FALSE)
  if (length(remaining_flat) > 0) {
    mcols(remaining_flat) <- mcols(affected_gr)[rep(seq_along(affected_gr), n_pieces), , drop = FALSE]
  }
  cat("  [timing] unlist + reattach metadata:", round(as.numeric(Sys.time() - t4, units = "secs"), 1), "s\n")
  
  n_fully_consumed <- sum(n_pieces == 0)
  bp_after <- sum(width(remaining_flat))
  cat("  Fully consumed by an enhancer (dropped entirely):", n_fully_consumed, "\n")
  cat("  Total bp removed from affected introns:", bp_before - bp_after, "\n")
  
  if (length(remaining_flat) == 0) {
    return(unaffected_gr)
  }
  
  c(unaffected_gr, remaining_flat)
}

message("\nSubtracting Enhancer regions from Intron territory (plus strand)...")
intron_plus_gr <- subtract_enhancer_from_introns(intron_plus_gr, enhancer_gr)
message("Subtracting Enhancer regions from Intron territory (minus strand)...")
intron_minus_gr <- subtract_enhancer_from_introns(intron_minus_gr, enhancer_gr)

##--------------------------------------------------------------## 11. RPKM computation via plyranges — bp-weighted overlap sum
## join_overlap_intersect() keeps ONLY the intersected (truly
## overlapping) portion of each feature/bedgraph-interval pair, so
## weighting by width here uses the TRUE overlap width, not the full
## width of either the feature or the bedgraph interval — this is the
## exact bp-weighted "sum(coverage x bp)" the RPKM formula calls for,
## by construction, regardless of how any other tool might compute it.
##
## Total feature length (the RPKM denominator) is computed separately,
## directly from the ORIGINAL un-intersected feature intervals — using
## the overlap result for this would silently undercount any truly
## zero-coverage bp within the feature, since genomecov -bg bedgraphs
## only report non-zero-coverage positions in the first place.

compute_rpkm <- function(features_gr, bedgraph_gr) {
  total_length <- tibble(
    feature_id = mcols(features_gr)$feature_id,
    total_bp   = width(features_gr)
  ) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::summarise(total_bp = sum(total_bp), .groups = "drop")
  
  overlaps <- plyranges::join_overlap_intersect(features_gr, bedgraph_gr)
  
  if (length(overlaps) == 0) {
    cov_sum <- tibble(feature_id = character(0), cov_bp_sum = numeric(0))
  } else {
    cov_sum <- tibble(
      feature_id = mcols(overlaps)$feature_id,
      cov_bp     = mcols(overlaps)$score * width(overlaps)
    ) %>%
      dplyr::group_by(feature_id) %>%
      dplyr::summarise(cov_bp_sum = sum(cov_bp), .groups = "drop")
  }
  
  total_length %>%
    dplyr::left_join(cov_sum, by = "feature_id") %>%
    dplyr::mutate(
      cov_bp_sum = tidyr::replace_na(cov_bp_sum, 0),
      rpkm       = cov_bp_sum / (total_bp / 1000)
    )
}

## Every feature type is counted and quantified at the INDIVIDUAL
## INTERVAL level — one row per promoter region, one row per enhancer
## region, one row per exon, one row per intron, one row per TES
## window. None of them are aggregated into a per-gene sum before RPKM
## is computed. Each interval gets a UNIQUE feature_id (not gene_id),
## with gene_id carried as separate metadata for the expression-bin
## join later — this is the SAME pattern Promoter/Enhancer already
## used (their own feature_id + nearest_gene_id), now applied
## consistently to Exon and Intron too instead of the gene-level
## aggregation used previously.
##
## No strand-aware AGGREGATION step is needed here (unlike the earlier
## compute_rpkm_stranded_combined() this replaces) — plus and minus
## intervals never share a feature_id under this unique-per-row scheme,
## so their RPKM results are simply stacked, not summed together.
compute_rpkm_stranded_individual <- function(plus_gr, minus_gr, plus_bgr, minus_bgr) {
  dplyr::bind_rows(
    compute_rpkm(plus_gr,  plus_bgr),
    compute_rpkm(minus_gr, minus_bgr)
  )
}

## Assign a unique feature_id per individual exon/intron interval,
## preserving the original gene_id as separate metadata for later
## joining. intron_plus_gr/intron_minus_gr here are already the
## enhancer-subtracted versions from section 10c.
exon_plus_gr_indiv    <- exon_plus_gr
exon_minus_gr_indiv   <- exon_minus_gr
intron_plus_gr_indiv  <- intron_plus_gr
intron_minus_gr_indiv <- intron_minus_gr

mcols(exon_plus_gr_indiv)$gene_id     <- mcols(exon_plus_gr_indiv)$feature_id
mcols(exon_plus_gr_indiv)$feature_id    <- paste0("exon_plus_",    seq_along(exon_plus_gr_indiv))
mcols(exon_minus_gr_indiv)$gene_id    <- mcols(exon_minus_gr_indiv)$feature_id
mcols(exon_minus_gr_indiv)$feature_id   <- paste0("exon_minus_",   seq_along(exon_minus_gr_indiv))

## Intron uses the STABLE interval_id assigned before section 10c's
## subtraction, NOT a fresh per-current-row ID — fragments split from
## the same original intron share this ID, so compute_rpkm()'s
## group_by(feature_id) correctly re-merges them (summing coverage and
## length across fragments) before computing RPKM, instead of treating
## each small fragment as its own independent feature.
mcols(intron_plus_gr_indiv)$gene_id    <- mcols(intron_plus_gr_indiv)$feature_id
mcols(intron_plus_gr_indiv)$feature_id   <- mcols(intron_plus_gr_indiv)$interval_id
mcols(intron_minus_gr_indiv)$gene_id   <- mcols(intron_minus_gr_indiv)$feature_id
mcols(intron_minus_gr_indiv)$feature_id  <- mcols(intron_minus_gr_indiv)$interval_id

## Lookup tables for joining gene_id (and, for introns, each interval's
## own length) back after RPKM computation strips extra metadata.
## Intron specifically must SUM fragment widths per feature_id here —
## multiple fragments now share one feature_id, so a naive per-row
## width would create duplicate feature_id entries and multiply rows
## on the later join.
exon_gene_lookup <- dplyr::bind_rows(
  tibble(feature_id = mcols(exon_plus_gr_indiv)$feature_id,
         gene_id    = mcols(exon_plus_gr_indiv)$gene_id),
  tibble(feature_id = mcols(exon_minus_gr_indiv)$feature_id,
         gene_id    = mcols(exon_minus_gr_indiv)$gene_id)
)
intron_gene_lookup <- dplyr::bind_rows(
  tibble(feature_id = mcols(intron_plus_gr_indiv)$feature_id,
         gene_id    = mcols(intron_plus_gr_indiv)$gene_id,
         frag_bp    = width(intron_plus_gr_indiv)),
  tibble(feature_id = mcols(intron_minus_gr_indiv)$feature_id,
         gene_id    = mcols(intron_minus_gr_indiv)$gene_id,
         frag_bp    = width(intron_minus_gr_indiv))
) %>%
  dplyr::group_by(feature_id, gene_id) %>%
  dplyr::summarise(intron_bp = sum(frag_bp), .groups = "drop")

message("\nComputing RPKM for all 5 feature types x 3 samples (interval-level throughout)...")

rpkm_results <- list()
for (s in BEDGRAPH_SAMPLES) {
  message("  ", s, "...")
  rpkm_results[[s]] <- list(
    promoter = compute_rpkm(promoter_gr, bedgraphs_unstranded[[s]]),
    enhancer = compute_rpkm(enhancer_gr, bedgraphs_unstranded[[s]]),
    tes      = compute_rpkm(tes_gr,      bedgraphs_unstranded[[s]]),
    exon     = compute_rpkm_stranded_individual(exon_plus_gr_indiv,   exon_minus_gr_indiv,
                                                bedgraphs[[s]]$plus, bedgraphs[[s]]$minus),
    intron   = compute_rpkm_stranded_individual(intron_plus_gr_indiv, intron_minus_gr_indiv,
                                                bedgraphs[[s]]$plus, bedgraphs[[s]]$minus)
  )
}

##--------------------------------------------------------------## 12. Average P2AB replicates, compute log2FC, anchor on WT_Input
## Same convention as ripseq_pol2_enrichment_analysis.R: PSEUDOCOUNT=1,
## P2AB_Rep2/Rep3 averaged at the RPKM level before the ratio, gene set
## anchored on WT_Input (baseline) > 0.

PSEUDOCOUNT <- 1

message("\nComputing log2FC(P2AB/Input) per feature type...")

feature_types <- c("promoter","enhancer","tes","exon","intron")

log2fc_by_type <- setNames(lapply(feature_types, function(ft) {
  p2ab_mean <- dplyr::inner_join(
    rpkm_results$P2AB_Rep2[[ft]] %>% dplyr::select(feature_id, rpkm_rep2 = rpkm),
    rpkm_results$P2AB_Rep3[[ft]] %>% dplyr::select(feature_id, rpkm_rep3 = rpkm),
    by = "feature_id"
  ) %>%
    dplyr::mutate(P2AB_rpkm = rowMeans(cbind(rpkm_rep2, rpkm_rep3), na.rm = TRUE)) %>%
    dplyr::select(feature_id, P2AB_rpkm)
  
  input_rpkm <- rpkm_results$WT_Input[[ft]] %>%
    dplyr::select(feature_id, WT_Input_rpkm = rpkm)
  
  combined <- dplyr::inner_join(p2ab_mean, input_rpkm, by = "feature_id") %>%
    dplyr::filter(WT_Input_rpkm > 0) %>%   # baseline-anchored, per established convention
    dplyr::mutate(
      log2fc      = log2((P2AB_rpkm + PSEUDOCOUNT) / (WT_Input_rpkm + PSEUDOCOUNT)),
      feature_type = ft
    )
  
  combined
}), feature_types)

for (ft in feature_types) {
  cat("  ", ft, ": ", nrow(log2fc_by_type[[ft]]), " features (WT_Input-anchored)\n", sep = "")
}

## Attach gene_id: Promoter/Enhancer via nearest-gene lookup (region-
## level features, gene is "nearest," not structural); Exon/Intron via
## their own gene_id lookup (structural — each interval genuinely
## belongs to one gene, now that feature_id is a unique per-interval
## ID rather than the gene_id itself); TES stays as feature_id==gene_id
## since it's still structurally one window per gene.
log2fc_by_type$promoter <- log2fc_by_type$promoter %>%
  dplyr::left_join(promoter_gene_lookup, by = "feature_id") %>%
  dplyr::rename(gene_id = nearest_gene_id)
log2fc_by_type$enhancer <- log2fc_by_type$enhancer %>%
  dplyr::left_join(enhancer_gene_lookup, by = "feature_id") %>%
  dplyr::rename(gene_id = nearest_gene_id)
log2fc_by_type$tes    <- log2fc_by_type$tes    %>% dplyr::mutate(gene_id = feature_id)
log2fc_by_type$exon   <- log2fc_by_type$exon   %>% dplyr::left_join(exon_gene_lookup,   by = "feature_id")
log2fc_by_type$intron <- log2fc_by_type$intron %>% dplyr::left_join(intron_gene_lookup, by = "feature_id")

log2fc_combined <- dplyr::bind_rows(log2fc_by_type) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id))

cat("\nCombined log2FC table:", nrow(log2fc_combined), "rows across all feature types\n")
cat("Rows with no gene_id (nearest-gene lookup failed):",
    sum(is.na(log2fc_combined$gene_id)), "\n")

##--------------------------------------------------------------## 13. Expression bins (WT TPM tertiles) and gene-length bins

EXPR_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
EXPR_CATS <- c("Low/Unexpressed","Medium","High")

message("\nLoading WT TPM and computing expression tertiles...")
tpm_df <- read_tsv(EXPR_FILE, show_col_types = FALSE) %>%
  dplyr::select(gene_id, TPM) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id))

expressed_tpm <- tpm_df %>% dplyr::filter(!is.na(TPM), TPM > 0) %>% dplyr::pull(TPM)
tpm_medium <- quantile(expressed_tpm, 1/3, na.rm = TRUE)
tpm_high   <- quantile(expressed_tpm, 2/3, na.rm = TRUE)

tpm_df <- tpm_df %>%
  dplyr::mutate(
    expr_bin = dplyr::case_when(
      is.na(TPM) | TPM == 0 ~ "Low/Unexpressed",
      TPM <= tpm_medium     ~ "Low/Unexpressed",
      TPM <= tpm_high       ~ "Medium",
      TRUE                  ~ "High"
    ),
    expr_bin = factor(expr_bin, levels = EXPR_CATS)
  )

## Gene-length bins — total gene BODY span (TSS to TES), tertile-binned.
## Distinct from intron length; reuses genes_gr already built in section 6.
GENE_LENGTH_CATS <- c("Short","Average","Long")

message("Computing gene-length tertiles...")
gene_length_df <- tibble(
  gene_id    = sub("\\..*", "", genes_gr$gene_id),
  gene_bp    = width(genes_gr)
)
len_medium <- quantile(gene_length_df$gene_bp, 1/3, na.rm = TRUE)
len_high   <- quantile(gene_length_df$gene_bp, 2/3, na.rm = TRUE)
gene_length_df <- gene_length_df %>%
  dplyr::mutate(
    length_bin = dplyr::case_when(
      gene_bp <= len_medium ~ "Short",
      gene_bp <= len_high   ~ "Average",
      TRUE                  ~ "Long"
    ),
    length_bin = factor(length_bin, levels = GENE_LENGTH_CATS)
  )

cat("Gene length tertile breaks: Short <=", round(len_medium), "bp | Average <=",
    round(len_high), "bp | Long >", round(len_high), "bp\n")

log2fc_combined <- log2fc_combined %>%
  dplyr::left_join(tpm_df %>% dplyr::select(gene_id, expr_bin), by = "gene_id") %>%
  dplyr::left_join(gene_length_df %>% dplyr::select(gene_id, length_bin), by = "gene_id")

cat("\nRows matched to an expression bin:", sum(!is.na(log2fc_combined$expr_bin)),
    "/", nrow(log2fc_combined), "\n")
cat("Rows matched to a length bin:      ", sum(!is.na(log2fc_combined$length_bin)),
    "/", nrow(log2fc_combined), "\n")

##--------------------------------------------------------------## 14. Boxplots — log2FC by feature type, faceted by expression, then length

message("\nBuilding boxplots...")

FEATURE_TYPE_LEVELS <- c("Promoter","Enhancer","Exon","Intron","TES")
FEATURE_TYPE_COLORS <- c(
  Promoter = "#882255",
  Enhancer = "#44AA99",
  Exon     = "#88CCEE",
  Intron   = "#CC6677",
  TES      = "#DDCC77"
)

plot_df <- log2fc_combined %>%
  dplyr::mutate(
    feature_type = dplyr::recode(feature_type,
                                 promoter = "Promoter", enhancer = "Enhancer",
                                 exon = "Exon", intron = "Intron", tes = "TES"),
    feature_type = factor(feature_type, levels = FEATURE_TYPE_LEVELS)
  )

## 14a. Faceted by expression (Low/Unexpressed, Medium, High)
p_by_expr <- plot_df %>%
  dplyr::filter(!is.na(expr_bin)) %>%
  ggplot(aes(x = feature_type, y = log2fc, color = feature_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.7, width = 0.6) +
  scale_color_manual(values = FEATURE_TYPE_COLORS, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ expr_bin, nrow = 1) +
  labs(
    x = NULL, y = "log2(P2AB / WT_Input)",
    title    = "Pol2 RIP-seq enrichment by feature type, stratified by expression",
    subtitle = "WT_Input-anchored (baseline > 0)  |  P2AB = mean of Rep2/Rep3"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_by_expr)

## 14b. Faceted by gene length (Short, Average, Long)
p_by_length <- plot_df %>%
  dplyr::filter(!is.na(length_bin)) %>%
  ggplot(aes(x = feature_type, y = log2fc, color = feature_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.7, width = 0.6) +
  scale_color_manual(values = FEATURE_TYPE_COLORS, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ length_bin, nrow = 1) +
  labs(
    x = NULL, y = "log2(P2AB / WT_Input)",
    title    = "Pol2 RIP-seq enrichment by feature type, stratified by gene length",
    subtitle = "WT_Input-anchored (baseline > 0)  |  P2AB = mean of Rep2/Rep3  |  gene length = TSS-to-TES span"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_by_length)

message("\nFeature-level enrichment analysis complete.")
message("Key objects: log2fc_combined, rpkm_results")

##--------------------------------------------------------------## 15. Final mirror plot — Low vs High expression only
## Top: boxplots (same as 14a, restricted to Low/Unexpressed + High).
## Bottom: bar chart of n features per type/bin, same colors, same
## facet columns, mirrored via scale_y_reverse() so it reads as a
## reflection of the boxplot panel above it. X-axis category labels
## appear once, on the bottom panel only.

message("\nBuilding final mirror plot (Low vs High expression)...")

mirror_df <- plot_df %>%
  dplyr::filter(!is.na(expr_bin), expr_bin %in% c("Low/Unexpressed","High")) %>%
  dplyr::mutate(expr_bin = factor(expr_bin, levels = c("Low/Unexpressed","High")))

p_mirror_top <- ggplot(mirror_df, aes(x = feature_type, y = log2fc, color = feature_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.7, width = 0.6) +
  scale_color_manual(values = FEATURE_TYPE_COLORS, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ expr_bin, nrow = 1) +
  labs(x = NULL, y = "log2(P2AB / WT_Input)") +
  base_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(b = 2)
  )

mirror_counts <- mirror_df %>%
  dplyr::count(expr_bin, feature_type, .drop = FALSE)

p_mirror_bottom <- ggplot(mirror_counts, aes(x = feature_type, y = n, fill = feature_type)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3, color = "black") +
  scale_fill_manual(values = FEATURE_TYPE_COLORS, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  scale_y_reverse(labels = scales::comma, expand = expansion(mult = c(0.05, 0.15))) +
  facet_wrap(~ expr_bin, nrow = 1) +
  labs(x = NULL, y = "n features") +
  base_theme +
  theme(
    strip.text        = element_blank(),
    strip.background  = element_blank(),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    plot.margin       = margin(t = 2)
  )

mirror_final <- (p_mirror_top / p_mirror_bottom) +
  patchwork::plot_layout(heights = c(3, 1.2)) +
  patchwork::plot_annotation(
    title    = "Pol2 RIP-seq enrichment by feature type \u2014 Low vs High expression",
    subtitle = "WT_Input-anchored (baseline > 0)  |  P2AB = mean of Rep2/Rep3  |  bars below = n features per box",
    theme    = theme(plot.title = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(color = "grey40", size = 11))
  )

print(mirror_final)

message("Mirror plot complete. Key object: mirror_final")

##--------------------------------------------------------------## 16. Second mirror plot — introns by length, Low vs High expression
## No new RPKM computation needed here anymore — section 11/12 now
## computes Intron at the interval level already (feature_id per
## intron, gene_id + intron_bp attached via intron_gene_lookup), so
## this just filters log2fc_combined down to Intron rows and adds
## length-tertile binning on top.

message("\nBuilding intron-length mirror plot...")

INTRON_LENGTH_CATS   <- c("Short","Medium","Long")
INTRON_LENGTH_COLORS <- c(Short = "#fdd49e", Medium = "#fc8d59", Long = "#7f0000")

intron_only <- log2fc_combined %>% dplyr::filter(feature_type == "intron")

intron_len_medium <- quantile(intron_only$intron_bp, 1/3, na.rm = TRUE)
intron_len_high   <- quantile(intron_only$intron_bp, 2/3, na.rm = TRUE)

intron_only <- intron_only %>%
  dplyr::mutate(
    intron_length_bin = dplyr::case_when(
      intron_bp <= intron_len_medium ~ "Short",
      intron_bp <= intron_len_high   ~ "Medium",
      TRUE                           ~ "Long"
    ),
    intron_length_bin = factor(intron_length_bin, levels = INTRON_LENGTH_CATS)
  )

cat("Intron length tertile breaks: Short <=", round(intron_len_medium), "bp | Medium <=",
    round(intron_len_high), "bp | Long >", round(intron_len_high), "bp\n")

## Same mirror-plot structure as section 15, categories swapped to
## intron_length_bin instead of feature_type
mirror_intron_df <- intron_only %>%
  dplyr::filter(!is.na(expr_bin), expr_bin %in% c("Low/Unexpressed","High"),
                !is.na(intron_length_bin)) %>%
  dplyr::mutate(expr_bin = factor(expr_bin, levels = c("Low/Unexpressed","High")))

p_intron_mirror_top <- ggplot(mirror_intron_df, aes(x = intron_length_bin, y = log2fc, color = intron_length_bin)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.7, width = 0.6) +
  scale_color_manual(values = INTRON_LENGTH_COLORS, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ expr_bin, nrow = 1) +
  labs(x = NULL, y = "log2(P2AB / WT_Input)") +
  base_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(b = 2)
  )

intron_mirror_counts <- mirror_intron_df %>%
  dplyr::count(expr_bin, intron_length_bin, .drop = FALSE)

p_intron_mirror_bottom <- ggplot(intron_mirror_counts, aes(x = intron_length_bin, y = n, fill = intron_length_bin)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3, color = "black") +
  scale_fill_manual(values = INTRON_LENGTH_COLORS, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  scale_y_reverse(labels = scales::comma, expand = expansion(mult = c(0.05, 0.15))) +
  facet_wrap(~ expr_bin, nrow = 1) +
  labs(x = NULL, y = "n introns") +
  base_theme +
  theme(
    strip.text        = element_blank(),
    strip.background  = element_blank(),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    plot.margin       = margin(t = 2)
  )

intron_mirror_final <- (p_intron_mirror_top / p_intron_mirror_bottom) +
  patchwork::plot_layout(heights = c(3, 1.2)) +
  patchwork::plot_annotation(
    title    = "Pol2 RIP-seq enrichment within introns, by intron length \u2014 Low vs High expression",
    subtitle = paste0(
      "WT_Input-anchored (baseline > 0)  |  P2AB = mean of Rep2/Rep3  |  ",
      "bars below = n introns per box  |  per-intron, not gene-aggregated"
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(color = "grey40", size = 11))
  )

print(intron_mirror_final)

message("Intron-length mirror plot complete. Key object: intron_mirror_final, intron_only")