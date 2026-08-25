##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Shared functions for loop anchor cCRE classification. Source
# this file at the top of any analysis script:
#   source("EP_analysis_functions.R")
#
# Two classifiers provided:
#   classify_anchors_priority() -- Promoter > Enhancer, ties broken by first hit
#   classify_anchors_midpoint() -- closest cCRE center to anchor midpoint wins
#
# Both exclude CTCF. Pass gr_wt_ep (CTCF-excluded GRanges) as ccre_gr.
#
# Figure output: none directly -- shared functions, sourced by
# other scripts in this folder and adapted by
# RIPseq_analysis/rip_seq_intron_dog_heatmap.R and
# rip_seq_logfc_analysis.R for their own BED-building step.
##################################################################

library(GenomicRanges)
library(dplyr)

# ── Shared constants ───────────────────────────────────────────────────────────

ANCHOR_SLOP_BP <- 2500L

PAIR_ORDER <- c("Promoter","Proximal-Enhancer",
                "Distal-Enhancer","Open-Chromatin","Other","Unannotated")

TYPE_PRIORITY <- c(
  "Promoter"           = 1L,
  "Proximal-Enhancer"  = 2L,
  "Distal-Enhancer"    = 3L,
  "Open-Chromatin"     = 4L,
  "Other"              = 5L
)

GROUP_MAP <- c(
  "Promoter — Promoter"                   = "P-P",
  "Promoter — Proximal-Enhancer"          = "E-P",
  "Promoter — Distal-Enhancer"            = "E-P",
  "Proximal-Enhancer — Proximal-Enhancer" = "E-E",
  "Proximal-Enhancer — Distal-Enhancer"   = "E-E",
  "Distal-Enhancer — Distal-Enhancer"     = "E-E"
)

GROUP_LEVELS <- c("P-P","E-P","E-E")
GROUP_COLORS <- c("P-P"="#882255","E-P"="#44AA99","E-E"="#4393c3")

PE_PAIR_TYPES <- c(
  "Promoter — Promoter",
  "Promoter — Proximal-Enhancer",
  "Promoter — Distal-Enhancer",
  "Proximal-Enhancer — Proximal-Enhancer",
  "Proximal-Enhancer — Distal-Enhancer",
  "Distal-Enhancer — Distal-Enhancer"
)

SINGLETON_LEVELS <- c("P-Gene","pE-Gene","dE-Gene")
SINGLETON_COLORS <- c("P-Gene"="#b2182b","pE-Gene"="#f4a582","dE-Gene"="#d6604d")

# ── cCRE loaders ──────────────────────────────────────────────────────────────

load_ccre_bed <- function(path, label=NULL) {
  df <- read.table(path, header=FALSE, sep="\t",
                   stringsAsFactors=FALSE, fill=TRUE,
                   col.names=c("chr","start","end","accession",
                               "score","strand","thickStart","thickEnd",
                               "rgb","ccre_class","activity"))
  df <- df %>%
    filter(!is.na(ccre_class), ccre_class != "Low-DNase") %>%
    dplyr::mutate(
      ccre_type = dplyr::case_when(
        ccre_class == "PLS"     ~ "Promoter",
        ccre_class == "pELS"    ~ "Proximal-Enhancer",
        ccre_class == "dELS"    ~ "Distal-Enhancer",
        ccre_class == "CA-CTCF" ~ "CTCF",
        ccre_class == "CA-only" ~ "Open-Chromatin",
        TRUE                    ~ "Other"
      ),
      source = if (!is.null(label)) label else basename(path)
    )
  message("Loaded ", nrow(df), " cCREs from ", basename(path))
  print(table(df$ccre_type))
  df
}

make_ccre_gr <- function(df) {
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end),
          accession=df$accession, ccre_type=df$ccre_type,
          ccre_class=df$ccre_class)
}

# ── Shared internal helpers ────────────────────────────────────────────────────

.fix_chr <- function(gr, ref_gr) {
  has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
  has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
  if (has_chr_ref && !has_chr_gr)
    seqlevels(gr) <- paste0("chr", seqlevels(gr))
  else if (!has_chr_ref && has_chr_gr)
    seqlevels(gr) <- gsub("^chr","",seqlevels(gr))
  gr
}

.make_anchor_gr <- function(loops_df, anchor_num, slop_bp) {
  if (anchor_num == 1) {
    chr <- loops_df$chr1; s <- loops_df$start1; e <- loops_df$end1
  } else {
    chr <- loops_df$chr2; s <- loops_df$start2; e <- loops_df$end2
  }
  GRanges(chr, IRanges(pmax(1L, s - slop_bp), e + slop_bp))
}

.make_pair <- function(a, b) {
  a_rank <- match(a, PAIR_ORDER, nomatch=99L)
  b_rank <- match(b, PAIR_ORDER, nomatch=99L)
  if (a_rank <= b_rank) paste(a, b, sep=" — ")
  else                  paste(b, a, sep=" — ")
}

# ── Classifier 1: Priority-based ──────────────────────────────────────────────
# Rule: Promoter > Proximal-Enhancer > Distal-Enhancer > Open-Chromatin
# When multiple cCREs overlap, the highest-priority type wins.
# CTCF excluded — pass gr_wt_ep (ccre_type != "CTCF") as ccre_gr.

classify_anchors_priority <- function(loops_df, ccre_gr, ccre_label="",
                                      slop_bp=ANCHOR_SLOP_BP) {
  classify_one <- function(anchor_gr) {
    anchor_gr_fixed <- .fix_chr(anchor_gr, ccre_gr)
    hits   <- findOverlaps(anchor_gr_fixed, ccre_gr)
    result <- rep("Unannotated", length(anchor_gr))
    if (length(hits) > 0) {
      hits_df <- data.frame(
        anchor_idx = queryHits(hits),
        ccre_type  = ccre_gr$ccre_type[subjectHits(hits)],
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(priority = TYPE_PRIORITY[ccre_type]) %>%
        dplyr::filter(!is.na(priority)) %>%
        dplyr::group_by(anchor_idx) %>%
        dplyr::filter(priority == min(priority)) %>%
        dplyr::slice_head(n=1) %>%
        dplyr::ungroup()
      if (nrow(hits_df) > 0)
        result[hits_df$anchor_idx] <- hits_df$ccre_type
    }
    result
  }
  
  message("classify_anchors_priority: ", nrow(loops_df), " loops against ",
          ccre_label, if (slop_bp > 0) paste0(" (slop=", slop_bp, "bp)") else "")
  
  a1 <- classify_one(.make_anchor_gr(loops_df, 1, slop_bp))
  a2 <- classify_one(.make_anchor_gr(loops_df, 2, slop_bp))
  
  loops_df %>%
    dplyr::mutate(
      anchor1_ccre   = a1,
      anchor2_ccre   = a2,
      loop_ccre_type = mapply(.make_pair, a1, a2, USE.NAMES=FALSE),
      ccre_source    = ccre_label
    )
}

# ── Classifier 2: Midpoint-based ──────────────────────────────────────────────
# Rule: when multiple cCREs overlap an anchor, the one whose center is
# closest to the anchor midpoint wins. No priority hierarchy.
# CTCF excluded — pass gr_wt_ep (ccre_type != "CTCF") as ccre_gr.

classify_anchors_midpoint <- function(loops_df, ccre_gr, ccre_label="",
                                      slop_bp=ANCHOR_SLOP_BP) {
  classify_one <- function(anchor_gr, midpoints) {
    anchor_gr_fixed <- .fix_chr(anchor_gr, ccre_gr)
    hits   <- findOverlaps(anchor_gr_fixed, ccre_gr)
    result <- rep("Unannotated", length(anchor_gr))
    if (length(hits) > 0) {
      hits_df <- data.frame(
        anchor_idx  = queryHits(hits),
        ccre_type   = ccre_gr$ccre_type[subjectHits(hits)],
        ccre_center = (start(ccre_gr)[subjectHits(hits)] +
                         end(ccre_gr)[subjectHits(hits)]) / 2,
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(
          anchor_mid  = midpoints[anchor_idx],
          dist_to_mid = abs(ccre_center - anchor_mid)
        ) %>%
        dplyr::filter(!is.na(ccre_type), ccre_type != "Other") %>%
        dplyr::group_by(anchor_idx) %>%
        dplyr::slice_min(dist_to_mid, n=1, with_ties=FALSE) %>%
        dplyr::ungroup()
      if (nrow(hits_df) > 0)
        result[hits_df$anchor_idx] <- hits_df$ccre_type
    }
    result
  }
  
  mid1 <- (loops_df$start1 + loops_df$end1) / 2
  mid2 <- (loops_df$start2 + loops_df$end2) / 2
  
  message("classify_anchors_midpoint: ", nrow(loops_df), " loops against ",
          ccre_label, if (slop_bp > 0) paste0(" (slop=", slop_bp, "bp)") else "")
  
  a1 <- classify_one(.make_anchor_gr(loops_df, 1, slop_bp), mid1)
  a2 <- classify_one(.make_anchor_gr(loops_df, 2, slop_bp), mid2)
  
  loops_df %>%
    dplyr::mutate(
      anchor1_ccre   = a1,
      anchor2_ccre   = a2,
      loop_ccre_type = mapply(.make_pair, a1, a2, USE.NAMES=FALSE),
      ccre_source    = ccre_label
    )
}

plot_apa_grid_generic <- function(apa_agg, bin_levels, bin_colors,
                                  title="APA", palette=ACTIVE_PALETTE,
                                  corner_size=3) {
  
  # Normalize each matrix by its own corner background
  # Result: 1.0 = background level, higher = loop enrichment
  normalize_to_background <- function(mat) {
    if (is.null(mat)) return(NULL)
    n <- nrow(mat)
    corners <- c(
      as.vector(mat[1:corner_size, 1:corner_size]),
      as.vector(mat[1:corner_size, (n-corner_size+1):n]),
      as.vector(mat[(n-corner_size+1):n, 1:corner_size]),
      as.vector(mat[(n-corner_size+1):n, (n-corner_size+1):n])
    )
    corner_mean <- mean(corners[!is.na(corners)], na.rm=TRUE)
    if (is.na(corner_mean) || corner_mean == 0) return(mat)
    mat / corner_mean
  }
  
  # Normalize all matrices
  apa_norm <- lapply(apa_agg, function(cond_list)
    lapply(cond_list, normalize_to_background)
  )
  
  # Global color scale across all normalized matrices
  all_vals   <- unlist(lapply(apa_norm, function(cl)
    unlist(lapply(cl, function(m) if (!is.null(m)) as.vector(m)))))
  global_max <- quantile(all_vals, 0.98, na.rm=TRUE)
  # Floor at 1.0 so background is always at the bottom of the scale
  global_min <- 1.0
  
  make_panel <- function(mat_norm) {
    if (is.null(mat_norm)) return(ggplot()+theme_void()+
                                    theme(plot.background=element_rect(fill="gray95",color=NA),
                                          plot.margin=margin(2,2,2,2)))
    n   <- nrow(mat_norm)
    mid <- ceiling(n/2)
    score <- mat_norm[mid, mid]  # fold-over-background at loop position
    
    df <- expand.grid(row=1:n, col=1:n) %>%
      dplyr::mutate(
        value = as.vector(mat_norm),
        x     = col - mid,
        y     = -(row - mid)
      )
    
    p <- ggplot(df, aes(x=x, y=y,
                        fill=pmin(pmax(value, global_min), global_max))) +
      geom_tile() +
      make_fill_scale(palette,
                      limits  = c(global_min, global_max),
                      name    = "Fold over\nbackground") +
      annotate("rect",
               xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5,
               fill=NA, color="black", linewidth=0.6) +
      coord_fixed() +
      scale_x_continuous(
        breaks = c(-BUFFER, 0, BUFFER),
        labels = c(paste0("-",BUFFER*RESOLUTION/1000,"kb"),
                   "loop",
                   paste0("+",BUFFER*RESOLUTION/1000,"kb"))
      ) +
      scale_y_continuous(
        breaks = c(-BUFFER, 0, BUFFER),
        labels = c(paste0("-",BUFFER*RESOLUTION/1000,"kb"),
                   "loop",
                   paste0("+",BUFFER*RESOLUTION/1000,"kb"))
      ) +
      labs(x=NULL, y=NULL) +
      theme_classic(base_size=8) +
      theme(axis.text=element_text(size=6),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            legend.position="none",
            plot.margin=margin(2,2,2,2))
    
    if (!is.na(score))
      p <- p + annotate("text",
                        x=-BUFFER+0.5, y=BUFFER-0.5,
                        label=paste0("APA=", round(score, 2)),
                        size=2.5, hjust=0, fontface="bold", color="white")
    p
  }
  
  col_headers <- lapply(bin_levels, function(b) {
    ggplot() +
      annotate("text", x=0.5, y=0.5, label=b,
               color=bin_colors[[b]], size=3.5, fontface="bold") +
      theme_void() +
      theme(plot.margin=margin(0,2,0,2))
  })
  
  plot_list <- c(list(ggplot()+theme_void()), col_headers)
  for (n in names(CONDITIONS)) {
    row_panels <- lapply(bin_levels, function(b) make_panel(apa_norm[[n]][[b]]))
    label_p <- ggplot() +
      annotate("text", x=0.5, y=0.5,
               label=CONDITIONS[[n]]$label,
               color=CONDITIONS[[n]]$color,
               size=3.5, fontface="bold", angle=90) +
      theme_void()
    plot_list <- c(plot_list, list(label_p), row_panels)
  }
  
  main_grid <- wrap_plots(
    plot_list,
    ncol    = length(bin_levels) + 1,
    widths  = c(0.12, rep(1, length(bin_levels))),
    heights = c(0.15, rep(1, length(CONDITIONS)))
  )
  
  # Legend using normalized scale
  legend_df <- data.frame(
    x    = 1:10, y = 1:10,
    fill = seq(global_min, global_max, length.out=10)
  )
  legend_plot <- ggplot(legend_df, aes(x=x, y=y, fill=fill)) +
    geom_tile() +
    make_fill_scale(palette,
                    limits = c(global_min, global_max),
                    name   = "Fold over\nbackground") +
    theme_void() +
    theme(legend.position    = "right",
          legend.key.height  = unit(1.5,"cm"),
          legend.key.width   = unit(0.3,"cm"),
          legend.text        = element_text(size=8),
          legend.title       = element_text(size=9))
  
  composite <- cowplot::plot_grid(
    main_grid, cowplot::get_legend(legend_plot),
    ncol=2, rel_widths=c(1, 0.07)
  )
  
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=13,
                        x=0, hjust=0, y=0.75) +
    cowplot::draw_label(
      paste0("Background-normalised (fold over corner mean)  |  ",
             RESOLUTION/1000, "kb  |  ±", BUFFER, " bins  |  ",
             "APA = fold-enrichment at loop pixel  |  palette: ", palette),
      size=8, color="gray40", x=0, hjust=0, y=0.25
    )
  
  final <- cowplot::plot_grid(title_grob, composite,
                              ncol=1, rel_heights=c(0.08,1))
  print(final)
  invisible(final)
}

# ══════════════════════════════════════════════════════════════════════════════
# classify_anchors_narrowpeak
# ChIP-seq + ATAC-seq based anchor classifier
# Add to EP_analysis_functions.R and source() from analysis scripts
#
# Rules:
#   1. H3K4me3 overlap → Promoter (ATAC not required)
#   2. ATAC + H3K27ac (no H3K4me3), within prox_kb of any H3K4me3 peak
#      → Proximal-Enhancer
#   3. ATAC + H3K27ac (no H3K4me3), beyond prox_kb → Distal-Enhancer
#   4. ATAC only (no H3K4me3, no H3K27ac) → Open-Chromatin
#   5. None → Unannotated
#
# Disambiguation: when H3K4me3 and H3K27ac both overlap the same anchor,
# the peak whose center is closest to the anchor midpoint wins.
# ATAC filter only applied to H3K27ac calls.
#
# CTCF: annotated independently as anchor1_ctcf / anchor2_ctcf (TRUE/FALSE)
# ══════════════════════════════════════════════════════════════════════════════

load_narrowpeak <- function(path) {
  if (!file.exists(path)) stop("narrowPeak not found: ", path)
  df <- read.table(path, header=FALSE, sep="\t",
                   stringsAsFactors=FALSE, fill=TRUE,
                   col.names=c("chr","start","end","name","score",
                               "strand","signalValue","pValue","qValue","peak"))
  GRanges(seqnames=df$chr, ranges=IRanges(df$start+1L, df$end),
          name=df$name, score=df$score)
}

# ══════════════════════════════════════════════════════════════════════════════
# classify_anchors_narrowpeak (v2 — TSS-geography first)
#
# Geography (ChIPseeker) gates classification; marks confirm type.
# An H3K4me3 peak 500kb from any gene is NOT a Promoter.
#
# Classification logic:
#
#  Geographic context (from ChIPseeker::annotatePeak):
#    PROMOTER_ZONE  = Promoter (default ±3kb of TSS)
#    GENIC_ZONE     = 5'UTR, 3'UTR, Exon, Intron, Downstream
#    INTERGENIC     = Distal Intergenic
#
#  Mark rules applied within each zone:
#
#    PROMOTER_ZONE + H3K4me3                    → Promoter
#    PROMOTER_ZONE + enhancer marks (no K4me3)  → Proximal-Enhancer
#    PROMOTER_ZONE + ATAC only                  → Open-Chromatin
#
#    GENIC_ZONE    + enhancer marks             → Proximal-Enhancer
#    GENIC_ZONE    + ATAC only                  → Open-Chromatin
#
#    INTERGENIC    + enhancer marks             → Distal-Enhancer
#    INTERGENIC    + ATAC only                  → Open-Chromatin
#
#    Any zone, no ATAC, no marks               → Unannotated
#
#  Enhancer marks = (H3K27ac + ATAC) OR (H3K4me1 + ATAC)
#  CTCF = independent boolean columns (anchor1_ctcf / anchor2_ctcf)
#
# Required packages: ChIPseeker, GenomicFeatures
#
# Building a TxDb for hg38 (GENCODE):
#   library(GenomicFeatures)
#   txdb <- makeTxDbFromGFF("gencode.v45.annotation.gtf.gz", format="gtf")
#   saveDb(txdb, "gencode_v45_hg38.sqlite")   # save so you don't rebuild each time
#   txdb <- loadDb("gencode_v45_hg38.sqlite")  # reload in future sessions
#
# Or from UCSC (faster, no GTF needed):
#   library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#   txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# ══════════════════════════════════════════════════════════════════════════════

classify_anchors_narrowpeak <- function(loops_df,
                                        atac_gr,
                                        h3k4me3_gr,
                                        h3k27ac_gr,
                                        txdb,
                                        h3k4me1_gr  = NULL,
                                        ctcf_gr     = NULL,
                                        tss_region  = c(-3000, 3000),
                                        ccre_label  = "",
                                        slop_bp     = 2500L) {
  
  stopifnot(requireNamespace("ChIPseeker", quietly=TRUE))
  
  # ── helpers ──────────────────────────────────────────────────────────────
  fix_gr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if      (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr","",seqlevels(gr))
    suppressWarnings(gr)
  }
  
  make_anchor_windows <- function(chr, start, end) {
    GRanges(chr, IRanges(pmax(1L, start - slop_bp), end + slop_bp))
  }
  
  # Map ChIPseeker annotation string to zone
  chipseeker_to_zone <- function(annotation) {
    dplyr::case_when(
      grepl("^Promoter",   annotation) ~ "PROMOTER",
      grepl("UTR|Exon|Intron|Downstream", annotation) ~ "GENIC",
      TRUE                             ~ "INTERGENIC"
    )
  }
  
  # Run ChIPseeker on an anchor GRanges, return zone per anchor
  get_geographic_zone <- function(anchor_gr, txdb, tss_region) {
    # ChIPseeker needs standard seqlevels
    suppressWarnings({
      anno <- ChIPseeker::annotatePeak(
        anchor_gr,
        tssRegion             = tss_region,
        TxDb                  = txdb,
        level                 = "transcript",
        assignGenomicAnnotation = TRUE,
        genomicAnnotationPriority = c("Promoter","5UTR","3UTR",
                                      "Exon","Intron","Downstream","Intergenic"),
        verbose               = FALSE
      )
    })
    chipseeker_to_zone(as.data.frame(anno)$annotation)
  }
  
  # Classify a single set of anchors given zones + mark overlaps
  classify_one <- function(anchor_gr, midpoints, zone,
                           atac_f, h3k4me3_f, h3k27ac_f, h3k4me1_f) {
    n <- length(anchor_gr)
    
    has_atac    <- overlapsAny(anchor_gr, atac_f)
    has_h3k4me3 <- overlapsAny(anchor_gr, h3k4me3_f)
    
    # H3K27ac — pick closest peak center to anchor midpoint (disambiguation)
    hits_ac     <- findOverlaps(anchor_gr, h3k27ac_f)
    has_h3k27ac <- rep(FALSE, n)
    if (length(hits_ac) > 0) {
      ac_df <- data.frame(
        anchor_idx  = queryHits(hits_ac),
        peak_center = (start(h3k27ac_f)[subjectHits(hits_ac)] +
                         end(h3k27ac_f)[subjectHits(hits_ac)]) / 2
      ) %>%
        dplyr::mutate(dist = abs(peak_center - midpoints[anchor_idx])) %>%
        dplyr::group_by(anchor_idx) %>%
        dplyr::slice_min(dist, n=1, with_ties=FALSE) %>%
        dplyr::ungroup()
      has_h3k27ac[ac_df$anchor_idx] <- TRUE
    }
    
    has_h3k4me1 <- if (!is.null(h3k4me1_f)) overlapsAny(anchor_gr, h3k4me1_f) else rep(FALSE, n)
    
    # Enhancer marks: H3K27ac+ATAC or H3K4me1+ATAC
    has_enhancer_marks <- (has_atac & has_h3k27ac) | (has_atac & has_h3k4me1)
    
    dplyr::case_when(
      # ── Promoter zone ──────────────────────────────────────────────────────
      # TSS-proximal + H3K4me3 → Promoter
      zone == "PROMOTER" & has_h3k4me3
      ~ "Promoter",
      
      # TSS-proximal + enhancer marks (no H3K4me3) → Proximal-Enhancer
      # (active regulatory element at promoter but lacking K4me3 — e.g. poised TSS)
      zone == "PROMOTER" & has_enhancer_marks
      ~ "Proximal-Enhancer",
      
      # TSS-proximal + ATAC only → Open-Chromatin (accessible but unmarked)
      zone == "PROMOTER" & has_atac
      ~ "Open-Chromatin",
      
      # ── Genic zone (UTR / Exon / Intron / Downstream) ────────────────────
      # Genic + enhancer marks → Proximal-Enhancer (intronic/genic enhancer)
      zone == "GENIC" & has_enhancer_marks
      ~ "Proximal-Enhancer",
      
      # Genic + ATAC only → Open-Chromatin
      zone == "GENIC" & has_atac
      ~ "Open-Chromatin",
      
      # ── Intergenic zone ────────────────────────────────────────────────────
      # Intergenic + enhancer marks → Distal-Enhancer
      zone == "INTERGENIC" & has_enhancer_marks
      ~ "Distal-Enhancer",
      
      # Intergenic + ATAC only → Open-Chromatin
      zone == "INTERGENIC" & has_atac
      ~ "Open-Chromatin",
      
      # ── Fallback ───────────────────────────────────────────────────────────
      TRUE ~ "Unannotated"
    )
  }
  
  # ── Build anchor windows ─────────────────────────────────────────────────
  anchor1_gr <- make_anchor_windows(loops_df$chr1, loops_df$start1, loops_df$end1)
  anchor2_gr <- make_anchor_windows(loops_df$chr2, loops_df$start2, loops_df$end2)
  mid1 <- (loops_df$start1 + loops_df$end1) / 2
  mid2 <- (loops_df$start2 + loops_df$end2) / 2
  
  # ── Fix chr naming ────────────────────────────────────────────────────────
  atac_f    <- fix_gr(atac_gr,    anchor1_gr)
  h3k4me3_f <- fix_gr(h3k4me3_gr, anchor1_gr)
  h3k27ac_f <- fix_gr(h3k27ac_gr, anchor1_gr)
  h3k4me1_f <- if (!is.null(h3k4me1_gr)) fix_gr(h3k4me1_gr, anchor1_gr) else NULL
  
  # ── Geographic annotation via ChIPseeker ─────────────────────────────────
  message("classify_anchors_narrowpeak: ", nrow(loops_df), " loops  |  ", ccre_label)
  message("  slop=", slop_bp, "bp  |  TSS region: [",
          tss_region[1], ", ", tss_region[2], "]")
  message("  Running ChIPseeker::annotatePeak on anchor1...")
  zone1 <- get_geographic_zone(anchor1_gr, txdb, tss_region)
  message("  Running ChIPseeker::annotatePeak on anchor2...")
  zone2 <- get_geographic_zone(anchor2_gr, txdb, tss_region)
  message("  Zone1: "); print(table(zone1))
  message("  Zone2: "); print(table(zone2))
  
  # ── Classify ──────────────────────────────────────────────────────────────
  a1 <- classify_one(anchor1_gr, mid1, zone1, atac_f, h3k4me3_f, h3k27ac_f, h3k4me1_f)
  a2 <- classify_one(anchor2_gr, mid2, zone2, atac_f, h3k4me3_f, h3k27ac_f, h3k4me1_f)
  message("  Anchor1 type: "); print(table(a1))
  message("  Anchor2 type: "); print(table(a2))
  
  # ── CTCF (independent) ────────────────────────────────────────────────────
  ctcf1 <- ctcf2 <- rep(FALSE, nrow(loops_df))
  if (!is.null(ctcf_gr)) {
    ctcf_f <- fix_gr(ctcf_gr, anchor1_gr)
    ctcf1  <- overlapsAny(anchor1_gr, ctcf_f)
    ctcf2  <- overlapsAny(anchor2_gr, ctcf_f)
    message("  CTCF anchor1: ", sum(ctcf1), "  anchor2: ", sum(ctcf2))
  }
  
  # ── Pair type ─────────────────────────────────────────────────────────────
  PAIR_ORDER <- c("CTCF","Promoter","Proximal-Enhancer",
                  "Distal-Enhancer","Open-Chromatin","Other","Unannotated")
  make_pair <- function(a, b) {
    a_rank <- match(a, PAIR_ORDER, nomatch=99L)
    b_rank <- match(b, PAIR_ORDER, nomatch=99L)
    if (a_rank <= b_rank) paste(a, b, sep=" — ") else paste(b, a, sep=" — ")
  }
  
  loops_df %>%
    dplyr::mutate(
      anchor1_ccre   = a1,
      anchor2_ccre   = a2,
      loop_ccre_type = mapply(make_pair, a1, a2, USE.NAMES=FALSE),
      anchor1_ctcf   = ctcf1,
      anchor2_ctcf   = ctcf2,
      ccre_source    = ccre_label,
      anchor1_zone   = zone1,    # ChIPseeker geographic zone (diagnostic)
      anchor2_zone   = zone2
    )
}
message("EP_analysis_functions.R loaded: classify_anchors_priority(), classify_anchors_midpoint(), plot_apa_grid_generic(), characterize_anchors() ")