##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Hi-C loop + contact map visualization pipeline -- GENOVA-style
# contact pyramids + gene/cCRE track, ported from the ChIA-PET
# locus plotter (chia_loop_analysis/updated_CHIA_PET_loci_plotter.R).
#
# Key components:
#   1. make_pyramid_heatmap() -- mimics GENOVA::pyramid()/
#      hic_matrixplot() defaults: white -> orange -> red -> black
#      continuous scale (GENOVA's "whitered" palette); pyramid
#      cropped to a fraction of its full height (dist_thresh) to
#      show local structure; single shared Mb-formatted x-axis on
#      the bottom-most panel; slim off-area colour key.
#   2. make_gene_ccre_panel() -- gene-model + ENCODE cCRE panel
#      (build_gene_annotation(), lane packing, UTRs, strand arrows,
#      TSS/TES marks, cCRE colour blocks), reused as-is from the
#      ChIA-PET script so both scripts render loci identically.
#   3. Composite assembly via cowplot::plot_grid, matching how the
#      ChIA-PET script assembles panels + a single shared legend column.
#
# Figure output: Fig. 6A-B (probable).
##################################################################

library(strawr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(scales)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════

# ── Downsampling group switch ──────────────────────────────────────────────────
# "ds_wt_pol2" — WT + Pol2 matched to Pol2 depth (~971M reads)
# "ds_all"     — WT + Pol2 + ActD matched to ActD depth (~712M reads)
# "unsampled"  — undownsampled originals; depth scaling applied
DS_GROUP <- "ds_wt_pol2"   # ← change this to switch groups

BASE     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
LOOP_DIR <- file.path(BASE, "loops")

if (DS_GROUP == "ds_wt_pol2") {
  hic_wt_loops  <- file.path(LOOP_DIR, "ds_wt_pol2/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_deg_loops <- file.path(LOOP_DIR, "ds_wt_pol2/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file   <- file.path(BASE, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file <- file.path(BASE, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT  <- 1.0
  DEPTH_SCALE_DEG <- 1.0
  MIN_OBS_WT      <- 0
  MIN_OBS_DEG     <- 0
  
} else if (DS_GROUP == "ds_all") {
  hic_wt_loops  <- file.path(LOOP_DIR, "ds_all/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_deg_loops <- file.path(LOOP_DIR, "ds_all/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file   <- file.path(BASE, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file <- file.path(BASE, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT  <- 1.0
  DEPTH_SCALE_DEG <- 1.0
  MIN_OBS_WT      <- 0
  MIN_OBS_DEG     <- 0
  
} else if (DS_GROUP == "unsampled") {
  hic_wt_loops  <- file.path(LOOP_DIR, "unsampled/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_deg_loops <- file.path(LOOP_DIR, "unsampled/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file   <- file.path(BASE, "mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file <- file.path(BASE, "mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  wt_total        <- 639409109 + 377220354
  pol2_total      <- 374163054 + 249677343
  DEPTH_SCALE_WT  <- 1.0
  DEPTH_SCALE_DEG <- wt_total / pol2_total
  MIN_OBS_WT      <- 10
  MIN_OBS_DEG     <- 25
  
} else {
  stop("Unknown DS_GROUP: '", DS_GROUP, "'. Must be 'ds_wt_pol2', 'ds_all', or 'unsampled'.")
}

# Depth scaling only makes sense for the unsampled group — ds_wt_pol2 and
# ds_all are already read-depth-matched, so their .hic files should never be
# rescaled. This catches a future edit (or a new DS_GROUP branch) that
# accidentally sets a non-1 factor outside "unsampled".
if (DS_GROUP != "unsampled" && (DEPTH_SCALE_WT != 1.0 || DEPTH_SCALE_DEG != 1.0)) {
  stop("DEPTH_SCALE_WT/DEG must be 1.0 for DS_GROUP = '", DS_GROUP,
       "' — depth scaling is only valid for the 'unsampled' group. ",
       "Got WT = ", DEPTH_SCALE_WT, ", Pol2 = ", DEPTH_SCALE_DEG, ".")
}

message("DS_GROUP: ", DS_GROUP)
message("Depth scale WT: ", DEPTH_SCALE_WT, "  Pol2: ", round(DEPTH_SCALE_DEG, 3))

ENCODE_CCRE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/encode/wt_hct116.bed"

# GENOVA's default "whitered" ramp: white -> orange -> red -> black
GENOVA_WHITERED <- c("white", "orange", "red", "black")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Bézier arc helper
# ══════════════════════════════════════════════════════════════════════════════

bezier_arc <- function(x1, x2, height, n = 100) {
  t  <- seq(0, 1, length.out = n)
  cx <- (x1 + x2) / 2
  cy <- height
  x  <- (1 - t)^2 * x1 + 2 * (1 - t) * t * cx + t^2 * x2
  y  <- (1 - t)^2 * 0  + 2 * (1 - t) * t * cy + t^2 * 0
  data.frame(x = x, y = y)
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Chromosome name normalisation
# ══════════════════════════════════════════════════════════════════════════════

match_chrom_style <- function(chrom, file_chroms) {
  if (chrom %in% file_chroms) return(chrom)
  stripped <- gsub("^chr", "", chrom)
  if (stripped %in% file_chroms) return(stripped)
  prefixed <- paste0("chr", chrom)
  if (prefixed %in% file_chroms) return(prefixed)
  stop("Cannot match chromosome '", chrom,
       "' to file chromosomes: ", paste(head(file_chroms), collapse = ", "))
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Load HiCCUPS loop calls
# ══════════════════════════════════════════════════════════════════════════════

load_hiccups_bedpe <- function(path) {
  all_lines <- readLines(path)
  
  header_line <- all_lines[grepl("^#chr1", all_lines)][1]
  col_names   <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  
  data_lines <- all_lines[!grepl("^#", all_lines)]
  if (length(data_lines) == 0) stop("No data lines found in ", path)
  
  df <- read.table(text = paste(data_lines, collapse = "\n"),
                   header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE)
  
  if (ncol(df) == length(col_names)) colnames(df) <- col_names
  
  out <- data.frame(
    chr1     = as.character(df[[1]]),
    start1   = as.integer(df[[2]]),
    end1     = as.integer(df[[3]]),
    chr2     = as.character(df[[4]]),
    start2   = as.integer(df[[5]]),
    end2     = as.integer(df[[6]]),
    observed = as.numeric(df[["observed"]]),
    fdrDonut = as.numeric(df[["fdrDonut"]]),
    fdrBL    = as.numeric(df[["fdrBL"]]),
    fdrH     = as.numeric(df[["fdrH"]]),
    fdrV     = as.numeric(df[["fdrV"]])
  ) %>%
    filter(!is.na(start1), !is.na(start2))
  
  message("Loaded ", nrow(out), " loops from ", basename(path))
  out
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Load Hi-C contact matrix via strawr
# ══════════════════════════════════════════════════════════════════════════════

load_hic_matrix <- function(hic_path, chrom, start, end,
                            resolution = 10000, norm = "SCALE",
                            depth_scale = 1.0) {
  available_chroms <- strawr::readHicChroms(hic_path)$name
  file_chrom <- match_chrom_style(chrom, available_chroms)
  
  raw <- tryCatch(
    strawr::straw(norm, hic_path,
                  paste0(file_chrom, ":", start, ":", end),
                  paste0(file_chrom, ":", start, ":", end),
                  "BP", resolution),
    error = function(e) {
      stop("straw() failed for '", hic_path, "' chrom ", file_chrom,
           " at resolution ", resolution, " norm ", norm,
           "\nError: ", e$message,
           "\nIf this is a missing-normalization-vector error, that .hic file ",
           "doesn't have '", norm, "' vectors baked in at this resolution — ",
           "check with strawr::readHicNormTypes(hic_path) and either rerun with ",
           "an available norm or regenerate the .hic with that norm added.")
    }
  )
  
  if (nrow(raw) == 0) {
    message("load_hic_matrix('", hic_path, "'): straw() returned 0 rows for ",
            file_chrom, ":", start, "-", end, " at resolution ", resolution,
            " norm ", norm, " — returning NULL. (This is the root cause if a ",
            "downstream panel renders blank.)")
    return(NULL)
  }
  message("load_hic_matrix('", hic_path, "'): retrieved ", nrow(raw), " bins for ",
          file_chrom, ":", start, "-", end, " at resolution ", resolution)
  
  if (depth_scale != 1.0) raw$counts <- raw$counts * depth_scale
  
  upper <- raw
  lower <- raw %>%
    filter(x != y) %>%
    mutate(tmp = x, x = y, y = tmp) %>%
    dplyr::select(-tmp)
  
  bind_rows(upper, lower) %>% arrange(x, y)
}

subtract_hic_matrices <- function(mat_wt, mat_deg) {
  if (is.null(mat_wt) || is.null(mat_deg)) return(NULL)
  full_join(
    mat_wt  %>% dplyr::rename(counts_wt  = counts),
    mat_deg %>% dplyr::rename(counts_deg = counts),
    by = c("x", "y")
  ) %>%
    mutate(
      counts_wt  = replace_na(counts_wt,  0),
      counts_deg = replace_na(counts_deg, 0),
      diff       = counts_deg - counts_wt
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. GENOVA-style pyramid heatmap panel
#    Mirrors GENOVA::pyramid()/hic_matrixplot(): whitered ramp, cropped triangle
#    height, minimal frame, Mb axis only drawn when show_axis = TRUE.
# ══════════════════════════════════════════════════════════════════════════════

make_pyramid_heatmap <- function(mat_df, roi_start, roi_end,
                                 resolution, title = "",
                                 colour_lim = NULL,
                                 diverging = FALSE,
                                 dist_thresh = 1.0,
                                 colour_scale_thresh = NULL,
                                 colours = GENOVA_WHITERED,
                                 show_axis = FALSE,
                                 highlight_loops = NULL,
                                 highlight_color = "black") {
  # colour_lim mirrors GENOVA's hic_matrixplot()/pyramid() intensity switch:
  #   NULL          -> auto-scaled to the 98th percentile of the plotted region
  #                    (GENOVA does this internally when colour_lim is unset)
  #   single number -> treated as the upper bound (c(0, x) for a normal map,
  #                    c(-x, x) for a diverging difference map) — equivalent
  #                    to the old cut.off = x usage
  #   c(min, max)   -> used exactly as given, same as colour_lim = c(min, max)
  #                    in GENOVA
  # dist_thresh: fraction of the FULL pyramid height to draw (1.0 = the whole
  #   triangle, so no loop is ever visually clipped no matter how far apart
  #   its anchors are; < 1.0 crops the apex off if you want to zoom in on
  #   short-range structure instead).
  # colour_scale_thresh: fraction of the full pyramid height used ONLY for
  #   the auto colour_lim quantile calculation (when colour_lim = NULL).
  #   Defaults to dist_thresh (i.e. the auto-scale is computed from whatever
  #   is displayed — the old behavior). Hi-C signal decays sharply with
  #   distance, so at dist_thresh = 1.0 (full triangle) that quantile gets
  #   diluted by a large area of near-zero long-range pixels, washing out
  #   contrast where the real short-range differences live — pass a smaller
  #   colour_scale_thresh (e.g. 0.3) to keep the full triangle on screen
  #   while scaling colour off just the short/mid-range portion.
  # highlight_loops: optional data frame with x_centroid/y_centroid columns
  #   (genomic coordinates of each anchor's midpoint — same columns produced
  #   by match_loop_uniqueness()/classify_loops()). Each row gets an open
  #   circle drawn at its actual matrix pixel, using the same rx/ry rotation
  #   as the heatmap tiles — not just projected onto the x-axis like the arc
  #   panel below the heatmap.
  
  empty_panel <- function(reason) {
    message("make_pyramid_heatmap('", title, "'): rendering blank panel — ", reason)
    ggplot() + theme_void() +
      labs(title = title) +
      theme(plot.title = element_text(size = 12, face = "bold"))
  }
  
  if (is.null(mat_df)) return(empty_panel("mat_df is NULL (upstream load_hic_matrix() returned no data)"))
  if (nrow(mat_df) == 0) return(empty_panel("mat_df has 0 rows"))
  
  value_col <- if (diverging) "diff" else "counts"
  
  mat_roi <- mat_df %>%
    filter(x >= roi_start, x <= roi_end,
           y >= roi_start, y <= roi_end) %>%
    mutate(
      rx = (x + y) / 2,
      ry = (y - x) / 2
    ) %>%
    filter(ry >= 0)
  
  if (nrow(mat_roi) == 0) {
    return(empty_panel(paste0("0 of ", nrow(mat_df), " mat_df rows fall inside roi_start/roi_end (",
                              roi_start, "-", roi_end, ") — check the .hic file's bin ",
                              "coordinates actually cover this region at this resolution")))
  }
  
  # Crop the pyramid to a fraction of its full height. dist_thresh = 1.0
  # (default) draws the entire triangle up to its apex, so long-range loops
  # never get clipped; set it below 1.0 only if you deliberately want to
  # zoom in on short-range contacts and are OK cutting off longer loops.
  full_height <- (roi_end - roi_start) / 2
  y_cap       <- full_height * dist_thresh
  mat_roi     <- mat_roi %>% filter(ry <= y_cap)
  
  if (nrow(mat_roi) == 0) {
    return(empty_panel(paste0("dist_thresh = ", dist_thresh, " cropped every remaining row ",
                              "(y_cap = ", round(y_cap), ") — raise dist_thresh")))
  }
  
  # Same rx/ry rotation as the tiles, so a highlighted loop's circle lands on
  # its actual pixel in the matrix rather than an x-axis projection.
  highlight_pts <- NULL
  if (!is.null(highlight_loops) && nrow(highlight_loops) > 0) {
    highlight_pts <- highlight_loops %>%
      dplyr::mutate(rx = (x_centroid + y_centroid) / 2,
                    ry = (y_centroid - x_centroid) / 2) %>%
      filter(rx >= roi_start, rx <= roi_end, ry >= 0, ry <= y_cap)
  }
  
  # Subset used ONLY for the auto colour_lim quantile — separate from the
  # display crop above, so a full-triangle display (dist_thresh = 1.0) isn't
  # forced to also auto-scale off the whole triangle's long-range near-zero
  # noise. Defaults to dist_thresh, i.e. identical to mat_roi, when
  # colour_scale_thresh isn't set (old behavior, unchanged).
  if (is.null(colour_scale_thresh)) colour_scale_thresh <- dist_thresh
  scale_subset <- if (colour_scale_thresh < dist_thresh) {
    mat_roi %>% filter(ry <= full_height * colour_scale_thresh)
  } else {
    mat_roi
  }
  if (nrow(scale_subset) == 0) scale_subset <- mat_roi
  
  if (diverging) {
    if (is.null(colour_lim)) {
      hi <- quantile(abs(scale_subset$diff), 0.98, na.rm = TRUE)
      lo <- -hi
      message("make_pyramid_heatmap('", title, "'): auto colour_lim = c(",
              round(lo, 2), ", ", round(hi, 2), ") from ", nrow(scale_subset),
              " rows (colour_scale_thresh = ", colour_scale_thresh, ", dist_thresh = ",
              dist_thresh, ")")
    } else if (length(colour_lim) == 1) {
      hi <- abs(colour_lim); lo <- -hi
    } else {
      lo <- colour_lim[1]; hi <- colour_lim[2]
    }
    fill_scale <- scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0,
      limits = c(lo, hi), name = "\u0394 contacts", na.value = "white"
    )
    fill_vals <- pmax(pmin(mat_roi$diff, hi), lo)
  } else {
    if (is.null(colour_lim)) {
      lo <- 0
      hi <- quantile(scale_subset$counts, 0.98, na.rm = TRUE)
    } else if (length(colour_lim) == 1) {
      lo <- 0; hi <- colour_lim
    } else {
      lo <- colour_lim[1]; hi <- colour_lim[2]
    }
    fill_scale <- scale_fill_gradientn(
      colours = colours, limits = c(lo, hi),
      name = "Contacts", na.value = colours[1]
    )
    fill_vals <- pmin(pmax(mat_roi$counts, lo), hi)
  }
  
  x_axis <- if (show_axis) {
    scale_x_continuous(
      limits = c(roi_start, roi_end), expand = c(0, 0),
      labels = label_number(scale = 1e-6, suffix = " Mb", accuracy = 0.01),
      name = "Genomic position"
    )
  } else {
    scale_x_continuous(limits = c(roi_start, roi_end), expand = c(0, 0),
                       labels = NULL, name = NULL)
  }
  
  p <- ggplot(mat_roi) +
    geom_tile(aes(x = rx, y = ry, width = resolution, height = resolution,
                  fill = fill_vals)) +
    fill_scale
  
  if (!is.null(highlight_pts) && nrow(highlight_pts) > 0) {
    p <- p + geom_point(data = highlight_pts, aes(x = rx, y = ry),
                        shape = 21, size = 4, stroke = 1.1,
                        fill = NA, color = highlight_color)
  }
  
  p +
    x_axis +
    scale_y_continuous(limits = c(0, y_cap), expand = c(0, 0)) +
    coord_cartesian(xlim = c(roi_start, roi_end)) +
    labs(title = title, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line    = element_blank(),
      axis.text.x  = if (show_axis) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (show_axis) element_line(linewidth = 0.3) else element_blank(),
      plot.title   = element_text(size = 12, face = "bold"),
      plot.margin  = margin(4, 4, if (show_axis) 4 else 0, 4),
      legend.position   = "right",
      legend.key.width  = unit(0.3, "cm"),
      legend.key.height = unit(0.6, "cm"),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Loop classification (WT vs. degraded) — unchanged logic, kept for the
#    optional difference column
# ══════════════════════════════════════════════════════════════════════════════

classify_loops <- function(loops_wt, loops_deg, chrom,
                           roi_start, roi_end,
                           min_observed_wt, min_observed_deg, max_fdr, tol = 25000) {
  
  filter_loops <- function(df, min_observed) {
    file_chrom <- tryCatch(match_chrom_style(chrom, unique(df$chr1)),
                           error = function(e) NULL)
    if (is.null(file_chrom)) return(data.frame())
    df %>%
      filter(chr1 == file_chrom,
             start1 <= roi_end & end1 >= roi_start,
             start2 <= roi_end & end2 >= roi_start,
             observed >= min_observed,
             fdrDonut <= max_fdr, fdrBL <= max_fdr,
             fdrH <= max_fdr, fdrV <= max_fdr)
  }
  
  wt_roi  <- filter_loops(loops_wt,  min_observed_wt)
  deg_roi <- filter_loops(loops_deg, min_observed_deg)
  if (nrow(wt_roi) == 0 && nrow(deg_roi) == 0) return(data.frame())
  
  if (nrow(wt_roi) > 0 && nrow(deg_roi) > 0) {
    library(InteractionSet)
    gi_wt <- GInteractions(
      anchor1 = GRanges(wt_roi$chr1, IRanges(wt_roi$start1, wt_roi$end1)),
      anchor2 = GRanges(wt_roi$chr2, IRanges(wt_roi$start2, wt_roi$end2))
    )
    gi_deg <- GInteractions(
      anchor1 = GRanges(deg_roi$chr1, IRanges(deg_roi$start1, deg_roi$end1)),
      anchor2 = GRanges(deg_roi$chr2, IRanges(deg_roi$start2, deg_roi$end2))
    )
    gi_wt_exp <- GInteractions(
      anchor1 = trim(resize(anchors(gi_wt, "first"),
                            width(anchors(gi_wt, "first")) + 2 * tol, fix = "center")),
      anchor2 = trim(resize(anchors(gi_wt, "second"),
                            width(anchors(gi_wt, "second")) + 2 * tol, fix = "center"))
    )
    hits <- findOverlaps(gi_deg, gi_wt_exp, use.region = "both")
    deg_roi$loop_class <- "degraded_only"
    deg_roi$loop_class[unique(queryHits(hits))] <- "common"
  } else {
    deg_roi$loop_class <- "degraded_only"
  }
  deg_roi %>% dplyr::mutate(x_centroid = (start1 + end1) / 2, y_centroid = (start2 + end2) / 2)
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. Loop matching + arc panel
#    match_loop_uniqueness() filters loops to the ROI/thresholds, flags any
#    loop with no matching anchor pair in the other condition (is_unique),
#    and computes x_centroid/y_centroid (anchor midpoints) once — shared by
#    both the arc panel and the matrix-highlight circles in
#    make_pyramid_heatmap(), so there's a single source of truth for "where
#    is this loop" instead of two panels re-deriving it independently.
# ══════════════════════════════════════════════════════════════════════════════

match_loop_uniqueness <- function(loops_df, other_loops_df, chrom, roi_start, roi_end,
                                  min_observed, other_min_observed, max_fdr,
                                  loop_tol = 25000) {
  
  file_chrom <- tryCatch(match_chrom_style(chrom, unique(loops_df$chr1)),
                         error = function(e) NULL)
  if (is.null(file_chrom)) return(data.frame())
  
  loops_roi <- loops_df %>%
    filter(chr1 == file_chrom,
           start1 <= roi_end & end1 >= roi_start,
           start2 <= roi_end & end2 >= roi_start,
           observed >= min_observed,
           fdrDonut <= max_fdr, fdrBL <= max_fdr,
           fdrH <= max_fdr, fdrV <= max_fdr) %>%
    dplyr::mutate(x_centroid = (start1 + end1) / 2, y_centroid = (start2 + end2) / 2)
  
  if (nrow(loops_roi) == 0) return(loops_roi)
  
  loops_roi$is_unique <- FALSE
  if (!is.null(other_loops_df)) {
    other_chrom <- tryCatch(match_chrom_style(chrom, unique(other_loops_df$chr1)),
                            error = function(e) NULL)
    if (is.null(other_chrom)) {
      loops_roi$is_unique <- TRUE
    } else {
      other_roi <- other_loops_df %>%
        filter(chr1 == other_chrom,
               start1 <= roi_end & end1 >= roi_start,
               start2 <= roi_end & end2 >= roi_start,
               observed >= other_min_observed,
               fdrDonut <= max_fdr, fdrBL <= max_fdr,
               fdrH <= max_fdr, fdrV <= max_fdr)
      if (nrow(other_roi) == 0) {
        loops_roi$is_unique <- TRUE
      } else {
        library(InteractionSet)
        gi_this  <- GInteractions(
          anchor1 = GRanges(loops_roi$chr1, IRanges(loops_roi$start1, loops_roi$end1)),
          anchor2 = GRanges(loops_roi$chr2, IRanges(loops_roi$start2, loops_roi$end2))
        )
        gi_other <- GInteractions(
          anchor1 = GRanges(other_roi$chr1, IRanges(other_roi$start1, other_roi$end1)),
          anchor2 = GRanges(other_roi$chr2, IRanges(other_roi$start2, other_roi$end2))
        )
        gi_other_exp <- GInteractions(
          anchor1 = trim(resize(anchors(gi_other, "first"),
                                width(anchors(gi_other, "first")) + 2 * loop_tol, fix = "center")),
          anchor2 = trim(resize(anchors(gi_other, "second"),
                                width(anchors(gi_other, "second")) + 2 * loop_tol, fix = "center"))
        )
        hits <- findOverlaps(gi_this, gi_other_exp, use.region = "both")
        loops_roi$is_unique <- TRUE
        loops_roi$is_unique[unique(queryHits(hits))] <- FALSE
      }
    }
  }
  loops_roi
}

make_loop_arc_panel <- function(loops_roi, roi_start, roi_end,
                                low_col = "#6baed6", high_col = "#08306b",
                                legend_name = "Observed") {
  # loops_roi: output of match_loop_uniqueness() (or classify_loops() with
  # x_centroid/y_centroid added) — already filtered to the ROI/thresholds.
  
  if (is.null(loops_roi) || nrow(loops_roi) == 0) return(ggplot() + theme_void())
  
  arc_data <- bind_rows(lapply(seq_len(nrow(loops_roi)), function(i) {
    row <- loops_roi[i, ]
    df  <- bezier_arc(row$x_centroid, row$y_centroid, height = 1)
    df$loop_id  <- i
    df$observed <- row$observed
    df
  }))
  
  ggplot() +
    geom_path(data = arc_data, aes(x = x, y = y, group = loop_id, color = observed),
              linewidth = 0.6, alpha = 0.75) +
    scale_color_gradient(low = low_col, high = high_col, name = legend_name) +
    scale_x_continuous(limits = c(roi_start, roi_end), expand = c(0, 0),
                       labels = NULL, name = NULL) +
    # Reversed: arc endpoints sit at y = 0 (the baseline directly under the
    # heatmap), and the panel is flipped so the y = 0.5 apex renders toward
    # the BOTTOM of the panel — loops hang down from the matrix instead of
    # bulging up into it.
    scale_y_reverse(limits = c(1.2, 0), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(0, 4, 0, 4),
      legend.position = "right",
      legend.key.width = unit(0.3, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Gene model + ENCODE cCRE panel
#    Ported directly from the ChIA-PET script's build_gene_annotation() and its
#    Panel-2 rendering block, so both pipelines draw loci identically.
# ══════════════════════════════════════════════════════════════════════════════

CCRE_COLORS_FILL <- c(
  "5UTR"    = "#41ab5d",
  "3UTR"    = "#fd8d3c",
  "PLS"     = "#e31a1c",
  "pELS"    = "#ff7f00",
  "dELS"    = "#fdbf6f",
  "CA-CTCF" = "#1f78b4",
  "CA-TF"   = "#6a3d9a",
  "TF"      = "#b15928"
)
CCRE_LABELS_FILL <- c(
  "5UTR"    = "5' UTR",
  "3UTR"    = "3' UTR",
  "PLS"     = "Promoter (PLS)",
  "pELS"    = "Prox. enhancer",
  "dELS"    = "Dist. enhancer",
  "CA-CTCF" = "CTCF-bound",
  "CA-TF"   = "TF-bound",
  "TF"      = "TF-only"
)

build_gene_annotation <- function(chrom, roi_start, roi_end, txdb, orgdb) {
  
  roi_gr    <- GRanges(chrom, IRanges(roi_start, roi_end))
  all_genes <- suppressMessages(subsetByOverlaps(genes(txdb), roi_gr))
  if (length(all_genes) == 0) return(NULL)
  
  genes_df <- as.data.frame(all_genes, row.names = NULL)
  sym_map  <- suppressMessages(
    AnnotationDbi::select(orgdb, keys = genes_df$gene_id, keytype = "ENTREZID",
                          columns = c("SYMBOL", "GENENAME"))
  )
  genes_df <- left_join(genes_df, sym_map, by = c("gene_id" = "ENTREZID")) %>%
    dplyr::mutate(
      SYMBOL      = ifelse(is.na(SYMBOL), gene_id, SYMBOL),
      is_lncrna   = grepl("long non-coding|lncRNA|antisense|LINC", GENENAME, ignore.case = TRUE),
      gene_type   = ifelse(is_lncrna, "lncRNA", "protein_coding"),
      disp_start  = pmax(start, roi_start),
      disp_end    = pmin(end,   roi_end),
      mid_x       = (disp_start + disp_end) / 2,
      strand_char = as.character(strand)
    )
  
  genes_df  <- genes_df %>% dplyr::arrange(start)
  lane_ends <- c()
  genes_df$lane <- 0L
  for (i in seq_len(nrow(genes_df))) {
    placed <- FALSE
    for (ln in seq_along(lane_ends)) {
      if (genes_df$start[i] > lane_ends[ln] + 5000) {
        genes_df$lane[i] <- ln
        lane_ends[ln]     <- genes_df$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      lane_ends <- c(lane_ends, genes_df$end[i])
      genes_df$lane[i] <- length(lane_ends)
    }
  }
  genes_df$y <- -(genes_df$lane - 1L) * 0.65
  gene_y_map <- setNames(genes_df$y, genes_df$gene_id)
  
  tx_by_gene  <- suppressMessages(transcriptsBy(txdb, by = "gene"))
  exons_by_tx <- suppressMessages(exonsBy(txdb, by = "tx"))
  utr5_by_tx  <- suppressMessages(fiveUTRsByTranscript(txdb))
  utr3_by_tx  <- suppressMessages(threeUTRsByTranscript(txdb))
  tx_lengths  <- suppressMessages(
    transcriptLengths(txdb, with.utr5_len = TRUE, with.utr3_len = TRUE)
  )
  
  exons_df <- data.frame(); introns_df <- data.frame()
  utrs_df  <- data.frame(); arrows_df  <- data.frame()
  chrom_bare <- gsub("^chr", "", chrom)
  
  for (gid in genes_df$gene_id) {
    if (!gid %in% names(tx_by_gene)) next
    gene_row    <- genes_df[genes_df$gene_id == gid, ][1, ]
    gene_left   <- gene_row$start
    gene_right  <- gene_row$end
    gene_strand <- gene_row$strand_char
    gene_type   <- gene_row$gene_type
    y_base      <- gene_y_map[gid]
    
    tx_sub <- tx_lengths[tx_lengths$gene_id == gid, ]
    if (nrow(tx_sub) == 0) next
    canon_name <- tx_sub$tx_name[which.max(tx_sub$tx_len)]
    all_tx_gr  <- tx_by_gene[[gid]]
    tx_id_map  <- setNames(all_tx_gr$tx_id, all_tx_gr$tx_name)
    canon_id   <- tx_id_map[canon_name]
    if (is.na(canon_id)) next
    canon_key  <- as.character(canon_id)
    if (!canon_key %in% names(exons_by_tx)) next
    
    raw_exons <- as.data.frame(exons_by_tx[[canon_key]], row.names = NULL) %>%
      filter(gsub("^chr", "", as.character(seqnames)) == chrom_bare,
             start <= gene_right, end >= gene_left) %>%
      dplyr::mutate(start = pmax(start, gene_left), end = pmin(end, gene_right),
                    gene_id = gid, gene_type = gene_type, y_base = y_base) %>%
      filter(end > start)
    if (nrow(raw_exons) == 0) next
    exons_df <- bind_rows(exons_df, raw_exons)
    
    get_utr_df <- function(utr_list, key, utype) {
      if (!key %in% names(utr_list)) return(NULL)
      as.data.frame(utr_list[[key]], row.names = NULL) %>%
        dplyr::mutate(utr_type = utype, gene_id = gid, y_base = y_base)
    }
    utrs_df <- bind_rows(utrs_df,
                         get_utr_df(utr5_by_tx, canon_key, "5UTR"),
                         get_utr_df(utr3_by_tx, canon_key, "3UTR"))
    
    if (nrow(raw_exons) >= 2) {
      ex_ord  <- raw_exons %>% dplyr::arrange(start)
      introns <- data.frame(
        intron_start = unname(ex_ord$end[-nrow(ex_ord)]),
        intron_end   = unname(ex_ord$start[-1]),
        gene_id = gid, gene_type = gene_type, y_base = y_base,
        row.names = NULL
      ) %>% filter(intron_end > intron_start)
      introns_df <- bind_rows(introns_df, introns)
    }
    
    gene_span <- gene_right - gene_left
    n_arr <- max(1L, floor(gene_span / 60000))
    arr_x <- seq(gene_left + gene_span / (2 * n_arr),
                 gene_right - gene_span / (2 * n_arr), length.out = n_arr)
    
    in_intron_fn <- function(ax) {
      if (nrow(introns_df) == 0) return(FALSE)
      sub <- introns_df[introns_df$gene_id == gid, ]
      any(sub$intron_start <= ax & sub$intron_end >= ax)
    }
    
    arrows_df <- bind_rows(arrows_df, data.frame(
      x = unname(arr_x), strand = gene_strand, gene_id = gid, y_base = y_base,
      in_intron = vapply(arr_x, in_intron_fn, logical(1)),
      row.names = NULL
    ))
  }
  
  list(genes_df = genes_df, exons_df = exons_df, introns_df = introns_df,
       utrs_df = utrs_df, arrows_df = arrows_df)
}

load_ccres <- function(ccre_path, chrom, roi_start, roi_end) {
  tryCatch({
    df <- read.table(ccre_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
      dplyr::select(chr = V1, start = V2, end = V3, id = V4, ccre_type = V10) %>%
      filter(chr == chrom, start <= roi_end, end >= roi_start,
             ccre_type %in% c("PLS", "pELS", "dELS", "CA-CTCF", "CA-TF", "TF"))
    if (nrow(df) == 0) NULL else df
  }, error = function(e) {
    message("  cCRE load error: ", e$message)
    NULL
  })
}

make_gene_ccre_panel <- function(chrom, roi_start, roi_end, txdb, orgdb,
                                 ccre_path = ENCODE_CCRE, show_axis = FALSE) {
  
  ann <- build_gene_annotation(chrom, roi_start, roi_end, txdb, orgdb)
  has_genes   <- !is.null(ann) && nrow(ann$genes_df)   > 0
  has_exons   <- !is.null(ann) && nrow(ann$exons_df)   > 0
  has_introns <- !is.null(ann) && nrow(ann$introns_df) > 0
  has_utrs    <- !is.null(ann) && nrow(ann$utrs_df)    > 0
  has_arrows  <- !is.null(ann) && nrow(ann$arrows_df)  > 0
  n_lanes     <- if (has_genes) max(ann$genes_df$lane) else 1L
  
  ccres     <- load_ccres(ccre_path, chrom, roi_start, roi_end)
  has_ccres <- !is.null(ccres) && nrow(ccres) > 0
  
  EXON_H <- 0.16; UTR_H <- 0.09; CCRE_Y <- 0.65; CCRE_H <- 0.18
  gene_y_top    <- if (has_ccres) CCRE_Y + CCRE_H + 0.15 else EXON_H + 0.38
  gene_y_bottom <- -(n_lanes - 1) * 0.65
  gene_y_min    <- gene_y_bottom - EXON_H - 0.55
  gene_y_max    <- gene_y_top
  
  p <- ggplot()
  
  if (has_ccres) {
    region_span <- roi_end - roi_start
    min_width   <- region_span * 0.003
    ccres_plot  <- ccres %>%
      dplyr::mutate(plot_end = pmax(end, start + min_width),
                    ccre_type = ifelse(ccre_type %in% names(CCRE_COLORS_FILL), ccre_type, "TF"))
    p <- p +
      annotate("segment", x = roi_start, xend = roi_end,
               y = CCRE_Y - CCRE_H - 0.06, yend = CCRE_Y - CCRE_H - 0.06,
               color = "gray88", linewidth = 0.3) +
      geom_rect(data = ccres_plot,
                aes(xmin = start, xmax = plot_end,
                    ymin = CCRE_Y - CCRE_H, ymax = CCRE_Y + CCRE_H, fill = ccre_type),
                color = NA, alpha = 0.9) +
      annotate("text", x = roi_start + (roi_end - roi_start) * 0.005, y = CCRE_Y,
               label = "cCREs", hjust = 0, vjust = 0.5, size = 2.4,
               color = "gray30", fontface = "bold")
  }
  
  if (has_introns) {
    p <- p + geom_segment(data = ann$introns_df,
                          aes(x = intron_start, xend = intron_end,
                              y = y_base, yend = y_base, linetype = gene_type),
                          color = "gray25", linewidth = 0.45)
  }
  
  if (has_arrows) {
    arr <- ann$arrows_df %>%
      dplyr::mutate(dx = ifelse(strand == "+", 8000, -8000),
                    x_from = x - abs(dx) / 2,
                    x_to   = x + abs(dx) / 2 * ifelse(strand == "+", 1, -1))
    p <- p + geom_segment(data = arr, aes(x = x_from, xend = x_to, y = y_base, yend = y_base),
                          color = "gray50", linewidth = 0.4,
                          arrow = arrow(length = unit(0.07, "cm"), type = "closed", ends = "last"))
  }
  
  if (has_utrs) {
    p <- p + geom_rect(data = ann$utrs_df,
                       aes(xmin = start, xmax = end,
                           ymin = y_base - UTR_H, ymax = y_base + UTR_H, fill = utr_type),
                       color = NA)
  }
  
  if (has_exons) {
    exon_colors <- setNames(
      ifelse(ann$genes_df$gene_type == "lncRNA", "#9ecae1", "#252525"),
      ann$genes_df$gene_id
    )
    p <- p + geom_rect(data = ann$exons_df,
                       aes(xmin = start, xmax = end,
                           ymin = y_base - EXON_H, ymax = y_base + EXON_H),
                       fill = exon_colors[ann$exons_df$gene_id], color = NA)
  }
  
  if (has_genes) {
    tss_tes_marks <- ann$genes_df %>%
      dplyr::mutate(tss = ifelse(strand_char == "+", start, end),
                    tes = ifelse(strand_char == "+", end, start)) %>%
      { bind_rows(
        dplyr::transmute(., x = tss, y = y - EXON_H - 0.07, yend = y + EXON_H + 0.07,
                         mark = "TSS"),
        dplyr::transmute(., x = tes, y = y - EXON_H - 0.07, yend = y + EXON_H + 0.07,
                         mark = "TES")
      ) }
    p <- p +
      geom_segment(data = tss_tes_marks,
                   aes(x = x, xend = x, y = y, yend = yend, color = mark,
                       linewidth = mark),
                   lineend = "butt")
  }
  
  if (has_genes) {
    p <- p +
      geom_text(data = ann$genes_df,
                aes(x = mid_x, y = y - EXON_H - 0.32, label = SYMBOL,
                    fontface = ifelse(gene_type == "lncRNA", "italic", "bold.italic")),
                size = 3.0, color = "gray15")
  }
  
  utr_types_present  <- if (has_utrs)  unique(ann$utrs_df$utr_type) else c()
  ccre_types_present <- if (has_ccres) unique(ccres$ccre_type)      else c()
  legend_breaks <- intersect(names(CCRE_COLORS_FILL), c(utr_types_present, ccre_types_present))
  legend_labels <- CCRE_LABELS_FILL[legend_breaks]
  
  x_axis <- if (show_axis) {
    scale_x_continuous(limits = c(roi_start, roi_end), expand = c(0, 0),
                       labels = label_number(scale = 1e-6, suffix = " Mb", accuracy = 0.01),
                       name = "Genomic position")
  } else {
    scale_x_continuous(limits = c(roi_start, roi_end), expand = c(0, 0), labels = NULL, name = NULL)
  }
  
  p +
    scale_linetype_manual(values = c(protein_coding = "solid", lncRNA = "dashed"),
                          name = "Gene type",
                          guide = guide_legend(override.aes = list(linewidth = 0.8))) +
    scale_fill_manual(values = CCRE_COLORS_FILL, breaks = legend_breaks,
                      labels = legend_labels, name = "Annotation", na.translate = FALSE) +
    scale_color_manual(values = c(TSS = "#1a9850", TES = "#d73027"),
                       labels = c(TSS = "TSS (5' end)", TES = "TES (3' end)"),
                       name = "Transcript end",
                       guide = guide_legend(override.aes = list(linewidth = 2))) +
    scale_linewidth_manual(values = c(TSS = 1.1, TES = 0.9), guide = "none") +
    x_axis +
    # Flipped to match the loop panel below it: this panel now sits directly
    # under the heatmap, and gets mirrored vertically (cCRE row/lane 1, which
    # was drawn near the top, now renders near the bottom of this panel,
    # closer to the loop track that hangs beneath it).
    scale_y_reverse(limits = c(gene_y_max, gene_y_min), expand = c(0, 0)) +
    theme_classic(base_size = 10) +
    theme(
      axis.line.y  = element_blank(), axis.ticks.y = element_blank(),
      axis.text.y  = element_blank(), axis.title.y = element_blank(),
      axis.text.x  = if (show_axis) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (show_axis) element_line(linewidth = 0.3) else element_blank(),
      legend.position = "right", legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      plot.margin = margin(0, 4, if (show_axis) 4 else 0, 4)
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Main composite plot — cowplot assembly, GENOVA-style pyramids + shared
#    gene/cCRE track, optional WT/Degraded/Difference columns
# ══════════════════════════════════════════════════════════════════════════════

plot_hic_region <- function(chrom, start, end,
                            hic_wt_path, hic_deg_path,
                            loops_wt_path, loops_deg_path,
                            txdb, orgdb,
                            label = NULL,
                            resolution = 10000, norm = "SCALE",
                            colour_lim = NULL,
                            colour_lim_diff = NULL,
                            colour_scale_thresh_diff = 0.3,
                            min_observed_wt = 10, min_observed_deg = 10,
                            depth_scale_wt = 1.0, depth_scale_deg = 1.0,
                            max_fdr = 0.1,
                            dist_thresh = 1.0,
                            show_difference = FALSE,
                            loop_tol = 25000,
                            ccre_path = ENCODE_CCRE,
                            highlight_lost_loops = TRUE) {
  # colour_lim / colour_lim_diff: same NULL / single-value / c(min,max) switch
  # as make_pyramid_heatmap()'s colour_lim (see that function for details).
  # colour_lim is shared between the Wildtype and Degraded columns so they're
  # directly comparable; colour_lim_diff controls the separate difference
  # panel's diverging scale when show_difference = TRUE.
  # colour_scale_thresh_diff: only affects auto-scaling (colour_lim_diff =
  # NULL) — restricts the auto colour_lim quantile calc to this fraction of
  # the full pyramid height even though the whole triangle (dist_thresh)
  # still gets drawn. Long-range Hi-C signal decays toward zero in both
  # conditions, so without this the auto-scale gets diluted by that
  # low-signal area and the difference map looks washed out. Set to NULL
  # (or equal to dist_thresh) to go back to scaling off the whole triangle.
  # min_observed_wt/deg and depth_scale_wt/deg let each condition carry its
  # own loop-calling threshold and read-depth normalization factor (e.g. the
  # MIN_OBS_WT/MIN_OBS_DEG/DEPTH_SCALE_WT/DEPTH_SCALE_DEG values from the
  # DS_GROUP parameter block).
  # dist_thresh defaults to 1.0 (the full pyramid) so no loop, however far
  # its anchors are apart, gets visually clipped — see make_pyramid_heatmap().
  # highlight_lost_loops: circle any loop's arc apex if it has no matching
  # anchor pair (within loop_tol bp) in the other condition — an open black
  # circle on the WT panel flags a WT-only ("lost") loop, and on the
  # Degraded panel flags a Degraded-only ("gained") loop.
  
  pad       <- (end - start) * 0.02
  roi_start <- start - pad
  roi_end   <- end   + pad
  
  loops_wt  <- load_hiccups_bedpe(loops_wt_path)
  loops_deg <- load_hiccups_bedpe(loops_deg_path)
  
  mat_wt  <- load_hic_matrix(hic_wt_path,  chrom, roi_start, roi_end, resolution, norm,
                             depth_scale = depth_scale_wt)
  mat_deg <- load_hic_matrix(hic_deg_path, chrom, roi_start, roi_end, resolution, norm,
                             depth_scale = depth_scale_deg)
  
  # Auto-scale once, shared across WT/Degraded, only if the caller didn't
  # already pin an explicit colour_lim.
  if (is.null(colour_lim) && !is.null(mat_wt) && !is.null(mat_deg)) {
    colour_lim <- quantile(c(mat_wt$counts, mat_deg$counts), 0.98, na.rm = TRUE)
  }
  
  build_column <- function(mat, loops, title, diverging = FALSE, low_col, high_col,
                           min_observed, other_loops, other_min_observed) {
    loops_matched <- match_loop_uniqueness(loops, other_loops, chrom, roi_start, roi_end,
                                           min_observed, other_min_observed, max_fdr,
                                           loop_tol = loop_tol)
    highlight_df <- if (highlight_lost_loops && nrow(loops_matched) > 0) {
      loops_matched %>% filter(is_unique)
    } else NULL
    
    p_hic   <- make_pyramid_heatmap(mat, roi_start, roi_end, resolution,
                                    title = title, colour_lim = colour_lim,
                                    diverging = diverging, dist_thresh = dist_thresh,
                                    show_axis = FALSE,
                                    highlight_loops = highlight_df)
    p_loops <- make_loop_arc_panel(loops_matched, roi_start, roi_end,
                                   low_col = low_col, high_col = high_col,
                                   legend_name = paste0("Observed (", title, ")"))
    p_genes <- make_gene_ccre_panel(chrom, roi_start, roi_end, txdb, orgdb,
                                    ccre_path = ccre_path, show_axis = TRUE)
    list(hic = p_hic, loops = p_loops, genes = p_genes)
  }
  
  col_wt  <- build_column(mat_wt,  loops_wt,  "Wildtype", FALSE, "#6baed6", "#08306b",
                          min_observed_wt,  other_loops = loops_deg, other_min_observed = min_observed_deg)
  col_deg <- build_column(mat_deg, loops_deg, "Degraded",  FALSE, "#fc9272", "#67000d",
                          min_observed_deg, other_loops = loops_wt,  other_min_observed = min_observed_wt)
  
  columns <- list(col_wt, col_deg)
  # Relative panel heights: heatmap tallest, loop arcs get more room than the
  # gene/cCRE track (matching the ChIA-PET script's 3 : ~1.65 ratio), not less.
  loop_panel_height <- 2.2
  gene_panel_range  <- 1.4
  
  if (show_difference) {
    if (is.null(mat_wt) || is.null(mat_deg)) {
      message("show_difference = TRUE but the ", 
              if (is.null(mat_wt)) "Wildtype " else "",
              if (is.null(mat_deg)) "Degraded " else "",
              "contact matrix returned zero bins for '", chrom, ":", start, "-", end,
              "' at resolution ", resolution, " (see the load_hic_matrix message ",
              "above). Rendering an empty difference column instead ",
              "of failing the whole plot — try a coarser resolution or a wider region.")
      p_diff_hic   <- ggplot() + theme_void() +
        labs(title = "Difference (deg \u2212 wt) \u2014 no data") +
        theme(plot.title = element_text(size = 12, face = "bold"))
      p_diff_loops <- ggplot() + theme_void()
      p_diff_genes <- make_gene_ccre_panel(chrom, roi_start, roi_end, txdb, orgdb,
                                           ccre_path = ccre_path, show_axis = TRUE)
      columns[[3]] <- list(hic = p_diff_hic, loops = p_diff_loops, genes = p_diff_genes)
    } else {
      diff_mat   <- subtract_hic_matrices(mat_wt, mat_deg)
      classified <- classify_loops(loops_wt, loops_deg, chrom, roi_start, roi_end,
                                   min_observed_wt, min_observed_deg, max_fdr, tol = loop_tol)
      p_diff_hic <- make_pyramid_heatmap(diff_mat, roi_start, roi_end, resolution,
                                         title = "Difference (deg \u2212 wt)",
                                         colour_lim = colour_lim_diff,
                                         colour_scale_thresh = colour_scale_thresh_diff,
                                         diverging = TRUE, dist_thresh = dist_thresh,
                                         show_axis = FALSE)
      p_diff_loops <- if (nrow(classified) > 0) {
        make_loop_arc_panel(classified, roi_start, roi_end,
                            low_col = "#2166ac", high_col = "#d6604d",
                            legend_name = "Observed (Degraded loop set)")
      } else ggplot() + theme_void()
      p_diff_genes <- make_gene_ccre_panel(chrom, roi_start, roi_end, txdb, orgdb,
                                           ccre_path = ccre_path, show_axis = TRUE)
      columns[[3]] <- list(hic = p_diff_hic, loops = p_diff_loops, genes = p_diff_genes)
    }
  }
  
  assemble_col <- function(col) {
    plot_grid(
      col$hic   + theme(legend.position = "none"),
      col$genes + theme(legend.position = "none"),
      col$loops + theme(legend.position = "none"),
      ncol = 1, align = "v", axis = "lr",
      rel_heights = c(2.4, gene_panel_range, loop_panel_height)
    )
  }
  
  panel_cols <- lapply(columns, assemble_col)
  
  # Legend column: WT and Degraded arcs use different color scales (see
  # make_loop_arc_panel's low_col/high_col), so a legend pulled from only one
  # column would mislabel — or, if that column happened to have zero loops
  # passing filters in this ROI, omit entirely — the other column's arc
  # colors. Pull from every column and skip any panel that has no legend
  # (empty ROI, or a get_legend/ggplot2 version mismatch) instead of letting
  # one missing legend silently drop from the assembled key.
  safe_legend <- function(p) {
    if (is.null(p)) return(NULL)
    tryCatch(get_legend(p + theme(legend.position = "right")), error = function(e) NULL)
  }
  
  legend_candidates <- list(
    safe_legend(col_wt$hic),
    safe_legend(col_wt$genes)
  )
  if (show_difference && length(columns) >= 3) {
    legend_candidates <- c(legend_candidates, list(safe_legend(columns[[3]]$hic)))
  }
  legend_candidates <- Filter(Negate(is.null), legend_candidates)
  
  legends <- if (length(legend_candidates) > 0) {
    plot_grid(plotlist = legend_candidates, ncol = 1, align = "v")
  } else {
    ggplot() + theme_void()
  }
  
  plot_title <- if (!is.null(label)) label else paste0(chrom, ":",
                                                       format(start, big.mark = ","), "-",
                                                       format(end,   big.mark = ","))
  
  body <- plot_grid(plotlist = c(panel_cols, list(legends)),
                    nrow = 1,
                    rel_widths = c(rep(1, length(panel_cols)), 0.32))
  
  title_grob <- ggdraw() +
    draw_label(plot_title, fontface = "bold", size = 14, x = 0.02, hjust = 0)
  
  final <- plot_grid(title_grob, body, ncol = 1, rel_heights = c(0.06, 1))
  print(final)
  invisible(final)
}



# ══════════════════════════════════════════════════════════════════════════════
# 10. Run
# ══════════════════════════════════════════════════════════════════════════════
interval=500000
plot_hic_region(
  chrom = "chr5", start = 130500953+interval, end = 132800685-interval,
  hic_wt_path = hic_wt_file, hic_deg_path = hic_pol2_file,
  loops_wt_path = hic_wt_loops, loops_deg_path = hic_deg_loops,
  txdb = txdb, orgdb = org.Hs.eg.db,
  label = "chr5 locus (zoomed out)",
  resolution = 10000, norm = "SCALE",
  colour_lim = c(0, 25),   
  colour_lim_diff = c(-10, 10),         # auto for the diverging difference panel
  min_observed_wt = MIN_OBS_WT, min_observed_deg = MIN_OBS_DEG,
  depth_scale_wt = DEPTH_SCALE_WT, depth_scale_deg = DEPTH_SCALE_DEG,
  max_fdr = 0.1,
  dist_thresh = 1.0, show_difference = TRUE, loop_tol = 25000
)

plot_hic_region(
  chrom = "chr5", start = 132618052, end = 132688395,
  hic_wt_path = hic_wt_file, hic_deg_path = hic_pol2_file,
  loops_wt_path = hic_wt_loops, loops_deg_path = hic_deg_loops,
  txdb = txdb, orgdb = org.Hs.eg.db,
  label = "chr5 locus (zoomed in)",
  resolution = 10000, norm = "SCALE",
  colour_lim = 40,                 # pin the WT/Degraded scale to 0–40 for this zoom level
  min_observed_wt = MIN_OBS_WT, min_observed_deg = MIN_OBS_DEG,
  depth_scale_wt = DEPTH_SCALE_WT, depth_scale_deg = DEPTH_SCALE_DEG,
  max_fdr = 0.1,
  dist_thresh = 1.0, show_difference = TRUE, loop_tol = 25000
)

plot_hic_region(
  chrom = "chr19", start = 44506844, end = 45008551,
  hic_wt_path = hic_wt_file, hic_deg_path = hic_pol2_file,
  loops_wt_path = hic_wt_loops, loops_deg_path = hic_deg_loops,
  txdb = txdb, orgdb = org.Hs.eg.db,
  label = "chr19 locus",
  resolution = 5000, norm = "SCALE",
  colour_lim = c(0, 50),            # explicit min/max, GENOVA colour_lim style
  min_observed_wt = MIN_OBS_WT, min_observed_deg = MIN_OBS_DEG,
  depth_scale_wt = DEPTH_SCALE_WT, depth_scale_deg = DEPTH_SCALE_DEG,
  max_fdr = 0.1,
  dist_thresh = 1.0, show_difference = FALSE
)

plot_hic_region(
  chrom = "chr22", start = 32000500, end = 33200000,
  hic_wt_path = hic_wt_file, hic_deg_path = hic_pol2_file,
  loops_wt_path = hic_wt_loops, loops_deg_path = hic_deg_loops,
  txdb = txdb, orgdb = org.Hs.eg.db,
  label = "chr19 locus",
  resolution = 5000, norm = "SCALE",
  colour_lim = c(0, 25),            # explicit min/max, GENOVA colour_lim style
  min_observed_wt = MIN_OBS_WT, min_observed_deg = MIN_OBS_DEG,
  depth_scale_wt = DEPTH_SCALE_WT, depth_scale_deg = DEPTH_SCALE_DEG,
  max_fdr = 0.1,
  dist_thresh = 1.0, show_difference = FALSE
)

#wt only lost	22	32650000	32675000	32662500	22	32975000	33000000	33000000	337500	7.933793544769287

# pol 2 only gained	3.9981333831460706	19	30220000	30230000	30225000	19	30650000	30660000	30655000	430000	0

plot_hic_region(
  chrom = "chr19", start = 30225000-1000000, end = 30655000+1000000,
  hic_wt_path = hic_wt_file, hic_deg_path = hic_pol2_file,
  loops_wt_path = hic_wt_loops, loops_deg_path = hic_deg_loops,
  txdb = txdb, orgdb = org.Hs.eg.db,
  label = "chr19 locus",
  resolution = 5000, norm = "SCALE",
  colour_lim = c(0, 10),            # explicit min/max, GENOVA colour_lim style
  min_observed_wt = MIN_OBS_WT, min_observed_deg = MIN_OBS_DEG,
  depth_scale_wt = DEPTH_SCALE_WT, depth_scale_deg = DEPTH_SCALE_DEG,
  max_fdr = 0.1,
  dist_thresh = 1.0, show_difference = FALSE
)

## THIS ONE IS GOOD USE
interval=250000
# wt_only	lost	-3.842638748	9	129325000	129350000	129337500	9	129600000	129625000	129612500	275000	13.346617698669434	0	exon+cCRE	exon+cCRE	Exon-Exon	NTMT1
plot_hic_region(
  chrom = "chr9", start = 129337500	-interval, end = 129612500+interval,
  hic_wt_path = hic_wt_file, hic_deg_path = hic_pol2_file,
  loops_wt_path = hic_wt_loops, loops_deg_path = hic_deg_loops,
  txdb = txdb, orgdb = org.Hs.eg.db,
  label = paste0("Chr9:",129337500	-interval, "-", 129612500+interval),
  resolution = 5000, norm = "SCALE",
  colour_lim = c(0, 10),   
  colour_lim_diff = c(-5, 5),          # explicit min/max, GENOVA colour_lim style
  min_observed_wt = MIN_OBS_WT, min_observed_deg = MIN_OBS_DEG,
  depth_scale_wt = DEPTH_SCALE_WT, depth_scale_deg = DEPTH_SCALE_DEG,
  max_fdr = 0.1,
  dist_thresh = 1.0, show_difference = T
)

## THIS ONE IS GOOD
interval=250000

# wt_only	lost	-3.159272915	22	32650000	32675000	32662500	22	32975000	33000000	33000000	337500	7.933793544769287	0	cCRE	cCRE	cCRE-cCRE	SYN3	Up	3.094559656
plot_hic_region(
  chrom = "chr22", start = 32662500	-interval, end = 33000000+interval,
  hic_wt_path = hic_wt_file, hic_deg_path = hic_pol2_file,
  loops_wt_path = hic_wt_loops, loops_deg_path = hic_deg_loops,
  txdb = txdb, orgdb = org.Hs.eg.db,
  label = paste0("Chr22:",32662500	-interval, "-", 32662500	+interval),
  resolution = 5000, norm = "SCALE",
  colour_lim = c(0, 10),   
  colour_lim_diff = c(-5, 5),        # explicit min/max, GENOVA colour_lim style
  min_observed_wt = MIN_OBS_WT, min_observed_deg = MIN_OBS_DEG,
  depth_scale_wt = DEPTH_SCALE_WT, depth_scale_deg = DEPTH_SCALE_DEG,
  max_fdr = 0.1,
  dist_thresh = 1.0, show_difference = T
)

## wt_only	lost	-2.716539643	7	101540000	101550000	101545000	7	101810000	101820000	101815000	270000	10.894538879394531	0.809621274	exon+cCRE	exon+cCRE	Exon-Exon

## THIS ONE IS GOOD
interval=250000

plot_hic_region(
  chrom = "chr7", start = 101545000	-interval, end = 101815000+interval,
  hic_wt_path = hic_wt_file, hic_deg_path = hic_pol2_file,
  loops_wt_path = hic_wt_loops, loops_deg_path = hic_deg_loops,
  txdb = txdb, orgdb = org.Hs.eg.db,
  label = "Chr7:101295000-101795000",
  resolution = 5000, norm = "SCALE",
  colour_lim = c(0, 10),   
  colour_lim_diff = c(-5, 5), # explicit min/max, GENOVA colour_lim style
  min_observed_wt = MIN_OBS_WT, min_observed_deg = MIN_OBS_DEG,
  depth_scale_wt = DEPTH_SCALE_WT, depth_scale_deg = DEPTH_SCALE_DEG,
  max_fdr = 0.1,
  dist_thresh = 1.0, show_difference = T
)
