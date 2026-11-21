##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Representative Hi-C contact maps with A/B compartment eigenvector
# tracks for Wildtype and POLR2A deg, for two chromosomes (reviewer request).
#
# One plot per chromosome: a square contact map with Wildtype in the lower-left
# triangle and POLR2A deg in the upper-right triangle. Each condition's A/B
# eigenvector track (A = positive, B = negative) sits outside the Mb tick marks
# on its own side: Wildtype on the left, POLR2A deg on top.
#
# Inputs (same as compartment_analysis.R, DS_GROUP "ds_wt_pol2"):
#   - Hi-C:        <BASE>/mega_ds_wt_pol2/<cond>/mega/aligned/inter_30.hic
#   - Eigenvectors <BASE>/compartments/ds_wt_pol2_<WT|Pol2>/..._compartments.bedgraph
# Contact matrices are read with strawr (SCALE-normalised), exactly as in
# hic_loop_track_plotter.R. The colour scale is shared between the two
# conditions of a chromosome so they are directly comparable.
# Plots are printed to the interactive device and also exported as one SVG per
# chromosome (EXPORT_DIR), with the contact-map panel exactly square.
##################################################################

library(strawr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(patchwork)

# ══════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════
DS_GROUP <- "ds_wt_pol2"          # "ds_wt_pol2" | "ds_all"  (as in compartment_analysis.R)

CHROMS   <- c("chr2", "chr17")    # two representative chromosomes
RES_MAP  <- 250000L               # contact map bin size (bp); eigenvectors are 500 kb
NORM     <- "SCALE"
TRACK_REL  <- 0.30               # track thickness relative to the square map (interactive plot)
COL_LOW_Q  <- 0.02                # colour limits = these quantiles of log10(contacts),
COL_HIGH_Q <- 0.995               #   pooled over both conditions of a chromosome

BASE_NEW  <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
COMP_BASE <- file.path(BASE_NEW, "compartments")

# SVG export
EXPORT_DIR <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/hi-c/scripts/plots"
MAP_IN     <- 7                   # contact-map panel side (inches); width = height
TRACK_IN   <- 1.2                 # eigenvector track thickness (inches)

if (DS_GROUP == "ds_wt_pol2") {
  CONDITIONS <- list(
    wt   = list(label = "Wildtype",
                hic      = file.path(BASE_NEW, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                comp_bed = file.path(COMP_BASE, "ds_wt_pol2_WT/ds_wt_pol2_WT_compartments.bedgraph")),
    pol2 = list(label = "POLR2A deg",
                hic      = file.path(BASE_NEW, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                comp_bed = file.path(COMP_BASE, "ds_wt_pol2_Pol2/ds_wt_pol2_Pol2_compartments.bedgraph"))
  )
} else if (DS_GROUP == "ds_all") {
  CONDITIONS <- list(
    wt   = list(label = "Wildtype",
                hic      = file.path(BASE_NEW, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                comp_bed = file.path(COMP_BASE, "ds_all_WT/ds_all_WT_compartments.bedgraph")),
    pol2 = list(label = "POLR2A deg",
                hic      = file.path(BASE_NEW, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                comp_bed = file.path(COMP_BASE, "ds_all_Pol2/ds_all_Pol2_compartments.bedgraph"))
  )
} else stop("DS_GROUP must be 'ds_wt_pol2' or 'ds_all'")

GENOVA_WHITERED <- c("white", "orange", "red", "black")   # as in hic_loop_track_plotter.R
COL_A <- "#D55E00"   # A compartment (positive eigenvector)
COL_B <- "#0072B2"   # B compartment (negative eigenvector)

# ══════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════
match_chrom_style <- function(chrom, file_chroms) {
  if (chrom %in% file_chroms) return(chrom)
  stripped <- gsub("^chr", "", chrom)
  if (stripped %in% file_chroms) return(stripped)
  prefixed <- paste0("chr", chrom)
  if (prefixed %in% file_chroms) return(prefixed)
  stop("Cannot match chromosome '", chrom, "' to file chromosomes: ",
       paste(head(file_chroms), collapse = ", "))
}

chrom_length <- function(hic_path, chrom) {
  ch <- strawr::readHicChroms(hic_path)
  ch$length[ch$name == match_chrom_style(chrom, ch$name)][1]
}

# Whole-chromosome SCALE-normalised matrix, mirrored to a full symmetric table
load_chrom_matrix <- function(hic_path, chrom, resolution, norm = NORM) {
  fc  <- match_chrom_style(chrom, strawr::readHicChroms(hic_path)$name)
  raw <- strawr::straw(norm, hic_path, fc, fc, "BP", resolution)
  if (nrow(raw) == 0) stop("straw() returned 0 rows for ", fc, " in ", hic_path)
  raw <- raw[is.finite(raw$counts), ]
  off <- raw[raw$x != raw$y, ]
  bind_rows(raw, data.frame(x = off$y, y = off$x, counts = off$counts))
}

load_eigen_chrom <- function(bed_path, chrom) {
  df <- read.table(bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                   col.names = c("chr", "start", "end", "eigen"))
  df$chr <- gsub("^chr", "", as.character(df$chr))
  df <- df[df$chr == gsub("^chr", "", chrom) & is.finite(df$eigen), ]
  if (nrow(df) == 0) stop("No eigenvector bins for ", chrom, " in ", bed_path)
  df
}

# ══════════════════════════════════════════════════════════════════
# Panels
# ══════════════════════════════════════════════════════════════════
theme_map <- function() {
  theme_classic(base_size = 11) +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# eigenvector track; side = "top" (position along x) or "left" (position along y, drawn rotated)
make_eigen_track <- function(eig, chr_len, ylim, side, label) {
  p <- ggplot(eig, aes(xmin = start, xmax = end, ymin = 0, ymax = eigen, fill = eigen > 0)) +
    geom_rect() +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    scale_fill_manual(values = c(`TRUE` = COL_A, `FALSE` = COL_B), guide = "none") +
    theme_map() +
    theme(axis.text.y = element_text(size = 8))
  if (side == "top") {
    p + scale_x_continuous(limits = c(0, chr_len), expand = c(0, 0)) +
      scale_y_continuous(limits = ylim, breaks = c(ylim[1], 0, ylim[2]),
                         labels = function(v) sprintf("%.2f", v)) +
      labs(x = NULL, y = label) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank())
  } else {
    # rotated: genomic position runs down the page (0 at top), A (positive) points left
    p + scale_x_reverse(limits = c(chr_len, 0), expand = c(0, 0)) +
      scale_y_reverse(limits = rev(ylim), breaks = c(ylim[2], 0, ylim[1]),
                      labels = function(v) sprintf("%.2f", v)) +
      coord_flip() +
      labs(x = label, y = NULL) +
      theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.line.y = element_blank())
  }
}

# square map: Wildtype in the lower-left triangle (y >= x), POLR2A deg in the upper-right (x > y);
# Mb tick marks on the top and left edges (where the tracks sit just outside them)
make_square_map <- function(mat_wt, mat_deg, chr_len, resolution, lims, legend = FALSE) {
  d_wt  <- mat_wt  %>% filter(y >= x)
  d_deg <- mat_deg %>% filter(x >  y)
  d <- bind_rows(d_wt, d_deg)
  d$v <- log10(d$counts)
  d$v[!is.finite(d$v)] <- NA
  p <- ggplot(d, aes(x = x + resolution / 2, y = y + resolution / 2, fill = v)) +
    geom_raster() +
    scale_fill_gradientn(colours = GENOVA_WHITERED, limits = lims, oob = squish,
                         na.value = "grey92", name = expression(log[10]~"contacts")) +
    scale_x_continuous(limits = c(0, chr_len), expand = c(0, 0), position = "top",
                       labels = function(v) v / 1e6) +
    scale_y_reverse(limits = c(chr_len, 0), expand = c(0, 0),
                    labels = function(v) v / 1e6) +
    labs(x = NULL, y = NULL) +
    theme_map() +
    theme(axis.text = element_text(size = 8),
          axis.line.x.top = element_line(), axis.line.x.bottom = element_blank(),
          axis.line.y.right = element_blank())
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

# One chromosome: returns the assembled plot + eigenvector concordance
build_chrom_figure <- function(chrom) {
  message("\n== ", chrom, " ==")
  chr_len <- chrom_length(CONDITIONS$wt$hic, chrom)
  mats    <- lapply(CONDITIONS, function(cd) load_chrom_matrix(cd$hic, chrom, RES_MAP))
  eigs    <- lapply(CONDITIONS, function(cd) load_eigen_chrom(cd$comp_bed, chrom))
  
  # shared colour limits (log10 contacts, both conditions) and shared eigenvector range
  allv <- log10(unlist(lapply(mats, function(m) m$counts)))
  allv <- allv[is.finite(allv)]
  lims <- as.numeric(quantile(allv, c(COL_LOW_Q, COL_HIGH_Q)))
  emax <- max(abs(unlist(lapply(eigs, function(e) e$eigen))))
  ylim <- c(-emax, emax)
  message("Colour limits (log10 contacts): ", paste(round(lims, 2), collapse = " to "))
  
  # eigenvector concordance between conditions (same 500 kb bins)
  j <- inner_join(eigs$wt, eigs$pol2, by = c("chr", "start", "end"), suffix = c("_wt", "_pol2"))
  r <- cor(j$eigen_wt, j$eigen_pol2)
  flip <- mean(sign(j$eigen_wt) != sign(j$eigen_pol2)) * 100
  message(sprintf("Eigenvector WT vs POLR2A deg: Pearson r = %.3f; %.1f%% of bins change A/B sign",
                  r, flip))
  
  trk_top  <- make_eigen_track(eigs$pol2, chr_len, ylim, "top",  CONDITIONS$pol2$label)
  trk_left <- make_eigen_track(eigs$wt,   chr_len, ylim, "left", CONDITIONS$wt$label)
  map      <- make_square_map(mats$wt, mats$pol2, chr_len, RES_MAP, lims)
  key      <- get_legend(make_square_map(mats$wt, mats$pol2, chr_len, RES_MAP, lims, legend = TRUE) +
                           theme(legend.position = "right"))
  
  body <- plot_grid(NULL, trk_top, trk_left, map, ncol = 2, align = "hv", axis = "tblr",
                    rel_widths = c(TRACK_REL, 1), rel_heights = c(TRACK_REL, 1))
  title_text <- paste0(chrom, " (Mb)   |   lower-left: ", CONDITIONS$wt$label, "   upper-right: ",
                       CONDITIONS$pol2$label, "   |   eigenvector: orange = A, blue = B")
  title <- ggdraw() + draw_label(title_text, fontface = "bold", size = 11, x = 0.02, hjust = 0)
  list(plot = plot_grid(title, plot_grid(body, key, nrow = 1, rel_widths = c(1, 0.12)),
                        ncol = 1, rel_heights = c(0.04, 1)),
       r = r, flip_pct = flip,
       parts = list(trk_top = trk_top, trk_left = trk_left, title_text = title_text,
                    map = make_square_map(mats$wt, mats$pol2, chr_len, RES_MAP, lims, legend = TRUE)))
}

# ══════════════════════════════════════════════════════════════════
# SVG export: one file per chromosome, contact-map panel exactly square
# ══════════════════════════════════════════════════════════════════
export_chrom_svg <- function(fig, chrom, outdir = EXPORT_DIR) {
  p <- wrap_plots(plot_spacer(), fig$parts$trk_top, fig$parts$trk_left, fig$parts$map, ncol = 2) +
    plot_layout(widths  = grid::unit(c(TRACK_IN, MAP_IN), "in"),
                heights = grid::unit(c(TRACK_IN, MAP_IN), "in"),
                guides  = "collect") +
    plot_annotation(title = fig$parts$title_text,
                    theme = theme(plot.title = element_text(face = "bold", size = 11)))
  # device size = fixed panels + axis/legend/title space, measured from the layout itself
  g <- patchwork::patchworkGrob(p)
  w <- grid::convertWidth(sum(g$widths),  "in", valueOnly = TRUE)
  h <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  f <- file.path(outdir, paste0(chrom, "_hic_eigenvector.svg"))
  svg(f, width = w, height = h)
  print(p)
  dev.off()
  message(sprintf("Wrote %s  (%.2f x %.2f in; map panel %.1f x %.1f in)", f, w, h, MAP_IN, MAP_IN))
  invisible(f)
}

# ══════════════════════════════════════════════════════════════════
# Run
# ══════════════════════════════════════════════════════════════════
run_contact_map_figure <- function() {
  inputs <- unlist(lapply(CONDITIONS, function(cd) c(cd$hic, cd$comp_bed)))
  missing <- inputs[!file.exists(inputs)]
  if (length(missing) > 0)
    stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "))
  
  if (!dir.exists(EXPORT_DIR)) stop("Export directory not found: ", EXPORT_DIR)
  
  graphics.off()
  figs <- list()
  for (ch in CHROMS) {
    figs[[ch]] <- build_chrom_figure(ch)
    print(figs[[ch]]$plot)          # interactive: one plot per chromosome
    export_chrom_svg(figs[[ch]], ch)  # SVG: map panel exactly square
  }
  invisible(list(per_chrom = lapply(figs, `[[`, "plot"),
                 eigen_cor = sapply(figs, `[[`, "r"),
                 eigen_flip_pct = sapply(figs, `[[`, "flip_pct")))
}

results_map <- run_contact_map_figure()
# re-draw one chromosome:  print(results_map$per_chrom[["chr2"]])