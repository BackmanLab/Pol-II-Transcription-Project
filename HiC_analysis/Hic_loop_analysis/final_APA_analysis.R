##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Aggregate Peak Analysis (APA) pipeline. HiCCUPS loops +
# narrowPeak anchor annotation, mariner-based contact extraction
# (NORM="NONE"), by pair type (P-P/E-P/E-E) and loop-length tertile.
#
# Figure output: Fig. 5B, 5E (probable); SI Fig. 8B-D (confirmed).
##################################################################
library(mariner)
library(strawr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(cowplot)
library(GenomicRanges)
library(InteractionSet)
library(RColorBrewer)
library(viridis)
source("~/Transcription Paper/Hic_loop_analysis/final_scripts/EP_analysis_functions.R")

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════

# ── Downsampling group switch ──────────────────────────────────────────────────
# "ds_wt_pol2" — WT + Pol2 matched to Pol2 depth (~971M reads); no ActD
# "ds_all"     — WT + Pol2 + ActD all matched to ActD depth (~712M reads)
# "unsampled"  — undownsampled originals; depth scaling applied
DS_GROUP <- "ds_wt_pol2"   # ← change this to switch groups

BASE     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
LOOP_DIR <- file.path(BASE, "loops")

if (DS_GROUP == "ds_wt_pol2") {
  # Downsampled: WT and Pol2 already matched — no depth scaling needed
  CONDITIONS <- list(
    wt   = list(
      label = "Wildtype",
      hic   = file.path(BASE, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
      loops = file.path(LOOP_DIR, "ds_wt_pol2/WT_HCT116_CTRL/mega/merged_loops.bedpe"),
      color = "#88CCEE"
    ),
    pol2 = list(
      label = "POLR2A deg",
      hic   = file.path(BASE, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
      loops = file.path(LOOP_DIR, "ds_wt_pol2/Pol2_6hrs_Aux/mega/merged_loops.bedpe"),
      color = "#CC6677"
    )
  )
  DEPTH_SCALE  <- list(wt = 1.0, pol2 = 1.0)
  MIN_OBSERVED <- list(wt = 0,  pol2 = 0)
  
} else if (DS_GROUP == "ds_all") {
  # Downsampled: WT, Pol2, ActD all matched — no depth scaling needed
  CONDITIONS <- list(
    wt   = list(
      label = "Wildtype",
      hic   = file.path(BASE, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
      loops = file.path(LOOP_DIR, "ds_all/WT_HCT116_CTRL/mega/merged_loops.bedpe"),
      color = "#88CCEE"
    ),
    pol2 = list(
      label = "POLR2A deg",
      hic   = file.path(BASE, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
      loops = file.path(LOOP_DIR, "ds_all/Pol2_6hrs_Aux/mega/merged_loops.bedpe"),
      color = "#CC6677"
    ),
    actd = list(
      label = "ActD 1hr",
      hic   = file.path(BASE, "mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic"),
      loops = file.path(LOOP_DIR, "ds_all/1hr_ActD/mega/merged_loops.bedpe"),
      color = "#44AA99"
    )
  )
  DEPTH_SCALE  <- list(wt = 1.0, pol2 = 1.0, actd = 1.0)
  MIN_OBSERVED <- list(wt = 10,  pol2 = 25,   actd = 25)
  
} else if (DS_GROUP == "unsampled") {
  # Unsampled: different depths — scale to WT
  wt_total   <- 639409109 + 377220354
  pol2_total <- 374163054 + 249677343
  actd_total <- 266069792 + 202465027
  CONDITIONS <- list(
    wt   = list(
      label = "Wildtype",
      hic   = file.path(BASE, "mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
      loops = file.path(LOOP_DIR, "unsampled/WT_HCT116_CTRL/mega/merged_loops.bedpe"),
      color = "#88CCEE"
    ),
    pol2 = list(
      label = "POLR2A deg",
      hic   = file.path(BASE, "mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
      loops = file.path(LOOP_DIR, "unsampled/Pol2_6hrs_Aux/mega/merged_loops.bedpe"),
      color = "#CC6677"
    ),
    actd = list(
      label = "ActD 1hr",
      hic   = file.path(BASE, "mega_unsampled/1hr_ActD/mega/aligned/inter_30.hic"),
      loops = file.path(LOOP_DIR, "unsampled/1hr_ActD/mega/merged_loops.bedpe"),
      color = "#44AA99"
    )
  )
  DEPTH_SCALE  <- list(wt = 1.0,
                       pol2 = wt_total / pol2_total,
                       actd = wt_total / actd_total)
  MIN_OBSERVED <- list(wt = 10, pol2 = 25, actd = 25)
  
} else {
  stop("Unknown DS_GROUP: '", DS_GROUP,
       "'. Must be 'ds_wt_pol2', 'ds_all', or 'unsampled'.")
}

# ── APA parameters ─────────────────────────────────────────────────────────────
RESOLUTION     <- 5000
BUFFER         <- 7
NORM           <- "NONE"
GENE_TOL_BP    <- 5000
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1

# ── File paths — expression, ChIP, ATAC ───────────────────────────────────────
EXPR_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
DESEQ_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
CHIP_DIR   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"
ATAC_DIR   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"

# ── cCRE group constants ───────────────────────────────────────────────────────
GROUP_LEVELS     <- c("P-P","E-P","E-E")
GROUP_COLORS     <- c("P-P"="#882255","E-P"="#44AA99","E-E"="#4393c3")
SINGLETON_LEVELS <- c("P-Gene","pE-Gene","dE-Gene")
SINGLETON_COLORS <- c("P-Gene"="#b2182b","pE-Gene"="#f4a582","dE-Gene"="#d6604d")
BINS       <- c("Low","Medium","High")
BIN_COLORS <- c("Low"="#fdae61","Medium"="#f46d43","High"="#d73027")

# ── Palette selectors ──────────────────────────────────────────────────────────
ACTIVE_PALETTE      <- "viridis"
ACTIVE_DIFF_PALETTE <- "PRGn"

# ══════════════════════════════════════════════════════════════════════════════
# Functions — palette helpers
# ══════════════════════════════════════════════════════════════════════════════
make_fill_scale <- function(palette_name, limits,
                            name = "Fold over\nbackground") {
  switch(palette_name,
         red        = scale_fill_gradient(low="white", high="#cc0000", limits=limits, name=name, oob=scales::squish),
         blue       = scale_fill_gradient(low="white", high="#1565c0", limits=limits, name=name, oob=scales::squish),
         purple     = scale_fill_gradient(low="white", high="#6a0dad", limits=limits, name=name, oob=scales::squish),
         orange_red = scale_fill_gradientn(colors=c("white","#fdae61","#f46d43","#d73027"), limits=limits, name=name, oob=scales::squish),
         viridis    = scale_fill_viridis_c(option="viridis", limits=limits, name=name, oob=scales::squish),
         magma      = scale_fill_viridis_c(option="magma",   limits=limits, name=name, oob=scales::squish),
         ylOrRd     = scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9,"YlOrRd"), limits=limits, name=name, oob=scales::squish),
         blues      = scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9,"Blues"),  limits=limits, name=name, oob=scales::squish),
         scale_fill_gradient(low="white", high="#cc0000", limits=limits, name=name, oob=scales::squish)
  )
}

make_diff_scale <- function(palette_name, limits, name = "log2\n(treat/WT)") {
  colors <- switch(palette_name,
                   PRGn       = RColorBrewer::brewer.pal(11, "PRGn"),
                   PuOr       = rev(RColorBrewer::brewer.pal(11, "PuOr")),
                   teal_coral = c("#00796b","#26a69a","#80cbc4","#e0f2f1","white",
                                  "#fce4ec","#ef9a9a","#e57373","#c62828"),
                   RdBu       = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                   rev(RColorBrewer::brewer.pal(11, "RdBu"))
  )
  scale_fill_gradientn(colors=colors, limits=c(-limits, limits),
                       name=name, oob=scales::squish)
}

# ══════════════════════════════════════════════════════════════════════════════
# Functions — loop loading
# ══════════════════════════════════════════════════════════════════════════════
load_hiccups_bedpe <- function(path) {
  all_lines   <- readLines(path)
  header_line <- all_lines[grepl("^#chr1", all_lines)][1]
  col_names   <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  data_lines  <- all_lines[!grepl("^#", all_lines)]
  df <- read.table(text=paste(data_lines, collapse="\n"),
                   header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
  if (ncol(df) == length(col_names)) colnames(df) <- col_names
  data.frame(
    chr1      = as.character(df[[1]]),
    start1    = as.integer(df[[2]]),
    end1      = as.integer(df[[3]]),
    chr2      = as.character(df[[4]]),
    start2    = as.integer(df[[5]]),
    end2      = as.integer(df[[6]]),
    observed  = as.numeric(df[["observed"]]),
    fdrDonut  = as.numeric(df[["fdrDonut"]]),
    fdrBL     = as.numeric(df[["fdrBL"]]),
    fdrH      = as.numeric(df[["fdrH"]]),
    fdrV      = as.numeric(df[["fdrV"]]),
    centroid1 = as.integer(df[["centroid1"]]),
    centroid2 = as.integer(df[["centroid2"]])
  ) %>% filter(!is.na(start1), !is.na(start2))
}

filter_loops <- function(df, min_observed = 10, max_fdr = 0.1) {
  df %>% filter(
    observed >= min_observed,
    fdrDonut <= max_fdr,
    fdrBL    <= max_fdr,
    fdrH     <= max_fdr,
    fdrV     <= max_fdr
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# Functions — gene linking and singletons
# ══════════════════════════════════════════════════════════════════════════════
link_loops_to_genes <- function(loops_df, gene_gr, tol_bp = GENE_TOL_BP) {
  loops_df <- loops_df %>% dplyr::mutate(loop_idx = row_number())
  fix_chr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if (has_chr_ref && !has_chr_gr)  seqlevels(gr) <- paste0("chr", seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    gr
  }
  anchor1_gr <- fix_chr(GRanges(loops_df$chr1, IRanges(loops_df$start1, loops_df$end1)), gene_gr)
  anchor2_gr <- fix_chr(GRanges(loops_df$chr2, IRanges(loops_df$start2, loops_df$end2)), gene_gr)
  gene_gr_exp <- trim(resize(gene_gr, width(gene_gr) + 2*tol_bp, fix="center"))
  hits1 <- findOverlaps(anchor1_gr, gene_gr_exp)
  hits2 <- findOverlaps(anchor2_gr, gene_gr_exp)
  all_pairs <- bind_rows(
    data.frame(loop_idx=queryHits(hits1), gene_idx=subjectHits(hits1)),
    data.frame(loop_idx=queryHits(hits2), gene_idx=subjectHits(hits2))
  ) %>% distinct(loop_idx, gene_idx) %>%
    dplyr::mutate(
      expr_bin  = gene_gr$expr_bin[gene_idx],
      de_status = gene_gr$de_status[gene_idx],
      TPM       = gene_gr$TPM[gene_idx],
      symbol    = gene_gr$symbol[gene_idx]
    )
  loop_annotations <- all_pairs %>%
    dplyr::group_by(loop_idx) %>%
    dplyr::slice_max(TPM, n=1, with_ties=FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(loop_idx, expr_bin, de_status, TPM, symbol)
  loops_df %>%
    left_join(loop_annotations, by="loop_idx") %>%
    dplyr::mutate(
      expr_bin  = replace_na(expr_bin,  "Unannotated"),
      de_status = replace_na(de_status, "Not DE"),
      symbol    = replace_na(symbol,    "Unknown")
    )
}

find_singleton_loops <- function(classified_loops, gene_gr, tol_bp = GENE_TOL_BP) {
  anchor1_gr <- GRanges(classified_loops$chr1, IRanges(classified_loops$start1, classified_loops$end1))
  anchor2_gr <- GRanges(classified_loops$chr2, IRanges(classified_loops$start2, classified_loops$end2))
  gene_gr_m    <- gene_gr
  has_chr_loop <- any(grepl("^chr", as.character(seqnames(anchor1_gr))))
  has_chr_gene <- any(grepl("^chr", as.character(seqnames(gene_gr))))
  if (has_chr_loop && !has_chr_gene)
    seqlevels(gene_gr_m) <- paste0("chr", seqlevels(gene_gr_m))
  else if (!has_chr_loop && has_chr_gene)
    seqlevels(gene_gr_m) <- gsub("^chr", "", seqlevels(gene_gr_m))
  gene_gr_exp <- trim(resize(gene_gr_m, width(gene_gr_m) + 2*tol_bp, fix="center"))
  hits1 <- findOverlaps(anchor1_gr, gene_gr_exp)
  hits2 <- findOverlaps(anchor2_gr, gene_gr_exp)
  df1 <- data.frame(loop_idx=queryHits(hits1), gene_idx=subjectHits(hits1),
                    other_ccre=classified_loops$anchor2_ccre[queryHits(hits1)],
                    stringsAsFactors=FALSE)
  df2 <- data.frame(loop_idx=queryHits(hits2), gene_idx=subjectHits(hits2),
                    other_ccre=classified_loops$anchor1_ccre[queryHits(hits2)],
                    stringsAsFactors=FALSE)
  bind_rows(df1, df2) %>%
    dplyr::mutate(
      de_status = gene_gr_m$de_status[gene_idx],
      symbol    = gene_gr_m$symbol[gene_idx]
    ) %>%
    filter(!is.na(other_ccre), other_ccre != "Unannotated") %>%
    dplyr::group_by(loop_idx, gene_idx) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      singleton_group = dplyr::case_when(
        other_ccre == "Promoter"           ~ "P-Gene",
        other_ccre == "Proximal-Enhancer"  ~ "pE-Gene",
        other_ccre == "Distal-Enhancer"    ~ "dE-Gene",
        TRUE                               ~ NA_character_
      )
    ) %>%
    filter(!is.na(singleton_group)) %>%
    dplyr::mutate(singleton_group = factor(singleton_group, levels=SINGLETON_LEVELS))
}

# ══════════════════════════════════════════════════════════════════════════════
# Functions — APA extraction
# ══════════════════════════════════════════════════════════════════════════════
loops_to_gi <- function(loops_df, resolution = RESOLUTION) {
  snap <- function(x) floor(x / resolution) * resolution
  pos1 <- snap(if ("centroid1" %in% colnames(loops_df) && !all(is.na(loops_df$centroid1)))
    loops_df$centroid1 else (loops_df$start1 + loops_df$end1) / 2)
  pos2 <- snap(if ("centroid2" %in% colnames(loops_df) && !all(is.na(loops_df$centroid2)))
    loops_df$centroid2 else (loops_df$start2 + loops_df$end2) / 2)
  gi <- GInteractions(
    anchor1         = GRanges(loops_df$chr1, IRanges(pos1, pos1 + resolution - 1)),
    anchor2         = GRanges(loops_df$chr2, IRanges(pos2, pos2 + resolution - 1)),
    expr_bin        = if ("expr_bin"        %in% colnames(loops_df)) loops_df$expr_bin        else "All",
    de_status       = if ("de_status"       %in% colnames(loops_df)) loops_df$de_status       else "All",
    group           = if ("group"           %in% colnames(loops_df)) loops_df$group           else "All",
    singleton_group = if ("singleton_group" %in% colnames(loops_df)) loops_df$singleton_group else "All"
  )
  GenomeInfoDb::seqlevelsStyle(gi) <- "ENSEMBL"
  gi
}

extract_apa_matrices_by <- function(loops_df, hic_file, group_col,
                                    group_levels, resolution = RESOLUTION,
                                    buffer = BUFFER, norm = NORM,
                                    min_loops = 10) {
  gi_all   <- loops_to_gi(loops_df, resolution)
  min_dist <- (2*buffer + 1) * resolution
  loop_dist <- abs(start(anchors(gi_all, "first")) - start(anchors(gi_all, "second")))
  gi_all   <- gi_all[loop_dist >= min_dist]
  message("Loops after diagonal filter: ", length(gi_all))
  results <- list()
  for (grp in group_levels) {
    gi_grp <- gi_all[mcols(gi_all)[[group_col]] == grp]
    if (length(gi_grp) < min_loops) {
      message("Skipping ", grp, " — ", length(gi_grp), " loops"); next
    }
    message("Extracting ", grp, " (", length(gi_grp), " loops)...")
    tryCatch({
      regions <- mariner::pixelsToMatrices(x=gi_grp, buffer=buffer)
      mats    <- mariner::pullHicMatrices(x=regions, files=hic_file,
                                          binSize=resolution, norm=norm,
                                          half="both", matrix="observed")
      results[[grp]] <- mats
    }, error = function(e) message("  Error: ", e$message))
  }
  results
}

extract_apa_matrices <- function(loops_annotated, hic_file,
                                 resolution = RESOLUTION, buffer = BUFFER,
                                 norm = NORM, min_loops = 10) {
  extract_apa_matrices_by(loops_annotated, hic_file,
                          group_col="expr_bin", group_levels=BINS,
                          resolution=resolution, buffer=buffer,
                          norm=norm, min_loops=min_loops)
}

aggregate_apa <- function(mat_list, scale_factor = 1.0) {
  lapply(mat_list, function(mats) {
    if (is.null(mats)) return(NULL)
    cnt      <- counts(mats)
    raw_mean <- apply(cnt[,,,1, drop=FALSE], c(1,2), mean, na.rm=TRUE)
    raw_mean * scale_factor
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# Functions — APA scoring and plotting
# ══════════════════════════════════════════════════════════════════════════════
apa_score <- function(mat, corner_size = 3) {
  if (is.null(mat)) return(NA_real_)
  n <- nrow(mat); mid <- ceiling(n / 2)
  central <- mat[mid, mid]
  corners <- c(as.vector(mat[1:corner_size, 1:corner_size]),
               as.vector(mat[1:corner_size, (n-corner_size+1):n]),
               as.vector(mat[(n-corner_size+1):n, 1:corner_size]),
               as.vector(mat[(n-corner_size+1):n, (n-corner_size+1):n]))
  corner_mean <- mean(corners[!is.na(corners)], na.rm=TRUE)
  if (is.na(corner_mean) || corner_mean == 0) return(NA_real_)
  central / corner_mean
}

get_group_n <- function(loops_df, group_col, group_levels, de_filter = NULL) {
  df <- as.data.frame(loops_df)
  if (!is.null(de_filter)) df <- df %>% filter(de_status == de_filter)
  bind_rows(lapply(group_levels, function(g) {
    grp_df <- df %>% filter(.data[[group_col]] == g)
    data.frame(bin=g, n_loops=nrow(grp_df),
               n_genes=dplyr::n_distinct(
                 grp_df$symbol[!is.na(grp_df$symbol) & grp_df$symbol != "Unknown"]))
  })) %>%
    dplyr::mutate(
      bin = factor(bin, levels=group_levels),
      axis_label = ifelse(n_genes > 0,
                          paste0(bin, "\n", n_loops, " loops / ", n_genes, " genes"),
                          paste0(bin, "\n", n_loops, " loops"))
    )
}

plot_apa_grid_generic <- function(apa_agg, bin_levels, bin_colors,
                                  title = "APA", palette = ACTIVE_PALETTE,
                                  corner_size = 3,
                                  conditions = names(CONDITIONS)) {
  normalize_to_bg <- function(mat) {
    if (is.null(mat)) return(NULL)
    n <- nrow(mat)
    corners <- c(as.vector(mat[1:corner_size, 1:corner_size]),
                 as.vector(mat[1:corner_size, (n-corner_size+1):n]),
                 as.vector(mat[(n-corner_size+1):n, 1:corner_size]),
                 as.vector(mat[(n-corner_size+1):n, (n-corner_size+1):n]))
    cm <- mean(corners[!is.na(corners)], na.rm=TRUE)
    if (is.na(cm) || cm == 0) return(mat)
    mat / cm
  }
  apa_norm   <- lapply(apa_agg, function(cl) lapply(cl, normalize_to_bg))
  all_vals   <- unlist(lapply(apa_norm, function(cl)
    unlist(lapply(cl, function(m) if (!is.null(m)) as.vector(m)))))
  global_max <- quantile(all_vals, 0.98, na.rm=TRUE)
  global_min <- 1.0
  
  make_panel <- function(mat_norm) {
    if (is.null(mat_norm)) return(ggplot() + theme_void() +
                                    theme(plot.background=element_rect(fill="gray95", color=NA), plot.margin=margin(2,2,2,2)))
    n <- nrow(mat_norm); mid <- ceiling(n/2); score <- mat_norm[mid, mid]
    df <- expand.grid(row=1:n, col=1:n) %>%
      dplyr::mutate(value=as.vector(mat_norm), x=col-mid, y=-(row-mid))
    p <- ggplot(df, aes(x=x, y=y, fill=pmin(pmax(value, global_min), global_max))) +
      geom_tile() +
      make_fill_scale(palette, limits=c(global_min, global_max), name="Fold over\nbackground") +
      annotate("rect", xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5,
               fill=NA, color="black", linewidth=0.6) +
      coord_fixed() +
      scale_x_continuous(breaks=c(-BUFFER, 0, BUFFER),
                         labels=c(paste0("-", BUFFER*RESOLUTION/1000, "kb"),
                                  "loop", paste0("+", BUFFER*RESOLUTION/1000, "kb"))) +
      scale_y_continuous(breaks=c(-BUFFER, 0, BUFFER),
                         labels=c(paste0("-", BUFFER*RESOLUTION/1000, "kb"),
                                  "loop", paste0("+", BUFFER*RESOLUTION/1000, "kb"))) +
      labs(x=NULL, y=NULL) + theme_classic(base_size=8) +
      theme(axis.text=element_text(size=9), axis.line=element_blank(),
            axis.ticks=element_blank(), legend.position="none",
            plot.margin=margin(2,2,2,2))
    if (!is.na(score))
      p <- p + annotate("text", x=-BUFFER+0.5, y=BUFFER-0.5,
                        label=paste0("APA=", round(score, 2)),
                        size=2.5, hjust=0, fontface="bold", color="white")
    p
  }
  col_headers <- lapply(bin_levels, function(b) {
    ggplot() + annotate("text", x=0.5, y=0.5, label=b,
                        color=bin_colors[[b]], size=3.5, fontface="bold") +
      theme_void() + theme(plot.margin=margin(0,2,0,2))
  })
  plot_list <- c(list(ggplot() + theme_void()), col_headers)
  for (n in conditions) {
    row_panels <- lapply(bin_levels, function(b) make_panel(apa_norm[[n]][[b]]))
    label_p <- ggplot() +
      annotate("text", x=0.5, y=0.5, label=CONDITIONS[[n]]$label,
               color=CONDITIONS[[n]]$color, size=3.5, fontface="bold", angle=90) +
      theme_void()
    plot_list <- c(plot_list, list(label_p), row_panels)
  }
  main_grid <- wrap_plots(plot_list, ncol=length(bin_levels)+1,
                          widths=c(0.12, rep(1, length(bin_levels))),
                          heights=c(0.15, rep(1, length(conditions))))
  legend_df <- data.frame(x=1:10, y=1:10,
                          fill=seq(global_min, global_max, length.out=10))
  legend_plot <- ggplot(legend_df, aes(x=x, y=y, fill=fill)) + geom_tile() +
    make_fill_scale(palette, limits=c(global_min, global_max), name="Fold over\nbackground") +
    theme_void() +
    theme(legend.position="right", legend.key.height=unit(1.5,"cm"),
          legend.key.width=unit(0.3,"cm"), legend.text=element_text(size=8),
          legend.title=element_text(size=9))
  composite <- cowplot::plot_grid(main_grid, cowplot::get_legend(legend_plot),
                                  ncol=2, rel_widths=c(1, 0.07))
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0, y=0.75) +
    cowplot::draw_label(paste0("Background-normalised  |  ", RESOLUTION/1000,
                               "kb  |  ±", BUFFER, " bins  |  ",
                               "APA=fold-enrichment at loop pixel  |  palette: ", palette),
                        size=8, color="gray40", x=0, hjust=0, y=0.25)
  final <- cowplot::plot_grid(title_grob, composite, ncol=1, rel_heights=c(0.08, 1))
  print(final); invisible(final)
}

plot_apa_grid <- function(apa_agg, bins = BINS,
                          title = "APA by expression level",
                          palette = ACTIVE_PALETTE) {
  all_vals   <- unlist(lapply(apa_agg, function(cl)
    unlist(lapply(cl, function(m) if (!is.null(m)) as.vector(m)))))
  global_max <- quantile(all_vals, 0.98, na.rm=TRUE)
  make_panel <- function(mat) {
    if (is.null(mat)) return(ggplot() + theme_void() +
                               theme(plot.background=element_rect(fill="gray95", color=NA), plot.margin=margin(2,2,2,2)))
    n <- nrow(mat); mid <- ceiling(n/2); score <- apa_score(mat)
    df <- expand.grid(row=1:n, col=1:n) %>%
      dplyr::mutate(value=as.vector(mat), x=col-mid, y=-(row-mid))
    p <- ggplot(df, aes(x=x, y=y, fill=pmin(pmax(value, 0), global_max))) +
      geom_tile() +
      make_fill_scale(palette, limits=c(0, global_max), name="Contacts") +
      annotate("rect", xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5,
               fill=NA, color="black", linewidth=0.6) +
      coord_fixed() +
      scale_x_continuous(breaks=c(-BUFFER, 0, BUFFER),
                         labels=c(paste0("-", BUFFER*RESOLUTION/1000, "kb"),
                                  "loop", paste0("+", BUFFER*RESOLUTION/1000, "kb"))) +
      scale_y_continuous(breaks=c(-BUFFER, 0, BUFFER),
                         labels=c(paste0("-", BUFFER*RESOLUTION/1000, "kb"),
                                  "loop", paste0("+", BUFFER*RESOLUTION/1000, "kb"))) +
      labs(x=NULL, y=NULL) + theme_classic(base_size=8) +
      theme(axis.text=element_text(size=6), axis.line=element_blank(),
            axis.ticks=element_blank(), legend.position="none",
            plot.margin=margin(2,2,2,2))
    if (!is.na(score))
      p <- p + annotate("text", x=-BUFFER+0.5, y=BUFFER-0.5,
                        label=paste0("APA=", round(score, 2)),
                        size=2.5, hjust=0, fontface="bold", color="white")
    p
  }
  col_headers <- lapply(bins, function(b) {
    ggplot() + annotate("text", x=0.5, y=0.5, label=b,
                        color=BIN_COLORS[b], size=3.5, fontface="bold") +
      theme_void() + theme(plot.margin=margin(0,2,0,2))
  })
  plot_list <- c(list(ggplot() + theme_void()), col_headers)
  for (n in names(CONDITIONS)) {
    row_panels <- lapply(bins, function(b) make_panel(apa_agg[[n]][[b]]))
    label_p <- ggplot() +
      annotate("text", x=0.5, y=0.5, label=CONDITIONS[[n]]$label,
               color=CONDITIONS[[n]]$color, size=3.5, fontface="bold", angle=90) +
      theme_void()
    plot_list <- c(plot_list, list(label_p), row_panels)
  }
  main_grid <- wrap_plots(plot_list, ncol=length(bins)+1,
                          widths=c(0.12, rep(1, length(bins))),
                          heights=c(0.15, rep(1, length(CONDITIONS))))
  legend_df <- data.frame(x=1:10, y=1:10, fill=seq(0, global_max, length.out=10))
  legend_plot <- ggplot(legend_df, aes(x=x, y=y, fill=fill)) + geom_tile() +
    make_fill_scale(palette, limits=c(0, global_max), name="Contacts") +
    theme_void() +
    theme(legend.position="right", legend.key.height=unit(1.5,"cm"),
          legend.key.width=unit(0.3,"cm"), legend.text=element_text(size=8),
          legend.title=element_text(size=9))
  composite <- cowplot::plot_grid(main_grid, cowplot::get_legend(legend_plot),
                                  ncol=2, rel_widths=c(1, 0.07))
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0, y=0.75) +
    cowplot::draw_label(paste0("Raw contacts  |  ", RESOLUTION/1000,
                               "kb  |  ±", BUFFER, " bins  |  palette: ", palette),
                        size=8, color="gray40", x=0, hjust=0, y=0.25)
  final <- cowplot::plot_grid(title_grob, composite, ncol=1, rel_heights=c(0.08, 1))
  print(final); invisible(final)
}

plot_apa_scores <- function(apa_agg, loops_annotated_cond,
                            de_filter = NULL, bins = BINS,
                            group_col = "expr_bin",
                            conditions = names(CONDITIONS),
                            title = "APA enrichment scores") {
  scores <- bind_rows(lapply(conditions, function(n) {
    bind_rows(lapply(bins, function(b) {
      data.frame(condition=CONDITIONS[[n]]$label, bin=b,
                 score=apa_score(apa_agg[[n]][[b]]))
    }))
  })) %>%
    dplyr::mutate(
      condition = factor(condition, levels=sapply(CONDITIONS[conditions], `[[`, "label")),
      bin       = factor(bin, levels=bins)
    )
  
  # ── x-axis labels ─────────────────────────────────────────────────────────
  if (is.list(loops_annotated_cond) && !is.data.frame(loops_annotated_cond)) {
    bin_n_list <- bind_rows(lapply(names(loops_annotated_cond), function(n) {
      df <- get_group_n(loops_annotated_cond[[n]], group_col=group_col,
                        group_levels=bins, de_filter=de_filter)
      df$cond_label <- CONDITIONS[[n]]$label
      df
    }))
    x_labels <- setNames(
      sapply(bins, function(b) {
        rows <- bin_n_list %>% filter(bin == b)
        paste0(b, "\n",
               paste0(rows$cond_label, ": ", rows$n_loops, " loops",
                      collapse="\n"))
      }),
      bins
    )
    bin_n <- bin_n_list
  } else {
    bin_n    <- get_group_n(loops_annotated_cond, group_col=group_col,
                            group_levels=bins, de_filter=de_filter)
    x_labels <- setNames(bin_n$axis_label, as.character(bin_n$bin))
  }
  
  treat_conds  <- conditions[conditions != "wt"]
  cond_colors  <- setNames(sapply(CONDITIONS[conditions], `[[`, "color"),
                           sapply(CONDITIONS[conditions], `[[`, "label"))
  treat_colors <- setNames(sapply(CONDITIONS[treat_conds], `[[`, "color"),
                           sapply(CONDITIONS[treat_conds], `[[`, "label"))
  
  p_raw <- scores %>% filter(!is.na(score)) %>%
    ggplot(aes(x=bin, y=score, color=condition, group=condition)) +
    geom_line(linewidth=1.1) + geom_point(size=3.5) +
    geom_text(aes(label=round(score, 2)), vjust=-0.9, size=2.8, show.legend=FALSE) +
    scale_color_manual(values=cond_colors, name=NULL) +
    scale_x_discrete(drop=FALSE, labels=x_labels) +
    scale_y_continuous(name="APA score (P2LL)", expand=expansion(mult=c(0.1, 0.2))) +
    labs(title=paste0(title, " — raw scores"), subtitle="P2LL", x=NULL) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold", size=13),
          plot.subtitle=element_text(size=9, color="gray40"),
          axis.text.x=element_text(size=9), legend.position="bottom")
  
  wt_scores <- scores %>% filter(condition == CONDITIONS$wt$label) %>%
    dplyr::select(bin, score_wt=score)
  fc_scores <- scores %>% filter(condition != CONDITIONS$wt$label) %>%
    left_join(wt_scores, by="bin") %>%
    dplyr::mutate(log2fc     = log2(score / score_wt),
                  pct_change = ((score - score_wt) / score_wt) * 100)
  
  p_fc <- fc_scores %>% filter(!is.na(log2fc)) %>%
    ggplot(aes(x=bin, y=log2fc, color=condition, group=condition)) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50", linewidth=0.6) +
    geom_line(linewidth=1.1) + geom_point(size=3.5) +
    geom_text(aes(label=round(log2fc, 2)), vjust=-0.9, size=2.8, show.legend=FALSE) +
    scale_color_manual(values=treat_colors, name=NULL) +
    scale_x_discrete(drop=FALSE, labels=x_labels) +
    scale_y_continuous(name="log2 fold change vs WT", expand=expansion(mult=c(0.2, 0.2))) +
    labs(title=paste0(title, " — log2FC vs WT"),
         subtitle="log2(treatment/WT)  |  dashed=no change", x=NULL) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold", size=13),
          plot.subtitle=element_text(size=9, color="gray40"),
          axis.text.x=element_text(size=9), legend.position="bottom")
  
  composite <- (p_raw / p_fc) + plot_layout(heights=c(1,1), guides="collect") &
    theme(legend.position="bottom")
  composite <- composite +
    plot_annotation(title=title,
                    theme=theme(plot.title=element_text(face="bold", size=14)))
  print(composite)
  invisible(list(composite=composite, raw=p_raw, fc=p_fc,
                 data=scores, fc_data=fc_scores, bin_n=bin_n))
}

plot_apa_diff <- function(apa_agg, bins = BINS,
                          title = "APA difference — treatment vs WT",
                          pseudocount = 0.01, palette = ACTIVE_DIFF_PALETTE,
                          corner_size = 3,
                          conditions = names(CONDITIONS)[names(CONDITIONS) != "wt"]) {
  norm_to_bg <- function(mat) {
    if (is.null(mat)) return(NULL)
    n <- nrow(mat)
    corners <- c(as.vector(mat[1:corner_size, 1:corner_size]),
                 as.vector(mat[1:corner_size, (n-corner_size+1):n]),
                 as.vector(mat[(n-corner_size+1):n, 1:corner_size]),
                 as.vector(mat[(n-corner_size+1):n, (n-corner_size+1):n]))
    cm <- mean(corners[!is.na(corners)], na.rm=TRUE)
    if (is.na(cm) || cm == 0) return(mat)
    mat / cm
  }
  compute_log2ratio <- function(treat_mat, wt_mat) {
    if (is.null(treat_mat) || is.null(wt_mat)) return(NULL)
    log2((norm_to_bg(treat_mat) + pseudocount) / (norm_to_bg(wt_mat) + pseudocount))
  }
  diffs <- lapply(conditions, function(n) {
    result <- lapply(bins, function(b)
      compute_log2ratio(apa_agg[[n]][[b]], apa_agg$wt[[b]]))
    names(result) <- bins; result
  }); names(diffs) <- conditions
  all_vals <- unlist(lapply(diffs, function(cond)
    unlist(lapply(cond, function(m) if (!is.null(m)) as.vector(m)))))
  limit <- quantile(abs(all_vals[!is.na(all_vals)]), 0.98)
  make_diff_panel <- function(mat) {
    if (is.null(mat)) return(ggplot() + theme_void() +
                               theme(plot.background=element_rect(fill="gray95", color=NA), plot.margin=margin(2,2,2,2)))
    n <- nrow(mat); mid <- ceiling(n/2)
    df <- expand.grid(row=1:n, col=1:n) %>%
      dplyr::mutate(value=as.vector(mat), x=col-mid, y=-(row-mid))
    p <- ggplot(df, aes(x=x, y=y, fill=pmin(pmax(value, -limit), limit))) +
      geom_tile() + make_diff_scale(palette, limits=limit) +
      annotate("rect", xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5,
               fill=NA, color="black", linewidth=0.6) +
      coord_fixed() +
      scale_x_continuous(breaks=c(-BUFFER, 0, BUFFER),
                         labels=c(paste0("-", BUFFER*RESOLUTION/1000, "kb"),
                                  "loop", paste0("+", BUFFER*RESOLUTION/1000, "kb"))) +
      scale_y_continuous(breaks=c(-BUFFER, 0, BUFFER),
                         labels=c(paste0("-", BUFFER*RESOLUTION/1000, "kb"),
                                  "loop", paste0("+", BUFFER*RESOLUTION/1000, "kb"))) +
      labs(x=NULL, y=NULL) + theme_classic(base_size=8) +
      theme(axis.text=element_text(size=9), axis.line=element_blank(),
            axis.ticks=element_blank(), legend.position="none",
            plot.margin=margin(2,2,2,2))
    score <- mat[mid, mid]
    if (!is.na(score))
      p <- p + annotate("text", x=-BUFFER+0.5, y=BUFFER-0.5,
                        label=paste0("Δlog2=", round(score, 2)),
                        size=2.5, hjust=0, fontface="bold")
    p
  }
  col_headers <- lapply(bins, function(b) {
    col <- if (b %in% names(BIN_COLORS)) BIN_COLORS[b] else
      if (b %in% names(GROUP_COLORS)) GROUP_COLORS[b] else
        if (b %in% names(SINGLETON_COLORS)) SINGLETON_COLORS[b] else "gray40"
    ggplot() + annotate("text", x=0.5, y=0.5, label=b,
                        color=col, size=3.5, fontface="bold") + theme_void()
  })
  plot_list <- c(list(ggplot() + theme_void()), col_headers)
  for (n in conditions) {
    row_panels <- lapply(bins, function(b) make_diff_panel(diffs[[n]][[b]]))
    label_p <- ggplot() +
      annotate("text", x=0.5, y=0.5,
               label=paste0(CONDITIONS[[n]]$label, "\nvs WT"),
               color=CONDITIONS[[n]]$color, size=3.5, fontface="bold", angle=90) +
      theme_void()
    plot_list <- c(plot_list, list(label_p), row_panels)
  }
  main_grid <- wrap_plots(plot_list, ncol=length(bins)+1,
                          widths=c(0.12, rep(1, length(bins))),
                          heights=c(0.15, rep(1, length(conditions))))
  legend_df <- data.frame(x=1:100, y=1, fill=seq(-limit, limit, length.out=100))
  legend_plot <- ggplot(legend_df, aes(x=x, y=y, fill=fill)) + geom_tile() +
    make_diff_scale(palette, limits=limit, name="Δlog2\nenrichment") +
    theme_void() +
    theme(legend.position="right", legend.text=element_text(size=8),
          legend.title=element_text(size=9, face="bold"))
  composite <- cowplot::plot_grid(main_grid, cowplot::get_legend(legend_plot),
                                  ncol=2, rel_widths=c(1, 0.1))
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0, y=0.75) +
    cowplot::draw_label(paste0("log2(treat/WT fold-enrichment)  |  ",
                               "Δlog2=central pixel  |  ",
                               RESOLUTION/1000, "kb  |  ±", BUFFER,
                               " bins  |  palette: ", palette),
                        size=8, color="gray40", x=0, hjust=0, y=0.25)
  final <- cowplot::plot_grid(title_grob, composite, ncol=1, rel_heights=c(0.08, 1))
  print(final); invisible(final)
}

# ══════════════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════════════
message("Loading expression data...")
expr      <- read_tsv(EXPR_FILE) %>% dplyr::select(gene_id, TPM)
deseq_res <- read_csv(DESEQ_FILE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
  left_join(expr, by="gene_id") %>%
  dplyr::mutate(
    de_status = dplyr::case_when(
      is.na(padj)                                               ~ "Filtered",
      padj >= PADJ_THRESHOLD                                    ~ "Not DE",
      padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD  ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD  ~ "Down",
      TRUE                                                      ~ "Not DE"
    ),
    de_status = factor(de_status, levels=c("Up","Not DE","Down","Filtered")),
    symbol    = dplyr::case_when(!is.na(symbol) & symbol != "" ~ symbol, TRUE ~ gene_id)
  )
message("DE breakdown:"); print(table(deseq_res$de_status))

expressed  <- deseq_res %>% filter(!is.na(TPM), TPM > 0) %>% pull(TPM)
tpm_medium <- quantile(expressed, 1/3, na.rm=TRUE)
tpm_high   <- quantile(expressed, 2/3, na.rm=TRUE)
deseq_res  <- deseq_res %>%
  dplyr::mutate(
    expr_bin = dplyr::case_when(
      is.na(TPM) | TPM == 0 ~ "Low",
      TPM <= tpm_medium     ~ "Low",
      TPM <= tpm_high       ~ "Medium",
      TRUE                  ~ "High"
    ),
    expr_bin = factor(expr_bin, levels=BINS)
  )
message("Expression bins:"); print(table(deseq_res$expr_bin))

gene_coords <- deseq_res %>%
  filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name)) %>%
  distinct(gene_id, .keep_all=TRUE)
gene_gr <- GRanges(
  seqnames = as.character(gene_coords$seq_name),
  ranges   = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end),
  symbol   = gene_coords$symbol,
  gene_id  = gene_coords$gene_id,
  de_status = as.character(gene_coords$de_status),
  expr_bin  = as.character(gene_coords$expr_bin),
  TPM       = gene_coords$TPM
)
active_gene_gr <- gene_gr[!is.na(gene_gr$TPM) & gene_gr$TPM > 0]
message("gene_gr: ", length(gene_gr), "  active (TPM>0): ", length(active_gene_gr))

# ── Loop loading ───────────────────────────────────────────────────────────────
message("Loading HiCCUPS loops (", DS_GROUP, ")...")
loops_raw <- lapply(names(CONDITIONS), function(n) {
  path <- CONDITIONS[[n]]$loops
  if (!file.exists(path)) stop("Loop file not found: ", path)
  filter_loops(load_hiccups_bedpe(path), min_observed=MIN_OBSERVED[[n]])
}); names(loops_raw) <- names(CONDITIONS)
invisible(lapply(names(loops_raw), function(n)
  message(CONDITIONS[[n]]$label, ": ", nrow(loops_raw[[n]]), " loops")))

message("Linking loops to genes...")
loops_annotated <- lapply(loops_raw, function(ldf) link_loops_to_genes(ldf, gene_gr))

# ── narrowPeak files ───────────────────────────────────────────────────────────
message("Loading narrowPeak files...")
atac_wt    <- load_narrowpeak(file.path(ATAC_DIR, "atac_wt.narrowPeak"))
h3k4me3_wt <- load_narrowpeak(file.path(CHIP_DIR, "H3K4me3.narrowPeak"))
h3k27ac_wt <- load_narrowpeak(file.path(CHIP_DIR, "h3k27ac.narrowPeak"))
h3k4me1_wt <- load_narrowpeak(file.path(CHIP_DIR, "H3K4me1.narrowPeak"))
ctcf_wt    <- load_narrowpeak(file.path(CHIP_DIR, "ctcf.narrowPeak"))
atac_ko    <- load_narrowpeak(file.path(ATAC_DIR, "atac_ko.narrowPeak"))
h3k4me3_ko <- load_narrowpeak(file.path(CHIP_DIR, "H3K4me3_ko.narrowPeak"))
h3k27ac_ko <- load_narrowpeak(file.path(CHIP_DIR, "h3k27ac_ko.narrowPeak"))
h3k4me1_ko <- load_narrowpeak(file.path(CHIP_DIR, "H3K4me1_ko.narrowPeak"))
ctcf_ko    <- load_narrowpeak(file.path(CHIP_DIR, "ctcf_ko.narrowPeak"))

# ── Expression-bin APA matrices ────────────────────────────────────────────────
message("Extracting APA matrices (expression bins)...")
apa_raw <- lapply(names(CONDITIONS), function(n) {
  message("--- ", CONDITIONS[[n]]$label, " ---")
  extract_apa_matrices(loops_annotated[[n]], CONDITIONS[[n]]$hic)
}); names(apa_raw) <- names(CONDITIONS)

apa_agg <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 1: All loops — expression bins
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 1: All loops (unmatched) ───────────────────────────────")
#apa.wtpol2.expr
plot_apa_grid_generic(apa_agg, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA by expression level — all loops")
plot_apa_scores(apa_agg,
                loops_annotated_cond = loops_annotated,
                group_col="expr_bin", de_filter=NULL,
                title="APA — all loops")
plot_apa_diff(apa_agg, bins=BINS,
              title="APA difference — all loops")


message("\n── Analysis 1 (matched): WT loop positions × all conditions ─────────")
apa_raw_1m <- lapply(names(CONDITIONS), function(n) {
  message("Pulling WT positions from ", CONDITIONS[[n]]$label, "...")
  extract_apa_matrices(loops_annotated$wt, CONDITIONS[[n]]$hic)
}); names(apa_raw_1m) <- names(CONDITIONS)
apa_agg_1m <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_1m[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
plot_apa_grid_generic(apa_agg_1m, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — all loops  |  WT positions in all conditions")
plot_apa_scores(apa_agg_1m, loops_annotated_cond=loops_annotated$wt,
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — all loops  |  WT positions matched")
plot_apa_diff(apa_agg_1m, bins=BINS,
              title="APA difference — all loops  |  WT positions, treatment vs WT")

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 1c: Exclusive loops — WT-only and pol2-only
#   WT-only loops:   scored at WT .hic and pol2 .hic at WT positions
#   pol2-only loops: scored at pol2 .hic and WT .hic at pol2 positions
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Finding WT-only and pol2-only exclusive loops ───────────────────")
EXCL_SLOP <- 5000L ## loop tolerance

make_gr_pair <- function(df) list(
  a1 = GRanges(df$chr1, IRanges(df$start1, df$end1)),
  a2 = GRanges(df$chr2, IRanges(df$start2, df$end2))
)
exp_gr <- function(gr) resize(gr, width(gr) + 2*EXCL_SLOP, fix="center")

gr_wt   <- make_gr_pair(loops_annotated$wt)
gr_pol2 <- make_gr_pair(loops_annotated$pol2)

hits1 <- findOverlaps(exp_gr(gr_wt$a1), exp_gr(gr_pol2$a1))
hits2 <- findOverlaps(exp_gr(gr_wt$a2), exp_gr(gr_pol2$a2))
shared_idx <- dplyr::inner_join(
  as.data.frame(hits1) %>% dplyr::rename(w=queryHits, p=subjectHits),
  as.data.frame(hits2) %>% dplyr::rename(w=queryHits, p=subjectHits),
  by=c("w","p")
) %>% dplyr::distinct(w, p)

wt_only_idx   <- setdiff(seq_len(nrow(loops_annotated$wt)),   shared_idx$w)
pol2_only_idx <- setdiff(seq_len(nrow(loops_annotated$pol2)), shared_idx$p)

wt_only_loops   <- loops_annotated$wt   %>% dplyr::slice(wt_only_idx)
pol2_only_loops <- loops_annotated$pol2 %>% dplyr::slice(pol2_only_idx)

message("Shared:          ", nrow(shared_idx))
message("WT-only loops:   ", nrow(wt_only_loops))
message("pol2-only loops: ", nrow(pol2_only_loops))

# ── WT-only: WT positions scored in both .hic files ──────────────────────────
message("\nExtracting APA — WT-only loops at WT positions in all conditions...")
#apa.wtpol2.expr.matched.wtloopsonly
apa_raw_wt_only <- lapply(names(CONDITIONS), function(n) {
  message("  ", CONDITIONS[[n]]$label, "...")
  extract_apa_matrices(wt_only_loops, CONDITIONS[[n]]$hic)
}); names(apa_raw_wt_only) <- names(CONDITIONS)

apa_agg_wt_only <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_wt_only[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)

plot_apa_grid_generic(apa_agg_wt_only, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — WT-only loops  |  WT positions scored in all conditions")
plot_apa_scores(apa_agg_wt_only, loops_annotated_cond=wt_only_loops,
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — WT-only loops  |  WT positions in all conditions")
plot_apa_diff(apa_agg_wt_only, bins=BINS,
              title="APA difference — WT-only loops  |  treatment vs WT at WT-only positions")

# ── pol2-only: pol2 positions scored in both .hic files ──────────────────────
message("\nExtracting APA — pol2-only loops at pol2 positions in all conditions...")
apa_raw_pol2_only <- lapply(names(CONDITIONS), function(n) {
  message("  ", CONDITIONS[[n]]$label, "...")
  extract_apa_matrices(pol2_only_loops, CONDITIONS[[n]]$hic)
}); names(apa_raw_pol2_only) <- names(CONDITIONS)

apa_agg_pol2_only <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_pol2_only[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)

plot_apa_grid_generic(apa_agg_pol2_only, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — pol2-only loops  |  pol2 positions scored in all conditions")
plot_apa_scores(apa_agg_pol2_only, loops_annotated_cond=pol2_only_loops,
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — pol2-only loops  |  pol2 positions in all conditions")
plot_apa_diff(apa_agg_pol2_only, bins=BINS,
              title="APA difference — pol2-only loops  |  treatment vs WT at pol2-only positions")

# ── Shared loops: scored in both .hic files at shared WT positions ────────────
message("\nExtracting APA — shared loops at WT positions in all conditions...")
shared_loops <- loops_annotated$wt %>% dplyr::slice(shared_idx$w)
message("Shared loops: ", nrow(shared_loops))

apa_raw_shared <- lapply(names(CONDITIONS), function(n) {
  message("  ", CONDITIONS[[n]]$label, "...")
  extract_apa_matrices(shared_loops, CONDITIONS[[n]]$hic)
}); names(apa_raw_shared) <- names(CONDITIONS)

apa_agg_shared <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_shared[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)

plot_apa_grid_generic(apa_agg_shared, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — shared loops  |  WT positions scored in all conditions")
plot_apa_scores(apa_agg_shared, loops_annotated_cond=shared_loops,
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — shared loops  |  WT positions in all conditions")
plot_apa_diff(apa_agg_shared, bins=BINS,
              title="APA difference — shared loops  |  treatment vs WT at shared positions")


# ══════════════════════════════════════════════════════════════════════════════
# 3x3 matched grid, single bin — rows 1&2 (WT-scored / Pol2-scored) share ONE
# fold-enrichment scale + ONE legend; row 3 (diff) keeps its own legend.
# ══════════════════════════════════════════════════════════════════════════════
plot_apa_matched_grid <- function(bin,
                                  apa_agg_wt_only,
                                  apa_agg_pol2_only,
                                  apa_agg_shared,
                                  bin_color    = "red",
                                  palette      = ACTIVE_PALETTE,
                                  diff_palette = ACTIVE_DIFF_PALETTE,
                                  corner_size  = 3,
                                  pseudocount  = 0.01,
                                  title        = NULL) {
  
  normalize_to_bg <- function(mat) {
    if (is.null(mat)) return(NULL)
    n <- nrow(mat)
    corners <- c(as.vector(mat[1:corner_size, 1:corner_size]),
                 as.vector(mat[1:corner_size, (n-corner_size+1):n]),
                 as.vector(mat[(n-corner_size+1):n, 1:corner_size]),
                 as.vector(mat[(n-corner_size+1):n, (n-corner_size+1):n]))
    cm <- mean(corners[!is.na(corners)], na.rm=TRUE)
    if (is.na(cm) || cm == 0) return(mat)
    mat / cm
  }
  
  axis_breaks <- c(-BUFFER, 0, BUFFER)
  axis_labs   <- c(paste0("-", BUFFER*RESOLUTION/1000, "kb"), "loop",
                   paste0("+", BUFFER*RESOLUTION/1000, "kb"))
  
  make_fold_panel <- function(mat_norm, global_min, global_max) {
    if (is.null(mat_norm)) return(ggplot() + theme_void() +
                                    theme(plot.background=element_rect(fill="gray95", color=NA)))
    n <- nrow(mat_norm); mid <- ceiling(n/2); score <- mat_norm[mid, mid]
    df <- expand.grid(row=1:n, col=1:n) %>%
      dplyr::mutate(value=as.vector(mat_norm), x=col-mid, y=-(row-mid))
    p <- ggplot(df, aes(x=x, y=y, fill=pmin(pmax(value, global_min), global_max))) +
      geom_tile() +
      make_fill_scale(palette, limits=c(global_min, global_max), name="Fold over\nbackground") +
      annotate("rect", xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5,
               fill=NA, color="black", linewidth=0.6) +
      coord_fixed() +
      scale_x_continuous(breaks=axis_breaks, labels=axis_labs) +
      scale_y_continuous(breaks=axis_breaks, labels=axis_labs) +
      labs(x=NULL, y=NULL) + theme_classic(base_size=8) +
      theme(axis.text=element_text(size=9), axis.line=element_blank(),
            axis.ticks=element_blank(), legend.position="none",
            plot.margin=margin(2,2,2,2))
    if (!is.na(score))
      p <- p + annotate("text", x=-BUFFER+0.5, y=BUFFER-0.5,
                        label=paste0("APA=", round(score, 2)),
                        size=2.5, hjust=0, fontface="bold", color="white")
    p
  }
  
  make_diff_panel <- function(mat_diff, limit) {
    if (is.null(mat_diff)) return(ggplot() + theme_void() +
                                    theme(plot.background=element_rect(fill="gray95", color=NA)))
    n <- nrow(mat_diff); mid <- ceiling(n/2); score <- mat_diff[mid, mid]
    df <- expand.grid(row=1:n, col=1:n) %>%
      dplyr::mutate(value=as.vector(mat_diff), x=col-mid, y=-(row-mid))
    p <- ggplot(df, aes(x=x, y=y, fill=pmin(pmax(value, -limit), limit))) +
      geom_tile() + make_diff_scale(diff_palette, limits=limit) +
      annotate("rect", xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5,
               fill=NA, color="black", linewidth=0.6) +
      coord_fixed() +
      scale_x_continuous(breaks=axis_breaks, labels=axis_labs) +
      scale_y_continuous(breaks=axis_breaks, labels=axis_labs) +
      labs(x=NULL, y=NULL) + theme_classic(base_size=8) +
      theme(axis.text=element_text(size=9), axis.line=element_blank(),
            axis.ticks=element_blank(), legend.position="none",
            plot.margin=margin(2,2,2,2))
    if (!is.na(score))
      p <- p + annotate("text", x=-BUFFER+0.5, y=BUFFER-0.5,
                        label=paste0("\u0394log2=", round(score, 2)),
                        size=2.5, hjust=0, fontface="bold")
    p
  }
  
  columns   <- list("WT-only loops"=apa_agg_wt_only,
                    "Pol2-only loops"=apa_agg_pol2_only,
                    "shared loops"=apa_agg_shared)
  col_names <- names(columns)
  
  wt_mats   <- lapply(columns, function(a) normalize_to_bg(a$wt[[bin]]))
  pol2_mats <- lapply(columns, function(a) normalize_to_bg(a$pol2[[bin]]))
  diff_mats <- setNames(lapply(col_names, function(n) {
    if (is.null(wt_mats[[n]]) || is.null(pol2_mats[[n]])) return(NULL)
    log2((pol2_mats[[n]] + pseudocount) / (wt_mats[[n]] + pseudocount))
  }), col_names)
  
  # ── Shared fold-enrichment scale across BOTH row 1 and row 2 ────────────────
  fold_vals <- unlist(lapply(c(wt_mats, pol2_mats), function(m) if (!is.null(m)) as.vector(m)))
  fold_max  <- quantile(fold_vals, 0.98, na.rm=TRUE)
  fold_min  <- 1.0
  
  diff_vals <- unlist(lapply(diff_mats, function(m) if (!is.null(m)) as.vector(m)))
  diff_lim  <- quantile(abs(diff_vals[!is.na(diff_vals)]), 0.98, na.rm=TRUE)
  
  header_plots  <- lapply(col_names, function(n)
    ggplot() + annotate("text", x=0.5, y=0.5, label=n, size=3.8, fontface="bold") + theme_void())
  bin_hdr_plots <- lapply(col_names, function(n)
    ggplot() + annotate("text", x=0.5, y=0.5, label=bin, color=bin_color,
                        size=3.5, fontface="bold") + theme_void())
  
  build_fold_legend <- function(gmin, gmax) cowplot::get_legend(
    ggplot(data.frame(x=1:10, y=1, fill=seq(gmin, gmax, length.out=10)),
           aes(x=x, y=y, fill=fill)) + geom_tile() +
      make_fill_scale(palette, limits=c(gmin, gmax), name="Fold over\nbackground") +
      theme(legend.position="right", legend.key.height=unit(1.4,"cm"),
            legend.key.width=unit(0.3,"cm"), legend.text=element_text(size=8),
            legend.title=element_text(size=9)))
  
  build_diff_legend <- function(limit) cowplot::get_legend(
    ggplot(data.frame(x=1:10, y=1, fill=seq(-limit, limit, length.out=10)),
           aes(x=x, y=y, fill=fill)) + geom_tile() +
      make_diff_scale(diff_palette, limits=limit, name="\u0394log2\nenrichment") +
      theme(legend.position="right", legend.key.height=unit(1.4,"cm"),
            legend.key.width=unit(0.3,"cm"), legend.text=element_text(size=8),
            legend.title=element_text(size=9)))
  
  row_wt_panels   <- lapply(col_names, function(n) make_fold_panel(wt_mats[[n]],   fold_min, fold_max))
  row_pol2_panels <- lapply(col_names, function(n) make_fold_panel(pol2_mats[[n]], fold_min, fold_max))
  row_diff_panels <- lapply(col_names, function(n) make_diff_panel(diff_mats[[n]], diff_lim))
  
  row_wt   <- cowplot::plot_grid(plotlist=row_wt_panels,   ncol=3)
  row_pol2 <- cowplot::plot_grid(plotlist=row_pol2_panels, ncol=3)
  row_diff <- cowplot::plot_grid(cowplot::plot_grid(plotlist=row_diff_panels, ncol=3),
                                 build_diff_legend(diff_lim), ncol=2, rel_widths=c(1, 0.15))
  
  # ── Rows 1+2 share one legend, vertically spanning both rows ───────────────
  rows_1_2_body <- cowplot::plot_grid(row_wt, row_pol2, ncol=1)
  rows_1_2      <- cowplot::plot_grid(rows_1_2_body, build_fold_legend(fold_min, fold_max),
                                      ncol=2, rel_widths=c(1, 0.15))
  
  header_row  <- cowplot::plot_grid(cowplot::plot_grid(plotlist=header_plots,  ncol=3),
                                    ggplot()+theme_void(), ncol=2, rel_widths=c(1, 0.15))
  bin_hdr_row <- cowplot::plot_grid(cowplot::plot_grid(plotlist=bin_hdr_plots, ncol=3),
                                    ggplot()+theme_void(), ncol=2, rel_widths=c(1, 0.15))
  
  composite <- cowplot::plot_grid(header_row, bin_hdr_row, rows_1_2, row_diff,
                                  ncol=1, rel_heights=c(0.08, 0.08, 2, 1))
  
  if (!is.null(title)) {
    title_grob <- cowplot::ggdraw() +
      cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0)
    composite <- cowplot::plot_grid(title_grob, composite, ncol=1, rel_heights=c(0.06, 1))
  }
  
  print(composite); invisible(composite)
}

# ══════════════════════════════════════════════════════════════════════════════
# Wrapper: one 3x3 matched grid per expression bin (3 separate figures)
# ══════════════════════════════════════════════════════════════════════════════
plot_apa_matched_grid_by_expr <- function(bins = c("Low","Medium","High"),
                                          apa_agg_wt_only,
                                          apa_agg_pol2_only,
                                          apa_agg_shared,
                                          bin_colors  = BIN_COLORS,
                                          out_prefix  = "apa_matched_grid",
                                          out_dir     = ".",
                                          width       = 11,
                                          height      = 9,
                                          save_files  = TRUE) {
  
  results <- setNames(vector("list", length(bins)), bins)
  
  for (b in bins) {
    message("Building matched grid — ", b, " expression...")
    
    composite <- plot_apa_matched_grid(
      bin                = b,
      apa_agg_wt_only    = apa_agg_wt_only,
      apa_agg_pol2_only  = apa_agg_pol2_only,
      apa_agg_shared     = apa_agg_shared,
      bin_color          = bin_colors[[b]],
      title              = paste0("APA \u2014 ", b, "-expression loops: matched WT vs Pol2 scoring")
    )
    
    results[[b]] <- composite
    
    if (save_files) {
      out_path <- file.path(out_dir, paste0(out_prefix, "_", tolower(b), ".pdf"))
      ggsave(out_path, composite, width = width, height = height, limitsize = FALSE)
      message("  Saved: ", out_path)
    }
  }
  
  invisible(results)
}

# ── Usage — produces 3 separate 3x3 figures, one per bin ─────────────────────
grids <- plot_apa_matched_grid_by_expr(
  bins               = c("Low","Medium","High"),
  apa_agg_wt_only    = apa_agg_wt_only,
  apa_agg_pol2_only  = apa_agg_pol2_only,
  apa_agg_shared     = apa_agg_shared,
  out_prefix         = "apa_matched_grid",
  out_dir            = "."
)
# grids$Low, grids$Medium, grids$High hold the ggplot/cowplot objects if you want
# to inspect or re-save any of them individually

# ══════════════════════════════════════════════════════════════════════════════
# Trajectory plots (APA score across P-P/E-P/E-E + log2FC vs WT) for each
# loop-subset category — reuses plot_apa_scores() as-is, no new logic.
# ══════════════════════════════════════════════════════════════════════════════
plot_apa_trajectory_by_ccre_subset <- function(apa_agg_ccre_wt_only,   wt_only_ccre,
                                               apa_agg_ccre_pol2_only, pol2_only_ccre,
                                               apa_agg_ccre_shared,    shared_ccre_wt,
                                               bins       = GROUP_LEVELS,
                                               out_prefix = "apa_trajectory_ccre",
                                               out_dir    = ".",
                                               width      = 8, height = 9,
                                               save_files = TRUE) {
  
  subsets <- list(
    wt_only = list(
      agg   = apa_agg_ccre_wt_only,
      loops = wt_only_ccre,
      title = "APA scores \u2014 WT-only cCRE loops  |  WT positions in all conditions"
    ),
    pol2_only = list(
      agg   = apa_agg_ccre_pol2_only,
      loops = pol2_only_ccre,
      title = "APA scores \u2014 pol2-only cCRE loops  |  pol2 positions in all conditions"
    ),
    shared = list(
      agg   = apa_agg_ccre_shared,
      loops = shared_ccre_wt,
      title = "APA scores \u2014 shared cCRE loops  |  WT positions in all conditions"
    )
  )
  
  results <- setNames(vector("list", length(subsets)), names(subsets))
  
  for (n in names(subsets)) {
    s <- subsets[[n]]
    message("Building cCRE trajectory plot \u2014 ", n, "...")
    
    out <- plot_apa_scores(
      apa_agg              = s$agg,
      loops_annotated_cond = s$loops,
      group_col            = "group",
      bins                 = bins,
      de_filter            = NULL,
      title                = s$title
    )
    
    results[[n]] <- out
    
    if (save_files) {
      out_path <- file.path(out_dir, paste0(out_prefix, "_", n, ".pdf"))
      ggsave(out_path, out$composite, width = width, height = height, limitsize = FALSE)
      message("  Saved: ", out_path)
    }
  }
  
  invisible(results)
}

# ── Usage — mirrors plot_apa_trajectory_by_subset, keyed on cCRE group ───────
ccre_trajectories <- plot_apa_trajectory_by_ccre_subset(
  apa_agg_ccre_wt_only   = apa_agg_ccre_wt_only,   wt_only_ccre   = wt_only_ccre,
  apa_agg_ccre_pol2_only = apa_agg_ccre_pol2_only, pol2_only_ccre = pol2_only_ccre,
  apa_agg_ccre_shared    = apa_agg_ccre_shared,    shared_ccre_wt = shared_ccre_wt
)
# ccre_trajectories$wt_only$composite, $pol2_only$composite, $shared$composite
# ccre_trajectories$<name>$data / $fc_data for the raw tables

# ══════════════════════════════════════════════════════════════════════════════
# Trajectory plots (P2LL across expression bins + log2FC vs WT) for each
# loop-subset category — reuses plot_apa_scores() as-is, no new logic.
# ══════════════════════════════════════════════════════════════════════════════
plot_apa_trajectory_by_subset <- function(apa_agg_wt_only,   wt_only_loops,
                                          apa_agg_pol2_only, pol2_only_loops,
                                          apa_agg_shared,    shared_loops,
                                          bins       = BINS,
                                          out_prefix = "apa_trajectory",
                                          out_dir    = ".",
                                          width      = 8, height = 9,
                                          save_files = TRUE) {
  
  subsets <- list(
    wt_only = list(
      agg   = apa_agg_wt_only,
      loops = wt_only_loops,
      title = "APA scores \u2014 WT-only loops  |  WT positions in all conditions"
    ),
    pol2_only = list(
      agg   = apa_agg_pol2_only,
      loops = pol2_only_loops,
      title = "APA scores \u2014 pol2-only loops  |  pol2 positions in all conditions"
    ),
    shared = list(
      agg   = apa_agg_shared,
      loops = shared_loops,
      title = "APA scores \u2014 shared loops  |  WT positions in all conditions"
    )
  )
  
  results <- setNames(vector("list", length(subsets)), names(subsets))
  
  for (n in names(subsets)) {
    s <- subsets[[n]]
    message("Building trajectory plot \u2014 ", n, "...")
    
    out <- plot_apa_scores(
      apa_agg              = s$agg,
      loops_annotated_cond = s$loops,
      group_col            = "expr_bin",
      bins                 = bins,
      de_filter            = NULL,
      title                = s$title
    )
    
    results[[n]] <- out
    
    if (save_files) {
      out_path <- file.path(out_dir, paste0(out_prefix, "_", n, ".pdf"))
      ggsave(out_path, out$composite, width = width, height = height, limitsize = FALSE)
      message("  Saved: ", out_path)
    }
  }
  
  invisible(results)
}

# ── Usage — reuses the exact objects your matched-grid wrapper already uses ──
trajectories <- plot_apa_trajectory_by_subset(
  apa_agg_wt_only   = apa_agg_wt_only,   wt_only_loops   = wt_only_loops,
  apa_agg_pol2_only = apa_agg_pol2_only, pol2_only_loops = pol2_only_loops,
  apa_agg_shared    = apa_agg_shared,    shared_loops    = shared_loops
)
# trajectories$wt_only$composite, trajectories$pol2_only$composite,
# trajectories$shared$composite hold the individual patchwork objects,
# and trajectories$<name>$data / $fc_data give you the raw score / log2FC
# tables if you want them in a table instead of a plot.

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 2: DE gene subsets
# ══════════════════════════════════════════════════════════════════════════════
extract_apa_de <- function(loops_annotated_list, de_direction) {
  result <- lapply(names(CONDITIONS), function(n) {
    loops_de <- loops_annotated_list[[n]] %>% filter(de_status == de_direction)
    message(CONDITIONS[[n]]$label, " (", de_direction, "): ", nrow(loops_de), " loops")
    extract_apa_matrices(loops_de, CONDITIONS[[n]]$hic)
  }); setNames(result, names(CONDITIONS))
}

message("\n── Analysis 2a: Upregulated genes (unmatched) ───────────────────────")
# apa.wtpol2.up.deg
apa_raw_up <- extract_apa_de(loops_annotated, "Up")
apa_agg_up <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_up[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
plot_apa_grid_generic(apa_agg_up, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — upregulated genes")
plot_apa_scores(apa_agg_up,
                loops_annotated_cond = loops_annotated,
                group_col="expr_bin", de_filter="Up",
                title="APA — upregulated genes")


message("\n── Analysis 2a (matched): WT upregulated positions × all conditions ─")
wt_loops_up <- loops_annotated$wt %>% filter(de_status == "Up")
message("WT upregulated loops: ", nrow(wt_loops_up))
apa_raw_up_m <- lapply(names(CONDITIONS), function(n) {
  extract_apa_matrices(wt_loops_up, CONDITIONS[[n]]$hic)
}); names(apa_raw_up_m) <- names(CONDITIONS)
apa_agg_up_m <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_up_m[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
plot_apa_grid_generic(apa_agg_up_m, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — upregulated genes  |  WT positions matched")
plot_apa_scores(apa_agg_up_m, loops_annotated_cond=wt_loops_up,
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — upregulated genes  |  WT positions matched")
plot_apa_diff(apa_agg_up_m, bins=BINS,
              title="APA difference — upregulated genes  |  WT positions, treatment vs WT")

message("\n── Analysis 2b: Downregulated genes (unmatched) ─────────────────────")
apa_raw_down <- extract_apa_de(loops_annotated, "Down")
apa_agg_down <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_down[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
plot_apa_grid_generic(apa_agg_down, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — downregulated genes")
plot_apa_scores(apa_agg_down, loops_annotated_cond=loops_annotated,
                group_col="expr_bin", de_filter="Down", title="APA — downregulated genes")

message("\n── Analysis 2b (matched): WT downregulated positions × all conditions")
wt_loops_down <- loops_annotated$wt %>% filter(de_status == "Down")
message("WT downregulated loops: ", nrow(wt_loops_down))
apa_raw_down_m <- lapply(names(CONDITIONS), function(n) {
  extract_apa_matrices(wt_loops_down, CONDITIONS[[n]]$hic)
}); names(apa_raw_down_m) <- names(CONDITIONS)
apa_agg_down_m <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_down_m[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
plot_apa_grid_generic(apa_agg_down_m, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — downregulated genes  |  WT positions matched")
plot_apa_scores(apa_agg_down_m, loops_annotated_cond=wt_loops_down,
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — downregulated genes  |  WT positions matched")
plot_apa_diff(apa_agg_down_m, bins=BINS,
              title="APA difference — downregulated genes  |  WT positions, treatment vs WT")

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 2c: Exclusive DE loops
#   WT-only loops at Up/Down genes → extracted from WT .hic
#   pol2-only loops at Up/Down genes → extracted from pol2 .hic
#   Requires wt_only_loops and pol2_only_loops from Analysis 1c
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 2c: Exclusive DE loops ──────────────────────────────────")

wt_only_up    <- wt_only_loops   %>% filter(de_status == "Up")
wt_only_down  <- wt_only_loops   %>% filter(de_status == "Down")
pol2_only_up  <- pol2_only_loops %>% filter(de_status == "Up")
pol2_only_down <- pol2_only_loops %>% filter(de_status == "Down")

message("WT-only Up:     ", nrow(wt_only_up))
message("WT-only Down:   ", nrow(wt_only_down))
message("pol2-only Up:   ", nrow(pol2_only_up))
message("pol2-only Down: ", nrow(pol2_only_down))

# ── Extract APA — each subset at its own .hic ─────────────────────────────
message("\nExtracting APA — exclusive Up loops...")
apa_raw_excl_up <- list(
  wt   = extract_apa_matrices(wt_only_up,   CONDITIONS$wt$hic),
  pol2 = extract_apa_matrices(pol2_only_up, CONDITIONS$pol2$hic)
)
apa_agg_excl_up <- list(
  wt   = aggregate_apa(apa_raw_excl_up$wt,   DEPTH_SCALE$wt),
  pol2 = aggregate_apa(apa_raw_excl_up$pol2, DEPTH_SCALE$pol2)
)

message("\nExtracting APA — exclusive Down loops...")
apa_raw_excl_down <- list(
  wt   = extract_apa_matrices(wt_only_down,   CONDITIONS$wt$hic),
  pol2 = extract_apa_matrices(pol2_only_down, CONDITIONS$pol2$hic)
)
apa_agg_excl_down <- list(
  wt   = aggregate_apa(apa_raw_excl_down$wt,   DEPTH_SCALE$wt),
  pol2 = aggregate_apa(apa_raw_excl_down$pol2, DEPTH_SCALE$pol2)
)

# ── APA grid ─────────────────────────────────────────────────────────────────
plot_apa_grid_generic(apa_agg_excl_up, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — exclusive Up loops  |  WT-only (WT .hic) vs pol2-only (pol2 .hic)")
plot_apa_grid_generic(apa_agg_excl_down, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — exclusive Down loops  |  WT-only (WT .hic) vs pol2-only (pol2 .hic)")

# ── Score plot ────────────────────────────────────────────────────────────────
plot_apa_scores(apa_agg_excl_up,
                loops_annotated_cond = list(wt=wt_only_up, pol2=pol2_only_up),
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — exclusive Up loops")
plot_apa_scores(apa_agg_excl_down,
                loops_annotated_cond = list(wt=wt_only_down, pol2=pol2_only_down),
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — exclusive Down loops")

# ── Difference plot ───────────────────────────────────────────────────────────
plot_apa_diff(apa_agg_excl_up, bins=BINS,
              title="APA difference — exclusive Up loops  |  pol2-only vs WT-only")
plot_apa_diff(apa_agg_excl_down, bins=BINS,
              title="APA difference — exclusive Down loops  |  pol2-only vs WT-only")


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 3: P-P / E-P / E-E — condition-matched narrowPeak annotation
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 3: P-P / E-P / E-E — narrowPeak + ATAC ──────────────────")
message("Classifying WT loops (WT peaks)...")

library(AnnotationDbi)
txdb <- loadDb("/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/GRCh38.115.sqlite")

wt_np_class <- classify_anchors_narrowpeak(
  loops_df   = loops_raw$wt,
  atac_gr    = atac_wt,
  h3k4me3_gr = h3k4me3_wt,
  h3k27ac_gr = h3k27ac_wt,
  h3k4me1_gr = h3k4me1_wt,
  ctcf_gr    = ctcf_wt,
  txdb       = txdb,
  tss_region = c(-3000, 3000),
  ccre_label = "WT_narrowpeak"
)

pol2_np_class <- classify_anchors_narrowpeak(
  loops_df   = loops_raw$pol2,
  atac_gr    = atac_ko,
  h3k4me3_gr = h3k4me3_ko,
  h3k27ac_gr = h3k27ac_ko,
  h3k4me1_gr = h3k4me1_ko,
  ctcf_gr    = ctcf_ko,
  txdb       = txdb,
  tss_region = c(-3000, 3000),
  ccre_label = "KO_narrowpeak"
)

np_classified <- list(wt=wt_np_class, pol2=pol2_np_class)

# If actd present, classify with WT peaks
if ("actd" %in% names(CONDITIONS)) {
  message("Classifying actd loops (WT peaks — no ActD ChIP available)...")
  actd_np_class <- classify_anchors_narrowpeak(
    loops_df   = loops_raw$actd,
    atac_gr    = atac_wt,
    h3k4me3_gr = h3k4me3_wt,
    h3k27ac_gr = h3k27ac_wt,
    h3k4me1_gr = h3k4me1_wt,
    ctcf_gr    = ctcf_wt,
    txdb       = txdb,
    tss_region = c(-3000, 3000),
    ccre_label = "WT_narrowpeak_actd"
  )
  np_classified$actd <- actd_np_class
}

np_grouped <- lapply(names(np_classified), function(n) {
  np_classified[[n]] %>%
    dplyr::mutate(
      group = dplyr::case_when(
        loop_ccre_type %in% names(GROUP_MAP) ~ GROUP_MAP[loop_ccre_type],
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group)) %>%
    dplyr::mutate(group = factor(group, levels=GROUP_LEVELS))
}); names(np_grouped) <- names(np_classified)

invisible(lapply(names(np_grouped), function(n) {
  message(CONDITIONS[[n]]$label, " group breakdown:")
  print(table(np_grouped[[n]]$group))
}))

message("\n── Analysis 3 (unmatched): each condition's own loops ───────────────")
apa_raw_ccre <- lapply(names(CONDITIONS), function(n) {
  message("Extracting ", CONDITIONS[[n]]$label, " P-P/E-P/E-E APA...")
  extract_apa_matrices_by(np_grouped[[n]], CONDITIONS[[n]]$hic,
                          group_col="group", group_levels=GROUP_LEVELS)
}); names(apa_raw_ccre) <- names(CONDITIONS)
apa_agg_ccre <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_ccre[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
#apa.ccre.wtpol2
plot_apa_grid_generic(apa_agg_ccre, bin_levels=GROUP_LEVELS,
                      bin_colors=as.list(GROUP_COLORS),
                      title="APA — P-P/E-P/E-E  |  narrowPeak+ATAC  |  unmatched")
plot_apa_scores(apa_agg_ccre, loops_annotated_cond=np_grouped,
                group_col="group", bins=GROUP_LEVELS,
                title="APA scores — P-P/E-P/E-E (unmatched)")

message("\n── Analysis 3 (matched): WT P-P/E-P/E-E positions × all conditions ──")
apa_raw_ccre_m <- lapply(names(CONDITIONS), function(n) {
  message("Pulling WT P-P/E-P/E-E positions from ", CONDITIONS[[n]]$label, "...")
  extract_apa_matrices_by(np_grouped$wt, CONDITIONS[[n]]$hic,
                          group_col="group", group_levels=GROUP_LEVELS)
}); names(apa_raw_ccre_m) <- names(CONDITIONS)
apa_agg_ccre_m <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_ccre_m[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
# apa.ccre.pol2wt.matched
plot_apa_grid_generic(apa_agg_ccre_m, bin_levels=GROUP_LEVELS,
                      bin_colors=as.list(GROUP_COLORS),
                      title="APA — P-P/E-P/E-E  |  WT positions matched")
plot_apa_scores(apa_agg_ccre_m, loops_annotated_cond=np_grouped,
                group_col="group", bins=GROUP_LEVELS,
                title="APA scores — P-P/E-P/E-E  |  WT positions matched")
plot_apa_diff(apa_agg_ccre_m, bins=GROUP_LEVELS,
              title="APA difference — P-P/E-P/E-E  |  WT positions, treatment vs WT")

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 3c: Exclusive cCRE loops — WT-only and pol2-only P-P/E-P/E-E
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 3c: Exclusive cCRE loops ────────────────────────────────")
EXCL_SLOP=5000

# Find shared/exclusive within the grouped (P-P/E-P/E-E) loop sets
make_gr_pair <- function(df) list(
  a1 = GRanges(df$chr1, IRanges(df$start1, df$end1)),
  a2 = GRanges(df$chr2, IRanges(df$start2, df$end2))
)
exp_gr  <- function(gr) resize(gr, width(gr) + 2*EXCL_SLOP, fix="center")

gr_wt_ccre   <- make_gr_pair(np_grouped$wt)
gr_pol2_ccre <- make_gr_pair(np_grouped$pol2)

hits1_ccre <- findOverlaps(exp_gr(gr_wt_ccre$a1), exp_gr(gr_pol2_ccre$a1))
hits2_ccre <- findOverlaps(exp_gr(gr_wt_ccre$a2), exp_gr(gr_pol2_ccre$a2))
shared_ccre <- dplyr::inner_join(
  as.data.frame(hits1_ccre) %>% dplyr::rename(w=queryHits, p=subjectHits),
  as.data.frame(hits2_ccre) %>% dplyr::rename(w=queryHits, p=subjectHits),
  by=c("w","p")
) %>% dplyr::distinct(w, p)

wt_only_ccre   <- np_grouped$wt   %>% dplyr::slice(setdiff(seq_len(nrow(np_grouped$wt)),   shared_ccre$w))
pol2_only_ccre <- np_grouped$pol2 %>% dplyr::slice(setdiff(seq_len(nrow(np_grouped$pol2)), shared_ccre$p))

message("Shared cCRE loops:    ", nrow(shared_ccre))
message("WT-only cCRE loops:   ", nrow(wt_only_ccre))
message("pol2-only cCRE loops: ", nrow(pol2_only_ccre))
message("WT-only group breakdown:");   print(table(wt_only_ccre$group))
message("pol2-only group breakdown:"); print(table(pol2_only_ccre$group))

# ── Extract APA — each exclusive set at its own .hic ─────────────────────────
message("\nExtracting APA — WT-only cCRE loops...")
#apa.ccre.wtpol2

apa_raw_ccre_wt_only <- list(
  wt   = extract_apa_matrices_by(wt_only_ccre, CONDITIONS$wt$hic,
                                 group_col="group", group_levels=GROUP_LEVELS),
  pol2 = extract_apa_matrices_by(wt_only_ccre, CONDITIONS$pol2$hic,
                                 group_col="group", group_levels=GROUP_LEVELS)
)
apa_agg_ccre_wt_only <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_ccre_wt_only[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)

message("\nExtracting APA — pol2-only cCRE loops...")
apa_raw_ccre_pol2_only <- list(
  wt   = extract_apa_matrices_by(pol2_only_ccre, CONDITIONS$wt$hic,
                                 group_col="group", group_levels=GROUP_LEVELS),
  pol2 = extract_apa_matrices_by(pol2_only_ccre, CONDITIONS$pol2$hic,
                                 group_col="group", group_levels=GROUP_LEVELS)
)
apa_agg_ccre_pol2_only <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_ccre_pol2_only[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)

# ── APA grid ──────────────────────────────────────────────────────────────────
#apa.ccre.wtpol2.wtloopsonly
plot_apa_grid_generic(apa_agg_ccre_wt_only, bin_levels=GROUP_LEVELS,
                      bin_colors=as.list(GROUP_COLORS),
                      title="APA — WT-only cCRE loops  |  WT positions in all conditions")
#apa.ccre.wtpol2.pol2loopsonly
plot_apa_grid_generic(apa_agg_ccre_pol2_only, bin_levels=GROUP_LEVELS,
                      bin_colors=as.list(GROUP_COLORS),
                      title="APA — pol2-only cCRE loops  |  pol2 positions in all conditions")

# ── Score plot ────────────────────────────────────────────────────────────────
plot_apa_scores(apa_agg_ccre_wt_only,
                loops_annotated_cond = list(wt=wt_only_ccre, pol2=wt_only_ccre),
                group_col="group", bins=GROUP_LEVELS,
                title="APA scores — WT-only cCRE loops  |  WT positions in all conditions")
plot_apa_scores(apa_agg_ccre_pol2_only,
                loops_annotated_cond = list(wt=pol2_only_ccre, pol2=pol2_only_ccre),
                group_col="group", bins=GROUP_LEVELS,
                title="APA scores — pol2-only cCRE loops  |  pol2 positions in all conditions")

# ── Difference plot ───────────────────────────────────────────────────────────
plot_apa_diff(apa_agg_ccre_wt_only, bins=GROUP_LEVELS,
              title="APA difference — WT-only cCRE loops  |  pol2 vs WT at WT positions")
plot_apa_diff(apa_agg_ccre_pol2_only, bins=GROUP_LEVELS,
              title="APA difference — pol2-only cCRE loops  |  pol2 vs WT at pol2 positions")


# ══════════════════════════════════════════════════════════════════════════════
# Missing piece from Analysis 3c: shared cCRE loops APA (WT positions,
# scored in both .hic files) — mirrors apa_agg_shared from Analysis 1c.
# Requires wt_only_ccre / pol2_only_ccre / shared_ccre already computed
# earlier in Analysis 3c.
# ══════════════════════════════════════════════════════════════════════════════
shared_ccre_wt <- np_grouped$wt %>% dplyr::slice(shared_ccre$w)
message("Shared cCRE loops (WT positions): ", nrow(shared_ccre_wt))

apa_raw_ccre_shared <- list(
  wt   = extract_apa_matrices_by(shared_ccre_wt, CONDITIONS$wt$hic,
                                 group_col="group", group_levels=GROUP_LEVELS),
  pol2 = extract_apa_matrices_by(shared_ccre_wt, CONDITIONS$pol2$hic,
                                 group_col="group", group_levels=GROUP_LEVELS)
)
apa_agg_ccre_shared <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_ccre_shared[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)

# ══════════════════════════════════════════════════════════════════════════════
# Wrapper: one 3x3 matched grid per cCRE loop type (P-P / E-P / E-E)
#   Same column structure as the expression-bin version:
#   Col1 = WT-only loops, Col2 = Pol2-only loops, Col3 = shared loops
#   Row1 = WT-scored, Row2 = Pol2-scored, Row3 = diff — one legend per row
# ══════════════════════════════════════════════════════════════════════════════
plot_apa_matched_grid_by_ccre_group <- function(groups      = GROUP_LEVELS,
                                                apa_agg_ccre_wt_only,
                                                apa_agg_ccre_pol2_only,
                                                apa_agg_ccre_shared,
                                                group_colors = GROUP_COLORS,
                                                out_prefix   = "apa_matched_grid_ccre",
                                                out_dir      = ".",
                                                width        = 11,
                                                height       = 9,
                                                save_files   = TRUE) {
  
  results <- setNames(vector("list", length(groups)), groups)
  
  for (g in groups) {
    message("Building matched grid \u2014 ", g, " loops...")
    
    composite <- plot_apa_matched_grid(
      bin                = g,
      apa_agg_wt_only    = apa_agg_ccre_wt_only,
      apa_agg_pol2_only  = apa_agg_ccre_pol2_only,
      apa_agg_shared     = apa_agg_ccre_shared,
      bin_color          = group_colors[[g]],
      title              = paste0("APA \u2014 ", g, " loops: matched WT vs Pol2 scoring")
    )
    
    results[[g]] <- composite
    
    if (save_files) {
      safe_name <- gsub("-", "", g)
      out_path  <- file.path(out_dir, paste0(out_prefix, "_", tolower(safe_name), ".pdf"))
      ggsave(out_path, composite, width = width, height = height, limitsize = FALSE)
      message("  Saved: ", out_path)
    }
  }
  
  invisible(results)
}

# ── Usage — produces 3 separate 3x3 figures, one per loop type ───────────────
ccre_grids <- plot_apa_matched_grid_by_ccre_group(
  groups                 = GROUP_LEVELS,
  apa_agg_ccre_wt_only   = apa_agg_ccre_wt_only,
  apa_agg_ccre_pol2_only = apa_agg_ccre_pol2_only,
  apa_agg_ccre_shared    = apa_agg_ccre_shared,
  out_prefix             = "apa_matched_grid_ccre",
  out_dir                = "."
)
# ccre_grids$`P-P`, ccre_grids$`E-P`, ccre_grids$`E-E` hold the composites

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 4: Singleton APA
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 4: Singleton P-Gene / pE-Gene / dE-Gene APA ────────────")
singleton_loops <- lapply(names(CONDITIONS), function(n) {
  message("Finding singletons: ", CONDITIONS[[n]]$label, "...")
  sing <- find_singleton_loops(np_classified[[n]], active_gene_gr)
  loop_coords <- loops_raw[[n]] %>%
    dplyr::mutate(loop_idx = row_number()) %>%
    dplyr::select(loop_idx, chr1, start1, end1, chr2, start2, end2,
                  observed, fdrDonut, fdrBL, fdrH, fdrV, centroid1, centroid2)
  result <- sing %>%
    dplyr::select(loop_idx, singleton_group, de_status, symbol) %>%
    distinct(loop_idx, singleton_group, .keep_all=TRUE) %>%
    left_join(loop_coords, by="loop_idx")
  message(CONDITIONS[[n]]$label, " singleton breakdown:")
  print(table(result$singleton_group))
  result
}); names(singleton_loops) <- names(CONDITIONS)

message("\n── Analysis 4 (unmatched) ───────────────────────────────────────────")
apa_raw_sing <- lapply(names(CONDITIONS), function(n) {
  message("Extracting ", CONDITIONS[[n]]$label, " singleton APA...")
  extract_apa_matrices_by(singleton_loops[[n]], CONDITIONS[[n]]$hic,
                          group_col="singleton_group", group_levels=SINGLETON_LEVELS)
}); names(apa_raw_sing) <- names(CONDITIONS)
apa_agg_sing <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_sing[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
plot_apa_grid_generic(apa_agg_sing, bin_levels=SINGLETON_LEVELS,
                      bin_colors=as.list(SINGLETON_COLORS),
                      title="APA — singleton  |  all active genes (unmatched)")
plot_apa_scores(apa_agg_sing, loops_annotated_cond=singleton_loops$wt,
                group_col="singleton_group", bins=SINGLETON_LEVELS,
                title="APA scores — singleton (unmatched)")

message("\n── Analysis 4 (matched): WT singleton positions × all conditions ────")
apa_raw_sing_m <- lapply(names(CONDITIONS), function(n) {
  message("Pulling WT singleton positions from ", CONDITIONS[[n]]$label, "...")
  extract_apa_matrices_by(singleton_loops$wt, CONDITIONS[[n]]$hic,
                          group_col="singleton_group", group_levels=SINGLETON_LEVELS)
}); names(apa_raw_sing_m) <- names(CONDITIONS)
apa_agg_sing_m <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_sing_m[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
plot_apa_grid_generic(apa_agg_sing_m, bin_levels=SINGLETON_LEVELS,
                      bin_colors=as.list(SINGLETON_COLORS),
                      title="APA — singleton  |  WT positions matched")
plot_apa_scores(apa_agg_sing_m, loops_annotated_cond=singleton_loops$wt,
                group_col="singleton_group", bins=SINGLETON_LEVELS,
                title="APA scores — singleton  |  WT positions matched")
plot_apa_diff(apa_agg_sing_m, bins=SINGLETON_LEVELS,
              title="APA difference — singleton  |  WT positions, treatment vs WT")

message("\n── Analysis 4b: Singleton by DE status ──────────────────────────────")
for (de_dir in c("Up", "Down")) {
  message("DE direction: ", de_dir)
  apa_raw_sing_de <- lapply(names(CONDITIONS), function(n) {
    extract_apa_matrices_by(
      singleton_loops[[n]] %>% filter(de_status == de_dir),
      CONDITIONS[[n]]$hic,
      group_col="singleton_group", group_levels=SINGLETON_LEVELS
    )
  }); names(apa_raw_sing_de) <- names(CONDITIONS)
  apa_agg_sing_de <- setNames(
    lapply(names(CONDITIONS), function(n)
      aggregate_apa(apa_raw_sing_de[[n]], DEPTH_SCALE[[n]])),
    names(CONDITIONS)
  )
  plot_apa_grid_generic(apa_agg_sing_de, bin_levels=SINGLETON_LEVELS,
                        bin_colors=as.list(SINGLETON_COLORS),
                        title=paste0("APA — singleton, ", de_dir, "regulated (unmatched)"))
  plot_apa_scores(apa_agg_sing_de,
                  loops_annotated_cond=singleton_loops$wt %>% filter(de_status == de_dir),
                  group_col="singleton_group", bins=SINGLETON_LEVELS,
                  title=paste0("APA scores — singleton, ", de_dir, "regulated (unmatched)"))
  wt_sing_de <- singleton_loops$wt %>% filter(de_status == de_dir)
  message("WT ", de_dir, " singleton loops: ", nrow(wt_sing_de))
  apa_raw_sing_de_m <- lapply(names(CONDITIONS), function(n) {
    extract_apa_matrices_by(wt_sing_de, CONDITIONS[[n]]$hic,
                            group_col="singleton_group", group_levels=SINGLETON_LEVELS)
  }); names(apa_raw_sing_de_m) <- names(CONDITIONS)
  apa_agg_sing_de_m <- setNames(
    lapply(names(CONDITIONS), function(n)
      aggregate_apa(apa_raw_sing_de_m[[n]], DEPTH_SCALE[[n]])),
    names(CONDITIONS)
  )
  plot_apa_grid_generic(apa_agg_sing_de_m, bin_levels=SINGLETON_LEVELS,
                        bin_colors=as.list(SINGLETON_COLORS),
                        title=paste0("APA — singleton, ", de_dir, "regulated  |  WT positions matched"))
  plot_apa_scores(apa_agg_sing_de_m, loops_annotated_cond=wt_sing_de,
                  group_col="singleton_group", bins=SINGLETON_LEVELS,
                  title=paste0("APA scores — singleton, ", de_dir, "regulated  |  WT positions matched"))
  plot_apa_diff(apa_agg_sing_de_m, bins=SINGLETON_LEVELS,
                title=paste0("APA difference — singleton, ", de_dir,
                             "regulated  |  WT positions, treatment vs WT"))
}



# ══════════════════════════════════════════════════════════════════════════════
# Analysis 5: WT loop positions × all conditions — expression bins (matched)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 5: WT positions × all conditions (expression bins) ──────")
apa_raw_wt_matched <- lapply(names(CONDITIONS), function(n) {
  message("Pulling WT loop positions from ", CONDITIONS[[n]]$label, " .hic...")
  extract_apa_matrices(loops_annotated$wt, CONDITIONS[[n]]$hic)
}); names(apa_raw_wt_matched) <- names(CONDITIONS)
apa_agg_wt_matched <- setNames(
  lapply(names(CONDITIONS), function(n)
    aggregate_apa(apa_raw_wt_matched[[n]], DEPTH_SCALE[[n]])),
  names(CONDITIONS)
)
plot_apa_grid_generic(apa_agg_wt_matched, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      title="APA — WT loop positions in all conditions")
plot_apa_scores(apa_agg_wt_matched, loops_annotated_cond=loops_annotated$wt,
                group_col="expr_bin", de_filter=NULL,
                title="APA scores — WT loop positions in all conditions")
plot_apa_diff(apa_agg_wt_matched, bins=BINS,
              title="APA difference — WT positions, treatment vs WT")

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 6: 2D-stratified APA — WT vs pol2 only
# Loop length × expression  |  Loop length × P-P/E-P/E-E  |  Expr × P-P/E-P/E-E
# Note: these analyses compare WT and Pol2 only, regardless of DS_GROUP
# ══════════════════════════════════════════════════════════════════════════════
LENGTH_LEVELS <- c("Short","Medium","Long")

add_length_bins <- function(loops_df, labels = LENGTH_LEVELS) {
  df <- loops_df %>%
    dplyr::mutate(
      loop_length = abs(
        ifelse(!is.na(centroid2), centroid2, (start2+end2)/2) -
          ifelse(!is.na(centroid1), centroid1, (start1+end1)/2)
      )
    )
  q <- quantile(df$loop_length, c(1/3, 2/3), na.rm=TRUE)
  breaks <- c(0, q[[1]], q[[2]], Inf)
  message(sprintf("Loop length breaks: Short <%.0fkb  |  Medium %.0f–%.0fkb  |  Long >%.0fkb",
                  q[[1]]/1000, q[[1]]/1000, q[[2]]/1000, q[[2]]/1000))
  df %>%
    dplyr::mutate(
      length_bin = factor(
        cut(loop_length, breaks=breaks, labels=labels, include.lowest=TRUE),
        levels=labels
      )
    )
}

# ── Helper: 2D matched APA comparing WT and Pol2 ──────────────────────────────
run_matched_2d_apa <- function(loops_df, outer_col, outer_levels,
                               inner_col, inner_levels, min_loops = 10) {
  lapply(outer_levels, function(outer_val) {
    loops_sub <- loops_df %>% filter(.data[[outer_col]] == outer_val)
    message(sprintf("  %s = %s: %d loops", outer_col, outer_val, nrow(loops_sub)))
    apa_wt_raw   <- extract_apa_matrices_by(loops_sub, CONDITIONS$wt$hic,
                                            group_col=inner_col,
                                            group_levels=inner_levels,
                                            min_loops=min_loops)
    apa_pol2_raw <- extract_apa_matrices_by(loops_sub, CONDITIONS$pol2$hic,
                                            group_col=inner_col,
                                            group_levels=inner_levels,
                                            min_loops=min_loops)
    list(
      wt   = aggregate_apa(apa_wt_raw,   DEPTH_SCALE$wt),
      pol2 = aggregate_apa(apa_pol2_raw, DEPTH_SCALE$pol2)
    )
  }) %>% setNames(outer_levels)
}

summarise_2d_scores <- function(apa_2d, outer_levels, inner_levels) {
  bind_rows(lapply(outer_levels, function(outer_val) {
    bind_rows(lapply(inner_levels, function(inner_val) {
      data.frame(
        outer    = outer_val,
        inner    = inner_val,
        apa_wt   = apa_score(apa_2d[[outer_val]]$wt[[inner_val]]),
        apa_pol2 = apa_score(apa_2d[[outer_val]]$pol2[[inner_val]])
      )
    }))
  })) %>%
    dplyr::mutate(
      log2fc     = log2(apa_pol2 / apa_wt),
      pct_change = ((apa_pol2 - apa_wt) / apa_wt) * 100
    )
}

message("\n── Analysis 6: preparing annotated WT loop sets ─────────────────────")
loops_wt_len <- add_length_bins(loops_annotated$wt)
message("Length bin distribution:"); print(table(loops_wt_len$length_bin))

np_wt_with_expr <- np_grouped$wt %>%
  left_join(
    loops_annotated$wt %>%
      dplyr::select(chr1, start1, end1, chr2, start2, end2,
                    expr_bin, de_status, TPM, symbol),
    by=c("chr1","start1","end1","chr2","start2","end2")
  ) %>%
  dplyr::mutate(
    expr_bin = replace_na(as.character(expr_bin), "Unannotated"),
    expr_bin = factor(expr_bin, levels=c(BINS, "Unannotated"))
  )

np_wt_with_len <- add_length_bins(np_grouped$wt) %>%
  left_join(
    loops_annotated$wt %>%
      dplyr::select(chr1, start1, end1, chr2, start2, end2,
                    expr_bin, de_status, TPM, symbol),
    by=c("chr1","start1","end1","chr2","start2","end2")
  ) %>%
  dplyr::mutate(
    expr_bin = replace_na(as.character(expr_bin), "Unannotated"),
    expr_bin = factor(expr_bin, levels=c(BINS, "Unannotated"))
  )

BUFFER <- 5  # tighter buffer for length-stratified analyses

message("\n── Analysis 6a: loop length × expression group ──────────────────────")
apa_6a <- run_matched_2d_apa(
  loops_df=loops_wt_len, outer_col="length_bin", outer_levels=LENGTH_LEVELS,
  inner_col="expr_bin",  inner_levels=BINS
)
scores_6a <- summarise_2d_scores(apa_6a, LENGTH_LEVELS, BINS)
message("APA scores — loop length × expression:"); print(scores_6a)
for (len in LENGTH_LEVELS) {
  loops_sub <- loops_wt_len %>% filter(length_bin == len)
  
  plot_apa_grid_generic(apa_agg=apa_6a[[len]], bin_levels=BINS,
                        bin_colors=as.list(BIN_COLORS), conditions=c("wt","pol2"),
                        title=paste0("APA — ", len, " loops  |  WT vs pol2  |  by expression"))
  
  plot_apa_scores(apa_agg=apa_6a[[len]], loops_annotated_cond=loops_sub,
                  group_col="expr_bin", bins=BINS, conditions=c("wt","pol2"),
                  title=paste0("APA scores — ", len, " loops  |  by expression"))
 
  plot_apa_diff(apa_agg=apa_6a[[len]], bins=BINS, conditions="pol2",
                title=paste0("APA diff — ", len, " loops  |  pol2 vs WT  |  by expression"))
}

message("\n── Analysis 6b: loop length × P-P/E-P/E-E ──────────────────────────")
apa_6b <- run_matched_2d_apa(
  loops_df=np_wt_with_len, outer_col="length_bin", outer_levels=LENGTH_LEVELS,
  inner_col="group",       inner_levels=GROUP_LEVELS
)
scores_6b <- summarise_2d_scores(apa_6b, LENGTH_LEVELS, GROUP_LEVELS)
message("APA scores — loop length × loop type:"); print(scores_6b)
LEN_COLORS <- setNames(c("#74add1","#abdda4","#f46d43"), LENGTH_LEVELS)
for (len in LENGTH_LEVELS) {
  loops_sub <- np_wt_with_len %>% filter(length_bin == len)
  #grid::grid.newpage()
  plot_apa_grid_generic(apa_agg=apa_6b[[len]], bin_levels=GROUP_LEVELS,
                        bin_colors=as.list(GROUP_COLORS), conditions=c("wt","pol2"),
                        title=paste0("APA — ", len, " loops  |  WT vs pol2  |  by loop type"))
  #grid::grid.newpage()
  plot_apa_scores(apa_agg=apa_6b[[len]], loops_annotated_cond=loops_sub,
                  group_col="group", bins=GROUP_LEVELS, conditions=c("wt","pol2"),
                  title=paste0("APA scores — ", len, " loops  |  by loop type"))
  #grid::grid.newpage()
  plot_apa_diff(apa_agg=apa_6b[[len]], bins=GROUP_LEVELS, conditions="pol2",
                title=paste0("APA diff — ", len, " loops  |  pol2 vs WT  |  by loop type"))
}

# PP / EP / EE by length individually
for (loop_type in GROUP_LEVELS) {
  safe_type <- gsub("-","", loop_type)
  apa_len_type <- list(
    wt   = setNames(lapply(LENGTH_LEVELS, function(len) apa_6b[[len]]$wt[[loop_type]]),   LENGTH_LEVELS),
    pol2 = setNames(lapply(LENGTH_LEVELS, function(len) apa_6b[[len]]$pol2[[loop_type]]), LENGTH_LEVELS)
  )
  #grid::grid.newpage()
  plot_apa_grid_generic(apa_agg=apa_len_type, bin_levels=LENGTH_LEVELS,
                        bin_colors=as.list(LEN_COLORS), conditions=c("wt","pol2"),
                        title=paste0("APA — ", loop_type, " loops  |  WT vs pol2  |  by loop length"))
  #grid::grid.newpage()
  plot_apa_scores(apa_agg=apa_len_type,
                  loops_annotated_cond=np_wt_with_len %>% filter(group == loop_type),
                  group_col="length_bin", bins=LENGTH_LEVELS, conditions=c("wt","pol2"),
                  title=paste0("APA scores — ", loop_type, " loops  |  by loop length"))
  #grid::grid.newpage()
  plot_apa_diff(apa_agg=apa_len_type, bins=LENGTH_LEVELS, conditions="pol2",
                title=paste0("APA diff — ", loop_type, " loops  |  pol2 vs WT  |  by loop length"))
}

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 6c: P-P/E-P/E-E by loop length — exclusive and shared
#   Mirrors Analysis 6b but for WT-only, pol2-only, and shared cCRE loops
#   Requires wt_only_ccre, pol2_only_ccre, shared_ccre from Analysis 3c
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 6c: cCRE loops by length — exclusive and shared ──────────")

# Shared loops — WT coordinates, both conditions
shared_ccre_wt <- np_grouped$wt %>% dplyr::slice(shared_ccre$w)

all_dists <- abs(as.integer((np_grouped$wt$start2   + np_grouped$wt$end2)  /2) -
                   as.integer((np_grouped$wt$start1   + np_grouped$wt$end1)  /2))

len_low  <- quantile(all_dists, 1/3, na.rm=TRUE)
len_high <- quantile(all_dists, 2/3, na.rm=TRUE)
LENGTH_LEVELS <- c("Short","Medium","Long")

add_length_bin <- function(df) {
  df %>%
    dplyr::mutate(
      loop_dist  = abs(as.integer((start2+end2)/2) - as.integer((start1+end1)/2)),
      length_bin = dplyr::case_when(
        loop_dist <= len_low  ~ "Short",
        loop_dist <= len_high ~ "Medium",
        TRUE                  ~ "Long"
      ),
      length_bin = factor(length_bin, levels=LENGTH_LEVELS)
    )
}

summary(all_dists / 1000)
quantile(all_dists / 1000, probs=c(0.1,0.25,0.33,0.5,0.67,0.75,0.9), na.rm=TRUE)

wt_only_len   <- add_length_bin(wt_only_ccre)
pol2_only_len <- add_length_bin(pol2_only_ccre)
shared_len    <- add_length_bin(shared_ccre_wt)

# Helper: build apa_6b-style structure for a loop set scored at two .hic files
build_len_ccre_apa <- function(loops_with_len, hic_a, hic_b,
                               label_a="wt", label_b="pol2",
                               min_loops=20L) {
  setNames(lapply(LENGTH_LEVELS, function(len) {
    loops_len <- loops_with_len %>% filter(length_bin == len)
    extract_bin <- function(hic) {
      if (nrow(loops_len) < min_loops) {
        message("  Skipping ", len, " (n=", nrow(loops_len), " < ", min_loops, ")")
        return(setNames(lapply(GROUP_LEVELS, function(g) NULL), GROUP_LEVELS))
      }
      # Check per group too
      result <- extract_apa_matrices_by(loops_len, hic,
                                        group_col="group", group_levels=GROUP_LEVELS)
      # NULL out any group with too few loops
      lapply(GROUP_LEVELS, function(g) {
        n_g <- sum(loops_len$group == g, na.rm=TRUE)
        if (n_g < min_loops) {
          message("  Skipping ", len, " — ", g, " (n=", n_g, " < ", min_loops, ")")
          NULL
        } else {
          result[[g]]
        }
      }) %>% setNames(GROUP_LEVELS)
    }
    list(wt=extract_bin(hic_a), pol2=extract_bin(hic_b)) %>%
      setNames(c(label_a, label_b))
  }), LENGTH_LEVELS)
}

message("\nExtracting APA — WT-only cCRE loops by length...")
apa_6c_wt_only <- build_len_ccre_apa(wt_only_len,
                                     CONDITIONS$wt$hic, CONDITIONS$pol2$hic)

message("\nExtracting APA — pol2-only cCRE loops by length...")
apa_6c_pol2_only <- build_len_ccre_apa(pol2_only_len,
                                       CONDITIONS$pol2$hic, CONDITIONS$wt$hic,
                                       label_a="pol2", label_b="wt")

message("\nExtracting APA — shared cCRE loops by length (WT positions)...")
apa_6c_shared <- build_len_ccre_apa(shared_len,
                                    CONDITIONS$wt$hic, CONDITIONS$pol2$hic)

plot_len_ccre <- function(apa_6c, loops_with_len, label_a, label_b, set_label,
                          min_loops=20L) {
  for (loop_type in GROUP_LEVELS) {
    matrices_a <- lapply(LENGTH_LEVELS, function(len)
      apa_6c[[len]][[label_a]][[loop_type]])
    matrices_b <- lapply(LENGTH_LEVELS, function(len)
      apa_6c[[len]][[label_b]][[loop_type]])
    
    if (any(sapply(c(matrices_a, matrices_b), is.null))) {
      message("Skipping ", loop_type, " in ", set_label,
              " — one or more length bins below minimum n=", min_loops)
      next
    }
    
    apa_len_type <- list(
      wt   = setNames(matrices_a, LENGTH_LEVELS),
      pol2 = setNames(matrices_b, LENGTH_LEVELS)
    )
    
    plot_apa_grid_generic(
      apa_agg=apa_len_type, bin_levels=LENGTH_LEVELS,
      bin_colors=as.list(LEN_COLORS), conditions=c("wt","pol2"),
      title=paste0("APA — ", loop_type, "  |  ", set_label, "  |  by loop length")
    )
    plot_apa_scores(
      apa_agg=apa_len_type,
      loops_annotated_cond=list(
        wt   = loops_with_len %>% filter(group == loop_type),
        pol2 = loops_with_len %>% filter(group == loop_type)
      ),
      group_col="length_bin", bins=LENGTH_LEVELS, conditions=c("wt","pol2"),
      title=paste0("APA scores — ", loop_type, "  |  ", set_label, "  |  by length")
    )
    plot_apa_diff(
      apa_agg=apa_len_type, bins=LENGTH_LEVELS, conditions="pol2",
      title=paste0("APA diff — ", loop_type, "  |  ", set_label, "  |  pol2 vs WT")
    )
  }
}

message("Length level names: ", paste(names(apa_6c_wt_only), collapse=", "))
message("Condition names:    ", paste(names(apa_6c_wt_only$Short), collapse=", "))
message("Group names:        ", paste(names(apa_6c_wt_only$Short$wt), collapse=", "))
message("GROUP_LEVELS:       ", paste(GROUP_LEVELS, collapse=", "))

message("\n── Plotting WT-only cCRE loops by length ─────────────────────────────")
plot_len_ccre(apa_6c_wt_only, wt_only_len,
              label_a="wt", label_b="pol2",
              set_label="WT-only loops  |  WT positions in both conditions")

message("\n── Plotting pol2-only cCRE loops by length ───────────────────────────")
plot_len_ccre(apa_6c_pol2_only, pol2_only_len,
              label_a="pol2", label_b="wt",
              set_label="pol2-only loops  |  pol2 positions in both conditions")

message("\n── Plotting shared cCRE loops by length ──────────────────────────────")
plot_len_ccre(apa_6c_shared, shared_len,
              label_a="wt", label_b="pol2",
              set_label="shared loops  |  WT positions in both conditions")

message("\n── Analysis 6c: expression group × P-P/E-P/E-E ─────────────────────")
apa_6c <- run_matched_2d_apa(
  loops_df=np_wt_with_expr %>% filter(expr_bin %in% BINS),
  outer_col="expr_bin",  outer_levels=BINS,
  inner_col="group",     inner_levels=GROUP_LEVELS
)
scores_6c <- summarise_2d_scores(apa_6c, BINS, GROUP_LEVELS)
message("APA scores — expression × loop type:"); print(scores_6c)
for (expr in BINS) {
  loops_sub <- np_wt_with_expr %>% filter(expr_bin == expr)
  grid::grid.newpage()
  plot_apa_grid_generic(apa_agg=apa_6c[[expr]], bin_levels=GROUP_LEVELS,
                        bin_colors=as.list(GROUP_COLORS), conditions=c("wt","pol2"),
                        title=paste0("APA — ", expr, " expression  |  WT vs pol2  |  by loop type"))
  grid::grid.newpage()
  plot_apa_scores(apa_agg=apa_6c[[expr]], loops_annotated_cond=loops_sub,
                  group_col="group", bins=GROUP_LEVELS, conditions=c("wt","pol2"),
                  title=paste0("APA scores — ", expr, " expression  |  by loop type"))
  grid::grid.newpage()
  plot_apa_diff(apa_agg=apa_6c[[expr]], bins=GROUP_LEVELS, conditions="pol2",
                title=paste0("APA diff — ", expr, " expression  |  pol2 vs WT  |  by loop type"))
}

message("\n── Analysis 6d: loop length × P-P/E-P/E-E — High expression only ──")
np_wt_high_len <- np_wt_with_len %>% filter(expr_bin == "High")
apa_6d <- run_matched_2d_apa(
  loops_df=np_wt_high_len, outer_col="length_bin", outer_levels=LENGTH_LEVELS,
  inner_col="group",        inner_levels=GROUP_LEVELS
)
scores_6d <- summarise_2d_scores(apa_6d, LENGTH_LEVELS, GROUP_LEVELS)
message("APA scores — length × loop type (High expression only):"); print(scores_6d)
for (len in LENGTH_LEVELS) {
  loops_sub <- np_wt_high_len %>% filter(length_bin == len)
  #grid::grid.newpage()
  plot_apa_grid_generic(apa_agg=apa_6d[[len]], bin_levels=GROUP_LEVELS,
                        bin_colors=as.list(GROUP_COLORS), conditions=c("wt","pol2"),
                        title=paste0("APA — ", len, " loops  |  High expr  |  WT vs pol2  |  by loop type"))
  #grid::grid.newpage()
  plot_apa_scores(apa_agg=apa_6d[[len]], loops_annotated_cond=loops_sub,
                  group_col="group", bins=GROUP_LEVELS, conditions=c("wt","pol2"),
                  title=paste0("APA scores — ", len, " loops  |  High expr  |  by loop type"))
  #grid::grid.newpage()
  plot_apa_diff(apa_agg=apa_6d[[len]], bins=GROUP_LEVELS, conditions="pol2",
                title=paste0("APA diff — ", len, " loops  |  High expr  |  pol2 vs WT  |  by loop type"))
}

message("\n── Analysis 6 summary ───────────────────────────────────────────────")
scores_6_all <- bind_rows(
  scores_6a %>% dplyr::mutate(analysis="Length × Expression"),
  scores_6b %>% dplyr::mutate(analysis="Length × Loop type"),
  scores_6c %>% dplyr::mutate(analysis="Expression × Loop type"),
  scores_6d %>% dplyr::mutate(analysis="Length × Loop type (High expr)")
) %>% dplyr::arrange(analysis, log2fc)
print(scores_6_all)

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 7: pE-pE vs pE-dE vs dE-dE loops
# ══════════════════════════════════════════════════════════════════════════════
BUFFER <- 7
message("\n── Analysis 7: pE-pE vs dE-dE loops (matched WT positions) ─────────")
EE_SUBTYPE_LEVELS <- c("pE-pE","pE-dE","dE-dE")
EE_SUBTYPE_COLORS <- c("pE-pE"="#74add1","pE-dE"="#abdda4","dE-dE"="#f46d43")
ee_subtyped <- np_classified$wt %>%
  filter(
    anchor1_ccre %in% c("Proximal-Enhancer","Distal-Enhancer"),
    anchor2_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")
  ) %>%
  dplyr::mutate(
    ee_subtype = dplyr::case_when(
      anchor1_ccre=="Proximal-Enhancer" & anchor2_ccre=="Proximal-Enhancer" ~ "pE-pE",
      anchor1_ccre=="Distal-Enhancer"   & anchor2_ccre=="Distal-Enhancer"   ~ "dE-dE",
      TRUE ~ "pE-dE"
    ),
    ee_subtype = factor(ee_subtype, levels=EE_SUBTYPE_LEVELS),
    group      = ee_subtype
  )
message("E-E subtype breakdown:"); print(table(ee_subtyped$ee_subtype))
apa_ee_wt_raw   <- extract_apa_matrices_by(ee_subtyped, CONDITIONS$wt$hic,
                                           group_col="group", group_levels=EE_SUBTYPE_LEVELS)
apa_ee_pol2_raw <- extract_apa_matrices_by(ee_subtyped, CONDITIONS$pol2$hic,
                                           group_col="group", group_levels=EE_SUBTYPE_LEVELS)
apa_agg_ee <- list(
  wt   = aggregate_apa(apa_ee_wt_raw,   DEPTH_SCALE$wt),
  pol2 = aggregate_apa(apa_ee_pol2_raw, DEPTH_SCALE$pol2)
)
scores_ee <- bind_rows(lapply(EE_SUBTYPE_LEVELS, function(g) {
  data.frame(subtype=g,
             apa_wt  =apa_score(apa_agg_ee$wt[[g]]),
             apa_pol2=apa_score(apa_agg_ee$pol2[[g]]))
})) %>%
  dplyr::mutate(log2fc=log2(apa_pol2/apa_wt),
                pct_change=((apa_pol2-apa_wt)/apa_wt)*100)
message("E-E subtype APA scores:"); print(scores_ee)
plot_apa_grid_generic(apa_agg=apa_agg_ee, bin_levels=EE_SUBTYPE_LEVELS,
                      bin_colors=as.list(EE_SUBTYPE_COLORS), conditions=c("wt","pol2"),
                      title="APA — pE-pE vs pE-dE vs dE-dE  |  WT positions in WT & pol2")
plot_apa_scores(apa_agg=apa_agg_ee, loops_annotated_cond=ee_subtyped,
                group_col="ee_subtype", bins=EE_SUBTYPE_LEVELS, conditions=c("wt","pol2"),
                title="APA scores — pE-pE vs pE-dE vs dE-dE  |  WT positions matched")
plot_apa_diff(apa_agg=apa_agg_ee, bins=EE_SUBTYPE_LEVELS, conditions="pol2",
              title="APA difference — pE-pE vs pE-dE vs dE-dE  |  pol2 vs WT")

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 8: P-P and E-P loops stratified by expression
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 8: P-P loops by expression group ────────────────────────")
pp_loops <- np_wt_with_expr %>%
  filter(group == "P-P", expr_bin %in% BINS) %>%
  dplyr::mutate(expr_bin = factor(expr_bin, levels=BINS))
message("P-P loops by expression bin:"); print(table(pp_loops$expr_bin))
apa_pp <- list(
  wt   = aggregate_apa(extract_apa_matrices_by(pp_loops, CONDITIONS$wt$hic,
                                               group_col="expr_bin", group_levels=BINS),
                       DEPTH_SCALE$wt),
  pol2 = aggregate_apa(extract_apa_matrices_by(pp_loops, CONDITIONS$pol2$hic,
                                               group_col="expr_bin", group_levels=BINS),
                       DEPTH_SCALE$pol2)
)
plot_apa_grid_generic(apa_agg=apa_pp, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      conditions=c("wt","pol2"),
                      title="APA — P-P loops  |  WT vs pol2  |  by expression")
plot_apa_scores(apa_agg=apa_pp, loops_annotated_cond=pp_loops,
                group_col="expr_bin", bins=BINS, conditions=c("wt","pol2"),
                title="APA scores — P-P loops  |  by expression")
plot_apa_diff(apa_agg=apa_pp, bins=BINS, conditions="pol2",
              title="APA difference — P-P loops  |  pol2 vs WT  |  by expression")

message("\n── Analysis 8b: E-P loops by expression group ───────────────────────")
ep_loops <- np_wt_with_expr %>%
  filter(group == "E-P", expr_bin %in% BINS) %>%
  dplyr::mutate(expr_bin = factor(expr_bin, levels=BINS))
message("E-P loops by expression bin:"); print(table(ep_loops$expr_bin))
apa_ep <- list(
  wt   = aggregate_apa(extract_apa_matrices_by(ep_loops, CONDITIONS$wt$hic,
                                               group_col="expr_bin", group_levels=BINS),
                       DEPTH_SCALE$wt),
  pol2 = aggregate_apa(extract_apa_matrices_by(ep_loops, CONDITIONS$pol2$hic,
                                               group_col="expr_bin", group_levels=BINS),
                       DEPTH_SCALE$pol2)
)
plot_apa_grid_generic(apa_agg=apa_ep, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      conditions=c("wt","pol2"),
                      title="APA — E-P loops  |  WT vs pol2  |  by expression")
plot_apa_scores(apa_agg=apa_ep, loops_annotated_cond=ep_loops,
                group_col="expr_bin", bins=BINS, conditions=c("wt","pol2"),
                title="APA scores — E-P loops  |  by expression")
plot_apa_diff(apa_agg=apa_ep, bins=BINS, conditions="pol2",
              title="APA difference — E-P loops  |  pol2 vs WT  |  by expression")

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 9: P-P, E-P, E-E and CTCF loop stratifications
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 9a: P-P loops by expression ─────────────────────────────")
apa_9_pp <- list(
  wt   = aggregate_apa(extract_apa_matrices_by(pp_loops, CONDITIONS$wt$hic,
                                               group_col="expr_bin", group_levels=BINS),
                       DEPTH_SCALE$wt),
  pol2 = aggregate_apa(extract_apa_matrices_by(pp_loops, CONDITIONS$pol2$hic,
                                               group_col="expr_bin", group_levels=BINS),
                       DEPTH_SCALE$pol2)
)
plot_apa_grid_generic(apa_agg=apa_9_pp, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      conditions=c("wt","pol2"),
                      title="APA — P-P loops  |  WT vs pol2  |  by expression")
plot_apa_scores(apa_agg=apa_9_pp, loops_annotated_cond=pp_loops,
                group_col="expr_bin", bins=BINS, conditions=c("wt","pol2"),
                title="APA scores — P-P loops  |  by expression")
plot_apa_diff(apa_agg=apa_9_pp, bins=BINS, conditions="pol2",
              title="APA diff — P-P loops  |  pol2 vs WT  |  by expression")

message("\n── Analysis 9b: E-P loops by expression ─────────────────────────────")
apa_9_ep <- list(
  wt   = aggregate_apa(extract_apa_matrices_by(ep_loops, CONDITIONS$wt$hic,
                                               group_col="expr_bin", group_levels=BINS),
                       DEPTH_SCALE$wt),
  pol2 = aggregate_apa(extract_apa_matrices_by(ep_loops, CONDITIONS$pol2$hic,
                                               group_col="expr_bin", group_levels=BINS),
                       DEPTH_SCALE$pol2)
)
plot_apa_grid_generic(apa_agg=apa_9_ep, bin_levels=BINS, bin_colors=as.list(BIN_COLORS),
                      conditions=c("wt","pol2"),
                      title="APA — E-P loops  |  WT vs pol2  |  by expression")
plot_apa_scores(apa_agg=apa_9_ep, loops_annotated_cond=ep_loops,
                group_col="expr_bin", bins=BINS, conditions=c("wt","pol2"),
                title="APA scores — E-P loops  |  by expression")
plot_apa_diff(apa_agg=apa_9_ep, bins=BINS, conditions="pol2",
              title="APA diff — E-P loops  |  pol2 vs WT  |  by expression")

message("\n── Analysis 9c: E-P loops by length ─────────────────────────────────")
ep_loops_len <- np_wt_with_len %>% filter(group == "E-P")
apa_9_ep_len <- list(
  wt   = aggregate_apa(extract_apa_matrices_by(ep_loops_len, CONDITIONS$wt$hic,
                                               group_col="length_bin", group_levels=LENGTH_LEVELS),
                       DEPTH_SCALE$wt),
  pol2 = aggregate_apa(extract_apa_matrices_by(ep_loops_len, CONDITIONS$pol2$hic,
                                               group_col="length_bin", group_levels=LENGTH_LEVELS),
                       DEPTH_SCALE$pol2)
)
plot_apa_grid_generic(apa_agg=apa_9_ep_len, bin_levels=LENGTH_LEVELS,
                      bin_colors=as.list(LEN_COLORS), conditions=c("wt","pol2"),
                      title="APA — E-P loops  |  WT vs pol2  |  by loop length")
plot_apa_scores(apa_agg=apa_9_ep_len, loops_annotated_cond=ep_loops_len,
                group_col="length_bin", bins=LENGTH_LEVELS, conditions=c("wt","pol2"),
                title="APA scores — E-P loops  |  by loop length")
plot_apa_diff(apa_agg=apa_9_ep_len, bins=LENGTH_LEVELS, conditions="pol2",
              title="APA diff — E-P loops  |  pol2 vs WT  |  by loop length")

message("\n── Analysis 9d: E-E loops by length ─────────────────────────────────")
ee_loops_len <- np_wt_with_len %>% filter(group == "E-E")
apa_9_ee_len <- list(
  wt   = aggregate_apa(extract_apa_matrices_by(ee_loops_len, CONDITIONS$wt$hic,
                                               group_col="length_bin", group_levels=LENGTH_LEVELS),
                       DEPTH_SCALE$wt),
  pol2 = aggregate_apa(extract_apa_matrices_by(ee_loops_len, CONDITIONS$pol2$hic,
                                               group_col="length_bin", group_levels=LENGTH_LEVELS),
                       DEPTH_SCALE$pol2)
)
plot_apa_grid_generic(apa_agg=apa_9_ee_len, bin_levels=LENGTH_LEVELS,
                      bin_colors=as.list(LEN_COLORS), conditions=c("wt","pol2"),
                      title="APA — E-E loops  |  WT vs pol2  |  by loop length")
plot_apa_scores(apa_agg=apa_9_ee_len, loops_annotated_cond=ee_loops_len,
                group_col="length_bin", bins=LENGTH_LEVELS, conditions=c("wt","pol2"),
                title="APA scores — E-E loops  |  by loop length")
plot_apa_diff(apa_agg=apa_9_ee_len, bins=LENGTH_LEVELS, conditions="pol2",
              title="APA diff — E-E loops  |  pol2 vs WT  |  by loop length")

message("\n── Analysis 9e: CTCF-CTCF, CTCF-P, CTCF-E loops ────────────────────")
CTCF_LEVELS <- c("CTCF-CTCF","CTCF-P","CTCF-E")
CTCF_COLORS <- c("CTCF-CTCF"="#5e4fa2","CTCF-P"="#3288bd","CTCF-E"="#66c2a5")
ctcf_loops <- np_classified$wt %>%
  filter(anchor1_ctcf | anchor2_ctcf) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      anchor1_ctcf & anchor2_ctcf ~ "CTCF-CTCF",
      (anchor1_ctcf & anchor1_ccre == "Promoter") |
        (anchor2_ctcf & anchor2_ccre == "Promoter") ~ "CTCF-P",
      (anchor1_ctcf & anchor1_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")) |
        (anchor2_ctcf & anchor2_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")) ~ "CTCF-E",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  dplyr::mutate(group = factor(group, levels=CTCF_LEVELS))
message("CTCF loop type breakdown:"); print(table(ctcf_loops$group))
apa_ctcf <- list(
  wt   = aggregate_apa(extract_apa_matrices_by(ctcf_loops, CONDITIONS$wt$hic,
                                               group_col="group", group_levels=CTCF_LEVELS),
                       DEPTH_SCALE$wt),
  pol2 = aggregate_apa(extract_apa_matrices_by(ctcf_loops, CONDITIONS$pol2$hic,
                                               group_col="group", group_levels=CTCF_LEVELS),
                       DEPTH_SCALE$pol2)
)
plot_apa_grid_generic(apa_agg=apa_ctcf, bin_levels=CTCF_LEVELS,
                      bin_colors=as.list(CTCF_COLORS), conditions=c("wt","pol2"),
                      title="APA — CTCF loop types  |  WT vs pol2")
plot_apa_scores(apa_agg=apa_ctcf, loops_annotated_cond=ctcf_loops,
                group_col="group", bins=CTCF_LEVELS, conditions=c("wt","pol2"),
                title="APA scores — CTCF loop types")
plot_apa_diff(apa_agg=apa_ctcf, bins=CTCF_LEVELS, conditions="pol2",
              title="APA diff — CTCF loop types  |  pol2 vs WT")

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 10: Unmatched comparisons — pol2 vs WT with WT-derived length breaks
# ══════════════════════════════════════════════════════════════════════════════
UNMATCHED_CONDITIONS <- c("wt","pol2")

wt_lengths <- abs(
  ifelse(!is.na(np_classified$wt$centroid2), np_classified$wt$centroid2,
         (np_classified$wt$start2 + np_classified$wt$end2)/2) -
    ifelse(!is.na(np_classified$wt$centroid1), np_classified$wt$centroid1,
           (np_classified$wt$start1 + np_classified$wt$end1)/2)
)
q_wt <- quantile(wt_lengths, c(1/3, 2/3), na.rm=TRUE)
WT_LENGTH_BREAKS <- c(0, q_wt[[1]], q_wt[[2]], Inf)
message(sprintf("WT length breaks: Short <%.0fkb  |  Medium %.0f–%.0fkb  |  Long >%.0fkb",
                q_wt[[1]]/1000, q_wt[[1]]/1000, q_wt[[2]]/1000, q_wt[[2]]/1000))

add_length_bins <- function(loops_df, breaks = WT_LENGTH_BREAKS, labels = LENGTH_LEVELS) {
  loops_df %>%
    dplyr::mutate(
      loop_length = abs(
        ifelse(!is.na(centroid2), centroid2, (start2+end2)/2) -
          ifelse(!is.na(centroid1), centroid1, (start1+end1)/2)
      ),
      length_bin = factor(
        cut(loop_length, breaks=breaks, labels=labels, include.lowest=TRUE),
        levels=labels
      )
    )
}

np_len <- lapply(names(np_classified), function(n) {
  add_length_bins(np_classified[[n]]) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        loop_ccre_type %in% names(GROUP_MAP) ~ GROUP_MAP[loop_ccre_type],
        TRUE ~ NA_character_
      ),
      group = factor(group, levels=GROUP_LEVELS)
    )
}); names(np_len) <- names(np_classified)

loops_wt_len   <- add_length_bins(loops_annotated$wt)
np_wt_with_len <- add_length_bins(np_grouped$wt) %>%
  left_join(
    loops_annotated$wt %>%
      dplyr::select(chr1, start1, end1, chr2, start2, end2,
                    expr_bin, de_status, TPM, symbol),
    by=c("chr1","start1","end1","chr2","start2","end2")
  ) %>%
  dplyr::mutate(
    expr_bin = replace_na(as.character(expr_bin), "Unannotated"),
    expr_bin = factor(expr_bin, levels=c(BINS, "Unannotated"))
  )

run_unmatched_by_length <- function(loop_type) {
  apa_raw <- lapply(UNMATCHED_CONDITIONS, function(n) {
    loops_sub <- np_len[[n]] %>%
      filter(group == loop_type) %>%
      dplyr::mutate(group = length_bin)
    message(sprintf("  %s %s: %d loops", CONDITIONS[[n]]$label, loop_type, nrow(loops_sub)))
    extract_apa_matrices_by(loops_sub, CONDITIONS[[n]]$hic,
                            group_col="group", group_levels=LENGTH_LEVELS)
  }); names(apa_raw) <- UNMATCHED_CONDITIONS
  setNames(lapply(UNMATCHED_CONDITIONS, function(n)
    aggregate_apa(apa_raw[[n]], DEPTH_SCALE[[n]])), UNMATCHED_CONDITIONS)
}

for (loop_type in GROUP_LEVELS) {
  label <- paste0("Analysis 10 — ", loop_type, " loops by length (unmatched)")
  message("\n── ", label, " ──")
  apa_unmatched <- run_unmatched_by_length(loop_type)
  plot_apa_grid_generic(apa_agg=apa_unmatched, bin_levels=LENGTH_LEVELS,
                        bin_colors=as.list(LEN_COLORS),
                        conditions=UNMATCHED_CONDITIONS,
                        title=paste0("APA — ", loop_type, " loops  |  unmatched  |  by loop length"))
  plot_apa_scores(apa_agg=apa_unmatched,
                  loops_annotated_cond=np_len$wt %>% filter(group == loop_type),
                  group_col="length_bin", bins=LENGTH_LEVELS,
                  conditions=UNMATCHED_CONDITIONS,
                  title=paste0("APA scores — ", loop_type, " loops  |  unmatched  |  by loop length"))
  if (loop_type != "P-P") {
    plot_apa_diff(apa_agg=apa_unmatched, bins=LENGTH_LEVELS,
                  conditions=UNMATCHED_CONDITIONS,
                  title=paste0("APA diff — ", loop_type, "  |  unmatched  |  by loop length"))
  }
}

message("\n── Analysis 10d: CTCF loop types (unmatched) ────────────────────────")
ctcf_per_cond <- lapply(UNMATCHED_CONDITIONS, function(n) {
  np_classified[[n]] %>%
    filter(anchor1_ctcf | anchor2_ctcf) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        anchor1_ctcf & anchor2_ctcf ~ "CTCF-CTCF",
        (anchor1_ctcf & anchor1_ccre == "Promoter") |
          (anchor2_ctcf & anchor2_ccre == "Promoter") ~ "CTCF-P",
        (anchor1_ctcf & anchor1_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")) |
          (anchor2_ctcf & anchor2_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")) ~ "CTCF-E",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group)) %>%
    dplyr::mutate(group = factor(group, levels=CTCF_LEVELS))
}); names(ctcf_per_cond) <- UNMATCHED_CONDITIONS

apa_raw_ctcf_u <- lapply(UNMATCHED_CONDITIONS, function(n) {
  extract_apa_matrices_by(ctcf_per_cond[[n]], CONDITIONS[[n]]$hic,
                          group_col="group", group_levels=CTCF_LEVELS)
}); names(apa_raw_ctcf_u) <- UNMATCHED_CONDITIONS
apa_10d <- setNames(lapply(UNMATCHED_CONDITIONS, function(n)
  aggregate_apa(apa_raw_ctcf_u[[n]], DEPTH_SCALE[[n]])), UNMATCHED_CONDITIONS)
plot_apa_grid_generic(apa_agg=apa_10d, bin_levels=CTCF_LEVELS,
                      bin_colors=as.list(CTCF_COLORS),
                      conditions=UNMATCHED_CONDITIONS,
                      title="APA — CTCF loop types  |  unmatched")
plot_apa_scores(apa_agg=apa_10d, loops_annotated_cond=ctcf_per_cond$wt,
                group_col="group", bins=CTCF_LEVELS,
                conditions=UNMATCHED_CONDITIONS,
                title="APA scores — CTCF loop types  |  unmatched")
