##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Scatter plots -- loop length vs. contact strength. Anchor
# classification via narrowPeak (ChIP-seq + ATAC), CTCF pair
# sub-type (CTCF-CTCF/CTCF-Promoter/CTCF-Enhancer), WT/POLR2A-degraded
# exclusive and shared loop sets.
#
# Figure output: SI Fig. 9A (probable); SI Fig. 9B-D (confirmed --
# exact CTCF-CTCF/CTCF-Promoter/CTCF-Enhancer category and color match).
##################################################################
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)
library(tidyr)
library(GenomicRanges)
library(viridis)
library(cowplot)
source("~/Transcription Paper/Hic_loop_analysis/final_scripts/EP_analysis_functions.R")

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════

# ── Downsampling group switch ──────────────────────────────────────────────────
DS_GROUP <- "ds_wt_pol2"   # ← "ds_wt_pol2" | "ds_all" | "unsampled"

BASE     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
LOOP_DIR <- file.path(BASE, "loops")

if (DS_GROUP == "ds_wt_pol2") {
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE",
                hic  =file.path(BASE,"mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                loops=file.path(LOOP_DIR,"ds_wt_pol2/WT_HCT116_CTRL/mega/merged_loops.bedpe")),
    pol2 = list(label="POLR2A deg", color="#CC6677",
                hic  =file.path(BASE,"mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                loops=file.path(LOOP_DIR,"ds_wt_pol2/Pol2_6hrs_Aux/mega/merged_loops.bedpe"))
  )
  DEPTH_SCALE  <- list(wt=1.0, pol2=1.0)
  MIN_OBSERVED <- list(wt=0,   pol2=0)
  
} else if (DS_GROUP == "ds_all") {
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE",
                hic  =file.path(BASE,"mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                loops=file.path(LOOP_DIR,"ds_all/WT_HCT116_CTRL/mega/merged_loops.bedpe")),
    pol2 = list(label="POLR2A deg", color="#CC6677",
                hic  =file.path(BASE,"mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                loops=file.path(LOOP_DIR,"ds_all/Pol2_6hrs_Aux/mega/merged_loops.bedpe")),
    actd = list(label="ActD 1hr",   color="#44AA99",
                hic  =file.path(BASE,"mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic"),
                loops=file.path(LOOP_DIR,"ds_all/1hr_ActD/mega/merged_loops.bedpe"))
  )
  DEPTH_SCALE  <- list(wt=1.0, pol2=1.0, actd=1.0)
  MIN_OBSERVED <- list(wt=0,   pol2=0,   actd=0)
  
} else if (DS_GROUP == "unsampled") {
  wt_total   <- 639409109 + 377220354
  pol2_total <- 374163054 + 249677343
  actd_total <- 266069792 + 202465027
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE",
                hic  =file.path(BASE,"mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                loops=file.path(LOOP_DIR,"unsampled/WT_HCT116_CTRL/mega/merged_loops.bedpe")),
    pol2 = list(label="POLR2A deg", color="#CC6677",
                hic  =file.path(BASE,"mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                loops=file.path(LOOP_DIR,"unsampled/Pol2_6hrs_Aux/mega/merged_loops.bedpe")),
    actd = list(label="ActD 1hr",   color="#44AA99",
                hic  =file.path(BASE,"mega_unsampled/1hr_ActD/mega/aligned/inter_30.hic"),
                loops=file.path(LOOP_DIR,"unsampled/1hr_ActD/mega/merged_loops.bedpe"))
  )
  DEPTH_SCALE  <- list(wt=1.0, pol2=wt_total/pol2_total, actd=wt_total/actd_total)
  MIN_OBSERVED <- list(wt=10,  pol2=25, actd=25)
  
} else {
  stop("Unknown DS_GROUP: '", DS_GROUP, "'")
}

message("DS_GROUP: ", DS_GROUP)

# ── Score mode switch ──────────────────────────────────────────────────────────
# "fold_enrichment" — observed / expectedDonut (controls for distance decay)
# "observed_scaled" — SCALE contacts depth-scaled to WT
SCORE_MODE <- "observed_scaled"   # ← change to switch

if (SCORE_MODE == "fold_enrichment") {
  SCORE_COL      <- "fold_enrichment"
  OBS_AXIS_LABEL <- "Loop fold enrichment (log10 obs/exp)"
  OBS_SUBTITLE   <- paste0("HiCCUPS obs/expectedDonut  |  ", DS_GROUP)
} else if (SCORE_MODE == "observed_scaled") {
  SCORE_COL      <- "observed_scaled"
  OBS_AXIS_LABEL <- "Observed contacts (log2, depth-scaled to WT)"
  OBS_SUBTITLE   <- paste0("HiCCUPS observed contacts  |  ", DS_GROUP)
} else {
  stop("Unknown SCORE_MODE: '", SCORE_MODE, "'")
}
message("SCORE_MODE: ", SCORE_MODE)

# ── Analysis parameters ────────────────────────────────────────────────────────
GENE_TOL_BP    <- 5000
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1
OBS_FILTER     <- 0

# ── File paths ─────────────────────────────────────────────────────────────────
EXPR_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
DESEQ_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
CHIP_DIR   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"
ATAC_DIR   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"

# ── Derived constants ──────────────────────────────────────────────────────────
COND_LEVELS   <- sapply(CONDITIONS, `[[`, "label")
COND_COLORS   <- setNames(sapply(CONDITIONS, `[[`, "color"), COND_LEVELS)
BINS          <- c("Low","Medium","High")
BIN_COLORS_CB <- c("Low"="#E69F00", "Medium"="#56B4E9", "High"="#009E73")
DE_COLORS_CB  <- c("Up"="#D55E00",  "Down"="#0072B2",  "Not DE"="#999999")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load expression data and DESeq2 results
# ══════════════════════════════════════════════════════════════════════════════
expr <- read_tsv(EXPR_FILE) %>% dplyr::select(gene_id, TPM)
deseq_res <- read_csv(DESEQ_FILE) %>%
  dplyr::mutate(
    gene_id   = sub("\\..*", "", gene_id),
    de_status = dplyr::case_when(
      is.na(padj)                                               ~ "Filtered",
      padj >= PADJ_THRESHOLD                                    ~ "Not DE",
      padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD  ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD  ~ "Down",
      TRUE                                                      ~ "Not DE"
    ),
    de_status = factor(de_status, levels=c("Up","Not DE","Down","Filtered")),
    symbol    = dplyr::case_when(!is.na(symbol) & symbol != "" ~ symbol,
                                 TRUE ~ gene_id)
  ) %>% left_join(expr, by="gene_id")

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

# ══════════════════════════════════════════════════════════════════════════════
# 2. Build gene GRanges
# ══════════════════════════════════════════════════════════════════════════════
gene_coords <- deseq_res %>%
  filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name)) %>%
  distinct(gene_id, .keep_all=TRUE)
gene_gr <- GRanges(
  seqnames = as.character(gene_coords$seq_name),
  ranges   = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end),
  symbol   = gene_coords$symbol,
  gene_id  = gene_coords$gene_id,
  de_status= as.character(gene_coords$de_status),
  expr_bin = as.character(gene_coords$expr_bin),
  TPM      = gene_coords$TPM
)
message("gene_gr: ", length(gene_gr), " genes")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Loop loader — HiCCUPS with fold enrichment
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
    chr1          = as.character(df[[1]]),
    start1        = as.integer(df[[2]]),
    end1          = as.integer(df[[3]]),
    chr2          = as.character(df[[4]]),
    start2        = as.integer(df[[5]]),
    end2          = as.integer(df[[6]]),
    observed      = as.numeric(df[["observed"]]),
    expectedDonut = as.numeric(df[["expectedDonut"]]),
    fdrDonut      = as.numeric(df[["fdrDonut"]]),
    fdrBL         = as.numeric(df[["fdrBL"]]),
    fdrH          = as.numeric(df[["fdrH"]]),
    fdrV          = as.numeric(df[["fdrV"]]),
    centroid1     = as.integer(df[["centroid1"]]),
    centroid2     = as.integer(df[["centroid2"]])
  ) %>% filter(!is.na(start1), !is.na(start2))
}

filter_loops <- function(df, min_observed=10, max_fdr=0.1) {
  df %>% filter(observed >= min_observed,
                fdrDonut <= max_fdr, fdrBL <= max_fdr,
                fdrH <= max_fdr, fdrV <= max_fdr)
}

message("\nLoading HiCCUPS loops (", DS_GROUP, ")...")
loops_raw <- lapply(names(CONDITIONS), function(n) {
  path <- CONDITIONS[[n]]$loops
  if (!file.exists(path)) stop("Loop file not found: ", path)
  filter_loops(load_hiccups_bedpe(path), min_observed=MIN_OBSERVED[[n]]) %>%
    dplyr::mutate(
      loop_idx        = row_number(),
      observed_scaled = observed * DEPTH_SCALE[[n]],
      fold_enrichment = observed / expectedDonut
    )
}); names(loops_raw) <- names(CONDITIONS)
invisible(lapply(names(loops_raw), function(n)
  message(CONDITIONS[[n]]$label, " loops: ", nrow(loops_raw[[n]]))))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Link loops to genes
# ══════════════════════════════════════════════════════════════════════════════
link_loops_to_genes <- function(loops_df, gene_gr, tol_bp=GENE_TOL_BP) {
  fix_chr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr","",seqlevels(gr))
    gr
  }
  anchor1_gr <- fix_chr(GRanges(loops_df$chr1, IRanges(loops_df$start1, loops_df$end1)), gene_gr)
  anchor2_gr <- fix_chr(GRanges(loops_df$chr2, IRanges(loops_df$start2, loops_df$end2)), gene_gr)
  gene_gr_exp <- if (tol_bp > 0)
    trim(resize(gene_gr, width(gene_gr) + 2*tol_bp, fix="center")) else gene_gr
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

message("\nLinking loops to genes...")
loops_annotated <- lapply(loops_raw, function(ldf) link_loops_to_genes(ldf, gene_gr))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Prepare scatter data
# ══════════════════════════════════════════════════════════════════════════════
prepare_scatter_data <- function(loops_annotated_list) {
  bind_rows(lapply(names(CONDITIONS), function(n) {
    loops_annotated_list[[n]] %>%
      as.data.frame() %>%
      filter(expr_bin %in% BINS) %>%
      dplyr::mutate(
        condition = CONDITIONS[[n]]$label,
        loop_kb   = abs(((start1+end1)/2) - ((start2+end2)/2)) / 1000,
        score     = .data[[SCORE_COL]],
        log2_obs  = log2(score + 1),
        expr_bin  = factor(expr_bin,   levels=BINS),
        de_status = factor(de_status,  levels=c("Up","Not DE","Down","Filtered")),
        condition = factor(condition,  levels=sapply(CONDITIONS, `[[`, "label"))
      ) %>%
      filter(!is.na(loop_kb), !is.na(score), score > OBS_FILTER) %>%
      dplyr::select(condition, loop_kb, observed, observed_scaled,
                    fold_enrichment, score, log2_obs,
                    expr_bin, de_status, symbol, TPM)
  }))
}

scatter_df <- prepare_scatter_data(loops_annotated)
message("Scatter data: ", nrow(scatter_df), " loop-condition pairs")
print(table(scatter_df$condition))


# ══════════════════════════════════════════════════════════════════════════════
# 7. Hexbin scatter with marginal boxplots — expression bins
# ══════════════════════════════════════════════════════════════════════════════
plot_loop_hexbin_box_tpm <- function(scatter_df, cap_kb=NULL, bins=60,
                                     title="Loop size vs contact strength") {
  plot_df <- scatter_df %>%
    filter(!is.na(loop_kb), !is.na(score), score > 0, expr_bin %in% BINS) %>%
    dplyr::mutate(log10_size=log10(loop_kb), log10_obs=log10(score),
                  expr_bin=factor(expr_bin, levels=BINS))
  if (!is.null(cap_kb)) plot_df <- plot_df %>% filter(loop_kb <= cap_kb)
  x_lims <- range(plot_df$log10_size, na.rm=TRUE)
  y_lims <- range(plot_df$log10_obs,  na.rm=TRUE)
  make_panel <- function(cond_label) {
    df       <- plot_df %>% filter(condition == cond_label)
    cond_key <- names(which(sapply(CONDITIONS, `[[`, "label") == cond_label))
    color    <- CONDITIONS[[cond_key]]$color
    med_x    <- median(df$log10_size, na.rm=TRUE)
    med_y    <- median(df$log10_obs,  na.rm=TRUE)
    whisker_stats <- df %>%
      dplyr::group_by(expr_bin) %>%
      dplyr::summarise(
        q3_size=quantile(log10_size,0.75,na.rm=TRUE), iqr_size=IQR(log10_size,na.rm=TRUE),
        max_size=max(log10_size,na.rm=TRUE), whisker_size=min(max_size, q3_size+1.5*iqr_size),
        med_size_raw=median(10^log10_size,na.rm=TRUE),
        q3_obs=quantile(log10_obs,0.75,na.rm=TRUE), iqr_obs=IQR(log10_obs,na.rm=TRUE),
        max_obs=max(log10_obs,na.rm=TRUE), whisker_obs=min(max_obs, q3_obs+1.5*iqr_obs),
        med_obs_raw=median(10^log10_obs,na.rm=TRUE), .groups="drop"
      ) %>%
      dplyr::mutate(expr_bin=factor(expr_bin, levels=BINS),
                    label_size=paste0(round(med_size_raw), " kb"),
                    label_obs=as.character(round(med_obs_raw, 2)))
    p_main <- ggplot(df, aes(x=log10_size, y=log10_obs)) +
      geom_hex(bins=bins) +
      scale_fill_viridis_c(option="inferno", name="Count",
                           trans="log10", limits=c(1,100), oob=scales::squish) +
      geom_vline(xintercept=med_x, linetype="dashed", color="#CC0000", linewidth=0.5) +
      geom_hline(yintercept=med_y, linetype="dashed", color="#CC0000", linewidth=0.5) +
      scale_x_continuous(name="Loop size (log10 kb)", limits=x_lims, expand=c(0.02,0)) +
      scale_y_continuous(name=OBS_AXIS_LABEL, limits=y_lims, expand=c(0.02,0)) +
      annotate("text", x=Inf, y=Inf, label=cond_label,
               hjust=1.1, vjust=1.5, size=4.5, fontface="bold", color=color) +
      theme_classic(base_size=11) +
      theme(legend.position=c(0.98,0.02), legend.justification=c(1,0),
            legend.key.size=unit(0.35,"cm"),
            legend.background=element_rect(fill=alpha("white",0.7),color=NA),
            plot.margin=margin(0,0,0,0))
    p_top <- ggplot(df, aes(x=log10_size, y=expr_bin, fill=expr_bin, color=expr_bin)) +
      geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
      geom_text(data=whisker_stats, aes(x=whisker_size, y=expr_bin, label=label_size),
                inherit.aes=FALSE, hjust=-0.2, vjust=0.5, size=2.5, color="gray30") +
      geom_vline(xintercept=med_x, linetype="dashed", color="#CC0000", linewidth=0.5) +
      scale_fill_manual(values=BIN_COLORS_CB, guide="none") +
      scale_color_manual(values=BIN_COLORS_CB, guide="none") +
      scale_x_continuous(limits=c(x_lims[1], x_lims[2]+0.8), expand=c(0.02,0)) +
      scale_y_discrete(limits=rev(BINS), labels=rev(BINS), expand=expansion(add=c(0.5,0.5))) +
      theme_void() +
      theme(axis.text.y=element_text(size=7, color="gray30", hjust=1), plot.margin=margin(2,0,0,0))
    p_right <- ggplot(df, aes(x=log10_obs, y=expr_bin, fill=expr_bin, color=expr_bin)) +
      geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
      geom_text(data=whisker_stats, aes(x=whisker_obs, y=expr_bin, label=label_obs),
                inherit.aes=FALSE, hjust=-0.2, vjust=0.5, size=2.5, color="gray30", angle=90) +
      geom_vline(xintercept=med_y, linetype="dashed", color="#CC0000", linewidth=0.5) +
      scale_fill_manual(values=BIN_COLORS_CB, guide="none") +
      scale_color_manual(values=BIN_COLORS_CB, guide="none") +
      scale_x_continuous(limits=c(y_lims[1], y_lims[2]+0.4), expand=c(0.02,0)) +
      scale_y_discrete(limits=rev(BINS), labels=rev(BINS), expand=expansion(add=c(0.5,0.5))) +
      coord_flip() + theme_void() + theme(plot.margin=margin(0,2,0,0))
    p_top + plot_spacer() + p_main + p_right +
      plot_layout(ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
  }
  panels    <- lapply(sapply(CONDITIONS, `[[`, "label"), make_panel)
  row       <- cowplot::plot_grid(plotlist=panels, nrow=1,
                                  labels=letters[seq_along(panels)], label_size=14)
  legend_df <- data.frame(x=1:3, y=1, expr_bin=factor(BINS, levels=BINS))
  legend_plot <- ggplot(legend_df, aes(x=x, y=y, fill=expr_bin, color=expr_bin)) +
    geom_point(size=4, shape=22, stroke=1) +
    scale_fill_manual(values=BIN_COLORS_CB, name="Expression bin") +
    scale_color_manual(values=BIN_COLORS_CB, name="Expression bin") +
    theme_void() + theme(legend.position="bottom")
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0, y=0.75) +
    cowplot::draw_label(paste0(OBS_SUBTITLE, "  |  marginal boxplots by expression bin"),
                        size=8, color="gray40", x=0, hjust=0, y=0.25)
  cowplot::plot_grid(title_grob, row, cowplot::get_legend(legend_plot),
                     ncol=1, rel_heights=c(0.08, 1, 0.07))
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Hexbin scatter with marginal boxplots — DE status
# ══════════════════════════════════════════════════════════════════════════════
plot_loop_hexbin_box_de_status <- function(scatter_df, cap_kb=NULL, bins=60,
                                           title="Loop size vs contact strength — DE genes") {
  DE_COLORS <- c("Up"="#D55E00", "Down"="#0072B2")
  DE_LEVELS <- c("Up","Down")
  plot_df <- scatter_df %>%
    filter(!is.na(loop_kb), !is.na(score), score > 0, de_status %in% DE_LEVELS) %>%
    dplyr::mutate(log10_size=log10(loop_kb), log10_obs=log10(score),
                  de_status=factor(de_status, levels=DE_LEVELS))
  if (!is.null(cap_kb)) plot_df <- plot_df %>% filter(loop_kb <= cap_kb)
  x_lims <- range(plot_df$log10_size, na.rm=TRUE)
  y_lims <- range(plot_df$log10_obs,  na.rm=TRUE)
  make_panel <- function(cond_label) {
    df       <- plot_df %>% filter(condition == cond_label)
    cond_key <- names(which(sapply(CONDITIONS, `[[`, "label") == cond_label))
    color    <- CONDITIONS[[cond_key]]$color
    med_x    <- median(df$log10_size, na.rm=TRUE)
    med_y    <- median(df$log10_obs,  na.rm=TRUE)
    whisker_stats <- df %>%
      dplyr::group_by(de_status) %>%
      dplyr::summarise(
        q3_size=quantile(log10_size,0.75,na.rm=TRUE), iqr_size=IQR(log10_size,na.rm=TRUE),
        max_size=max(log10_size,na.rm=TRUE), whisker_size=min(max_size, q3_size+1.5*iqr_size),
        med_size_raw=median(10^log10_size,na.rm=TRUE),
        q3_obs=quantile(log10_obs,0.75,na.rm=TRUE), iqr_obs=IQR(log10_obs,na.rm=TRUE),
        max_obs=max(log10_obs,na.rm=TRUE), whisker_obs=min(max_obs, q3_obs+1.5*iqr_obs),
        med_obs_raw=median(10^log10_obs,na.rm=TRUE), .groups="drop"
      ) %>%
      dplyr::mutate(de_status=factor(de_status, levels=DE_LEVELS),
                    label_size=paste0(round(med_size_raw), " kb"),
                    label_obs=as.character(round(med_obs_raw, 2)))
    p_main <- ggplot(df, aes(x=log10_size, y=log10_obs)) +
      geom_hex(bins=bins) +
      scale_fill_viridis_c(option="inferno", name="Count",
                           trans="log10", limits=c(1,100), oob=scales::squish) +
      geom_vline(xintercept=med_x, linetype="dashed", color="#CC0000", linewidth=0.5) +
      geom_hline(yintercept=med_y, linetype="dashed", color="#CC0000", linewidth=0.5) +
      scale_x_continuous(name="Loop size (log10 kb)", limits=x_lims, expand=c(0.02,0)) +
      scale_y_continuous(name=OBS_AXIS_LABEL, limits=y_lims, expand=c(0.02,0)) +
      annotate("text", x=Inf, y=Inf, label=cond_label,
               hjust=1.1, vjust=1.5, size=4.5, fontface="bold", color=color) +
      theme_classic(base_size=11) +
      theme(legend.position=c(0.98,0.02), legend.justification=c(1,0),
            legend.key.size=unit(0.35,"cm"),
            legend.background=element_rect(fill=alpha("white",0.7),color=NA),
            plot.margin=margin(0,0,0,0))
    p_top <- ggplot(df, aes(x=log10_size, y=de_status, fill=de_status, color=de_status)) +
      geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
      geom_text(data=whisker_stats, aes(x=whisker_size, y=de_status, label=label_size),
                inherit.aes=FALSE, hjust=-0.2, vjust=0.5, size=2.5, color="gray30") +
      geom_vline(xintercept=med_x, linetype="dashed", color="#CC0000", linewidth=0.5) +
      scale_fill_manual(values=DE_COLORS, guide="none") +
      scale_color_manual(values=DE_COLORS, guide="none") +
      scale_x_continuous(limits=c(x_lims[1], x_lims[2]+0.8), expand=c(0.02,0)) +
      scale_y_discrete(limits=rev(DE_LEVELS), labels=rev(DE_LEVELS),
                       expand=expansion(add=c(0.5,0.5))) +
      theme_void() +
      theme(axis.text.y=element_text(size=7, color="gray30", hjust=1), plot.margin=margin(2,0,0,0))
    p_right <- ggplot(df, aes(x=log10_obs, y=de_status, fill=de_status, color=de_status)) +
      geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
      geom_text(data=whisker_stats, aes(x=whisker_obs, y=de_status, label=label_obs),
                inherit.aes=FALSE, hjust=-0.2, vjust=0.5, size=2.5, color="gray30", angle=90) +
      geom_vline(xintercept=med_y, linetype="dashed", color="#CC0000", linewidth=0.5) +
      scale_fill_manual(values=DE_COLORS, guide="none") +
      scale_color_manual(values=DE_COLORS, guide="none") +
      scale_x_continuous(limits=c(y_lims[1], y_lims[2]+0.4), expand=c(0.02,0)) +
      scale_y_discrete(limits=rev(DE_LEVELS), labels=rev(DE_LEVELS),
                       expand=expansion(add=c(0.5,0.5))) +
      coord_flip() + theme_void() + theme(plot.margin=margin(0,2,0,0))
    p_top + plot_spacer() + p_main + p_right +
      plot_layout(ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
  }
  panels <- lapply(sapply(CONDITIONS, `[[`, "label"), make_panel)
  row    <- cowplot::plot_grid(plotlist=panels, nrow=1,
                               labels=letters[seq_along(panels)], label_size=14)
  de_legend_df <- data.frame(x=1:2, y=1, de_status=factor(DE_LEVELS, levels=DE_LEVELS))
  de_legend_plot <- ggplot(de_legend_df, aes(x=x, y=y, fill=de_status, color=de_status)) +
    geom_point(size=4, shape=22, stroke=1) +
    scale_fill_manual(values=DE_COLORS, name="DE status") +
    scale_color_manual(values=DE_COLORS, name="DE status") +
    theme_void() + theme(legend.position="bottom")
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0, y=0.75) +
    cowplot::draw_label(paste0(OBS_SUBTITLE, "  |  marginal boxplots by DE status"),
                        size=8, color="gray40", x=0, hjust=0, y=0.25)
  cowplot::plot_grid(title_grob, row, cowplot::get_legend(de_legend_plot),
                     ncol=1, rel_heights=c(0.08, 1, 0.07))
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Run expression and DE scatter plots
# ══════════════════════════════════════════════════════════════════════════════
p_hexbin_box_tpm <- plot_loop_hexbin_box_tpm(
  scatter_df, title=paste0("Loop size vs contact strength — ", DS_GROUP))
print(p_hexbin_box_tpm)

p_hexbin_box_de <- plot_loop_hexbin_box_de_status(
  scatter_df, title=paste0("Loop size vs contact strength — DE genes — ", DS_GROUP))
print(p_hexbin_box_de)

# ══════════════════════════════════════════════════════════════════════════════
# 9b. Strict DEG scatter — both anchors on a DE gene
# ══════════════════════════════════════════════════════════════════════════════
link_loops_to_deg_strict <- function(loops_df, gene_gr, tol_bp=GENE_TOL_BP) {
  fix_chr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr","",seqlevels(gr))
    gr
  }
  anchor1_gr <- fix_chr(GRanges(loops_df$chr1, IRanges(loops_df$start1, loops_df$end1)), gene_gr)
  anchor2_gr <- fix_chr(GRanges(loops_df$chr2, IRanges(loops_df$start2, loops_df$end2)), gene_gr)
  de_gene_gr <- gene_gr[gene_gr$de_status %in% c("Up","Down")]
  gene_gr_exp <- if (tol_bp > 0)
    trim(resize(de_gene_gr, width(de_gene_gr) + 2*tol_bp, fix="center")) else de_gene_gr
  hits1 <- findOverlaps(anchor1_gr, gene_gr_exp)
  hits2 <- findOverlaps(anchor2_gr, gene_gr_exp)
  both_de <- intersect(unique(queryHits(hits1)), unique(queryHits(hits2)))
  message("  Loops with both anchors on DE genes: ", length(both_de))
  if (length(both_de) == 0) return(loops_df %>%
                                     dplyr::mutate(expr_bin="Unannotated", de_status="Not DE", symbol="Unknown", TPM=NA))
  all_pairs <- bind_rows(
    data.frame(loop_idx=queryHits(hits1)[queryHits(hits1) %in% both_de],
               gene_idx=subjectHits(hits1)[queryHits(hits1) %in% both_de]),
    data.frame(loop_idx=queryHits(hits2)[queryHits(hits2) %in% both_de],
               gene_idx=subjectHits(hits2)[queryHits(hits2) %in% both_de])
  ) %>% distinct(loop_idx, gene_idx) %>%
    dplyr::mutate(expr_bin=de_gene_gr$expr_bin[gene_idx],
                  de_status=de_gene_gr$de_status[gene_idx],
                  TPM=de_gene_gr$TPM[gene_idx], symbol=de_gene_gr$symbol[gene_idx])
  loop_annotations <- all_pairs %>%
    dplyr::group_by(loop_idx) %>%
    dplyr::slice_max(TPM, n=1, with_ties=FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(loop_idx, expr_bin, de_status, TPM, symbol)
  loops_df %>%
    left_join(loop_annotations, by="loop_idx") %>%
    dplyr::mutate(expr_bin=replace_na(expr_bin,"Unannotated"),
                  de_status=replace_na(de_status,"Not DE"),
                  symbol=replace_na(symbol,"Unknown")) %>%
    filter(loop_idx %in% both_de)
}

message("\nLinking loops to DE genes — strict (both anchors)...")
loops_annotated_strict <- lapply(loops_raw, function(ldf)
  link_loops_to_deg_strict(ldf, gene_gr))
scatter_df_strict <- prepare_scatter_data(loops_annotated_strict)
message("Strict DEG scatter: ", nrow(scatter_df_strict), " loop-condition pairs")

p_hexbin_box_de_strict <- plot_loop_hexbin_box_de_status(
  scatter_df_strict,
  title=paste0("Loop size vs contact strength — both anchors on DE genes — ", DS_GROUP))
print(p_hexbin_box_de_strict)

# ══════════════════════════════════════════════════════════════════════════════
# 10. cCRE scatter — WT vs Pol2 only (narrowPeak + ATAC, condition-matched)
# ══════════════════════════════════════════════════════════════════════════════
stopifnot("wt" %in% names(CONDITIONS), "pol2" %in% names(CONDITIONS))

message("\nLoading narrowPeak files...")
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

library(AnnotationDbi)
txdb <- loadDb("/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/GRCh38.115.sqlite")

message("Classifying WT loops (WT peaks)...")
loops_np_wt <- classify_anchors_narrowpeak(
  loops_df    = loops_raw$wt,
  atac_gr     = atac_wt,
  h3k4me3_gr  = h3k4me3_wt,
  h3k27ac_gr  = h3k27ac_wt,
  h3k4me1_gr  = h3k4me1_wt,
  ctcf_gr     = ctcf_wt,
  txdb        = txdb,
  tss_region  = c(-3000, 3000),
  ccre_label  = "WT_narrowpeak"
)
message("Classifying pol2 loops (KO peaks)...")
loops_np_pol2 <- classify_anchors_narrowpeak(
  loops_df    = loops_raw$pol2,
  atac_gr     = atac_ko,
  h3k4me3_gr  = h3k4me3_ko,
  h3k27ac_gr  = h3k27ac_ko,
  h3k4me1_gr  = h3k4me1_ko,
  ctcf_gr     = ctcf_ko,
  txdb        = txdb,
  tss_region  = c(-3000, 3000),
  ccre_label  = "KO_narrowpeak"
)


NP_COLORS <- c(
  "Promoter"          = "#e31a1c",
  "Proximal-Enhancer" = "#ff7f00",
  "Distal-Enhancer"   = "#fdbf6f",
  "CTCF"              = "#1f78b4"
)
NP_LABELS <- c(
  "Promoter"          = "Promoter",
  "Proximal-Enhancer" = "Proximal enhancer",
  "Distal-Enhancer"   = "Distal enhancer",
  "CTCF"              = "CTCF-bound"
)
NP_PRIORITY <- c("Promoter"=1, "Proximal-Enhancer"=2, "Distal-Enhancer"=3, "CTCF"=4)

build_np_scatter_df <- function(classified, cond_name) {
  score_vals <- if (SCORE_MODE == "fold_enrichment") {
    classified$observed / classified$expectedDonut
  } else {
    classified$observed * DEPTH_SCALE[[cond_name]]
  }
  classified %>%
    dplyr::mutate(score_val = score_vals) %>%
    filter(!is.na(score_val), score_val > OBS_FILTER,
           !is.na(anchor1_ccre) | !is.na(anchor2_ccre)) %>%
    dplyr::mutate(
      condition  = CONDITIONS[[cond_name]]$label,
      loop_kb    = abs(((start1+end1)/2) - ((start2+end2)/2)) / 1000,
      log10_obs  = log10(score_val),
      log10_size = log10(loop_kb),
      a1_rank    = NP_PRIORITY[anchor1_ccre],
      a2_rank    = NP_PRIORITY[anchor2_ccre],
      ccre_class = dplyr::case_when(
        is.na(a1_rank) & is.na(a2_rank) ~ NA_character_,
        is.na(a1_rank)                  ~ anchor2_ccre,
        is.na(a2_rank)                  ~ anchor1_ccre,
        a1_rank <= a2_rank              ~ anchor1_ccre,
        TRUE                            ~ anchor2_ccre
      ),
      condition = factor(condition, levels=COND_LEVELS)
    ) %>%
    filter(!is.na(loop_kb), loop_kb > 0, !is.na(ccre_class)) %>%
    dplyr::select(condition, loop_kb, log10_size, log10_obs,
                  ccre_class, anchor1_ccre, anchor2_ccre)
}

scatter_enh_np <- bind_rows(
  build_np_scatter_df(loops_np_wt,   "wt"),
  build_np_scatter_df(loops_np_pol2, "pol2")
) %>% filter(ccre_class %in% c("Promoter","Proximal-Enhancer","Distal-Enhancer"))

build_np_ctcf_scatter_df <- function(classified, cond_name) {
  score_vals <- if (SCORE_MODE == "fold_enrichment") {
    classified$observed / classified$expectedDonut
  } else {
    classified$observed * DEPTH_SCALE[[cond_name]]
  }
  classified %>%
    dplyr::mutate(score_val = score_vals) %>%
    filter(!is.na(score_val), score_val > OBS_FILTER,
           anchor1_ctcf | anchor2_ctcf) %>%
    dplyr::mutate(
      condition  = CONDITIONS[[cond_name]]$label,
      loop_kb    = abs(((start1+end1)/2) - ((start2+end2)/2)) / 1000,
      log10_obs  = log10(score_val),
      log10_size = log10(loop_kb),
      ccre_class = "CTCF",                          # ← allows plot_np_scatter to work
      condition  = factor(condition, levels=COND_LEVELS)
    ) %>%
    filter(!is.na(loop_kb), loop_kb > 0) %>%
    dplyr::select(condition, loop_kb, log10_size, log10_obs,
                  ccre_class, anchor1_ccre, anchor2_ccre, anchor1_ctcf, anchor2_ctcf)
}

scatter_ctcf_np <- bind_rows(
  build_np_ctcf_scatter_df(loops_np_wt,   "wt"),
  build_np_ctcf_scatter_df(loops_np_pol2, "pol2")
)

all_data_np <- bind_rows(scatter_enh_np, scatter_ctcf_np)
x_lims_np <- c(range(all_data_np$log10_size, na.rm=TRUE)[1], log10(10000))
y_lims_np <- c(range(all_data_np$log10_obs, na.rm=TRUE)[1], 1.5)

plot_np_scatter <- function(scatter_df, active_classes, x_lims, y_lims,
                            title, subtitle, alpha=0.35, point_size=0.6) {
  plot_df <- scatter_df %>%
    filter(ccre_class %in% active_classes) %>%
    dplyr::mutate(ccre_class=factor(ccre_class, levels=active_classes))
  whisker_stats <- plot_df %>%
    dplyr::group_by(condition, ccre_class) %>%
    dplyr::summarise(
      q3_size=quantile(log10_size,0.75,na.rm=TRUE), iqr_size=IQR(log10_size,na.rm=TRUE),
      max_size=max(log10_size,na.rm=TRUE), whisker_size=min(max_size, q3_size+1.5*iqr_size),
      med_size_raw=median(10^log10_size,na.rm=TRUE),
      q3_obs=quantile(log10_obs,0.75,na.rm=TRUE), iqr_obs=IQR(log10_obs,na.rm=TRUE),
      max_obs=max(log10_obs,na.rm=TRUE), whisker_obs=min(max_obs, q3_obs+1.5*iqr_obs),
      med_obs_raw=median(10^log10_obs,na.rm=TRUE), .groups="drop"
    ) %>%
    dplyr::mutate(ccre_class=factor(ccre_class, levels=active_classes),
                  label_size=paste0(round(med_size_raw), " kb"),
                  label_obs=as.character(round(med_obs_raw, 2)))
  make_panel <- function(cond_label) {
    df         <- plot_df %>% filter(condition == cond_label)
    ws         <- whisker_stats %>% filter(condition == cond_label)
    cond_key   <- names(which(sapply(CONDITIONS, `[[`, "label") == cond_label))
    cond_color <- CONDITIONS[[cond_key]]$color
    med_x <- median(df$log10_size, na.rm=TRUE)
    med_y <- median(df$log10_obs,  na.rm=TRUE)
    class_counts <- df %>%
      dplyr::count(ccre_class) %>%
      dplyr::mutate(ccre_class=factor(ccre_class, levels=active_classes),
                    ann_label=paste0("n=", n)) %>%
      dplyr::arrange(ccre_class) %>%
      dplyr::mutate(ann_y=y_lims[2]-diff(y_lims)*0.12-
                      (seq_len(dplyr::n())-1)*diff(y_lims)*0.05)
    p_main <- ggplot(df, aes(x=log10_size, y=log10_obs, color=ccre_class)) +
      geom_vline(xintercept=med_x, linetype="dashed", color="gray60", linewidth=0.4) +
      geom_hline(yintercept=med_y, linetype="dashed", color="gray60", linewidth=0.4) +
      geom_point(alpha=alpha, size=point_size) +
      geom_text(data=class_counts,
                aes(x=x_lims[2], y=ann_y, label=ann_label, color=ccre_class),
                inherit.aes=FALSE, hjust=1, vjust=0.5, size=2.5, fontface="bold") +
      scale_color_manual(values=NP_COLORS[active_classes], name="Anchor type",
                         labels=NP_LABELS[active_classes]) +
      scale_x_continuous(limits=x_lims, expand=c(0.02,0),
                         breaks=log10(c(10,100,1000,10000,100000)),
                         labels=c("10 kb","100 kb","1 Mb","10 Mb","100 Mb"), name="Loop size") +
      scale_y_continuous(limits=y_lims, expand=c(0.02,0), name=OBS_AXIS_LABEL) +
      annotate("text", x=Inf, y=Inf, label=cond_label,
               hjust=1.1, vjust=1.5, size=4.5, fontface="bold", color=cond_color) +
      theme_classic(base_size=11) +
      theme(legend.position="none", plot.margin=margin(0,0,0,0))
    p_top <- ggplot(df, aes(x=log10_size, y=ccre_class, fill=ccre_class, color=ccre_class)) +
      geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
      geom_text(data=ws, aes(x=whisker_size, y=ccre_class, label=label_size),
                inherit.aes=FALSE, hjust=-0.2, vjust=0.5, size=2.5, color="gray30") +
      geom_vline(xintercept=med_x, linetype="dashed", color="gray60", linewidth=0.4) +
      scale_fill_manual(values=NP_COLORS[active_classes], guide="none") +
      scale_color_manual(values=NP_COLORS[active_classes], guide="none") +
      scale_x_continuous(limits=c(x_lims[1], x_lims[2]+0.8), expand=c(0.02,0)) +
      scale_y_discrete(limits=rev(active_classes), labels=rev(NP_LABELS[active_classes]),
                       expand=expansion(add=c(0.5,0.5))) +
      theme_void() +
      theme(axis.text.y=element_text(size=7, color="gray30", hjust=1), plot.margin=margin(2,0,0,0))
    p_right <- ggplot(df, aes(x=log10_obs, y=ccre_class, fill=ccre_class, color=ccre_class)) +
      geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
      geom_text(data=ws, aes(x=whisker_obs, y=ccre_class, label=label_obs),
                inherit.aes=FALSE, hjust=-0.2, vjust=0.5, size=2.5, color="gray30", angle=90) +
      geom_vline(xintercept=med_y, linetype="dashed", color="gray60", linewidth=0.4) +
      scale_fill_manual(values=NP_COLORS[active_classes], guide="none") +
      scale_color_manual(values=NP_COLORS[active_classes], guide="none") +
      scale_x_continuous(limits=c(y_lims[1], y_lims[2]+0.3), expand=c(0.02,0)) +
      scale_y_discrete(limits=rev(active_classes), labels=rev(NP_LABELS[active_classes]),
                       expand=expansion(add=c(0.5,0.5))) +
      coord_flip() + theme_void() + theme(plot.margin=margin(0,2,0,0))
    p_top + plot_spacer() + p_main + p_right +
      plot_layout(ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
  }
  wt_pol2_labels <- sapply(CONDITIONS[c("wt","pol2")], `[[`, "label")
  panels  <- lapply(wt_pol2_labels, make_panel)
  row     <- cowplot::plot_grid(plotlist=panels, nrow=1, labels=c("a","b"), label_size=14)
  legend_df <- data.frame(x=seq_along(active_classes), y=1,
                          ccre_class=factor(active_classes, levels=active_classes))
  legend_plot <- ggplot(legend_df, aes(x=x, y=y, fill=ccre_class, color=ccre_class)) +
    geom_point(size=4, shape=22, stroke=1) +
    scale_fill_manual(values=NP_COLORS[active_classes], name="Anchor type",
                      labels=NP_LABELS[active_classes]) +
    scale_color_manual(values=NP_COLORS[active_classes], name="Anchor type",
                       labels=NP_LABELS[active_classes]) +
    theme_void() + theme(legend.position="bottom", legend.text=element_text(size=10),
                         legend.title=element_text(size=11, face="bold"))
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0, y=0.75) +
    cowplot::draw_label(subtitle, size=8, color="gray40", x=0, hjust=0, y=0.25)
  cowplot::plot_grid(title_grob, row, cowplot::get_legend(legend_plot),
                     ncol=1, rel_heights=c(0.10, 1, 0.07))
}

p_enh_np <- plot_np_scatter(
  scatter_enh_np, c("Promoter","Proximal-Enhancer","Distal-Enhancer"),
  x_lims_np, y_lims_np,
  title="Loop size vs contact strength — Promoter & Enhancer loops",
  subtitle=paste0("Condition-matched narrowPeak + ATAC  |  ", DS_GROUP,
                  "  |  obs > ", OBS_FILTER)
)
p_ctcf_np <- plot_np_scatter(
  scatter_ctcf_np, "CTCF",
  x_lims_np, y_lims_np,
  title="Loop size vs contact strength — CTCF loops",
  subtitle=paste0("WT peaks for both conditions  |  ", DS_GROUP,
                  "  |  obs > ", OBS_FILTER)
)
print(p_enh_np)
print(p_ctcf_np)

# ══════════════════════════════════════════════════════════════════════════════
# 11. Shared vs exclusive loop scatter — WT only | Shared | pol2 only
#     Marginal boxplots by expression bin with median + n labels
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Finding shared / exclusive loops (WT vs pol2) ────────────────────")

find_shared_loops <- function(loops_a, loops_b, slop=5000L) {
  make_gr <- function(df) list(
    a1=GRanges(df$chr1, IRanges(df$start1, df$end1)),
    a2=GRanges(df$chr2, IRanges(df$start2, df$end2))
  )
  gr_a <- make_gr(loops_a); gr_b <- make_gr(loops_b)
  expand <- function(gr) resize(gr, width(gr) + 2*slop, fix="center")
  hits_a1 <- findOverlaps(expand(gr_a$a1), expand(gr_b$a1))
  hits_a2 <- findOverlaps(expand(gr_a$a2), expand(gr_b$a2))
  shared <- dplyr::inner_join(
    as.data.frame(hits_a1) %>% dplyr::rename(a_idx=queryHits, b_idx=subjectHits),
    as.data.frame(hits_a2) %>% dplyr::rename(a_idx=queryHits, b_idx=subjectHits),
    by=c("a_idx","b_idx")
  ) %>% dplyr::distinct(a_idx, b_idx)
  list(
    shared   = shared,
    a_only   = setdiff(seq_len(nrow(loops_a)), shared$a_idx),
    b_only   = setdiff(seq_len(nrow(loops_b)), shared$b_idx),
    n_shared = nrow(shared),
    n_a_only = nrow(loops_a) - dplyr::n_distinct(shared$a_idx),
    n_b_only = nrow(loops_b) - dplyr::n_distinct(shared$b_idx)
  )
}

loop_overlap <- find_shared_loops(loops_raw$wt, loops_raw$pol2, slop=5000L)
message("WT loops:   ", nrow(loops_raw$wt))
message("pol2 loops: ", nrow(loops_raw$pol2))
message("Shared:     ", loop_overlap$n_shared,
        " (", round(100*loop_overlap$n_shared/nrow(loops_raw$wt),1), "% of WT)")
message("WT only:    ", loop_overlap$n_a_only)
message("pol2 only:  ", loop_overlap$n_b_only)

build_subset_scatter <- function() {
  df_wt_only <- loops_annotated$wt %>% as.data.frame() %>%
    dplyr::slice(loop_overlap$a_only) %>%
    dplyr::mutate(panel="WT only",
                  loop_kb=abs(((start1+end1)/2)-((start2+end2)/2))/1000,
                  score=.data[[SCORE_COL]])
  df_pol2_only <- loops_annotated$pol2 %>% as.data.frame() %>%
    dplyr::slice(loop_overlap$b_only) %>%
    dplyr::mutate(panel="pol2 only",
                  loop_kb=abs(((start1+end1)/2)-((start2+end2)/2))/1000,
                  score=.data[[SCORE_COL]])
  df_shared <- loops_annotated$wt %>% as.data.frame() %>%
    dplyr::slice(loop_overlap$shared$a_idx) %>%
    dplyr::mutate(panel="Shared",
                  loop_kb=abs(((start1+end1)/2)-((start2+end2)/2))/1000,
                  score=.data[[SCORE_COL]])
  bind_rows(df_wt_only, df_pol2_only, df_shared) %>%
    filter(expr_bin %in% BINS, !is.na(loop_kb), !is.na(score), score > 0) %>%
    dplyr::mutate(
      log10_size  = log10(loop_kb),
      log10_score = log10(score),
      expr_bin    = factor(expr_bin, levels=BINS),
      panel       = factor(panel, levels=c("WT only","Shared","pol2 only"))
    )
}

subset_df <- build_subset_scatter()
message("Subset scatter rows: ", nrow(subset_df))
print(table(subset_df$panel, subset_df$expr_bin))

# Cap axes: x at 10 Mb (10,000 kb), y at 1
x_lims_sub <- c(range(subset_df$log10_size,  na.rm=TRUE)[1], log10(10000))
y_lims_sub <- c(range(subset_df$log10_score, na.rm=TRUE)[1], 1)

make_subset_panel <- function(panel_name) {
  df <- subset_df %>% filter(panel == panel_name)
  if (nrow(df) == 0) return(ggplot() + theme_void() +
                              labs(title=paste0(panel_name, " (no data)")))
  
  whisker_stats <- df %>%
    dplyr::group_by(expr_bin) %>%
    dplyr::summarise(
      n             = dplyr::n(),
      q3_size       = quantile(log10_size,  0.75, na.rm=TRUE),
      iqr_size      = IQR(log10_size, na.rm=TRUE),
      max_size      = max(log10_size, na.rm=TRUE),
      whisker_size  = min(max_size, q3_size + 1.5*iqr_size),
      med_size_raw  = median(10^log10_size, na.rm=TRUE),
      q3_score      = quantile(log10_score, 0.75, na.rm=TRUE),
      iqr_score     = IQR(log10_score, na.rm=TRUE),
      max_score     = min(max(log10_score, na.rm=TRUE), 1),
      whisker_score = min(max_score, q3_score + 1.5*iqr_score),
      med_score_raw = median(10^log10_score, na.rm=TRUE),
      .groups       = "drop"
    ) %>%
    dplyr::mutate(
      expr_bin    = factor(expr_bin, levels=BINS),
      label_size  = paste0(round(med_size_raw), " kb (n=", n, ")"),
      label_score = as.character(round(med_score_raw, 2))
    )
  
  panel_color <- switch(panel_name,
                        "WT only"   = CONDITIONS$wt$color,
                        "pol2 only" = CONDITIONS$pol2$color,
                        "gray35"
  )
  med_x <- median(df$log10_size,  na.rm=TRUE)
  med_y <- median(df$log10_score, na.rm=TRUE)
  
  # Main scatter
  p_main <- ggplot(df, aes(x=log10_size, y=log10_score, color=expr_bin)) +
    geom_point(alpha=0.25, size=0.5) +
    geom_vline(xintercept=med_x, linetype="dashed", color="gray50", linewidth=0.4) +
    geom_hline(yintercept=med_y, linetype="dashed", color="gray50", linewidth=0.4) +
    scale_color_manual(values=BIN_COLORS_CB, name="Expression") +
    scale_x_continuous(limits=x_lims_sub, expand=c(0.02,0),
                       breaks=log10(c(10,100,1000,10000)),
                       labels=c("10","100","1k","10k"), name="Loop size (kb)") +
    scale_y_continuous(limits=y_lims_sub, expand=c(0.02,0), name=OBS_AXIS_LABEL) +
    annotate("text", x=Inf, y=Inf, label=panel_name,
             hjust=1.1, vjust=1.5, size=4.5, fontface="bold", color=panel_color) +
    annotate("text", x=-Inf, y=-Inf, label=paste0("n=", nrow(df)),
             hjust=-0.1, vjust=-0.4, size=3.2, color="gray40") +
    theme_classic(base_size=11) +
    theme(legend.position="none", plot.margin=margin(0,0,0,0))
  
  # Top margin — horizontal boxplots, x=loop size, y=expr_bin
  p_top <- ggplot(df, aes(x=log10_size, y=expr_bin, fill=expr_bin, color=expr_bin)) +
    geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
    geom_text(data=whisker_stats,
              aes(x=whisker_size, y=expr_bin, label=label_size),
              inherit.aes=FALSE, hjust=-0.1, vjust=0.5, size=2.2, color="gray30") +
    geom_vline(xintercept=med_x, linetype="dashed", color="gray50", linewidth=0.4) +
    scale_fill_manual(values=BIN_COLORS_CB, guide="none") +
    scale_color_manual(values=BIN_COLORS_CB, guide="none") +
    scale_x_continuous(limits=c(x_lims_sub[1], x_lims_sub[2]+1.5), expand=c(0.02,0)) +
    scale_y_discrete(limits=rev(BINS), labels=rev(BINS),
                     expand=expansion(add=c(0.5,0.5))) +
    theme_void() +
    theme(axis.text.y=element_text(size=7, color="gray30", hjust=1),
          plot.margin=margin(2,0,0,0))
  
  # Right margin — vertical boxplots, x=expr_bin, y=score
  # No coord_flip: y-axis naturally aligns with p_main's y-axis
  # Labels sit above each whisker via vjust
  p_right <- ggplot(df, aes(x=expr_bin, y=log10_score, fill=expr_bin, color=expr_bin)) +
    geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
    geom_text(data=whisker_stats,
              aes(x=expr_bin, y=whisker_score, label=label_score),
              inherit.aes=FALSE, hjust=0.5, vjust=-0.4,
              size=2.2, color="gray30") +
    geom_hline(yintercept=med_y, linetype="dashed", color="gray50", linewidth=0.4) +
    scale_fill_manual(values=BIN_COLORS_CB, guide="none") +
    scale_color_manual(values=BIN_COLORS_CB, guide="none") +
    scale_y_continuous(limits=c(y_lims_sub[1], y_lims_sub[2]+0.15), expand=c(0.02,0)) +
    scale_x_discrete(limits=BINS, expand=expansion(add=c(0.5,0.5))) +
    theme_void() +
    theme(plot.margin=margin(0,2,0,0))
  
  p_top + plot_spacer() + p_main + p_right +
    plot_layout(ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
}

panel_wt_only   <- make_subset_panel("WT only")
panel_shared    <- make_subset_panel("Shared")
panel_pol2_only <- make_subset_panel("pol2 only")

row_subset <- cowplot::plot_grid(
  panel_wt_only, panel_shared, panel_pol2_only,
  nrow=1, labels=c("a","b","c"), label_size=14
)

legend_df_sub <- data.frame(x=1:3, y=1, expr_bin=factor(BINS, levels=BINS))
legend_sub <- ggplot(legend_df_sub, aes(x=x, y=y, fill=expr_bin, color=expr_bin)) +
  geom_point(size=4, shape=22, stroke=1) +
  scale_fill_manual(values=BIN_COLORS_CB, name="Expression bin") +
  scale_color_manual(values=BIN_COLORS_CB, name="Expression bin") +
  theme_void() + theme(legend.position="bottom")

title_sub <- cowplot::ggdraw() +
  cowplot::draw_label(
    paste0("Shared vs exclusive loops: WT only | Shared | pol2 only  |  ", DS_GROUP),
    fontface="bold", size=13, x=0, hjust=0, y=0.75) +
  cowplot::draw_label(
    paste0(OBS_SUBTITLE,
           "  |  Shared: ", loop_overlap$n_shared,
           "  |  WT only: ", loop_overlap$n_a_only,
           "  |  pol2 only: ", loop_overlap$n_b_only,
           "  |  marginal boxplots by expression bin"),
    size=8, color="gray40", x=0, hjust=0, y=0.25)

p_shared_exclusive <- cowplot::plot_grid(
  title_sub, row_subset, cowplot::get_legend(legend_sub),
  ncol=1, rel_heights=c(0.08, 1, 0.07)
)
print(p_shared_exclusive)

# ══════════════════════════════════════════════════════════════════════════════
# 12. Shared vs exclusive CTCF loop scatter — WT only | Shared | pol2 only
#     Marginal boxplots by CTCF pair type (CTCF-CTCF, CTCF-Promoter, CTCF-Enhancer)
#     Requires loops_np_wt and loops_np_pol2 from section 10
# ══════════════════════════════════════════════════════════════════════════════
message("\n── CTCF shared / exclusive loops (WT vs pol2) ───────────────────────")

CTCF_PAIR_TYPES  <- c("CTCF \u2014 CTCF", "CTCF \u2014 Promoter", "CTCF \u2014 Enhancer")
CTCF_PAIR_COLORS <- c(
  "CTCF \u2014 CTCF"     = "#40004B",
  "CTCF \u2014 Promoter" = "#762A83",
  "CTCF \u2014 Enhancer" = "#C2A5CF"
)

label_ctcf_pairs <- function(df) {
  df %>%
    dplyr::mutate(
      ctcf_pair_type = dplyr::case_when(
        anchor1_ctcf & anchor2_ctcf                                               ~ "CTCF \u2014 CTCF",
        anchor1_ctcf & anchor2_ccre == "Promoter"                                 ~ "CTCF \u2014 Promoter",
        anchor2_ctcf & anchor1_ccre == "Promoter"                                 ~ "CTCF \u2014 Promoter",
        anchor1_ctcf & anchor2_ccre %in% c("Proximal-Enhancer","Distal-Enhancer") ~ "CTCF \u2014 Enhancer",
        anchor2_ctcf & anchor1_ccre %in% c("Proximal-Enhancer","Distal-Enhancer") ~ "CTCF \u2014 Enhancer",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(ctcf_pair_type))
}

ctcf_wt_classified   <- label_ctcf_pairs(loops_np_wt)
ctcf_pol2_classified <- label_ctcf_pairs(loops_np_pol2)

message("WT CTCF loops:   ", nrow(ctcf_wt_classified))
message("pol2 CTCF loops: ", nrow(ctcf_pol2_classified))
print(table(ctcf_wt_classified$ctcf_pair_type))
print(table(ctcf_pol2_classified$ctcf_pair_type))

ctcf_overlap <- find_shared_loops(ctcf_wt_classified, ctcf_pol2_classified, slop=5000L)
message("CTCF shared:    ", ctcf_overlap$n_shared,
        " (", round(100*ctcf_overlap$n_shared/nrow(ctcf_wt_classified),1), "% of WT)")
message("CTCF WT only:   ", ctcf_overlap$n_a_only)
message("CTCF pol2 only: ", ctcf_overlap$n_b_only)

build_ctcf_subset_scatter <- function() {
  get_score <- function(df, cond_name) {
    if (SCORE_MODE == "fold_enrichment") df$observed / df$expectedDonut
    else df$observed * DEPTH_SCALE[[cond_name]]
  }
  df_wt_only <- ctcf_wt_classified %>%
    dplyr::slice(ctcf_overlap$a_only) %>%
    dplyr::mutate(panel="WT only",
                  loop_kb=abs(((start1+end1)/2)-((start2+end2)/2))/1000,
                  score=get_score(., "wt"))
  df_pol2_only <- ctcf_pol2_classified %>%
    dplyr::slice(ctcf_overlap$b_only) %>%
    dplyr::mutate(panel="pol2 only",
                  loop_kb=abs(((start1+end1)/2)-((start2+end2)/2))/1000,
                  score=get_score(., "pol2"))
  df_shared <- ctcf_wt_classified %>%
    dplyr::slice(ctcf_overlap$shared$a_idx) %>%
    dplyr::mutate(panel="Shared",
                  loop_kb=abs(((start1+end1)/2)-((start2+end2)/2))/1000,
                  score=get_score(., "wt"))
  bind_rows(df_wt_only, df_pol2_only, df_shared) %>%
    filter(!is.na(loop_kb), !is.na(score), score > 0,
           ctcf_pair_type %in% CTCF_PAIR_TYPES) %>%
    dplyr::mutate(
      log10_size     = log10(loop_kb),
      log10_score    = log10(score),
      ctcf_pair_type = factor(ctcf_pair_type, levels=CTCF_PAIR_TYPES),
      panel          = factor(panel, levels=c("WT only","Shared","pol2 only"))
    )
}

ctcf_subset_df <- build_ctcf_subset_scatter()
message("CTCF subset scatter rows: ", nrow(ctcf_subset_df))
print(table(ctcf_subset_df$panel, ctcf_subset_df$ctcf_pair_type))

x_lims_ctcf <- c(range(ctcf_subset_df$log10_size,  na.rm=TRUE)[1], log10(10000))
y_lims_ctcf <- c(range(ctcf_subset_df$log10_score, na.rm=TRUE)[1], 1)

make_ctcf_subset_panel <- function(panel_name) {
  df <- ctcf_subset_df %>% filter(panel == panel_name)
  if (nrow(df) == 0) return(ggplot() + theme_void() +
                              labs(title=paste0(panel_name, " (no data)")))
  
  whisker_stats <- df %>%
    dplyr::group_by(ctcf_pair_type) %>%
    dplyr::summarise(
      n             = dplyr::n(),
      q3_size       = quantile(log10_size,  0.75, na.rm=TRUE),
      iqr_size      = IQR(log10_size, na.rm=TRUE),
      max_size      = max(log10_size, na.rm=TRUE),
      whisker_size  = min(max_size, q3_size + 1.5*iqr_size),
      med_size_raw  = median(10^log10_size, na.rm=TRUE),
      q3_score      = quantile(log10_score, 0.75, na.rm=TRUE),
      iqr_score     = IQR(log10_score, na.rm=TRUE),
      max_score     = min(max(log10_score, na.rm=TRUE), 1),
      whisker_score = min(max_score, q3_score + 1.5*iqr_score),
      med_score_raw = median(10^log10_score, na.rm=TRUE),
      .groups       = "drop"
    ) %>%
    dplyr::mutate(
      ctcf_pair_type = factor(ctcf_pair_type, levels=CTCF_PAIR_TYPES),
      label_size     = paste0(round(med_size_raw), " kb (n=", n, ")"),
      label_score    = as.character(round(med_score_raw, 2))
    )
  
  panel_color <- switch(panel_name,
                        "WT only"   = CONDITIONS$wt$color,
                        "pol2 only" = CONDITIONS$pol2$color,
                        "gray35"
  )
  med_x <- median(df$log10_size,  na.rm=TRUE)
  med_y <- median(df$log10_score, na.rm=TRUE)
  
  p_main <- ggplot(df, aes(x=log10_size, y=log10_score, color=ctcf_pair_type)) +
    geom_point(alpha=0.3, size=0.5) +
    geom_vline(xintercept=med_x, linetype="dashed", color="gray50", linewidth=0.4) +
    geom_hline(yintercept=med_y, linetype="dashed", color="gray50", linewidth=0.4) +
    scale_color_manual(values=CTCF_PAIR_COLORS, name="CTCF pair type") +
    scale_x_continuous(limits=x_lims_ctcf, expand=c(0.02,0),
                       breaks=log10(c(10,100,1000,10000)),
                       labels=c("10","100","1k","10k"), name="Loop size (kb)") +
    scale_y_continuous(limits=y_lims_ctcf, expand=c(0.02,0), name=OBS_AXIS_LABEL) +
    annotate("text", x=Inf, y=Inf, label=panel_name,
             hjust=1.1, vjust=1.5, size=4.5, fontface="bold", color=panel_color) +
    annotate("text", x=-Inf, y=-Inf, label=paste0("n=", nrow(df)),
             hjust=-0.1, vjust=-0.4, size=3.2, color="gray40") +
    theme_classic(base_size=11) +
    theme(legend.position="none", plot.margin=margin(0,0,0,0))
  
  # Top margin — horizontal boxplots, x=loop size, y=ctcf_pair_type
  p_top <- ggplot(df, aes(x=log10_size, y=ctcf_pair_type,
                          fill=ctcf_pair_type, color=ctcf_pair_type)) +
    geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
    geom_text(data=whisker_stats,
              aes(x=whisker_size, y=ctcf_pair_type, label=label_size),
              inherit.aes=FALSE, hjust=-0.1, vjust=0.5, size=2.2, color="gray30") +
    geom_vline(xintercept=med_x, linetype="dashed", color="gray50", linewidth=0.4) +
    scale_fill_manual(values=CTCF_PAIR_COLORS, guide="none") +
    scale_color_manual(values=CTCF_PAIR_COLORS, guide="none") +
    scale_x_continuous(limits=c(x_lims_ctcf[1], x_lims_ctcf[2]+1.5), expand=c(0.02,0)) +
    scale_y_discrete(limits=rev(CTCF_PAIR_TYPES), labels=rev(CTCF_PAIR_TYPES),
                     expand=expansion(add=c(0.5,0.5))) +
    theme_void() +
    theme(axis.text.y=element_text(size=7, color="gray30", hjust=1),
          plot.margin=margin(2,0,0,0))
  
  # Right margin — vertical boxplots, x=ctcf_pair_type, y=score
  # No coord_flip: y-axis naturally aligns with p_main's y-axis
  p_right <- ggplot(df, aes(x=ctcf_pair_type, y=log10_score,
                            fill=ctcf_pair_type, color=ctcf_pair_type)) +
    geom_boxplot(outlier.shape=NA, alpha=0.6, linewidth=0.4, width=0.6) +
    geom_text(data=whisker_stats,
              aes(x=ctcf_pair_type, y=whisker_score, label=label_score),
              inherit.aes=FALSE, hjust=0.5, vjust=-0.4,
              size=2.2, color="gray30") +
    geom_hline(yintercept=med_y, linetype="dashed", color="gray50", linewidth=0.4) +
    scale_fill_manual(values=CTCF_PAIR_COLORS, guide="none") +
    scale_color_manual(values=CTCF_PAIR_COLORS, guide="none") +
    scale_y_continuous(limits=c(y_lims_ctcf[1], y_lims_ctcf[2]+0.15), expand=c(0.02,0)) +
    scale_x_discrete(limits=CTCF_PAIR_TYPES, expand=expansion(add=c(0.5,0.5))) +
    theme_void() +
    theme(plot.margin=margin(0,2,0,0))
  
  p_top + plot_spacer() + p_main + p_right +
    plot_layout(ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
}

panel_ctcf_wt   <- make_ctcf_subset_panel("WT only")
panel_ctcf_sh   <- make_ctcf_subset_panel("Shared")
panel_ctcf_pol2 <- make_ctcf_subset_panel("pol2 only")

row_ctcf <- cowplot::plot_grid(
  panel_ctcf_wt, panel_ctcf_sh, panel_ctcf_pol2,
  nrow=1, labels=c("a","b","c"), label_size=14
)

legend_ctcf_df <- data.frame(x=1:3, y=1,
                             ctcf_pair_type=factor(CTCF_PAIR_TYPES, levels=CTCF_PAIR_TYPES))
legend_ctcf <- ggplot(legend_ctcf_df,
                      aes(x=x, y=y, fill=ctcf_pair_type, color=ctcf_pair_type)) +
  geom_point(size=4, shape=22, stroke=1) +
  scale_fill_manual(values=CTCF_PAIR_COLORS, name="CTCF pair type") +
  scale_color_manual(values=CTCF_PAIR_COLORS, name="CTCF pair type") +
  theme_void() + theme(legend.position="bottom")

title_ctcf <- cowplot::ggdraw() +
  cowplot::draw_label(
    paste0("CTCF loops — shared vs exclusive: WT only | Shared | pol2 only  |  ", DS_GROUP),
    fontface="bold", size=13, x=0, hjust=0, y=0.75) +
  cowplot::draw_label(
    paste0(OBS_SUBTITLE,
           "  |  Shared: ", ctcf_overlap$n_shared,
           "  |  WT only: ", ctcf_overlap$n_a_only,
           "  |  pol2 only: ", ctcf_overlap$n_b_only,
           "  |  narrowPeak + ATAC  |  marginals by CTCF pair type"),
    size=8, color="gray40", x=0, hjust=0, y=0.25)

p_ctcf_shared_exclusive <- cowplot::plot_grid(
  title_ctcf, row_ctcf, cowplot::get_legend(legend_ctcf),
  ncol=1, rel_heights=c(0.08, 1, 0.07)
)
print(p_ctcf_shared_exclusive)
