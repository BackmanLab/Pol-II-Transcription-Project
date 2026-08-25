##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# HiCCUPS loop cCRE analysis. Anchor classification via narrowPeak
# (ChIP-seq + ATAC): loop anchor cCRE composition, and pure
# CTCF-CTCF loop classification (no Promoter/Enhancer co-occupancy).
#
# Figure output: SI Fig. 8E (probable, non-gene anchor
# classification); SI Fig. 9E (probable, pure CTCF Gained/Stable/Lost).
##################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(cowplot)
library(GenomicRanges)
library(scales)
library(ggdist)
library(ggridges)
library(ggExtra)
library(strawr)
source("~/Transcription Paper/Hic_loop_analysis/final_scripts/EP_analysis_functions.R")

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

message("DS_GROUP: ", DS_GROUP)
message("Depth scale WT: ", DEPTH_SCALE_WT, "  Pol2: ", round(DEPTH_SCALE_DEG, 3))

DESEQ_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
EXPR_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
CHIP_DIR   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"
ATAC_DIR   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"

# ── Score mode switch ──────────────────────────────────────────────────────────
# "fold_enrichment"  — observed / expectedDonut (controls for distance decay)
# "observed_scaled"  — SCALE-normalized contact counts (field standard for APA)
SCORE_MODE <- "fold_enrichment"   # ← change this to switch

if (SCORE_MODE == "fold_enrichment") {
  SCORE_COL   <- "fold_enrichment"
  SCORE_LABEL <- "Loop fold enrichment (log10 obs/exp)"
  CAP_VAL     <- 20     # typical HiCCUPS fold enrichment range: 1.5–20
} else if (SCORE_MODE == "observed_scaled") {
  SCORE_COL   <- "observed_scaled"
  SCORE_LABEL <- "SCALE-normalized contacts (log10)"
  CAP_VAL     <- 500
} else {
  stop("Unknown SCORE_MODE: '", SCORE_MODE, "'. Must be 'fold_enrichment' or 'observed_scaled'.")
}
message("SCORE_MODE: ", SCORE_MODE, "  →  score column: ", SCORE_COL)

PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1
MIN_N          <- 10

# ══════════════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════════════
type_colors_auxin <- c(
  "Promoter — Promoter"                   = "#b2182b",
  "Promoter — Proximal-Enhancer"          = "#f4a582",
  "Promoter — Distal-Enhancer"            = "#d6604d",
  "Proximal-Enhancer — Proximal-Enhancer" = "#92c5de",
  "Proximal-Enhancer — Distal-Enhancer"   = "#74B9D3",
  "Distal-Enhancer — Distal-Enhancer"     = "#4393c3",
  "CTCF — Promoter"                       = "#762A83",
  "CTCF — Proximal-Enhancer"              = "#9970AB",
  "CTCF — Distal-Enhancer"               = "#C2A5CF",
  "CTCF — CTCF"                          = "#40004B"
)
single_anchor_colors <- c(
  "Promoter"          = "#b2182b",
  "Proximal-Enhancer" = "#f4a582",
  "Distal-Enhancer"   = "#d6604d",
  "CTCF"              = "#762A83",
  "Open-Chromatin"    = "#8073AC"
)
GROUP_COLORS <- c("P-P"="#882255", "E-P"="#44AA99", "E-E"="#4393c3")
CTCF_GROUP_COLORS <- c(
  "CTCF — CTCF"     = "#40004B",
  "CTCF — Promoter" = "#762A83",
  "CTCF — Enhancer" = "#C2A5CF"
)
COND_COLORS <- c("Wildtype"="#88CCEE", "POLR2A deg"="#CC6677")
ctcf_np_pair_colors <- c(
  "CTCF — CTCF"     = "#40004B",
  "CTCF — Promoter" = "#762A83",
  "CTCF — Enhancer" = "#C2A5CF"
)
GROUP_MAP <- c(
  "Promoter — Promoter"                   = "P-P",
  "Promoter — Distal-Enhancer"            = "E-P",
  "Promoter — Proximal-Enhancer"          = "E-P",
  "Distal-Enhancer — Distal-Enhancer"     = "E-E",
  "Proximal-Enhancer — Proximal-Enhancer" = "E-E",
  "Proximal-Enhancer — Distal-Enhancer"   = "E-E"
)
CTCF_GROUP_MAP <- c(
  "CTCF — CTCF"              = "CTCF — CTCF",
  "CTCF — Promoter"          = "CTCF — Promoter",
  "CTCF — Proximal-Enhancer" = "CTCF — Enhancer",
  "CTCF — Distal-Enhancer"   = "CTCF — Enhancer",
  "CTCF — Enhancer"          = "CTCF — Enhancer"
)
interesting_types <- c(
  "Promoter — Distal-Enhancer",
  "Promoter — Proximal-Enhancer",
  "Promoter — Promoter",
  "Distal-Enhancer — Distal-Enhancer",
  "Proximal-Enhancer — Proximal-Enhancer",
  "CTCF — Promoter",
  "CTCF — Distal-Enhancer",
  "CTCF — CTCF"
)
ctcf_np_pair_types  <- c("CTCF — CTCF","CTCF — Promoter","CTCF — Enhancer")
ep_anchor_types     <- c("Promoter","Proximal-Enhancer","Distal-Enhancer")

ridge_order_matched <- c("P-P — Wildtype","P-P — POLR2A deg",
                         "E-P — Wildtype","E-P — POLR2A deg",
                         "E-E — Wildtype","E-E — POLR2A deg")
ridge_fill_m <- c(
  "P-P — Wildtype"   = scales::alpha(GROUP_COLORS[["P-P"]], 0.85),
  "P-P — POLR2A deg" = scales::alpha(GROUP_COLORS[["P-P"]], 0.35),
  "E-P — Wildtype"   = scales::alpha(GROUP_COLORS[["E-P"]], 0.85),
  "E-P — POLR2A deg" = scales::alpha(GROUP_COLORS[["E-P"]], 0.35),
  "E-E — Wildtype"   = scales::alpha(GROUP_COLORS[["E-E"]], 0.85),
  "E-E — POLR2A deg" = scales::alpha(GROUP_COLORS[["E-E"]], 0.35)
)
ridge_color_m <- c(
  "P-P — Wildtype"   = GROUP_COLORS[["P-P"]], "P-P — POLR2A deg" = GROUP_COLORS[["P-P"]],
  "E-P — Wildtype"   = GROUP_COLORS[["E-P"]], "E-P — POLR2A deg" = GROUP_COLORS[["E-P"]],
  "E-E — Wildtype"   = GROUP_COLORS[["E-E"]], "E-E — POLR2A deg" = GROUP_COLORS[["E-E"]]
)
ctcf_ridge_order <- as.vector(rbind(paste0(names(CTCF_GROUP_COLORS)," — Wildtype"),
                                    paste0(names(CTCF_GROUP_COLORS)," — POLR2A deg")))
ctcf_ridge_fill <- unlist(lapply(names(CTCF_GROUP_COLORS), function(g) {
  setNames(c(scales::alpha(CTCF_GROUP_COLORS[g],0.85), scales::alpha(CTCF_GROUP_COLORS[g],0.35)),
           c(paste0(g," — Wildtype"), paste0(g," — POLR2A deg")))
}))
ctcf_ridge_color <- unlist(lapply(names(CTCF_GROUP_COLORS), function(g) {
  setNames(c(CTCF_GROUP_COLORS[g], CTCF_GROUP_COLORS[g]),
           c(paste0(g," — Wildtype"), paste0(g," — POLR2A deg")))
}))

# ══════════════════════════════════════════════════════════════════════════════
# Functions
# ══════════════════════════════════════════════════════════════════════════════
load_hiccups_bedpe <- function(path) {
  all_lines   <- readLines(path)
  header_line <- all_lines[grepl("^#chr1", all_lines)][1]
  col_names   <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  data_lines  <- all_lines[!grepl("^#", all_lines)]
  if (length(data_lines) == 0) stop("No data lines in ", path)
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
    observed      = as.numeric(df[["observed"]]),
    expectedDonut = as.numeric(df[["expectedDonut"]]),
    fdrDonut  = as.numeric(df[["fdrDonut"]]),
    fdrBL     = as.numeric(df[["fdrBL"]]),
    fdrH      = as.numeric(df[["fdrH"]]),
    fdrV      = as.numeric(df[["fdrV"]]),
    centroid1 = as.integer(df[["centroid1"]]),
    centroid2 = as.integer(df[["centroid2"]])
  ) %>% filter(!is.na(start1), !is.na(start2))
}

filter_loops <- function(df, min_observed = 10, max_fdr = 0.1) {
  df %>% filter(observed >= min_observed,
                fdrDonut <= max_fdr, fdrBL <= max_fdr,
                fdrH     <= max_fdr, fdrV  <= max_fdr)
}

link_loops_to_de_genes <- function(classified_loops, gene_gr,
                                   condition_label, de_filter = TRUE) {
  anchor1_gr <- GRanges(classified_loops$chr1, IRanges(classified_loops$start1, classified_loops$end1))
  anchor2_gr <- GRanges(classified_loops$chr2, IRanges(classified_loops$start2, classified_loops$end2))
  gene_gr_m  <- gene_gr
  has_chr_loop <- any(grepl("^chr", as.character(seqnames(anchor1_gr))))
  has_chr_gene <- any(grepl("^chr", as.character(seqnames(gene_gr))))
  if (has_chr_loop && !has_chr_gene) seqlevels(gene_gr_m) <- paste0("chr",seqlevels(gene_gr_m))
  else if (!has_chr_loop && has_chr_gene) seqlevels(gene_gr_m) <- gsub("^chr","",seqlevels(gene_gr_m))
  hits1 <- findOverlaps(anchor1_gr, gene_gr_m)
  hits2 <- findOverlaps(anchor2_gr, gene_gr_m)
  all_pairs <- bind_rows(
    data.frame(loop_idx=queryHits(hits1), gene_idx=subjectHits(hits1)),
    data.frame(loop_idx=queryHits(hits2), gene_idx=subjectHits(hits2))
  ) %>%
    distinct(loop_idx, gene_idx) %>%
    dplyr::mutate(
      loop_ccre_type  = classified_loops$loop_ccre_type[loop_idx],
      fold_enrichment = classified_loops$fold_enrichment[loop_idx],
      observed_scaled = classified_loops$observed_scaled[loop_idx],
      symbol          = gene_gr_m$symbol[gene_idx],
      de_status       = as.character(gene_gr_m$de_status[gene_idx]),
      condition       = condition_label
    ) %>%
    filter(!is.na(loop_ccre_type))
  if (de_filter) all_pairs <- all_pairs %>% filter(de_status %in% c("Up","Down"))
  message(condition_label, ": ", nrow(all_pairs), " loop-gene pairs",
          if (de_filter) " at DE genes" else " at expressed genes (TPM>1)")
  all_pairs
}

add_group <- function(df, cond_label, gmap = GROUP_MAP) {
  df %>%
    dplyr::mutate(
      loop_dist = abs(as.integer((start2+end2)/2) - as.integer((start1+end1)/2)),
      group     = dplyr::case_when(loop_ccre_type %in% names(gmap) ~ gmap[loop_ccre_type],
                                   TRUE ~ NA_character_),
      condition = cond_label
    ) %>%
    filter(!is.na(group), loop_dist > 0) %>%
    dplyr::mutate(group=factor(group, levels=names(GROUP_COLORS)),
                  condition=factor(condition, levels=c("Wildtype","POLR2A deg")))
}

label_ctcf_np_loops <- function(np_classified) {
  np_classified %>%
    dplyr::mutate(
      ctcf_pair_type = dplyr::case_when(
        anchor1_ctcf & anchor2_ctcf                                                ~ "CTCF — CTCF",
        anchor1_ctcf & anchor2_ccre == "Promoter"                                  ~ "CTCF — Promoter",
        anchor2_ctcf & anchor1_ccre == "Promoter"                                  ~ "CTCF — Promoter",
        anchor1_ctcf & anchor2_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")  ~ "CTCF — Enhancer",
        anchor2_ctcf & anchor1_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")  ~ "CTCF — Enhancer",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(ctcf_pair_type))
}

check_group_counts <- function(np_classified, label = "") {
  message("\n── Loop pair counts: ", label, " ──────────────────────────────────────")
  message("Total loops: ", nrow(np_classified))
  message("\nAll pair types (top 15):")
  print(sort(table(np_classified$loop_ccre_type), decreasing=TRUE)[1:15])
  grouped <- np_classified %>%
    dplyr::mutate(group=dplyr::case_when(
      loop_ccre_type %in% names(GROUP_MAP) ~ GROUP_MAP[loop_ccre_type], TRUE ~ NA_character_))
  message("\nP-P / E-P / E-E counts:"); print(table(grouped$group, useNA="ifany"))
  grp_tbl <- table(grouped$group)
  message("\nP-P / E-P / E-E as % of classified loops:")
  print(round(100*grp_tbl/sum(grp_tbl,na.rm=TRUE),1))
}

make_group_plots <- function(combined, group_map, group_colors,
                             score_col    = "fold_enrichment",
                             score_label  = SCORE_LABEL,
                             title_prefix = "",
                             ds_label     = DS_GROUP) {
  bar_df <- combined %>%
    dplyr::count(condition, group) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(prop=n/sum(n), pct=paste0(round(prop*100,1),"%")) %>%
    dplyr::ungroup()
  totals <- bar_df %>% dplyr::group_by(condition) %>%
    dplyr::summarise(total=sum(n),.groups="drop")
  p_count <- ggplot(bar_df, aes(x=condition,y=n,fill=group)) +
    geom_bar(stat="identity",width=0.65,color="white",linewidth=0.3) +
    geom_text(aes(label=n),position=position_stack(vjust=0.5),color="white",size=3,fontface="bold") +
    geom_text(data=totals,aes(x=condition,y=total,label=paste0("n=",total)),
              inherit.aes=FALSE,vjust=-0.4,size=3,color="gray30") +
    scale_fill_manual(values=group_colors,name=NULL) +
    scale_y_continuous(name="Loop count",expand=expansion(mult=c(0,0.12))) +
    labs(x=NULL) + theme_classic(base_size=11) +
    theme(axis.text.x=element_text(size=11,face="bold",
                                   color=COND_COLORS[levels(combined$condition)]),
          legend.position="none")
  p_prop <- ggplot(bar_df, aes(x=condition,y=prop,fill=group)) +
    geom_bar(stat="identity",width=0.65,color="white",linewidth=0.3) +
    geom_text(aes(label=pct),position=position_stack(vjust=0.5),color="white",size=3,fontface="bold") +
    scale_fill_manual(values=group_colors,name=NULL) +
    scale_y_continuous(name="Proportion",labels=scales::percent,expand=expansion(mult=c(0,0.05))) +
    labs(x=NULL) + theme_classic(base_size=11) +
    theme(axis.text.x=element_text(size=11,face="bold",
                                   color=COND_COLORS[levels(combined$condition)]),
          legend.position="right",legend.key.size=unit(0.4,"cm"),legend.text=element_text(size=9))
  p_bar <- (p_count|p_prop) +
    plot_layout(guides="collect") &
    theme(legend.position="right") &
    plot_annotation(
      title    = paste0(title_prefix," loop counts — Wildtype vs POLR2A deg"),
      subtitle = paste0("HiCCUPS loops  |  WT × WT peaks  |  pol2 × KO peaks  |  ", ds_label),
      theme=theme(plot.title=element_text(face="bold",size=13),
                  plot.subtitle=element_text(size=8,color="gray40"))
    )
  group_levels <- names(group_colors)
  ridge_order  <- as.vector(rbind(paste0(group_levels," — Wildtype"),
                                  paste0(group_levels," — POLR2A deg")))
  n_labels <- combined %>% dplyr::count(group,condition) %>%
    dplyr::mutate(ridge_label=paste0(group," — ",condition),
                  axis_label=paste0(group," — ",condition," (n=",n,")"))
  label_lookup <- setNames(n_labels$axis_label, n_labels$ridge_label)
  ridge_df <- combined %>%
    dplyr::mutate(ridge_label=factor(paste0(group," — ",condition),levels=rev(ridge_order))) %>%
    filter(!is.na(ridge_label), !is.na(.data[[score_col]]), .data[[score_col]] > 0)
  median_df <- ridge_df %>%
    dplyr::group_by(ridge_label) %>%
    dplyr::summarise(med_log=median(log10(.data[[score_col]]),na.rm=TRUE),.groups="drop") %>%
    dplyr::mutate(y_num=as.numeric(ridge_label),med_label=sprintf("%.2f", 10^med_log - 1))
  ridge_fill <- unlist(lapply(group_levels, function(g) {
    setNames(c(scales::alpha(group_colors[[g]],0.85),scales::alpha(group_colors[[g]],0.35)),
             c(paste0(g," — Wildtype"),paste0(g," — POLR2A deg")))
  }))
  ridge_color <- unlist(lapply(group_levels, function(g) {
    setNames(c(group_colors[[g]],group_colors[[g]]),
             c(paste0(g," — Wildtype"),paste0(g," — POLR2A deg")))
  }))
  p_ridge <- ggplot(ridge_df,aes(x=log10(.data[[score_col]]),y=ridge_label,
                                 fill=ridge_label,color=ridge_label)) +
    ggridges::geom_density_ridges(scale=0.9,alpha=0.85,linewidth=0.6,quantile_lines=TRUE,quantiles=2) +
    geom_text(data=median_df,aes(x=med_log,y=y_num+0.15,label=med_label),
              color="white",size=2.8,fontface="bold",inherit.aes=FALSE) +
    scale_fill_manual(values=ridge_fill,guide="none") +
    scale_color_manual(values=ridge_color,guide="none") +
    scale_y_discrete(labels=function(x) label_lookup[x]) +
    scale_x_continuous(name=score_label) +
    labs(title=paste0(title_prefix," loop strength distributions"),
         subtitle="dark = Wildtype  |  light = POLR2A deg  |  line = median",y=NULL) +
    theme_ridges(font_size=11) +
    theme(plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=8,color="gray40"),axis.text.y=element_text(size=10))
  make_scatter <- function(df, cond_label) {
    ggplot(df,aes(x=log10(loop_dist),y=log10(.data[[score_col]]),color=group)) +
      geom_point(alpha=0.25,size=0.6) +
      scale_color_manual(values=group_colors,name=NULL) +
      scale_x_continuous(name="Loop distance (log10 bp)",limits=c(4.5,8)) +
      scale_y_continuous(name=score_label) +
      labs(title=cond_label) + theme_classic(base_size=11) +
      theme(plot.title=element_text(face="bold",size=11,color=COND_COLORS[cond_label]),
            legend.position="bottom")
  }
  s_wt   <- make_scatter(combined %>% filter(condition=="Wildtype"),   "Wildtype")
  s_pol2 <- make_scatter(combined %>% filter(condition=="POLR2A deg"), "POLR2A deg")
  m_wt   <- ggExtra::ggMarginal(s_wt,   type="boxplot",groupColour=TRUE,groupFill=TRUE,size=6,alpha=0.5,outlier.shape=NA)
  m_pol2 <- ggExtra::ggMarginal(s_pol2, type="boxplot",groupColour=TRUE,groupFill=TRUE,size=6,alpha=0.5,outlier.shape=NA)
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(paste0(title_prefix," — loop distance vs strength"),
                        fontface="bold",size=13,x=0.02,hjust=0,y=0.75) +
    cowplot::draw_label(paste0("HiCCUPS loops  |  WT × WT peaks  |  pol2 × KO peaks  |  ", ds_label),
                        size=8,color="gray40",x=0.02,hjust=0,y=0.25)
  p_scatter <- cowplot::plot_grid(title_grob,
                                  cowplot::plot_grid(m_wt,m_pol2,ncol=2,labels=c("a","b")),
                                  ncol=1,rel_heights=c(0.08,1))
  list(bar=p_bar, ridge=p_ridge, scatter=p_scatter)
}

plot_unmatched_ccre <- function(wt_classified, deg_classified,
                                title, cap_val = CAP_VAL, min_n = MIN_N) {
  ep_types <- interesting_types[!grepl("CTCF", interesting_types)]
  plot_df <- bind_rows(
    wt_classified  %>% dplyr::mutate(condition="Wildtype"),
    deg_classified %>% dplyr::mutate(condition="POLR2A deg")
  ) %>%
    filter(loop_ccre_type %in% ep_types) %>%
    dplyr::mutate(condition=factor(condition,levels=c("Wildtype","POLR2A deg")),
                  loop_ccre_type=factor(loop_ccre_type,levels=ep_types),
                  log_score=pmin(log10(.data[[SCORE_COL]]), log10(CAP_VAL))) %>%
    dplyr::group_by(loop_ccre_type,condition) %>%
    dplyr::filter(dplyr::n() >= min_n) %>%
    dplyr::ungroup()
  if (nrow(plot_df)==0) { message("No data to plot"); return(NULL) }
  types_present <- levels(droplevels(plot_df$loop_ccre_type))
  x_map <- bind_rows(lapply(seq_along(types_present), function(i) {
    base <- (i-1)*1.6
    data.frame(loop_ccre_type=types_present[i],
               condition=c("Wildtype","POLR2A deg"),x_pos=c(base+1.0,base+1.45))
  }))
  plot_df <- plot_df %>% left_join(x_map, by=c("loop_ccre_type","condition"))
  x_labels <- plot_df %>%
    dplyr::group_by(x_pos,condition,loop_ccre_type) %>%
    dplyr::summarise(n=dplyr::n(),.groups="drop") %>%
    dplyr::mutate(label=paste0(condition,"\n(n=",n,")")) %>%
    dplyr::arrange(x_pos)
  band_df <- bind_rows(lapply(seq_along(types_present), function(i) {
    base <- (i-1)*1.6
    data.frame(loop_ccre_type=types_present[i],xmin=base+0.75,xmax=base+1.70,mid=base+1.225)
  }))
  get_pval <- function(ltype) {
    wt  <- plot_df %>% filter(loop_ccre_type==ltype,condition=="Wildtype") %>% pull(log_score)
    deg <- plot_df %>% filter(loop_ccre_type==ltype,condition=="POLR2A deg") %>% pull(log_score)
    if (length(wt)<3||length(deg)<3) return("n/a")
    p <- wilcox.test(wt,deg,paired=FALSE)$p.value
    dplyr::case_when(p<0.0001~"****",p<0.001~"***",p<0.01~"**",p<0.05~"*",TRUE~"ns")
  }
  y_max <- max(plot_df$log_score,na.rm=TRUE)
  p <- ggplot(plot_df,aes(x=x_pos,y=log_score,fill=condition,color=condition,group=x_pos))
  for (i in seq_len(nrow(band_df))) {
    ltype <- band_df$loop_ccre_type[i]
    p <- p +
      annotate("rect",xmin=band_df$xmin[i],xmax=band_df$xmax[i],ymin=-Inf,ymax=Inf,
               fill=type_colors_auxin[ltype],alpha=0.12) +
      annotate("text",x=band_df$mid[i],y=y_max*1.08,label=ltype,
               size=2.8,color=type_colors_auxin[ltype],fontface="bold") +
      annotate("text",x=band_df$mid[i],y=y_max*1.03,label=get_pval(ltype),size=4,color="gray30")
  }
  p +
    geom_boxplot(width=0.1,outlier.shape=NA,color="gray20",alpha=0.9,linewidth=0.5) +
    ggdist::stat_halfeye(aes(thickness=after_stat(pdf*n)),side="right",scale=0.18,
                         .width=0,point_colour=NA,alpha=0.7,normalize="groups",
                         position=position_nudge(x=0.12)) +
    scale_fill_manual(values=COND_COLORS,name=NULL) +
    scale_color_manual(values=c("Wildtype"="#4A90C4","POLR2A deg"="#AA4455"),name=NULL,guide="none") +
    scale_x_continuous(breaks=x_labels$x_pos,labels=x_labels$label,
                       limits=c(0.6,max(x_labels$x_pos)+0.6),expand=c(0,0)) +
    scale_y_continuous(name=SCORE_LABEL,
                       limits=c(NA,y_max*1.15),expand=expansion(mult=c(0.02,0))) +
    labs(title=title,x=NULL) + theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold",size=13),
          axis.text.x=element_text(size=7,angle=45,hjust=1,vjust=1),
          legend.position="bottom",plot.margin=margin(8,8,40,8))
}

plot_ccre_stacked <- function(all_pairs, de_direction=NULL, title,
                              min_n=5, min_label_prop=0.02) {
  plot_df <- all_pairs %>%
    { if (!is.null(de_direction)) filter(.,de_status==de_direction) else . } %>%
    filter(loop_ccre_type %in% interesting_types) %>%
    distinct(loop_idx,loop_ccre_type,condition) %>%
    dplyr::count(loop_ccre_type,condition,name="n_loops") %>%
    tidyr::complete(loop_ccre_type=interesting_types,
                    condition=levels(all_pairs$condition),fill=list(n_loops=0)) %>%
    dplyr::group_by(loop_ccre_type) %>% dplyr::filter(max(n_loops)>=min_n) %>% dplyr::ungroup() %>%
    dplyr::mutate(loop_ccre_type=factor(loop_ccre_type,levels=interesting_types),
                  condition=factor(condition,levels=c("Wildtype","POLR2A deg"))) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(prop=n_loops/sum(n_loops),
                  pct=ifelse(prop>=min_label_prop,paste0(round(prop*100,1),"%"),"")) %>%
    dplyr::ungroup()
  if (nrow(plot_df)==0) return(ggplot()+theme_void()+labs(title=title))
  totals <- plot_df %>% dplyr::group_by(condition) %>% dplyr::summarise(total=sum(n_loops),.groups="drop")
  p_count <- ggplot(plot_df,aes(x=condition,y=n_loops,fill=loop_ccre_type)) +
    geom_col(width=0.65,color="white",linewidth=0.3) +
    geom_text(aes(label=ifelse(n_loops>=min_n,n_loops,"")),
              position=position_stack(vjust=0.5),color="white",size=3,fontface="bold") +
    geom_text(data=totals,aes(x=condition,y=total,label=paste0("n=",total)),
              inherit.aes=FALSE,vjust=-0.4,size=3,color="gray30") +
    scale_fill_manual(values=type_colors_auxin,name="Loop type",drop=FALSE) +
    scale_y_continuous(name="Loop count",expand=expansion(mult=c(0,0.1))) +
    labs(x=NULL) + theme_classic(base_size=11) +
    theme(axis.text.x=element_text(size=11,face="bold",
                                   color=COND_COLORS[c("Wildtype","POLR2A deg")]),
          legend.position="none",plot.margin=margin(4,4,4,4))
  p_prop <- ggplot(plot_df,aes(x=condition,y=prop,fill=loop_ccre_type)) +
    geom_col(width=0.65,color="white",linewidth=0.3) +
    geom_text(aes(label=pct),position=position_stack(vjust=0.5),color="white",size=3,fontface="bold") +
    scale_fill_manual(values=type_colors_auxin,name="Loop type",drop=FALSE) +
    scale_y_continuous(name="Proportion",labels=scales::percent,expand=expansion(mult=c(0,0.05))) +
    labs(x=NULL) + theme_classic(base_size=11) +
    theme(axis.text.x=element_text(size=11,face="bold",
                                   color=COND_COLORS[c("Wildtype","POLR2A deg")]),
          legend.position="right",legend.key.size=unit(0.4,"cm"),
          legend.text=element_text(size=9),legend.title=element_text(size=10,face="bold"),
          plot.margin=margin(4,4,4,4))
  (p_count|p_prop) +
    plot_layout(guides="collect") &
    theme(legend.position="right") &
    plot_annotation(
      title=title,
      subtitle=paste0(if (!is.null(de_direction)) paste0(de_direction,"regulated genes  |  ")
                      else "All expressed genes (TPM>1)  |  ",
                      "Left=count  |  Right=proportion  |  narrowPeak + ATAC  |  ", DS_GROUP),
      theme=theme(plot.title=element_text(face="bold",size=13),
                  plot.subtitle=element_text(size=8,color="gray40"))
    )
}

make_ccre_count_plot_paired <- function(all_pairs, de_direction, title, min_n=5) {
  cond_levels <- levels(all_pairs$condition)
  plot_df <- all_pairs %>%
    { if (!is.null(de_direction)) filter(.,de_status==de_direction) else . } %>%
    filter(loop_ccre_type %in% interesting_types) %>%
    distinct(loop_idx,loop_ccre_type,condition) %>%
    dplyr::count(loop_ccre_type,condition,name="n_loops") %>%
    tidyr::complete(loop_ccre_type=interesting_types,condition=cond_levels,fill=list(n_loops=0)) %>%
    dplyr::group_by(loop_ccre_type) %>% dplyr::filter(max(n_loops)>=min_n) %>% dplyr::ungroup() %>%
    dplyr::mutate(loop_ccre_type=factor(loop_ccre_type,levels=interesting_types),
                  condition=factor(condition, levels=cond_levels),
                  fill_key=factor(paste0(condition,"__",loop_ccre_type),
                                  levels=c(paste0(cond_levels[1],"__",interesting_types),
                                           paste0(cond_levels[2],"__",interesting_types))))
  if (nrow(plot_df)==0) return(ggplot()+theme_void()+labs(title=title))
  type_order <- plot_df %>% dplyr::group_by(loop_ccre_type) %>%
    dplyr::summarise(total=sum(n_loops),.groups="drop") %>% arrange(desc(total)) %>%
    dplyr::mutate(x_label=gsub(" — ","\n",as.character(loop_ccre_type)))
  plot_df <- plot_df %>%
    dplyr::mutate(x_label=factor(gsub(" — ","\n",as.character(loop_ccre_type)),levels=type_order$x_label))
  fill_values <- c(
    setNames(rep("gray70",length(interesting_types)),paste0(cond_levels[1],"__",interesting_types)),
    setNames(type_colors_auxin,paste0(cond_levels[2],"__",names(type_colors_auxin)))
  )
  ggplot(plot_df,aes(x=x_label,y=n_loops,fill=fill_key,group=condition)) +
    geom_col(position=position_dodge(width=0.7),width=0.65,color="white",linewidth=0.3,alpha=0.92) +
    geom_text(aes(label=ifelse(n_loops>0,n_loops,"")),
              position=position_dodge(width=0.7),vjust=-0.4,size=2.8,color="gray30") +
    scale_fill_manual(values=fill_values,
                      breaks=c(paste0(cond_levels[1],"__Promoter — Distal-Enhancer"),
                               paste0(cond_levels[2],"__Promoter — Distal-Enhancer")),
                      labels=c("WT (grey)","POLR2A deg (colored)"),name=NULL) +
    scale_y_continuous(name="Number of loops",expand=expansion(mult=c(0,0.12))) +
    labs(title=title,
         subtitle=paste0("Left=Wildtype (grey)  |  Right=POLR2A deg (colored)  |  ",
                         "narrowPeak + ATAC  |  padj<",PADJ_THRESHOLD,"  |lfc|>",LFC_THRESHOLD),
         x=NULL) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=9,color="gray40"),
          axis.text.x=element_text(size=9),legend.position="bottom",plot.margin=margin(8,8,5,8))
}

link_loops_singleton_de <- function(classified_loops, gene_gr, condition_label) {
  anchor1_gr <- GRanges(classified_loops$chr1, IRanges(classified_loops$start1, classified_loops$end1))
  anchor2_gr <- GRanges(classified_loops$chr2, IRanges(classified_loops$start2, classified_loops$end2))
  gene_gr_m  <- gene_gr
  has_chr_loop <- any(grepl("^chr",as.character(seqnames(anchor1_gr))))
  has_chr_gene <- any(grepl("^chr",as.character(seqnames(gene_gr))))
  if (has_chr_loop && !has_chr_gene) seqlevels(gene_gr_m) <- paste0("chr",seqlevels(gene_gr_m))
  else if (!has_chr_loop && has_chr_gene) seqlevels(gene_gr_m) <- gsub("^chr","",seqlevels(gene_gr_m))
  hits1 <- findOverlaps(anchor1_gr,gene_gr_m)
  hits2 <- findOverlaps(anchor2_gr,gene_gr_m)
  df1 <- data.frame(loop_idx=queryHits(hits1),gene_idx=subjectHits(hits1),
                    other_anchor_ccre=classified_loops$anchor2_ccre[queryHits(hits1)],
                    stringsAsFactors=FALSE)
  df2 <- data.frame(loop_idx=queryHits(hits2),gene_idx=subjectHits(hits2),
                    other_anchor_ccre=classified_loops$anchor1_ccre[queryHits(hits2)],
                    stringsAsFactors=FALSE)
  result <- bind_rows(df1,df2) %>%
    dplyr::mutate(de_status=gene_gr_m$de_status[gene_idx],
                  symbol=gene_gr_m$symbol[gene_idx],
                  fold_enrichment=classified_loops$fold_enrichment[loop_idx],
                  observed_scaled=classified_loops$observed_scaled[loop_idx],
                  loop_ccre_type=classified_loops$loop_ccre_type[loop_idx],
                  condition=condition_label) %>%
    filter(de_status %in% c("Up","Down"),
           !is.na(other_anchor_ccre), other_anchor_ccre != "Unannotated") %>%
    dplyr::group_by(loop_idx,gene_idx) %>% dplyr::filter(dplyr::n()==1) %>% dplyr::ungroup()
  message(condition_label,": ",nrow(result)," singleton loop-gene pairs at DE genes")
  message("  Up: ",sum(result$de_status=="Up"),"  Down: ",sum(result$de_status=="Down"))
  print(table(result$other_anchor_ccre))
  result
}

make_singleton_count_plot <- function(all_pairs, de_direction, ccre_types, title, min_n=5) {
  cond_levels <- levels(all_pairs$condition)
  plot_df <- all_pairs %>%
    filter(de_status==de_direction,other_anchor_ccre %in% ccre_types) %>%
    distinct(loop_idx,gene_idx,other_anchor_ccre,condition) %>%
    dplyr::count(other_anchor_ccre,condition,name="n_loops") %>%
    tidyr::complete(other_anchor_ccre=ccre_types,condition=cond_levels,fill=list(n_loops=0)) %>%
    dplyr::group_by(other_anchor_ccre) %>% dplyr::filter(max(n_loops)>=min_n) %>% dplyr::ungroup() %>%
    dplyr::mutate(other_anchor_ccre=factor(other_anchor_ccre,levels=ccre_types),
                  condition=factor(condition, levels=cond_levels),
                  fill_key=factor(paste0(condition,"__",other_anchor_ccre),
                                  levels=c(paste0(cond_levels[1],"__",ccre_types),
                                           paste0(cond_levels[2],"__",ccre_types))))
  if (nrow(plot_df)==0) return(ggplot()+theme_void()+labs(title=title))
  type_order <- plot_df %>% dplyr::group_by(other_anchor_ccre) %>%
    dplyr::summarise(total=sum(n_loops),.groups="drop") %>% arrange(desc(total))
  plot_df <- plot_df %>%
    dplyr::mutate(other_anchor_ccre=factor(other_anchor_ccre,levels=type_order$other_anchor_ccre))
  fill_values <- c(
    setNames(rep("gray70",length(ccre_types)),paste0(cond_levels[1],"__",ccre_types)),
    setNames(single_anchor_colors[ccre_types],paste0(cond_levels[2],"__",ccre_types))
  )
  ggplot(plot_df,aes(x=other_anchor_ccre,y=n_loops,fill=fill_key,group=condition)) +
    geom_col(position=position_dodge(width=0.7),width=0.65,color="white",linewidth=0.3,alpha=0.92) +
    geom_text(aes(label=ifelse(n_loops>0,n_loops,"")),
              position=position_dodge(width=0.7),vjust=-0.6,size=2.8,color="gray30") +
    scale_fill_manual(values=fill_values,
                      breaks=c(paste0(cond_levels[1],"__",ccre_types[1]),
                               paste0(cond_levels[2],"__",ccre_types[1])),
                      labels=c("WT (grey)","POLR2A deg (colored)"),name=NULL) +
    scale_y_continuous(name="Number of loops",expand=expansion(mult=c(0,0.20))) +
    labs(title=title,
         subtitle=paste0("Left=Wildtype (grey)  |  Right=POLR2A deg (colored)  |  ",
                         "One anchor in DE gene body  |  narrowPeak + ATAC  |  ",
                         "padj<",PADJ_THRESHOLD,"  |lfc|>",LFC_THRESHOLD),x=NULL) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=9,color="gray40"),
          axis.text.x=element_text(size=10),legend.position="bottom",plot.margin=margin(8,8,5,8))
}

collapse_ctcf_types <- function(df) {
  df %>%
    dplyr::mutate(
      loop_ccre_type = dplyr::case_when(
        loop_ccre_type %in% c("CTCF — Unannotated","CTCF — Other","CTCF — Open-Chromatin") ~ "CTCF — Other",
        TRUE ~ loop_ccre_type
      )
    )
}

# CTCF-specific singleton plot — one anchor in DE gene, other anchor is CTCF
# Shows breakdown by CTCF loop sub-type (CTCF-CTCF, CTCF-Promoter, CTCF-Enhancer, CTCF-Other)
make_ctcf_singleton_plot <- function(singleton_df, de_direction, title, min_n=5) {
  ctcf_types_here <- c("CTCF — CTCF","CTCF — Promoter",
                       "CTCF — Proximal-Enhancer","CTCF — Distal-Enhancer","CTCF — Other")
  ctcf_colors_here <- c(
    "CTCF — CTCF"              = "#40004B",
    "CTCF — Promoter"          = "#762A83",
    "CTCF — Proximal-Enhancer" = "#9970AB",
    "CTCF — Distal-Enhancer"   = "#C2A5CF",
    "CTCF — Other"             = "#E7D4E8"
  )
  cond_levels <- levels(singleton_df$condition)
  plot_df <- singleton_df %>%
    filter(de_status==de_direction, other_anchor_ccre=="CTCF",
           loop_ccre_type %in% ctcf_types_here) %>%
    distinct(loop_idx, gene_idx, loop_ccre_type, condition) %>%
    dplyr::count(loop_ccre_type, condition, name="n_loops") %>%
    tidyr::complete(loop_ccre_type=ctcf_types_here, condition=cond_levels,
                    fill=list(n_loops=0)) %>%
    dplyr::group_by(loop_ccre_type) %>%
    dplyr::filter(max(n_loops) >= min_n) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(loop_ccre_type=factor(loop_ccre_type, levels=ctcf_types_here),
                  condition=factor(condition, levels=cond_levels),
                  fill_key=factor(paste0(condition,"__",loop_ccre_type),
                                  levels=c(paste0(cond_levels[1],"__",ctcf_types_here),
                                           paste0(cond_levels[2],"__",ctcf_types_here))))
  if (nrow(plot_df)==0) return(ggplot()+theme_void()+labs(title=title))
  fill_values <- c(
    setNames(rep("gray70",length(ctcf_types_here)), paste0(cond_levels[1],"__",ctcf_types_here)),
    setNames(ctcf_colors_here[ctcf_types_here],     paste0(cond_levels[2],"__",ctcf_types_here))
  )
  ggplot(plot_df, aes(x=loop_ccre_type, y=n_loops, fill=fill_key, group=condition)) +
    geom_col(position=position_dodge(width=0.7), width=0.65, color="white", linewidth=0.3, alpha=0.92) +
    geom_text(aes(label=ifelse(n_loops>0,n_loops,"")),
              position=position_dodge(width=0.7), vjust=-0.6, size=2.8, color="gray30") +
    scale_fill_manual(values=fill_values,
                      breaks=c(paste0(cond_levels[1],"__",ctcf_types_here[1]),
                               paste0(cond_levels[2],"__",ctcf_types_here[1])),
                      labels=c("WT (grey)","POLR2A deg (colored)"), name=NULL) +
    scale_x_discrete(labels=function(x) gsub(" — ","\n",x)) +
    scale_y_continuous(name="Number of loops", expand=expansion(mult=c(0,0.20))) +
    labs(title=title,
         subtitle=paste0("Left=Wildtype (grey)  |  Right=POLR2A deg (colored)  |  ",
                         "One anchor in DE gene body  |  CTCF at other anchor  |  ",
                         "narrowPeak + ATAC  |  padj<",PADJ_THRESHOLD,"  |lfc|>",LFC_THRESHOLD),
         x=NULL) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=9,color="gray40"),
          axis.text.x=element_text(size=9), legend.position="bottom",
          plot.margin=margin(8,8,5,8))
}

# Extract observed contacts at loop anchor positions from a .hic file via strawr
# Used to compare WT loop positions scored in the pol2 .hic (matched ridgeline)
extract_loop_contacts <- function(loops_df, hic_file, norm="SCALE", resolution=5000) {
  snap <- function(x) as.integer(floor(x / resolution) * resolution)
  results <- rep(0, nrow(loops_df))
  c1_pos     <- snap(as.integer((loops_df$start1 + loops_df$end1) / 2))
  c2_pos     <- snap(as.integer((loops_df$start2 + loops_df$end2) / 2))
  chr1_clean <- gsub("^chr", "", loops_df$chr1)
  chr2_clean <- gsub("^chr", "", loops_df$chr2)
  chr_pairs  <- unique(paste(chr1_clean, chr2_clean, sep="|"))
  message("Extracting contacts from: ", basename(hic_file))
  message("  Normalization: ", norm, "  Resolution: ", resolution/1000, "kb  Loops: ", nrow(loops_df))
  for (cp in chr_pairs) {
    parts  <- strsplit(cp, "\\|")[[1]]
    c1_chr <- parts[1]; c2_chr <- parts[2]
    idx    <- which(chr1_clean == c1_chr & chr2_clean == c2_chr)
    # Query only the bounding region of loops on this chr — avoids loading full matrix
    c1_min <- min(c1_pos[idx]); c1_max <- max(c1_pos[idx]) + resolution
    c2_min <- min(c2_pos[idx]); c2_max <- max(c2_pos[idx]) + resolution
    loc1 <- paste0(c1_chr, ":", c1_min, ":", c1_max)
    loc2 <- paste0(c2_chr, ":", c2_min, ":", c2_max)
    tryCatch({
      mat <- strawr::straw(norm, hic_file, loc1, loc2, "BP", resolution)
      if (nrow(mat) == 0) { results[idx] <- 0; next }
      lookup <- setNames(mat$counts, paste(mat$x, mat$y, sep="|"))
      keys   <- paste(pmin(c1_pos[idx], c2_pos[idx]),
                      pmax(c1_pos[idx], c2_pos[idx]), sep="|")
      vals   <- lookup[keys]
      results[idx] <- ifelse(is.na(vals), 0, vals)
    }, error=function(e) {
      message("  Warning: ", c1_chr, " — ", e$message)
    })
  }
  n_hit <- sum(results > 0)
  message("  Signal: ", n_hit, " / ", nrow(loops_df), " loops  (",
          round(100*n_hit/nrow(loops_df),1), "%)")
  results
}

# ══════════════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Loading expression data ───────────────────────────────────────────")
expr <- read_tsv(EXPR_FILE) %>% dplyr::select(gene_id, TPM)
deseq_res <- read_csv(DESEQ_FILE) %>%
  left_join(expr, by="gene_id") %>%
  dplyr::mutate(
    gene_id   = sub("\\..*","",gene_id),
    de_status = dplyr::case_when(
      is.na(padj)                                               ~ "Filtered",
      padj >= PADJ_THRESHOLD                                    ~ "Not DE",
      padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD  ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD  ~ "Down",
      TRUE                                                      ~ "Not DE"
    ),
    de_status = factor(de_status,levels=c("Up","Not DE","Down","Filtered")),
    symbol    = dplyr::case_when(!is.na(symbol)&symbol!="" ~ symbol, TRUE ~ gene_id)
  )
message("DE breakdown:"); print(table(deseq_res$de_status))

gene_coords <- deseq_res %>%
  filter(!is.na(gene_seq_start),!is.na(gene_seq_end),!is.na(seq_name)) %>%
  distinct(gene_id,.keep_all=TRUE)
gene_gr <- GRanges(
  seqnames=as.character(gene_coords$seq_name),
  ranges=IRanges(gene_coords$gene_seq_start,gene_coords$gene_seq_end),
  symbol=gene_coords$symbol, gene_id=gene_coords$gene_id,
  de_status=as.character(gene_coords$de_status), TPM=gene_coords$TPM
)
gene_gr_expr <- gene_gr[!is.na(gene_gr$TPM) & gene_gr$TPM > 1]
message("All genes: ",length(gene_gr),"  Expressed (TPM>1): ",length(gene_gr_expr))

message("\n── Loading HiCCUPS loops (", DS_GROUP, ") ───────────────────────────")
loops_wt  <- filter_loops(load_hiccups_bedpe(hic_wt_loops),  min_observed=MIN_OBS_WT) %>%
  dplyr::mutate(loop_idx=row_number(),
                loop_dist=abs(as.integer((start2+end2)/2)-as.integer((start1+end1)/2)),
                observed_scaled  = observed * DEPTH_SCALE_WT,
                fold_enrichment  = observed / expectedDonut)
loops_pol2 <- filter_loops(load_hiccups_bedpe(hic_deg_loops), min_observed=MIN_OBS_DEG) %>%
  dplyr::mutate(loop_idx=row_number(),
                loop_dist=abs(as.integer((start2+end2)/2)-as.integer((start1+end1)/2)),
                observed_scaled  = observed * DEPTH_SCALE_DEG,
                fold_enrichment  = observed / expectedDonut)
message("WT loops:       ", nrow(loops_wt))
message("POLR2A deg loops: ", nrow(loops_pol2))

message("\n── Loading narrowPeak files ──────────────────────────────────────────")
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

message("\n── Loading TxDb ──────────────────────────────────────────────────────")
# Build once, save, reload on future runs (build takes ~5–10 min):
# library(GenomicFeatures)
# txdb <- makeTxDbFromGFF(
#   "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf",
#   format="gtf")
# saveDb(txdb, "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/GRCh38.115.sqlite")
library(GenomicFeatures)
txdb <- loadDb("/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/GRCh38.115.sqlite")
message("TxDb loaded: ", basename("/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/GRCh38.115.sqlite"))

# ══════════════════════════════════════════════════════════════════════════════
# Classify loop anchors
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Classifying anchors — narrowPeak + ATAC + TSS geography ──────────")
wt_np_class <- classify_anchors_narrowpeak(
  loops_df   = loops_wt,
  atac_gr    = atac_wt,
  h3k4me3_gr = h3k4me3_wt,
  h3k27ac_gr = h3k27ac_wt,
  txdb       = txdb,
  h3k4me1_gr = h3k4me1_wt,
  ctcf_gr    = ctcf_wt,
  tss_region = c(-3000, 3000),
  ccre_label = "WT_narrowpeak"
)
pol2_np_class <- classify_anchors_narrowpeak(
  loops_df   = loops_pol2,
  atac_gr    = atac_ko,
  h3k4me3_gr = h3k4me3_ko,
  h3k27ac_gr = h3k27ac_ko,
  txdb       = txdb,
  h3k4me1_gr = h3k4me1_ko,
  ctcf_gr    = ctcf_ko,
  tss_region = c(-3000, 3000),
  ccre_label = "KO_narrowpeak"
)
check_group_counts(wt_np_class,   "WT narrowPeak")
check_group_counts(pol2_np_class, "pol2 KO narrowPeak")

wt_np_grouped   <- add_group(wt_np_class,   "Wildtype")
pol2_np_grouped <- add_group(pol2_np_class, "POLR2A deg")
message("WT EP/EE/PP loops:   ", nrow(wt_np_grouped))
message("pol2 EP/EE/PP loops: ", nrow(pol2_np_grouped))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 1: Global unmatched comparison
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 1: Global unmatched comparison ───────────────────────────")
print(plot_unmatched_ccre(wt_np_class, pol2_np_class,
                          paste0("Hi-C loop strength: WT vs POLR2A deg — narrowPeak + ATAC  |  ", DS_GROUP)))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 2: DE gene loop counts
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 2: DE gene loop counts ──────────────────────────────────")
pairs_wt_np   <- link_loops_to_de_genes(wt_np_class,   gene_gr, "Wildtype")
pairs_pol2_np <- link_loops_to_de_genes(pol2_np_class, gene_gr, "POLR2A deg")
all_de_pairs_np <- bind_rows(pairs_wt_np, pairs_pol2_np) %>%
  dplyr::mutate(condition=factor(condition,levels=c("Wildtype","POLR2A deg")))

(make_ccre_count_plot_paired(all_de_pairs_np,"Up",
                             "cCRE loop types at upregulated genes — narrowPeak + ATAC") /
    make_ccre_count_plot_paired(all_de_pairs_np,"Down",
                                "cCRE loop types at downregulated genes — narrowPeak + ATAC")) +
  plot_layout(heights=c(1,1)) &
  plot_annotation(title=paste0("HiCCUPS loop counts at DE genes — ", DS_GROUP),
                  theme=theme(plot.title=element_text(face="bold",size=14)))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 2b: Singleton anchor analysis at DE genes
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 2b: Singleton DE gene analysis ───────────────────────────")
singleton_wt_np   <- link_loops_singleton_de(wt_np_class,   gene_gr, "Wildtype")
singleton_pol2_np <- link_loops_singleton_de(pol2_np_class, gene_gr, "POLR2A deg")
all_singleton_np  <- bind_rows(singleton_wt_np, singleton_pol2_np) %>%
  dplyr::mutate(condition=factor(condition,levels=c("Wildtype","POLR2A deg")))

(make_singleton_count_plot(all_singleton_np,"Up",  ep_anchor_types,
                           "Non-gene anchor type at upregulated genes — narrowPeak + ATAC") /
    make_singleton_count_plot(all_singleton_np,"Down",ep_anchor_types,
                              "Non-gene anchor type at downregulated genes — narrowPeak + ATAC")) +
  plot_annotation(title=paste0("Singleton: narrowPeak + ATAC — WT vs POLR2A deg  |  ", DS_GROUP),
                  theme=theme(plot.title=element_text(face="bold",size=14)))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 2c: Global expressed genes (TPM > 1)
# ══════════════════════════════════════════════════════════════════════════════

pairs_wt_expr   <- link_loops_to_de_genes(wt_np_class,   gene_gr_expr, "Wildtype",   de_filter=FALSE)
pairs_pol2_expr <- link_loops_to_de_genes(pol2_np_class, gene_gr_expr, "POLR2A deg", de_filter=FALSE)
all_expr_pairs  <- bind_rows(pairs_wt_expr,pairs_pol2_expr) %>%
  dplyr::mutate(condition=factor(condition,levels=c("Wildtype","POLR2A deg")))
print(make_ccre_count_plot_paired(all_expr_pairs,NULL,
                                  "cCRE loop types at expressed genes (TPM>1) — narrowPeak + ATAC"))
print(plot_ccre_stacked(all_expr_pairs,NULL,
                        "cCRE loop types at expressed genes (TPM>1) — narrowPeak + ATAC"))


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 3: CTCF loop pair types — condition-matched
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 4: CTCF loop pair types ─────────────────────────────────")
wt_ctcf_np_pairs   <- label_ctcf_np_loops(wt_np_class)
pol2_ctcf_np_pairs <- label_ctcf_np_loops(pol2_np_class)
message("WT CTCF pair breakdown:");   print(table(wt_ctcf_np_pairs$ctcf_pair_type))
message("pol2 CTCF pair breakdown:"); print(table(pol2_ctcf_np_pairs$ctcf_pair_type))

plot_df_ctcf <- bind_rows(
  wt_ctcf_np_pairs   %>% dplyr::mutate(condition="Wildtype"),
  pol2_ctcf_np_pairs %>% dplyr::mutate(condition="POLR2A deg")
) %>%
  dplyr::mutate(condition=factor(condition,levels=c("Wildtype","POLR2A deg")),
                ctcf_pair_type=factor(ctcf_pair_type,levels=ctcf_np_pair_types)) %>%
  dplyr::count(ctcf_pair_type,condition) %>%
  tidyr::complete(ctcf_pair_type=ctcf_np_pair_types,condition=c("Wildtype","POLR2A deg"),
                  fill=list(n=0)) %>%
  dplyr::mutate(
    condition  = factor(condition, levels=c("Wildtype","POLR2A deg")),   # re-factor after complete
    fill_key   = factor(paste0(condition,"__",ctcf_pair_type),
                        levels=c(paste0("Wildtype__",ctcf_np_pair_types),
                                 paste0("POLR2A deg__",ctcf_np_pair_types)))
  )

fill_values_ctcf <- c(
  setNames(rep("gray70",length(ctcf_np_pair_types)),paste0("Wildtype__",ctcf_np_pair_types)),
  setNames(ctcf_np_pair_colors,paste0("POLR2A deg__",ctcf_np_pair_types))
)

print(ggplot(plot_df_ctcf,aes(x=ctcf_pair_type,y=n,fill=fill_key,group=condition)) +
        geom_col(position=position_dodge(width=0.7),width=0.65,color="white",linewidth=0.3,alpha=0.92) +
        geom_text(aes(label=ifelse(n>0,n,"")),position=position_dodge(width=0.7),
                  vjust=-0.4,size=3,color="gray30") +
        scale_fill_manual(values=fill_values_ctcf,
                          breaks=c("Wildtype__CTCF — CTCF","POLR2A deg__CTCF — CTCF"),
                          labels=c("WT (grey)","POLR2A deg (colored)"),name=NULL) +
        scale_x_discrete(labels=function(x) gsub(" — ","\n",x)) +
        scale_y_continuous(name="Number of loops",expand=expansion(mult=c(0,0.12))) +
        labs(title=paste0("CTCF loop pair types — narrowPeak + ATAC  |  ", DS_GROUP),
             subtitle="Condition-matched peaks  |  WT loops × WT CTCF  |  pol2 loops × KO CTCF",
             x=NULL) +
        theme_classic(base_size=12) +
        theme(plot.title=element_text(face="bold",size=13),
              plot.subtitle=element_text(size=8,color="gray40"),
              axis.text.x=element_text(size=10),legend.position="bottom"))

## THIS PLOTS CTCF HICCUPS OBSERVED NOT STRAWR

ctcf_np_combined <- bind_rows(
  wt_ctcf_np_pairs %>%
    dplyr::mutate(group=dplyr::case_when(ctcf_pair_type %in% names(CTCF_GROUP_MAP) ~
                                           CTCF_GROUP_MAP[ctcf_pair_type],TRUE~NA_character_),
                  condition="Wildtype"),
  pol2_ctcf_np_pairs %>%
    dplyr::mutate(group=dplyr::case_when(ctcf_pair_type %in% names(CTCF_GROUP_MAP) ~
                                           CTCF_GROUP_MAP[ctcf_pair_type],TRUE~NA_character_),
                  condition="POLR2A deg")
) %>%
  filter(!is.na(group),!is.na(.data[[SCORE_COL]]),.data[[SCORE_COL]]>0) %>%
  dplyr::mutate(condition=factor(condition,levels=c("Wildtype","POLR2A deg")),
                group=factor(group,levels=names(CTCF_GROUP_COLORS)),
                ridge_label=factor(paste0(group," — ",condition),levels=rev(ctcf_ridge_order)))

n_ctcf_np <- ctcf_np_combined %>% dplyr::count(group,condition) %>%
  dplyr::mutate(ridge_label=paste0(group," — ",condition),
                axis_label=paste0(group," — ",condition," (n=",n,")"))
ctcf_np_label_lookup <- setNames(n_ctcf_np$axis_label, n_ctcf_np$ridge_label)
ctcf_np_combined
print(ggplot(ctcf_np_combined,aes(x=log10(.data[[SCORE_COL]]),y=ridge_label,
                                  fill=ridge_label,color=ridge_label)) +
        ggridges::geom_density_ridges(scale=0.9,alpha=0.85,linewidth=0.6,quantile_lines=TRUE,quantiles=2) +
        scale_fill_manual(values=ctcf_ridge_fill,   guide="none") +
        scale_color_manual(values=ctcf_ridge_color, guide="none") +
        scale_y_discrete(labels=function(x) ctcf_np_label_lookup[x]) +
        scale_x_continuous(name=SCORE_LABEL,
                           limits=c(0, log10(CAP_VAL)))+
        labs(title=paste0("CTCF loop strength — narrowPeak + ATAC — WT vs POLR2A deg  |  ", DS_GROUP),
             subtitle="Condition-matched peaks  |  dark=WT  |  light=pol2  |  line=median",
             y=NULL) +
        theme_ridges(font_size=11) +
        theme(plot.title=element_text(face="bold",size=13),
              plot.subtitle=element_text(size=8,color="gray40"),
              axis.text.y=element_text(size=10)))
# ══════════════════════════════════════════════════════════════════════════════
# Analysis 4b: WT-only CTCF loop count bar
# Count of WT loops by CTCF group type (CTCF-CTCF, CTCF-Promoter, CTCF-Enhancer)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 4b: WT CTCF loop count bar ──────────────────────────────")
wt_ctcf_grouped <- wt_ctcf_np_pairs %>%
  dplyr::mutate(
    group = dplyr::case_when(
      ctcf_pair_type %in% names(CTCF_GROUP_MAP) ~ CTCF_GROUP_MAP[ctcf_pair_type],
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  dplyr::mutate(group=factor(group, levels=names(CTCF_GROUP_COLORS)))

bar_ctcf_wt <- wt_ctcf_grouped %>%
  dplyr::count(group) %>%
  dplyr::mutate(prop=n/sum(n), pct=paste0(round(prop*100,1),"%"))

print(ggplot(bar_ctcf_wt, aes(x=group, y=n, fill=group)) +
        geom_bar(stat="identity", width=0.65, color="white", linewidth=0.3) +
        geom_text(aes(label=paste0(n,"\n(",pct,")")), vjust=-0.3, size=3.5, color="gray30") +
        scale_fill_manual(values=CTCF_GROUP_COLORS, guide="none") +
        scale_y_continuous(name="Loop count", expand=expansion(mult=c(0,0.15))) +
        labs(title=paste0("CTCF loop counts — Wildtype  |  ", DS_GROUP),
             subtitle="narrowPeak + ATAC  |  WT HiCCUPS loops × WT CTCF peaks",
             x=NULL) +
        theme_classic(base_size=12) +
        theme(plot.title=element_text(face="bold",size=13),
              plot.subtitle=element_text(size=9,color="gray40"),
              axis.text.x=element_text(size=11,face="bold",
                                       color=CTCF_GROUP_COLORS[levels(wt_ctcf_grouped$group)])))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Matched ridgeline analyses — shared/exclusive loops
# Uses strawr matrix="oe" for both conditions → directly comparable obs/exp
# Replaces Analyses 3b, 3c, 4c with proper apples-to-apples comparison
# Adds shared/exclusive ridgelines for CTCF and EP/EE/PP loops
# ══════════════════════════════════════════════════════════════════════════════

# ── Updated extract_loop_contacts: now supports matrix="oe" ──────────────────
extract_loop_contacts <- function(loops_df, hic_file, norm="SCALE",
                                  resolution=5000, matrix="oe") {
  snap <- function(x) as.integer(floor(x / resolution) * resolution)
  results   <- rep(NA_real_, nrow(loops_df))
  c1_pos    <- snap(as.integer((loops_df$start1 + loops_df$end1) / 2))
  c2_pos    <- snap(as.integer((loops_df$start2 + loops_df$end2) / 2))
  chr1_clean <- gsub("^chr", "", loops_df$chr1)
  chr2_clean <- gsub("^chr", "", loops_df$chr2)
  chr_pairs  <- unique(paste(chr1_clean, chr2_clean, sep="|"))
  message("Extracting ", matrix, " from: ", basename(hic_file),
          "  norm=", norm, "  res=", resolution/1000, "kb  n=", nrow(loops_df))
  for (cp in chr_pairs) {
    parts  <- strsplit(cp, "\\|")[[1]]
    c1_chr <- parts[1]; c2_chr <- parts[2]
    idx    <- which(chr1_clean == c1_chr & chr2_clean == c2_chr)
    c1_min <- min(c1_pos[idx]); c1_max <- max(c1_pos[idx]) + resolution
    c2_min <- min(c2_pos[idx]); c2_max <- max(c2_pos[idx]) + resolution
    loc1   <- paste0(c1_chr, ":", c1_min, ":", c1_max)
    loc2   <- paste0(c2_chr, ":", c2_min, ":", c2_max)
    tryCatch({
      mat <- strawr::straw(norm, hic_file, loc1, loc2, "BP", resolution,
                           matrix=matrix)
      if (nrow(mat) == 0) { results[idx] <- NA; next }
      lookup <- setNames(mat$counts,
                         paste(pmin(mat$x, mat$y), pmax(mat$x, mat$y), sep="|"))
      keys   <- paste(pmin(c1_pos[idx], c2_pos[idx]),
                      pmax(c1_pos[idx], c2_pos[idx]), sep="|")
      vals   <- lookup[keys]
      results[idx] <- ifelse(is.na(vals), NA, vals)
    }, error=function(e) {
      message("  Warning: ", c1_chr, " — ", e$message)
    })
  }
  n_hit <- sum(!is.na(results) & results > 0, na.rm=TRUE)
  message("  Signal: ", n_hit, " / ", nrow(loops_df),
          " (", round(100*n_hit/nrow(loops_df),1), "%)")
  results
}

# ── Generic loop overlap function — finds shared loops between two sets ───────
find_shared_loops <- function(loops_a, loops_b, slop=5000L) {
  make_gr <- function(df) {
    list(
      a1 = GRanges(df$chr1, IRanges(df$start1, df$end1)),
      a2 = GRanges(df$chr2, IRanges(df$start2, df$end2))
    )
  }
  gr_a <- make_gr(loops_a); gr_b <- make_gr(loops_b)
  expand <- function(gr) resize(gr, width(gr) + 2*slop, fix="center")
  hits_a1 <- findOverlaps(expand(gr_a$a1), expand(gr_b$a1))
  hits_a2 <- findOverlaps(expand(gr_a$a2), expand(gr_b$a2))
  shared <- inner_join(
    as.data.frame(hits_a1) %>% dplyr::rename(a_idx=queryHits, b_idx=subjectHits),
    as.data.frame(hits_a2) %>% dplyr::rename(a_idx=queryHits, b_idx=subjectHits),
    by = c("a_idx","b_idx")
  ) %>% distinct(a_idx, b_idx)
  list(
    shared    = shared,
    a_only    = setdiff(seq_len(nrow(loops_a)), shared$a_idx),
    b_only    = setdiff(seq_len(nrow(loops_b)), shared$b_idx),
    n_shared  = nrow(shared),
    n_a_only  = nrow(loops_a) - n_distinct(shared$a_idx),
    n_b_only  = nrow(loops_b) - n_distinct(shared$b_idx)
  )
}

# ── Generic matched ridgeline builder ────────────────────────────────────────
# loops_a: loop set whose positions are used as anchors (scored in both .hic)
# label_a/b: condition labels
# group_col: column in loops_a to use for grouping ridge facets
# groups: named color vector (names = group values)
# ridge_order: factor levels for ridge_label (controls y-axis order)
# ridge_fill/color: named color vectors for fill/outline
make_matched_ridge <- function(loops_a,
                               hic_file_a, hic_file_b,
                               label_a = "Wildtype",
                               label_b = "POLR2A deg",
                               group_col,
                               groups,
                               ridge_order,
                               ridge_fill,
                               ridge_color,
                               norm       = "SCALE",
                               resolution = 5000,
                               title      = "",
                               subtitle   = "",
                               cap_val    = CAP_VAL,
                               min_n      = 10L) {
  message("Extracting obs/exp — ", label_a, " .hic at ", label_a, " loop positions...")
  oe_a <- extract_loop_contacts(loops_a, hic_file_a, norm=norm,
                                resolution=resolution, matrix="oe")
  message("Extracting obs/exp — ", label_b, " .hic at ", label_a, " loop positions...")
  oe_b <- extract_loop_contacts(loops_a, hic_file_b, norm=norm,
                                resolution=resolution, matrix="oe")
  
  ridge_df <- loops_a %>%
    dplyr::mutate(oe_a = oe_a, oe_b = oe_b,
                  group = .data[[group_col]]) %>%
    filter(!is.na(group)) %>%
    dplyr::select(group, oe_a, oe_b) %>%
    tidyr::pivot_longer(cols=c(oe_a, oe_b),
                        names_to="score_cond", values_to="score") %>%
    dplyr::mutate(
      score_cond  = factor(score_cond,
                           levels=c("oe_a","oe_b"),
                           labels=c(label_a, label_b)),
      ridge_label = factor(paste0(group," — ",score_cond),
                           levels=rev(ridge_order)),
      score       = replace_na(score, 0)
    ) %>%
    filter(!is.na(ridge_label))
  
  # Drop groups with too few loops
  group_counts <- ridge_df %>%
    filter(score_cond == label_a) %>%
    dplyr::count(group)
  keep_groups <- group_counts %>% filter(n >= min_n) %>% pull(group)
  ridge_df <- ridge_df %>% filter(group %in% keep_groups)
  if (nrow(ridge_df) == 0) {
    message("No groups with n >= ", min_n, " — skipping plot")
    return(invisible(NULL))
  }
  
  n_labels <- ridge_df %>% dplyr::count(group, score_cond) %>%
    dplyr::mutate(ridge_label = paste0(group," — ",score_cond),
                  axis_label  = paste0(group," — ",score_cond," (n=",n,")"))
  label_lookup <- setNames(n_labels$axis_label, n_labels$ridge_label)
  
  zero_pct <- ridge_df %>%
    filter(score_cond == label_b) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(pct_zero=round(100*mean(score==0),1), .groups="drop")
  message("\n% positions with 0 obs/exp in ", label_b, ":")
  print(zero_pct)
  
  median_df <- ridge_df %>%
    dplyr::group_by(ridge_label) %>%
    dplyr::summarise(med=median(log10(score+1),na.rm=TRUE), .groups="drop") %>%
    dplyr::mutate(y_num=as.numeric(ridge_label),
                  med_label=sprintf("%.2f", 10^med - 1))
  
  ggplot(ridge_df, aes(x=log10(score+1), y=ridge_label,
                       fill=ridge_label, color=ridge_label)) +
    ggridges::geom_density_ridges(scale=0.9, alpha=0.85, linewidth=0.6,
                                  quantile_lines=TRUE, quantiles=2) +
    geom_text(data=median_df, aes(x=med, y=y_num+0.15, label=med_label),
              color="white", size=2.8, fontface="bold", inherit.aes=FALSE) +
    scale_fill_manual(values=ridge_fill,   guide="none") +
    scale_color_manual(values=ridge_color, guide="none") +
    scale_y_discrete(labels=function(x) label_lookup[x]) +
    scale_x_continuous(name="Obs/exp (log10(n+1))",
                       limits=c(0, log10(cap_val+1))) +
    labs(title=title, subtitle=subtitle, y=NULL) +
    theme_ridges(font_size=11) +
    theme(plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=8,color="gray40"),
          axis.text.y=element_text(size=10))
}

# ══════════════════════════════════════════════════════════════════════════════
# Find shared / exclusive loops — CTCF
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Loop overlap: CTCF ────────────────────────────────────────────────")
ctcf_overlap <- find_shared_loops(wt_ctcf_np_pairs, pol2_ctcf_np_pairs, slop=5000L)
message("WT CTCF:           ", nrow(wt_ctcf_np_pairs))
message("pol2 CTCF:         ", nrow(pol2_ctcf_np_pairs))
message("Shared:            ", ctcf_overlap$n_shared,
        " (", round(100*ctcf_overlap$n_shared/nrow(wt_ctcf_np_pairs),1), "% of WT)")
message("WT only:           ", ctcf_overlap$n_a_only)
message("pol2 only:         ", ctcf_overlap$n_b_only)

# ── Find shared / exclusive loops — EP/EE/PP ─────────────────────────────────
message("\n── Loop overlap: EP/EE/PP ────────────────────────────────────────────")
ep_overlap <- find_shared_loops(wt_np_grouped, pol2_np_grouped, slop=5000L)
message("WT EP/EE/PP:       ", nrow(wt_np_grouped))
message("pol2 EP/EE/PP:     ", nrow(pol2_np_grouped))
message("Shared:            ", ep_overlap$n_shared,
        " (", round(100*ep_overlap$n_shared/nrow(wt_np_grouped),1), "% of WT)")
message("WT only:           ", ep_overlap$n_a_only)
message("pol2 only:         ", ep_overlap$n_b_only)

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 3b (updated): EP/EE/PP shared loops — obs/exp from both .hic files
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 3b: EP/EE/PP shared loops matched ridgeline ─────────────")


print(make_matched_ridge(
  loops_a     = wt_np_grouped[ep_overlap$shared$a_idx, ],
  hic_file_a  = hic_wt_file,
  hic_file_b  = hic_pol2_file,
  label_a     = "Wildtype",
  label_b     = "POLR2A deg",
  group_col   = "group",
  groups      = GROUP_COLORS,
  ridge_order = ridge_order_matched,
  ridge_fill  = ridge_fill_m,
  ridge_color = ridge_color_m,
  title       = paste0("Shared EP/EE/PP loops: WT vs pol2 obs/exp  |  ", DS_GROUP),
  subtitle    = paste0("Same loop positions  |  strawr SCALE oe  |  ",
                       ep_overlap$n_shared, " shared loops  |  dark=WT  |  light=pol2")
))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 3c (updated): EP/EE/PP WT-only loops — obs/exp from both .hic files
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 3c: EP/EE/PP WT-only loops ──────────────────────────────")
print(make_matched_ridge(
  loops_a     = wt_np_grouped[ep_overlap$a_only, ],
  hic_file_a  = hic_wt_file,
  hic_file_b  = hic_pol2_file,
  label_a     = "Wildtype",
  label_b     = "POLR2A deg",
  group_col   = "group",
  groups      = GROUP_COLORS,
  ridge_order = ridge_order_matched,
  ridge_fill  = ridge_fill_m,
  ridge_color = ridge_color_m,
  title       = paste0("WT-only EP/EE/PP loops: obs/exp at WT positions  |  ", DS_GROUP),
  subtitle    = paste0(ep_overlap$n_a_only, " WT-exclusive loops  |  ",
                       "strawr SCALE oe  |  dark=WT  |  light=pol2 at same positions")
))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 3d (new): EP/EE/PP pol2-only loops — obs/exp from both .hic files
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 3d: EP/EE/PP pol2-only loops ────────────────────────────")
print(make_matched_ridge(
  loops_a     = pol2_np_grouped[ep_overlap$b_only, ],
  hic_file_a  = hic_pol2_file,
  hic_file_b  = hic_wt_file,
  label_a     = "POLR2A deg",
  label_b     = "Wildtype",
  group_col   = "group",
  groups      = GROUP_COLORS,
  ridge_order = ridge_order_matched,
  ridge_fill  = ridge_fill_m,
  ridge_color = ridge_color_m,
  title       = paste0("pol2-only EP/EE/PP loops: obs/exp at pol2 positions  |  ", DS_GROUP),
  subtitle    = paste0(ep_overlap$n_b_only, " pol2-exclusive loops  |  ",
                       "strawr SCALE oe  |  dark=pol2  |  light=WT at same positions")
))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 4c (updated): CTCF shared loops — obs/exp from both .hic files
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 4c: CTCF shared loops matched ridgeline ─────────────────")
print(make_matched_ridge(
  loops_a     = wt_ctcf_np_pairs[ctcf_overlap$shared$a_idx, ],
  hic_file_a  = hic_wt_file,
  hic_file_b  = hic_pol2_file,
  label_a     = "Wildtype",
  label_b     = "POLR2A deg",
  group_col   = "ctcf_pair_type",
  groups      = ctcf_np_pair_colors,
  ridge_order = ctcf_ridge_order,
  ridge_fill  = ctcf_ridge_fill,
  ridge_color = ctcf_ridge_color,
  title       = paste0("Shared CTCF loops: WT vs pol2 obs/exp  |  ", DS_GROUP),
  subtitle    = paste0(ctcf_overlap$n_shared, " shared CTCF loops  |  ",
                       "strawr SCALE oe  |  dark=WT  |  light=pol2")
))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 4d (new): CTCF WT-only loops — obs/exp from both .hic files
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 4d: CTCF WT-only loops ──────────────────────────────────")
print(make_matched_ridge(
  loops_a     = wt_ctcf_np_pairs[ctcf_overlap$a_only, ],
  hic_file_a  = hic_wt_file,
  hic_file_b  = hic_pol2_file,
  label_a     = "Wildtype",
  label_b     = "POLR2A deg",
  group_col   = "ctcf_pair_type",
  groups      = ctcf_np_pair_colors,
  ridge_order = ctcf_ridge_order,
  ridge_fill  = ctcf_ridge_fill,
  ridge_color = ctcf_ridge_color,
  title       = paste0("WT-only CTCF loops: obs/exp at WT positions  |  ", DS_GROUP),
  subtitle    = paste0(ctcf_overlap$n_a_only, " WT-exclusive CTCF loops  |  ",
                       "strawr SCALE oe  |  dark=WT  |  light=pol2 at same positions")
))

# ══════════════════════════════════════════════════════════════════════════════
# Analysis 4e (new): CTCF pol2-only loops — obs/exp from both .hic files
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Analysis 4e: CTCF pol2-only loops ────────────────────────────────")
print(make_matched_ridge(
  loops_a     = pol2_ctcf_np_pairs[ctcf_overlap$b_only, ],
  hic_file_a  = hic_pol2_file,
  hic_file_b  = hic_wt_file,
  label_a     = "POLR2A deg",
  label_b     = "Wildtype",
  group_col   = "ctcf_pair_type",
  groups      = ctcf_np_pair_colors,
  ridge_order = ctcf_ridge_order,
  ridge_fill  = ctcf_ridge_fill,
  ridge_color = ctcf_ridge_color,
  title       = paste0("pol2-only CTCF loops: obs/exp at pol2 positions  |  ", DS_GROUP),
  subtitle    = paste0(ctcf_overlap$n_b_only, " pol2-exclusive CTCF loops  |  ",
                       "strawr SCALE oe  |  dark=pol2  |  light=WT at same positions")
))

# ══════════════════════════════════════════════════════════════════════════════
# 16. Pure CTCF-CTCF shared loops — WT vs pol2 ridgeline
#     Shared only, both anchors CTCF, no Promoter/Enhancer at either anchor.
#     Tests whether globally diminished contacts reflect E-P loop loss
#     or a true reduction at structural CTCF-CTCF loops.
# ══════════════════════════════════════════════════════════════════════════════
make_matched_ridge <- function(classified_df, group_colors, group_levels,
                               group_var, title, subtitle) {
  ridge_order <- as.vector(rbind(
    paste0(group_levels, " \u2014 Gained"),
    paste0(group_levels, " \u2014 Stable"),
    paste0(group_levels, " \u2014 Lost")
  ))
  ridge_fill <- unlist(lapply(group_levels, function(g) {
    col <- group_colors[[g]]
    setNames(c(scales::alpha(col,0.85), scales::alpha(col,0.55), scales::alpha(col,0.25)),
             paste0(g, " \u2014 ", c("Gained","Stable","Lost")))
  }))
  ridge_color <- unlist(lapply(group_levels, function(g) {
    setNames(rep(group_colors[[g]], 3),
             paste0(g, " \u2014 ", c("Gained","Stable","Lost")))
  }))
  base_df <- classified_df %>%
    filter(!is.na(.data[[group_var]]), .data[[group_var]] %in% group_levels) %>%
    dplyr::mutate(grp=factor(.data[[group_var]], levels=group_levels)) %>%
    dplyr::select(grp, loop_class, wt_contacts, pol2_contacts) %>%
    tidyr::pivot_longer(cols=c(wt_contacts, pol2_contacts),
                        names_to="cond", values_to="contacts") %>%
    filter(!is.na(contacts), contacts > 0) %>%
    dplyr::mutate(
      cond_label  = ifelse(cond == "wt_contacts", "WT", "POLR2A deg"),
      ridge_label = factor(paste0(grp, " \u2014 ", loop_class),
                           levels=rev(ridge_order))
    ) %>%
    filter(!is.na(ridge_label))
  wt_df   <- base_df %>% filter(cond_label == "WT")
  pol2_df <- base_df %>% filter(cond_label == "POLR2A deg")
  n_labels <- base_df %>%
    dplyr::group_by(ridge_label, cond_label) %>%
    dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
    tidyr::pivot_wider(names_from=cond_label, values_from=n, values_fill=0) %>%
    dplyr::mutate(axis_label=paste0(ridge_label,
                                    " (WT n=", WT, " / pol2 n=", `POLR2A deg`, ")"))
  label_lookup <- setNames(n_labels$axis_label, n_labels$ridge_label)
  pol2_medians <- pol2_df %>%
    dplyr::group_by(ridge_label) %>%
    dplyr::summarise(
      med_log = median(log10(contacts), na.rm=TRUE),
      y_num   = as.numeric(ridge_label[1]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(med_label=sprintf("%.2f", 10^med_log))
  ggplot(mapping=aes(x=log10(contacts), y=ridge_label)) +
    ggridges::geom_density_ridges(
      data=wt_df, fill=scales::alpha("#3d3d3d",0.60),
      color=scales::alpha("black",0.8), linewidth=0.7, linetype="dashed",
      scale=0.9, quantile_lines=TRUE, quantiles=2
    ) +
    ggridges::geom_density_ridges(
      data=pol2_df, aes(fill=ridge_label, color=ridge_label),
      alpha=0.70, linewidth=0.6, scale=0.9, quantile_lines=TRUE, quantiles=2
    ) +
    geom_text(data=pol2_medians,
              aes(x=med_log, y=y_num+0.15, label=med_label),
              inherit.aes=FALSE, color="white", size=2.5, fontface="bold") +
    scale_fill_manual(values=ridge_fill,   guide="none") +
    scale_color_manual(values=ridge_color, guide="none") +
    scale_y_discrete(labels=function(x) label_lookup[x]) +
    scale_x_continuous(name=SCORE_LABEL) +
    labs(title=title, subtitle=subtitle, y=NULL,
         caption="Charcoal dashed = WT  |  Color = POLR2A deg") +
    theme_ridges(font_size=11) +
    theme(plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=8,color="gray40"),
          plot.caption=element_text(size=8,color="gray40"),
          axis.text.y=element_text(size=9))
}

pure_ctcf_shared <- all_annotated %>%
  filter(
    source == "shared",
    anchor1_ctcf & anchor2_ctcf,
    !anchor1_ccre %in% c("Promoter","Proximal-Enhancer","Distal-Enhancer"),
    !anchor2_ccre %in% c("Promoter","Proximal-Enhancer","Distal-Enhancer")
  )

message("Pure CTCF-CTCF shared loops: ", nrow(pure_ctcf_shared))
message("anchor1_ccre:"); print(table(pure_ctcf_shared$anchor1_ccre))
message("anchor2_ccre:"); print(table(pure_ctcf_shared$anchor2_ccre))
message("loop_class:");   print(table(pure_ctcf_shared$loop_class))

print(make_matched_ridge(
  classified_df = pure_ctcf_shared,
  group_colors  = CLASS_COLORS,
  group_levels  = c("Gained","Stable","Lost"),
  group_var     = "loop_class",
  title         = paste0("Pure CTCF-CTCF shared loops — WT vs pol2  |  ", DS_GROUP),
  subtitle       = paste0(
    "Both anchors CTCF only — no Promoter/Enhancer co-occupancy  |  ",
    "shared loops only  |  n=", nrow(pure_ctcf_shared), "  |  ",
    "charcoal dashed=WT  color=POLR2A deg"
  )
))
