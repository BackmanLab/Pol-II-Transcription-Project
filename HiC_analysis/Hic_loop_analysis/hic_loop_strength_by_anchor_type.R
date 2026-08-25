##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Wildtype Hi-C loop strength by anchor type
# (Exon-Exon / Exon-cCRE / Exon-Intron / cCRE-cCRE / Intron-Intron).
#
# Extracted from Hic_loop_analysis/hic_general_analysis.R -- this is
# the WT panel of that script's already-active
# plot_hic_anchor_raincloud_wt_pol2() call (source:
# classify_hic_loop_anchors_wt() + plot_hic_anchor_raincloud_single()).
# Isolated into its own script for clarity; not new analysis logic.
##################################################################

library(strawr)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(ggdist)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ══════════════════════════════════════════════════════════════════════════════
# 0. Downsampling group switch — same convention as the APA / LFC pipelines
# ══════════════════════════════════════════════════════════════════════════════
DS_GROUP <- "ds_wt_pol2"   # <- "ds_wt_pol2" | "ds_all" | "unsampled"

BASE     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
LOOP_DIR <- file.path(BASE, "loops")

if (DS_GROUP == "ds_wt_pol2") {
  wt_loops   <- file.path(LOOP_DIR, "ds_wt_pol2/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  pol2_loops <- file.path(LOOP_DIR, "ds_wt_pol2/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  wt_hic     <- file.path(BASE, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic   <- file.path(BASE, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT <- 1.0; DEPTH_SCALE_POL2 <- 1.0; MIN_OBS_WT <- 0; MIN_OBS_POL2 <- 0
} else if (DS_GROUP == "ds_all") {
  wt_loops   <- file.path(LOOP_DIR, "ds_all/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  pol2_loops <- file.path(LOOP_DIR, "ds_all/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  wt_hic     <- file.path(BASE, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic   <- file.path(BASE, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT <- 1.0; DEPTH_SCALE_POL2 <- 1.0; MIN_OBS_WT <- 0; MIN_OBS_POL2 <- 0
} else if (DS_GROUP == "unsampled") {
  wt_loops   <- file.path(LOOP_DIR, "unsampled/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  pol2_loops <- file.path(LOOP_DIR, "unsampled/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  wt_hic     <- file.path(BASE, "mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic   <- file.path(BASE, "mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT <- 1.0; DEPTH_SCALE_POL2 <- 1.0; MIN_OBS_WT <- 0; MIN_OBS_POL2 <- 0
} else {
  stop("Unknown DS_GROUP: '", DS_GROUP, "'. Must be 'ds_wt_pol2', 'ds_all', or 'unsampled'.")
}
message("DS_GROUP: ", DS_GROUP)

ENCODE_CCRE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/encode/wt_hct116.bed"

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load HiCCUPS loop calls
# ══════════════════════════════════════════════════════════════════════════════

load_hiccups_bedpe <- function(path) {
  all_lines   <- readLines(path)
  header_line <- all_lines[grepl("^#chr1", all_lines)][1]
  col_names   <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  data_lines  <- all_lines[!grepl("^#", all_lines)]
  if (length(data_lines) == 0) stop("No data lines found in ", path)

  df <- read.table(text=paste(data_lines, collapse="\n"),
                   header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
  if (ncol(df) == length(col_names)) colnames(df) <- col_names

  out <- data.frame(
    chr1=as.character(df[[1]]), start1=as.integer(df[[2]]), end1=as.integer(df[[3]]),
    chr2=as.character(df[[4]]), start2=as.integer(df[[5]]), end2=as.integer(df[[6]]),
    observed=as.numeric(df[["observed"]]), fdrDonut=as.numeric(df[["fdrDonut"]]),
    fdrBL=as.numeric(df[["fdrBL"]]), fdrH=as.numeric(df[["fdrH"]]), fdrV=as.numeric(df[["fdrV"]]),
    centroid1=as.integer(df[["centroid1"]]), centroid2=as.integer(df[["centroid2"]])
  ) %>% filter(!is.na(start1), !is.na(start2))
  message("Loaded ", nrow(out), " loops from ", basename(path))
  out
}

filter_loops <- function(df, min_observed=10, max_fdr=0.1) {
  df %>% filter(observed >= min_observed, fdrDonut <= max_fdr, fdrBL <= max_fdr,
                fdrH <= max_fdr, fdrV <= max_fdr)
}

find_shared_loops <- function(loops_a, loops_b, slop=25000L) {
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
  list(shared=shared, a_only=setdiff(seq_len(nrow(loops_a)), shared$a_idx),
       b_only=setdiff(seq_len(nrow(loops_b)), shared$b_idx))
}

build_union_loops <- function(loops_wt, loops_deg, min_observed_wt=MIN_OBS_WT,
                              min_observed_deg=MIN_OBS_POL2, max_fdr=0.1, slop=25000L) {
  wt_f  <- filter_loops(loops_wt,  min_observed=min_observed_wt,  max_fdr=max_fdr) %>%
    dplyr::mutate(loop_dist=abs(centroid2-centroid1))
  deg_f <- filter_loops(loops_deg, min_observed=min_observed_deg, max_fdr=max_fdr) %>%
    dplyr::mutate(loop_dist=abs(centroid2-centroid1))
  overlap <- find_shared_loops(wt_f, deg_f, slop=slop)
  shared_wt <- wt_f %>% dplyr::slice(overlap$shared$a_idx) %>% dplyr::mutate(source="shared")
  wt_only_loops <- wt_f %>% dplyr::slice(overlap$a_only) %>% dplyr::mutate(source="wt_only")
  pol2_only_loops <- deg_f %>% dplyr::slice(overlap$b_only) %>% dplyr::mutate(source="pol2_only")
  bind_rows(shared_wt, wt_only_loops, pol2_only_loops) %>% dplyr::mutate(loop_idx=row_number())
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Real SCALE-normalized contact extraction via strawr
# ══════════════════════════════════════════════════════════════════════════════

STRAW_NORM       <- "SCALE"
STRAW_RESOLUTION <- 5000L

extract_loop_contacts <- function(loops_df, hic_file, norm=STRAW_NORM,
                                  resolution=STRAW_RESOLUTION, matrix="observed") {
  snap    <- function(x) as.integer(floor(x / resolution) * resolution)
  results <- rep(NA_real_, nrow(loops_df))
  has_centroids <- all(c("centroid1","centroid2") %in% colnames(loops_df)) &&
    !any(is.na(loops_df$centroid1)) && !any(is.na(loops_df$centroid2))
  if (has_centroids) {
    c1_pos <- snap(as.integer(loops_df$centroid1)); c2_pos <- snap(as.integer(loops_df$centroid2))
  } else {
    c1_pos <- snap(as.integer((loops_df$start1 + loops_df$end1) / 2))
    c2_pos <- snap(as.integer((loops_df$start2 + loops_df$end2) / 2))
  }
  chr1_clean <- gsub("^chr", "", loops_df$chr1); chr2_clean <- gsub("^chr", "", loops_df$chr2)
  chr_pairs  <- unique(paste(chr1_clean, chr2_clean, sep="|"))
  for (cp in chr_pairs) {
    parts  <- strsplit(cp, "\\|")[[1]]; c1_chr <- parts[1]; c2_chr <- parts[2]
    idx    <- which(chr1_clean == c1_chr & chr2_clean == c2_chr)
    c1_min <- min(c1_pos[idx]); c1_max <- max(c1_pos[idx]) + resolution
    c2_min <- min(c2_pos[idx]); c2_max <- max(c2_pos[idx]) + resolution
    loc1   <- paste0(c1_chr, ":", c1_min, ":", c1_max); loc2 <- paste0(c2_chr, ":", c2_min, ":", c2_max)
    tryCatch({
      mat <- strawr::straw(norm, hic_file, loc1, loc2, "BP", resolution, matrix=matrix)
      if (nrow(mat) == 0) { results[idx] <- NA; next }
      lookup <- setNames(mat$counts, paste(pmin(mat$x, mat$y), pmax(mat$x, mat$y), sep="|"))
      keys   <- paste(pmin(c1_pos[idx], c2_pos[idx]), pmax(c1_pos[idx], c2_pos[idx]), sep="|")
      vals   <- lookup[keys]
      results[idx] <- ifelse(is.na(vals), NA, vals)
    }, error=function(e) message("  Warning: ", c1_chr, " - ", e$message))
  }
  results
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Build union loop set with SCALE-normalized contacts
# ══════════════════════════════════════════════════════════════════════════════

loops_wt_df  <- load_hiccups_bedpe(wt_loops)
loops_deg_df <- load_hiccups_bedpe(pol2_loops)

PSEUDOCOUNT <- 1

matched <- build_union_loops(loops_wt_df, loops_deg_df, min_observed_wt=MIN_OBS_WT,
                             min_observed_deg=MIN_OBS_POL2, max_fdr=0.1, slop=25000L)

scale_wt_raw  <- extract_loop_contacts(matched, wt_hic,   norm=STRAW_NORM, resolution=STRAW_RESOLUTION)
scale_deg_raw <- extract_loop_contacts(matched, pol2_hic, norm=STRAW_NORM, resolution=STRAW_RESOLUTION)

matched <- matched %>%
  dplyr::mutate(
    obs_wt=tidyr::replace_na(scale_wt_raw, 0) * DEPTH_SCALE_WT,
    obs_deg=tidyr::replace_na(scale_deg_raw, 0) * DEPTH_SCALE_POL2,
    scale_wt=obs_wt, scale_deg=obs_deg
  )

# ══════════════════════════════════════════════════════════════════════════════
# 4. Anchor classification: exon / intron / cCRE, WT-only (single condition)
# ══════════════════════════════════════════════════════════════════════════════

load_ccre_gr <- function(ccre_path) {
  df <- read.table(ccre_path, header=FALSE, sep="\t", stringsAsFactors=FALSE) %>%
    dplyr::select(chr=V1, start=V2, end=V3, id=V4, ccre_type=V10) %>%
    filter(ccre_type %in% c("PLS", "pELS", "dELS", "CA-CTCF"))
  GRanges(seqnames=df$chr, ranges=IRanges(df$start, df$end), ccre_type=df$ccre_type)
}

message("Building exon/intron ranges from txdb...")
all_exons   <- reduce(exons(txdb))
all_genes   <- genes(txdb)
all_introns <- GenomicRanges::setdiff(all_genes, all_exons)
all_introns <- all_introns[!overlapsAny(all_introns, all_exons)]
exon_gr   <- all_exons
intron_gr <- all_introns
ccre_gr   <- load_ccre_gr(ENCODE_CCRE)

classify_hic_loop_anchors_wt <- function(hic_loops, exon_gr, intron_gr, ccre_gr=NULL,
                                         anchor_window=STRAW_RESOLUTION) {
  hic_loops <- hic_loops %>% dplyr::mutate(loop_idx=row_number())
  fix_chr <- function(gr, reference_gr) {
    ref_chroms <- as.character(seqnames(reference_gr)); has_chr_ref <- any(grepl("^chr", ref_chroms))
    gr_chroms  <- as.character(seqnames(gr)); has_chr_gr <- any(grepl("^chr", gr_chroms))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    gr
  }
  half_w <- anchor_window / 2
  anchor1_gr <- GRanges(seqnames=hic_loops$chr1, ranges=IRanges(hic_loops$centroid1-half_w, hic_loops$centroid1+half_w))
  anchor2_gr <- GRanges(seqnames=hic_loops$chr2, ranges=IRanges(hic_loops$centroid2-half_w, hic_loops$centroid2+half_w))
  anchor1_gr <- fix_chr(anchor1_gr, exon_gr); anchor2_gr <- fix_chr(anchor2_gr, exon_gr)
  a1_exon <- overlapsAny(anchor1_gr, exon_gr); a1_intron <- overlapsAny(anchor1_gr, intron_gr)
  a2_exon <- overlapsAny(anchor2_gr, exon_gr); a2_intron <- overlapsAny(anchor2_gr, intron_gr)
  if (!is.null(ccre_gr)) {
    ccre_gr_fixed <- fix_chr(ccre_gr, exon_gr)
    a1_ccre <- overlapsAny(anchor1_gr, ccre_gr_fixed); a2_ccre <- overlapsAny(anchor2_gr, ccre_gr_fixed)
  } else { a1_ccre <- rep(FALSE, length(anchor1_gr)); a2_ccre <- rep(FALSE, length(anchor2_gr)) }
  classify_anchor <- function(is_exon, is_intron, is_ccre) {
    dplyr::case_when(
      is_exon & is_ccre ~ "exon+cCRE", is_exon & !is_ccre ~ "exon",
      !is_exon & is_ccre ~ "cCRE", !is_exon & !is_ccre & is_intron ~ "intron",
      TRUE ~ "intergenic")
  }
  anchor1_class <- classify_anchor(a1_exon, a1_intron, a1_ccre)
  anchor2_class <- classify_anchor(a2_exon, a2_intron, a2_ccre)
  loop_type <- dplyr::case_when(
    grepl("exon", anchor1_class) & grepl("exon", anchor2_class) ~ "Exon-Exon",
    (grepl("exon", anchor1_class) & anchor2_class == "cCRE") |
      (anchor1_class == "cCRE" & grepl("exon", anchor2_class)) ~ "Exon-cCRE",
    (grepl("exon", anchor1_class) & anchor2_class == "intron") |
      (anchor1_class == "intron" & grepl("exon", anchor2_class)) ~ "Exon-Intron",
    anchor1_class == "cCRE" & anchor2_class == "cCRE" ~ "cCRE-cCRE",
    anchor1_class == "intron" & anchor2_class == "intron" ~ "Intron-Intron",
    TRUE ~ "Other")
  hic_loops %>% dplyr::mutate(
    anchor1_class=anchor1_class, anchor2_class=anchor2_class,
    loop_anchor_type=factor(loop_type, levels=c("Exon-Exon","Exon-cCRE","Exon-Intron",
                                                 "cCRE-cCRE","Intron-Intron","Other")))
}

hic_classified <- classify_hic_loop_anchors_wt(matched, exon_gr, intron_gr, ccre_gr=ccre_gr)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Plot — Wildtype panel only (Fig. 5g)
# ══════════════════════════════════════════════════════════════════════════════

ANCHOR_TYPE_COLORS <- c("Exon-Exon"="#d6604d", "Exon-cCRE"="#984ea3", "Exon-Intron"="#4DAF4A",
                        "cCRE-cCRE"="#ff7f00", "Intron-Intron"="#2166ac")
ANCHOR_TYPE_LEVELS <- names(ANCHOR_TYPE_COLORS)

plot_hic_anchor_raincloud_single <- function(classified_loops, value_col, source_filter, panel_title,
                                             group_levels=ANCHOR_TYPE_LEVELS,
                                             group_colors=ANCHOR_TYPE_COLORS, y_limits=NULL) {
  plot_df <- classified_loops %>%
    filter(loop_anchor_type %in% group_levels, source %in% source_filter,
           !is.na(.data[[value_col]]), .data[[value_col]] > 0) %>%
    dplyr::mutate(log10_observed=pmin(log10(.data[[value_col]] + 1), log10(1001)),
                  loop_anchor_type=factor(loop_anchor_type, levels=group_levels))
  counts <- plot_df %>% dplyr::count(loop_anchor_type) %>%
    dplyr::mutate(loop_anchor_type=factor(loop_anchor_type, levels=group_levels),
                  label=paste0(loop_anchor_type, "\n(n=", n, ")")) %>%
    dplyr::arrange(loop_anchor_type)
  plot_df <- plot_df %>% left_join(counts %>% dplyr::select(loop_anchor_type, label), by="loop_anchor_type") %>%
    dplyr::mutate(label=factor(label, levels=counts$label))
  ggplot(plot_df, aes(x=label, y=log10_observed, fill=loop_anchor_type, color=loop_anchor_type)) +
    ggdist::stat_halfeye(aes(thickness=after_stat(pdf * n)), side="right", scale=0.4, width=0.6,
                         .width=0, point_colour=NA, alpha=0.75, normalize="groups",
                         position=position_nudge(x=0.05)) +
    geom_boxplot(width=0.12, outlier.shape=NA, color="gray20", alpha=0.9, linewidth=0.6,
                position=position_nudge(x=-0.07)) +
    scale_fill_manual(values=group_colors, guide="none") +
    scale_color_manual(values=group_colors, guide="none") +
    scale_y_continuous(name=paste0("SCALE contacts (", STRAW_RESOLUTION/1000, "kb, log10)"),
                       limits=if (!is.null(y_limits)) y_limits else c(NA, log10(1001)*1.05),
                       breaks=c(0,1,2,3), labels=c("1","10","100","1000"),
                       expand=expansion(mult=c(0.02,0))) +
    scale_x_discrete(expand=expansion(add=c(0.3, 0.5))) +
    labs(title=panel_title, x=NULL) +
    theme_classic(base_size=13) +
    theme(plot.title=element_text(face="bold", size=14),
          axis.text.x=element_text(size=10, angle=30, hjust=1, vjust=1),
          plot.margin=margin(8, 16, 8, 8), aspect.ratio=1.4)
}

fig5g <- plot_hic_anchor_raincloud_single(
  hic_classified, value_col="scale_wt", source_filter=c("shared", "wt_only"),
  panel_title="Wildtype Hi-C"
)
print(fig5g)
