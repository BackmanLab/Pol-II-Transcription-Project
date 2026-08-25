##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# ChIA-PET locus plotter -- full self-contained pipeline for
# multi-track genome browser-style locus figures. Panel order:
# ChIA-PET arcs / gene annotation + cCREs / POLR2A CUT&Tag / Total
# RNA-seq / EU-seq.
#
# Two analyses:
#   1. Loop arcs from public POLR2A ChIA-PET, colored by contact
#      strength, with gene and cCRE annotation tracks
#   2. Signal tracks (CUT&Tag, Total RNA-seq, EU-seq) at the same
#      locus, WT vs. POLR2A-KO
#
# Figure output: Fig. 10a-b (MMP7/MMP20 and RAD50/IL13/IL4 loci)
##################################################################

library(InteractionSet)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
library(AnnotationDbi)
library(scales)

# ── File paths ─────────────────────────────────────────────────────────────────

CHIA_DIR  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/chia"
CHIA_REP1 <- file.path(CHIA_DIR, "pol2_rep1.bedpe")
CHIA_REP2 <- file.path(CHIA_DIR, "pol2_rep2.bedpe")

TOTAL_RNA_WT_BW <- "/projects/b1042/BackmanLab/Lucas/101823_rnaseq/coverage/cpm/WT.cpm.bigWig"
TOTAL_RNA_KO_BW <- "/projects/b1042/BackmanLab/Lucas/101823_rnaseq/coverage/cpm/POLR2A_8hrs_aux.cpm.bigWig"
EU_SEQ_WT_BW    <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/coverage/subtracted/WT.subtracted.bigWig"
EU_SEQ_KO_BW    <- "/projects/b1042/BackmanLab/Lucas/090124_euseq/coverage/subtracted/POLR2A_8hrs_aux.subtracted.bigWig"
ENCODE_CCRE     <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/encode/wt_hct116.bed"

# ENCODE standard cCRE colors
CCRE_COLORS <- c(
  "PLS"     = "#ff0000",   # Promoter-like signature
  "pELS"    = "#ffa500",   # Proximal enhancer-like
  "dELS"    = "#ffcd00",   # Distal enhancer-like
  "CA-CTCF" = "#00b0f0",   # CTCF-bound
  "CA-TF"   = "#be28e5",   # TF-bound
  "CA"      = "#06da93",   # Chromatin accessible only
  "TF"      = "#be28e5"    # TF-only
)

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load and merge ChIA-PET replicates
# ══════════════════════════════════════════════════════════════════════════════

load_chiapet_bedpe <- function(path) {
  df <- read.table(path, header=FALSE, sep="\t",
                   comment.char="#",
                   col.names=c("chr1","start1","end1",
                               "chr2","start2","end2","score"))
  data.frame(
    chr1   = df$chr1,
    start1 = as.integer(df$start1),
    end1   = as.integer(df$end1),
    chr2   = df$chr2,
    start2 = as.integer(df$start2),
    end2   = as.integer(df$end2),
    score  = as.numeric(df$score)
  )
}

merge_chiapet_loops <- function(rep1_df, rep2_df, tol=5000) {
  
  make_gi <- function(df) {
    GInteractions(
      anchor1 = GRanges(df$chr1, IRanges(df$start1, df$end1)),
      anchor2 = GRanges(df$chr2, IRanges(df$start2, df$end2)),
      score   = df$score
    )
  }
  
  gi1 <- make_gi(rep1_df)
  gi2 <- make_gi(rep2_df)
  
  gi1_expanded <- GInteractions(
    anchor1 = trim(resize(anchors(gi1,"first"),
                          width(anchors(gi1,"first"))  + 2*tol, fix="center")),
    anchor2 = trim(resize(anchors(gi1,"second"),
                          width(anchors(gi1,"second")) + 2*tol, fix="center")),
    score   = gi1$score
  )
  
  hits <- findOverlaps(gi1_expanded, gi2, use.region="both")
  
  matched <- data.frame(
    idx1 = queryHits(hits),
    idx2 = subjectHits(hits)
  ) %>%
    dplyr::mutate(
      chr1   = as.character(seqnames(anchors(gi1,"first")))[idx1],
      start1 = as.integer((start(anchors(gi1,"first"))[idx1]  +
                             start(anchors(gi2,"first"))[idx2])  / 2),
      end1   = as.integer((end(anchors(gi1,"first"))[idx1]    +
                             end(anchors(gi2,"first"))[idx2])    / 2),
      chr2   = as.character(seqnames(anchors(gi1,"second")))[idx1],
      start2 = as.integer((start(anchors(gi1,"second"))[idx1] +
                             start(anchors(gi2,"second"))[idx2]) / 2),
      end2   = as.integer((end(anchors(gi1,"second"))[idx1]   +
                             end(anchors(gi2,"second"))[idx2])   / 2),
      score  = (gi1$score[idx1] + gi2$score[idx2]) / 2,
      source = "both_replicates"
    ) %>%
    dplyr::select(chr1, start1, end1, chr2, start2, end2, score, source)
  
  rep1_only <- rep1_df[setdiff(seq_along(gi1), queryHits(hits)), ] %>%
    dplyr::mutate(source="rep1_only") %>%
    dplyr::select(chr1, start1, end1, chr2, start2, end2, score, source)
  
  rep2_only <- rep2_df[setdiff(seq_along(gi2), subjectHits(hits)), ] %>%
    dplyr::mutate(source="rep2_only") %>%
    dplyr::select(chr1, start1, end1, chr2, start2, end2, score, source)
  
  bind_rows(matched, rep1_only, rep2_only)
}

message("Loading ChIA-PET replicates...")
rep1_df <- load_chiapet_bedpe(CHIA_REP1)
rep2_df <- load_chiapet_bedpe(CHIA_REP2)
message("  Rep1: ", nrow(rep1_df), " loops")
message("  Rep2: ", nrow(rep2_df), " loops")

merged_loops <- merge_chiapet_loops(rep1_df, rep2_df, tol=5000)
message("Merged loop source breakdown:")
print(table(merged_loops$source))
message("Score distribution:")
print(quantile(merged_loops$score, probs=c(0, 0.25, 0.5, 0.75, 0.9, 1)))

loops <- merged_loops %>% filter(score >= 2)
message(nrow(loops), " loops retained after PET >= 2 filter")
print(table(loops$source))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Helpers
# ══════════════════════════════════════════════════════════════════════════════

bezier_arc <- function(x1, x2, height, n=120) {
  t  <- seq(0, 1, length.out=n)
  cx <- (x1 + x2) / 2
  cy <- height
  x  <- (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2
  y  <- (1-t)^2 * 0  + 2*(1-t)*t * cy + t^2 * 0
  data.frame(x=x, y=y)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Gene annotation builder
# ══════════════════════════════════════════════════════════════════════════════

build_gene_annotation <- function(chrom, roi_start, roi_end, txdb, orgdb) {
  
  roi_gr    <- GRanges(chrom, IRanges(roi_start, roi_end))
  all_genes <- suppressMessages(subsetByOverlaps(genes(txdb), roi_gr))
  
  if (length(all_genes) == 0) {
    message("  No genes found in region")
    return(NULL)
  }
  
  genes_df <- as.data.frame(all_genes, row.names=NULL)
  
  sym_map <- suppressMessages(
    AnnotationDbi::select(orgdb,
                          keys    = genes_df$gene_id,
                          keytype = "ENTREZID",
                          columns = c("SYMBOL","GENENAME"))
  )
  genes_df <- left_join(genes_df, sym_map, by=c("gene_id"="ENTREZID")) %>%
    dplyr::mutate(
      SYMBOL      = ifelse(is.na(SYMBOL), gene_id, SYMBOL),
      is_lncrna   = grepl("long non-coding|lncRNA|antisense|LINC",
                          GENENAME, ignore.case=TRUE),
      gene_type   = ifelse(is_lncrna, "lncRNA", "protein_coding"),
      disp_start  = pmax(start, roi_start),
      disp_end    = pmin(end,   roi_end),
      mid_x       = (disp_start + disp_end) / 2,
      strand_char = as.character(strand)
    )
  
  # Greedy lane packing
  genes_df  <- genes_df %>% dplyr::arrange(start)
  lane_ends <- c()
  genes_df$lane <- 0L
  for (i in seq_len(nrow(genes_df))) {
    placed <- FALSE
    for (ln in seq_along(lane_ends)) {
      if (genes_df$start[i] > lane_ends[ln] + 5000) {
        genes_df$lane[i] <- ln
        lane_ends[ln]     <- genes_df$end[i]
        placed            <- TRUE
        break
      }
    }
    if (!placed) {
      lane_ends        <- c(lane_ends, genes_df$end[i])
      genes_df$lane[i] <- length(lane_ends)
    }
  }
  genes_df$y <- -(genes_df$lane - 1L) * 0.65
  
  gene_y_map <- setNames(genes_df$y, genes_df$gene_id)
  
  tx_by_gene  <- suppressMessages(transcriptsBy(txdb, by="gene"))
  exons_by_tx <- suppressMessages(exonsBy(txdb, by="tx"))
  utr5_by_tx  <- suppressMessages(fiveUTRsByTranscript(txdb))
  utr3_by_tx  <- suppressMessages(threeUTRsByTranscript(txdb))
  tx_lengths  <- suppressMessages(
    transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
  )
  
  exons_df   <- data.frame()
  introns_df <- data.frame()
  utrs_df    <- data.frame()
  arrows_df  <- data.frame()
  chrom_bare <- gsub("^chr","", chrom)
  
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
    raw_exons <- as.data.frame(exons_by_tx[[canon_key]], row.names=NULL) %>%
      filter(gsub("^chr","", as.character(seqnames)) == chrom_bare,
             start <= gene_right, end >= gene_left) %>%
      dplyr::mutate(start     = pmax(start, gene_left),
                    end       = pmin(end,   gene_right),
                    gene_id   = gid,
                    gene_type = gene_type,
                    y_base    = y_base) %>%
      filter(end > start)
    
    if (nrow(raw_exons) == 0) next
    exons_df <- bind_rows(exons_df, raw_exons)
    
    get_utr_df <- function(utr_list, key, utype) {
      if (!key %in% names(utr_list)) return(NULL)
      as.data.frame(utr_list[[key]], row.names=NULL) %>%
        dplyr::mutate(utr_type=utype, gene_id=gid, y_base=y_base)
    }
    utrs_df <- bind_rows(utrs_df,
                         get_utr_df(utr5_by_tx, canon_key, "5UTR"),
                         get_utr_df(utr3_by_tx, canon_key, "3UTR"))
    
    if (nrow(raw_exons) >= 2) {
      ex_ord <- raw_exons %>% dplyr::arrange(start)
      introns <- data.frame(
        intron_start = ex_ord$end[-nrow(ex_ord)],
        intron_end   = ex_ord$start[-1],
        gene_id      = gid,
        gene_type    = gene_type,
        y_base       = y_base
      ) %>% filter(intron_end > intron_start)
      introns_df <- bind_rows(introns_df, introns)
    }
    
    gene_span <- gene_right - gene_left
    n_arr     <- max(1L, floor(gene_span / 60000))
    arr_x     <- seq(gene_left  + gene_span / (2*n_arr),
                     gene_right - gene_span / (2*n_arr),
                     length.out = n_arr)
    
    in_intron_fn <- function(ax) {
      if (nrow(introns_df) == 0) return(FALSE)
      sub <- introns_df[introns_df$gene_id == gid, ]
      any(sub$intron_start <= ax & sub$intron_end >= ax)
    }
    
    arrows_df <- bind_rows(arrows_df, data.frame(
      x         = arr_x,
      strand    = gene_strand,
      gene_id   = gid,
      y_base    = y_base,
      in_intron = vapply(arr_x, in_intron_fn, logical(1))
    ))
  }
  
  list(genes_df=genes_df, exons_df=exons_df, introns_df=introns_df,
       utrs_df=utrs_df, arrows_df=arrows_df)
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. BigWig importer
# ══════════════════════════════════════════════════════════════════════════════

import_bigwig_region <- function(bw_path, chrom, start, end, n_bins=500) {
  
  if (is.null(bw_path) || !file.exists(bw_path)) {
    message("  BigWig not found: ", bw_path)
    return(NULL)
  }
  
  roi_gr     <- GRanges(chrom, IRanges(start, end))
  chrom_bare <- gsub("^chr","", chrom)
  
  tryCatch({
    vals <- rtracklayer::import(bw_path, format="BigWig", which=roi_gr)
    if (length(vals) == 0) return(NULL)
    
    vals_df <- as.data.frame(vals, row.names=NULL) %>%
      filter(gsub("^chr","", as.character(seqnames)) == chrom_bare)
    if (nrow(vals_df) == 0) return(NULL)
    
    bin_width <- (end - start) / n_bins
    bins <- data.frame(
      bin_mid   = seq(start + bin_width/2, end,           length.out=n_bins),
      bin_start = seq(start,               end-bin_width, length.out=n_bins),
      bin_end   = seq(start + bin_width,   end,           length.out=n_bins),
      score     = 0
    )
    
    for (i in seq_len(nrow(vals_df))) {
      v_s <- vals_df$start[i]; v_e <- vals_df$end[i]; v_v <- vals_df$score[i]
      ov  <- bins$bin_start < v_e & bins$bin_end > v_s
      if (!any(ov)) next
      overlap_len    <- pmin(bins$bin_end[ov], v_e) -
        pmax(bins$bin_start[ov], v_s)
      bins$score[ov] <- bins$score[ov] + v_v * overlap_len / bin_width
    }
    bins
    
  }, error=function(e) {
    message("  BigWig error (", basename(bw_path), "): ", e$message)
    NULL
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Main plotter
# ══════════════════════════════════════════════════════════════════════════════

plot_chiapet_region <- function(
    chrom, start, end,
    loops_df,
    txdb, orgdb,
    total_rna_wt_bw = TOTAL_RNA_WT_BW,
    total_rna_ko_bw = TOTAL_RNA_KO_BW,
    eu_seq_wt_bw    = EU_SEQ_WT_BW,
    eu_seq_ko_bw    = EU_SEQ_KO_BW,
    pol2_wt_bw      = "/projects/b1042/BackmanLab/Lucas/eric_splicing/cut/Binned_Coverage/WT_P2AB.bw",
    ccre_path       = ENCODE_CCRE,
    label           = NULL,
    min_pet         = 2,
    n_signal_bins   = 400,
    source_filter   = NULL
) {
  
  pad       <- (end - start) * 0.05
  roi_start <- start - pad
  roi_end   <- end   + pad
  
  message("Region: ", chrom, ":", format(roi_start, big.mark=","),
          " - ", format(roi_end, big.mark=","))
  
  # ── Both anchors must fall within ROI ─────────────────────────────────────
  loops_roi <- loops_df %>%
    filter(
      chr1   == chrom,
      chr2   == chrom,
      score  >= min_pet,
      start1 >= roi_start,
      end1   <= roi_end,
      start2 >= roi_start,
      end2   <= roi_end
    )
  
  if (!is.null(source_filter) && "source" %in% names(loops_roi)) {
    loops_roi <- loops_roi %>% filter(source %in% source_filter)
  }
  
  message(nrow(loops_roi), " loops with both anchors in ROI (PET >= ",
          min_pet, ")")
  if ("source" %in% names(loops_roi)) print(table(loops_roi$source))
  
  if (nrow(loops_roi) == 0) {
    message("No loops found — try lowering min_pet, check chrom format, ",
            "or widen the region coordinates")
    return(invisible(NULL))
  }
  
  loops_roi <- loops_roi %>%
    dplyr::mutate(
      x1 = (start1 + end1) / 2,
      x2 = (start2 + end2) / 2
    )
  
  max_pet       <- max(loops_roi$score)
  min_pet_shown <- min(loops_roi$score)
  
  # ── Arc building ───────────────────────────────────────────────────────────
  # Each loop produces two outputs from the same filtered pass: the bezier
  # curve (still anchored at the window midpoints, for the arc's shape) and
  # an anchor-width span pair — one segment per anchor, drawn at baseline —
  # so the plot shows the true called window instead of implying a single
  # precise coordinate the data doesn't actually have.
  arc_list    <- list()
  anchor_list <- list()
  for (i in seq_len(nrow(loops_roi))) {
    row <- loops_roi[i, ]
    if (abs(row$x2 - row$x1) < 100) next
    scaled     <- ((row$score - min_pet_shown) /
                     max(max_pet - min_pet_shown, 1e-6)) ^ 0.5
    span       <- abs(row$x2 - row$x1)
    max_height <- min(span / 2, 10)
    arc_height <- max_height * (0.20 + 0.80 * scaled)
    df         <- bezier_arc(row$x1, row$x2, arc_height)
    df$loop_id <- i
    df$score   <- row$score
    df$source  <- if ("source" %in% names(row)) row$source else "unknown"
    arc_list[[length(arc_list) + 1]] <- df
    
    anchor_list[[length(anchor_list) + 1]] <- data.frame(
      xmin    = c(row$start1, row$start2),
      xmax    = c(row$end1,   row$end2),
      loop_id = i,
      score   = row$score
    )
  }
  arc_data     <- bind_rows(arc_list)
  anchor_spans <- bind_rows(anchor_list)
  arc_y_max <- max(arc_data$y, na.rm=TRUE)
  message(nrow(arc_data) / 120, " arcs rendered, max height: ",
          round(arc_y_max, 2))
  
  # ── Gene annotation ────────────────────────────────────────────────────────
  message("Building gene annotation...")
  ann <- build_gene_annotation(chrom, roi_start, roi_end, txdb, orgdb)
  
  has_genes   <- !is.null(ann) && nrow(ann$genes_df)   > 0
  has_exons   <- !is.null(ann) && nrow(ann$exons_df)   > 0
  has_introns <- !is.null(ann) && nrow(ann$introns_df) > 0
  has_utrs    <- !is.null(ann) && nrow(ann$utrs_df)    > 0
  has_arrows  <- !is.null(ann) && nrow(ann$arrows_df)  > 0
  n_lanes     <- if (has_genes) max(ann$genes_df$lane) else 1L
  
  # ── Load cCREs — informative types only ───────────────────────────────────
  message("Loading cCREs...")
  ccres <- tryCatch({
    df <- read.table(ccre_path, header=FALSE, sep="\t",
                     stringsAsFactors=FALSE) %>%
      dplyr::select(chr=V1, start=V2, end=V3, id=V4, ccre_type=V10) %>%
      filter(
        chr      == chrom,
        start    <= roi_end,
        end      >= roi_start,
        ccre_type %in% c("PLS","pELS","dELS","CA-CTCF","CA-TF","TF")
      )
    if (nrow(df) == 0) {
      message("  No informative cCREs in region"); NULL
    } else {
      message("  ", nrow(df), " cCREs in region")
      print(table(df$ccre_type))
      df
    }
  }, error=function(e) {
    message("  cCRE load error: ", e$message); NULL
  })
  has_ccres <- !is.null(ccres) && nrow(ccres) > 0
  
  # ── BigWig signals ─────────────────────────────────────────────────────────
  message("Importing bigWig signals...")
  sig_pol2_wt  <- import_bigwig_region(pol2_wt_bw,      chrom,
                                       roi_start, roi_end, n_signal_bins)
  sig_total_wt <- import_bigwig_region(total_rna_wt_bw, chrom,
                                       roi_start, roi_end, n_signal_bins)
  sig_total_ko <- import_bigwig_region(total_rna_ko_bw, chrom,
                                       roi_start, roi_end, n_signal_bins)
  sig_eu_wt    <- import_bigwig_region(eu_seq_wt_bw,    chrom,
                                       roi_start, roi_end, n_signal_bins)
  sig_eu_ko    <- import_bigwig_region(eu_seq_ko_bw,    chrom,
                                       roi_start, roi_end, n_signal_bins)
  
  pol2_ymax <- if (!is.null(sig_pol2_wt)) {
    max(sig_pol2_wt$score, na.rm=TRUE) * 1.1
  } else 0.01
  
  total_ymax <- max(
    if (!is.null(sig_total_wt)) max(sig_total_wt$score, na.rm=TRUE) else 0,
    if (!is.null(sig_total_ko)) max(sig_total_ko$score, na.rm=TRUE) else 0
  ) * 1.1
  
  eu_ymax <- max(
    if (!is.null(sig_eu_wt)) max(sig_eu_wt$score, na.rm=TRUE) else 0,
    if (!is.null(sig_eu_ko)) max(sig_eu_ko$score, na.rm=TRUE) else 0
  ) * 1.1
  
  # ── x limits ──────────────────────────────────────────────────────────────
  x_lim <- c(roi_start, roi_end)
  
  plot_title <- if (!is.null(label)) {
    paste0("POLR2A ChIA-PET — ", label)
  } else {
    paste0(chrom, ":", format(start, big.mark=","),
           "–", format(end, big.mark=","))
  }
  
  # ── cCRE and fill scale definitions ───────────────────────────────────────
  CCRE_COLORS <- c(
    "PLS"     = "#e31a1c",
    "pELS"    = "#ff7f00",
    "dELS"    = "#fdbf6f",
    "CA-CTCF" = "#1f78b4",
    "CA-TF"   = "#6a3d9a",
    "TF"      = "#b15928"
  )
  
  all_fill_values <- c(
    "5UTR"    = "#41ab5d",
    "3UTR"    = "#fd8d3c",
    "PLS"     = "#e31a1c",
    "pELS"    = "#ff7f00",
    "dELS"    = "#fdbf6f",
    "CA-CTCF" = "#1f78b4",
    "CA-TF"   = "#6a3d9a",
    "TF"      = "#b15928"
  )
  
  all_fill_labels <- c(
    "5UTR"    = "5' UTR",
    "3UTR"    = "3' UTR",
    "PLS"     = "Promoter (PLS)",
    "pELS"    = "Prox. enhancer",
    "dELS"    = "Dist. enhancer",
    "CA-CTCF" = "CTCF-bound",
    "CA-TF"   = "TF-bound",
    "TF"      = "TF-only"
  )
  
  # ══════════════════════════════════════════════════════════════════════════
  # Panel 1: ChIA-PET arcs
  # ══════════════════════════════════════════════════════════════════════════
  p_loops <- ggplot() +
    geom_path(data=arc_data,
              aes(x=x, y=y, group=loop_id, color=score),
              linewidth=0.6, alpha=0.75) +
    # Anchor windows, drawn as thick baseline segments under each arc foot —
    # this is the actual called peak width (V2:V3 / V5:V6 in the bedpe), not
    # a single coordinate. Same color mapping as the arcs so it reads as
    # "the wide version of that arc's foot" rather than a separate feature.
    geom_segment(data=anchor_spans,
                 aes(x=xmin, xend=xmax, y=0, yend=0, color=score),
                 linewidth=3, lineend="butt", alpha=0.85) +
    scale_color_gradient2(
      low="#6baed6", mid="#2171b5", high="#08306b",
      midpoint=(max_pet + min_pet_shown) / 2,
      name="PET count"
    ) +
    scale_x_continuous(limits=x_lim, expand=c(0, 0),
                       labels=NULL, name=NULL) +
    scale_y_continuous(limits=c(0, arc_y_max * 1.08), expand=c(0, 0)) +
    labs(title=plot_title, y=NULL) +
    theme_classic(base_size=11) +
    theme(
      axis.line.x       = element_blank(),
      axis.ticks.x      = element_blank(),
      axis.text.x       = element_blank(),
      axis.title.x      = element_blank(),
      axis.line.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.text.y       = element_blank(),
      plot.title        = element_text(face="bold", size=12),
      legend.position   = "right",
      legend.key.height = unit(0.7,"cm"),
      plot.margin       = margin(4, 6, 0, 6)
    )
  
  # ══════════════════════════════════════════════════════════════════════════
  # Panel 2: Gene annotation + cCREs
  # ══════════════════════════════════════════════════════════════════════════
  EXON_H <- 0.16
  UTR_H  <- 0.09
  CCRE_Y <- 0.65
  CCRE_H <- 0.18
  
  gene_y_top       <- if (has_ccres) CCRE_Y + CCRE_H + 0.15 else EXON_H + 0.38
  gene_y_bottom    <- -(n_lanes - 1) * 0.65
  gene_y_min       <- gene_y_bottom - EXON_H - 0.55
  gene_y_max       <- gene_y_top
  gene_panel_range <- gene_y_max - gene_y_min
  
  p_genes <- ggplot()
  
  # ── cCRE row ──────────────────────────────────────────────────────────────
  if (has_ccres) {
    region_span <- roi_end - roi_start
    min_width   <- region_span * 0.003
    ccres_plot  <- ccres %>%
      dplyr::mutate(
        plot_end  = pmax(end, start + min_width),
        ccre_type = ifelse(ccre_type %in% names(CCRE_COLORS),
                           ccre_type, "TF")
      )
    p_genes <- p_genes +
      annotate("segment",
               x=roi_start, xend=roi_end,
               y=CCRE_Y - CCRE_H - 0.06, yend=CCRE_Y - CCRE_H - 0.06,
               color="gray88", linewidth=0.3) +
      geom_rect(data=ccres_plot,
                aes(xmin=start, xmax=plot_end,
                    ymin=CCRE_Y - CCRE_H, ymax=CCRE_Y + CCRE_H,
                    fill=ccre_type),
                color=NA, alpha=0.9) +
      annotate("text",
               x=roi_start + (roi_end - roi_start) * 0.005, y=CCRE_Y,
               label="cCREs", hjust=0, vjust=0.5,
               size=2.4, color="gray30", fontface="bold")
  }
  
  # ── Intron lines ──────────────────────────────────────────────────────────
  if (has_introns) {
    p_genes <- p_genes +
      geom_segment(data=ann$introns_df,
                   aes(x=intron_start, xend=intron_end,
                       y=y_base, yend=y_base, linetype=gene_type),
                   color="gray25", linewidth=0.45)
  }
  
  # ── Strand arrows ─────────────────────────────────────────────────────────
  if (has_arrows) {
    arr <- ann$arrows_df %>%
      dplyr::mutate(
        dx     = ifelse(strand=="+",  8000, -8000),
        x_from = x - abs(dx)/2,
        x_to   = x + abs(dx)/2 * ifelse(strand=="+", 1, -1)
      )
    p_genes <- p_genes +
      geom_segment(data=arr,
                   aes(x=x_from, xend=x_to, y=y_base, yend=y_base),
                   color="gray50", linewidth=0.4,
                   arrow=arrow(length=unit(0.07,"cm"),
                               type="closed", ends="last"))
  }
  
  # ── UTR blocks ────────────────────────────────────────────────────────────
  if (has_utrs) {
    p_genes <- p_genes +
      geom_rect(data=ann$utrs_df,
                aes(xmin=start, xmax=end,
                    ymin=y_base - UTR_H, ymax=y_base + UTR_H,
                    fill=utr_type),
                color=NA)
  }
  
  # ── CDS exon blocks ───────────────────────────────────────────────────────
  if (has_exons) {
    exon_colors <- setNames(
      ifelse(ann$genes_df$gene_type=="lncRNA", "#9ecae1", "#252525"),
      ann$genes_df$gene_id
    )
    p_genes <- p_genes +
      geom_rect(data=ann$exons_df,
                aes(xmin=start, xmax=end,
                    ymin=y_base - EXON_H, ymax=y_base + EXON_H),
                fill=exon_colors[ann$exons_df$gene_id],
                color=NA)
  }
  
  # ── TSS / TES / 5' / 3' / gene name labels ────────────────────────────────
  if (has_genes) {
    tss_tes <- ann$genes_df %>%
      dplyr::mutate(
        tss = ifelse(strand_char=="+", start, end),
        tes = ifelse(strand_char=="+", end,   start)
      )
    p_genes <- p_genes +
      geom_segment(data=tss_tes,
                   aes(x=tss, xend=tss,
                       y=y - EXON_H - 0.07, yend=y + EXON_H + 0.07),
                   color="#1a9850", linewidth=1.1) +
      geom_text(data=tss_tes,
                aes(x=tss, y=y + EXON_H + 0.14, label="TSS"),
                size=2.3, color="#1a9850", fontface="bold", hjust=0.5) +
      geom_segment(data=tss_tes,
                   aes(x=tes, xend=tes,
                       y=y - EXON_H - 0.07, yend=y + EXON_H + 0.07),
                   color="#d73027", linewidth=0.9, linetype="dashed") +
      geom_text(data=tss_tes,
                aes(x=tes, y=y + EXON_H + 0.14, label="TES"),
                size=2.3, color="#d73027", fontface="bold", hjust=0.5) +
      geom_text(data=tss_tes,
                aes(x=tss, y=y - EXON_H - 0.18, label="5'"),
                size=2.5, color="gray40", hjust=0.5) +
      geom_text(data=tss_tes,
                aes(x=tes, y=y - EXON_H - 0.18, label="3'"),
                size=2.5, color="gray40", hjust=0.5) +
      geom_text(data=ann$genes_df,
                aes(x=mid_x, y=y - EXON_H - 0.32, label=SYMBOL,
                    fontface=ifelse(gene_type=="lncRNA",
                                    "italic","bold.italic")),
                size=3.0, color="gray15")
  }
  
  # ── Combined fill scale ────────────────────────────────────────────────────
  utr_types_present  <- if (has_utrs)   unique(ann$utrs_df$utr_type) else c()
  ccre_types_present <- if (has_ccres)  unique(ccres$ccre_type)      else c()
  
  legend_breaks <- intersect(
    names(all_fill_values),
    c(utr_types_present, ccre_types_present)
  )
  legend_labels <- all_fill_labels[legend_breaks]
  
  p_genes <- p_genes +
    scale_linetype_manual(
      values = c(protein_coding="solid", lncRNA="dashed"),
      name   = "Gene type",
      guide  = guide_legend(override.aes=list(linewidth=0.8))
    ) +
    scale_fill_manual(
      values      = all_fill_values,
      breaks      = legend_breaks,
      labels      = legend_labels,
      name        = "Annotation",
      na.translate= FALSE
    ) +
    scale_x_continuous(limits=x_lim, expand=c(0, 0),
                       labels=NULL, name=NULL) +
    scale_y_continuous(limits=c(gene_y_min, gene_y_max), expand=c(0, 0)) +
    theme_classic(base_size=10) +
    theme(
      axis.line.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.x      = element_blank(),
      axis.title.x     = element_blank(),
      axis.line.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      axis.title.y     = element_blank(),
      legend.position  = "right",
      legend.text      = element_text(size=8),
      legend.title     = element_text(size=9, face="bold"),
      plot.margin      = margin(0, 6, 0, 6)
    )
  
  # ══════════════════════════════════════════════════════════════════════════
  # Signal track helper
  # ══════════════════════════════════════════════════════════════════════════
  make_signal_panel <- function(sig_wt, sig_ko=NULL,
                                fill_color, ko_color=NULL,
                                y_max, y_label,
                                bottom_panel=FALSE) {
    p <- ggplot()
    
    if (!is.null(sig_wt)) {
      p <- p +
        geom_area(data=sig_wt, aes(x=bin_mid, y=score),
                  fill=fill_color, alpha=0.45, color=NA) +
        geom_line(data=sig_wt, aes(x=bin_mid, y=score, color="WT"),
                  linewidth=0.4)
    }
    if (!is.null(sig_ko) && !is.null(ko_color)) {
      p <- p +
        geom_line(data=sig_ko, aes(x=bin_mid, y=score, color="POLR2A deg"),
                  linewidth=0.55)
    }
    
    # Build color scale — only include conditions present
    cond_colors <- c(WT=fill_color)
    if (!is.null(sig_ko) && !is.null(ko_color)) {
      cond_colors <- c(cond_colors, "POLR2A deg"=ko_color)
    }
    
    p <- p +
      scale_color_manual(values=cond_colors, name=NULL) +
      scale_x_continuous(
        limits = x_lim,
        expand = c(0, 0),
        labels = if (bottom_panel) {
          label_number(scale=1e-6, suffix=" Mb", accuracy=0.01)
        } else {
          NULL
        },
        name = if (bottom_panel) "Genomic position" else NULL
      ) +
      scale_y_continuous(
        name   = y_label,
        limits = c(0, max(y_max, 0.01)),
        expand = expansion(mult=c(0, 0.05))
      ) +
      theme_classic(base_size=10) +
      theme(
        axis.title.y    = element_text(size=8, color="gray30", angle=90),
        axis.text.y     = element_text(size=7),
        legend.position = "right",
        legend.text     = element_text(size=8),
        legend.key.size = unit(0.4,"cm"),
        plot.margin     = margin(0, 6, if (bottom_panel) 4 else 0, 6)
      )
    
    if (!bottom_panel) {
      p <- p + theme(
        axis.line.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.x = element_blank()
      )
    } else {
      p <- p + theme(
        axis.text.x  = element_text(size=8),
        axis.title.x = element_text(size=9)
      )
    }
    p
  }
  
  # ── Panel 3: POLR2A CUT&TAG (WT only — single condition) ─────────────────
  p_pol2 <- make_signal_panel(
    sig_wt       = sig_pol2_wt,
    sig_ko       = NULL,          # single condition track
    fill_color   = "#7b2d8b",     # purple — distinct from RNA tracks
    ko_color     = NULL,
    y_max        = pol2_ymax,
    y_label      = "POLR2A\nCUT&TAG",
    bottom_panel = FALSE
  )
  
  # ── Panel 4: Total RNA-seq ────────────────────────────────────────────────
  p_total_rna <- make_signal_panel(
    sig_wt       = sig_total_wt,
    sig_ko       = sig_total_ko,
    fill_color   = "#2166ac",
    ko_color     = "#d73027",
    y_max        = total_ymax,
    y_label      = "Total RNA\n(CPM)",
    bottom_panel = FALSE
  )
  
  # ── Panel 5: EU-seq ───────────────────────────────────────────────────────
  p_eu_seq <- make_signal_panel(
    sig_wt       = sig_eu_wt,
    sig_ko       = sig_eu_ko,
    fill_color   = "#1b7837",
    ko_color     = "#d73027",
    y_max        = eu_ymax,
    y_label      = "EU-seq\n(subtracted)",
    bottom_panel = TRUE
  )
  
  # ══════════════════════════════════════════════════════════════════════════
  # Assemble with cowplot
  # ══════════════════════════════════════════════════════════════════════════
  leg_loops <- get_legend(p_loops     + theme(legend.position="right"))
  leg_genes <- get_legend(p_genes     + theme(legend.position="right"))
  leg_pol2  <- get_legend(p_pol2      + theme(legend.position="right"))
  leg_total <- get_legend(p_total_rna + theme(legend.position="right"))
  leg_eu    <- get_legend(p_eu_seq    + theme(legend.position="right"))
  
  all_legends <- plot_grid(
    leg_loops, leg_genes, leg_pol2, leg_total, leg_eu,
    ncol=1, align="v"
  )
  
  panels <- plot_grid(
    p_loops     + theme(legend.position="none"),
    p_genes     + theme(legend.position="none"),
    p_pol2      + theme(legend.position="none"),
    p_total_rna + theme(legend.position="none"),
    p_eu_seq    + theme(legend.position="none"),
    ncol        = 1,
    align       = "v",
    axis        = "lr",
    rel_heights = c(
      3,                           # loops
      gene_panel_range * 0.55,     # genes + cCREs
      1.2,                         # POLR2A CUT&TAG — slightly smaller,
      # single condition so less visual weight
      1.5,                         # total RNA
      1.5                          # EU-seq
    )
  )
  
  final <- plot_grid(
    panels, all_legends,
    ncol       = 2,
    rel_widths = c(1, 0.18)
  )
  
  print(final)
  invisible(final)
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Example calls
# ══════════════════════════════════════════════════════════════════════════════

## TITAN
plot_chiapet_region("chr2", 178530000, 178825000, loops,
                    txdb, org.Hs.eg.db, label="TTN locus")

plot_chiapet_region("chr5", 132550822, 132685508, loops,
                    txdb, org.Hs.eg.db)

plot_chiapet_region("chr8", 127698577, 127754815, loops,
                    txdb, org.Hs.eg.db, label="MYC locus")

plot_chiapet_region("chr5", 173890089, 174110000, loops,
                    txdb, org.Hs.eg.db, label="C5orf47 locus")

plot_chiapet_region("chr11", 102519233, 102630000, loops,
                    txdb, org.Hs.eg.db, label="MMP7 locus")
