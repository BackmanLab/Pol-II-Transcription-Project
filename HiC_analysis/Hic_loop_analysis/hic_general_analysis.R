##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Loop strength vs. DEG status (MA-style plot). Uses the same
# DS_GROUP path switch and loop-loading architecture as the APA
# pipeline and the Gained/Stable/Lost LFC pipeline.
#
# Also contains the anchor-type classification/raincloud logic
# extracted into HiC_analysis/Hic_loop_analysis/hic_loop_strength_
# by_anchor_type.R (Fig. 5G).
#
# Figure output: no confirmed figure for the MA-plot portion itself.
##################################################################

# Tol's "Muted" palette
sample_colors <- c(
  H3K27ac    = "#CC6677",  ## rose
  H3K4me3 = "#CC6677",  ## rose
  H3K27ac_ko = "#DDCC77",  ## sand
  H3K9me3    = "#332288",  ## indigo
  H3K9me3_ko = "#88CCEE",  ## cyan
  POL2       = "#44AA99"   ## teal
)

# Okabe-Ito Divergent Pair
# Okabe-Ito pair (#0072B2 / #D55E00)

library(strawr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(GenomicRanges)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readr)
library(tidyr)
library(ggdist)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ══════════════════════════════════════════════════════════════════════════════
# 0. Downsampling group switch — same convention as the APA / LFC pipelines
# ══════════════════════════════════════════════════════════════════════════════
DS_GROUP <- "ds_wt_pol2"   # ← "ds_wt_pol2" | "ds_all" | "unsampled"

BASE     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
LOOP_DIR <- file.path(BASE, "loops")

if (DS_GROUP == "ds_wt_pol2") {
  wt_loops   <- file.path(LOOP_DIR, "ds_wt_pol2/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  pol2_loops <- file.path(LOOP_DIR, "ds_wt_pol2/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  wt_hic     <- file.path(BASE, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic   <- file.path(BASE, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT   <- 1.0
  DEPTH_SCALE_POL2 <- 1.0
  MIN_OBS_WT       <- 0
  MIN_OBS_POL2     <- 0
  
} else if (DS_GROUP == "ds_all") {
  wt_loops   <- file.path(LOOP_DIR, "ds_all/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  pol2_loops <- file.path(LOOP_DIR, "ds_all/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  wt_hic     <- file.path(BASE, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic   <- file.path(BASE, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT   <- 1.0
  DEPTH_SCALE_POL2 <- 1.0
  MIN_OBS_WT       <- 0
  MIN_OBS_POL2     <- 0
  
} else if (DS_GROUP == "unsampled") {
  wt_loops   <- file.path(LOOP_DIR, "unsampled/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  pol2_loops <- file.path(LOOP_DIR, "unsampled/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  wt_hic     <- file.path(BASE, "mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic   <- file.path(BASE, "mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  wt_total   <- 639409109 + 377220354
  pol2_total <- 374163054 + 249677343
  DEPTH_SCALE_WT   <- 1.0
  DEPTH_SCALE_POL2 <- wt_total / pol2_total
  MIN_OBS_WT       <- 10
  MIN_OBS_POL2     <- 25
  
} else {
  stop("Unknown DS_GROUP: '", DS_GROUP, "'. Must be 'ds_wt_pol2', 'ds_all', or 'unsampled'.")
}
message("DS_GROUP: ", DS_GROUP)
message("  wt_loops:   ", wt_loops)
message("  pol2_loops: ", pol2_loops)
message("  wt_hic:     ", wt_hic)
message("  pol2_hic:   ", pol2_hic)
message("  depth scale WT=", round(DEPTH_SCALE_WT,3), "  pol2=", round(DEPTH_SCALE_POL2,3))

ENCODE_CCRE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/encode/wt_hct116.bed"

# ══════════════════════════════════════════════════════════════════════════════
# 1. Bézier arc helper
# ══════════════════════════════════════════════════════════════════════════════

bezier_arc <- function(x1, x2, height, n=100) {
  t  <- seq(0, 1, length.out=n)
  cx <- (x1 + x2) / 2
  cy <- height
  x  <- (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2
  y  <- (1-t)^2 * 0  + 2*(1-t)*t * cy + t^2 * 0
  data.frame(x=x, y=y)
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
       "' to file chromosomes: ", paste(head(file_chroms), collapse=", "))
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Load HiCCUPS loop calls — shared architecture with the APA / LFC pipelines
# ══════════════════════════════════════════════════════════════════════════════

load_hiccups_bedpe <- function(path) {
  all_lines   <- readLines(path)
  header_line <- all_lines[grepl("^#chr1", all_lines)][1]
  col_names   <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  data_lines  <- all_lines[!grepl("^#", all_lines)]
  if (length(data_lines) == 0) stop("No data lines found in ", path)
  
  df <- read.table(text=paste(data_lines, collapse="\n"),
                   header=FALSE, sep="\t",
                   stringsAsFactors=FALSE, fill=TRUE)
  
  if (ncol(df) == length(col_names)) colnames(df) <- col_names
  
  out <- data.frame(
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
  ) %>%
    filter(!is.na(start1), !is.na(start2))
  
  message("Loaded ", nrow(out), " loops from ", basename(path))
  message("Chromosome style in file: ",
          paste(head(unique(out$chr1)), collapse=", "))
  out
}

filter_loops <- function(df, min_observed=10, max_fdr=0.1) {
  df %>% filter(observed >= min_observed,
                fdrDonut <= max_fdr, fdrBL <= max_fdr,
                fdrH     <= max_fdr, fdrV  <= max_fdr)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3b. Real SCALE-normalized contact extraction via strawr
#     (the HiCCUPS "observed" column has NO normalization applied — this is
#     the same extract_loop_contacts() used in the Gained/Stable/Lost LFC
#     pipeline, kept identical across both scripts)
# ══════════════════════════════════════════════════════════════════════════════

STRAW_NORM       <- "SCALE"
STRAW_RESOLUTION <- 5000L   # bp — should match the resolution loops were called at

extract_loop_contacts <- function(loops_df, hic_file,
                                  norm      = STRAW_NORM,
                                  resolution = STRAW_RESOLUTION,
                                  matrix    = "observed") {
  snap    <- function(x) as.integer(floor(x / resolution) * resolution)
  results <- rep(NA_real_, nrow(loops_df))
  
  # Use the loop's own reported centroid (the actual called peak position —
  # WT's own centroid for shared/wt_only rows, Pol2's own centroid for
  # pol2_only rows, since that's the only bedpe that called a loop there at
  # all) rather than recomputing (start+end)/2 from the raw anchor span,
  # which can be tens of kb off the true peak due to HiCCUPS resolution-
  # merged anchor reporting. Falls back to the anchor midpoint only if this
  # loops_df doesn't carry centroid columns at all.
  has_centroids <- all(c("centroid1","centroid2") %in% colnames(loops_df)) &&
    !any(is.na(loops_df$centroid1)) && !any(is.na(loops_df$centroid2))
  if (has_centroids) {
    message("  Position source: centroid1/centroid2 (as reported by whichever ",
            "bedpe called this loop)")
    c1_pos <- snap(as.integer(loops_df$centroid1))
    c2_pos <- snap(as.integer(loops_df$centroid2))
  } else {
    message("  Position source: FALLBACK to (start+end)/2 — centroid1/centroid2 ",
            "missing or NA on this loops_df; positions may be off from the true peak")
    c1_pos <- snap(as.integer((loops_df$start1 + loops_df$end1) / 2))
    c2_pos <- snap(as.integer((loops_df$start2 + loops_df$end2) / 2))
  }
  
  chr1_clean <- gsub("^chr", "", loops_df$chr1)
  chr2_clean <- gsub("^chr", "", loops_df$chr2)
  chr_pairs  <- unique(paste(chr1_clean, chr2_clean, sep="|"))
  message("Extracting ", matrix, " from ", basename(hic_file),
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
      mat <- strawr::straw(norm, hic_file, loc1, loc2, "BP", resolution, matrix=matrix)
      if (nrow(mat) == 0) { results[idx] <- NA; next }
      lookup <- setNames(mat$counts,
                         paste(pmin(mat$x, mat$y), pmax(mat$x, mat$y), sep="|"))
      keys         <- paste(pmin(c1_pos[idx], c2_pos[idx]),
                            pmax(c1_pos[idx], c2_pos[idx]), sep="|")
      vals         <- lookup[keys]
      results[idx] <- ifelse(is.na(vals), NA, vals)
    }, error=function(e) message("  Warning: ", c1_chr, " — ", e$message))
  }
  n_hit <- sum(!is.na(results) & results > 0, na.rm=TRUE)
  message("  Signal: ", n_hit, " / ", nrow(loops_df),
          " (", round(100*n_hit/nrow(loops_df),1), "%)")
  results
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Load DESeq2 results
# ══════════════════════════════════════════════════════════════════════════════

DESEQ_FILE     <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1

deseq_res <- read_csv(DESEQ_FILE) %>%
  dplyr::mutate(
    gene_id   = sub("\\..*", "", gene_id),
    de_status = dplyr::case_when(
      is.na(padj)                                              ~ "Filtered",
      padj >= PADJ_THRESHOLD                                   ~ "Not DE",
      padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD ~ "Down",
      TRUE                                                     ~ "Not DE"
    ),
    de_status = factor(de_status,
                       levels = c("Up", "Not DE", "Down", "Filtered"))
  )

# ══════════════════════════════════════════════════════════════════════════════
# 5. Build GRanges from DESeq2 gene coordinates
# seq_name in your data has no "chr" prefix (e.g. "5" not "chr5")
# so we normalise to match the loop files
# ══════════════════════════════════════════════════════════════════════════════

deseq_gr <- deseq_res %>%
  filter(!is.na(gene_seq_start), !is.na(gene_seq_end),
         !is.na(seq_name)) %>%
  dplyr::mutate(
    seqname_clean = as.character(seq_name)
  )

gene_gr <- GRanges(
  seqnames    = deseq_gr$seqname_clean,
  ranges      = IRanges(deseq_gr$gene_seq_start, deseq_gr$gene_seq_end),
  symbol      = deseq_gr$symbol,
  log2fc_gene = deseq_gr$log2FoldChange,
  padj        = deseq_gr$padj,
  de_status   = as.character(deseq_gr$de_status)
)

# ══════════════════════════════════════════════════════════════════════════════
# 6. Build the UNION loop set across conditions (shared / wt_only / pol2_only)
#    — same architecture as the Gained/Stable/Lost LFC pipeline's
#    find_shared_loops() + union_loops, replacing the old intersection-only
#    get_matched_loops(). That old version queried findOverlaps(gi_deg, ...),
#    which meant a loop HiCCUPS called ONLY in WT (i.e. genuinely lost in
#    Pol2) could never generate a hit and was silently dropped before ever
#    reaching a plot — structurally incapable of showing loop loss. This
#    version keeps wt_only and pol2_only loops explicitly, tagged via a
#    `source` column, with real SCALE contacts extracted at those exact
#    positions in BOTH .hic files below (section 6b).
# ══════════════════════════════════════════════════════════════════════════════

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
  list(
    shared   = shared,
    a_only   = setdiff(seq_len(nrow(loops_a)), shared$a_idx),
    b_only   = setdiff(seq_len(nrow(loops_b)), shared$b_idx),
    n_shared = nrow(shared),
    n_a_only = nrow(loops_a) - dplyr::n_distinct(shared$a_idx),
    n_b_only = nrow(loops_b) - dplyr::n_distinct(shared$b_idx)
  )
}

build_union_loops <- function(loops_wt, loops_deg,
                              min_observed_wt  = MIN_OBS_WT,
                              min_observed_deg = MIN_OBS_POL2,
                              max_fdr          = 0.1,
                              slop             = 25000L) {
  
  wt_f  <- filter_loops(loops_wt,  min_observed=min_observed_wt,  max_fdr=max_fdr) %>%
    dplyr::mutate(loop_dist=abs(centroid2-centroid1))
  deg_f <- filter_loops(loops_deg, min_observed=min_observed_deg, max_fdr=max_fdr) %>%
    dplyr::mutate(loop_dist=abs(centroid2-centroid1))
  
  message("WT loops (filtered):   ", nrow(wt_f))
  message("Pol2 loops (filtered): ", nrow(deg_f))
  
  overlap <- find_shared_loops(wt_f, deg_f, slop=slop)
  message("Shared:    ", overlap$n_shared,
          "  (", round(100*overlap$n_shared/nrow(wt_f),1), "% of WT)")
  message("WT only:   ", overlap$n_a_only)
  message("pol2 only: ", overlap$n_b_only)
  
  shared_wt <- wt_f %>%
    dplyr::slice(overlap$shared$a_idx) %>%
    dplyr::mutate(source="shared")
  
  wt_only_loops <- wt_f %>%
    dplyr::slice(overlap$a_only) %>%
    dplyr::mutate(source="wt_only")
  
  pol2_only_loops <- deg_f %>%
    dplyr::slice(overlap$b_only) %>%
    dplyr::mutate(source="pol2_only")
  
  union_loops <- bind_rows(shared_wt, wt_only_loops, pol2_only_loops) %>%
    dplyr::mutate(loop_idx=row_number())
  
  message("Union loop set: ", nrow(union_loops),
          "  (shared=", sum(union_loops$source=="shared"),
          "  wt_only=", sum(union_loops$source=="wt_only"),
          "  pol2_only=", sum(union_loops$source=="pol2_only"), ")")
  
  union_loops
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. Overlap loop anchors with DE genes
# A loop is linked to a gene if EITHER anchor overlaps the gene body
# Returns one row per loop-gene pair
# ══════════════════════════════════════════════════════════════════════════════

link_loops_to_genes <- function(matched_loops, gene_gr) {
  
  anchor1_gr <- GRanges(
    seqnames = matched_loops$chr1,
    ranges   = IRanges(matched_loops$start1, matched_loops$end1)
  )
  anchor2_gr <- GRanges(
    seqnames = matched_loops$chr2,
    ranges   = IRanges(matched_loops$start2, matched_loops$end2)
  )
  
  hits1 <- findOverlaps(anchor1_gr, gene_gr)
  hits2 <- findOverlaps(anchor2_gr, gene_gr)
  
  pairs1 <- data.frame(loop_idx=queryHits(hits1), gene_idx=subjectHits(hits1), anchor="anchor1")
  pairs2 <- data.frame(loop_idx=queryHits(hits2), gene_idx=subjectHits(hits2), anchor="anchor2")
  
  all_pairs <- bind_rows(pairs1, pairs2) %>%
    distinct(loop_idx, gene_idx, .keep_all=TRUE)
  
  if (nrow(all_pairs) == 0) {
    message("No loop-gene overlaps found")
    return(NULL)
  }
  
  result <- all_pairs %>%
    dplyr::mutate(
      log2fc_loop = matched_loops$log2fc_loop[loop_idx],
      loop_class  = matched_loops$loop_class[loop_idx],
      loop_dist   = matched_loops$loop_dist[loop_idx],
      obs_wt      = matched_loops$obs_wt[loop_idx],
      obs_deg     = matched_loops$obs_deg[loop_idx],
      symbol      = gene_gr$symbol[gene_idx],
      log2fc_gene = gene_gr$log2fc_gene[gene_idx],
      padj        = gene_gr$padj[gene_idx],
      de_status   = gene_gr$de_status[gene_idx]
    ) %>%
    filter(de_status %in% c("Up", "Down"))
  
  message(nrow(result), " loop-gene pairs found (",
          sum(result$de_status=="Up"), " upregulated, ",
          sum(result$de_status=="Down"), " downregulated genes)")
  result
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Plot: loop log2FC vs gene log2FC, split by DE direction
# ══════════════════════════════════════════════════════════════════════════════

plot_loop_gene_integration <- function(loop_gene_pairs) {
  
  if (is.null(loop_gene_pairs) || nrow(loop_gene_pairs) == 0) {
    stop("No loop-gene pairs to plot")
  }
  
  color_map <- c(
    "lost"      = "#2166ac",
    "gained"    = "#d6604d",
    "unchanged" = "gray70"
  )
  
  make_panel <- function(df, title) {
    n_lost      <- sum(df$loop_class=="lost")
    n_gained    <- sum(df$loop_class=="gained")
    n_unchanged <- sum(df$loop_class=="unchanged")
    
    ggplot(df, aes(x=log2fc_loop, y=log2fc_gene, color=loop_class)) +
      geom_point(alpha=0.6, size=1.8) +
      geom_vline(xintercept=c(-0.5, 0.5),
                 linetype="dashed", color="gray60", linewidth=0.4) +
      geom_hline(yintercept=0,
                 linetype="solid",  color="gray40", linewidth=0.4) +
      geom_hline(yintercept=c(-LFC_THRESHOLD, LFC_THRESHOLD),
                 linetype="dashed", color="gray60", linewidth=0.4) +
      ggrepel::geom_text_repel(
        data    = df %>%
          arrange(desc(abs(log2fc_gene))) %>%
          slice_head(n=15),
        aes(label=symbol),
        size    = 2.8,
        color   = "black",
        max.overlaps = 20,
        segment.color = "gray70"
      ) +
      scale_color_manual(
        values = color_map,
        labels = c(
          "lost"      = paste0("Lost (n=",      n_lost,      ")"),
          "gained"    = paste0("Gained (n=",    n_gained,    ")"),
          "unchanged" = paste0("Unchanged (n=", n_unchanged, ")")
        ),
        name = "Loop change"
      ) +
      labs(
        title    = title,
        subtitle = paste0(DS_GROUP, "  |  n = ", nrow(df), " loop-gene pairs"),
        x        = "Loop log2FC (degraded / wildtype)",
        y        = "Gene log2FC (DESeq2)"
      ) +
      theme_classic(base_size=12) +
      theme(
        plot.title    = element_text(face="bold", size=13),
        plot.subtitle = element_text(size=9, color="gray40"),
        legend.position = "bottom"
      )
  }
  
  up_pairs   <- loop_gene_pairs %>% filter(de_status=="Up")
  down_pairs <- loop_gene_pairs %>% filter(de_status=="Down")
  
  p_up   <- make_panel(up_pairs,   "Upregulated genes")
  p_down <- make_panel(down_pairs, "Downregulated genes")
  
  composite <- (p_up | p_down) +
    plot_annotation(
      title = "Loop contact changes at DE genes after POLR2A degradation",
      theme = theme(plot.title=element_text(face="bold", size=14))
    )
  
  print(composite)
  invisible(list(up=p_up, down=p_down, data=loop_gene_pairs))
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Run the full integration
# ══════════════════════════════════════════════════════════════════════════════

loops_wt_df  <- load_hiccups_bedpe(wt_loops)
loops_deg_df <- load_hiccups_bedpe(pol2_loops)

PSEUDOCOUNT <- 1  # added before log2 to avoid log(0) — matches LFC pipeline convention

matched <- build_union_loops(loops_wt_df, loops_deg_df,
                             min_observed_wt  = MIN_OBS_WT,
                             min_observed_deg = MIN_OBS_POL2,
                             max_fdr          = 0.1,
                             slop             = 25000L)

message("\n── Extracting SCALE-normalized contacts at ALL union positions ──────")
message("(this now includes wt_only and pol2_only loops, not just shared ones —")
message(" so a WT-called loop with no Pol2 signal will show near-zero scale_deg")
message(" instead of being dropped before it ever reaches a plot)")
scale_wt_raw  <- extract_loop_contacts(matched, wt_hic,   norm=STRAW_NORM, resolution=STRAW_RESOLUTION)
scale_deg_raw <- extract_loop_contacts(matched, pol2_hic, norm=STRAW_NORM, resolution=STRAW_RESOLUTION)

matched <- matched %>%
  dplyr::mutate(
    obs_wt      = tidyr::replace_na(scale_wt_raw,  0) * DEPTH_SCALE_WT,
    obs_deg     = tidyr::replace_na(scale_deg_raw, 0) * DEPTH_SCALE_POL2,
    scale_wt    = obs_wt,   # alias — the raincloud panel functions read scale_wt/scale_deg
    scale_deg   = obs_deg,
    log2fc_loop = log2((obs_deg + PSEUDOCOUNT) / (obs_wt + PSEUDOCOUNT)),
    loop_class  = dplyr::case_when(
      log2fc_loop >  0.5 ~ "gained",
      log2fc_loop < -0.5 ~ "lost",
      TRUE               ~ "unchanged"
    )
  )

message("Union loops: ", nrow(matched), "  |  DS_GROUP=", DS_GROUP)
message("By source:"); print(table(matched$source))
message("By class:");  print(table(matched$loop_class))
message("Source x class cross-tab (this is the sanity check — wt_only should")
message("skew heavily toward 'lost', pol2_only heavily toward 'gained'):")
print(table(matched$source, matched$loop_class))

loop_gene_pairs <- link_loops_to_genes(matched, gene_gr)

result <- plot_loop_gene_integration(loop_gene_pairs)

# Sanity checks
head(deseq_res)
colnames(deseq_res)


#########################################################################################
#  Updated ChIA-PET anchor classifier — now includes cCRE category
#  (unaffected by the DS_GROUP change — uses its own `loops` object)
#########################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# Updated ChIA-PET anchor classifier — matches Hi-C classification logic
# Each anchor checked independently against exon, intron, and cCRE
# An anchor landing on a cCRE inside an intron is classified as cCRE
# ══════════════════════════════════════════════════════════════════════════════

classify_loop_anchors <- function(chiapet_loops, exon_gr, intron_gr,
                                  ccre_gr=NULL, anchor_window=STRAW_RESOLUTION) {
  
  chiapet_loops <- chiapet_loops %>%
    dplyr::mutate(loop_idx=row_number())
  
  fix_chr <- function(gr, reference_gr) {
    ref_chroms  <- as.character(seqnames(reference_gr))
    has_chr_ref <- any(grepl("^chr", ref_chroms))
    gr_chroms   <- as.character(seqnames(gr))
    has_chr_gr  <- any(grepl("^chr", gr_chroms))
    if (has_chr_ref && !has_chr_gr) {
      seqlevels(gr) <- paste0("chr", seqlevels(gr))
    } else if (!has_chr_ref && has_chr_gr) {
      seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    }
    gr
  }
  
  # Use a tight, resolution-sized window centered on each anchor's own
  # midpoint — NOT the raw start1:end1 / start2:end2 span, which (as with
  # the Hi-C HiCCUPS classifier) can be much wider than the loop's actual
  # called position and would otherwise act as an unintended large
  # tolerance against small features like exons. Falls back to the anchor
  # midpoint (start+end)/2 if this loop source doesn't carry its own
  # centroid1/centroid2 columns (ChIA-PET loop callers vary in this).
  half_w <- anchor_window / 2
  center1 <- if ("centroid1" %in% colnames(chiapet_loops))
    chiapet_loops$centroid1 else (chiapet_loops$start1 + chiapet_loops$end1) / 2
  center2 <- if ("centroid2" %in% colnames(chiapet_loops))
    chiapet_loops$centroid2 else (chiapet_loops$start2 + chiapet_loops$end2) / 2
  message("Anchor classification window: center \u00b1", half_w,
          "bp (", anchor_window, "bp total) \u2014 NOT the raw anchor span")
  
  anchor1_gr <- GRanges(
    seqnames = chiapet_loops$chr1,
    ranges   = IRanges(center1 - half_w, center1 + half_w)
  )
  anchor2_gr <- GRanges(
    seqnames = chiapet_loops$chr2,
    ranges   = IRanges(center2 - half_w, center2 + half_w)
  )
  
  anchor1_gr <- fix_chr(anchor1_gr, exon_gr)
  anchor2_gr <- fix_chr(anchor2_gr, exon_gr)
  
  a1_exon   <- overlapsAny(anchor1_gr, exon_gr)
  a1_intron <- overlapsAny(anchor1_gr, intron_gr)
  a2_exon   <- overlapsAny(anchor2_gr, exon_gr)
  a2_intron <- overlapsAny(anchor2_gr, intron_gr)
  
  if (!is.null(ccre_gr)) {
    ccre_gr_fixed <- fix_chr(ccre_gr, exon_gr)
    a1_ccre <- overlapsAny(anchor1_gr, ccre_gr_fixed)
    a2_ccre <- overlapsAny(anchor2_gr, ccre_gr_fixed)
  } else {
    a1_ccre <- rep(FALSE, length(anchor1_gr))
    a2_ccre <- rep(FALSE, length(anchor2_gr))
  }
  
  classify_anchor <- function(is_exon, is_intron, is_ccre) {
    dplyr::case_when(
      is_exon  & is_ccre                  ~ "exon+cCRE",
      is_exon  & !is_ccre                 ~ "exon",
      !is_exon & is_ccre                  ~ "cCRE",
      !is_exon & !is_ccre & is_intron     ~ "intron",
      TRUE                                ~ "intergenic"
    )
  }
  
  anchor1_class <- classify_anchor(a1_exon, a1_intron, a1_ccre)
  anchor2_class <- classify_anchor(a2_exon, a2_intron, a2_ccre)
  
  message("Anchor1 class breakdown:")
  print(table(anchor1_class))
  message("Anchor2 class breakdown:")
  print(table(anchor2_class))
  
  loop_type <- dplyr::case_when(
    grepl("exon", anchor1_class) & grepl("exon", anchor2_class)
    ~ "Exon-Exon",
    (grepl("exon", anchor1_class) & anchor2_class == "cCRE") |
      (anchor1_class == "cCRE" & grepl("exon", anchor2_class))
    ~ "Exon-cCRE",
    (grepl("exon", anchor1_class) & anchor2_class == "intron") |
      (anchor1_class == "intron" & grepl("exon", anchor2_class))
    ~ "Exon-Intron",
    anchor1_class == "cCRE" & anchor2_class == "cCRE"
    ~ "cCRE-cCRE",
    anchor1_class == "intron" & anchor2_class == "intron"
    ~ "Intron-Intron",
    TRUE ~ "Other"
  )
  
  chiapet_loops %>%
    dplyr::mutate(
      anchor1_class    = anchor1_class,
      anchor2_class    = anchor2_class,
      loop_anchor_type = factor(
        loop_type,
        levels = c("Exon-Exon",
                   "Exon-cCRE",
                   "Exon-Intron",
                   "cCRE-cCRE",
                   "Intron-Intron",
                   "Other")
      )
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# Plot function — unchanged
# ══════════════════════════════════════════════════════════════════════════════

plot_loop_strength_raincloud <- function(classified_loops) {
  
  plot_df <- as.data.frame(classified_loops) %>%
    filter(loop_anchor_type %in% c("Exon-Exon",
                                   "Exon-cCRE",
                                   "Exon-Intron",
                                   "cCRE-cCRE",
                                   "Intron-Intron")) %>%
    droplevels() %>%
    dplyr::mutate(
      log10_score_capped = pmin(log10(score + 1), log10(101)),
      loop_anchor_type   = factor(loop_anchor_type,
                                  levels=c("Exon-Exon",
                                           "Exon-cCRE",
                                           "Exon-Intron",
                                           "cCRE-cCRE",
                                           "Intron-Intron"))
    )
  
  counts <- plot_df %>%
    dplyr::count(loop_anchor_type) %>%
    dplyr::mutate(
      loop_anchor_type = factor(loop_anchor_type,
                                levels=c("Exon-Exon",
                                         "Exon-cCRE",
                                         "Exon-Intron",
                                         "cCRE-cCRE",
                                         "Intron-Intron")),
      label = paste0(loop_anchor_type, "\n(n=", n, ")")
    ) %>%
    dplyr::arrange(loop_anchor_type)
  
  plot_df <- plot_df %>%
    left_join(counts %>% dplyr::select(loop_anchor_type, label),
              by="loop_anchor_type") %>%
    dplyr::mutate(label=factor(label, levels=counts$label))
  
  colors <- c(
    "Exon-Exon"     = "#d6604d",
    "Exon-cCRE"     = "#984ea3",
    "Exon-Intron"   = "#4DAF4A",
    "cCRE-cCRE"     = "#ff7f00",
    "Intron-Intron" = "#2166ac"
  )
  
  ggplot(plot_df, aes(x=label, y=log10_score_capped,
                      fill=loop_anchor_type,
                      color=loop_anchor_type)) +
    
    ggdist::stat_halfeye(
      aes(thickness = after_stat(pdf * n)),
      side         = "right",
      scale        = 0.4,
      width        = 0.6,
      .width       = 0,
      point_colour = NA,
      alpha        = 0.75,
      normalize    = "groups",
      position     = position_nudge(x=0.05)
    ) +
    
    geom_boxplot(
      width         = 0.12,
      outlier.shape = NA,
      color         = "gray20",
      alpha         = 0.9,
      linewidth     = 0.6,
      position      = position_nudge(x=-0.07)
    ) +
    
    scale_fill_manual(values=colors,  guide="none") +
    scale_color_manual(values=colors, guide="none") +
    
    scale_y_continuous(
      name   = "Loop strength (PET count)",
      limits = c(NA, 2.15),
      breaks = c(log10(2), log10(5), 1, log10(20), log10(50), 2),
      labels = c("2", "5", "10", "20", "50", "100"),
      expand = expansion(mult=c(0.02, 0))
    ) +
    
    scale_x_discrete(expand=expansion(add=c(0.3, 0.5))) +
    
    labs(
      title = "POLR2A ChIA-PET loop strength by anchor type",
      x     = NULL
    ) +
    
    theme_classic(base_size=13) +
    theme(
      plot.title   = element_text(face="bold", size=14),
      axis.text.x  = element_text(size=10, angle=30, hjust=1, vjust=1),
      plot.margin  = margin(8, 16, 8, 8),
      aspect.ratio = 1.4
    )
}

# ── Run ───────────────────────────────────────────────────────────────────────
# NOTE: this section expects `loops`, `exon_gr`, `intron_gr`, `ccre_gr` to
# already exist from your ChIA-PET loading step — unrelated to the Hi-C
# DS_GROUP switch above.
#
# chiapet_classified <- classify_loop_anchors(loops, exon_gr, intron_gr,
#                                             ccre_gr=ccre_gr)
#
# cat("\nLoop anchor type breakdown:\n")
# print(table(chiapet_classified$loop_anchor_type))
#
# plot_loop_strength_raincloud(chiapet_classified)


#########################################################################################
#  HIC loop summarization (intronic, exonic, etc) — with DEGs
#########################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# Load cCREs for anchor classification
# ══════════════════════════════════════════════════════════════════════════════

load_ccre_gr <- function(ccre_path) {
  df <- read.table(ccre_path, header=FALSE, sep="\t",
                   stringsAsFactors=FALSE) %>%
    dplyr::select(chr=V1, start=V2, end=V3, id=V4, ccre_type=V10) %>%
    filter(ccre_type %in% c("PLS", "pELS", "dELS", "CA-CTCF"))
  
  message("Loaded ", nrow(df), " cCRE elements:")
  print(table(df$ccre_type))
  
  GRanges(
    seqnames = df$chr,
    ranges   = IRanges(df$start, df$end),
    ccre_type = df$ccre_type
  )
}

# Load once and pass to classify function
if (!exists("ccre_gr")) {
  ccre_gr <- load_ccre_gr(ENCODE_CCRE)
}

# ══════════════════════════════════════════════════════════════════════════════
# Updated anchor classifier — now includes cCRE category
# ══════════════════════════════════════════════════════════════════════════════

classify_hic_loop_anchors <- function(matched_loops, exon_gr, intron_gr,
                                      ccre_gr=NULL, anchor_window=STRAW_RESOLUTION) {
  
  fix_chr <- function(gr, reference_gr) {
    ref_chroms  <- as.character(seqnames(reference_gr))
    has_chr_ref <- any(grepl("^chr", ref_chroms))
    gr_chroms   <- as.character(seqnames(gr))
    has_chr_gr  <- any(grepl("^chr", gr_chroms))
    if (has_chr_ref && !has_chr_gr) {
      seqlevels(gr) <- paste0("chr", seqlevels(gr))
    } else if (!has_chr_ref && has_chr_gr) {
      seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    }
    gr
  }
  
  classify_anchor <- function(gr) {
    gr        <- fix_chr(gr, exon_gr)
    is_exon   <- overlapsAny(gr, exon_gr)
    is_intron <- overlapsAny(gr, intron_gr)
    
    is_ccre <- if (!is.null(ccre_gr)) {
      ccre_gr_fixed <- fix_chr(ccre_gr, gr)
      overlapsAny(gr, ccre_gr_fixed)
    } else {
      rep(FALSE, length(gr))
    }
    
    dplyr::case_when(
      is_exon               ~ "exon",
      is_intron & !is_exon  ~ "intron",
      is_ccre               ~ "cCRE",
      TRUE                  ~ "intergenic"
    )
  }
  
  # Tight, resolution-sized window centered on each anchor's own centroid —
  # NOT the raw start1:end1 / start2:end2 span (see classify_hic_loop_anchors_wt
  # for the full rationale: merged_loops.bedpe anchor spans can be much wider
  # than a single bin due to HiCCUPS resolution-merging, which would otherwise
  # act as an unintended large tolerance against small features like exons).
  half_w <- anchor_window / 2
  message("Anchor classification window: centroid \u00b1", half_w,
          "bp (", anchor_window, "bp total) \u2014 NOT the raw bedpe anchor span")
  
  anchor1_gr <- GRanges(
    seqnames = matched_loops$chr1,
    ranges   = IRanges(matched_loops$centroid1 - half_w, matched_loops$centroid1 + half_w)
  )
  anchor2_gr <- GRanges(
    seqnames = matched_loops$chr2,
    ranges   = IRanges(matched_loops$centroid2 - half_w, matched_loops$centroid2 + half_w)
  )
  
  matched_loops %>%
    dplyr::mutate(
      anchor1_class    = classify_anchor(anchor1_gr),
      anchor2_class    = classify_anchor(anchor2_gr),
      loop_anchor_type = dplyr::case_when(
        anchor1_class == "exon"   & anchor2_class == "exon"   ~ "Exon-Exon",
        anchor1_class == "intron" & anchor2_class == "intron" ~ "Intron-Intron",
        (anchor1_class == "exon"  & anchor2_class == "intron") |
          (anchor1_class == "intron" & anchor2_class == "exon") ~ "Exon-Intron",
        (anchor1_class == "exon"  & anchor2_class == "cCRE") |
          (anchor1_class == "cCRE"  & anchor2_class == "exon") ~ "Exon-cCRE",
        TRUE ~ "Other"
      ),
      loop_anchor_type = factor(loop_anchor_type,
                                levels=c("Exon-Exon",
                                         "Exon-Intron",
                                         "Exon-cCRE",
                                         "Intron-Intron"))
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# Updated data prep — now includes Exon-cCRE group
# ══════════════════════════════════════════════════════════════════════════════

prepare_hic_raincloud_data <- function(classified_matched) {
  classified_matched %>%
    filter(loop_anchor_type %in% c("Exon-Exon",
                                   "Exon-Intron",
                                   "Exon-cCRE",
                                   "Intron-Intron")) %>%
    dplyr::select(loop_anchor_type, obs_wt, obs_deg) %>%
    filter(!is.na(obs_wt),  obs_wt  > 0,
           !is.na(obs_deg), obs_deg > 0) %>%
    tidyr::pivot_longer(
      cols      = c(obs_wt, obs_deg),
      names_to  = "condition",
      values_to = "observed"
    ) %>%
    dplyr::mutate(
      condition = dplyr::case_when(
        condition == "obs_wt"  ~ "Wildtype",
        condition == "obs_deg" ~ "Degraded"
      ),
      condition      = factor(condition, levels=c("Wildtype", "Degraded")),
      log10_observed = pmin(log10(observed + 1), log10(1001)),
      x_pos = dplyr::case_when(
        loop_anchor_type=="Exon-Exon"     & condition=="Wildtype"  ~ 1.0,
        loop_anchor_type=="Exon-Exon"     & condition=="Degraded"  ~ 1.4,
        loop_anchor_type=="Exon-Intron"   & condition=="Wildtype"  ~ 2.4,
        loop_anchor_type=="Exon-Intron"   & condition=="Degraded"  ~ 2.8,
        loop_anchor_type=="Exon-cCRE"     & condition=="Wildtype"  ~ 3.8,
        loop_anchor_type=="Exon-cCRE"     & condition=="Degraded"  ~ 4.2,
        loop_anchor_type=="Intron-Intron" & condition=="Wildtype"  ~ 5.2,
        loop_anchor_type=="Intron-Intron" & condition=="Degraded"  ~ 5.6
      )
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# Updated plot function — four anchor type groups
# ══════════════════════════════════════════════════════════════════════════════

plot_hic_loop_strength_raincloud <- function(classified_matched,
                                             exon_gr, intron_gr,
                                             ccre_gr=NULL) {
  
  classified <- classify_hic_loop_anchors(classified_matched,
                                          exon_gr, intron_gr,
                                          ccre_gr=ccre_gr)
  
  message("Hi-C loop anchor type breakdown:")
  print(table(classified$loop_anchor_type))
  
  plot_df <- prepare_hic_raincloud_data(classified)
  
  anchor_colors <- c(
    "Exon-Exon"     = "#CC6677",
    "Exon-Intron"   = "#44AA99",
    "Exon-cCRE"     = "#984ea3",
    "Intron-Intron" = "#2166ac"
  )
  
  get_pval <- function(anchor) {
    wt  <- plot_df %>%
      filter(loop_anchor_type==anchor, condition=="Wildtype") %>%
      pull(log10_observed)
    deg <- plot_df %>%
      filter(loop_anchor_type==anchor, condition=="Degraded") %>%
      pull(log10_observed)
    if (length(wt) != length(deg) || length(wt) < 3) return("ns")
    p <- wilcox.test(wt, deg, paired=TRUE)$p.value
    dplyr::case_when(
      p < 0.0001 ~ "****",
      p < 0.001  ~ "***",
      p < 0.01   ~ "**",
      p < 0.05   ~ "*",
      TRUE       ~ "ns"
    )
  }
  
  pval_ee  <- get_pval("Exon-Exon")
  pval_ei  <- get_pval("Exon-Intron")
  pval_ec  <- get_pval("Exon-cCRE")
  pval_ii  <- get_pval("Intron-Intron")
  
  n_per_group <- plot_df %>% dplyr::count(x_pos)
  
  counts <- plot_df %>%
    dplyr::distinct(x_pos, condition, loop_anchor_type) %>%
    dplyr::left_join(n_per_group, by="x_pos") %>%
    dplyr::mutate(x_label = paste0(condition, "\n(n=", n, ")"))
  
  ggplot(plot_df, aes(x=x_pos, y=log10_observed,
                      fill=condition, color=condition)) +
    
    annotate("rect", xmin=0.75, xmax=1.65,
             ymin=-Inf, ymax=Inf,
             fill=anchor_colors["Exon-Exon"], alpha=0.1) +
    annotate("rect", xmin=2.15, xmax=3.05,
             ymin=-Inf, ymax=Inf,
             fill=anchor_colors["Exon-Intron"], alpha=0.1) +
    annotate("rect", xmin=3.55, xmax=4.45,
             ymin=-Inf, ymax=Inf,
             fill=anchor_colors["Exon-cCRE"], alpha=0.1) +
    annotate("rect", xmin=4.95, xmax=5.85,
             ymin=-Inf, ymax=Inf,
             fill=anchor_colors["Intron-Intron"], alpha=0.1) +
    
    annotate("text", x=1.2,  y=3.25, label="Exon-Exon",
             size=3.5, color=anchor_colors["Exon-Exon"], fontface="bold") +
    annotate("text", x=2.6,  y=3.25, label="Exon-Intron",
             size=3.5, color=anchor_colors["Exon-Intron"], fontface="bold") +
    annotate("text", x=4.0,  y=3.25, label="Exon-cCRE",
             size=3.5, color=anchor_colors["Exon-cCRE"], fontface="bold") +
    annotate("text", x=5.4,  y=3.25, label="Intron-Intron",
             size=3.5, color=anchor_colors["Intron-Intron"], fontface="bold") +
    
    annotate("text", x=1.2,  y=3.10, label=pval_ee, size=4, color="gray30") +
    annotate("text", x=2.6,  y=3.10, label=pval_ei, size=4, color="gray30") +
    annotate("text", x=4.0,  y=3.10, label=pval_ec, size=4, color="gray30") +
    annotate("text", x=5.4,  y=3.10, label=pval_ii, size=4, color="gray30") +
    
    ggdist::stat_halfeye(
      aes(thickness = after_stat(pdf * n),
          group     = interaction(loop_anchor_type, condition)),
      side         = "right",
      scale        = 0.28,
      width        = 0.45,
      .width       = 0,
      point_colour = NA,
      alpha        = 0.7,
      normalize    = "groups",
      position     = position_nudge(x=0.03)
    ) +
    
    geom_boxplot(
      aes(group = interaction(loop_anchor_type, condition)),
      width         = 0.08,
      outlier.shape = NA,
      color         = "gray20",
      alpha         = 0.9,
      linewidth     = 0.5,
      position      = position_nudge(x=-0.05)
    ) +
    
    scale_fill_manual(
      values = c("Wildtype"="#0072B2", "Degraded"="#D55E00"),
      name   = NULL
    ) +
    scale_color_manual(
      values = c("Wildtype"="#0072B2", "Degraded"="#D55E00"),
      name   = NULL,
      guide  = "none"
    ) +
    
    scale_x_continuous(
      breaks = counts$x_pos,
      labels = counts$x_label,
      limits = c(0.55, 6.4),
      expand = c(0, 0),
      oob    = scales::squish
    ) +
    
    scale_y_continuous(
      name   = "Observed contacts (KR normalised)",
      limits = c(NA, 3.35),
      breaks = c(0, 1, 2, 3),
      labels = c("1", "10", "100", "1000"),
      expand = expansion(mult=c(0.02, 0))
    ) +
    
    labs(
      title    = "Hi-C loop contact strength by anchor type",
      subtitle = paste0(DS_GROUP, "  |  Wildtype vs POLR2A degraded  |  paired Wilcoxon test"),
      x        = NULL
    ) +
    
    theme_classic(base_size=13) +
    theme(
      plot.title      = element_text(face="bold", size=14),
      plot.subtitle   = element_text(size=9, color="gray40"),
      axis.text.x     = element_text(size=9, angle=37, hjust=1, vjust=1),
      legend.position = "bottom"
    )
}

# ── Run ───────────────────────────────────────────────────────────────────────
if (!exists("exon_gr") || !exists("intron_gr")) {
  message("Building exon/intron ranges from txdb...")
  all_exons   <- reduce(exons(txdb))
  all_genes   <- genes(txdb)
  all_introns <- GenomicRanges::setdiff(all_genes, all_exons)
  all_introns <- all_introns[!overlapsAny(all_introns, all_exons)]
  exon_gr   <- all_exons
  intron_gr <- all_introns
  message("Exon ranges:   ", length(exon_gr))
  message("Intron ranges: ", length(intron_gr))
}

plot_hic_loop_strength_raincloud(matched, exon_gr, intron_gr, ccre_gr=ccre_gr)


#########################################################################################
#  HIC loop summarization (intronic, exonic, etc) — no DEGs, WT only
#########################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 1. Build exon and intron GRanges (reuse from above if already built)
# ══════════════════════════════════════════════════════════════════════════════

build_exon_intron_ranges <- function(txdb) {
  message("Building exon ranges...")
  all_exons <- reduce(exons(txdb))
  
  message("Building intron ranges...")
  all_genes   <- genes(txdb)
  all_introns <- GenomicRanges::setdiff(all_genes, all_exons)
  all_introns <- all_introns[!overlapsAny(all_introns, all_exons)]
  
  message("Exon ranges:   ", length(all_exons))
  message("Intron ranges: ", length(all_introns))
  
  list(exons=all_exons, introns=all_introns)
}

if (!exists("exon_gr") || !exists("intron_gr")) {
  ranges    <- build_exon_intron_ranges(txdb)
  exon_gr   <- ranges$exons
  intron_gr <- ranges$introns
}
if (!exists("ccre_gr")) {
  ccre_gr <- load_ccre_gr(ENCODE_CCRE)
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Classify Hi-C loop anchors (WT-only variant, kept distinct from the
#    DEG-integrated classify_hic_loop_anchors() defined above so this
#    section can be run independently)
#    Each anchor is checked independently against exon, intron, and cCRE
#    An anchor landing on a cCRE inside an intron is classified as cCRE
# ══════════════════════════════════════════════════════════════════════════════

classify_hic_loop_anchors_wt <- function(hic_loops, exon_gr, intron_gr,
                                         ccre_gr=NULL,
                                         anchor_window = STRAW_RESOLUTION) {
  
  hic_loops <- hic_loops %>%
    dplyr::mutate(loop_idx=row_number())
  
  fix_chr <- function(gr, reference_gr) {
    ref_chroms  <- as.character(seqnames(reference_gr))
    has_chr_ref <- any(grepl("^chr", ref_chroms))
    gr_chroms   <- as.character(seqnames(gr))
    has_chr_gr  <- any(grepl("^chr", gr_chroms))
    if (has_chr_ref && !has_chr_gr) {
      seqlevels(gr) <- paste0("chr", seqlevels(gr))
    } else if (!has_chr_ref && has_chr_gr) {
      seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    }
    gr
  }
  
  # Use a tight, resolution-sized window centered on each anchor's own
  # centroid — NOT the raw start1:end1 / start2:end2 span reported in
  # merged_loops.bedpe. That raw span can be much wider than a single bin,
  # since HiCCUPS merged loop calls often report the union footprint across
  # every resolution the loop was called at (5/10/25kb), which would
  # otherwise act as a large, invisible tolerance when checking overlap
  # against small features like exons (often only ~100-300bp).
  half_w <- anchor_window / 2
  message("Anchor classification window: centroid \u00b1", half_w,
          "bp (", anchor_window, "bp total) \u2014 NOT the raw bedpe anchor span")
  
  anchor1_gr <- GRanges(
    seqnames = hic_loops$chr1,
    ranges   = IRanges(hic_loops$centroid1 - half_w, hic_loops$centroid1 + half_w)
  )
  anchor2_gr <- GRanges(
    seqnames = hic_loops$chr2,
    ranges   = IRanges(hic_loops$centroid2 - half_w, hic_loops$centroid2 + half_w)
  )
  
  anchor1_gr <- fix_chr(anchor1_gr, exon_gr)
  anchor2_gr <- fix_chr(anchor2_gr, exon_gr)
  
  a1_exon   <- overlapsAny(anchor1_gr, exon_gr)
  a1_intron <- overlapsAny(anchor1_gr, intron_gr)
  a2_exon   <- overlapsAny(anchor2_gr, exon_gr)
  a2_intron <- overlapsAny(anchor2_gr, intron_gr)
  
  if (!is.null(ccre_gr)) {
    ccre_gr_fixed <- fix_chr(ccre_gr, exon_gr)
    a1_ccre <- overlapsAny(anchor1_gr, ccre_gr_fixed)
    a2_ccre <- overlapsAny(anchor2_gr, ccre_gr_fixed)
  } else {
    a1_ccre <- rep(FALSE, length(anchor1_gr))
    a2_ccre <- rep(FALSE, length(anchor2_gr))
  }
  
  classify_anchor <- function(is_exon, is_intron, is_ccre) {
    dplyr::case_when(
      is_exon  & is_ccre                  ~ "exon+cCRE",
      is_exon  & !is_ccre                 ~ "exon",
      !is_exon & is_ccre                  ~ "cCRE",
      !is_exon & !is_ccre & is_intron     ~ "intron",
      TRUE                                ~ "intergenic"
    )
  }
  
  anchor1_class <- classify_anchor(a1_exon, a1_intron, a1_ccre)
  anchor2_class <- classify_anchor(a2_exon, a2_intron, a2_ccre)
  
  message("Anchor1 class breakdown:")
  print(table(anchor1_class))
  message("Anchor2 class breakdown:")
  print(table(anchor2_class))
  message("Raw cCRE overlaps — anchor1: ", sum(a1_ccre),
          "  anchor2: ", sum(a2_ccre),
          "  both: ", sum(a1_ccre & a2_ccre))
  message("Full anchor class cross-tabulation:")
  print(table(anchor1_class, anchor2_class))
  
  loop_type <- dplyr::case_when(
    grepl("exon", anchor1_class) & grepl("exon", anchor2_class)
    ~ "Exon-Exon",
    (grepl("exon", anchor1_class) & anchor2_class == "cCRE") |
      (anchor1_class == "cCRE" & grepl("exon", anchor2_class))
    ~ "Exon-cCRE",
    (grepl("exon", anchor1_class) & anchor2_class == "intron") |
      (anchor1_class == "intron" & grepl("exon", anchor2_class))
    ~ "Exon-Intron",
    anchor1_class == "cCRE" & anchor2_class == "cCRE"
    ~ "cCRE-cCRE",
    anchor1_class == "intron" & anchor2_class == "intron"
    ~ "Intron-Intron",
    TRUE ~ "Other"
  )
  
  out <- hic_loops %>%
    dplyr::mutate(
      anchor1_class    = anchor1_class,
      anchor2_class    = anchor2_class,
      loop_anchor_type = factor(
        loop_type,
        levels = c("Exon-Exon",
                   "Exon-cCRE",
                   "Exon-Intron",
                   "cCRE-cCRE",
                   "Intron-Intron",
                   "Other")
      )
    )
  
  message("Final loop_anchor_type breakdown:")
  print(table(out$loop_anchor_type))
  
  out
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Plot — WT and Pol2 as two INDEPENDENT panels (not paired/dodged), WT on
#    the left, Pol2 on the right. Each panel is filtered and counted on its
#    own value column only (no joint scale_wt>0 & scale_deg>0 requirement),
#    so the n per anchor type can genuinely differ between panels. Colors
#    are by anchor type — same palette as the original single-condition
#    raincloud — since these are two separate plots, not a WT-vs-Pol2
#    comparison within one plot.
# ══════════════════════════════════════════════════════════════════════════════

ANCHOR_TYPE_COLORS <- c(
  "Exon-Exon"     = "#d6604d",
  "Exon-cCRE"     = "#984ea3",
  "Exon-Intron"   = "#4DAF4A",
  "cCRE-cCRE"     = "#ff7f00",
  "Intron-Intron" = "#2166ac"
)
ANCHOR_TYPE_LEVELS <- names(ANCHOR_TYPE_COLORS)

prepare_hic_anchor_plot_data <- function(classified_loops, value_col, source_filter,
                                         group_levels = ANCHOR_TYPE_LEVELS) {
  classified_loops %>%
    filter(
      loop_anchor_type %in% group_levels,
      source %in% source_filter,          # ← actual loop membership, not background presence
      !is.na(.data[[value_col]]),
      .data[[value_col]] > 0
    ) %>%
    dplyr::mutate(
      log10_observed   = pmin(log10(.data[[value_col]] + 1), log10(1001)),
      loop_anchor_type = factor(loop_anchor_type, levels=group_levels)
    )
}

plot_hic_anchor_raincloud_single <- function(classified_loops, value_col, source_filter, panel_title,
                                             group_levels = ANCHOR_TYPE_LEVELS,
                                             group_colors = ANCHOR_TYPE_COLORS,
                                             y_limits     = NULL) {
  
  plot_df <- prepare_hic_anchor_plot_data(classified_loops, value_col, source_filter, group_levels)
  
  message(panel_title, " — Hi-C loop anchor type breakdown ",
          "(source restricted to: ", paste(source_filter, collapse=", "), "):")
  print(table(plot_df$loop_anchor_type))
  
  counts <- plot_df %>%
    dplyr::count(loop_anchor_type) %>%
    dplyr::mutate(
      loop_anchor_type = factor(loop_anchor_type, levels=group_levels),
      label = paste0(loop_anchor_type, "\n(n=", n, ")")
    ) %>%
    dplyr::arrange(loop_anchor_type)
  
  plot_df <- plot_df %>%
    left_join(counts %>% dplyr::select(loop_anchor_type, label), by="loop_anchor_type") %>%
    dplyr::mutate(label=factor(label, levels=counts$label))
  
  y_scale <- scale_y_continuous(
    name   = paste0("SCALE contacts (", STRAW_RESOLUTION/1000, "kb, log10)"),
    limits = if (!is.null(y_limits)) y_limits else c(NA, log10(1001)*1.05),
    breaks = c(0,1,2,3), labels=c("1","10","100","1000"),
    expand = expansion(mult=c(0.02,0))
  )
  
  ggplot(plot_df, aes(x=label, y=log10_observed,
                      fill=loop_anchor_type, color=loop_anchor_type)) +
    
    ggdist::stat_halfeye(
      aes(thickness = after_stat(pdf * n)),
      side="right", scale=0.4, width=0.6, .width=0, point_colour=NA,
      alpha=0.75, normalize="groups", position=position_nudge(x=0.05)
    ) +
    
    geom_boxplot(
      width=0.12, outlier.shape=NA, color="gray20", alpha=0.9,
      linewidth=0.6, position=position_nudge(x=-0.07)
    ) +
    
    scale_fill_manual(values=group_colors,  guide="none") +
    scale_color_manual(values=group_colors, guide="none") +
    
    y_scale +
    scale_x_discrete(expand=expansion(add=c(0.3, 0.5))) +
    
    labs(title=panel_title, x=NULL) +
    
    theme_classic(base_size=13) +
    theme(
      plot.title   = element_text(face="bold", size=14),
      axis.text.x  = element_text(size=10, angle=30, hjust=1, vjust=1),
      plot.margin  = margin(8, 16, 8, 8),
      aspect.ratio = 1.4
    )
}

plot_hic_anchor_raincloud_wt_pol2 <- function(classified_loops,
                                              group_levels = ANCHOR_TYPE_LEVELS,
                                              group_colors = ANCHOR_TYPE_COLORS,
                                              share_y_axis = TRUE) {
  
  # Each panel is restricted to loops HiCCUPS actually called in that
  # condition (source membership) — NOT "does this position have any
  # nonzero background SCALE signal," which is true almost everywhere in
  # real Hi-C data and was previously making both panels converge on
  # ~the same full union regardless of which condition actually has the loop.
  wt_sources  <- c("shared", "wt_only")
  deg_sources <- c("shared", "pol2_only")
  
  y_limits <- NULL
  if (share_y_axis) {
    wt_vals  <- classified_loops$scale_wt[classified_loops$source %in% wt_sources &
                                            !is.na(classified_loops$scale_wt) & classified_loops$scale_wt > 0]
    deg_vals <- classified_loops$scale_deg[classified_loops$source %in% deg_sources &
                                             !is.na(classified_loops$scale_deg) & classified_loops$scale_deg > 0]
    y_max    <- log10(max(c(wt_vals, deg_vals), na.rm=TRUE) + 1) * 1.05
    y_limits <- c(NA, y_max)
  }
  
  p_wt <- plot_hic_anchor_raincloud_single(
    classified_loops, value_col="scale_wt", source_filter=wt_sources, panel_title="Wildtype",
    group_levels=group_levels, group_colors=group_colors, y_limits=y_limits
  )
  p_deg <- plot_hic_anchor_raincloud_single(
    classified_loops, value_col="scale_deg", source_filter=deg_sources, panel_title="POLR2A deg",
    group_levels=group_levels, group_colors=group_colors, y_limits=y_limits
  )
  
  composite <- cowplot::plot_grid(p_wt, p_deg, ncol=2)
  
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label("Hi-C loop contact strength by anchor type",
                        fontface="bold", size=15, x=0, hjust=0, y=0.7) +
    cowplot::draw_label(paste0(DS_GROUP, "  |  SCALE-normalized contacts via strawr  |  ",
                               "each panel = that condition's own HiCCUPS-called loops  |  ",
                               "WT = left, Pol2 = right"),
                        size=9, color="gray40", x=0, hjust=0, y=0.25)
  
  final <- cowplot::plot_grid(title_grob, composite, ncol=1, rel_heights=c(0.12, 1))
  print(final)
  invisible(list(wt=p_wt, deg=p_deg, composite=final))
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Run
# ══════════════════════════════════════════════════════════════════════════════

hic_classified <- classify_hic_loop_anchors_wt(
  matched, exon_gr, intron_gr, ccre_gr=ccre_gr
)

message("\nFinal loop anchor type breakdown:")
print(table(hic_classified$loop_anchor_type))

plot_hic_anchor_raincloud_wt_pol2(hic_classified)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Export top differential loops for Hi-C locus plotting — fully annotated
#     Pulls the 100 most-lost and 100 most-gained loops by log2fc_loop and
#     writes a CSV annotated with: anchor coordinates + centroids, cCRE/exon/
#     intron anchor classification (from hic_classified, so it uses the same
#     tight centroid-window fix as the raincloud plots above), and the
#     overlapping gene(s) + DE status — enough to hand individual examples to
#     a locus-plotting tool (Juicebox, pyGenomeTracks, cooltools, etc.) and
#     know immediately what's biologically at that position.
#     Placed here (end of script) because it needs hic_classified, which
#     isn't built until exon_gr/intron_gr/ccre_gr exist.
# ══════════════════════════════════════════════════════════════════════════════

GENE_TOL_BP <- 5000L  # padding around the gene body, same convention as the APA/LFC pipelines

# ── Gene-linking variant that includes ALL genes (not just DE), using the
#    same tight centroid-window anchor definition as the anchor classifiers
#    above, so a wide raw bedpe anchor span can't spuriously link to a
#    nearby gene the loop isn't actually near. ────────────────────────────────
link_loops_to_genes_full <- function(loops_df, gene_gr,
                                     anchor_window = STRAW_RESOLUTION,
                                     tol_bp        = GENE_TOL_BP) {
  
  fix_chr <- function(gr, reference_gr) {
    ref_chroms  <- as.character(seqnames(reference_gr))
    has_chr_ref <- any(grepl("^chr", ref_chroms))
    gr_chroms   <- as.character(seqnames(gr))
    has_chr_gr  <- any(grepl("^chr", gr_chroms))
    if (has_chr_ref && !has_chr_gr) {
      seqlevels(gr) <- paste0("chr", seqlevels(gr))
    } else if (!has_chr_ref && has_chr_gr) {
      seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    }
    gr
  }
  
  half_w <- anchor_window / 2
  anchor1_gr <- GRanges(loops_df$chr1, IRanges(loops_df$centroid1 - half_w, loops_df$centroid1 + half_w))
  anchor2_gr <- GRanges(loops_df$chr2, IRanges(loops_df$centroid2 - half_w, loops_df$centroid2 + half_w))
  anchor1_gr <- fix_chr(anchor1_gr, gene_gr)
  anchor2_gr <- fix_chr(anchor2_gr, gene_gr)
  
  gene_gr_exp <- trim(resize(gene_gr, width(gene_gr) + 2*tol_bp, fix="center"))
  
  hits1 <- findOverlaps(anchor1_gr, gene_gr_exp)
  hits2 <- findOverlaps(anchor2_gr, gene_gr_exp)
  
  pairs <- bind_rows(
    data.frame(loop_idx=loops_df$loop_idx[queryHits(hits1)], gene_idx=subjectHits(hits1)),
    data.frame(loop_idx=loops_df$loop_idx[queryHits(hits2)], gene_idx=subjectHits(hits2))
  ) %>%
    dplyr::distinct(loop_idx, gene_idx) %>%
    dplyr::mutate(
      symbol      = gene_gr$symbol[gene_idx],
      log2fc_gene = gene_gr$log2fc_gene[gene_idx],
      padj        = gene_gr$padj[gene_idx],
      de_status   = as.character(gene_gr$de_status[gene_idx])
    )
  
  if (nrow(pairs) == 0) {
    message("No loop-gene overlaps found at all — returning loops_df unannotated.")
    return(loops_df %>% dplyr::mutate(
      gene_symbol=NA_character_, gene_de_status="No gene",
      gene_log2fc=NA_real_, gene_padj=NA_real_,
      all_overlapping_genes="", n_overlapping_genes=0L
    ))
  }
  
  # One representative gene per loop: DE (Up/Down) genes ranked ahead of
  # Not DE / Filtered, tie-broken by smallest padj (no TPM available in
  # this script's gene_gr to rank by expression instead).
  primary_gene <- pairs %>%
    dplyr::mutate(de_priority = dplyr::case_when(de_status %in% c("Up","Down") ~ 1L, TRUE ~ 2L)) %>%
    dplyr::arrange(loop_idx, de_priority, padj) %>%
    dplyr::group_by(loop_idx) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(loop_idx, gene_symbol=symbol, gene_de_status=de_status,
                  gene_log2fc=log2fc_gene, gene_padj=padj)
  
  all_genes <- pairs %>%
    dplyr::group_by(loop_idx) %>%
    dplyr::summarise(
      all_overlapping_genes = paste(unique(symbol), collapse=","),
      n_overlapping_genes   = dplyr::n_distinct(symbol),
      .groups = "drop"
    )
  
  loops_df %>%
    dplyr::left_join(primary_gene, by="loop_idx") %>%
    dplyr::left_join(all_genes,    by="loop_idx") %>%
    dplyr::mutate(
      gene_symbol           = tidyr::replace_na(gene_symbol, "None"),
      gene_de_status        = tidyr::replace_na(gene_de_status, "No gene"),
      all_overlapping_genes = tidyr::replace_na(all_overlapping_genes, ""),
      n_overlapping_genes   = tidyr::replace_na(n_overlapping_genes, 0L)
    )
}

message("\n── Annotating union loops with overlapping genes (all, not just DE) ──")
matched_gene_annotated <- link_loops_to_genes_full(matched, gene_gr)
message("Loops with at least one overlapping gene: ",
        sum(matched_gene_annotated$gene_symbol != "None"), " / ", nrow(matched_gene_annotated))

# ── Combine cCRE/exon/intron classification (hic_classified) with the gene
#    annotation, then pull the top differential loops for export ────────────
N_TOP_DIFFERENTIAL <- 100
PLOT_PADDING_BP     <- 50000L  # padding around the loop span for locus plots

if (any(matched$chr1 != matched$chr2)) {
  message("NOTE: found inter-chromosomal rows in `matched` — plot_region below ",
          "assumes intra-chromosomal loops and uses chr1 only; check these rows.")
}

fully_annotated_loops <- hic_classified %>%
  dplyr::select(loop_idx, anchor1_class, anchor2_class, loop_anchor_type) %>%
  dplyr::left_join(matched_gene_annotated, ., by="loop_idx") %>%
  dplyr::filter(!is.na(log2fc_loop), loop_class %in% c("gained", "lost")) %>%
  dplyr::mutate(
    plot_chr    = chr1,
    plot_start  = pmax(0, pmin(start1, start2) - PLOT_PADDING_BP),
    plot_end    = pmax(end1, end2) + PLOT_PADDING_BP,
    plot_region = paste0(plot_chr, ":", plot_start, "-", plot_end),
    loop_label  = paste0(chr1, ":", start1, "-", end1, " x ", chr2, ":", start2, "-", end2)
  ) %>%
  dplyr::select(
    loop_idx, source, loop_class, log2fc_loop,
    chr1, start1, end1, centroid1,
    chr2, start2, end2, centroid2,
    loop_dist, obs_wt, obs_deg,
    anchor1_class, anchor2_class, loop_anchor_type,
    gene_symbol, gene_de_status, gene_log2fc, gene_padj,
    all_overlapping_genes, n_overlapping_genes,
    loop_label, plot_region
  )

top_lost <- fully_annotated_loops %>%
  filter(loop_class == "lost") %>%
  dplyr::arrange(log2fc_loop) %>%              # most negative first = most lost
  dplyr::slice_head(n=N_TOP_DIFFERENTIAL)

top_gained <- fully_annotated_loops %>%
  filter(loop_class == "gained") %>%
  dplyr::arrange(dplyr::desc(log2fc_loop)) %>% # most positive first = most gained
  dplyr::slice_head(n=N_TOP_DIFFERENTIAL)

top_export <- bind_rows(top_lost, top_gained) %>%
  dplyr::arrange(dplyr::desc(abs(log2fc_loop)))

message("Top differential loops for export: ",
        nrow(top_lost), " lost, ", nrow(top_gained), " gained")
message("Of these, with a DE gene overlap: ",
        sum(top_export$gene_de_status %in% c("Up","Down")), " / ", nrow(top_export))

export_path <- paste0("top_differential_loops_annotated_", DS_GROUP, ".csv")
readr::write_csv(top_export, export_path)
message("Wrote: ", export_path)