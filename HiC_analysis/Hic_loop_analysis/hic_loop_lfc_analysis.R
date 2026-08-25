##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# HiCCUPS loop LFC analysis -- Gained/Stable/Lost classification.
# Contacts extracted via strawr at all loop positions (both
# conditions). Anchor classification via narrowPeak (ChIP-seq +
# ATAC) + TxDb.
#
# Figure output: Fig. 5A, 5C, 5D (probable); SI Fig. 8A, 8F (probable).
##################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(cowplot)
library(GenomicRanges)
library(AnnotationDbi)
library(scales)
library(ggdist)
library(ggridges)
library(strawr)
source("~/Transcription Paper/Hic_loop_analysis/final_scripts/EP_analysis_functions.R")

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════

# ── Downsampling group switch ──────────────────────────────────────────────────
DS_GROUP <- "ds_wt_pol2"   # ← "ds_wt_pol2" | "ds_all" | "unsampled"

BASE     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
LOOP_DIR <- file.path(BASE, "loops")

if (DS_GROUP == "ds_wt_pol2") {
  hic_wt_loops   <- file.path(LOOP_DIR, "ds_wt_pol2/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_pol2_loops <- file.path(LOOP_DIR, "ds_wt_pol2/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file    <- file.path(BASE, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file  <- file.path(BASE, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT  <- 1.0
  DEPTH_SCALE_POL2 <- 1.0
  MIN_OBS_WT       <- 0
  MIN_OBS_POL2     <- 0
  
} else if (DS_GROUP == "ds_all") {
  hic_wt_loops   <- file.path(LOOP_DIR, "ds_all/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_pol2_loops <- file.path(LOOP_DIR, "ds_all/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file    <- file.path(BASE, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file  <- file.path(BASE, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT  <- 1.0
  DEPTH_SCALE_POL2 <- 1.0
  MIN_OBS_WT       <- 0
  MIN_OBS_POL2     <- 0
  
} else if (DS_GROUP == "unsampled") {
  hic_wt_loops   <- file.path(LOOP_DIR, "unsampled/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_pol2_loops <- file.path(LOOP_DIR, "unsampled/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file    <- file.path(BASE, "mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file  <- file.path(BASE, "mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  wt_total         <- 639409109 + 377220354
  pol2_total       <- 374163054 + 249677343
  DEPTH_SCALE_WT  <- 1.0
  DEPTH_SCALE_POL2 <- wt_total / pol2_total
  MIN_OBS_WT       <- 10
  MIN_OBS_POL2     <- 25
  
} else {
  stop("Unknown DS_GROUP: '", DS_GROUP, "'")
}
message("DS_GROUP: ", DS_GROUP)

# ── strawr parameters ──────────────────────────────────────────────────────────
STRAW_NORM       <- "SCALE"
STRAW_RESOLUTION <- 5000L    # bp — must match the resolution used for loop calling

# ── LFC / classification thresholds ───────────────────────────────────────────
LOG2FC_THRESHOLD   <- 0.5
MIN_WT_CONTACTS    <- 5      # minimum depth-scaled WT contacts to call "Lost"
MIN_POL2_CONTACTS  <- 5      # minimum depth-scaled pol2 contacts to call "Gained"
PSEUDOCOUNT        <- 1      # added before log2 to avoid log(0)
ANCHOR_SLOP_BP     <- 5000L  # slop for shared-loop matching

# ── Expression / DE thresholds ────────────────────────────────────────────────
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1
GENE_TOL_BP    <- 5000
BINS           <- c("Low","Medium","High")

# ── File paths ─────────────────────────────────────────────────────────────────
EXPR_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
DESEQ_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
TXDB_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/GRCh38.115.sqlite"
CHIP_DIR   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"
ATAC_DIR   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks"

# ══════════════════════════════════════════════════════════════════════════════
# Constants — colors
# ══════════════════════════════════════════════════════════════════════════════
COND_COLORS  <- c("Wildtype"="#88CCEE", "POLR2A deg"="#CC6677")
CLASS_COLORS <- c("Gained"="#882255", "Stable"="#cccccc", "Lost"="#4393c3")
BIN_COLORS   <- c("Low"="#fdae61", "Medium"="#f46d43", "High"="#d73027",
                  "Unannotated"="#cccccc")
DE_COLORS    <- c("Up"="#D55E00", "Not DE"="#999999",
                  "Down"="#0072B2", "Filtered"="#cccccc")
GROUP_COLORS <- c("P-P"="#882255", "E-P"="#44AA99", "E-E"="#4393c3")
GROUP_MAP <- c(
  "Promoter \u2014 Promoter"                   = "P-P",
  "Promoter \u2014 Distal-Enhancer"            = "E-P",
  "Promoter \u2014 Proximal-Enhancer"          = "E-P",
  "Distal-Enhancer \u2014 Distal-Enhancer"     = "E-E",
  "Proximal-Enhancer \u2014 Proximal-Enhancer" = "E-E",
  "Proximal-Enhancer \u2014 Distal-Enhancer"   = "E-E"
)
CTCF_NP_PAIR_TYPES  <- c("CTCF \u2014 CTCF", "CTCF \u2014 Promoter", "CTCF \u2014 Enhancer")
CTCF_NP_PAIR_COLORS <- c(
  "CTCF \u2014 CTCF"     = "#40004B",
  "CTCF \u2014 Promoter" = "#762A83",
  "CTCF \u2014 Enhancer" = "#C2A5CF"
)
SCORE_LABEL <- paste0("SCALE contacts (log10, depth-scaled to WT)  |  ", STRAW_NORM)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Expression + DESeq2 data
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Loading expression data ─────────────────────────────────────────────")
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
message("DE breakdown:"); print(table(deseq_res$de_status))

gene_coords <- deseq_res %>%
  filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name)) %>%
  distinct(gene_id, .keep_all=TRUE)
gene_gr <- GRanges(
  seqnames  = as.character(gene_coords$seq_name),
  ranges    = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end),
  symbol    = gene_coords$symbol,
  gene_id   = gene_coords$gene_id,
  de_status = as.character(gene_coords$de_status),
  expr_bin  = as.character(gene_coords$expr_bin),
  TPM       = gene_coords$TPM
)
message("gene_gr: ", length(gene_gr), " genes")

# ══════════════════════════════════════════════════════════════════════════════
# 2. HiCCUPS loop loaders + shared-loop finder
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
                fdrH     <= max_fdr, fdrV  <= max_fdr)
}

find_shared_loops <- function(loops_a, loops_b, slop=ANCHOR_SLOP_BP) {
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

# ══════════════════════════════════════════════════════════════════════════════
# 3. Contact extraction via strawr
#    Queries SCALE contacts from a .hic file at each loop centroid pair.
#    Centroids are snapped to the nearest STRAW_RESOLUTION bin.
#    Returns a numeric vector of contact values, length == nrow(loops_df).
# ══════════════════════════════════════════════════════════════════════════════

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
      mat <- strawr::straw(norm, hic_file, loc1, loc2, "BP", resolution,
                           matrix=matrix)
      if (nrow(mat) == 0) { results[idx] <- NA; next }
      lookup <- setNames(mat$counts,
                         paste(pmin(mat$x, mat$y), pmax(mat$x, mat$y), sep="|"))
      keys         <- paste(pmin(c1_pos[idx], c2_pos[idx]),
                            pmax(c1_pos[idx], c2_pos[idx]), sep="|")
      vals         <- lookup[keys]
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


# ══════════════════════════════════════════════════════════════════════════════
# 4. Load HiCCUPS loops and build union set
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Loading HiCCUPS loops (", DS_GROUP, ") ─────────────────────────────")
loops_wt   <- filter_loops(load_hiccups_bedpe(hic_wt_loops),   MIN_OBS_WT) %>%
  dplyr::mutate(loop_idx=row_number(),
                loop_dist=abs(centroid2-centroid1))
loops_pol2 <- filter_loops(load_hiccups_bedpe(hic_pol2_loops), MIN_OBS_POL2) %>%
  dplyr::mutate(loop_idx=row_number(),
                loop_dist=abs(centroid2-centroid1))
message("WT loops:    ", nrow(loops_wt))
message("pol2 loops:  ", nrow(loops_pol2))

# Find shared / exclusive
message("\n── Finding shared / exclusive loops ────────────────────────────────────")
loop_overlap <- find_shared_loops(loops_wt, loops_pol2)
message("Shared:    ", loop_overlap$n_shared,
        " (", round(100*loop_overlap$n_shared/nrow(loops_wt),1), "% of WT)")
message("WT only:   ", loop_overlap$n_a_only)
message("pol2 only: ", loop_overlap$n_b_only)

# Build union set:
#   shared  → use WT coordinates (representative position for both)
#   wt_only → WT coordinates
#   pol2_only → pol2 coordinates
# "source" column records origin so we know how to interpret results
shared_wt <- loops_wt %>%
  dplyr::slice(loop_overlap$shared$a_idx) %>%
  dplyr::mutate(source="shared")

wt_only_loops <- loops_wt %>%
  dplyr::slice(loop_overlap$a_only) %>%
  dplyr::mutate(source="wt_only")

pol2_only_loops <- loops_pol2 %>%
  dplyr::slice(loop_overlap$b_only) %>%
  dplyr::mutate(source="pol2_only")

union_loops <- bind_rows(shared_wt, wt_only_loops, pol2_only_loops) %>%
  dplyr::mutate(loop_idx=row_number())

message("Union loop set: ", nrow(union_loops),
        "  (shared=", sum(union_loops$source=="shared"),
        "  wt_only=", sum(union_loops$source=="wt_only"),
        "  pol2_only=", sum(union_loops$source=="pol2_only"), ")")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Extract SCALE contacts from both .hic files at all union positions
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Extracting contacts via strawr ──────────────────────────────────────")
wt_contacts_raw   <- extract_loop_contacts(union_loops, hic_wt_file)
pol2_contacts_raw <- extract_loop_contacts(union_loops, hic_pol2_file)

# Depth-scale pol2 contacts to WT sequencing depth
wt_contacts_scaled   <- replace_na(wt_contacts_raw,   0) * DEPTH_SCALE_WT
pol2_contacts_scaled <- replace_na(pol2_contacts_raw, 0) * DEPTH_SCALE_POL2

message("WT contacts  — median: ",   round(median(wt_contacts_scaled),1),
        "  zero: ", sum(wt_contacts_scaled == 0))
message("pol2 contacts — median: ",  round(median(pol2_contacts_scaled),1),
        "  zero: ", sum(pol2_contacts_scaled == 0))

# ══════════════════════════════════════════════════════════════════════════════
# 6. Classify Gained / Stable / Lost
#    LFC = log2((pol2 + PSEUDOCOUNT) / (wt + PSEUDOCOUNT))
#    Uses depth-scaled SCALE contacts for both conditions at all positions.
# ══════════════════════════════════════════════════════════════════════════════
union_loops <- union_loops %>%
  dplyr::mutate(
    wt_contacts   = wt_contacts_scaled,
    pol2_contacts = pol2_contacts_scaled,
    log2fc        = log2((pol2_contacts + PSEUDOCOUNT) /
                           (wt_contacts + PSEUDOCOUNT)),
    loop_class    = dplyr::case_when(
      wt_contacts   >= MIN_WT_CONTACTS   & log2fc < -LOG2FC_THRESHOLD ~ "Lost",
      pol2_contacts >= MIN_POL2_CONTACTS & log2fc >  LOG2FC_THRESHOLD ~ "Gained",
      TRUE                                                             ~ "Stable"
    ),
    loop_class = factor(loop_class, levels=c("Gained","Stable","Lost"))
  )

message("\nLoop class breakdown:")
print(table(union_loops$loop_class, union_loops$source))

# log2FC distribution
print(ggplot(union_loops, aes(x=log2fc, fill=source)) +
        geom_histogram(bins=100, alpha=0.8, position="identity") +
        geom_vline(xintercept=c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
                   linetype="dashed", color=c("#4393c3","#882255"), linewidth=0.8) +
        scale_x_continuous(limits=c(-4,4)) +
        scale_fill_manual(values=c("shared"="#2166ac","wt_only"="#4dac26","pol2_only"="#d01c8b"),
                          name="Loop origin") +
        labs(title="log2FC distribution — union HiCCUPS loops",
             subtitle=paste0("log2((pol2 SCALE + ", PSEUDOCOUNT, ") / (WT SCALE + ",
                             PSEUDOCOUNT, "))  |  thresholds: \u00b1", LOG2FC_THRESHOLD,
                             "  |  ", DS_GROUP),
             x="log2FC (POLR2A deg / WT)", y="Count") +
        theme_classic(base_size=12) +
        theme(plot.title=element_text(face="bold",size=13)))

label_ctcf_pairs <- function(df) {
  df %>%
    dplyr::mutate(
      ctcf_pair_type = dplyr::case_when(
        anchor1_ctcf & anchor2_ctcf                                                ~ "CTCF \u2014 CTCF",
        anchor1_ctcf & anchor2_ccre == "Promoter"                                  ~ "CTCF \u2014 Promoter",
        anchor2_ctcf & anchor1_ccre == "Promoter"                                  ~ "CTCF \u2014 Promoter",
        anchor1_ctcf & anchor2_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")  ~ "CTCF \u2014 Enhancer",
        anchor2_ctcf & anchor1_ccre %in% c("Proximal-Enhancer","Distal-Enhancer")  ~ "CTCF \u2014 Enhancer",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(ctcf_pair_type))
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. Load narrowPeak files + TxDb, classify anchors (condition-matched)
#    WT peaks  → shared + wt_only loops (Lost / Stable)
#    KO peaks  → pol2_only loops (Gained)
#    Bind back together so all loops share the same column schema.
# ══════════════════════════════════════════════════════════════════════════════


message("\n── Loading TxDb ────────────────────────────────────────────────────────")
txdb <- loadDb(TXDB_FILE)

message("\n── Loading narrowPeak files ────────────────────────────────────────────")
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

# Split union loops by source before classifying
loops_wt_ref   <- union_loops %>% filter(source %in% c("shared", "wt_only"))
loops_ko_ref   <- union_loops %>% filter(source == "pol2_only")

message("\n── Classifying WT-reference loops (shared + wt_only) against WT peaks ──")
classified_wt_ref <- classify_anchors_narrowpeak(
  loops_df   = loops_wt_ref,
  atac_gr    = atac_wt,
  h3k4me3_gr = h3k4me3_wt,
  h3k27ac_gr = h3k27ac_wt,
  h3k4me1_gr = h3k4me1_wt,
  ctcf_gr    = ctcf_wt,
  txdb       = txdb,
  tss_region = c(-3000, 3000),
  ccre_label = "WT_narrowpeak"
)

message("\n── Classifying Gained loops (pol2_only) against KO peaks ───────────────")
classified_ko_ref <- classify_anchors_narrowpeak(
  loops_df   = loops_ko_ref,
  atac_gr    = atac_ko,
  h3k4me3_gr = h3k4me3_ko,
  h3k27ac_gr = h3k27ac_ko,
  h3k4me1_gr = h3k4me1_ko,
  ctcf_gr    = ctcf_ko,
  txdb       = txdb,
  tss_region = c(-3000, 3000),
  ccre_label = "KO_narrowpeak"
)

all_loops_classified <- bind_rows(classified_wt_ref, classified_ko_ref) %>%
  dplyr::arrange(loop_idx)

message("Combined classified loops: ", nrow(all_loops_classified))
message("Loop type breakdown:")
print(sort(table(all_loops_classified$loop_ccre_type), decreasing=TRUE)[1:10])
# ══════════════════════════════════════════════════════════════════════════════
# 8. Gene annotation
# ══════════════════════════════════════════════════════════════════════════════
link_loops_to_genes <- function(loops_df, gene_gr, tol_bp=GENE_TOL_BP) {
  fix_chr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr",seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr","",seqlevels(gr))
    gr
  }
  anchor1_gr <- fix_chr(GRanges(loops_df$chr1, IRanges(loops_df$start1, loops_df$end1)), gene_gr)
  anchor2_gr <- fix_chr(GRanges(loops_df$chr2, IRanges(loops_df$start2, loops_df$end2)), gene_gr)
  gene_gr_exp <- trim(resize(gene_gr, width(gene_gr)+2*tol_bp, fix="center"))
  hits1 <- findOverlaps(anchor1_gr, gene_gr_exp)
  hits2 <- findOverlaps(anchor2_gr, gene_gr_exp)
  all_pairs <- bind_rows(
    data.frame(loop_idx=queryHits(hits1), gene_idx=subjectHits(hits1)),
    data.frame(loop_idx=queryHits(hits2), gene_idx=subjectHits(hits2))
  ) %>%
    distinct(loop_idx, gene_idx) %>%
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
      expr_bin  = replace_na(as.character(expr_bin),  "Unannotated"),
      de_status = replace_na(as.character(de_status), "Not DE"),
      symbol    = replace_na(symbol, "Unknown")
    )
}

message("\n── Linking loops to genes ──────────────────────────────────────────────")
all_annotated <- link_loops_to_genes(all_loops_classified, gene_gr)

# ══════════════════════════════════════════════════════════════════════════════
# 9. Helper: stacked bar (count + proportion)
# ══════════════════════════════════════════════════════════════════════════════
plot_class_bars <- function(df, fill_var, fill_colors, fill_levels,
                            title, subtitle,
                            min_label_n=10, min_label_prop=0.02) {
  plot_df <- df %>%
    dplyr::mutate(
      loop_class = factor(loop_class, levels=c("Gained","Stable","Lost")),
      fill_val   = factor(.data[[fill_var]], levels=fill_levels)
    ) %>%
    filter(!is.na(fill_val)) %>%
    dplyr::count(loop_class, fill_val) %>%
    dplyr::group_by(loop_class) %>%
    dplyr::mutate(
      prop        = n / sum(n),
      pct_label   = ifelse(prop >= min_label_prop, paste0(round(prop*100,1),"%"), ""),
      count_label = ifelse(n >= min_label_n, as.character(n), "")
    ) %>%
    dplyr::ungroup()
  totals <- plot_df %>%
    dplyr::group_by(loop_class) %>%
    dplyr::summarise(total=sum(n), .groups="drop")
  
  p_count <- ggplot(plot_df, aes(x=loop_class, y=n, fill=fill_val)) +
    geom_bar(stat="identity", width=0.65, color="white", linewidth=0.3) +
    geom_text(aes(label=count_label), position=position_stack(vjust=0.5),
              color="white", size=3, fontface="bold") +
    geom_text(data=totals, aes(x=loop_class,y=total,label=paste0("n=",total)),
              inherit.aes=FALSE, vjust=-0.4, size=3, color="gray30") +
    scale_fill_manual(values=fill_colors, name=NULL, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    scale_y_continuous(name="Loop count", expand=expansion(mult=c(0,0.12))) +
    labs(x=NULL) + theme_classic(base_size=11) +
    theme(axis.text.x=element_text(size=11, face="bold",
                                   color=CLASS_COLORS[c("Gained","Stable","Lost")]),
          legend.position="none", plot.margin=margin(4,4,4,4))
  
  p_prop <- ggplot(plot_df, aes(x=loop_class, y=prop, fill=fill_val)) +
    geom_bar(stat="identity", width=0.65, color="white", linewidth=0.3) +
    geom_text(aes(label=pct_label), position=position_stack(vjust=0.5),
              color="white", size=3, fontface="bold") +
    scale_fill_manual(values=fill_colors, name=NULL, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    scale_y_continuous(name="Proportion", labels=scales::percent,
                       expand=expansion(mult=c(0,0.05))) +
    labs(x=NULL) + theme_classic(base_size=11) +
    theme(axis.text.x=element_text(size=11, face="bold",
                                   color=CLASS_COLORS[c("Gained","Stable","Lost")]),
          legend.position="right", legend.key.size=unit(0.4,"cm"),
          legend.text=element_text(size=9), plot.margin=margin(4,4,4,4))
  
  composite <- (p_count | p_prop) +
    plot_layout(guides="collect") & theme(legend.position="right")
  
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0, y=0.75) +
    cowplot::draw_label(subtitle, size=8, color="gray40", x=0, hjust=0, y=0.25)
  
  final <- cowplot::plot_grid(title_grob, composite, ncol=1, rel_heights=c(0.12,1))
  print(final); invisible(final)
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. Plot 1: Loop class × EP/PP pair type (P-P / E-P / E-E)
# ══════════════════════════════════════════════════════════════════════════════
ep_loops <- all_annotated %>%
  dplyr::mutate(
    group = dplyr::case_when(
      loop_ccre_type %in% names(GROUP_MAP) ~ GROUP_MAP[loop_ccre_type],
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  dplyr::mutate(group=factor(group, levels=names(GROUP_COLORS)))

plot_class_bars(
  df=ep_loops, fill_var="group",
  fill_colors=GROUP_COLORS, fill_levels=names(GROUP_COLORS),
  title="HiCCUPS loop LFC — Promoter & Enhancer loops",
  subtitle=paste0(DS_GROUP, "  |  narrowPeak + ATAC  |  WT peaks reference  |  ",
                  "SCALE contacts  |  log2FC \u00b1", LOG2FC_THRESHOLD,
                  "  |  slop=", ANCHOR_SLOP_BP/1000, "kb")
)

# ══════════════════════════════════════════════════════════════════════════════
# 11. Plot 2: Loop class × CTCF pair type
# ══════════════════════════════════════════════════════════════════════════════
ctcf_loops <- label_ctcf_pairs(all_annotated)
message("CTCF pair type × class:"); print(table(ctcf_loops$loop_class, ctcf_loops$ctcf_pair_type))

plot_class_bars(
  df=ctcf_loops, fill_var="ctcf_pair_type",
  fill_colors=CTCF_NP_PAIR_COLORS, fill_levels=CTCF_NP_PAIR_TYPES,
  title="HiCCUPS loop LFC — CTCF loops",
  subtitle=paste0(DS_GROUP, "  |  narrowPeak + ATAC  |  WT peaks reference  |  ",
                  "SCALE contacts  |  log2FC \u00b1", LOG2FC_THRESHOLD)
)

# ══════════════════════════════════════════════════════════════════════════════
# 12. Plot 3: Loop class × expression bin
# ══════════════════════════════════════════════════════════════════════════════
plot_class_bars(
  df=all_annotated %>%
    dplyr::mutate(expr_bin=factor(expr_bin, levels=c("Low","Medium","High","Unannotated"))),
  fill_var="expr_bin", fill_colors=BIN_COLORS,
  fill_levels=c("Low","Medium","High","Unannotated"),
  title="HiCCUPS loop LFC — by expression level",
  subtitle=paste0(DS_GROUP, "  |  all loops  |  WT TPM tertiles  |  \u00b1",
                  GENE_TOL_BP/1000, "kb gene tolerance")
)

# ══════════════════════════════════════════════════════════════════════════════
# 13. Plot 4: Loop class × DE status
# ══════════════════════════════════════════════════════════════════════════════
plot_class_bars(
  df=all_annotated %>%
    dplyr::mutate(de_status=factor(de_status, levels=c("Up","Not DE","Down","Filtered"))),
  fill_var="de_status", fill_colors=DE_COLORS,
  fill_levels=c("Up","Not DE","Down","Filtered"),
  title="HiCCUPS loop LFC — by DE status",
  subtitle=paste0(DS_GROUP, "  |  all loops  |  padj < ", PADJ_THRESHOLD,
                  "  |  |log2FC| > ", LFC_THRESHOLD)
)

# ══════════════════════════════════════════════════════════════════════════════
# Loop length (Short/Medium/Long) stratified by expression — Low vs High,
# side by side. Count-only bars (mirrors the LEFT panel of plot_class_bars()).
#   Global tertile breaks computed once across ALL loops, so "Short"/"Medium"/
#   "Long" mean the same absolute size range in both panels.
# ══════════════════════════════════════════════════════════════════════════════
# ── Add a y_max parameter to the count-bar helper ────────────────────────────
plot_class_bars_count_only <- function(df, fill_var, fill_colors, fill_levels,
                                       title, subtitle, min_label_n=10,
                                       y_max=NULL) {
  plot_df <- df %>%
    dplyr::mutate(
      loop_class = factor(loop_class, levels=c("Gained","Stable","Lost")),
      fill_val   = factor(.data[[fill_var]], levels=fill_levels)
    ) %>%
    filter(!is.na(fill_val)) %>%
    dplyr::count(loop_class, fill_val) %>%
    dplyr::group_by(loop_class) %>%
    dplyr::mutate(count_label = ifelse(n >= min_label_n, as.character(n), "")) %>%
    dplyr::ungroup()
  totals <- plot_df %>%
    dplyr::group_by(loop_class) %>%
    dplyr::summarise(total=sum(n), .groups="drop")
  
  y_scale <- if (!is.null(y_max)) {
    scale_y_continuous(name="Loop count", limits=c(0, y_max),
                       expand=expansion(mult=c(0, 0.02)))
  } else {
    scale_y_continuous(name="Loop count", expand=expansion(mult=c(0, 0.12)))
  }
  
  p <- ggplot(plot_df, aes(x=loop_class, y=n, fill=fill_val)) +
    geom_bar(stat="identity", width=0.65, color="white", linewidth=0.3) +
    geom_text(aes(label=count_label), position=position_stack(vjust=0.5),
              color="white", size=3, fontface="bold") +
    geom_text(data=totals, aes(x=loop_class,y=total,label=paste0("n=",total)),
              inherit.aes=FALSE, vjust=-0.4, size=3, color="gray30") +
    scale_fill_manual(values=fill_colors, name=NULL, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    y_scale +
    labs(x=NULL) + theme_classic(base_size=11) +
    theme(axis.text.x=element_text(size=11, face="bold",
                                   color=CLASS_COLORS[c("Gained","Stable","Lost")]),
          legend.position="none", plot.margin=margin(4,4,4,4))
  
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface="bold", size=12, x=0, hjust=0, y=0.75) +
    cowplot::draw_label(subtitle, size=7.5, color="gray40", x=0, hjust=0, y=0.25)
  
  cowplot::plot_grid(title_grob, p, ncol=1, rel_heights=c(0.2, 1))
}

# ── Missing data prep this block depends on, reconstructed: ─────────────────
#    LENGTH_TIER_COLORS was never defined anywhere in this script, and
#    LENGTH_TIER_LEVELS was only defined ~500 lines later (too late for a
#    top-to-bottom script) — low_len_df/high_len_df/shared_legend were never
#    built at all. This mirrors the "Low vs High expression by loop length"
#    feature from the sibling MA-loop script, adapted to this script's
#    all_annotated/loop_dist/expr_bin columns.
LENGTH_TIER_LEVELS <- c("Short","Medium","Long")
LENGTH_TIER_COLORS <- c("Short"="#fdd49e", "Medium"="#fc8d59", "Long"="#7f0000")

length_tertiles <- quantile(all_annotated$loop_dist, c(1/3, 2/3), na.rm=TRUE)
message(sprintf("Loop length tertiles (all loops): Short <%.0fkb | Medium %.0f-%.0fkb | Long >%.0fkb",
                length_tertiles[[1]]/1000, length_tertiles[[1]]/1000,
                length_tertiles[[2]]/1000, length_tertiles[[2]]/1000))

all_annotated_len <- all_annotated %>%
  dplyr::filter(!is.na(loop_dist), loop_dist > 0) %>%
  dplyr::mutate(
    length_bin = dplyr::case_when(
      loop_dist <= length_tertiles[[1]] ~ "Short",
      loop_dist <= length_tertiles[[2]] ~ "Medium",
      TRUE                              ~ "Long"
    ),
    length_bin = factor(length_bin, levels=LENGTH_TIER_LEVELS)
  )

low_len_df  <- all_annotated_len %>% dplyr::filter(expr_bin == "Low")
high_len_df <- all_annotated_len %>% dplyr::filter(expr_bin == "High")

legend_source <- ggplot(data.frame(fill_val=factor(LENGTH_TIER_LEVELS, levels=LENGTH_TIER_LEVELS)),
                        aes(x=1, y=1, fill=fill_val)) +
  geom_col() +
  scale_fill_manual(values=LENGTH_TIER_COLORS, name=NULL) +
  theme(legend.position="bottom", legend.text=element_text(size=10))
shared_legend <- cowplot::get_legend(legend_source)

# ── Compute one shared max across BOTH subsets before building either panel ──
shared_y_max <- max(
  low_len_df  %>% dplyr::count(loop_class) %>% dplyr::pull(n),
  high_len_df %>% dplyr::count(loop_class) %>% dplyr::pull(n)
) * 1.08   # small headroom for the "n=" total label above the tallest bar

p_low <- plot_class_bars_count_only(
  df=low_len_df, fill_var="length_bin",
  fill_colors=LENGTH_TIER_COLORS, fill_levels=LENGTH_TIER_LEVELS,
  title="Low expression",
  subtitle=paste0(DS_GROUP, "  |  loop length tertiles, all loops  |  \u00b1",
                  GENE_TOL_BP/1000, "kb gene tolerance"),
  y_max=shared_y_max
)

p_high <- plot_class_bars_count_only(
  df=high_len_df, fill_var="length_bin",
  fill_colors=LENGTH_TIER_COLORS, fill_levels=LENGTH_TIER_LEVELS,
  title="High expression",
  subtitle=paste0(DS_GROUP, "  |  loop length tertiles, all loops  |  \u00b1",
                  GENE_TOL_BP/1000, "kb gene tolerance"),
  y_max=shared_y_max
)

# ── legend + composite assembly unchanged from before ────────────────────────
composite <- cowplot::plot_grid(p_low, p_high, ncol=2)
overall_title <- cowplot::ggdraw() +
  cowplot::draw_label("HiCCUPS loop LFC \u2014 by loop length, Low vs High expression",
                      fontface="bold", size=14, x=0, hjust=0)
final <- cowplot::plot_grid(overall_title, composite, shared_legend,
                            ncol=1, rel_heights=c(0.08, 1, 0.08))
print(final)

# ══════════════════════════════════════════════════════════════════════════════
# 14. Plot 5: Matched ridgeline — WT vs pol2 SCALE contacts by loop class × type
#     All union loops: extracted contacts give genuine WT and pol2 values
#     even for exclusive loops.
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
      data      = wt_df,
      fill      = scales::alpha("#3d3d3d", 0.60),
      color     = scales::alpha("black",   0.8),
      #alpha     = 1,
      linewidth = 0.7,
      linetype  = "dashed",
      scale     = 0.9,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    ggridges::geom_density_ridges(
      data      = pol2_df,
      aes(fill=ridge_label, color=ridge_label),
      alpha     = 0.70,
      linewidth = 0.6,
      scale     = 0.9,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    geom_text(
      data        = pol2_medians,
      aes(x=med_log, y=y_num + 0.15, label=med_label),
      inherit.aes = FALSE,
      color       = "white",
      size        = 2.5,
      fontface    = "bold"
    ) +
    scale_fill_manual(values=ridge_fill, guide="none") +
    scale_color_manual(values=ridge_color, guide="none") +
    scale_y_discrete(labels=function(x) label_lookup[x]) +
    scale_x_continuous(name=SCORE_LABEL) +
    labs(title=title, subtitle=subtitle, y=NULL,
         caption="Charcoal dashed = WT  |  Color = POLR2A deg") +
    theme_ridges(font_size=11) +
    theme(
      plot.title    = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=8, color="gray40"),
      plot.caption  = element_text(size=8, color="gray40"),
      axis.text.y   = element_text(size=9)
    )
}

# EP/PP ridgeline
ep_annotated <- all_annotated %>%
  dplyr::mutate(
    group = dplyr::case_when(
      loop_ccre_type %in% names(GROUP_MAP) ~ GROUP_MAP[loop_ccre_type],
      TRUE ~ NA_character_
    )
  )

print(make_matched_ridge(
  classified_df = ep_annotated,
  group_colors  = GROUP_COLORS,
  group_levels  = names(GROUP_COLORS),
  group_var     = "group",
  title         = paste0("EP/PP matched ridgeline  |  ", DS_GROUP),
  subtitle      = paste0("All union loop positions  |  SCALE contacts extracted via strawr  |  ",
                         "solid=WT  dashed=pol2  |  narrowPeak + ATAC  |  WT peaks ref")
))

print(make_matched_ridge(
  classified_df = all_annotated %>% filter(expr_bin %in% BINS),
  group_colors  = c("Low"="#fdae61", "Medium"="#f46d43", "High"="#d73027"),
  group_levels  = BINS,
  group_var     = "expr_bin",
  title         = paste0("Expression matched ridgeline  |  ", DS_GROUP),
  subtitle      = paste0("All union loop positions  |  SCALE contacts extracted via strawr  |  ",
                         "charcoal dashed=WT  color=POLR2A deg  |  WT TPM tertiles")
))

# CTCF ridgeline
ctcf_annotated <- label_ctcf_pairs(all_annotated)

print(make_matched_ridge(
  classified_df = ctcf_annotated,
  group_colors  = CTCF_NP_PAIR_COLORS,
  group_levels  = CTCF_NP_PAIR_TYPES,
  group_var     = "ctcf_pair_type",
  title         = paste0("CTCF matched ridgeline  |  ", DS_GROUP),
  subtitle      = paste0("All union loop positions  |  SCALE contacts extracted via strawr  |  ",
                         "solid=WT  dashed=pol2  |  narrowPeak + ATAC  |  WT peaks ref")
))

# ══════════════════════════════════════════════════════════════════════════════
# 14b. Plot 2b: Gained vs Lost mirror bar — anchor type composition
#      Back-to-back proportional bars, Stable excluded.
#      Gained classified against KO peaks; Lost against WT peaks.
# ══════════════════════════════════════════════════════════════════════════════
plot_gained_lost_mirror <- function(df, fill_var, fill_colors, fill_levels,
                                    title, subtitle, min_label_prop=0.03) {
  plot_df <- df %>%
    filter(loop_class %in% c("Gained","Lost")) %>%
    dplyr::mutate(
      loop_class = factor(loop_class, levels=c("Gained","Lost")),
      fill_val   = factor(.data[[fill_var]], levels=fill_levels)
    ) %>%
    filter(!is.na(fill_val)) %>%
    dplyr::count(loop_class, fill_val) %>%
    dplyr::group_by(loop_class) %>%
    dplyr::mutate(
      prop      = n / sum(n),
      # Mirror: Gained goes right (+), Lost goes left (-)
      prop_dir  = ifelse(loop_class == "Gained", prop, -prop),
      pct_label = ifelse(prop >= min_label_prop, paste0(round(prop*100,1),"%"), "")
    ) %>%
    dplyr::ungroup()
  
  totals <- plot_df %>%
    dplyr::group_by(loop_class) %>%
    dplyr::summarise(total=sum(n), .groups="drop") %>%
    dplyr::mutate(
      x_pos = ifelse(loop_class == "Gained", 0.52, -0.52),
      label = paste0(loop_class, "\n(n=", total, ")")
    )
  
  x_lim <- max(abs(plot_df$prop_dir), na.rm=TRUE)
  x_lim <- ceiling(x_lim * 10) / 10 + 0.05
  
  ggplot(plot_df, aes(x=prop_dir, y=fill_val, fill=fill_val)) +
    geom_col(data=plot_df %>% filter(loop_class=="Gained"),
             aes(x=prop_dir), width=0.7, color="white", linewidth=0.3) +
    geom_col(data=plot_df %>% filter(loop_class=="Lost"),
             aes(x=prop_dir), width=0.7, color="white", linewidth=0.3) +
    geom_text(data=plot_df %>% filter(loop_class=="Gained"),
              aes(x=prop_dir, label=pct_label),
              hjust=-0.1, size=3, color="gray20") +
    geom_text(data=plot_df %>% filter(loop_class=="Lost"),
              aes(x=prop_dir, label=pct_label),
              hjust=1.1, size=3, color="gray20") +
    geom_vline(xintercept=0, color="gray30", linewidth=0.6) +
    geom_text(data=totals,
              aes(x=x_pos, y=Inf, label=label),
              inherit.aes=FALSE, vjust=1.4, size=3.5, fontface="bold",
              color=CLASS_COLORS[as.character(totals$loop_class)]) +
    scale_fill_manual(values=fill_colors, name=NULL, drop=FALSE) +
    scale_x_continuous(
      limits = c(-x_lim, x_lim),
      breaks = seq(-x_lim, x_lim, by=0.1),
      labels = function(x) scales::percent(abs(x)),
      name   = "\u2190 Lost    Proportion    Gained \u2192"
    ) +
    scale_y_discrete(limits=rev(fill_levels)) +
    labs(title=title, subtitle=subtitle, y=NULL) +
    theme_classic(base_size=11) +
    theme(
      plot.title    = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=8, color="gray40"),
      legend.position = "right",
      legend.key.size = unit(0.4,"cm"),
      legend.text     = element_text(size=9),
      panel.grid.major.x = element_line(color="gray90", linewidth=0.3),
      axis.line.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# EP/PP mirror
ep_mirror <- all_annotated %>%
  dplyr::mutate(
    fill_val = dplyr::case_when(
      loop_ccre_type %in% names(GROUP_MAP) ~ GROUP_MAP[loop_ccre_type],
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(fill_val)) %>%
  dplyr::rename(group=fill_val)

print(plot_gained_lost_mirror(
  df          = ep_mirror,
  fill_var    = "group",
  fill_colors = GROUP_COLORS,
  fill_levels = names(GROUP_COLORS),
  title       = "Gained vs Lost — Promoter & Enhancer anchor composition",
  subtitle    = paste0(DS_GROUP, "  |  Gained=KO peaks  Lost=WT peaks  |  ",
                       "log2FC \u00b1", LOG2FC_THRESHOLD, "  |  Stable excluded")
))

# CTCF mirror
ctcf_mirror <- label_ctcf_pairs(all_annotated)

print(plot_gained_lost_mirror(
  df          = ctcf_mirror,
  fill_var    = "ctcf_pair_type",
  fill_colors = CTCF_NP_PAIR_COLORS,
  fill_levels = CTCF_NP_PAIR_TYPES,
  title       = "Gained vs Lost — CTCF pair type composition",
  subtitle    = paste0(DS_GROUP, "  |  Gained=KO peaks  Lost=WT peaks  |  ",
                       "log2FC \u00b1", LOG2FC_THRESHOLD, "  |  Stable excluded")
))

# ══════════════════════════════════════════════════════════════════════════════
# 15. Loop length by class — ridgeline + stacked bar
# ══════════════════════════════════════════════════════════════════════════════
LENGTH_BIN_COLORS <- c(
  "< 100 kb"     = "#fdd49e",
  "100–500 kb"   = "#fc8d59",
  "500 kb–1 Mb"  = "#d7301f",
  "> 1 Mb"       = "#7f0000"
)
LENGTH_BIN_LEVELS <- c("< 100 kb", "100–500 kb", "500 kb–1 Mb", "> 1 Mb")

loop_length_df <- all_annotated %>%
  filter(!is.na(loop_dist), loop_dist > 0) %>%
  dplyr::mutate(
    length_bin = dplyr::case_when(
      loop_dist <  100000 ~ "< 100 kb",
      loop_dist <  500000 ~ "100\u2013500 kb",
      loop_dist < 1000000 ~ "500 kb\u20131 Mb",
      TRUE                ~ "> 1 Mb"
    ),
    length_bin = factor(length_bin, levels=LENGTH_BIN_LEVELS)
  )

# Stacked bar — proportion of each length bin within each loop class
plot_class_bars(
  df          = loop_length_df,
  fill_var    = "length_bin",
  fill_colors = LENGTH_BIN_COLORS,
  fill_levels = LENGTH_BIN_LEVELS,
  title       = "HiCCUPS loop LFC — by loop length",
  subtitle    = paste0(DS_GROUP, "  |  all union loops  |  log2FC \u00b1", LOG2FC_THRESHOLD)
)

# Ridgeline — one distribution per class, loop size on x-axis
length_medians <- loop_length_df %>%
  dplyr::group_by(loop_class) %>%
  dplyr::summarise(
    med_log = median(log10(loop_dist), na.rm=TRUE),
    med_kb  = median(loop_dist / 1000, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_num     = as.numeric(factor(loop_class, levels=c("Gained","Stable","Lost"))),
    med_label = paste0(round(med_kb), " kb")
  )

n_per_class <- loop_length_df %>%
  dplyr::count(loop_class) %>%
  dplyr::mutate(axis_label=paste0(loop_class, " (n=", n, ")"))
label_lookup_length <- setNames(n_per_class$axis_label, as.character(n_per_class$loop_class))

print(ggplot(loop_length_df,
             aes(x=log10(loop_dist), y=loop_class,
                 fill=loop_class, color=loop_class)) +
        ggridges::geom_density_ridges(scale=0.9, alpha=0.75, linewidth=0.6,
                                      quantile_lines=TRUE, quantiles=2) +
        geom_text(data=length_medians,
                  aes(x=med_log, y=y_num + 0.15, label=med_label),
                  inherit.aes=FALSE, color="white", size=2.8, fontface="bold") +
        scale_fill_manual(values=CLASS_COLORS,  guide="none") +
        scale_color_manual(values=CLASS_COLORS, guide="none") +
        scale_y_discrete(labels=function(x) label_lookup_length[x]) +
        scale_x_continuous(
          name   = "Loop size",
          breaks = log10(c(10000, 50000, 100000, 500000, 1000000, 5000000)),
          labels = c("10 kb","50 kb","100 kb","500 kb","1 Mb","5 Mb")
        ) +
        labs(
          title    = paste0("Loop length distribution by class  |  ", DS_GROUP),
          subtitle = paste0("All union loops  |  log2FC \u00b1", LOG2FC_THRESHOLD,
                            "  |  white label = median"),
          y        = NULL,
          caption  = "Gained = pol2-origin coordinates  |  Stable/Lost = WT-origin coordinates"
        ) +
        theme_ridges(font_size=11) +
        theme(
          plot.title    = element_text(face="bold", size=13),
          plot.subtitle = element_text(size=8,  color="gray40"),
          plot.caption  = element_text(size=8,  color="gray40"),
          axis.text.y   = element_text(size=11, face="bold",
                                       color=CLASS_COLORS[c("Lost","Stable","Gained")])
        ))


# ══════════════════════════════════════════════════════════════════════════════
# Gained/Lost ridgelines by DE status, split by loop source (pol2_only / shared)
#   4 ridge rows: Gained-Up, Lost-Up, Gained-Down, Lost-Down
#   Each row = colored ridge (pol2 contacts, colored by class) overlaid on
#              charcoal dashed ridge (WT contacts, same loop positions)
# ══════════════════════════════════════════════════════════════════════════════
plot_gained_lost_ridge_by_de <- function(df,
                                         source_filter,
                                         title,
                                         subtitle,
                                         wt_color = "#3d3d3d") {
  
  base <- df %>%
    dplyr::filter(source == source_filter,
                  loop_class %in% c("Gained","Lost"),
                  de_status %in% c("Up","Down")) %>%
    dplyr::mutate(
      category = factor(paste0(loop_class, " \u2014 ", de_status),
                        levels = c("Gained \u2014 Up", "Lost \u2014 Up",
                                   "Gained \u2014 Down", "Lost \u2014 Down"))
    ) %>%
    dplyr::filter(!is.na(category))
  
  if (nrow(base) == 0) {
    message("No loops matched source='", source_filter,
            "' with class Gained/Lost x de_status Up/Down.")
    return(invisible(NULL))
  }
  
  message(source_filter, " category counts:")
  print(table(base$category))
  
  ridge_levels <- c("Gained \u2014 Up", "Lost \u2014 Up",
                    "Gained \u2014 Down", "Lost \u2014 Down")
  
  ridge_fill <- c(
    "Gained \u2014 Up"   = scales::alpha(CLASS_COLORS[["Gained"]], 0.75),
    "Lost \u2014 Up"     = scales::alpha(CLASS_COLORS[["Lost"]],   0.75),
    "Gained \u2014 Down" = scales::alpha(CLASS_COLORS[["Gained"]], 0.75),
    "Lost \u2014 Down"   = scales::alpha(CLASS_COLORS[["Lost"]],   0.75)
  )
  ridge_color <- c(
    "Gained \u2014 Up"   = CLASS_COLORS[["Gained"]],
    "Lost \u2014 Up"     = CLASS_COLORS[["Lost"]],
    "Gained \u2014 Down" = CLASS_COLORS[["Gained"]],
    "Lost \u2014 Down"   = CLASS_COLORS[["Lost"]]
  )
  
  long_df <- base %>%
    dplyr::select(category, wt_contacts, pol2_contacts) %>%
    tidyr::pivot_longer(cols=c(wt_contacts, pol2_contacts),
                        names_to="cond", values_to="contacts") %>%
    dplyr::filter(!is.na(contacts), contacts > 0) %>%
    dplyr::mutate(
      cond_label  = ifelse(cond == "wt_contacts", "WT", "POLR2A deg"),
      ridge_label = factor(as.character(category), levels=rev(ridge_levels))
    ) %>%
    dplyr::filter(!is.na(ridge_label))
  
  wt_df   <- long_df %>% dplyr::filter(cond_label == "WT")
  pol2_df <- long_df %>% dplyr::filter(cond_label == "POLR2A deg")
  
  n_labels <- long_df %>%
    dplyr::group_by(ridge_label, cond_label) %>%
    dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
    tidyr::pivot_wider(names_from=cond_label, values_from=n, values_fill=0) %>%
    dplyr::mutate(axis_label=paste0(ridge_label,
                                    " (WT n=", WT, " / pol2 n=", `POLR2A deg`, ")"))
  label_lookup <- setNames(n_labels$axis_label, as.character(n_labels$ridge_label))
  
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
      data      = wt_df,
      fill      = scales::alpha(wt_color, 0.60),
      color     = scales::alpha("black",  0.8),
      linewidth = 0.7,
      linetype  = "dashed",
      scale     = 0.9,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    ggridges::geom_density_ridges(
      data      = pol2_df,
      aes(fill=ridge_label, color=ridge_label),
      alpha     = 0.70,
      linewidth = 0.6,
      scale     = 0.9,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    geom_text(
      data        = pol2_medians,
      aes(x=med_log, y=y_num + 0.15, label=med_label),
      inherit.aes = FALSE,
      color       = "white",
      size        = 2.5,
      fontface    = "bold"
    ) +
    scale_fill_manual(values=ridge_fill, guide="none") +
    scale_color_manual(values=ridge_color, guide="none") +
    scale_y_discrete(labels=function(x) label_lookup[x]) +
    scale_x_continuous(name=SCORE_LABEL) +
    labs(title=title, subtitle=subtitle, y=NULL,
         caption="Charcoal dashed = WT  |  Color = POLR2A deg") +
    theme_ridges(font_size=11) +
    theme(
      plot.title    = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=8, color="gray40"),
      plot.caption  = element_text(size=8, color="gray40"),
      axis.text.y   = element_text(size=9)
    )
}

# ── Plot 1: pol2-only loops ───────────────────────────────────────────────────
print(plot_gained_lost_ridge_by_de(
  df            = all_annotated,
  source_filter = "pol2_only",
  title         = "Pol2-only loops \u2014 Gained/Lost by DE status",
  subtitle      = paste0(DS_GROUP, "  |  source=pol2_only  |  SCALE contacts extracted via strawr  |  ",
                         "solid=pol2, dashed=WT  |  log2FC \u00b1", LOG2FC_THRESHOLD)
))

# ── Plot 2: shared loops ──────────────────────────────────────────────────────
print(plot_gained_lost_ridge_by_de(
  df            = all_annotated,
  source_filter = "shared",
  title         = "Shared loops \u2014 Gained/Lost by DE status",
  subtitle      = paste0(DS_GROUP, "  |  source=shared  |  SCALE contacts extracted via strawr  |  ",
                         "solid=pol2, dashed=WT  |  log2FC \u00b1", LOG2FC_THRESHOLD)
))

# ── Plot 2: shared loops ──────────────────────────────────────────────────────
print(plot_gained_lost_ridge_by_de(
  df            = all_annotated,
  source_filter = "wt_only",
  title         = "WT only loops \u2014 Gained/Lost by DE status",
  subtitle      = paste0(DS_GROUP, "  |  source=shared  |  SCALE contacts extracted via strawr  |  ",
                         "solid=pol2, dashed=WT  |  log2FC \u00b1", LOG2FC_THRESHOLD)
))

# ══════════════════════════════════════════════════════════════════════════════
# Length-tier setup for P-P / E-P ridgelines
# ══════════════════════════════════════════════════════════════════════════════
LENGTH_TIER_LEVELS <- c("Short","Medium","Long")

add_length_tertiles <- function(df, dist_col="loop_dist", labels=LENGTH_TIER_LEVELS) {
  q <- quantile(df[[dist_col]], c(1/3, 2/3), na.rm=TRUE)
  message(sprintf("Loop length tertiles: Short <%.0fkb  |  Medium %.0f\u2013%.0fkb  |  Long >%.0fkb",
                  q[[1]]/1000, q[[1]]/1000, q[[2]]/1000, q[[2]]/1000))
  df %>%
    dplyr::mutate(
      length_tier = dplyr::case_when(
        .data[[dist_col]] <= q[[1]] ~ labels[1],
        .data[[dist_col]] <= q[[2]] ~ labels[2],
        TRUE                       ~ labels[3]
      ),
      length_tier = factor(length_tier, levels=labels)
    )
}

pp_ep_base <- all_annotated %>%
  dplyr::mutate(
    group = dplyr::case_when(
      loop_ccre_type %in% names(GROUP_MAP) ~ GROUP_MAP[loop_ccre_type],
      TRUE ~ NA_character_
    )
  ) %>%
  filter(group %in% c("P-P","E-P"))

pp_ep_annotated <- pp_ep_base %>%
  filter(!is.na(loop_dist), loop_dist > 0) %>%
  add_length_tertiles(dist_col="loop_dist")

message("P-P/E-P loops by length tier:")
print(table(pp_ep_annotated$group, pp_ep_annotated$length_tier))

for (len in LENGTH_TIER_LEVELS) {
  sub_df <- pp_ep_annotated %>% filter(length_tier == len)
  print(make_matched_ridge(
    classified_df = sub_df,
    group_colors  = GROUP_COLORS,
    group_levels  = c("P-P","E-P"),
    group_var     = "group",
    title         = paste0("EP/PP matched ridgeline \u2014 ", len, " loops  |  ", DS_GROUP),
    subtitle      = paste0("Loop length tier: ", len, "  |  All union loop positions  |  ",
                           "SCALE contacts extracted via strawr  |  solid=WT dashed=pol2")
  ))
}

# ══════════════════════════════════════════════════════════════════════════════
# Combined ridgeline: group + length_tier + loop_class, all stacked in one plot
#   Color = group (P-P/E-P), alpha = loop_class (same convention as make_matched_ridge)
#   Row order nests group > length_tier > loop_class
# ══════════════════════════════════════════════════════════════════════════════
make_matched_ridge_length_combined <- function(classified_df,
                                               group_colors,
                                               group_levels  = c("P-P","E-P"),
                                               length_levels = LENGTH_TIER_LEVELS,
                                               title, subtitle) {
  
  class_levels <- c("Gained","Stable","Lost")
  class_alphas <- c(Gained=0.85, Stable=0.55, Lost=0.25)
  
  combo_levels <- c()
  ridge_fill   <- c()
  ridge_color  <- c()
  for (g in group_levels) {
    col <- group_colors[[g]]
    for (l in length_levels) {
      for (cl in class_levels) {
        lbl <- paste0(g, " ", l, " \u2014 ", cl)
        combo_levels <- c(combo_levels, lbl)
        ridge_fill[lbl]  <- scales::alpha(col, class_alphas[[cl]])
        ridge_color[lbl] <- col
      }
    }
  }
  
  base_df <- classified_df %>%
    filter(!is.na(group), group %in% group_levels,
           !is.na(length_tier), length_tier %in% length_levels) %>%
    dplyr::select(group, length_tier, loop_class, wt_contacts, pol2_contacts) %>%
    tidyr::pivot_longer(cols=c(wt_contacts, pol2_contacts),
                        names_to="cond", values_to="contacts") %>%
    filter(!is.na(contacts), contacts > 0) %>%
    dplyr::mutate(
      cond_label  = ifelse(cond == "wt_contacts", "WT", "POLR2A deg"),
      ridge_label = factor(paste0(group, " ", length_tier, " \u2014 ", loop_class),
                           levels=rev(combo_levels))
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
  label_lookup <- setNames(n_labels$axis_label, as.character(n_labels$ridge_label))
  
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
      data=wt_df, fill=scales::alpha("#3d3d3d", 0.60), color=scales::alpha("black", 0.8),
      linewidth=0.7, linetype="dashed", scale=0.9, quantile_lines=TRUE, quantiles=2
    ) +
    ggridges::geom_density_ridges(
      data=pol2_df, aes(fill=ridge_label, color=ridge_label),
      alpha=0.70, linewidth=0.6, scale=0.9, quantile_lines=TRUE, quantiles=2
    ) +
    geom_text(
      data=pol2_medians, aes(x=med_log, y=y_num+0.15, label=med_label),
      inherit.aes=FALSE, color="white", size=2.3, fontface="bold"
    ) +
    scale_fill_manual(values=ridge_fill, guide="none") +
    scale_color_manual(values=ridge_color, guide="none") +
    scale_y_discrete(labels=function(x) label_lookup[x]) +
    scale_x_continuous(name=SCORE_LABEL) +
    labs(title=title, subtitle=subtitle, y=NULL,
         caption="Charcoal dashed = WT  |  Color = POLR2A deg") +
    theme_ridges(font_size=10) +
    theme(
      plot.title=element_text(face="bold", size=13),
      plot.subtitle=element_text(size=8, color="gray40"),
      plot.caption=element_text(size=8, color="gray40"),
      axis.text.y=element_text(size=8)
    )
}

combined_ridge <- make_matched_ridge_length_combined(
  classified_df = pp_ep_annotated,
  group_colors  = GROUP_COLORS,
  title         = paste0("EP/PP matched ridgeline by loop length  |  ", DS_GROUP),
  subtitle      = "All union loop positions  |  SCALE contacts extracted via strawr  |  solid=WT dashed=pol2"
)
print(combined_ridge)

# ── Missing data prep this needs, reconstructed: ─────────────────────────────
#    combo_levels only existed as a LOCAL variable inside
#    make_matched_ridge_length_combined() above — out of scope here — and
#    pp_ep_annotated never had a row_combo column added. Both are needed by
#    the mirror_matched_len call below (row_var="row_combo").
combo_levels <- as.vector(outer(c("P-P","E-P"), LENGTH_TIER_LEVELS,
                                function(g,l) paste0(g, " ", l)))
pp_ep_annotated <- pp_ep_annotated %>%
  dplyr::mutate(row_combo = factor(paste0(group, " ", length_tier), levels=combo_levels))

# ══════════════════════════════════════════════════════════════════════════════
# Mirrored ridge, matched-style: Gained (up) / Lost (down), each side showing
# BOTH its own pol2 (colored) and WT (charcoal dashed) distribution — same
# fill/alpha/linetype convention as make_matched_ridge(), median line drawn
# INSIDE each density (terminating at the curve's own height), not as an
# external tick.
# ══════════════════════════════════════════════════════════════════════════════
plot_gained_lost_mirror_matched <- function(classified_df,
                                            row_var,
                                            row_levels,
                                            title, subtitle,
                                            gained_color = CLASS_COLORS[["Gained"]],
                                            lost_color   = CLASS_COLORS[["Lost"]],
                                            wt_color     = "#3d3d3d",
                                            row_height   = 0.42,
                                            min_n        = 10,
                                            n_points     = 512) {
  
  df <- classified_df %>%
    dplyr::filter(loop_class %in% c("Gained","Lost")) %>%
    dplyr::mutate(row_label = factor(.data[[row_var]], levels=row_levels)) %>%
    dplyr::filter(!is.na(row_label))
  
  row_y <- setNames(seq_along(row_levels), row_levels)
  
  build_density <- function(vals) {
    vals <- vals[!is.na(vals) & vals > 0]
    if (length(vals) < min_n) return(NULL)
    stats::density(log10(vals), n=n_points)
  }
  
  # ── Pass 1: build every raw (un-normalized) density, tagged by row/class ────
  records <- list()
  for (rl in row_levels) {
    for (cl in c("Gained","Lost")) {
      sub    <- df %>% dplyr::filter(row_label == rl, loop_class == cl)
      d_pol2 <- build_density(sub$pol2_contacts)
      d_wt   <- build_density(sub$wt_contacts)
      if (is.null(d_pol2) || is.null(d_wt)) {
        message("Skipping ", rl, " \u2014 ", cl,
                " (pol2 n=", sum(!is.na(sub$pol2_contacts) & sub$pol2_contacts>0),
                ", wt n=",   sum(!is.na(sub$wt_contacts)   & sub$wt_contacts>0),
                ") \u2014 below min_n=", min_n)
        next
      }
      records[[paste(rl, cl)]] <- list(
        row=rl, class=cl,
        pol2_med = median(log10(sub$pol2_contacts[sub$pol2_contacts>0]), na.rm=TRUE),
        wt_med   = median(log10(sub$wt_contacts[sub$wt_contacts>0]),     na.rm=TRUE),
        d_pol2=d_pol2, d_wt=d_wt
      )
    }
  }
  if (length(records) == 0) stop("No categories passed min_n=", min_n)
  
  # ── Two GLOBAL scale factors (mirrors the two separate geom_density_ridges
  #    layers in make_matched_ridge — one shared scale for ALL pol2 curves,
  #    one shared scale for ALL WT curves) ─────────────────────────────────────
  pol2_scale <- row_height / max(sapply(records, function(r) max(r$d_pol2$y)))
  wt_scale   <- row_height / max(sapply(records, function(r) max(r$d_wt$y)))
  interp_h   <- function(d, x0, scale) approx(d$x, d$y, xout=x0)$y * scale
  
  ribbon_rows <- list(); seg_rows <- list(); label_rows <- list()
  for (key in names(records)) {
    r        <- records[[key]]
    baseline <- row_y[[r$row]]
    dir_sign <- if (r$class == "Gained") 1 else -1
    
    ribbon_rows[[paste0(key,"_pol2")]] <- data.frame(
      row=r$row, class=r$class, layer="pol2", x=r$d_pol2$x,
      ymin = baseline + pmin(0, dir_sign * r$d_pol2$y * pol2_scale),
      ymax = baseline + pmax(0, dir_sign * r$d_pol2$y * pol2_scale))
    ribbon_rows[[paste0(key,"_wt")]] <- data.frame(
      row=r$row, class=r$class, layer="wt", x=r$d_wt$x,
      ymin = baseline + pmin(0, dir_sign * r$d_wt$y * wt_scale),
      ymax = baseline + pmax(0, dir_sign * r$d_wt$y * wt_scale))
    
    pol2_h <- interp_h(r$d_pol2, r$pol2_med, pol2_scale)
    wt_h   <- interp_h(r$d_wt,   r$wt_med,   wt_scale)
    
    seg_rows[[paste0(key,"_pol2")]] <- data.frame(
      row=r$row, class=r$class, layer="pol2",
      x=r$pol2_med, y=baseline, yend=baseline + dir_sign * pol2_h)
    seg_rows[[paste0(key,"_wt")]] <- data.frame(
      row=r$row, class=r$class, layer="wt",
      x=r$wt_med, y=baseline, yend=baseline + dir_sign * wt_h)
    
    label_rows[[key]] <- data.frame(
      row=r$row, class=r$class,
      x=r$pol2_med, y=baseline + dir_sign * (pol2_h + 0.04),
      label=sprintf("%.1f", 10^r$pol2_med))
  }
  
  ribbons <- dplyr::bind_rows(ribbon_rows)
  segs    <- dplyr::bind_rows(seg_rows)
  labels  <- dplyr::bind_rows(label_rows)
  
  n_labels <- df %>%
    dplyr::count(row_label, loop_class, .drop=FALSE) %>%
    tidyr::pivot_wider(names_from=loop_class, values_from=n, values_fill=0) %>%
    dplyr::mutate(axis_label = paste0(row_label,
                                      "  (Gained n=", Gained, " / Lost n=", Lost, ")"))
  label_lookup <- setNames(n_labels$axis_label, as.character(n_labels$row_label))
  class_colors <- c(Gained=gained_color, Lost=lost_color)
  
  ggplot() +
    geom_hline(yintercept=row_y, color="gray85", linewidth=0.4) +
    # WT — same alpha/linetype convention as make_matched_ridge's WT layer
    geom_ribbon(data=ribbons %>% dplyr::filter(layer=="wt"),
                aes(x=x, ymin=ymin, ymax=ymax, group=interaction(row,class)),
                fill=scales::alpha(wt_color, 0.60),
                color=scales::alpha("black", 0.8),
                linetype="dashed", linewidth=0.7) +
    # pol2 — same alpha/linewidth convention as make_matched_ridge's pol2 layer
    geom_ribbon(data=ribbons %>% dplyr::filter(layer=="pol2"),
                aes(x=x, ymin=ymin, ymax=ymax, fill=class, color=class,
                    group=interaction(row,class)),
                alpha=0.70, linewidth=0.6) +
    # median lines — drawn INSIDE the density, terminating at its own height
    geom_segment(data=segs %>% dplyr::filter(layer=="wt"),
                 aes(x=x, xend=x, y=y, yend=yend),
                 color=scales::alpha("black",0.8), linetype="dashed", linewidth=0.7) +
    geom_segment(data=segs %>% dplyr::filter(layer=="pol2"),
                 aes(x=x, xend=x, y=y, yend=yend, color=class), linewidth=0.6) +
    geom_label(data=labels, aes(x=x, y=y, label=label),
               fill="white", alpha=0.85, color="black", fontface="bold",
               size=2.5, label.size=0, label.padding=unit(0.12,"lines")) +
    scale_fill_manual(values=class_colors, guide="none") +
    scale_color_manual(values=class_colors, guide="none") +
    scale_y_continuous(breaks=row_y,
                       labels=function(x) label_lookup[names(row_y)[match(x, row_y)]]) +
    scale_x_continuous(name=SCORE_LABEL) +
    labs(title=title, subtitle=subtitle, y=NULL,
         caption="\u2191 Gained / \u2193 Lost = pol2 contacts (colored)  |  gray dashed = matching WT distribution") +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold", size=13),
          plot.subtitle=element_text(size=8, color="gray40"),
          plot.caption=element_text(size=8, color="gray40"),
          axis.text.y=element_text(size=9),
          axis.line.y=element_blank(), axis.ticks.y=element_blank())
}

print(plot_gained_lost_mirror_matched(
  classified_df = pp_ep_base,
  row_var       = "group",
  row_levels    = c("P-P","E-P"),
  title         = paste0("Gained vs Lost \u2014 P-P / E-P  |  ", DS_GROUP),
  subtitle      = "pol2 SCALE contacts + matching WT, log10  |  all union loop positions"
))

mirror_matched_len <- plot_gained_lost_mirror_matched(
  classified_df = pp_ep_annotated,
  row_var       = "row_combo",
  row_levels    = combo_levels,
  title         = paste0("Gained vs Lost by loop length \u2014 P-P / E-P  |  ", DS_GROUP),
  subtitle      = "pol2 SCALE contacts + matching WT, log10  |  all union loop positions"
)
print(mirror_matched_len)
