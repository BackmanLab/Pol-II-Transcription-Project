##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Hi-C loop analyses stratified by EU-seq (nascent) transcription change instead of total RNA-seq
# DESeq2 status: tertiles of the EU-seq ABSOLUTE LOSS in scaled RPKM, WT_level * (1 - 2^LFC) with
# WT_level = 2^mean_wt (log2(RPKM + 1) scale of eu_results), i.e. 2^mean_wt - 2^mean_polr2a. It is
# tertiled as eu_abs = 2^mean_polr2a - 2^mean_wt (negative = transcription lost), so T1 = lowest =
# largest absolute loss. The Ref stratum (|EU-seq LFC| <= REF_LFC) is still defined on the LFC.
# Loop calling, contact extraction (SCALE, 5 kb), Gained / Stable / Lost classification
# (|log2FC| > 0.5 and >= 5 contacts) and APA follow hic_loop_lfc_analysis.R and
# final_APA_analysis.R; only the per-gene stratifier differs (eu_results of
# EU_seq_zscore_analysis.R). A loop takes the stratum of the gene overlapping its anchors
# (highest WT TPM among EU-seq genes); loops with no EU-seq gene at an anchor are left out.
#
# Absolute reference: loops on genes with |EU-seq LFC| <= REF_LFC form an extra column to the left
# of T1 in every stratified plot (these loops also stay in their own tertile).
#
# Scoring (APA panels and scoring figures), each as log2(POLR2A deg / WT); more negative = more
# contact lost; bars are 95% bootstrap CIs over loops:
#   - APA-style 3x3 patch  (3x3 center patch over the four 3x3 corner blocks)
#   - APA-style 5x5 patch  (5x5 center patch over the four 3x3 corner blocks)
#   - whole 15x15 grid     (summed contacts, no background correction)
#   - center pixel         (no background correction)
#
# Outputs (SVGs in PLOT_DIR; every plot is also printed):
#   1  log2FC loop distribution (histogram and kernel density, each with boxplots) and Lost | Gained
#      ridges per stratum, for WT-only / shared / Pol2-only loops
#   2  APA of loops with an anchor on a gene body per stratum, for WT-only / shared / Pol2-only loops
#   3  APA of WT loops by stratum (Ref, T1-T3): WT, POLR2A deg and difference, with the APA
#      3x3 / 5x5 / whole-grid / center-pixel numbers on the panels; scoring figure with CIs
#   4  WT loops split into short and long loops (median loop distance): log2FC density + boxplots,
#      Lost | Gained ridges, APA and scoring figure per distance class
#   5  WT loops split into two halves by steady-state expression (median WT TPM of the EU-seq
#      genes): APA and scoring figure per expression half
#   6  non-gene anchor type of loops with one anchor in a gene body, per stratum
#   7  Gained vs Lost loop counts per stratum: loops with a promoter or enhancer anchor, and all
#      loops at genes
#
# Needs: eu_results in the session (or EU_RESULTS_FILE), EP_analysis_functions.R, ChIPseeker,
# GenomicFeatures, strawr, mariner, InteractionSet, ggridges, cowplot, rtracklayer.
# Anchors of both Wildtype and POLR2A-degraded loops are typed against the untreated (DMSO)
# ATAC, H3K4me3, H3K27ac and H3K4me1 narrowPeaks (loop_chromatin_regression/data).
##################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(cowplot)
library(GenomicRanges)
library(scales)
library(ggridges)
library(strawr)
library(mariner)
library(InteractionSet)
library(RColorBrewer)
library(viridis)

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════

# ── Data source switch ─────────────────────────────────────────────────────────
DATA_SOURCE <- "hic"       # ← "hic" (our Hi-C and its HiCCUPS loops) | "microc" (ENCODE Micro-C and its HiCCUPS loops)
DS_GROUP    <- "ds_wt_pol2"   # Hi-C only: "ds_wt_pol2" | "ds_all" | "unsampled"

BASE     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
LOOP_DIR <- file.path(BASE, "loops")

# Micro-C (ENCODE POLR2A-degron HCT116): WT = ENCSR533VUS, POLR2A deg = ENCSR218URU.
# Each loop bedpe is derived_from the .hic it is paired with (ENCODE file records):
#   ENCFF247KKX <- ENCFF534TYD (WT),  ENCFF452UUK <- ENCFF302CUI (POLR2A deg)
MICROC_DIR            <- "/projects/b1042/BackmanLab/Lucas/microc_encode"
MICROC_WT_LOOPS       <- file.path(MICROC_DIR, "ENCFF247KKX.bedpe.gz")   # HiCCUPS loops, WT
MICROC_POL2_LOOPS     <- file.path(MICROC_DIR, "ENCFF452UUK.bedpe.gz")   # HiCCUPS loops, POLR2A deg
MICROC_DEPTH_SCALE_WT   <- 1.0   # multiplies WT contacts
MICROC_DEPTH_SCALE_POL2 <- 1.0   # multiplies POLR2A deg contacts (set to WT total / POLR2A deg total contacts to depth-match)
MICROC_MIN_OBS_WT     <- 0       # minimum HiCCUPS "observed" to keep a loop (the FDR filters still apply)
MICROC_MIN_OBS_POL2   <- 0

if (DATA_SOURCE == "microc") {
  hic_wt_loops   <- MICROC_WT_LOOPS
  hic_pol2_loops <- MICROC_POL2_LOOPS
  hic_wt_file    <- file.path(MICROC_DIR, "ENCFF534TYD.hic")
  hic_pol2_file  <- file.path(MICROC_DIR, "ENCFF302CUI.hic")
  DEPTH_SCALE_WT   <- MICROC_DEPTH_SCALE_WT
  DEPTH_SCALE_POL2 <- MICROC_DEPTH_SCALE_POL2
  MIN_OBS_WT       <- MICROC_MIN_OBS_WT
  MIN_OBS_POL2     <- MICROC_MIN_OBS_POL2
} else if (DATA_SOURCE != "hic") {
  stop("Unknown DATA_SOURCE: '", DATA_SOURCE, "' (use \"hic\" or \"microc\")")
} else if (DS_GROUP == "ds_wt_pol2") {
  hic_wt_loops   <- file.path(LOOP_DIR, "ds_wt_pol2/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_pol2_loops <- file.path(LOOP_DIR, "ds_wt_pol2/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file    <- file.path(BASE, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file  <- file.path(BASE, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT   <- 1.0
  DEPTH_SCALE_POL2 <- 1.0
  MIN_OBS_WT       <- 0
  MIN_OBS_POL2     <- 0
} else if (DS_GROUP == "ds_all") {
  hic_wt_loops   <- file.path(LOOP_DIR, "ds_all/WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_pol2_loops <- file.path(LOOP_DIR, "ds_all/Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file    <- file.path(BASE, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file  <- file.path(BASE, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  DEPTH_SCALE_WT   <- 1.0
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
  DEPTH_SCALE_WT   <- 1.0
  DEPTH_SCALE_POL2 <- wt_total / pol2_total
  MIN_OBS_WT       <- 10
  MIN_OBS_POL2     <- 25
} else {
  stop("Unknown DS_GROUP: '", DS_GROUP, "'")
}
# label (subtitles / messages) and file-name tag (SVGs) of the chosen data, so Hi-C and Micro-C outputs never overwrite each other
DATA_LABEL <- if (DATA_SOURCE == "microc") "Micro-C" else DS_GROUP
DATA_TAG   <- paste0(if (DATA_SOURCE == "microc") "microc" else paste0("hic-", gsub("_", "", sub("^ds_", "", DS_GROUP))), "-absloss")
message("Data source: ", DATA_SOURCE, " (", DATA_LABEL, ")")
message("  loops WT:   ", hic_wt_loops, "\n  loops deg:  ", hic_pol2_loops,
        "\n  .hic WT:    ", hic_wt_file, "\n  .hic deg:   ", hic_pol2_file)

CONDITIONS   <- list(wt   = list(label="Wildtype",   hic=hic_wt_file),
                     pol2 = list(label="POLR2A deg", hic=hic_pol2_file))
DEPTH_SCALE  <- list(wt=DEPTH_SCALE_WT, pol2=DEPTH_SCALE_POL2)

# ── strawr / LFC / classification (hic_loop_lfc_analysis.R) ────────────────────
STRAW_NORM         <- "SCALE"
STRAW_RESOLUTION   <- 5000L
LOG2FC_THRESHOLD   <- 0.5
MIN_WT_CONTACTS    <- 5
MIN_POL2_CONTACTS  <- 5
PSEUDOCOUNT        <- 1
ANCHOR_SLOP_BP     <- 5000L
GENE_TOL_BP        <- 10000     # loop-to-gene overlap tolerance: gene body +/- 10 kb
GENEBODY_TOL_BP    <- 5000      # same, for the gene-body APA section: gene body +/- 5 kb
# ── APA (final_APA_analysis.R) ─────────────────────────────────────────────────
RESOLUTION <- 5000
BUFFER     <- 7
NORM       <- "NONE"
ACTIVE_PALETTE      <- "viridis"
ACTIVE_DIFF_PALETTE <- "PRGn"

# ── EU-seq stratifier ──────────────────────────────────────────────────────────
LFC_THRESHOLD <- 1              # |LFC| of eu_results$de_status (EU_seq_zscore_analysis.R)
STATUS_LEVELS <- c("Up", "Not DE", "Down")
TERT_LEVELS   <- c("T1 (lowest)", "T2 (middle)", "T3 (highest)")
EU_RESULTS_FILE <- NULL         # optional csv of eu_results, used only if eu_results is not in the session


# ── File paths ─────────────────────────────────────────────────────────────────
HG38_REF      <- "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference"
GTF_FILE      <- file.path(HG38_REF, "Homo_sapiens.GRCh38.115.chr.gtf")   # Ensembl 115
EXPR_FILE     <- file.path(HG38_REF, "ENCFF435PHM.tsv")                    # ENCODE HCT116 RSEM (gene_id, TPM, ...)
PEAK_DIR      <- "/projects/b1042/BackmanLab/Lucas/loop_chromatin_regression/data"
EP_FUNCS_FILE <- "~/Transcription Paper/Hic_loop_analysis/final_scripts/EP_analysis_functions.R"
PLOT_DIR      <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/euseq_stratified_loops"   # SVGs are saved here (created if missing)

NP_FILES <- c(atac    = file.path(PEAK_DIR, "ATAC_DMSO.narrowPeak"),
              h3k4me3 = file.path(PEAK_DIR, "H3K4me3_DMSO.narrowPeak"),
              h3k27ac = file.path(PEAK_DIR, "H3K27ac_DMSO.narrowPeak"),
              h3k4me1 = file.path(PEAK_DIR, "H3K4me1_DMSO.narrowPeak"))
NP_TSS_REGION <- c(-3000, 3000)

check_files <- function(paths) {
  miss <- paths[!file.exists(path.expand(paths))]
  if (length(miss) > 0) stop("Input not found:\n  ", paste(miss, collapse="\n  "), call. = FALSE)
  invisible(TRUE)
}
check_files(c(hic_wt_loops, hic_pol2_loops, hic_wt_file, hic_pol2_file, EXPR_FILE, GTF_FILE, EP_FUNCS_FILE, NP_FILES))
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(PLOT_DIR)) stop("PLOT_DIR could not be created: ", PLOT_DIR)

# ── Constants — colors ─────────────────────────────────────────────────────────
CLASS_COLORS  <- c("Gained"="#882255", "Stable"="#cccccc", "Lost"="#4393c3")
SOURCE_COLORS <- c("shared"="#2166ac", "wt_only"="#4dac26", "pol2_only"="#d01c8b")
DIST_COLORS   <- c(short="#E69F00", long="#5b5b5b")
single_anchor_colors <- c(
  "Promoter"          = "#b2182b",
  "Proximal-Enhancer" = "#f4a582",
  "Distal-Enhancer"   = "#d6604d"
)
ep_anchor_types <- c("Promoter", "Proximal-Enhancer", "Distal-Enhancer")
GENEBODY_GROUP  <- "Gene body"
SCORE_LABEL   <- paste0("SCALE contacts (log10, depth-scaled to WT)  |  ", STRAW_NORM)
RIDGE_SETS    <- c("wt_only"="WT-only loops", "shared"="Shared loops", "pol2_only"="Pol2-only loops")

# ── Analysis settings ──────────────────────────────────────────────────────────
WT_SRC       <- "wt_only"     # WT loops for APA / scoring; c("wt_only", "shared") for the full WT loop set
REF_LFC      <- 0.25          # absolute reference: loops on genes with |EU-seq LFC| <= this
REF_LEVEL    <- paste0("Ref |LFC| <= ", REF_LFC)
STRAT_COL    <- "eu_tertile"
STRAT_LEVELS <- c(REF_LEVEL, TERT_LEVELS)       # the reference column sits left of T1
STRAT_LABEL  <- "EU-seq absolute-loss tertile"
EXPR_LEVELS  <- c("Expr low (bottom half)", "Expr high (top half)")
EXPR_TAGS    <- c("low", "high")
EXPR_COLORS  <- setNames(c("#9ecae1", "#08519c"), EXPR_LEVELS)
N_BOOT       <- 500           # bootstrap resamples (over loops) for the scoring CIs

# ══════════════════════════════════════════════════════════════════════════════
# Functions — palettes (final_APA_analysis.R)
# ══════════════════════════════════════════════════════════════════════════════
make_fill_scale <- function(palette_name, limits, name = "Fold over\nbackground") {
  switch(palette_name,
         viridis = scale_fill_viridis_c(option="viridis", limits=limits, name=name, oob=scales::squish),
         scale_fill_gradient(low="white", high="#cc0000", limits=limits, name=name, oob=scales::squish))
}
make_diff_scale <- function(palette_name, limits, name = "log2\n(treat/WT)") {
  colors <- switch(palette_name,
                   PRGn = RColorBrewer::brewer.pal(11, "PRGn"),
                   rev(RColorBrewer::brewer.pal(11, "RdBu")))
  scale_fill_gradientn(colors=colors, limits=c(-limits, limits),
                       name=name, oob=scales::squish)
}

# ══════════════════════════════════════════════════════════════════════════════
# Functions — loop loading, shared loops, contact extraction (hic_loop_lfc_analysis.R)
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
  need <- c("observed", "expectedDonut", "fdrDonut", "fdrBL", "fdrH", "fdrV", "centroid1", "centroid2")
  if (!all(need %in% colnames(df)))
    stop(basename(path), " lacks HiCCUPS columns: ", paste(setdiff(need, colnames(df)), collapse=", "),
         "\n  columns found: ", paste(colnames(df), collapse=", "))
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

# chromosome names stored in a .hic file (loop tables and .hic files can use "1" or "chr1")
hic_chrom_names <- function(hic_file) {
  ch <- as.character(strawr::readHicChroms(hic_file)$name)
  ch[!toupper(ch) %in% "ALL"]
}
hic_chr_name <- function(x, hic_chroms) {
  cand <- cbind(x, paste0("chr", x))
  out  <- apply(cand, 1, function(r) { m <- r[r %in% hic_chroms]; if (length(m)) m[1] else NA_character_ })
  out
}

extract_loop_contacts <- function(loops_df, hic_file,
                                  norm       = STRAW_NORM,
                                  resolution = STRAW_RESOLUTION,
                                  matrix     = "observed") {
  snap    <- function(x) as.integer(floor(x / resolution) * resolution)
  results <- rep(NA_real_, nrow(loops_df))
  has_centroids <- all(c("centroid1","centroid2") %in% colnames(loops_df)) &&
    !any(is.na(loops_df$centroid1)) && !any(is.na(loops_df$centroid2))
  if (!has_centroids) stop("centroid1/centroid2 missing or NA on loops_df")
  c1_pos <- snap(as.integer(loops_df$centroid1))
  c2_pos <- snap(as.integer(loops_df$centroid2))
  chr1_clean <- gsub("^chr", "", loops_df$chr1)
  chr2_clean <- gsub("^chr", "", loops_df$chr2)
  chr_pairs  <- unique(paste(chr1_clean, chr2_clean, sep="|"))
  hic_chroms <- hic_chrom_names(hic_file)
  message("  .hic chromosome names: ", paste(head(hic_chroms, 4), collapse=", "), ", ...")
  hic1 <- hic_chr_name(unique(chr1_clean), hic_chroms); names(hic1) <- unique(chr1_clean)
  hic2 <- hic_chr_name(unique(chr2_clean), hic_chroms); names(hic2) <- unique(chr2_clean)
  if (anyNA(hic1) || anyNA(hic2))
    stop("Loop chromosomes not found in ", basename(hic_file), ": ",
         paste(unique(c(names(hic1)[is.na(hic1)], names(hic2)[is.na(hic2)])), collapse=", "),
         "  (.hic has: ", paste(hic_chroms, collapse=", "), ")")
  message("Extracting ", matrix, " from ", basename(hic_file),
          "  norm=", norm, "  res=", resolution/1000, "kb  n=", nrow(loops_df))
  for (cp in chr_pairs) {
    parts  <- strsplit(cp, "\\|")[[1]]
    c1_chr <- parts[1]; c2_chr <- parts[2]
    idx    <- which(chr1_clean == c1_chr & chr2_clean == c2_chr)
    c1_min <- min(c1_pos[idx]); c1_max <- max(c1_pos[idx]) + resolution
    c2_min <- min(c2_pos[idx]); c2_max <- max(c2_pos[idx]) + resolution
    loc1   <- paste0(hic1[[c1_chr]], ":", c1_min, ":", c1_max)
    loc2   <- paste0(hic2[[c2_chr]], ":", c2_min, ":", c2_max)
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
# Functions — EU-seq stratifier and loop-to-gene linking
# ══════════════════════════════════════════════════════════════════════════════
build_eu_genes <- function(eu_in) {
  need <- c("gene_id", "lfc", "de_status", "mean_wt", "mean_polr2a")
  if (!all(need %in% colnames(eu_in)))
    stop("eu_results is missing column(s): ", paste(setdiff(need, colnames(eu_in)), collapse=", "))
  eu <- eu_in %>%
    dplyr::transmute(gene_id   = sub("\\..*", "", as.character(gene_id)),
                     eu_lfc    = as.numeric(lfc),
                     eu_abs    = 2^as.numeric(mean_polr2a) - 2^as.numeric(mean_wt),   # scaled RPKM; negative = lost
                     eu_chk    = as.numeric(mean_polr2a) - as.numeric(mean_wt),
                     eu_wt     = as.numeric(mean_wt),       # log2(RPKM + 1), WT
                     eu_pol2   = as.numeric(mean_polr2a),   # log2(RPKM + 1), POLR2A deg
                     eu_status = factor(as.character(de_status), levels=STATUS_LEVELS)) %>%
    dplyr::filter(!is.na(eu_lfc), !is.na(eu_abs))
  # self-checks: one row per gene, known status labels, status matches the |LFC| threshold
  if (anyDuplicated(eu$gene_id) > 0) stop("eu_results has duplicated gene_id")
  if (!isTRUE(all.equal(eu$eu_lfc, eu$eu_chk))) stop("eu_results$lfc differs from mean_polr2a - mean_wt")
  # absolute loss written as WT_level * (1 - 2^LFC) must equal 2^mean_wt - 2^mean_polr2a
  wt_lin <- 2^as.numeric(eu_in$mean_wt[match(eu$gene_id, sub("\\..*", "", as.character(eu_in$gene_id)))])
  if (!isTRUE(all.equal(wt_lin * (1 - 2^eu$eu_lfc), -eu$eu_abs))) stop("absolute loss does not equal WT_level * (1 - 2^LFC)")
  eu$eu_chk <- NULL
  if (anyNA(eu$eu_status)) stop("eu_results$de_status has values outside: ", paste(STATUS_LEVELS, collapse=", "))
  chk <- dplyr::case_when(eu$eu_lfc >  LFC_THRESHOLD ~ "Up",
                          eu$eu_lfc < -LFC_THRESHOLD ~ "Down",
                          TRUE                       ~ "Not DE")
  if (!all(chk == as.character(eu$eu_status)))
    stop("eu_results$de_status does not match |lfc| > ", LFC_THRESHOLD)
  cuts <- stats::quantile(eu$eu_abs, c(1/3, 2/3), na.rm=TRUE)
  eu <- eu %>%
    dplyr::mutate(eu_tertile = factor(dplyr::case_when(eu_abs <= cuts[[1]] ~ TERT_LEVELS[1],
                                                       eu_abs <= cuts[[2]] ~ TERT_LEVELS[2],
                                                       TRUE                ~ TERT_LEVELS[3]),
                                      levels=TERT_LEVELS))
  share <- as.numeric(table(eu$eu_tertile)) / nrow(eu)
  if (max(abs(share - 1/3)) > 0.02) stop("EU-seq absolute-loss tertiles are not balanced: ", paste(round(share, 3), collapse=" / "))
  attr(eu, "tertile_cuts") <- cuts
  eu
}

# Loop -> gene: any gene overlapping an anchor (gene body +/- tol_bp); the gene with the
# highest TPM is kept (repo rule). Loops with no EU-seq gene get NA strata.
link_loops_to_eu_genes <- function(loops_df, gene_gr, tol_bp = GENE_TOL_BP) {
  loops_df$.lk <- seq_len(nrow(loops_df))
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
    data.frame(.lk=queryHits(hits1), gene_idx=subjectHits(hits1)),
    data.frame(.lk=queryHits(hits2), gene_idx=subjectHits(hits2))
  ) %>%
    distinct(.lk, gene_idx) %>%
    dplyr::mutate(eu_status  = gene_gr$eu_status[gene_idx],
                  eu_tertile = gene_gr$eu_tertile[gene_idx],
                  eu_lfc     = gene_gr$eu_lfc[gene_idx],
                  eu_abs     = gene_gr$eu_abs[gene_idx],
                  TPM        = gene_gr$TPM[gene_idx],
                  symbol     = gene_gr$symbol[gene_idx])
  ann <- all_pairs %>%
    dplyr::group_by(.lk) %>%
    dplyr::slice_max(TPM, n=1, with_ties=FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(.lk, eu_status, eu_tertile, eu_lfc, eu_abs, TPM, symbol)
  out <- loops_df %>% left_join(ann, by=".lk") %>% dplyr::select(-.lk)
  out$eu_status  <- factor(out$eu_status,  levels=STATUS_LEVELS)
  out$eu_tertile <- factor(out$eu_tertile, levels=TERT_LEVELS)
  # self-check: every loop kept exactly once
  stopifnot(nrow(out) == nrow(loops_df))
  out
}

# One anchor in a gene body (no tolerance, as FINAL_EP_loop_analysis_hiccups.R), other anchor
# a classified element; loops with the same gene at both anchors are dropped.
link_loops_singleton_eu <- function(classified_loops, gene_gr, condition_label) {
  anchor1_gr <- GRanges(classified_loops$chr1, IRanges(classified_loops$start1, classified_loops$end1))
  anchor2_gr <- GRanges(classified_loops$chr2, IRanges(classified_loops$start2, classified_loops$end2))
  gene_gr_m  <- gene_gr
  has_chr_loop <- any(grepl("^chr", as.character(seqnames(anchor1_gr))))
  has_chr_gene <- any(grepl("^chr", as.character(seqnames(gene_gr))))
  if (has_chr_loop && !has_chr_gene) seqlevels(gene_gr_m) <- paste0("chr", seqlevels(gene_gr_m))
  else if (!has_chr_loop && has_chr_gene) seqlevels(gene_gr_m) <- gsub("^chr", "", seqlevels(gene_gr_m))
  hits1 <- findOverlaps(anchor1_gr, gene_gr_m)
  hits2 <- findOverlaps(anchor2_gr, gene_gr_m)
  df1 <- data.frame(loop_idx=queryHits(hits1), gene_idx=subjectHits(hits1),
                    other_anchor_ccre=classified_loops$anchor2_ccre[queryHits(hits1)],
                    stringsAsFactors=FALSE)
  df2 <- data.frame(loop_idx=queryHits(hits2), gene_idx=subjectHits(hits2),
                    other_anchor_ccre=classified_loops$anchor1_ccre[queryHits(hits2)],
                    stringsAsFactors=FALSE)
  result <- bind_rows(df1, df2) %>%
    dplyr::mutate(eu_status  = gene_gr_m$eu_status[gene_idx],
                  eu_tertile = gene_gr_m$eu_tertile[gene_idx],
                  eu_lfc     = gene_gr_m$eu_lfc[gene_idx],
                  eu_abs     = gene_gr_m$eu_abs[gene_idx],
                  symbol     = gene_gr_m$symbol[gene_idx],
                  condition  = condition_label) %>%
    filter(!is.na(other_anchor_ccre), other_anchor_ccre != "Unannotated") %>%
    dplyr::group_by(loop_idx, gene_idx) %>% dplyr::filter(dplyr::n() == 1) %>% dplyr::ungroup()
  message(condition_label, ": ", nrow(result), " singleton loop-gene pairs at EU-seq genes")
  result
}

# ══════════════════════════════════════════════════════════════════════════════
# Functions — APA (final_APA_analysis.R)
# ══════════════════════════════════════════════════════════════════════════════
loops_to_gi <- function(loops_df, resolution = RESOLUTION, style = "ENSEMBL") {
  snap <- function(x) floor(x / resolution) * resolution
  pos1 <- snap(loops_df$centroid1)
  pos2 <- snap(loops_df$centroid2)
  gi <- GInteractions(
    anchor1    = GRanges(loops_df$chr1, IRanges(pos1, pos1 + resolution - 1)),
    anchor2    = GRanges(loops_df$chr2, IRanges(pos2, pos2 + resolution - 1)),
    eu_status  = as.character(loops_df$eu_status),
    eu_tertile = as.character(loops_df$eu_tertile),
    group      = as.character(loops_df$group)
  )
  GenomeInfoDb::seqlevelsStyle(gi) <- style
  gi
}

extract_apa_matrices_by <- function(loops_df, hic_file, group_col, group_levels,
                                    resolution = RESOLUTION, buffer = BUFFER,
                                    norm = NORM, min_loops = 10) {
  hic_style <- if (any(grepl("^chr", hic_chrom_names(hic_file)))) "UCSC" else "ENSEMBL"
  gi_all   <- loops_to_gi(loops_df, resolution, style = hic_style)
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

aggregate_apa <- function(mat_list, scale_factor = 1.0) {
  lapply(mat_list, function(mats) {
    if (is.null(mats)) return(NULL)
    cnt      <- counts(mats)
    raw_mean <- apply(cnt[,,,1, drop=FALSE], c(1,2), mean, na.rm=TRUE)
    raw_mean * scale_factor
  })
}

# APA of one loop type split by stratum: each stratum's loops at the Wildtype loop
# positions, pulled from the Wildtype and the POLR2A-degraded .hic.
run_apa_strata <- function(loops_df, group_value, strat_col, strat_levels) {
  sub <- loops_df %>% filter(group == group_value, !is.na(.data[[strat_col]]))
  message(group_value, " loops with an EU-seq stratum: ", nrow(sub))
  raw <- lapply(CONDITIONS, function(cd)
    extract_apa_matrices_by(sub, cd$hic, group_col=strat_col, group_levels=strat_levels))
  n_loops <- vapply(raw$wt, function(m) dim(counts(m))[3], numeric(1))
  # self-check: the same loops are pulled from both .hic files
  stopifnot(identical(n_loops, vapply(raw$pol2, function(m) dim(counts(m))[3], numeric(1))))
  list(wt   = aggregate_apa(raw$wt,   DEPTH_SCALE$wt),
       pol2 = aggregate_apa(raw$pol2, DEPTH_SCALE$pol2),
       n    = n_loops)
}

# 3-row grid per loop type (matched-grid layout of final_APA_analysis.R): columns are the strata;
# row 1 Wildtype, row 2 POLR2A-degraded (one shared fold-enrichment scale), row 3 difference.
# Panel numbers: APA = patch / corner background (center pixel, 3x3, 5x5); difference panels give
# log2(POLR2A deg / WT) of APA 3x3, APA 5x5, the whole 15x15 grid and the center pixel.
plot_apa_strata_grid <- function(apa, group_value, strat_levels,
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
  empty_panel <- function() ggplot() + theme_void() +
    theme(plot.background=element_rect(fill="gray95", color=NA))
  make_fold_panel <- function(mat_norm, global_min, global_max, raw) {
    if (is.null(mat_norm)) return(empty_panel())
    n <- nrow(mat_norm); mid <- ceiling(n/2)
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
    mt <- apa_pixel_metrics(raw, corner_size)      # APA = patch / corner background
    p + annotate("text", x=-BUFFER+0.5, y=BUFFER+0.4, vjust=1, hjust=0, size=2.3,
                 fontface="bold", color="white",
                 label=sprintf("APA center=%.2f\nAPA 3x3=%.2f\nAPA 5x5=%.2f",
                               mt[["center"]]/mt[["bg"]], mt[["p3"]]/mt[["bg"]], mt[["p5"]]/mt[["bg"]]))
  }
  make_diff_panel <- function(mat_diff, limit, raw_wt, raw_pol2) {
    if (is.null(mat_diff)) return(empty_panel())
    n <- nrow(mat_diff); mid <- ceiling(n/2)
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
    d <- apa_delta(raw_wt, raw_pol2, corner_size)   # log2(POLR2A deg / WT) of the reported metrics
    p + annotate("text", x=-BUFFER+0.5, y=BUFFER+0.4, vjust=1, hjust=0, size=2.3, fontface="bold",
                 label=sprintf("Δlog2 APA 3x3=%.2f\nΔlog2 APA 5x5=%.2f\nΔlog2 grid=%.2f\nΔlog2 center=%.2f",
                               d[["apa3"]], d[["apa5"]], d[["grid"]], d[["center"]]))
  }
  
  col_names <- strat_levels
  wt_mats   <- setNames(lapply(col_names, function(n) normalize_to_bg(apa$wt[[n]])),   col_names)
  pol2_mats <- setNames(lapply(col_names, function(n) normalize_to_bg(apa$pol2[[n]])), col_names)
  diff_mats <- setNames(lapply(col_names, function(n) {
    if (is.null(wt_mats[[n]]) || is.null(pol2_mats[[n]])) return(NULL)
    log2((pol2_mats[[n]] + pseudocount) / (wt_mats[[n]] + pseudocount))
  }), col_names)
  
  fold_vals <- unlist(lapply(c(wt_mats, pol2_mats), function(m) if (!is.null(m)) as.vector(m)))
  if (length(fold_vals) == 0) { message("No APA matrices for ", group_value); return(invisible(NULL)) }
  fold_max  <- quantile(fold_vals, 0.98, na.rm=TRUE)
  fold_min  <- 1.0
  diff_vals <- unlist(lapply(diff_mats, function(m) if (!is.null(m)) as.vector(m)))
  diff_lim  <- quantile(abs(diff_vals[!is.na(diff_vals)]), 0.98, na.rm=TRUE)
  
  n_lab <- ifelse(col_names %in% names(apa$n),
                  paste0("N = ", apa$n[col_names]), "N < 10")
  header_plots  <- lapply(seq_along(col_names), function(i)
    ggplot() + annotate("text", x=0.5, y=0.5, label=paste0(col_names[i], ", ", n_lab[i]),
                        size=3.8, fontface="bold") + theme_void())
  bin_hdr_plots <- lapply(col_names, function(n)
    ggplot() + annotate("text", x=0.5, y=0.5, label=group_value, color=bin_color,
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
      make_diff_scale(diff_palette, limits=limit, name="Δlog2\nenrichment") +
      theme(legend.position="right", legend.key.height=unit(1.4,"cm"),
            legend.key.width=unit(0.3,"cm"), legend.text=element_text(size=8),
            legend.title=element_text(size=9)))
  
  nc <- length(col_names)
  row_wt_panels   <- lapply(col_names, function(n) make_fold_panel(wt_mats[[n]],   fold_min, fold_max, apa$wt[[n]]))
  row_pol2_panels <- lapply(col_names, function(n) make_fold_panel(pol2_mats[[n]], fold_min, fold_max, apa$pol2[[n]]))
  row_diff_panels <- lapply(col_names, function(n) make_diff_panel(diff_mats[[n]], diff_lim, apa$wt[[n]], apa$pol2[[n]]))
  row_wt   <- cowplot::plot_grid(plotlist=row_wt_panels,   ncol=nc)
  row_pol2 <- cowplot::plot_grid(plotlist=row_pol2_panels, ncol=nc)
  row_diff <- cowplot::plot_grid(cowplot::plot_grid(plotlist=row_diff_panels, ncol=nc),
                                 build_diff_legend(diff_lim), ncol=2, rel_widths=c(1, 0.15))
  rows_1_2_body <- cowplot::plot_grid(row_wt, row_pol2, ncol=1)
  rows_1_2      <- cowplot::plot_grid(rows_1_2_body, build_fold_legend(fold_min, fold_max),
                                      ncol=2, rel_widths=c(1, 0.15))
  header_row  <- cowplot::plot_grid(cowplot::plot_grid(plotlist=header_plots,  ncol=nc),
                                    ggplot()+theme_void(), ncol=2, rel_widths=c(1, 0.15))
  bin_hdr_row <- cowplot::plot_grid(cowplot::plot_grid(plotlist=bin_hdr_plots, ncol=nc),
                                    ggplot()+theme_void(), ncol=2, rel_widths=c(1, 0.15))
  composite <- cowplot::plot_grid(header_row, bin_hdr_row, rows_1_2, row_diff,
                                  ncol=1, rel_heights=c(0.08, 0.08, 2, 1))
  if (!is.null(title)) {
    title_grob <- cowplot::ggdraw() +
      cowplot::draw_label(title, fontface="bold", size=13, x=0, hjust=0)
    composite <- cowplot::plot_grid(title_grob, composite, ncol=1, rel_heights=c(0.06, 1))
  }
  composite
}

# ══════════════════════════════════════════════════════════════════════════════
# Functions — absolute reference, APA pixel metrics, scoring
# ══════════════════════════════════════════════════════════════════════════════
# Absolute reference: loops on genes with |EU-seq LFC| <= REF_LFC are copied into an extra stratum
# (REF_LEVEL) next to the tertile loops; the copies keep every other column. eu_tertile becomes character.
add_ref <- function(df) {
  ref <- df %>% dplyr::filter(!is.na(eu_lfc), abs(eu_lfc) <= REF_LFC)
  dplyr::bind_rows(
    ref %>% dplyr::mutate(eu_tertile = REF_LEVEL),
    df  %>% dplyr::filter(!is.na(eu_tertile)) %>% dplyr::mutate(eu_tertile = as.character(eu_tertile)))
}

# Pixel metrics of a mean APA matrix (counts): center pixel, 3x3 and 5x5 patch means, whole-grid
# mean, and the corner background (mean of the four corner_size x corner_size blocks, as the APA panels)
apa_pixel_metrics <- function(m, corner_size = 3) {
  n <- nrow(m); mid <- ceiling(n / 2); cs <- corner_size
  patch <- function(h) mean(m[(mid - h):(mid + h), (mid - h):(mid + h)], na.rm = TRUE)
  bg <- mean(c(m[1:cs, 1:cs], m[1:cs, (n-cs+1):n], m[(n-cs+1):n, 1:cs], m[(n-cs+1):n, (n-cs+1):n]), na.rm = TRUE)
  c(center = m[mid, mid], p3 = patch(1), p5 = patch(2), grid = mean(m, na.rm = TRUE), bg = bg)
}
# log2(POLR2A deg / WT) of the four reported metrics for one pair of mean APA matrices
apa_delta <- function(wt, pol2, corner_size = 3) {
  w <- apa_pixel_metrics(wt, corner_size); p <- apa_pixel_metrics(pol2, corner_size)
  c(apa3   = log2((p[["p3"]] / p[["bg"]]) / (w[["p3"]] / w[["bg"]])),
    apa5   = log2((p[["p5"]] / p[["bg"]]) / (w[["p5"]] / w[["bg"]])),
    grid   = log2(p[["grid"]] / w[["grid"]]),
    center = log2(p[["center"]] / w[["center"]]))
}

# Per-loop APA matrices (loops x 225 pixels), pulled once from each .hic and depth-scaled as the APA panels
get_loop_matrices <- function(loops_df) {
  out <- lapply(names(CONDITIONS), function(nm) {
    m   <- extract_apa_matrices_by(loops_df, CONDITIONS[[nm]]$hic, "group", "WT loops")[["WT loops"]]
    cnt <- counts(m)[, , , 1, drop = FALSE]
    stopifnot(dim(cnt)[3] == nrow(loops_df))      # no loop lost to the diagonal filter, same order
    matrix(aperm(cnt[, , , 1], c(3, 1, 2)), nrow = dim(cnt)[3]) * DEPTH_SCALE[[nm]]
  })
  setNames(out, names(CONDITIONS))
}

SC_NB <- 2 * BUFFER + 1; SC_MID <- BUFFER + 1; SC_CS <- 3
SC_PIX <- matrix(seq_len(SC_NB * SC_NB), SC_NB, SC_NB)
sc_patch  <- function(h) as.vector(SC_PIX[(SC_MID - h):(SC_MID + h), (SC_MID - h):(SC_MID + h)])
sc_corner <- c(SC_PIX[1:SC_CS, 1:SC_CS], SC_PIX[1:SC_CS, (SC_NB-SC_CS+1):SC_NB],
               SC_PIX[(SC_NB-SC_CS+1):SC_NB, 1:SC_CS], SC_PIX[(SC_NB-SC_CS+1):SC_NB, (SC_NB-SC_CS+1):SC_NB])
# per-loop sums over a pixel set; W = WT, P = POLR2A deg
scoring_inputs <- function(pl) {
  rs <- function(M, idx) rowSums(M[, idx, drop = FALSE], na.rm = TRUE)
  list(W1 = rs(pl$wt, sc_patch(0)), P1 = rs(pl$pol2, sc_patch(0)),
       W3 = rs(pl$wt, sc_patch(1)), P3 = rs(pl$pol2, sc_patch(1)),
       W5 = rs(pl$wt, sc_patch(2)), P5 = rs(pl$pol2, sc_patch(2)),
       CW = rs(pl$wt, sc_corner),   CP = rs(pl$pol2, sc_corner),
       AW = rs(pl$wt, seq_len(SC_NB * SC_NB)), AP = rs(pl$pol2, seq_len(SC_NB * SC_NB)))
}
lr     <- function(p, w) log2(sum(p) / sum(w))
apa_lr <- function(p, cp, w, cw) log2((sum(p) / sum(cp)) / (sum(w) / sum(cw)))   # patch over background, per file
# the four reported metrics; each takes the per-loop sums L and a vector of loop indices i
SCORE_METRICS <- list(
  "APA-style 3x3 patch (corner background)" = function(L, i) apa_lr(L$P3[i], L$CP[i], L$W3[i], L$CW[i]),
  "APA-style 5x5 patch (corner background)" = function(L, i) apa_lr(L$P5[i], L$CP[i], L$W5[i], L$CW[i]),
  "Whole grid (15x15), raw"                 = function(L, i) lr(L$AP[i], L$AW[i]),
  "Center pixel, raw"                       = function(L, i) lr(L$P1[i], L$W1[i]))

# Loop index sets per stratum (reference + tertiles) of the loop table d, optionally within a mask
make_groups <- function(d, mask = rep(TRUE, nrow(d))) {
  g <- c(list(which(mask & !is.na(d$eu_lfc) & abs(d$eu_lfc) <= REF_LFC)),
         lapply(TERT_LEVELS, function(lv) which(mask & d$eu_tertile == lv)))
  setNames(g, STRAT_LEVELS)
}

# Bootstrap over loops for every metric x group; also T1 - T3 with its CI. groups with < min_n loops are NA.
run_scoring <- function(L, groups, metrics = SCORE_METRICS, n_boot = N_BOOT, min_n = 10,
                        diff_pair = TERT_LEVELS[c(1, 3)]) {
  ok <- lengths(groups) >= min_n
  rows <- list(); drows <- list()
  for (m in names(metrics)) {
    f   <- metrics[[m]]
    est <- lapply(groups[ok], function(i) f(L, i))
    bs  <- lapply(groups[ok], function(i) replicate(n_boot, f(L, i[sample.int(length(i), replace = TRUE)])))
    pick <- function(fun) vapply(names(groups), function(g) if (ok[[g]]) fun(g) else NA_real_, numeric(1))
    rows[[m]] <- data.frame(metric = m, group = names(groups), n = unname(lengths(groups)),
                            estimate = pick(function(g) est[[g]]),
                            lo = pick(function(g) quantile(bs[[g]], 0.025, names = FALSE)),
                            hi = pick(function(g) quantile(bs[[g]], 0.975, names = FALSE)),
                            row.names = NULL)
    if (all(ok[diff_pair])) {
      d <- bs[[diff_pair[1]]] - bs[[diff_pair[2]]]
      drows[[m]] <- data.frame(metric = m, diff = est[[diff_pair[1]]] - est[[diff_pair[2]]],
                               lo = quantile(d, 0.025, names = FALSE), hi = quantile(d, 0.975, names = FALSE))
    }
  }
  list(long = dplyr::bind_rows(rows), diff = if (length(drows)) dplyr::bind_rows(drows) else NULL)
}

print_scoring <- function(sc, header) {
  message(header, "\n  log2(POLR2A deg / WT), estimate [95% bootstrap CI over loops]; T1 - T3 < 0 means T1 loses more")
  tab <- sc$long %>%
    dplyr::mutate(txt = ifelse(is.na(estimate), "n < 10", sprintf("%.2f [%.2f, %.2f]", estimate, lo, hi))) %>%
    dplyr::select(metric, group, txt) %>%
    tidyr::pivot_wider(names_from = group, values_from = txt)
  if (!is.null(sc$diff))
    tab <- dplyr::left_join(tab, sc$diff %>% dplyr::transmute(metric, `T1 - T3` = sprintf("%.2f [%.2f, %.2f]", diff, lo, hi)),
                            by = "metric")
  print(as.data.frame(tab), row.names = FALSE)
  message("  loops per stratum: ", paste(unique(paste0(sc$long$group, " = ", sc$long$n)), collapse = ", "))
}

# Scoring figure: one facet per metric; x = stratum (Ref, T1-T3). split_colors: color + dodge by a `split` column
plot_scoring <- function(long_df, title, subtitle, split_colors = NULL) {
  d <- long_df %>% dplyr::mutate(metric = factor(metric, levels = names(SCORE_METRICS)),
                                 group  = factor(group, levels = STRAT_LEVELS))
  p <- ggplot(d, aes(x = group, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60")
  if (is.null(split_colors)) {
    nn  <- d %>% dplyr::filter(metric == levels(metric)[1]) %>% dplyr::distinct(group, n)
    lab <- setNames(paste0(nn$group, "\nN = ", nn$n), as.character(nn$group))
    p <- p + geom_pointrange(aes(ymin = lo, ymax = hi), color = "#882255", size = 0.4) +
      scale_x_discrete(labels = lab, drop = FALSE)
  } else {
    d$split <- factor(d$split, levels = names(split_colors))
    p <- ggplot(d, aes(x = group, y = estimate, color = split)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.4) +
      scale_color_manual(values = split_colors, name = NULL) + scale_x_discrete(drop = FALSE)
  }
  p + facet_wrap(~ metric, nrow = 1, scales = "free_y", labeller = label_wrap_gen(24)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "log2(POLR2A deg / WT)") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8, color = "gray40"),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
          axis.text.x = element_text(size = 7), legend.position = "bottom")
}

# ══════════════════════════════════════════════════════════════════════════════
# Functions — log2FC distribution, Lost | Gained ridges, Gained / Lost bars, non-gene anchor types, SVG export
# ══════════════════════════════════════════════════════════════════════════════
# Log2FC loop distribution per stratum with a boxplot strip under each panel. type = "histogram"
# (counts) or "density" (partially transparent kernel densities, one per group). group_col /
# group_levels / group_colors choose the groups drawn within each stratum (default: loop source).
plot_loop_lfc_by_stratum_box <- function(df, strat_col, strat_levels, strat_label,
                                         type = c("histogram", "density"),
                                         group_col    = "source",
                                         group_levels = c("pol2_only", "shared", "wt_only"),
                                         group_colors = SOURCE_COLORS) {
  type <- match.arg(type)
  src_levels <- group_levels
  d <- df %>% filter(!is.na(.data[[strat_col]]), !is.na(log2fc))
  d$stratum <- factor(as.character(d[[strat_col]]), levels=strat_levels)
  d$source  <- factor(d[[group_col]], levels=src_levels)
  n_tbl <- d %>% dplyr::count(stratum, .drop=FALSE)
  lab   <- setNames(paste0(n_tbl$stratum, "\nN = ", n_tbl$n, " loops"), as.character(n_tbl$stratum))
  
  # numbers attached to each box: N and median log2FC per source within each stratum
  num_tbl <- d %>%
    dplyr::group_by(stratum, source, .drop=FALSE) %>%
    dplyr::summarise(n=dplyr::n(), med=median(log2fc), .groups="drop") %>%
    dplyr::filter(n > 0) %>%
    dplyr::mutate(num_label=sprintf("N=%d, median=%.2f", n, med))
  
  top_theme <- theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold", size=13),
          plot.subtitle=element_text(size=8, color="gray40"),
          strip.background=element_blank(), strip.text=element_text(face="bold"),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())
  top_labs <- labs(title="Log2FC Loop Distribution",
                   subtitle=paste0("Log2((POLR2A deg + ", PSEUDOCOUNT, ")/(WT + ", PSEUDOCOUNT,
                                   "))  |  Threshold = ", LOG2FC_THRESHOLD, "  |  loops by ", strat_label),
                   x=NULL, y=if (type == "histogram") "Count" else "Density")
  thr_lines <- list(
    geom_vline(xintercept=-LOG2FC_THRESHOLD, linetype="dashed", color="#4393c3", linewidth=0.8),
    geom_vline(xintercept= LOG2FC_THRESHOLD, linetype="dashed", color="#882255", linewidth=0.8))
  p_hist <- if (type == "histogram") {
    ggplot(d, aes(x=log2fc, fill=source)) +
      geom_histogram(bins=100, alpha=0.8, position="identity") + thr_lines +
      facet_wrap(~ stratum, nrow=1, scales="free_y", labeller=labeller(stratum=lab)) +
      scale_x_continuous(limits=c(-4, 4)) +
      scale_fill_manual(values=group_colors, breaks=src_levels, name=NULL) +
      top_labs + top_theme
  } else {
    ggplot(d, aes(x=log2fc, fill=source, color=source)) +
      geom_density(alpha=0.35, linewidth=0.8) + thr_lines +
      facet_wrap(~ stratum, nrow=1, labeller=labeller(stratum=lab)) +
      coord_cartesian(xlim=c(-4, 4)) +
      scale_fill_manual(values=group_colors, breaks=src_levels, name=NULL) +
      scale_color_manual(values=group_colors, breaks=src_levels, name=NULL) +
      top_labs + top_theme
  }
  
  p_box <- ggplot(d, aes(x=log2fc, y=source, color=source, fill=source)) +
    geom_boxplot(alpha=0.35, linewidth=0.5, outlier.size=0.3, outlier.alpha=0.25, width=0.6) +
    geom_vline(xintercept=-LOG2FC_THRESHOLD, linetype="dashed", color="#4393c3", linewidth=0.5) +
    geom_vline(xintercept= LOG2FC_THRESHOLD, linetype="dashed", color="#882255", linewidth=0.5) +
    geom_text(data=num_tbl, aes(x=3.95, y=source, label=num_label, color=source),
              inherit.aes=FALSE, hjust=1, vjust=-1.1, size=2.6, show.legend=FALSE) +
    facet_wrap(~ stratum, nrow=1) +
    scale_color_manual(values=group_colors, guide="none") +
    scale_fill_manual(values=group_colors, guide="none") +
    scale_y_discrete(limits=src_levels, drop=FALSE) +
    coord_cartesian(xlim=c(-4, 4)) +
    labs(x="log2FC (POLR2A deg / WT)", y=NULL) +
    theme_classic(base_size=12) +
    theme(strip.background=element_blank(), strip.text=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())
  
  legend <- cowplot::get_legend(p_hist)
  stacked <- cowplot::plot_grid(p_hist + theme(legend.position="none"), p_box,
                                ncol=1, align="v", axis="lr", rel_heights=c(5, 1.8))
  cowplot::plot_grid(stacked, legend, nrow=1, rel_widths=c(1, 0.16))
}

# F: Lost (left) vs Gained (right) contact distributions per stratum for one loop source.
# Wildtype dashed, POLR2A-degraded solid. The label of each ridge gives the number of loops
# and their median log2FC (POLR2A deg / WT), not a contact value.
plot_lost_gained_ridge <- function(df, source_filter, strat_col, strat_levels,
                                   title, subtitle, wt_color = "#3d3d3d") {
  class_levels <- c("Lost", "Gained")
  key_levels   <- as.vector(rbind(paste0("Lost — ",   strat_levels),
                                  paste0("Gained — ", strat_levels)))
  base <- df %>%
    dplyr::filter(source == source_filter, loop_class %in% class_levels,
                  !is.na(.data[[strat_col]]), !is.na(log2fc)) %>%
    dplyr::mutate(class_grp  = factor(as.character(loop_class), levels=class_levels),
                  ridge_key  = factor(paste0(loop_class, " — ", .data[[strat_col]]), levels=rev(key_levels)))
  if (nrow(base) == 0) {
    message("No loops matched source='", source_filter, "' with class Gained/Lost x ", strat_col)
    return(invisible(NULL))
  }
  message(source_filter, " Lost / Gained loops by ", strat_col, ":"); print(table(base$class_grp, base[[strat_col]]))
  long_df <- base %>%
    dplyr::select(class_grp, ridge_key, wt_contacts, pol2_contacts) %>%
    tidyr::pivot_longer(cols=c(wt_contacts, pol2_contacts), names_to="cond", values_to="contacts") %>%
    dplyr::filter(!is.na(contacts), contacts > 0) %>%
    dplyr::mutate(cond_label = ifelse(cond == "wt_contacts", "WT", "POLR2A deg"))
  cnt <- long_df %>% dplyr::group_by(ridge_key, cond_label) %>% dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
    tidyr::pivot_wider(names_from=cond_label, values_from=n, values_fill=0)
  # a ridge needs >= 3 loops in each condition
  keep <- as.character(cnt$ridge_key[cnt$WT >= 3 & cnt$`POLR2A deg` >= 3])
  dropped <- setdiff(as.character(cnt$ridge_key), keep)
  if (length(dropped) > 0) message("  Dropped (< 3 loops): ", paste(dropped, collapse="; "))
  long_df <- long_df %>% dplyr::filter(as.character(ridge_key) %in% keep)
  if (nrow(long_df) == 0) return(invisible(NULL))
  lab_tbl <- base %>%
    dplyr::group_by(ridge_key) %>%
    dplyr::summarise(n=dplyr::n(), med_lfc=median(log2fc), .groups="drop") %>%
    dplyr::filter(as.character(ridge_key) %in% keep) %>%
    dplyr::mutate(axis_label=sprintf("%s (n=%d, median LFC=%.2f)",
                                     sub("^[^—]*— ", "", as.character(ridge_key)), n, med_lfc))
  label_lookup <- setNames(lab_tbl$axis_label, as.character(lab_tbl$ridge_key))
  class_cols <- CLASS_COLORS[class_levels]
  wt_df   <- long_df %>% dplyr::filter(cond_label == "WT")
  pol2_df <- long_df %>% dplyr::filter(cond_label == "POLR2A deg")
  ggplot(mapping=aes(x=log10(contacts), y=ridge_key)) +
    ggridges::geom_density_ridges(
      data=wt_df, fill=scales::alpha(wt_color, 0.60), color=scales::alpha("black", 0.8),
      linewidth=0.7, linetype="dashed", scale=0.9, quantile_lines=TRUE, quantiles=2) +
    ggridges::geom_density_ridges(
      data=pol2_df, aes(fill=class_grp, color=class_grp),
      alpha=0.70, linewidth=0.6, scale=0.9, quantile_lines=TRUE, quantiles=2) +
    facet_wrap(~ class_grp, nrow=1, scales="free_y") +
    scale_fill_manual(values=scales::alpha(class_cols, 0.75), guide="none") +
    scale_color_manual(values=class_cols, guide="none") +
    scale_y_discrete(labels=function(x) label_lookup[x]) +
    scale_x_continuous(name=SCORE_LABEL) +
    labs(title=title, subtitle=subtitle, y=NULL,
         caption="Charcoal dashed = Wildtype  |  Color = POLR2A deg  |  Vertical lines = median contact of each distribution  |  Axis labels: n loops and median log2FC") +
    theme_ridges(font_size=11) +
    theme(plot.title    = element_text(face="bold", size=13),
          plot.subtitle = element_text(size=8, color="gray40"),
          plot.caption  = element_text(size=8, color="gray40"),
          axis.text.y   = element_text(size=9),
          strip.text    = element_text(face="bold"))
}

# Gained vs Lost loop counts per stratum (grouped bars, Lost left / Gained right)
make_gain_loss_bar <- function(df, strat_col, strat_levels, title, subtitle) {
  class_levels <- c("Lost", "Gained")
  plot_df <- df %>%
    dplyr::filter(loop_class %in% class_levels, !is.na(.data[[strat_col]])) %>%
    dplyr::mutate(stratum    = factor(as.character(.data[[strat_col]]), levels=strat_levels),
                  class_grp  = factor(as.character(loop_class), levels=class_levels)) %>%
    dplyr::count(stratum, class_grp, name="n_loops") %>%
    tidyr::complete(stratum=factor(strat_levels, levels=strat_levels),
                    class_grp=factor(class_levels, levels=class_levels), fill=list(n_loops=0))
  message(title); print(tidyr::pivot_wider(plot_df, names_from=class_grp, values_from=n_loops))
  ggplot(plot_df, aes(x=stratum, y=n_loops, fill=class_grp)) +
    geom_col(position=position_dodge(width=0.75), width=0.7, color="white", linewidth=0.3, alpha=0.92) +
    geom_text(aes(label=n_loops), position=position_dodge(width=0.75), vjust=-0.6,
              size=3, color="gray30") +
    scale_fill_manual(values=CLASS_COLORS[class_levels], name=NULL) +
    scale_y_continuous(name="Number of loops", expand=expansion(mult=c(0, 0.15))) +
    labs(title=title, subtitle=subtitle, x=NULL) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold", size=13),
          plot.subtitle=element_text(size=9, color="gray40"),
          legend.position="bottom")
}

# E: non-gene anchor type counts (Wildtype loops vs POLR2A-degraded loops) at one stratum
make_singleton_count_plot_eu <- function(all_pairs, strat_col, stratum, ccre_types,
                                         title, strat_label, min_n = 5) {
  cond_levels <- levels(all_pairs$condition)
  plot_df <- all_pairs %>%
    filter(.data[[strat_col]] == stratum, other_anchor_ccre %in% ccre_types) %>%
    distinct(loop_idx, gene_idx, other_anchor_ccre, condition) %>%
    dplyr::count(other_anchor_ccre, condition, name="n_loops") %>%
    tidyr::complete(other_anchor_ccre=ccre_types, condition=cond_levels, fill=list(n_loops=0)) %>%
    dplyr::group_by(other_anchor_ccre) %>% dplyr::filter(max(n_loops) >= min_n) %>% dplyr::ungroup() %>%
    dplyr::mutate(other_anchor_ccre=factor(other_anchor_ccre, levels=ccre_types),
                  condition=factor(condition, levels=cond_levels),
                  fill_key=factor(paste0(condition, "__", other_anchor_ccre),
                                  levels=c(paste0(cond_levels[1], "__", ccre_types),
                                           paste0(cond_levels[2], "__", ccre_types))))
  if (nrow(plot_df) == 0) return(ggplot() + theme_void() + labs(title=title))
  type_order <- plot_df %>% dplyr::group_by(other_anchor_ccre) %>%
    dplyr::summarise(total=sum(n_loops), .groups="drop") %>% arrange(desc(total))
  plot_df <- plot_df %>%
    dplyr::mutate(other_anchor_ccre=factor(other_anchor_ccre, levels=type_order$other_anchor_ccre))
  fill_values <- c(
    setNames(rep("gray70", length(ccre_types)), paste0(cond_levels[1], "__", ccre_types)),
    setNames(single_anchor_colors[ccre_types],  paste0(cond_levels[2], "__", ccre_types))
  )
  ggplot(plot_df, aes(x=other_anchor_ccre, y=n_loops, fill=fill_key, group=condition)) +
    geom_col(position=position_dodge(width=0.7), width=0.65, color="white", linewidth=0.3, alpha=0.92) +
    geom_text(aes(label=ifelse(n_loops > 0, n_loops, "")),
              position=position_dodge(width=0.7), vjust=-0.6, size=2.8, color="gray30") +
    scale_fill_manual(values=fill_values,
                      breaks=c(paste0(cond_levels[1], "__", ccre_types[1]),
                               paste0(cond_levels[2], "__", ccre_types[1])),
                      labels=c("WT (grey)", "POLR2A deg (colored)"), name=NULL) +
    scale_y_continuous(name="Number of loops", expand=expansion(mult=c(0, 0.20))) +
    labs(title=title,
         subtitle=paste0("Left=Wildtype (grey)  |  Right=POLR2A deg (colored)  |  ",
                         "One anchor in gene body  |  narrowPeak + ATAC  |  ", strat_label),
         x=NULL) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold", size=13),
          plot.subtitle=element_text(size=9, color="gray40"),
          axis.text.x=element_text(size=10), legend.position="bottom",
          plot.margin=margin(8, 8, 5, 8))
}

# ── SVG export: every plot is printed and also saved in PLOT_DIR ──────────────
save_plot_svg <- function(p, strat_name, panel, w, h) {
  if (is.null(p)) return(invisible(NULL))
  print(p)
  f <- file.path(PLOT_DIR, paste0(DATA_TAG, "_euseq-", strat_name, "_", panel, ".svg"))
  svg(f, width=w, height=h); print(p); dev.off()
  message("  saved ", f)
  invisible(f)
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. EU-seq stratifier -> gene ranges
# ══════════════════════════════════════════════════════════════════════════════
message("\n── EU-seq stratifier ───────────────────────────────────────────────────")
if (exists("eu_results", inherits = TRUE)) {
  eu_in <- get("eu_results", inherits = TRUE)
  message("Using eu_results from this session (", nrow(eu_in), " genes)")
} else if (!is.null(EU_RESULTS_FILE)) {
  check_files(EU_RESULTS_FILE)
  eu_in <- read_csv(EU_RESULTS_FILE, show_col_types = FALSE)
  message("Read eu_results from ", EU_RESULTS_FILE, " (", nrow(eu_in), " genes)")
} else {
  stop("eu_results is not in this session: run EU_seq_zscore_analysis.R up to the eu_results ",
       "table, or set EU_RESULTS_FILE to a csv of it.")
}
eu_genes <- build_eu_genes(eu_in)
message("EU-seq genes: ", nrow(eu_genes))
message("Status counts:");  print(table(eu_genes$eu_status))
message("Tertile cut points (EU-seq absolute change, scaled RPKM; negative = lost): ",
        paste(round(attr(eu_genes, "tertile_cuts"), 3), collapse=" | "))
message("Tertile counts:"); print(table(eu_genes$eu_tertile))

# Gene-level EU-seq signal per tertile. Absolute level = each tertile's share of the TOTAL EU-seq signal of all EU-seq genes
# (scaled RPKM = 2^x - 1 of the log2(RPKM + 1) means in eu_results), as % of the total WT signal, for WT and POLR2A deg
# (WT shares sum to 100%; POLR2A deg shares sum to the % of total transcription that remains). 95% bootstrap CI over genes.
gene_tot <- eu_genes %>%
  dplyr::mutate(stratum = factor(eu_tertile, levels = TERT_LEVELS), wt_rpkm = 2^eu_wt - 1, resid_rpkm = 2^eu_pol2 - 1)
stopifnot(all(gene_tot$wt_rpkm > 0), !anyNA(gene_tot$stratum))
share_fun <- function(d) {
  tot <- sum(d$wt_rpkm)
  data.frame(stratum = factor(TERT_LEVELS, levels = TERT_LEVELS),
             Wildtype     = vapply(TERT_LEVELS, function(lv) 100 * sum(d$wt_rpkm[d$stratum == lv])    / tot, numeric(1)),
             `POLR2A deg` = vapply(TERT_LEVELS, function(lv) 100 * sum(d$resid_rpkm[d$stratum == lv]) / tot, numeric(1)),
             check.names = FALSE, row.names = NULL)
}
set.seed(1)
share_est <- share_fun(gene_tot) %>% tidyr::pivot_longer(c(Wildtype, `POLR2A deg`), names_to = "series", values_to = "est")
share_bs  <- dplyr::bind_rows(lapply(seq_len(N_BOOT), function(b) share_fun(gene_tot[sample.int(nrow(gene_tot), replace = TRUE), ]))) %>%
  tidyr::pivot_longer(c(Wildtype, `POLR2A deg`), names_to = "series", values_to = "v") %>%
  dplyr::group_by(stratum, series) %>%
  dplyr::summarise(lo = quantile(v, 0.025, names = FALSE), hi = quantile(v, 0.975, names = FALSE), .groups = "drop")
gene_share <- dplyr::left_join(share_est, share_bs, by = c("stratum", "series")) %>%
  dplyr::mutate(series = factor(series, levels = c("Wildtype", "POLR2A deg")))
stopifnot(abs(sum(gene_share$est[gene_share$series == "Wildtype"]) - 100) < 1e-6)
message("Share of total WT EU-seq signal (all EU-seq genes), % (95% bootstrap CI over genes):"); print(as.data.frame(gene_share), row.names = FALSE)
message("Total EU-seq signal remaining after degradation: ", round(sum(gene_share$est[gene_share$series == "POLR2A deg"]), 1), "% of WT")
gene_lv <- gene_tot %>% dplyr::group_by(stratum) %>%
  dplyr::summarise(N = dplyr::n(), Wildtype = median(wt_rpkm), `POLR2A deg` = median(resid_rpkm), LFC = median(eu_lfc), .groups = "drop")
message("Median EU-seq level (scaled RPKM) and LFC per tertile (genes):"); print(as.data.frame(gene_lv), row.names = FALSE)
P0_NAME <- "Median EU-seq level (scaled RPKM)"
P1_NAME <- "Absolute residual EU-seq signal (% of total WT signal, all genes)"
P2_NAME <- "Median EU-seq LFC (POLR2A deg / WT)"
gene_lv_long <- dplyr::bind_rows(
  gene_lv %>% tidyr::pivot_longer(c(Wildtype, `POLR2A deg`), names_to = "series", values_to = "value") %>%
    dplyr::transmute(stratum, series, value, panel = P0_NAME),
  gene_share %>% dplyr::filter(series == "POLR2A deg") %>% dplyr::transmute(stratum, series, value = est, panel = P1_NAME),
  gene_lv %>% dplyr::transmute(stratum, series = "LFC", value = LFC, panel = P2_NAME)) %>%
  dplyr::mutate(series = factor(series, levels = c("Wildtype", "POLR2A deg", "LFC")),
                panel  = factor(panel, levels = c(P0_NAME, P1_NAME, P2_NAME)))
gene_lab <- setNames(paste0(gene_lv$stratum, "\nN = ", gene_lv$N, " genes"), as.character(gene_lv$stratum))
p_gene_lv <- ggplot(gene_lv_long, aes(x = stratum, y = value, fill = series)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(data = gene_share %>% dplyr::filter(series == "POLR2A deg") %>%
                  dplyr::mutate(panel = factor(P1_NAME, levels = c(P0_NAME, P1_NAME, P2_NAME))),
                aes(x = stratum, ymin = lo, ymax = hi), width = 0.25, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, color = "gray40") +
  scale_fill_manual(values = c(Wildtype = "gray60", `POLR2A deg` = "#882255", LFC = "#4393c3"), name = NULL, drop = TRUE) +
  scale_x_discrete(labels = gene_lab) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  labs(title = paste0("EU-seq signal and LFC of genes by ", STRAT_LABEL),
       subtitle = paste0(nrow(eu_genes), " EU-seq genes  |  left: median level per gene  |  middle: POLR2A deg signal of the tertile as % of the total WT EU-seq signal, ",
                         "error bars = 95% bootstrap CI over genes  |  right: median LFC per tertile"),
       x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8, color = "gray40"),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 8), legend.position = "bottom")
save_plot_svg(p_gene_lv, "tertile", "gene-level_euseq-level-lfc-bars", w = 16, h = 5)

expr <- read_tsv(EXPR_FILE, show_col_types = FALSE) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%
  dplyr::select(gene_id, TPM) %>% distinct(gene_id, .keep_all = TRUE)
# gene coordinates and symbols from the same Ensembl GTF the EU-seq script used (all Ensembl genes)
if (!requireNamespace("rtracklayer", quietly = TRUE)) stop("Package not installed: rtracklayer")
gtf_genes <- as.data.frame(rtracklayer::import(GTF_FILE, feature.type = "gene"))
need_cols <- c("gene_id", "gene_name", "seqnames", "start", "end")
if (!all(need_cols %in% colnames(gtf_genes)))
  stop("GTF gene records are missing column(s): ", paste(setdiff(need_cols, colnames(gtf_genes)), collapse=", "))
gene_tbl <- gtf_genes %>%
  dplyr::transmute(gene_id        = sub("\\..*", "", gene_id),
                   symbol         = gene_name,
                   gene_seq_start = start,
                   gene_seq_end   = end,
                   seq_name       = as.character(seqnames))
gene_coords <- gene_tbl %>%
  dplyr::mutate(symbol = dplyr::case_when(!is.na(symbol) & symbol != "" ~ symbol, TRUE ~ gene_id)) %>%
  left_join(expr, by = "gene_id") %>%
  filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name)) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  inner_join(eu_genes, by = "gene_id")
match_frac <- nrow(gene_coords) / nrow(eu_genes)
message("EU-seq genes with coordinates: ", nrow(gene_coords), " / ", nrow(eu_genes),
        " (", round(100*match_frac, 1), "%)")
if (match_frac < 0.9) {
  message("Example EU-seq gene_ids: ", paste(head(eu_genes$gene_id, 3), collapse = ", "))
  message("Example GTF gene_ids:    ", paste(head(gene_tbl$gene_id, 3), collapse = ", "))
  message("Unmatched EU-seq gene_ids: ", paste(head(setdiff(eu_genes$gene_id, gene_tbl$gene_id), 5), collapse = ", "))
  stop("Fewer than 90% of EU-seq genes matched gene coordinates; check gene_id format")
}
gene_gr <- GRanges(
  seqnames   = as.character(gene_coords$seq_name),
  ranges     = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end),
  symbol     = gene_coords$symbol,
  gene_id    = gene_coords$gene_id,
  eu_status  = as.character(gene_coords$eu_status),
  eu_tertile = as.character(gene_coords$eu_tertile),
  eu_lfc     = gene_coords$eu_lfc,
  eu_abs     = gene_coords$eu_abs,
  TPM        = gene_coords$TPM
)
message("gene_gr (EU-seq genes): ", length(gene_gr))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Loops, union set, contacts, Gained / Stable / Lost (hic_loop_lfc_analysis.R)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Loading HiCCUPS loops (", DATA_LABEL, ") ─────────────────────────────")
loops_wt <- filter_loops(load_hiccups_bedpe(hic_wt_loops), MIN_OBS_WT) %>%
  dplyr::mutate(loop_idx=row_number(), loop_dist=abs(centroid2-centroid1))
loops_pol2 <- filter_loops(load_hiccups_bedpe(hic_pol2_loops), MIN_OBS_POL2) %>%
  dplyr::mutate(loop_idx=row_number(), loop_dist=abs(centroid2-centroid1))
message("WT loops: ", nrow(loops_wt), "   POLR2A deg loops: ", nrow(loops_pol2))

loop_overlap <- find_shared_loops(loops_wt, loops_pol2)
message("Shared: ", loop_overlap$n_shared, "  WT only: ", loop_overlap$n_a_only,
        "  pol2 only: ", loop_overlap$n_b_only)
union_loops <- bind_rows(
  loops_wt   %>% dplyr::slice(loop_overlap$shared$a_idx) %>% dplyr::mutate(source="shared"),
  loops_wt   %>% dplyr::slice(loop_overlap$a_only)       %>% dplyr::mutate(source="wt_only"),
  loops_pol2 %>% dplyr::slice(loop_overlap$b_only)       %>% dplyr::mutate(source="pol2_only")
) %>% dplyr::mutate(loop_idx=row_number())

wt_contacts_scaled   <- replace_na(extract_loop_contacts(union_loops, hic_wt_file),   0) * DEPTH_SCALE_WT
pol2_contacts_scaled <- replace_na(extract_loop_contacts(union_loops, hic_pol2_file), 0) * DEPTH_SCALE_POL2
message("WT contacts   — median: ", round(median(wt_contacts_scaled), 1),
        "  zero: ", sum(wt_contacts_scaled == 0), " / ", length(wt_contacts_scaled))
message("POLR2A contacts — median: ", round(median(pol2_contacts_scaled), 1),
        "  zero: ", sum(pol2_contacts_scaled == 0), " / ", length(pol2_contacts_scaled))
# self-check: the two .hic files must give real, different contact values at the loops
if (mean(wt_contacts_scaled > 0) < 0.5 || mean(pol2_contacts_scaled > 0) < 0.5)
  stop("Fewer than 50% of loops have contacts in the .hic files; see the 'Warning:' / 'Signal:' lines above ",
       "(norm = ", STRAW_NORM, ", res = ", STRAW_RESOLUTION, ", files = ", hic_wt_file, ", ", hic_pol2_file, ")")
if (identical(wt_contacts_scaled, pol2_contacts_scaled))
  stop("WT and POLR2A contacts are identical at every loop; check hic_wt_file / hic_pol2_file")
union_loops <- union_loops %>%
  dplyr::mutate(
    wt_contacts   = wt_contacts_scaled,
    pol2_contacts = pol2_contacts_scaled,
    log2fc        = log2((pol2_contacts + PSEUDOCOUNT) / (wt_contacts + PSEUDOCOUNT)),
    loop_class    = dplyr::case_when(
      wt_contacts   >= MIN_WT_CONTACTS   & log2fc < -LOG2FC_THRESHOLD ~ "Lost",
      pol2_contacts >= MIN_POL2_CONTACTS & log2fc >  LOG2FC_THRESHOLD ~ "Gained",
      TRUE                                                             ~ "Stable"),
    loop_class    = factor(loop_class, levels=c("Gained","Stable","Lost"))
  )
message("Loop class breakdown:"); print(table(union_loops$loop_class, union_loops$source))

union_eu <- link_loops_to_eu_genes(union_loops, gene_gr)
message("Union loops with an EU-seq gene at an anchor: ", sum(!is.na(union_eu$eu_status)),
        " / ", nrow(union_eu))
message("Loops by source x EU-seq status:");  print(table(union_eu$source, union_eu$eu_status, useNA="ifany"))
message("Loops by source x EU-seq tertile:"); print(table(union_eu$source, union_eu$eu_tertile, useNA="ifany"))
# self-check: status and tertile are assigned together (both NA or both set)
stopifnot(identical(is.na(union_eu$eu_status), is.na(union_eu$eu_tertile)))


# ══════════════════════════════════════════════════════════════════════════════
# 3. Log2FC distributions and Lost | Gained ridges per stratum (Ref + tertiles)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Log2FC distributions and Lost / Gained ridges ───────────────────────")
message("EU-seq genes with |LFC| <= ", REF_LFC, ": ", sum(abs(gene_gr$eu_lfc) <= REF_LFC, na.rm = TRUE), " / ", length(gene_gr))
union_ref <- add_ref(union_eu)
message("Loops by source x stratum:"); print(table(union_ref$source, factor(union_ref$eu_tertile, levels = STRAT_LEVELS)))

# Absolute residual transcription vs loop loss by tertile:
#   absolute residual EU-seq signal = POLR2A deg signal of the tertile as % of the total WT EU-seq signal of all EU-seq genes
#   (95% bootstrap CI over genes) and median log2FC of WT loops (95% bootstrap CI over loops)
gene_res_s <- gene_share %>% dplyr::filter(series == "POLR2A deg") %>%
  dplyr::left_join(gene_lv %>% dplyr::select(stratum, N_genes = N), by = "stratum") %>%
  dplyr::transmute(stratum, N_genes, resid = est, rlo = lo, rhi = hi)
set.seed(1)
loop_res_s <- union_eu %>%
  dplyr::filter(source %in% WT_SRC, !is.na(log2fc), !is.na(eu_tertile)) %>%
  dplyr::mutate(stratum = factor(eu_tertile, levels = TERT_LEVELS)) %>%
  dplyr::group_by(stratum) %>%
  dplyr::summarise(N_loops = dplyr::n(), lfc = median(log2fc),
                   lo = quantile(replicate(N_BOOT, median(sample(log2fc, replace = TRUE))), 0.025, names = FALSE),
                   hi = quantile(replicate(N_BOOT, median(sample(log2fc, replace = TRUE))), 0.975, names = FALSE),
                   .groups = "drop")
res_tab <- dplyr::left_join(gene_res_s, loop_res_s, by = "stratum")
message("Absolute residual EU-seq signal (genes) and loop log2FC (", paste(WT_SRC, collapse = " + "), " loops) per tertile:")
print(as.data.frame(res_tab), row.names = FALSE)

res_lab_g <- setNames(paste0(res_tab$stratum, "\nN = ", res_tab$N_genes, " genes"), as.character(res_tab$stratum))
res_lab_l <- setNames(paste0(res_tab$stratum, "\nN = ", res_tab$N_loops, " loops"), as.character(res_tab$stratum))
res_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 10), axis.text.x = element_text(size = 8))
p_res_lvl <- ggplot(res_tab, aes(x = stratum, y = resid)) +
  geom_col(fill = "#882255", width = 0.7) +
  geom_errorbar(aes(ymin = rlo, ymax = rhi), width = 0.25) +
  scale_x_discrete(labels = res_lab_g) +
  labs(title = "Absolute residual EU-seq signal", x = NULL,
       y = "% of total WT EU-seq signal remaining (95% CI)") + res_theme
p_res_loop <- ggplot(res_tab, aes(x = stratum, y = lfc)) +
  geom_col(fill = "gray45", width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
  geom_hline(yintercept = 0, color = "gray30") +
  scale_x_discrete(labels = res_lab_l) +
  labs(title = paste0("Loop contact change (", paste(WT_SRC, collapse = " + "), " loops)"), x = NULL,
       y = "Loop log2FC, median (bars = 95% CI)") + res_theme
save_plot_svg(cowplot::plot_grid(p_res_lvl, p_res_loop, nrow = 1),
              "tertile", "residual-transcription-vs-loop-loss", w = 10, h = 5)

save_plot_svg(plot_loop_lfc_by_stratum_box(union_ref, STRAT_COL, STRAT_LEVELS, STRAT_LABEL, type = "histogram"),
              "tertile", "log2fc-histogram-box", w = 14, h = 6)
save_plot_svg(plot_loop_lfc_by_stratum_box(union_ref, STRAT_COL, STRAT_LEVELS, STRAT_LABEL, type = "density"),
              "tertile", "log2fc-density-box", w = 14, h = 6)

# WT loops (WT-only + shared), not split by stratum: log2FC of all WT loops next to WT loops linked to an EU-seq gene
lfc_sets <- c("All WT loops", "WT loops at EU-seq genes")
wt_all <- union_loops %>% dplyr::filter(source %in% c("wt_only", "shared"))
wt_eu  <- union_eu    %>% dplyr::filter(source %in% c("wt_only", "shared"), !is.na(eu_status))
lfc_cmp <- dplyr::bind_rows(data.frame(set = lfc_sets[1], log2fc = wt_all$log2fc),
                            data.frame(set = lfc_sets[2], log2fc = wt_eu$log2fc))
cmp_stats <- lfc_cmp %>% dplyr::group_by(set) %>%
  dplyr::summarise(n = dplyr::n(), med = median(log2fc, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(set = factor(set, levels = lfc_sets)) %>% dplyr::arrange(set) %>%
  dplyr::mutate(lab = factor(sprintf("%s\nN = %d loops, median = %.2f", set, n, med),
                             levels = sprintf("%s\nN = %d loops, median = %.2f", set, n, med)))
message("WT loops log2FC (all vs at EU-seq genes):"); print(as.data.frame(cmp_stats %>% dplyr::select(-lab)), row.names = FALSE)
lfc_cmp <- dplyr::left_join(lfc_cmp %>% dplyr::mutate(set = factor(set, levels = lfc_sets)),
                            cmp_stats %>% dplyr::select(set, lab), by = "set")
p_lfc_cmp <- ggplot(lfc_cmp, aes(x = log2fc)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "gray60", color = NA, alpha = 0.8) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(data = cmp_stats, aes(xintercept = med), color = "#882255", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ lab, nrow = 1) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(title = "Log2FC of WT loops",
       subtitle = paste0(DATA_LABEL, "  |  WT-only + shared loops  |  Log2((POLR2A deg + ", PSEUDOCOUNT, ")/(WT + ", PSEUDOCOUNT,
                         "))  |  solid line = 0, dashed line = median"),
       x = "log2FC (POLR2A deg / WT)", y = "Density") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8, color = "gray40"),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9))
save_plot_svg(p_lfc_cmp, "unstratified", "wt-loops_log2fc-all-vs-eu-genes", w = 10, h = 5)

# WT-only loops, log2FC density of every stratum overlaid in one panel (shared and Pol2-only loops left out)
wt_ov <- union_ref %>% dplyr::filter(source == "wt_only", !is.na(log2fc), !is.na(eu_tertile)) %>%
  dplyr::mutate(stratum = factor(eu_tertile, levels = STRAT_LEVELS))
n_ov <- table(wt_ov$stratum)
ov_labels <- setNames(paste0(names(n_ov), "  (N = ", as.integer(n_ov), ")"), names(n_ov))
ov_colors <- setNames(c("#7f7f7f", "#882255", "#DDAA33", "#004488"), STRAT_LEVELS)
p_wt_overlay <- ggplot(wt_ov, aes(x = log2fc, color = stratum, fill = stratum)) +
  geom_density(alpha = 0.12, linewidth = 0.9) +
  geom_vline(xintercept = -LOG2FC_THRESHOLD, linetype = "dashed", color = "#4393c3", linewidth = 0.6) +
  geom_vline(xintercept =  LOG2FC_THRESHOLD, linetype = "dashed", color = "#882255", linewidth = 0.6) +
  scale_color_manual(values = ov_colors, labels = ov_labels, name = STRAT_LABEL) +
  scale_fill_manual(values = ov_colors, labels = ov_labels, name = STRAT_LABEL) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(title = "Log2FC of WT-only loops by EU-seq stratum",
       subtitle = paste0(DATA_LABEL, "  |  Log2((POLR2A deg + ", PSEUDOCOUNT, ")/(WT + ", PSEUDOCOUNT, "))  |  Threshold = ", LOG2FC_THRESHOLD),
       x = "log2FC (POLR2A deg / WT)", y = "Density") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8, color = "gray40"))
save_plot_svg(p_wt_overlay, "tertile", "wt-only_log2fc-density-overlay", w = 8, h = 5.5)

ridge_plots <- lapply(names(RIDGE_SETS), function(src)
  plot_lost_gained_ridge(
    df = union_ref, source_filter = src, strat_col = STRAT_COL, strat_levels = STRAT_LEVELS,
    title = paste0(RIDGE_SETS[[src]], " — Lost vs Gained loops by ", STRAT_LABEL),
    subtitle = paste0(DATA_LABEL, "  |  source=", src, "  |  |log2FC| > ", LOG2FC_THRESHOLD,
                      " and >= ", MIN_WT_CONTACTS, " contacts  |  SCALE contacts extracted via strawr")))
ridge_plots <- Filter(Negate(is.null), ridge_plots)
if (length(ridge_plots) > 0)
  save_plot_svg(cowplot::plot_grid(plotlist = ridge_plots, ncol = 1),
                "tertile", "lost-gained-ridges", w = 14, h = 5.5 * length(ridge_plots))

# ══════════════════════════════════════════════════════════════════════════════
# 4. APA of loops with an anchor on a gene body, per stratum, for WT-only / shared / Pol2-only
#    loops (gene body +/- GENEBODY_TOL_BP, no cCRE / zone filter)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Gene-body APA ───────────────────────────────────────────────────────")
union_gene <- union_loops %>%
  dplyr::mutate(group = GENEBODY_GROUP) %>%
  link_loops_to_eu_genes(gene_gr, tol_bp = GENEBODY_TOL_BP)
message("Loops with an anchor within ", GENEBODY_TOL_BP / 1000, " kb of an EU-seq gene body: ", sum(!is.na(union_gene$eu_status)), " / ", nrow(union_gene))
print(table(union_gene$source, union_gene$eu_tertile, useNA = "ifany"))
stopifnot(sum(!is.na(union_gene$eu_status)) > 0)
for (src in names(RIDGE_SETS)) {
  wg <- add_ref(union_gene %>% dplyr::filter(source == src))
  message("APA: ", RIDGE_SETS[[src]], " by ", STRAT_LABEL)
  apa <- run_apa_strata(wg, GENEBODY_GROUP, STRAT_COL, STRAT_LEVELS)
  save_plot_svg(
    plot_apa_strata_grid(apa, GENEBODY_GROUP, STRAT_LEVELS, bin_color = "#6a6a6a",
                         title = paste0("APA — ", RIDGE_SETS[[src]], " with an anchor on a gene body by ",
                                        STRAT_LABEL, ": WT vs POLR2A deg scoring")),
    "tertile", paste0("apa-genebody_", gsub("_", "-", src)), w = 14, h = 9)
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. WT loops: APA by stratum (Ref + tertiles) and scoring
# ══════════════════════════════════════════════════════════════════════════════
message("\n── WT loops: APA and scoring ───────────────────────────────────────────")
# short / long cut = median loop distance of the WT loops; expression halves = median WT TPM of the EU-seq genes
dist_cut <- median(union_eu$loop_dist[union_eu$source %in% WT_SRC])
gene_tpm <- ifelse(is.na(gene_gr$TPM), 0, gene_gr$TPM)
expr_cut <- median(gene_tpm)
if (!(expr_cut > 0)) stop("Median WT TPM of the EU-seq genes is ", expr_cut, "; the two expression halves would be unbalanced")
to_expr_half <- function(tpm) factor(ifelse(is.na(tpm) | tpm <= expr_cut, EXPR_LEVELS[1], EXPR_LEVELS[2]), levels = EXPR_LEVELS)
message("Short / long cut: ", round(dist_cut / 1000), " kb   |   expression cut (median WT TPM, EU-seq genes): ", round(expr_cut, 3))

wt_df <- union_eu %>%
  dplyr::filter(source %in% WT_SRC, !is.na(eu_tertile)) %>%
  dplyr::mutate(eu_tertile = as.character(eu_tertile),
                group      = "WT loops",
                dist_class = factor(ifelse(loop_dist <= dist_cut, "short", "long"), levels = names(DIST_COLORS)),
                expr_half  = to_expr_half(TPM))
message("WT loops with an EU-seq tertile: ", nrow(wt_df)); print(table(wt_df$eu_tertile))
n_ref <- sum(!is.na(wt_df$eu_lfc) & abs(wt_df$eu_lfc) <= REF_LFC)
message("WT loops in the reference stratum (|LFC| <= ", REF_LFC, "): ", n_ref)
if (n_ref < 10) warning("Fewer than 10 WT loops sit on genes with |LFC| <= ", REF_LFC, "; the reference column will be empty (widen REF_LFC)")

apa_main <- run_apa_strata(add_ref(wt_df), "WT loops", STRAT_COL, STRAT_LEVELS)
save_plot_svg(
  plot_apa_strata_grid(apa_main, "WT loops", STRAT_LEVELS, bin_color = "#6a6a6a",
                       title = paste0("APA — WT loops by ", STRAT_LABEL, " (reference |EU-seq LFC| <= ", REF_LFC,
                                      "): WT vs POLR2A deg scoring")),
  "tertile", "wt-apa", w = 14, h = 9)

# per-loop matrices of the loops that survive the APA diagonal filter (same loops as the APA panels)
wt_p <- wt_df[abs(floor(wt_df$centroid1 / RESOLUTION) * RESOLUTION -
                    floor(wt_df$centroid2 / RESOLUTION) * RESOLUTION) >= (2 * BUFFER + 1) * RESOLUTION, ]
message("Loops for scoring (after the diagonal filter): ", nrow(wt_p))
pl <- get_loop_matrices(wt_p)
L  <- scoring_inputs(pl)

set.seed(1)
sc_main <- run_scoring(L, make_groups(wt_p))
print_scoring(sc_main, "WT loops by stratum:")
# self-check: the scoring reproduces the numbers printed on the APA panels
est_of <- function(sc, lv) vapply(names(SCORE_METRICS), function(m) sc$long$estimate[sc$long$group == lv & sc$long$metric == m], numeric(1))
chk <- vapply(TERT_LEVELS, function(lv)
  max(abs(unname(apa_delta(apa_main$wt[[lv]], apa_main$pol2[[lv]])) - unname(est_of(sc_main, lv)))), numeric(1))
stopifnot(all(chk < 0.05))
save_plot_svg(
  plot_scoring(sc_main$long, "Loop contact change by EU-seq stratum, by scoring method (WT loops)",
               paste0("log2(POLR2A deg / WT); points = estimate, bars = 95% bootstrap CI over loops (", N_BOOT,
                      " resamples); T1 = most-decreased tertile; Ref = genes with |EU-seq LFC| <= ", REF_LFC)),
  "tertile", "wt-apa_scoring", w = 14, h = 4.8)

# ══════════════════════════════════════════════════════════════════════════════
# 6. WT loops split into short and long loops (median loop distance)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Short vs long WT loops ──────────────────────────────────────────────")
sl_df <- add_ref(union_eu %>%
                   dplyr::filter(source %in% WT_SRC) %>%
                   dplyr::mutate(source = ifelse(loop_dist <= dist_cut, "short", "long")))   # the ridge function filters on `source`
message("Loops by distance class x stratum:"); print(table(sl_df$source, factor(sl_df$eu_tertile, levels = STRAT_LEVELS)))
stopifnot(all(sl_df$loop_dist[sl_df$source == "short"] <= dist_cut), all(sl_df$loop_dist[sl_df$source == "long"] > dist_cut))

save_plot_svg(
  plot_loop_lfc_by_stratum_box(sl_df, STRAT_COL, STRAT_LEVELS, STRAT_LABEL, type = "density",
                               group_levels = names(DIST_COLORS), group_colors = DIST_COLORS),
  "tertile", "wt-short-long_log2fc-density-box", w = 14, h = 6)

sl_ridges <- lapply(names(DIST_COLORS), function(cl)
  plot_lost_gained_ridge(
    df = sl_df, source_filter = cl, strat_col = STRAT_COL, strat_levels = STRAT_LEVELS,
    title = paste0("WT ", cl, " loops — Lost vs Gained loops by ", STRAT_LABEL),
    subtitle = paste0(DATA_LABEL, "  |  ", cl, " = ", ifelse(cl == "short", "<=", ">"), " ", round(dist_cut / 1000),
                      " kb (median WT loop distance)  |  |log2FC| > ", LOG2FC_THRESHOLD, " and >= ",
                      MIN_WT_CONTACTS, " contacts")))
sl_ridges <- Filter(Negate(is.null), sl_ridges)
if (length(sl_ridges) > 0)
  save_plot_svg(cowplot::plot_grid(plotlist = sl_ridges, ncol = 1),
                "tertile", "wt-short-long_lost-gained-ridges", w = 14, h = 5.5 * length(sl_ridges))

for (cl in names(DIST_COLORS)) {
  message("APA: WT ", cl, " loops by ", STRAT_LABEL)
  apa_sl <- run_apa_strata(add_ref(wt_df %>% dplyr::mutate(group = as.character(dist_class))), cl, STRAT_COL, STRAT_LEVELS)
  save_plot_svg(
    plot_apa_strata_grid(apa_sl, cl, STRAT_LEVELS, bin_color = DIST_COLORS[[cl]],
                         title = paste0("APA — WT ", cl, " loops (", ifelse(cl == "short", "<=", ">"), " ",
                                        round(dist_cut / 1000), " kb) by ", STRAT_LABEL, ": WT vs POLR2A deg scoring")),
    "tertile", paste0("wt-apa_", cl), w = 14, h = 9)
}

sc_sl <- lapply(setNames(names(DIST_COLORS), names(DIST_COLORS)), function(cl) {
  r <- run_scoring(L, make_groups(wt_p, wt_p$dist_class == cl))
  print_scoring(r, paste0("WT ", cl, " loops by stratum:"))
  r$long$split <- cl
  r
})
save_plot_svg(
  plot_scoring(dplyr::bind_rows(lapply(sc_sl, `[[`, "long")),
               "Loop contact change by EU-seq stratum, short vs long WT loops",
               paste0("log2(POLR2A deg / WT); short <= ", round(dist_cut / 1000), " kb < long (median WT loop distance); ",
                      "bars = 95% bootstrap CI over loops"),
               split_colors = DIST_COLORS),
  "tertile", "wt-short-long_scoring", w = 14, h = 5.2)

# ══════════════════════════════════════════════════════════════════════════════
# 7. WT loops split into two halves by steady-state expression
#    (WT TPM of the loop's linked gene vs the median WT TPM of the EU-seq genes)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Expression halves ───────────────────────────────────────────────────")
message("EU-seq genes, EU-seq tertile x expression half:")
print(table(EU = gene_gr$eu_tertile, Expr = to_expr_half(gene_tpm)))
message("WT loops, EU-seq tertile x expression half:")
print(table(EU = wt_df$eu_tertile, Expr = wt_df$expr_half))

for (k in seq_along(EXPR_LEVELS)) {
  ex <- EXPR_LEVELS[k]
  sub <- wt_df %>% dplyr::filter(expr_half == ex)
  message("APA: ", ex, " (", nrow(sub), " loops)")
  apa_ex <- run_apa_strata(add_ref(sub), "WT loops", STRAT_COL, STRAT_LEVELS)
  save_plot_svg(
    plot_apa_strata_grid(apa_ex, "WT loops", STRAT_LEVELS, bin_color = unname(EXPR_COLORS[k]),
                         title = paste0("APA — WT loops at ", ex, " genes by ", STRAT_LABEL, ": WT vs POLR2A deg scoring")),
    "tertile", paste0("wt-apa_expr-", EXPR_TAGS[k]), w = 14, h = 9)
}

sc_ex <- lapply(setNames(EXPR_LEVELS, EXPR_LEVELS), function(ex) {
  r <- run_scoring(L, make_groups(wt_p, wt_p$expr_half == ex))
  print_scoring(r, paste0("WT loops, ", ex, ", by stratum:"))
  r$long$split <- ex
  r
})
save_plot_svg(
  plot_scoring(dplyr::bind_rows(lapply(sc_ex, `[[`, "long")),
               "Loop contact change by EU-seq stratum within expression halves (WT loops)",
               paste0("log2(POLR2A deg / WT); expression halves split at the median WT TPM of the EU-seq genes (",
                      round(expr_cut, 2), " TPM); bars = 95% bootstrap CI over loops"),
               split_colors = EXPR_COLORS),
  "tertile", "wt-expr-halves_scoring", w = 14, h = 5.2)

# ── Scatter: WT loops, EU-seq LFC (x) vs loop LFC (y), raw log2 values, one panel per WT expression tertile ──
# Expression tertiles = WT TPM tertiles over the EU-seq genes; a loop takes the TPM of its linked gene.
EXPR3_LEVELS <- c("Expr T1 (lowest)", "Expr T2 (middle)", "Expr T3 (highest)")
expr3_cuts   <- quantile(gene_tpm, c(1/3, 2/3))
stopifnot(expr3_cuts[[1]] < expr3_cuts[[2]])
message("WT TPM tertile cuts over ", length(gene_gr), " EU-seq genes: ", paste(round(expr3_cuts, 3), collapse = " | "))
to_expr3 <- function(tpm) factor(dplyr::case_when(
  is.na(tpm) | tpm <= expr3_cuts[[1]] ~ EXPR3_LEVELS[1],
  tpm <= expr3_cuts[[2]]              ~ EXPR3_LEVELS[2],
  TRUE                                ~ EXPR3_LEVELS[3]), levels = EXPR3_LEVELS)
sc_df <- wt_df %>%
  dplyr::filter(!is.na(eu_abs), !is.na(log2fc)) %>%
  dplyr::mutate(expr3 = to_expr3(TPM))
message("Genes (EU-seq tertile x WT expression tertile):")
print(table(EU = gene_gr$eu_tertile, Expr = to_expr3(gene_tpm)))
message("WT loops plotted (EU-seq tertile x WT expression tertile):")
print(table(EU = sc_df$eu_tertile, Expr = sc_df$expr3))
sc_stats <- sc_df %>%
  dplyr::group_by(expr3) %>%
  dplyr::summarise(n   = dplyr::n(),
                   rho = unname(suppressWarnings(cor.test(eu_abs, log2fc, method = "spearman", exact = FALSE)$estimate)),
                   p   = suppressWarnings(cor.test(eu_abs, log2fc, method = "spearman", exact = FALSE)$p.value),
                   .groups = "drop") %>%
  dplyr::mutate(lab = sprintf("%s\nN = %d loops\nSpearman rho = %.3f, p = %s", expr3, n, rho, signif(p, 2)))
message("EU-seq absolute change vs loop LFC, WT loops by expression tertile:"); print(as.data.frame(sc_stats %>% dplyr::select(-lab)), row.names = FALSE)
sc_df <- dplyr::left_join(sc_df, sc_stats %>% dplyr::select(expr3, lab), by = "expr3")
p_scatter <- ggplot(sc_df, aes(x = eu_abs, y = log2fc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 0.6, alpha = 0.3, color = "gray30") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#882255", fill = "#882255", alpha = 0.2) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  facet_wrap(~ lab, nrow = 1) +
  labs(title = "WT loops: EU-seq absolute change vs loop LFC, by WT expression tertile",
       subtitle = paste0(DATA_LABEL, "  |  loop LFC = log2((POLR2A deg + ", PSEUDOCOUNT, ") / (WT + ", PSEUDOCOUNT, ")) of the loop pixel  |  ",
                         "expression tertile cuts (WT TPM) ", paste(round(expr3_cuts, 2), collapse = " | ")),
       x = "EU-seq absolute change of the linked gene (scaled RPKM, POLR2A deg - WT; pseudo-log axis)", y = "Loop LFC, log2(POLR2A deg / WT)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8, color = "gray40"),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9))
save_plot_svg(p_scatter, "tertile", "wt-scatter_euseq-lfc-vs-loop-lfc_expr-tertiles", w = 14, h = 5.5)

# ══════════════════════════════════════════════════════════════════════════════
# 8. Anchor classification (narrowPeak + ATAC + TSS geography; repo classifier)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Classifying loop anchors ────────────────────────────────────────────")
for (pkg in c("ChIPseeker", "GenomicFeatures")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package not installed: ", pkg)
}
ep <- new.env()                    # keeps the repo's constants out of this session
sys.source(path.expand(EP_FUNCS_FILE), envir = ep)
message("Building TxDb from ", GTF_FILE)
np_txdb  <- GenomicFeatures::makeTxDbFromGFF(GTF_FILE, format = "gtf")
np_peaks <- lapply(NP_FILES, ep$load_narrowpeak)

# loop coordinates written with the chromosome style of the TxDb (ChIPseeker needs a match)
tx_chr  <- any(grepl("^chr", GenomeInfoDb::seqlevels(np_txdb)))
chr_fmt <- function(x) if (tx_chr) paste0("chr", sub("^chr", "", x)) else sub("^chr", "", x)
classify_loops <- function(loops_df, label) {
  np_in <- loops_df
  np_in$chr1 <- chr_fmt(loops_df$chr1)
  np_in$chr2 <- chr_fmt(loops_df$chr2)
  res <- ep$classify_anchors_narrowpeak(
    loops_df = np_in, atac_gr = np_peaks$atac, h3k4me3_gr = np_peaks$h3k4me3,
    h3k27ac_gr = np_peaks$h3k27ac, h3k4me1_gr = np_peaks$h3k4me1,
    txdb = np_txdb, tss_region = NP_TSS_REGION, ccre_label = label)
  stopifnot(nrow(res) == nrow(loops_df))
  res$chr1 <- loops_df$chr1          # restore the .hic / HiCCUPS chromosome naming for APA
  res$chr2 <- loops_df$chr2
  res
}
wt_np_class   <- classify_loops(loops_wt,   "WT_narrowpeak")
pol2_np_class <- classify_loops(loops_pol2, "POLR2A_deg_narrowpeak")
# self-checks: row-preserving, expected columns, and a chromosome-style mismatch would leave
# no promoter anchors
stopifnot(nrow(wt_np_class) == nrow(loops_wt), nrow(pol2_np_class) == nrow(loops_pol2),
          all(c("anchor1_ccre", "anchor2_ccre", "loop_ccre_type", "anchor1_zone", "anchor2_zone") %in%
                colnames(wt_np_class)),
          sum(c(wt_np_class$anchor1_zone, wt_np_class$anchor2_zone) == "PROMOTER") > 0)
message("Wildtype pair types (top 12):")
print(head(sort(table(wt_np_class$loop_ccre_type), decreasing = TRUE), 12))

# union loops carry the anchor classes of the loop set they came from (shared / WT-only loops
# from the Wildtype classification, Pol2-only loops from the POLR2A-degraded classification)
union_cls <- dplyr::bind_rows(wt_np_class[loop_overlap$shared$a_idx, ],
                              wt_np_class[loop_overlap$a_only, ],
                              pol2_np_class[loop_overlap$b_only, ])
stopifnot(nrow(union_cls) == nrow(union_eu),
          identical(union_cls$start1, union_eu$start1), identical(union_cls$end2, union_eu$end2),
          identical(as.character(union_cls$chr1), as.character(union_eu$chr1)))
union_eu$pe_loop <- union_cls$anchor1_ccre %in% ep_anchor_types | union_cls$anchor2_ccre %in% ep_anchor_types
message("Loops with a promoter or enhancer anchor: ", sum(union_eu$pe_loop), " / ", nrow(union_eu))

# ══════════════════════════════════════════════════════════════════════════════
# 9. Non-gene anchor type at genes and Gained vs Lost counts per stratum (Ref + tertiles)
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Non-gene anchor types ───────────────────────────────────────────────")
all_singleton <- bind_rows(link_loops_singleton_eu(wt_np_class,   gene_gr, "Wildtype"),
                           link_loops_singleton_eu(pol2_np_class, gene_gr, "POLR2A deg")) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Wildtype", "POLR2A deg")))
all_singleton <- add_ref(all_singleton)
e_plots <- lapply(STRAT_LEVELS, function(lv)
  make_singleton_count_plot_eu(all_singleton, STRAT_COL, lv, ep_anchor_types,
                               title = paste0("Non-gene anchor type at EU-seq ", lv, " genes"),
                               strat_label = STRAT_LABEL))
save_plot_svg(cowplot::plot_grid(plotlist = e_plots, ncol = 1),
              "tertile", "nongene-anchor-type", w = 8, h = 4.5 * length(e_plots))

union_ref <- add_ref(union_eu)      # again, now with the pe_loop column
gl_sub <- paste0(DATA_LABEL, "  |  Lost / Gained: |log2FC| > ", LOG2FC_THRESHOLD, " and >= ",
                 MIN_WT_CONTACTS, " contacts  |  loops by ", STRAT_LABEL)
save_plot_svg(make_gain_loss_bar(union_ref %>% filter(pe_loop), STRAT_COL, STRAT_LEVELS,
                                 "Lost vs Gained loops with a promoter or enhancer anchor", gl_sub),
              "tertile", "gained-lost-bar_promoter-enhancer", w = 8, h = 4.5)
save_plot_svg(make_gain_loss_bar(union_ref, STRAT_COL, STRAT_LEVELS,
                                 "Lost vs Gained loops at all genes", gl_sub),
              "tertile", "gained-lost-bar_all-genes", w = 8, h = 4.5)

message("\nDone.")