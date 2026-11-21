##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Shared main script for the loop loss vs. Pol II occupancy vs. histone marks
# analyses. Holds the parameters, the functions and the feature build that the
# regression scripts use; it runs no regression itself.
#
# Loops are the union of the Wildtype and POLR2A-degraded HiCCUPS calls
# (shared and Wildtype-only loops at Wildtype positions, POLR2A-degraded-only
# loops at their own positions). Per loop it builds: Wildtype and
# POLR2A-degraded contacts, log2FC, Lost/Gained class, anchor signal from
# ENCODE POLR2A ChIP-seq (ENCSR000EUU, untreated HCT116) and the Table 5
# Mint-ChIP-seq (POLR2A-degron HCT116, DMSO vs +5-Ph-IAA), and the Wildtype
# TPM bin of the highest-expressed gene at either anchor.
#
# Assay switch: ASSAY = "hic" (default; HiCCUPS loops and mega Hi-C maps) or
# "microc" (ENCODE POLR2A-degron HCT116 Micro-C, ENCSR533VUS / ENCSR218URU).
# Set it before running a regression script (ASSAY <- "microc") or export
# LOOP_ASSAY=microc. Only the input paths and the contact-extraction mode differ.
#
# Per loop it also records the linked anchor gene (anchor_gene_id) and, when the
# ENCODE cCRE file is available, the promoter / enhancer type of each anchor
# (loop_type: P-P, E-P, E-E, Other).
#
# Objects left in the session for the regression scripts:
#   features   all union loops, one row per loop
#   wt_reg     Wildtype loops with Wildtype contacts >= MIN_WT_CONTACTS
#   run_lm_models(), run_logit_models(), compare_lm_formulas(), plot_model_results(),
#   MARK_TERMS, NUM_TERMS, ...
#
# Sourced by:
#   fork_pooled_curves.R     pooled loop change across transcription and histone-mark change
#
# Adapted from: Hic_loop_analysis/hic_loop_lfc_analysis.R
#   (loop union, Lost/Gained rule, find_shared_loops) and
#   Hic_loop_analysis/final_APA_analysis.R (gene linking, expression bins)
##################################################################
library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(strawr)
library(ggplot2)

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════
# ── Assay switch ──────────────────────────────────────────────────────────────
# "hic" (default) or "microc". Set ASSAY before running a script, or export
# LOOP_ASSAY. Everything outside the two branches below is shared.
if (!exists("ASSAY")) ASSAY <- Sys.getenv("LOOP_ASSAY", unset = "hic")
stopifnot(ASSAY %in% c("hic", "microc"))

PROJ     <- "/projects/b1042/BackmanLab/Lucas/loop_chromatin_regression"
MC_DIR   <- "/projects/b1042/BackmanLab/Lucas/microc_encode"
N_CORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2"))

if (ASSAY == "hic") {
  BASE     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
  LOOP_DIR <- file.path(BASE, "loops")
  
  DS_GROUP  <- "ds_wt_pol2"   # loop calls:           "ds_wt_pol2" | "ds_all"
  HIC_GROUP <- "ds_wt_pol2"   # contact maps (.hic):  "ds_wt_pol2" | "ds_all"
  stopifnot(DS_GROUP  %in% c("ds_wt_pol2", "ds_all"),
            HIC_GROUP %in% c("ds_wt_pol2", "ds_all"))
  if (DS_GROUP != HIC_GROUP) {
    message("CHECK: loop calls (", DS_GROUP, ") and contact maps (", HIC_GROUP,
            ") come from different sets")
  }
  
  hic_wt_loops   <- file.path(LOOP_DIR, DS_GROUP, "WT_HCT116_CTRL/mega/merged_loops.bedpe")
  hic_pol2_loops <- file.path(LOOP_DIR, DS_GROUP, "Pol2_6hrs_Aux/mega/merged_loops.bedpe")
  hic_wt_file    <- file.path(BASE, paste0("mega_", HIC_GROUP), "WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  hic_pol2_file  <- file.path(BASE, paste0("mega_", HIC_GROUP), "Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  HIC_BYTES      <- NULL                 # no size check for the lab's own maps
  CONTACT_MODE   <- "block"              # one strawr query per chromosome
  OUT_DIR        <- file.path(PROJ, "out")
} else {
  # ENCODE POLR2A-degron HCT116 Micro-C: control ENCSR533VUS, 6 h 5-Ph-IAA ENCSR218URU.
  # ASSUMED pairing of loop calls and matrices (not yet checked against the
  # ENCODE derived_from field): ENCFF954EMJ <-> ENCFF534TYD, ENCFF604OYY <-> ENCFF302CUI.
  DS_GROUP  <- "ENCODE Micro-C loops"
  HIC_GROUP <- "ENCODE Micro-C maps"
  hic_wt_loops   <- file.path(MC_DIR, "ENCFF954EMJ.bedpe.gz")
  hic_pol2_loops <- file.path(MC_DIR, "ENCFF604OYY.bedpe.gz")
  hic_wt_file    <- file.path(MC_DIR, "ENCFF534TYD.hic")
  hic_pol2_file  <- file.path(MC_DIR, "ENCFF302CUI.hic")
  HIC_BYTES      <- c(wt = 27301466053, pol2 = 33785386023)   # ENCODE file sizes
  CONTACT_MODE   <- "per_loop"           # one strawr query per loop pixel (cached)
  OUT_DIR        <- file.path(PROJ, "out", "microc")
}
message("Assay: ", ASSAY, "  |  loop calls: ", DS_GROUP, "  |  contact maps: ", HIC_GROUP)
MIN_OBS_WT    <- 0
MIN_OBS_POL2  <- 0

# ── Promoter / enhancer annotation of anchors (ENCODE Registry-V4 cCREs, HCT116 Core) ──
# Required for ASSAY = "microc"; for "hic" it is used only if the file exists.
CCRE_FILE <- file.path(MC_DIR, "anchor_functional_analysis",
                       "ENCFF431JDU_ENCFF964OOU_ENCFF787LMI_ENCFF388PVO.bed")
CCRE_REQUIRED <- (ASSAY == "microc")
CCRE_SLOP_BP  <- 2500L       # anchor +/- this distance may overlap a cCRE (EP_analysis_functions.R)

# ── strawr / LFC parameters (same as hic_loop_lfc_analysis.R) ─────────────────
STRAW_NORM        <- "SCALE"
STRAW_RESOLUTION  <- 5000L
MIN_WT_CONTACTS   <- 5       # minimum Wildtype contacts to call "Lost"
MIN_POL2_CONTACTS <- 5       # minimum POLR2A-degraded contacts to call "Gained"
LOG2FC_THRESHOLD  <- 0.5
PSEUDOCOUNT       <- 1
ANCHOR_SLOP_BP    <- 5000L   # slop for shared-loop matching

# ── Expression (Wildtype TPM) ─────────────────────────────────────────────────
REF_DIR     <- "/projects/b1042/BackmanLab/Lucas/101326_hg38_rnaseq/reference"
GTF_FILE    <- file.path(REF_DIR, "Homo_sapiens.GRCh38.115.chr.gtf")
EXPR_FILE   <- file.path(REF_DIR, "ENCFF435PHM.tsv")
GENE_TOL_BP <- 5000L         # gene body +/- this distance may overlap an anchor

# ── Signal tracks (one bigWig per mark per condition) ─────────────────────────
DATA_DIR <- file.path(PROJ, "data")
MARKS    <- c("H3K27ac", "H3K4me3", "H3K9me3")
BW_FILES <- c(
  Pol2 = file.path(DATA_DIR, "POLR2A_ENCODE.bw"),
  setNames(file.path(DATA_DIR, paste0(MARKS, "_DMSO.bw")), paste0(MARKS, "_DMSO")),
  setNames(file.path(DATA_DIR, paste0(MARKS, "_IAA.bw")),  paste0(MARKS, "_IAA"))
)

# ── Cross-validation ──────────────────────────────────────────────────────────
CV_FOLDS <- 10L
CV_SEED  <- 1L

# ── SVG export of plots: set PLOT_DIR in the regression script before it sources this file ──
# (NULL = do not save). Every saved plot has the same size (inches).
if (!exists("PLOT_DIR")) PLOT_DIR <- NULL
PLOT_W <- 8   # standard SVG size (inches)
PLOT_H <- 6
if (!is.null(PLOT_DIR) && !dir.exists(PLOT_DIR)) stop("PLOT_DIR not found: ", PLOT_DIR)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Functions
# ══════════════════════════════════════════════════════════════════════════════
check_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "))
  }
  invisible(TRUE)
}

# ── From hic_loop_lfc_analysis.R (load_hiccups_bedpe adds a column-name check) ─
load_hiccups_bedpe <- function(path) {
  all_lines   <- readLines(path)
  header_line <- all_lines[grepl("^#chr1", all_lines)][1]
  col_names   <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  data_lines  <- all_lines[!grepl("^#", all_lines)]
  if (length(data_lines) == 0) stop("No data lines in ", path)
  df <- read.table(text=paste(data_lines, collapse="\n"),
                   header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
  if (ncol(df) == length(col_names)) colnames(df) <- col_names
  need <- c("observed", "expectedDonut", "fdrDonut", "fdrBL", "fdrH", "fdrV",
            "centroid1", "centroid2")
  if (!all(need %in% colnames(df))) {
    stop(basename(path), ": expected bedpe columns not found (", ncol(df), " data columns vs ",
         length(col_names), " header names); missing: ",
         paste(setdiff(need, colnames(df)), collapse = ", "))
  }
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

check_loops <- function(df, label) {
  stopifnot(nrow(df) > 0,
            !anyNA(df[, c("chr1", "start1", "end1", "chr2", "start2", "end2",
                          "centroid1", "centroid2")]),
            all(df$chr1 == df$chr2),
            all(df$loop_dist > 0))
  message(label, " loops: ", nrow(df))
}

extract_loop_contacts <- function(loops_df, hic_file,
                                  norm      = STRAW_NORM,
                                  resolution = STRAW_RESOLUTION,
                                  matrix    = "observed") {
  snap    <- function(x) as.integer(floor(x / resolution) * resolution)
  results <- rep(NA_real_, nrow(loops_df))
  
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
  
  # Chromosome names exactly as stored in this .hic ("chr1" or "1")
  hic_chroms <- strawr::readHicChroms(hic_file)$name
  hic_name   <- setNames(hic_chroms, sub("^chr", "", hic_chroms))
  absent     <- setdiff(unique(c(chr1_clean, chr2_clean)), names(hic_name))
  if (length(absent) > 0) {
    stop(basename(hic_file), ": loop chromosome(s) not in .hic: ",
         paste(absent, collapse = ", "))
  }
  message("Extracting ", matrix, " from ", basename(hic_file),
          "  norm=", norm, "  res=", resolution/1000, "kb  n=", nrow(loops_df))
  for (cp in chr_pairs) {
    parts  <- strsplit(cp, "\\|")[[1]]
    c1_chr <- parts[1]; c2_chr <- parts[2]
    idx    <- which(chr1_clean == c1_chr & chr2_clean == c2_chr)
    c1_min <- min(c1_pos[idx]); c1_max <- max(c1_pos[idx]) + resolution
    c2_min <- min(c2_pos[idx]); c2_max <- max(c2_pos[idx]) + resolution
    loc1   <- paste0(hic_name[[c1_chr]], ":", c1_min, ":", c1_max)
    loc2   <- paste0(hic_name[[c2_chr]], ":", c2_min, ":", c2_max)
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

# ── Micro-C: contact at the pixel holding both loop centroids ─────────────────
# One strawr call per small 2D tile of loop pixels (not one per pixel): every strawr
# call re-reads the .hic index and normalisation vectors, which dominates the time
# on a 27-34 GB Micro-C file. Pixels are grouped by chromosome and by the
# `tile_bp` block holding each anchor, and each tile is read as one rectangle.
# Pixel with no record -> 0; a failed query -> NA (counted and reported).
# Results are cached per .hic file in `cache_file` (keyed by chromosome:start:end of
# the pixel, saved after every chunk of tiles), so an interrupted run resumes.
# The first chunk is checked against single-pixel queries (the previous method).
extract_contacts_per_loop <- function(loops_df, hic_file, cache_file,
                                      norm = STRAW_NORM, resolution = STRAW_RESOLUTION,
                                      n_cores = N_CORES, tile_bp = 2000000L, chunk = 50L) {
  hic_chroms <- strawr::readHicChroms(hic_file)$name
  hic_name   <- setNames(hic_chroms, sub("^chr", "", hic_chroms))
  chr_clean  <- sub("^chr", "", loops_df$chr1)
  absent     <- setdiff(unique(chr_clean), names(hic_name))
  if (length(absent) > 0) {
    stop(basename(hic_file), ": loop chromosome(s) not in .hic: ", paste(absent, collapse = ", "))
  }
  snap <- function(x) as.integer(floor(x / resolution) * resolution)
  lo   <- snap(pmin(loops_df$centroid1, loops_df$centroid2))
  hi   <- snap(pmax(loops_df$centroid1, loops_df$centroid2))
  chr  <- unname(hic_name[chr_clean])
  key  <- paste(chr, lo, hi, sep = ":")
  
  meta  <- list(file = basename(hic_file), size = file.size(hic_file),
                norm = norm, resolution = resolution)
  cache <- c(meta, list(vals = numeric(0)))
  if (file.exists(cache_file)) {
    old <- readRDS(cache_file)
    if (identical(old[names(meta)], meta)) {
      cache <- old
    } else {
      message("  Cache ", basename(cache_file), " is for a different file/norm/resolution; ignored")
    }
  }
  need <- setdiff(unique(key), names(cache$vals))
  message("Extracting contacts from ", basename(hic_file), "  norm=", norm, "  res=",
          resolution / 1000, "kb  pixels=", length(unique(key)), "  cached=",
          length(unique(key)) - length(need), "  cores=", n_cores)
  if (length(need) > 0) {
    first <- match(need, key)
    tile  <- paste(chr[first], lo[first] %/% tile_bp, hi[first] %/% tile_bp, sep = ":")
    tiles <- unname(split(seq_along(need), tile))
    message("  ", length(need), " new pixels in ", length(tiles), " tiles (",
            tile_bp / 1e6, " Mb x ", tile_bp / 1e6, " Mb)")
    
    # Contact at one pixel by a single-pixel query (check only)
    pixel_query <- function(i) {
      mat <- strawr::straw(norm, hic_file,
                           paste0(chr[i], ":", lo[i], ":", lo[i] + resolution - 1L),
                           paste0(chr[i], ":", hi[i], ":", hi[i] + resolution - 1L), "BP", resolution)
      v <- mat$counts[mat$x == lo[i] & mat$y == hi[i]]
      if (length(v) == 0L) 0 else if (is.finite(v[1])) as.numeric(v[1]) else 0
    }
    # Contacts at all pixels of one tile from one rectangle query
    tile_query <- function(ix) {
      i <- first[ix]
      tryCatch({
        mat <- strawr::straw(norm, hic_file,
                             paste0(chr[i[1]], ":", min(lo[i]), ":", max(lo[i]) + resolution - 1L),
                             paste0(chr[i[1]], ":", min(hi[i]), ":", max(hi[i]) + resolution - 1L),
                             "BP", resolution)
        look <- setNames(mat$counts, paste(pmin(mat$x, mat$y), pmax(mat$x, mat$y), sep = ":"))
        v <- unname(look[paste(pmin(lo[i], hi[i]), pmax(lo[i], hi[i]), sep = ":")])
        v[!is.finite(v)] <- 0
        v
      }, error = function(e) rep(NA_real_, length(ix)))
    }
    
    groups <- split(seq_along(tiles), ceiling(seq_along(tiles) / chunk))
    n_fail <- 0L; n_done <- 0L
    t0 <- Sys.time()
    for (g in seq_along(groups)) {
      tl  <- tiles[groups[[g]]]
      res <- parallel::mclapply(tl, tile_query, mc.cores = n_cores)
      ok  <- vapply(res, function(z) is.numeric(z), logical(1))
      ix  <- unlist(tl[ok]); v <- unlist(res[ok])
      n_fail <- n_fail + sum(lengths(tl[!ok])) + sum(is.na(v))
      good <- !is.na(v)
      if (g == 1L) {
        # Check on the live file: tile values equal single-pixel queries
        pick <- ix[good]
        if (length(pick) > 0) {
          pick <- pick[order(-v[good])][unique(round(seq(1, length(pick), length.out = min(6, length(pick)))))]
        }
        for (k in pick) {
          single <- pixel_query(first[k])
          tilev  <- v[match(k, ix)]
          if (!isTRUE(all.equal(single, tilev))) {
            stop("Tile query and single-pixel query differ at ", need[k], ": ", tilev, " vs ", single)
          }
        }
        message("  Check: tile queries equal single-pixel queries at ", length(pick), " pixels")
      }
      cache$vals <- c(cache$vals, setNames(v[good], need[ix[good]]))
      saveRDS(cache, cache_file)
      n_done <- n_done + length(unlist(tl))
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message("  ", n_done, " / ", length(need), " new pixels queried (", n_fail, " failed)  ",
              round(el, 1), " min elapsed, ~", round(el / n_done * (length(need) - n_done), 1),
              " min left")
    }
    if (n_fail > 0) message("  CHECK: ", n_fail, " pixel queries failed; they are re-queried on the next run")
  }
  out <- unname(cache$vals[key])
  n_hit <- sum(!is.na(out) & out > 0)
  message("  Signal: ", n_hit, " / ", length(out), " (", round(100 * n_hit / length(out), 1), "%)")
  out
}

# One entry point: block mode for the lab's Hi-C maps, per-loop mode for Micro-C
get_contacts <- function(loops_df, hic_file, tag) {
  if (CONTACT_MODE == "per_loop") {
    extract_contacts_per_loop(loops_df, hic_file,
                              file.path(OUT_DIR, paste0("contacts_cache_", tag, ".rds")))
  } else {
    extract_loop_contacts(loops_df, hic_file)
  }
}

# ── Width-weighted mean bigWig signal over each anchor ────────────────────────
# Positions with no bigWig interval count as 0. Returns one value per anchor.
anchor_signal <- function(bw_path, chr, s_pos, e_pos) {
  bwf       <- rtracklayer::BigWigFile(bw_path)
  bw_levels <- GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(bwf))
  chr_clean <- sub("^chr", "", chr)
  chr_bw    <- if (any(grepl("^chr", bw_levels))) paste0("chr", chr_clean) else chr_clean
  absent    <- setdiff(unique(chr_bw), bw_levels)
  if (length(absent) > 0) {
    stop(basename(bw_path), ": loop chromosome(s) not in bigWig: ",
         paste(absent, collapse = ", "))
  }
  anch <- GRanges(chr_bw, IRanges(s_pos + 1L, e_pos))
  win  <- reduce(anch)
  sig  <- rtracklayer::import(bw_path, format = "bigWig", which = win)
  out  <- numeric(length(anch))
  if (length(sig) == 0) return(out)
  score <- sig$score
  score[is.na(score)] <- 0
  hits  <- findOverlaps(anch, sig)
  q     <- queryHits(hits)
  s     <- subjectHits(hits)
  ov_bp <- pmin(end(anch)[q], end(sig)[s]) - pmax(start(anch)[q], start(sig)[s]) + 1L
  num   <- rowsum(ov_bp * score[s], q)
  out[as.integer(rownames(num))] <- num[, 1]
  out / width(anch)
}

# Pseudocount = 1st percentile of the positive values (printed below)
pc_for <- function(v) {
  v <- v[v > 0]
  stopifnot(length(v) > 0)
  unname(quantile(v, 0.01))
}

# ── Highest-TPM gene within GENE_TOL_BP of either anchor (final_APA_analysis.R) ──
link_loops_to_genes <- function(loops_df, gene_gr, tol_bp = GENE_TOL_BP) {
  chr_of     <- function(x) paste0("chr", sub("^chr", "", x))
  anchor1_gr <- GRanges(chr_of(loops_df$chr1), IRanges(loops_df$start1, loops_df$end1))
  anchor2_gr <- GRanges(chr_of(loops_df$chr2), IRanges(loops_df$start2, loops_df$end2))
  gene_gr_exp <- trim(resize(gene_gr, width(gene_gr) + 2*tol_bp, fix="center"))
  hits1 <- findOverlaps(anchor1_gr, gene_gr_exp)
  hits2 <- findOverlaps(anchor2_gr, gene_gr_exp)
  pairs <- dplyr::bind_rows(
    data.frame(loop_idx = queryHits(hits1), gene_idx = subjectHits(hits1)),
    data.frame(loop_idx = queryHits(hits2), gene_idx = subjectHits(hits2))
  ) %>%
    dplyr::distinct(loop_idx, gene_idx) %>%
    dplyr::mutate(anchor_gene_tpm = gene_gr$TPM[gene_idx],
                  expr_bin        = gene_gr$expr_bin[gene_idx],
                  anchor_gene_id  = gene_gr$gene_id[gene_idx])
  ann <- pairs %>%
    dplyr::group_by(loop_idx) %>%
    dplyr::slice_max(anchor_gene_tpm, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(loop_idx, anchor_gene_tpm, expr_bin, anchor_gene_id)
  out <- loops_df
  out$loop_idx <- seq_len(nrow(out))
  out <- dplyr::left_join(out, ann, by = "loop_idx")
  out$expr_bin <- tidyr::replace_na(out$expr_bin, "Unannotated")
  out$loop_idx <- NULL
  out
}

# ── ENCODE cCREs -> promoter / enhancer anchors (priority rule of EP_analysis_functions.R) ──
# Keeps PLS (Promoter), pELS (Proximal-Enhancer), dELS (Distal-Enhancer); CTCF-only
# and other classes are dropped, as in classify_anchors_priority(). The class
# column is found by its values (a ",CTCF-bound" suffix is ignored).
load_ccre_gr <- function(path) {
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", comment.char = "", fill = TRUE)
  message("cCRE file: ", nrow(df), " rows x ", ncol(df), " columns; first row: ",
          paste(unlist(df[1, ]), collapse = " | "))
  if (!is.numeric(df$V2) || !is.numeric(df$V3)) stop("cCRE file: V2/V3 are not numeric (header line?)")
  if (!any(grepl("^chr", df$V1))) stop("cCRE file: V1 holds no chr-prefixed chromosome names")
  types <- c("PLS", "pELS", "dELS", "CA-CTCF")
  class_cols <- names(df)[vapply(df, function(x) is.character(x) &&
                                   any(sub(",.*", "", x) %in% types), logical(1))]
  if (length(class_cols) != 1) {
    stop("cCRE file: expected exactly one cCRE class column, found: ",
         paste(class_cols, collapse = ", "))
  }
  cls <- sub(",.*", "", df[[class_cols]])
  message("cCRE classes in file:")
  print(table(cls))
  type <- c(PLS = "Promoter", pELS = "Proximal-Enhancer", dELS = "Distal-Enhancer")[cls]
  keep <- !is.na(type)
  stopifnot(sum(keep) > 0)
  GRanges(df$V1[keep], IRanges(df$V2[keep] + 1L, df$V3[keep]),   # BED 0-based -> 1-based
          ccre_type = unname(type[keep]))
}

# Per anchor: highest-priority overlapping cCRE (Promoter > Proximal- > Distal-Enhancer),
# else "Unannotated". Per loop: P-P, E-P (promoter + any enhancer), E-E, else Other.
annotate_anchor_types <- function(loops_df, ccre_gr, slop_bp = CCRE_SLOP_BP) {
  prio   <- c("Promoter" = 1L, "Proximal-Enhancer" = 2L, "Distal-Enhancer" = 3L)
  chr_of <- function(x) paste0("chr", sub("^chr", "", x))
  classify_one <- function(chr, s, e) {
    anchor <- GRanges(chr_of(chr), IRanges(pmax(1L, s - slop_bp), e + slop_bp))
    hits   <- findOverlaps(anchor, ccre_gr)
    res    <- rep("Unannotated", length(anchor))
    if (length(hits) > 0) {
      best <- tapply(prio[ccre_gr$ccre_type[subjectHits(hits)]], queryHits(hits), min)
      res[as.integer(names(best))] <- names(prio)[as.integer(best)]
    }
    res
  }
  a1  <- classify_one(loops_df$chr1, loops_df$start1, loops_df$end1)
  a2  <- classify_one(loops_df$chr2, loops_df$start2, loops_df$end2)
  enh <- c("Proximal-Enhancer", "Distal-Enhancer")
  loops_df$anchor1_ccre <- a1
  loops_df$anchor2_ccre <- a2
  loops_df$loop_type <- dplyr::case_when(
    a1 == "Promoter" & a2 == "Promoter"                             ~ "P-P",
    (a1 == "Promoter" & a2 %in% enh) | (a2 == "Promoter" & a1 %in% enh) ~ "E-P",
    a1 %in% enh & a2 %in% enh                                       ~ "E-E",
    TRUE                                                            ~ "Other"
  )
  loops_df$loop_type <- factor(loops_df$loop_type, levels = c("P-P", "E-P", "E-E", "Other"))
  loops_df
}

# Custom nested linear models on log2fc: z-scores `z_terms`, then CV R2 and AIC per formula.
# Returns the fits so the caller can run anova() on any pair.
compare_lm_formulas <- function(reg, formulas, tag, z_terms) {
  message("\n══ ", tag, "  (n = ", nrow(reg), " loops) ══════════════════════════")
  stopifnot(all(c("log2fc", z_terms) %in% colnames(reg)),
            all(is.finite(as.matrix(reg[, c("log2fc", z_terms)]))))
  reg[z_terms] <- lapply(reg[z_terms], function(x) as.numeric(scale(x)))
  fits <- lapply(formulas, lm, data = reg)
  set.seed(CV_SEED)
  folds <- sample(rep_len(seq_len(CV_FOLDS), nrow(reg)))
  tbl <- data.frame(
    model   = names(fits),
    n_terms = sapply(fits, function(f) length(coef(f))),
    cv_R2   = sapply(formulas, cv_r2, data = reg, folds = folds),
    AIC     = sapply(fits, AIC),
    row.names = NULL
  )
  tbl$dAIC_vs_first <- tbl$AIC - tbl$AIC[1]
  message("Model comparison (", CV_FOLDS, "-fold CV R2, AIC)")
  print(tbl)
  invisible(list(model_tbl = tbl, fits = fits, data = reg))
}

vif_table <- function(fit) {
  X <- model.matrix(fit)[, -1, drop = FALSE]
  vifs <- sapply(seq_len(ncol(X)), function(j) {
    r2 <- summary(lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared
    1 / (1 - r2)
  })
  data.frame(term = colnames(X), VIF = round(vifs, 2), row.names = NULL)
}

# Nested model set. Base = covariates only; the others add Pol II, marks, or both.
model_formulas <- function(resp) {
  base_rhs <- "log2_wt_contacts + log2_loop_dist"
  mark_rhs <- paste(MARK_TERMS, collapse = " + ")
  list(
    Base            = as.formula(paste(resp, "~", base_rhs)),
    Pol2            = as.formula(paste(resp, "~", base_rhs, "+ pol2_occ")),
    Marks           = as.formula(paste(resp, "~", base_rhs, "+", mark_rhs)),
    Pol2_plus_Marks = as.formula(paste(resp, "~", base_rhs, "+ pol2_occ +", mark_rhs))
  )
}

# Variance partition from the cross-validated metric of the four models
partition_tbl <- function(cv, metric) {
  inc_pol2  <- unname(cv["Pol2"]            - cv["Base"])
  inc_marks <- unname(cv["Marks"]           - cv["Base"])
  inc_both  <- unname(cv["Pol2_plus_Marks"] - cv["Base"])
  out <- data.frame(
    component = c("Pol II unique (added to marks)",
                  "Histone marks unique (added to Pol II)",
                  "Shared Pol II / marks",
                  "Total beyond Base covariates"),
    delta     = c(inc_both - inc_marks, inc_both - inc_pol2,
                  inc_pol2 + inc_marks - inc_both, inc_both)
  )
  colnames(out)[2] <- paste0("delta_cv_", metric)
  out
}

cv_r2 <- function(formula, data, folds) {
  y    <- data[[all.vars(formula)[1]]]
  pred <- rep(NA_real_, nrow(data))
  for (k in sort(unique(folds))) {
    test       <- folds == k
    fit        <- lm(formula, data = data[!test, ])
    pred[test] <- predict(fit, newdata = data[test, ])
  }
  stopifnot(!anyNA(pred))
  1 - sum((y - pred)^2) / sum((y - mean(y))^2)
}

# AUC = probability that a random positive outranks a random negative
auc <- function(y, p) {
  r  <- rank(p)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cv_auc <- function(formula, data, folds) {
  y    <- data[[all.vars(formula)[1]]]
  pred <- rep(NA_real_, nrow(data))
  for (k in sort(unique(folds))) {
    test       <- folds == k
    fit        <- glm(formula, family = binomial, data = data[!test, ])
    pred[test] <- predict(fit, newdata = data[test, ], type = "response")
  }
  stopifnot(!anyNA(pred))
  auc(y, pred)
}

# ── Linear models for log2FC (Wildtype loop regression and expression fork) ────────────────────────────────
run_lm_models <- function(reg, tag) {
  message("\n══ ", tag, "  (n = ", nrow(reg), " loops) ══════════════════════════")
  stopifnot(nrow(reg) > 20 * length(NUM_TERMS),
            all(c("log2fc", NUM_TERMS) %in% colnames(reg)),
            all(is.finite(as.matrix(reg[, c("log2fc", NUM_TERMS)]))))
  # z-score every numeric predictor so coefficients are comparable
  reg[NUM_TERMS] <- lapply(reg[NUM_TERMS], function(x) as.numeric(scale(x)))
  stopifnot(all(is.finite(as.matrix(reg[, NUM_TERMS]))))
  
  formulas <- model_formulas("log2fc")
  fits     <- lapply(formulas, lm, data = reg)
  
  set.seed(CV_SEED)
  folds <- sample(rep_len(seq_len(CV_FOLDS), nrow(reg)))
  
  model_tbl <- data.frame(
    model   = names(fits),
    n_terms = sapply(fits, function(f) length(coef(f))),
    adj_R2  = sapply(fits, function(f) summary(f)$adj.r.squared),
    cv_R2   = sapply(formulas, cv_r2, data = reg, folds = folds),
    AIC     = sapply(fits, AIC),
    row.names = NULL
  )
  model_tbl$dAIC_vs_Base <- model_tbl$AIC - model_tbl$AIC[model_tbl$model == "Base"]
  message("Model comparison (", CV_FOLDS, "-fold CV R2, AIC)")
  print(model_tbl)
  
  partition <- partition_tbl(setNames(model_tbl$cv_R2, model_tbl$model), "R2")
  message("Variance partition (delta CV R2 over Base)")
  print(partition)
  
  message("Nested F test: does Pol II add to marks?  (Marks vs Pol2_plus_Marks)")
  print(anova(fits$Marks, fits$Pol2_plus_Marks))
  message("Nested F test: do marks add to Pol II?  (Pol2 vs Pol2_plus_Marks)")
  print(anova(fits$Pol2, fits$Pol2_plus_Marks))
  
  message("Collinearity (VIF, Pol2_plus_Marks model)")
  vt <- vif_table(fits$Pol2_plus_Marks)
  print(vt[vt$term %in% NUM_TERMS, ])
  
  message("Standardized coefficients (Pol2_plus_Marks model)")
  print(round(summary(fits$Pol2_plus_Marks)$coefficients, 4))
  
  invisible(list(model_tbl = model_tbl, partition = partition, fits = fits))
}

# ── Logistic models for Lost (1) vs Gained (0) (Lost vs Gained fork) ───────────────────────
run_logit_models <- function(reg, tag) {
  message("\n══ ", tag, "  (Lost = ", sum(reg$lost == 1),
          ", Gained = ", sum(reg$lost == 0), ") ══════════════════════════")
  stopifnot(all(reg$lost %in% c(0, 1)),
            min(table(reg$lost)) >= 2 * CV_FOLDS,
            all(c("lost", NUM_TERMS) %in% colnames(reg)),
            all(is.finite(as.matrix(reg[, NUM_TERMS]))))
  reg[NUM_TERMS] <- lapply(reg[NUM_TERMS], function(x) as.numeric(scale(x)))
  stopifnot(all(is.finite(as.matrix(reg[, NUM_TERMS]))))
  
  formulas <- model_formulas("lost")
  fits     <- lapply(formulas, glm, family = binomial, data = reg)
  
  # folds assigned within each class so every fold holds both classes
  set.seed(CV_SEED)
  folds <- integer(nrow(reg))
  for (cl in c(0, 1)) {
    idx        <- which(reg$lost == cl)
    folds[idx] <- sample(rep_len(seq_len(CV_FOLDS), length(idx)))
  }
  
  model_tbl <- data.frame(
    model   = names(fits),
    n_terms = sapply(fits, function(f) length(coef(f))),
    cv_AUC  = sapply(formulas, cv_auc, data = reg, folds = folds),
    AIC     = sapply(fits, AIC),
    row.names = NULL
  )
  model_tbl$dAIC_vs_Base <- model_tbl$AIC - model_tbl$AIC[model_tbl$model == "Base"]
  message("Model comparison (", CV_FOLDS, "-fold CV AUC, AIC)")
  print(model_tbl)
  
  partition <- partition_tbl(setNames(model_tbl$cv_AUC, model_tbl$model), "AUC")
  message("Variance partition (delta CV AUC over Base)")
  print(partition)
  
  message("Nested likelihood-ratio test: does Pol II add to marks?  (Marks vs Pol2_plus_Marks)")
  print(anova(fits$Marks, fits$Pol2_plus_Marks, test = "Chisq"))
  message("Nested likelihood-ratio test: do marks add to Pol II?  (Pol2 vs Pol2_plus_Marks)")
  print(anova(fits$Pol2, fits$Pol2_plus_Marks, test = "Chisq"))
  
  message("Collinearity (VIF, Pol2_plus_Marks model)")
  vt <- vif_table(fits$Pol2_plus_Marks)
  print(vt[vt$term %in% NUM_TERMS, ])
  
  message("Standardized coefficients (Pol2_plus_Marks model; positive = more likely Lost)")
  print(round(summary(fits$Pol2_plus_Marks)$coefficients, 4))
  
  invisible(list(model_tbl = model_tbl, partition = partition, fits = fits))
}

# ── Save one plot as SVG in PLOT_DIR (no-op when PLOT_DIR is NULL) ────────────
# File name: <assay and loop/map set>_<analysis tag>_<plot type>.svg (words shortened by PLOT_ABBR)
PLOT_ABBR <- c("POLR2A-degraded" = "deg", "Wildtype" = "wt", "log2FC" = "lfc",
               "comparison" = "comp", "coefficients" = "coef",
               "ds_wt_pol2" = "wtpol2", "ds_all" = "all")
plot_run_tag <- function() {
  if (ASSAY == "microc") "microc"
  else if (DS_GROUP == HIC_GROUP) paste0("hic-", DS_GROUP)
  else paste0("hic-loops-", DS_GROUP, "-maps-", HIC_GROUP)
}

save_plot_svg <- function(p, tag, type, w = PLOT_W, h = PLOT_H) {
  if (is.null(PLOT_DIR)) return(invisible(NULL))
  if (!dir.exists(PLOT_DIR)) stop("PLOT_DIR not found: ", PLOT_DIR)
  nm <- paste(plot_run_tag(), tag, type, sep = "_")
  for (k in names(PLOT_ABBR)) nm <- gsub(k, PLOT_ABBR[[k]], nm, ignore.case = TRUE)
  nm <- gsub("^-+|-+$", "", gsub("[^A-Za-z0-9._-]+", "-", nm))
  f  <- file.path(PLOT_DIR, paste0(nm, ".svg"))
  svg(f, width = w, height = h)
  print(p)
  dev.off()
  message("  saved ", f)
  invisible(f)
}

# ── Plots for one set of model results, printed to the graphics device ────────
# 1) cross-validated metric per model  2) variance partition  3) coefficients
plot_model_results <- function(res, tag) {
  mt     <- res$model_tbl
  metric <- intersect(c("cv_R2", "cv_AUC"), colnames(mt))
  stopifnot(length(metric) == 1,
            all(c("Base", "Pol2", "Marks", "Pol2_plus_Marks") %in% mt$model))
  metric_lab <- c(cv_R2 = "Cross-validated R-squared", cv_AUC = "Cross-validated AUC")[[metric]]
  n_lab      <- paste0("N = ", format(nrow(res$fits$Base$model), big.mark = ","), " loops")
  mdl_levels <- c("Base", "Pol2", "Marks", "Pol2_plus_Marks")
  mdl_labels <- c("Base covariates", "Pol II", "Histone marks", "Pol II + marks")
  col_main   <- "#3b6b8f"
  col_ink    <- "#333333"
  th <- theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40"),
          axis.title.x = element_blank())
  
  # 1) model comparison
  d1 <- data.frame(model = factor(mt$model, levels = mdl_levels, labels = mdl_labels),
                   value = mt[[metric]])
  p_cv <- ggplot(d1, aes(model, value)) +
    { if (metric == "cv_R2") geom_hline(yintercept = 0, colour = "grey60") } +
    { if (metric == "cv_R2") geom_segment(aes(xend = model, yend = 0), colour = col_main, linewidth = 0.6) } +
    { if (metric == "cv_AUC") geom_hline(yintercept = 0.5, colour = "grey60", linetype = "dashed") } +
    geom_point(colour = col_main, size = 3.5) +
    geom_text(aes(label = sprintf("%.3f", value)), vjust = -1.1, colour = col_ink, size = 3.6) +
    { if (metric == "cv_AUC") scale_y_continuous(limits = c(0.5, 1), expand = expansion(mult = c(0.02, 0.05)))
      else scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) } +
    labs(title = paste0(tag, ": model comparison"), subtitle = n_lab, y = metric_lab) + th
  print(p_cv)
  save_plot_svg(p_cv, tag, "comp")
  
  # 2) variance partition (gain over the Base covariates)
  part <- res$partition[1:3, ]
  d2 <- data.frame(component = factor(c("Pol II unique", "Histone marks unique", "Shared"),
                                      levels = c("Pol II unique", "Histone marks unique", "Shared")),
                   value = part[[2]])
  p_part <- ggplot(d2, aes(component, value)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_col(fill = col_main, width = 0.55) +
    geom_text(aes(label = sprintf("%.4f", value), vjust = ifelse(value >= 0, -0.5, 1.4)),
              colour = col_ink, size = 3.6) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    labs(title = paste0(tag, ": what each predictor set adds"), subtitle = n_lab,
         y = paste0("Gain in ", tolower(substr(metric_lab, 1, 1)), substring(metric_lab, 2),
                    " over Base covariates")) + th
  print(p_part)
  save_plot_svg(p_part, tag, "gain")
  
  # 3) standardized coefficients, Pol II + marks model
  cf <- summary(res$fits$Pol2_plus_Marks)$coefficients
  stopifnot(all(c("Estimate", "Std. Error") %in% colnames(cf)))
  cf <- cf[rownames(cf) != "(Intercept)", , drop = FALSE]
  term <- rownames(cf)
  d3 <- data.frame(term = term, est = cf[, "Estimate"], se = cf[, "Std. Error"])
  d3$group <- ifelse(grepl("^log2_", d3$term), "Covariates",
                     ifelse(d3$term == "pol2_occ", "Pol II", "Histone marks"))
  d3$group <- factor(d3$group, levels = c("Covariates", "Pol II", "Histone marks"))
  d3$label <- c(log2_wt_contacts = "log2 Wildtype contacts", log2_loop_dist = "log2 loop distance",
                pol2_occ = "Pol II occupancy")[d3$term]
  is_mark <- is.na(d3$label)
  d3$label[is_mark] <- sub("_base$", " baseline", sub("_delta$", " change", d3$term[is_mark]))
  d3$label <- factor(d3$label, levels = rev(d3$label))
  p_coef <- ggplot(d3, aes(est, label)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_linerange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se),
                   orientation = "y", colour = col_main, linewidth = 0.6) +
    geom_point(colour = col_main, size = 3) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = paste0(tag, ": standardized coefficients"),
         subtitle = paste0(n_lab, "; Pol II + marks model, 95% confidence intervals"),
         x = "Coefficient per 1 SD of predictor", y = NULL) +
    th + theme(panel.grid.major.x = element_line(colour = "grey90"),
               panel.grid.major.y = element_blank(),
               axis.title.x = element_text(),
               strip.placement = "outside",
               strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"))
  print(p_coef)
  save_plot_svg(p_coef, tag, "coef")
  
  invisible(list(p_cv = p_cv, p_part = p_part, p_coef = p_coef))
}

# ══════════════════════════════════════════════════════════════════════════════
# Stage A. Loop union, contacts, log2FC, Lost/Gained, anchor signal, expression
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Input checks ────────────────────────────────────────────────────────")
use_ccre <- CCRE_REQUIRED || file.exists(CCRE_FILE)
check_files(c(hic_wt_loops, hic_pol2_loops, hic_wt_file, hic_pol2_file,
              EXPR_FILE, GTF_FILE, BW_FILES, if (use_ccre) CCRE_FILE))
if (!use_ccre) message("cCRE file not found; loop_type will not be built: ", CCRE_FILE)
if (!is.null(HIC_BYTES)) {
  hic_paths <- c(wt = hic_wt_file, pol2 = hic_pol2_file)
  for (n in names(HIC_BYTES)) {
    stopifnot(file.size(hic_paths[[n]]) == HIC_BYTES[[n]],
              STRAW_RESOLUTION %in% strawr::readHicBpResolutions(hic_paths[[n]]))
  }
  message("Micro-C .hic files complete; ", STRAW_RESOLUTION / 1000, " kb resolution present")
}

message("\n── Loading HiCCUPS loops ───────────────────────────────────────────────")
loops_wt <- filter_loops(load_hiccups_bedpe(hic_wt_loops), MIN_OBS_WT) %>%
  dplyr::mutate(loop_dist = abs(centroid2 - centroid1))
loops_pol2 <- filter_loops(load_hiccups_bedpe(hic_pol2_loops), MIN_OBS_POL2) %>%
  dplyr::mutate(loop_dist = abs(centroid2 - centroid1))
check_loops(loops_wt,   "Wildtype")
check_loops(loops_pol2, "POLR2A-degraded")

message("\n── Shared / exclusive loops, union set ─────────────────────────────────")
ov <- find_shared_loops(loops_wt, loops_pol2)
n_dup <- nrow(ov$shared) - dplyr::n_distinct(ov$shared$a_idx)
message("Shared: ", ov$n_shared, "  Wildtype only: ", ov$n_a_only,
        "  POLR2A-degraded only: ", ov$n_b_only,
        "  (Wildtype loops matching >1 degraded loop: ", n_dup, ", counted once)")
union_loops <- dplyr::bind_rows(
  loops_wt   %>% dplyr::slice(unique(ov$shared$a_idx)) %>% dplyr::mutate(source = "shared"),
  loops_wt   %>% dplyr::slice(ov$a_only)               %>% dplyr::mutate(source = "wt_only"),
  loops_pol2 %>% dplyr::slice(ov$b_only)               %>% dplyr::mutate(source = "pol2_only")
)
stopifnot(sum(union_loops$source != "pol2_only") == nrow(loops_wt),
          sum(union_loops$source == "pol2_only") == ov$n_b_only)
union_loops$anchor_bp <- ((union_loops$end1 - union_loops$start1) +
                            (union_loops$end2 - union_loops$start2)) / 2
message("Union loops: ", nrow(union_loops))
message("Mean anchor width (bp), quantiles: ",
        paste(names(quantile(union_loops$anchor_bp, c(0, .25, .5, .75, 1))),
              round(quantile(union_loops$anchor_bp, c(0, .25, .5, .75, 1))), collapse = "  "))

message("\n── Contacts via strawr ─────────────────────────────────────────────────")
union_loops$wt_contacts   <- replace_na(get_contacts(union_loops, hic_wt_file,   "wt"),   0)
union_loops$pol2_contacts <- replace_na(get_contacts(union_loops, hic_pol2_file, "pol2"), 0)
frac_wt   <- mean(union_loops$wt_contacts[union_loops$source != "pol2_only"] > 0)
frac_pol2 <- mean(union_loops$pol2_contacts[union_loops$source != "wt_only"] > 0)
message("Wildtype / shared loops with Wildtype contacts > 0: ", round(100 * frac_wt, 1), "%")
message("Shared / degraded-only loops with degraded contacts > 0: ", round(100 * frac_pol2, 1), "%")
stopifnot(frac_wt > 0.5, frac_pol2 > 0.5)

union_loops <- union_loops %>%
  dplyr::mutate(
    log2fc = log2((pol2_contacts + PSEUDOCOUNT) / (wt_contacts + PSEUDOCOUNT)),
    loop_class = dplyr::case_when(
      wt_contacts   >= MIN_WT_CONTACTS   & log2fc < -LOG2FC_THRESHOLD ~ "Lost",
      pol2_contacts >= MIN_POL2_CONTACTS & log2fc >  LOG2FC_THRESHOLD ~ "Gained",
      TRUE                                                             ~ "Stable"
    ),
    loop_class = factor(loop_class, levels = c("Gained", "Stable", "Lost"))
  )
message("Loop class by source:")
print(table(union_loops$loop_class, union_loops$source))

message("\n── Anchor signal from bigWigs ──────────────────────────────────────────")
loop_raw <- lapply(names(BW_FILES), function(nm) {
  message("  ", nm, " : ", BW_FILES[[nm]])
  a1 <- anchor_signal(BW_FILES[[nm]], union_loops$chr1, union_loops$start1, union_loops$end1)
  a2 <- anchor_signal(BW_FILES[[nm]], union_loops$chr2, union_loops$start2, union_loops$end2)
  (a1 + a2) / 2
})
names(loop_raw) <- names(BW_FILES)
sig_raw <- as.data.frame(loop_raw)

# Structure self-checks on the live data
stopifnot(nrow(sig_raw) == nrow(union_loops),
          all(is.finite(as.matrix(sig_raw))))
track_summary <- data.frame(
  track        = names(sig_raw),
  median       = sapply(sig_raw, median),
  max          = sapply(sig_raw, max),
  frac_nonzero = sapply(sig_raw, function(x) mean(x > 0)),
  row.names    = NULL
)
print(track_summary)
stopifnot(all(track_summary$frac_nonzero > 0))
if (any(track_summary$frac_nonzero < 0.9)) {
  message("CHECK: tracks with <90% of loops having nonzero anchor signal: ",
          paste(track_summary$track[track_summary$frac_nonzero < 0.9], collapse = ", "))
}
for (m in MARKS) {
  d <- sig_raw[[paste0(m, "_DMSO")]]
  i <- sig_raw[[paste0(m, "_IAA")]]
  stopifnot(!isTRUE(all.equal(d, i)))
  message("  ", m, " DMSO vs IAA loop-level correlation: ", round(cor(d, i), 3))
}

pc_pol2  <- pc_for(sig_raw$Pol2)
pc_marks <- sapply(MARKS, function(m) {
  pc_for(c(sig_raw[[paste0(m, "_DMSO")]], sig_raw[[paste0(m, "_IAA")]]))
})
message("Pseudocounts (1st percentile of positive loop-level values):")
print(c(Pol2 = pc_pol2, pc_marks))

features <- union_loops
features$pol2_occ <- log2(sig_raw$Pol2 + pc_pol2)
for (m in MARKS) {
  d <- sig_raw[[paste0(m, "_DMSO")]]
  i <- sig_raw[[paste0(m, "_IAA")]]
  features[[paste0(m, "_base")]]  <- log2(d + pc_marks[[m]])
  features[[paste0(m, "_delta")]] <- log2((i + pc_marks[[m]]) / (d + pc_marks[[m]]))
}
MARK_TERMS <- c(paste0(MARKS, "_base"), paste0(MARKS, "_delta"))
NUM_TERMS  <- c("log2_wt_contacts", "log2_loop_dist", "pol2_occ", MARK_TERMS)
stopifnot(all(c("pol2_occ", MARK_TERMS) %in% colnames(features)),
          all(is.finite(as.matrix(features[, c("pol2_occ", MARK_TERMS)]))))

message("\n── Wildtype expression and anchor genes ────────────────────────────────")
expr_raw <- read_tsv(EXPR_FILE, show_col_types = FALSE)
message("Expression table columns: ", paste(colnames(expr_raw), collapse = ", "))
stopifnot(all(c("gene_id", "TPM") %in% colnames(expr_raw)))
expr <- expr_raw %>%
  dplyr::select(gene_id, TPM) %>%
  dplyr::filter(!grepl("_PAR_Y", gene_id)) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id))
stopifnot(!anyDuplicated(expr$gene_id))

gene_gr <- rtracklayer::import(GTF_FILE, format = "gtf", feature.type = "gene")
stopifnot(length(gene_gr) > 0,
          "gene_id" %in% colnames(mcols(gene_gr)),
          any(grepl("^chr", seqlevels(gene_gr))))
gene_gr$gene_id <- sub("\\..*", "", gene_gr$gene_id)
gene_gr$TPM     <- expr$TPM[match(gene_gr$gene_id, expr$gene_id)]
n_expr_pos      <- sum(expr$TPM > 0, na.rm = TRUE)
n_expr_in_gtf   <- sum(expr$TPM > 0 & expr$gene_id %in% gene_gr$gene_id, na.rm = TRUE)
message("Expressed genes (TPM > 0) in the table: ", n_expr_pos,
        "  found in the GTF: ", n_expr_in_gtf,
        " (", round(100 * n_expr_in_gtf / n_expr_pos, 1), "%)")
stopifnot(n_expr_in_gtf / n_expr_pos > 0.5)
gene_gr <- gene_gr[!is.na(gene_gr$TPM)]

# Tertiles of expressed genes; TPM == 0 counts as Low (same rule as final_APA_analysis.R)
expressed  <- gene_gr$TPM[gene_gr$TPM > 0]
tpm_medium <- quantile(expressed, 1/3)
tpm_high   <- quantile(expressed, 2/3)
message("TPM thresholds: Low <= ", round(tpm_medium, 2), "  Medium <= ", round(tpm_high, 2),
        "  High > ", round(tpm_high, 2))
gene_gr$expr_bin <- dplyr::case_when(
  gene_gr$TPM <= tpm_medium ~ "Low",
  gene_gr$TPM <= tpm_high   ~ "Medium",
  TRUE                      ~ "High"
)
features <- link_loops_to_genes(features, gene_gr)
stopifnot(nrow(features) == nrow(union_loops))
message("Loops by anchor-gene expression bin:")
print(table(features$expr_bin, features$source))

if (use_ccre) {
  message("\n── Promoter / enhancer anchors (cCRE) ──────────────────────────────────")
  ccre_gr  <- load_ccre_gr(CCRE_FILE)
  features <- annotate_anchor_types(features, ccre_gr)
  stopifnot(nrow(features) == nrow(union_loops))
  message("Loops by type (P-P, E-P, E-E, Other) and source:")
  print(table(features$loop_type, features$source))
  message("Anchors by cCRE class (both anchors pooled):")
  print(table(c(features$anchor1_ccre, features$anchor2_ccre)))
} else {
  features$loop_type <- factor(NA_character_, levels = c("P-P", "E-P", "E-E", "Other"))
}

write_tsv(features, file.path(OUT_DIR, "union_loops_features.tsv"))
message("Wrote ", nrow(features), " loops to ",
        file.path(OUT_DIR, "union_loops_features.tsv"))
FEATURES_ASSAY <- ASSAY      # the regression scripts compare this with ASSAY before reusing `features`

# ══════════════════════════════════════════════════════════════════════════════
# Wildtype loop regression set
# ══════════════════════════════════════════════════════════════════════════════
wt_reg <- features %>%
  dplyr::filter(source %in% c("shared", "wt_only"), wt_contacts >= MIN_WT_CONTACTS) %>%
  dplyr::mutate(
    log2_wt_contacts = log2(wt_contacts + PSEUDOCOUNT),
    log2_loop_dist   = log2(loop_dist)
  )
message("\nWildtype loops with WT contacts >= ", MIN_WT_CONTACTS, ": ", nrow(wt_reg),
        " / ", sum(features$source %in% c("shared", "wt_only")))
stopifnot(nrow(wt_reg) > 0)