##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Micro-C loop strength by anchor type, ENCODE POLR2A-degron HCT116
# (Wildtype = ENCSR533VUS, ENCFF534TYD.hic; POLR2A-degraded = ENCSR218URU,
# ENCFF302CUI.hic).
#
# Reads the per-loop anchor annotation exported by
# loop_anchor_analysis.R (Section 3b), extracts the SCALE-normalised
# contact at the 5 kb pixel containing both HiCCUPS loop centroids from each
# condition's own .hic file, and plots loop strength by loop anchor type
# (raincloud as in Hic_loop_analysis/hic_loop_strength_by_anchor_type.R).
# Loops with small anchors only. The black line on each box is the bootstrap
# 95% CI of the median (loops resampled with replacement, N_BOOT resamples).
#
# Contacts are read from the .hic files in tiles: each query returns the block of
# the contact matrix that contains a group of loops (TILE_BP x TILE_BP), not one
# query per loop and not a whole chromosome. A random subset of loops is also
# queried one at a time and compared with the tiled values.
#
# Output (ANALYSIS_DIR): microc_loop_contacts_<condition>.txt (per-condition cache,
# reused if present), microc_loop_contacts.txt (all loops with contacts), and
# figures/microc_loop_strength_by_anchor_type.svg. The plot prints to the live
# device and is saved as SVG.
#
# Run interactively on Quest (needs strawr and ggdist).
##################################################################
suppressPackageStartupMessages({
  library(strawr)
  library(dplyr)
  library(ggplot2)
  library(ggdist)
})

# ══════════════════════════════════════════════════════════════════════════════
# 0. Parameters
# ══════════════════════════════════════════════════════════════════════════════
BASE         <- "/projects/b1042/BackmanLab/Lucas/microc_encode"
ANALYSIS_DIR <- file.path(BASE, "anchor_functional_analysis")
OUT_DIR      <- file.path(ANALYSIS_DIR, "figures")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

LOOP_ANNOT_PATH <- file.path(ANALYSIS_DIR, "microc_loop_anchor_annotations.txt")
CONTACTS_PATH   <- file.path(ANALYSIS_DIR, "microc_loop_contacts.txt")

HIC_FILES <- c(wt   = file.path(BASE, "ENCFF534TYD.hic"),     # ENCSR533VUS
               pol2 = file.path(BASE, "ENCFF302CUI.hic"))     # ENCSR218URU
HIC_BYTES <- c(wt=27301466053, pol2=33785386023)              # ENCODE file sizes
COND_DISPLAY <- c(wt="Wildtype", pol2="POLR2A-degraded")
COND_LEVELS  <- unname(COND_DISPLAY)                          # control first / left

SMALL_ANCHOR_MAX <- 2000L                                     # bp; as in the annotation script
SMALL_LABEL      <- paste0("Small (<=", SMALL_ANCHOR_MAX / 1000, " kb)")
STRAW_NORM       <- "SCALE"
STRAW_RESOLUTION <- 5000L
TILE_BP          <- 5000000L                                  # tile edge, multiple of STRAW_RESOLUTION
N_CHECK_LOOPS    <- 100L                                      # loops re-queried one at a time as a check
slurm_cores      <- Sys.getenv(c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"))
slurm_cores      <- slurm_cores[nzchar(slurm_cores)]
N_CORES          <- if (length(slurm_cores)) as.integer(slurm_cores[1]) else parallel::detectCores()   # allocated cores
N_BOOT           <- 2000L                                     # bootstrap resamples for the median CIs
SEED             <- 1L

ANCHOR_TYPE_COLORS <- c("Exon-Exon"="#d6604d", "Exon-cCRE"="#984ea3", "Exon-Intron"="#4DAF4A",
                        "cCRE-cCRE"="#ff7f00", "Intron-Intron"="#2166ac")   # as in the Hi-C loop-strength figure
LOOP_TYPES   <- c(names(ANCHOR_TYPE_COLORS), "Other")                       # as in the annotation file
OBS_CI_COLOR <- "black"

save_svg <- function(plot_obj, path, width, height) {
  use_svglite <- requireNamespace("svglite", quietly=TRUE)
  if (!use_svglite) message("svglite not installed; saving with grDevices::svg()")
  dev <- if (use_svglite) svglite::svglite else grDevices::svg
  ggplot2::ggsave(path, plot_obj, device=dev, width=width, height=height, units="in")
  stopifnot(file.exists(path), file.size(path) > 0)
  stopifnot(any(grepl("<svg", readLines(path, n=5, warn=FALSE))))
  message("Saved: ", path)
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Loops with anchor annotations (Section 3b export of the annotation script)
# ══════════════════════════════════════════════════════════════════════════════
stopifnot(file.exists(LOOP_ANNOT_PATH))
annot <- read.delim(LOOP_ANNOT_PATH, stringsAsFactors=FALSE)
need  <- c("condition", "anchor_width_class", "chr1", "start1", "end1", "chr2", "start2", "end2",
           "loop_length_bp", "loop_type", "centroid1", "centroid2")
if (!all(need %in% names(annot)))
  stop("Annotation file lacks columns: ", paste(setdiff(need, names(annot)), collapse=", "),
       " (re-run Section 3b of loop_anchor_analysis.R)")
stopifnot(all(annot$condition %in% COND_LEVELS), all(annot$chr1 == annot$chr2),
          !anyNA(annot$centroid1), !anyNA(annot$centroid2), SMALL_LABEL %in% annot$anchor_width_class,
          all(annot$loop_type %in% LOOP_TYPES))
loops_small <- annot %>% dplyr::filter(anchor_width_class == SMALL_LABEL)
message("Loops with small anchors: ", nrow(loops_small), " of ", nrow(annot))
print(table(loops_small$condition, loops_small$loop_type))

# ══════════════════════════════════════════════════════════════════════════════
# 2. SCALE-normalised contact at the pixel containing both centroids, read in tiles
# ══════════════════════════════════════════════════════════════════════════════
for (n in names(HIC_FILES)) {
  stopifnot(file.exists(HIC_FILES[[n]]), file.size(HIC_FILES[[n]]) == HIC_BYTES[[n]])
  stopifnot(STRAW_RESOLUTION %in% strawr::readHicBpResolutions(HIC_FILES[[n]]))
}
message("\nMicro-C .hic files complete; ", STRAW_RESOLUTION / 1000, " kb resolution present")
stopifnot(TILE_BP %% STRAW_RESOLUTION == 0)

# One loop, one query (reference implementation used for the check)
loop_contact_single <- function(chr, lo, hi, hic_file, norm, resolution) {
  mat <- strawr::straw(norm, hic_file,
                       paste0(chr, ":", lo, ":", lo + resolution - 1L),
                       paste0(chr, ":", hi, ":", hi + resolution - 1L), "BP", resolution)
  v <- mat$counts[pmin(mat$x, mat$y) == lo & pmax(mat$x, mat$y) == hi]
  if (length(v) == 1L) v else NA_real_
}

extract_loop_contacts <- function(loops_df, hic_file, norm=STRAW_NORM, resolution=STRAW_RESOLUTION,
                                  tile_bp=TILE_BP, n_cores=N_CORES, n_check=N_CHECK_LOOPS) {
  chroms  <- strawr::readHicChroms(hic_file)
  use_chr <- any(grepl("^chr", chroms$name))
  to_hic  <- function(x) if (use_chr) paste0("chr", gsub("^chr", "", x)) else gsub("^chr", "", x)
  chr     <- to_hic(loops_df$chr1)
  stopifnot(all(chr %in% chroms$name))
  chr_len <- setNames(chroms$length, chroms$name)
  snap    <- function(x) as.integer(floor(x / resolution) * resolution)
  lo      <- snap(pmin(loops_df$centroid1, loops_df$centroid2))
  hi      <- snap(pmax(loops_df$centroid1, loops_df$centroid2))
  ti      <- lo %/% as.integer(tile_bp)
  tj      <- hi %/% as.integer(tile_bp)
  stopifnot(all(ti <= tj))
  n_bins  <- as.integer(tile_bp %/% resolution)
  tiles   <- split(seq_along(lo), paste(chr, ti, tj, sep="|"))
  
  # one tile: the block [ti, tj] of the chromosome's contact matrix, matched to its loops by bin
  query_tile <- function(ids) {
    tryCatch({
      k  <- ids[1]
      x0 <- ti[k] * as.integer(tile_bp); y0 <- tj[k] * as.integer(tile_bp)
      x1 <- min(x0 + as.integer(tile_bp) - 1L, as.integer(chr_len[[chr[k]]]) - 1L)
      y1 <- min(y0 + as.integer(tile_bp) - 1L, as.integer(chr_len[[chr[k]]]) - 1L)
      mat <- strawr::straw(norm, hic_file, paste0(chr[k], ":", x0, ":", x1),
                           paste0(chr[k], ":", y0, ":", y1), "BP", resolution)
      ia  <- (pmin(mat$x, mat$y) - x0) %/% resolution
      ib  <- (pmax(mat$x, mat$y) - y0) %/% resolution
      ok  <- ia >= 0 & ia < n_bins & ib >= 0 & ib < n_bins
      pix <- ia[ok] * n_bins + ib[ok]
      key <- ((lo[ids] - x0) %/% resolution) * n_bins + (hi[ids] - y0) %/% resolution
      list(ids=ids, vals=mat$counts[ok][match(key, pix)], err=NA_character_)
    }, error=function(e) list(ids=ids, vals=rep(NA_real_, length(ids)), err=conditionMessage(e)))
  }
  
  out     <- rep(NA_real_, length(lo))
  n_fail  <- 0L
  chunks  <- split(seq_along(tiles), ceiling(seq_along(tiles) / (n_cores * 4L)))
  t0      <- Sys.time()
  done    <- 0L
  message("  ", length(lo), " loops in ", length(tiles), " tiles of ", tile_bp / 1e6, " Mb; ", n_cores, " cores")
  for (ch in chunks) {
    res <- parallel::mclapply(tiles[ch], query_tile, mc.cores=n_cores, mc.preschedule=FALSE)
    for (r in res) {
      if (!is.list(r) || is.null(r$ids)) { n_fail <- n_fail + 1L; next }
      out[r$ids] <- r$vals
      if (!is.na(r$err)) { n_fail <- n_fail + 1L; message("  tile failed: ", r$err) }
    }
    done <- done + length(ch)
    el   <- as.numeric(difftime(Sys.time(), t0, units="mins"))
    message(sprintf("  tiles %d / %d; %.1f min elapsed; about %.1f min remaining",
                    done, length(tiles), el, el * (length(tiles) - done) / done))
  }
  out[is.nan(out)] <- NA_real_
  stopifnot(length(out) == length(lo))
  message("  failed tiles: ", n_fail, "; loops with a contact value: ", sum(!is.na(out)), " / ", length(out))
  
  # Check on the live file: loops queried one at a time must give the tiled value
  set.seed(SEED)
  chk <- sample.int(length(lo), min(n_check, length(lo)))
  t1  <- Sys.time()
  ref <- vapply(chk, function(i) tryCatch(loop_contact_single(chr[i], lo[i], hi[i], hic_file, norm, resolution),
                                          error=function(e) NA_real_), numeric(1))
  ref[is.nan(ref)] <- NA_real_
  message(sprintf("  single-loop check: %d loops, %.2f s per loop", length(chk),
                  as.numeric(difftime(Sys.time(), t1, units="secs")) / length(chk)))
  stopifnot(identical(is.na(ref), is.na(out[chk])), isTRUE(all.equal(ref[!is.na(ref)], out[chk][!is.na(ref)])))
  out
}

contact_list <- list()
for (cn in names(COND_DISPLAY)) {
  cond  <- COND_DISPLAY[[cn]]
  L     <- loops_small %>% dplyr::filter(condition == cond)
  cache <- file.path(ANALYSIS_DIR, paste0("microc_loop_contacts_", cn, ".txt"))
  keep  <- c("chr1", "start1", "start2", "centroid1", "centroid2")
  if (file.exists(cache)) {
    cached <- read.delim(cache, stringsAsFactors=FALSE)
    if (all(c("norm", "resolution_bp") %in% names(cached)) && nrow(cached) == nrow(L) &&
        all(cached[, keep] == L[, keep]) && all(cached$norm == STRAW_NORM) &&
        all(cached$resolution_bp == STRAW_RESOLUTION)) {
      message("\n", cond, ": using cached contacts ", basename(cache))
      contact_list[[cond]] <- cached
      next
    }
  }
  message("\n", cond, ": extracting ", nrow(L), " loop contacts from ", basename(HIC_FILES[[cn]]), "...")
  L$contacts      <- extract_loop_contacts(L, HIC_FILES[[cn]])
  L$norm          <- STRAW_NORM
  L$resolution_bp <- STRAW_RESOLUTION
  stopifnot(mean(!is.na(L$contacts)) > 0.9)
  write.table(L, cache, sep="\t", quote=FALSE, row.names=FALSE)
  contact_list[[cond]] <- L
}
contacts_all <- dplyr::bind_rows(contact_list)
stopifnot(nrow(contacts_all) == nrow(loops_small))
write.table(contacts_all, CONTACTS_PATH, sep="\t", quote=FALSE, row.names=FALSE)
message("\nSaved: ", CONTACTS_PATH, " (", nrow(contacts_all), " loops)")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Loop strength by anchor type
# ══════════════════════════════════════════════════════════════════════════════
strength_df <- contacts_all %>%
  dplyr::filter(loop_type %in% names(ANCHOR_TYPE_COLORS), !is.na(contacts), contacts > 0) %>%
  dplyr::mutate(cond=factor(condition, levels=COND_LEVELS),
                loop_anchor_type=factor(loop_type, levels=names(ANCHOR_TYPE_COLORS)),
                log10_observed=pmin(log10(contacts + 1), log10(1001)))
counts <- strength_df %>% dplyr::count(cond, loop_anchor_type) %>%
  dplyr::mutate(label=paste0(loop_anchor_type, "\n(N=", format(n, big.mark=",", trim=TRUE), ")"))
strength_df <- strength_df %>%
  dplyr::left_join(counts %>% dplyr::select(cond, loop_anchor_type, label), by=c("cond", "loop_anchor_type"))
strength_df$label <- factor(strength_df$label, levels=unique(counts$label[order(counts$cond, counts$loop_anchor_type)]))

# Bootstrap 95% CI of the median (loops resampled with replacement)
boot_median_ci <- function(x, n_boot=N_BOOT)
  unname(quantile(replicate(n_boot, median(sample(x, length(x), replace=TRUE))), c(0.025, 0.975)))
set.seed(SEED)
median_ci <- strength_df %>% dplyr::group_by(cond, loop_anchor_type, label) %>%
  dplyr::summarise(N=dplyr::n(), median_log10=median(log10_observed),
                   ci=list(boot_median_ci(log10_observed)), .groups="drop") %>%
  dplyr::mutate(ci_lo=vapply(ci, `[`, numeric(1), 1), ci_hi=vapply(ci, `[`, numeric(1), 2))
stopifnot(all(median_ci$ci_lo <= median_ci$median_log10 + 1e-9), all(median_ci$ci_hi >= median_ci$median_log10 - 1e-9))
message("\nLoop strength (", STRAW_NORM, " contacts) by anchor type; median and bootstrap 95% CI:")
print(median_ci %>% dplyr::transmute(cond, loop_anchor_type, N,
                                     median_contacts=round(10^median_log10 - 1, 2),
                                     ci_lo_contacts=round(10^ci_lo - 1, 2), ci_hi_contacts=round(10^ci_hi - 1, 2)))

fig_strength <- ggplot(strength_df, aes(x=label, y=log10_observed, fill=loop_anchor_type, color=loop_anchor_type)) +
  ggdist::stat_halfeye(aes(thickness=after_stat(pdf * n)), side="right", scale=0.4, width=0.6,
                       .width=0, point_colour=NA, alpha=0.75, normalize="groups",
                       position=position_nudge(x=0.05)) +
  geom_boxplot(width=0.12, outlier.shape=NA, color="gray20", alpha=0.9, linewidth=0.6,
               position=position_nudge(x=-0.07)) +
  geom_linerange(data=median_ci, aes(x=label, ymin=ci_lo, ymax=ci_hi), inherit.aes=FALSE,
                 color=OBS_CI_COLOR, linewidth=1.2, position=position_nudge(x=-0.07)) +
  facet_wrap(~cond, scales="free_x") +
  scale_fill_manual(values=ANCHOR_TYPE_COLORS, guide="none") +
  scale_color_manual(values=ANCHOR_TYPE_COLORS, guide="none") +
  scale_y_continuous(name=paste0(STRAW_NORM, " contacts (", STRAW_RESOLUTION / 1000, "kb, log10)"),
                     limits=c(NA, log10(1001) * 1.05), breaks=c(0, 1, 2, 3),
                     labels=c("1", "10", "100", "1000"), expand=expansion(mult=c(0.02, 0))) +
  scale_x_discrete(expand=expansion(add=c(0.3, 0.5))) +
  labs(title=paste0("Micro-C loop strength by anchor type (anchors <=", SMALL_ANCHOR_MAX / 1000, " kb)"),
       subtitle="Box: median and quartiles. Black line: bootstrap 95% CI of the median (over loops).", x=NULL) +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=14), plot.subtitle=element_text(size=9),
        axis.text.x=element_text(size=9, angle=30, hjust=1, vjust=1),
        strip.background=element_blank(), strip.text=element_text(face="bold"),
        plot.margin=margin(8, 16, 8, 8))
print(fig_strength)
save_svg(fig_strength, file.path(OUT_DIR, "microc_loop_strength_by_anchor_type.svg"), width=11, height=6)