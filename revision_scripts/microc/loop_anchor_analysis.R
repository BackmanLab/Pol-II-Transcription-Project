##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# ENCODE POLR2A-degron HCT116 Micro-C: do loop anchors fall on functional
# elements (exons and cCREs)?
# Wildtype = ENCSR533VUS (untreated); POLR2A-degraded = ENCSR218URU
# (6 h 5-Ph-IAA). HiCCUPS loop calls as released by ENCODE.
#
# 1. Anchor classes (Exon+cCRE / Exon / cCRE / Intron / Intergenic) and loop
#    types (Exon-Exon / Exon-cCRE / Exon-Intron / cCRE-cCRE / Intron-Intron /
#    Other), split by anchor width (small <= SMALL_ANCHOR_MAX vs larger).
#    cCREs inside introns count as cCRE anchors; "Intron" = gene body with no
#    exon and no cCRE.
# 2. Background model: each loop is moved as a unit to a random position on its
#    own chromosome (same anchor width, same anchor separation, same bin grid;
#    N_PERM permutations). Overlap is pre-computed once for every grid bin, so
#    each permutation is a lookup.
#    Statistics:
#      A  anchors overlapping an exon
#      B  exon overlap among anchors overlapping a gene body
#      C  loops with both anchors on exons
#      D  anchors overlapping a cCRE
#      E  anchors on a functional element (exon or cCRE)
#      F  loops with both anchors functional
#      and the fraction of loops in each loop type.
#    Empirical one-sided p = (1 + #(null >= observed)) / (N_PERM + 1)
#    (p_depleted uses <=).
# 3. Loop positions with anchor annotations are exported (Section 3b) for
#    loop_anchor_contact_strength_analysis.R, which extracts the SCALE-normalised
#    Micro-C contacts at the loop centroid pixel and plots loop strength by
#    anchor type (run it after this script, interactively).
# 4. Confidence intervals: every bar plot shows a bootstrap 95% CI of the observed
#    value (loops resampled with replacement, N_BOOT resamples; both anchors of a
#    loop stay together) next to the random-placement mean and 95% range.
#
# Annotation: UCSC hg38 knownGene TxDb (strand-collapsed, loop anchors are
# unstranded); cCREs from the ENCODE Registry-V4 bed.
#
# Run interactively on Quest. Plots print to the live device and are
# also saved as SVG to <ANALYSIS_DIR>/figures/.
##################################################################
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(dplyr)
  library(ggplot2)
})
library(GenomicFeatures)                       # ANNOTATION_SOURCE
library(TxDb.Hsapiens.UCSC.hg38.knownGene)     # ANNOTATION_SOURCE

# ══════════════════════════════════════════════════════════════════════════════
# 0. Parameters
# ══════════════════════════════════════════════════════════════════════════════
BASE         <- "/projects/b1042/BackmanLab/Lucas/microc_encode"
ANALYSIS_DIR <- file.path(BASE, "anchor_functional_analysis")
OUT_DIR      <- file.path(ANALYSIS_DIR, "figures")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

LOOP_FILES <- c(
  wt   = file.path(BASE, "ENCFF954EMJ.bedpe.gz"),   # ENCSR533VUS
  pol2 = file.path(BASE, "ENCFF604OYY.bedpe.gz")    # ENCSR218URU
)
COND_DISPLAY <- c(wt="Wildtype", pol2="POLR2A-degraded")
COND_LEVELS  <- unname(COND_DISPLAY)                 # control first / left
COND_COLORS  <- setNames(c("#88CCEE", "#CC6677"), COND_LEVELS)

SMALL_ANCHOR_MAX <- 2000L    # bp; anchors <= this are "small"
MAX_FDR          <- 0.1
MIN_OBSERVED     <- 0
N_PERM           <- 10000L
N_BOOT           <- 2000L    # bootstrap resamples of loops for the observed-value CIs
SEED             <- 1L
STD_CHROMS       <- paste0("chr", c(1:22, "X"))

WIDTH_LEVELS <- c(paste0("Small (<=", SMALL_ANCHOR_MAX / 1000, " kb)"),
                  paste0("Larger (>",  SMALL_ANCHOR_MAX / 1000, " kb)"))

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load HiCCUPS loop calls (bedpe, as in Hic_loop_analysis scripts)
# ══════════════════════════════════════════════════════════════════════════════
load_hiccups_bedpe <- function(path) {
  all_lines   <- readLines(path)
  header_line <- all_lines[grepl("^#chr1", all_lines)][1]
  stopifnot(!is.na(header_line))
  col_names   <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  data_lines  <- all_lines[!grepl("^#", all_lines)]
  if (length(data_lines) == 0) stop("No data lines found in ", path)
  
  df <- read.table(text=paste(data_lines, collapse="\n"), header=FALSE, sep="\t",
                   stringsAsFactors=FALSE, fill=TRUE, quote="", comment.char="")
  stopifnot(ncol(df) == length(col_names))
  colnames(df) <- col_names
  need <- c("chr1", "x1", "x2", "chr2", "y1", "y2", "observed",
            "fdrBL", "fdrDonut", "fdrH", "fdrV", "centroid1", "centroid2")
  stopifnot(all(need %in% colnames(df)))
  
  out <- data.frame(
    chr1=as.character(df$chr1), x1=as.integer(df$x1), x2=as.integer(df$x2),
    chr2=as.character(df$chr2), y1=as.integer(df$y1), y2=as.integer(df$y2),
    observed=as.numeric(df$observed), fdrDonut=as.numeric(df$fdrDonut),
    fdrBL=as.numeric(df$fdrBL), fdrH=as.numeric(df$fdrH), fdrV=as.numeric(df$fdrV),
    centroid1=as.integer(df$centroid1), centroid2=as.integer(df$centroid2),
    stringsAsFactors=FALSE
  ) %>% dplyr::filter(!is.na(x1), !is.na(y1))
  message("Loaded ", nrow(out), " loops from ", basename(path))
  out
}

filter_loops <- function(df, min_observed=MIN_OBSERVED, max_fdr=MAX_FDR) {
  df %>% dplyr::filter(observed >= min_observed, fdrDonut <= max_fdr, fdrBL <= max_fdr,
                       fdrH <= max_fdr, fdrV <= max_fdr)
}

loops_list <- lapply(names(LOOP_FILES), function(n) {
  raw <- load_hiccups_bedpe(LOOP_FILES[[n]])
  flt <- filter_loops(raw)
  message("  ", COND_DISPLAY[[n]], ": ", nrow(flt), " of ", nrow(raw), " loops pass filters")
  flt <- flt %>% dplyr::filter(chr1 %in% STD_CHROMS, chr2 %in% STD_CHROMS)
  flt$cond <- COND_DISPLAY[[n]]
  flt
})
loops_all <- dplyr::bind_rows(loops_list)

# Structure checks on the live loop tables
stopifnot(all(loops_all$chr1 == loops_all$chr2))                    # intrachromosomal
stopifnot(all((loops_all$x2 - loops_all$x1) == (loops_all$y2 - loops_all$y1)))  # equal-width anchors
stopifnot(!anyNA(loops_all$centroid1), !anyNA(loops_all$centroid2))
cen_dist <- pmax(pmax(loops_all$x1 - loops_all$centroid1, loops_all$centroid1 - loops_all$x2, 0),
                 pmax(loops_all$y1 - loops_all$centroid2, loops_all$centroid2 - loops_all$y2, 0))
message("Loops with both centroids inside their anchors: ", round(100 * mean(cen_dist == 0), 1), "%")
message("Centroid distance outside anchor (bp), quantiles:")
print(quantile(cen_dist, c(0.5, 0.8, 0.9, 0.95, 0.99, 1)))
loops_all$width       <- loops_all$x2 - loops_all$x1
loops_all$width_class <- ifelse(loops_all$width <= SMALL_ANCHOR_MAX, WIDTH_LEVELS[1], WIDTH_LEVELS[2])
loops_all$cond        <- factor(loops_all$cond, levels=COND_LEVELS)
loops_all$width_class <- factor(loops_all$width_class, levels=WIDTH_LEVELS)
stopifnot(!anyNA(loops_all$cond))
message("\nLoops by anchor width (bp) and condition:")
print(table(loops_all$width, loops_all$cond))
stopifnot(any(loops_all$width_class == WIDTH_LEVELS[1]))            # small anchors exist

# ══════════════════════════════════════════════════════════════════════════════
# 2. Exon / intron / gene-body ranges (UCSC hg38 knownGene, as in the
#    anchor-type scripts), strand-collapsed
# ══════════════════════════════════════════════════════════════════════════════
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene      # ANNOTATION_SOURCE

exon_gr <- GenomicRanges::reduce(exons(txdb), ignore.strand=TRUE)
# Gene bodies = transcript spans: exons of transcripts without an Entrez Gene ID
# are in exons(txdb) but not in genes(txdb)
gene_gr <- GenomicRanges::reduce(transcripts(txdb), ignore.strand=TRUE)
exon_gr <- keepSeqlevels(exon_gr, intersect(seqlevels(exon_gr), STD_CHROMS), pruning.mode="coarse")
gene_gr <- keepSeqlevels(gene_gr, intersect(seqlevels(gene_gr), STD_CHROMS), pruning.mode="coarse")
intron_gr <- GenomicRanges::setdiff(gene_gr, exon_gr)

seqlens <- seqlengths(txdb)[STD_CHROMS]
stopifnot(!anyNA(seqlens))
stopifnot(all(overlapsAny(exon_gr, gene_gr, ignore.strand=TRUE)))   # exons lie inside gene bodies
message("\nAnnotation: ", length(exon_gr), " merged exon ranges, ",
        length(gene_gr), " merged gene bodies, ", length(intron_gr), " intron ranges")
message("Exon bp / gene-body bp / genome bp: ",
        sum(as.numeric(width(exon_gr))), " / ", sum(as.numeric(width(gene_gr))), " / ",
        sum(as.numeric(seqlens)))

# Fast interval-overlap test against a set of disjoint ranges (global coordinates)
make_overlap_fn <- function(gr, chroms, chrom_lengths) {
  gr      <- GenomicRanges::reduce(gr, ignore.strand=TRUE)
  offsets <- setNames(c(0, cumsum(as.numeric(chrom_lengths)))[seq_along(chroms)], chroms)
  chr_gr  <- as.character(seqnames(gr))
  r_start <- unname(offsets[chr_gr]) + start(gr)
  r_end   <- unname(offsets[chr_gr]) + end(gr)
  o       <- order(r_start)
  r_start <- r_start[o]; r_end <- r_end[o]
  stopifnot(all(r_start[-1] > r_end[-length(r_end)]))
  function(chr, start, end) {
    g_start <- unname(offsets[chr]) + start
    g_end   <- unname(offsets[chr]) + end
    idx     <- findInterval(g_end, r_start)
    idx > 0L & r_end[pmax(idx, 1L)] >= g_start
  }
}
overlap_exon   <- make_overlap_fn(exon_gr,   STD_CHROMS, seqlens)
overlap_gene   <- make_overlap_fn(gene_gr,   STD_CHROMS, seqlens)
overlap_intron <- make_overlap_fn(intron_gr, STD_CHROMS, seqlens)

make_anchor_set <- function(L) {
  # bedpe is 0-based half-open -> 1-based closed
  list(chr=c(L$chr1, L$chr2), start=c(L$x1, L$y1) + 1L, end=c(L$x2, L$y2), n_loops=nrow(L))
}

# ══════════════════════════════════════════════════════════════════════════════
# 2b. cCREs, functional-element definition (exon OR cCRE) and loop-type classes
#     Registry-V4 cCRE bed; cCRE classes and anchor/loop classification as in
#     hic_loop_strength_by_anchor_type.R. A cCRE inside an intron counts as a
#     cCRE (positive) anchor, so "Intron" means gene body with no exon and no cCRE.
# ══════════════════════════════════════════════════════════════════════════════
CCRE_PATH  <- file.path(ANALYSIS_DIR, "ENCFF431JDU_ENCFF964OOU_ENCFF787LMI_ENCFF388PVO.bed")   # Registry-V4
stopifnot(file.exists(CCRE_PATH))
CCRE_TYPES <- c("PLS", "pELS", "dELS", "CA-CTCF")

CLASS_ORDER <- c("Exon+cCRE", "Exon", "cCRE", "Intron", "Intergenic")      # integer codes 1..5
LOOP_TYPES  <- c("Exon-Exon", "Exon-cCRE", "Exon-Intron", "cCRE-cCRE", "Intron-Intron", "Other")
STAT_LABELS <- c(A="Anchors overlapping an exon",
                 B="Exon overlap among gene-body anchors",
                 C="Loops with both anchors on exons",
                 D="Anchors overlapping a cCRE",
                 E="Anchors on a functional element (exon or cCRE)",
                 F="Loops with both anchors functional",
                 setNames(paste0("Loops: ", LOOP_TYPES), LOOP_TYPES))
STAT_IDS    <- names(STAT_LABELS)

load_ccre_gr <- function(path, types) {
  df <- read.table(path, header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char="", fill=TRUE)
  message("cCRE file: ", nrow(df), " rows x ", ncol(df), " columns; first row: ", paste(unlist(df[1, ]), collapse=" | "))
  if (!is.numeric(df$V2) || !is.numeric(df$V3)) stop("V2/V3 are not numeric (header line in file?)")
  if (!any(df$V1 %in% STD_CHROMS)) stop("V1 holds no chr-prefixed chromosome names")
  # class column = the one column whose values start with a cCRE class (value before any comma)
  class_cols <- names(df)[vapply(df, function(x) is.character(x) && any(sub(",.*", "", x) %in% types), logical(1))]
  if (length(class_cols) != 1) stop("Expected exactly one cCRE class column, found: ", paste(class_cols, collapse=", "))
  message("cCRE class column: ", class_cols)
  print(table(df[[class_cols]]))
  df$ccre_type <- sub(",.*", "", df[[class_cols]])
  df <- df[df$ccre_type %in% types & df$V1 %in% STD_CHROMS, ]
  stopifnot(nrow(df) > 0)
  gr <- GRanges(df$V1, IRanges(df$V2 + 1L, df$V3), ccre_type=df$ccre_type)   # BED 0-based -> 1-based
  stopifnot(all(end(gr) <= seqlens[as.character(seqnames(gr))]))        # same genome build as loops/txdb
  gr
}
ccre_all  <- load_ccre_gr(CCRE_PATH, CCRE_TYPES)
ccre_gr   <- GenomicRanges::reduce(ccre_all, ignore.strand=TRUE)
func_gr   <- GenomicRanges::reduce(c(granges(exon_gr), granges(ccre_gr)), ignore.strand=TRUE)
overlap_ccre <- make_overlap_fn(ccre_gr, STD_CHROMS, seqlens)

genome_bp <- sum(as.numeric(seqlens))
cover <- data.frame(
  element=c("Exon", "cCRE", "Exon or cCRE", "Gene body"),
  bp=c(sum(as.numeric(width(exon_gr))), sum(as.numeric(width(ccre_gr))),
       sum(as.numeric(width(func_gr))), sum(as.numeric(width(gene_gr)))))
cover$pct_of_genome <- round(100 * cover$bp / genome_bp, 2)
message("\ncCRE types kept: ", paste(CCRE_TYPES, collapse=", "), "; ", length(ccre_all),
        " cCREs (", length(ccre_gr), " after merging)")
message("Genome coverage of annotation classes (denominator: chr1-22,X incl. gaps):")
print(cover)

anchor_hits <- function(chr, start, end) {
  list(exon=overlap_exon(chr, start, end), gene=overlap_gene(chr, start, end),
       ccre=overlap_ccre(chr, start, end))
}

# Priority classes: Exon+cCRE (1) > Exon (2) > cCRE (3) > Intron (4) > Intergenic (5)
anchor_class_code <- function(hit_exon, hit_gene, hit_ccre) {
  ifelse(hit_exon & hit_ccre, 1L, ifelse(hit_exon, 2L, ifelse(hit_ccre, 3L, ifelse(hit_gene, 4L, 5L))))
}
# Loop type from the two anchor class codes (same rules as classify_hic_loop_anchors_wt)
loop_type_code <- function(c1, c2) {
  e1 <- c1 <= 2L; e2 <- c2 <= 2L; k1 <- c1 == 3L; k2 <- c2 == 3L; i1 <- c1 == 4L; i2 <- c2 == 4L
  out <- rep(6L, length(c1))                                  # Other
  out[i1 & i2]              <- 5L                             # Intron-Intron
  out[k1 & k2]              <- 4L                             # cCRE-cCRE
  out[(e1 & i2) | (i1 & e2)] <- 3L                            # Exon-Intron
  out[(e1 & k2) | (k1 & e2)] <- 2L                            # Exon-cCRE
  out[e1 & e2]              <- 1L                             # Exon-Exon
  out
}
# Anchors are ordered: first n_loops = anchor 1, next n_loops = anchor 2
anchor_stats <- function(hit_exon, hit_gene, hit_ccre, n_loops) {
  first  <- seq_len(n_loops)
  second <- n_loops + first
  func   <- hit_exon | hit_ccre
  code   <- anchor_class_code(hit_exon, hit_gene, hit_ccre)
  ltype  <- loop_type_code(code[first], code[second])
  c(A=mean(hit_exon),
    B=sum(hit_exon) / sum(hit_gene),
    C=mean(hit_exon[first] & hit_exon[second]),
    D=mean(hit_ccre),
    E=mean(func),
    F=mean(func[first] & func[second]),
    setNames(tabulate(ltype, nbins=length(LOOP_TYPES)) / n_loops, LOOP_TYPES))
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Item 1 — observed anchor annotation
# ══════════════════════════════════════════════════════════════════════════════
obs_sets <- list()
for (cn in COND_LEVELS) {
  for (wc in WIDTH_LEVELS) {
    L <- loops_all %>% dplyr::filter(cond == cn, width_class == wc)
    if (nrow(L) == 0) next
    av <- make_anchor_set(L)
    h  <- anchor_hits(av$chr, av$start, av$end)
    code  <- anchor_class_code(h$exon, h$gene, h$ccre)
    ltype <- loop_type_code(code[seq_len(av$n_loops)], code[av$n_loops + seq_len(av$n_loops)])
    obs_sets[[paste(cn, wc, sep="|")]] <- list(
      cond=cn, width_class=wc, loops=L, av=av, hit=h, code=code, loop_type=ltype,
      stats=anchor_stats(h$exon, h$gene, h$ccre, av$n_loops))
  }
}

# Structure self-checks on the live anchors
for (key in names(obs_sets)) {
  s  <- obs_sets[[key]]
  gr <- GRanges(s$av$chr, IRanges(s$av$start, s$av$end))
  stopifnot(identical(s$hit$exon, as.vector(suppressWarnings(overlapsAny(gr, exon_gr, ignore.strand=TRUE)))))
  stopifnot(identical(s$hit$gene, as.vector(suppressWarnings(overlapsAny(gr, gene_gr, ignore.strand=TRUE)))))
  stopifnot(identical(s$hit$ccre, as.vector(suppressWarnings(overlapsAny(gr, ccre_all, ignore.strand=TRUE)))))
  stopifnot(identical(s$hit$exon | s$hit$ccre,
                      as.vector(suppressWarnings(overlapsAny(gr, func_gr, ignore.strand=TRUE)))))
  stopifnot(all((s$hit$exon | overlap_intron(s$av$chr, s$av$start, s$av$end)) == s$hit$gene))
  stopifnot(all(s$av$chr %in% STD_CHROMS), all(s$av$end <= seqlens[s$av$chr]))
  # Loop types agree with the string-based rules of classify_hic_loop_anchors_wt
  cls   <- c("exon+cCRE", "exon", "cCRE", "intron", "intergenic")[s$code]
  a1    <- cls[seq_len(s$av$n_loops)]; a2 <- cls[s$av$n_loops + seq_len(s$av$n_loops)]
  ref   <- dplyr::case_when(
    grepl("exon", a1) & grepl("exon", a2) ~ "Exon-Exon",
    (grepl("exon", a1) & a2 == "cCRE") | (a1 == "cCRE" & grepl("exon", a2)) ~ "Exon-cCRE",
    (grepl("exon", a1) & a2 == "intron") | (a1 == "intron" & grepl("exon", a2)) ~ "Exon-Intron",
    a1 == "cCRE" & a2 == "cCRE" ~ "cCRE-cCRE",
    a1 == "intron" & a2 == "intron" ~ "Intron-Intron",
    TRUE ~ "Other")
  stopifnot(identical(LOOP_TYPES[s$loop_type], ref))
  stopifnot(abs(sum(s$stats[LOOP_TYPES]) - 1) < 1e-9)
}
message("Structure checks passed for ", length(obs_sets), " condition x width sets")

anchor_summary <- dplyr::bind_rows(lapply(obs_sets, function(s) {
  n   <- length(s$code)
  out <- data.frame(cond=s$cond, width_class=s$width_class, n_loops=s$av$n_loops, n_anchors=n,
                    stringsAsFactors=FALSE)
  for (k in seq_along(CLASS_ORDER)) out[[paste0("pct_", CLASS_ORDER[k])]] <- 100 * mean(s$code == k)
  out$pct_functional <- 100 * mean(s$code <= 3L)
  stopifnot(abs(out$pct_functional - 100 * unname(s$stats["E"])) < 1e-9)   # matches statistic E
  out
}))
rownames(anchor_summary) <- NULL
message("\n── Item 1: observed anchor classes (%) ─────────────────────────────")
print(anchor_summary %>% dplyr::mutate(dplyr::across(dplyr::starts_with("pct_"), ~round(.x, 1))))

loop_type_summary <- dplyr::bind_rows(lapply(obs_sets, function(s)
  data.frame(cond=s$cond, width_class=s$width_class, loop_type=LOOP_TYPES,
             n=tabulate(s$loop_type, nbins=length(LOOP_TYPES)),
             pct=100 * tabulate(s$loop_type, nbins=length(LOOP_TYPES)) / s$av$n_loops,
             stringsAsFactors=FALSE)))
message("\nObserved loops by anchor-type pair:")
print(loop_type_summary %>% dplyr::mutate(pct=round(pct, 1)))

# ══════════════════════════════════════════════════════════════════════════════
# 3b. Export: loop positions with anchor annotations (tab-delimited, one row per loop,
#     both conditions and anchor-width classes). Coordinates as in the ENCODE bedpe
#     (0-based, half-open).
# ══════════════════════════════════════════════════════════════════════════════
LOOP_ANNOT_PATH <- file.path(ANALYSIS_DIR, "microc_loop_anchor_annotations.txt")
loop_annot <- dplyr::bind_rows(lapply(obs_sets, function(s) {
  n <- s$av$n_loops
  f <- seq_len(n)
  data.frame(condition=as.character(s$cond), anchor_width_class=as.character(s$width_class),
             chr1=s$loops$chr1, start1=s$loops$x1, end1=s$loops$x2,
             chr2=s$loops$chr2, start2=s$loops$y1, end2=s$loops$y2,
             centroid1=s$loops$centroid1, centroid2=s$loops$centroid2,
             loop_length_bp=s$loops$y1 - s$loops$x1,
             loop_type=LOOP_TYPES[s$loop_type],
             anchor1_class=CLASS_ORDER[s$code[f]], anchor2_class=CLASS_ORDER[s$code[n + f]],
             anchor1_exon=s$hit$exon[f], anchor1_gene_body=s$hit$gene[f], anchor1_ccre=s$hit$ccre[f],
             anchor2_exon=s$hit$exon[n + f], anchor2_gene_body=s$hit$gene[n + f], anchor2_ccre=s$hit$ccre[n + f],
             stringsAsFactors=FALSE)
}))
loop_annot <- loop_annot[order(match(loop_annot$condition, COND_LEVELS), match(loop_annot$chr1, STD_CHROMS),
                               loop_annot$start1, loop_annot$start2), ]
rownames(loop_annot) <- NULL

# Structure checks on the live table
stopifnot(nrow(loop_annot) == nrow(loops_all), !anyNA(loop_annot))
stopifnot(all(loop_annot$end1 - loop_annot$start1 == loop_annot$end2 - loop_annot$start2))
chk_n <- loop_annot %>% dplyr::count(condition, anchor_width_class, loop_type) %>%
  dplyr::inner_join(loop_type_summary %>% dplyr::filter(n > 0) %>%
                      dplyr::transmute(condition=as.character(cond), anchor_width_class=as.character(width_class),
                                       loop_type, n_ref=n),
                    by=c("condition", "anchor_width_class", "loop_type"))
stopifnot(nrow(chk_n) == sum(loop_type_summary$n > 0), all(chk_n$n == chk_n$n_ref))   # same counts as the loop-type table

write.table(loop_annot, LOOP_ANNOT_PATH, sep="\t", quote=FALSE, row.names=FALSE)
re_read <- read.delim(LOOP_ANNOT_PATH, stringsAsFactors=FALSE)
stopifnot(nrow(re_read) == nrow(loop_annot), identical(names(re_read), names(loop_annot)),
          all(re_read$start1 == loop_annot$start1), all(re_read$loop_type == loop_annot$loop_type))
message("\nSaved: ", LOOP_ANNOT_PATH, " (", nrow(loop_annot), " loops x ", ncol(loop_annot), " columns)")
print(utils::head(loop_annot, 3))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Item 2 — background model. Each loop is moved as a unit to a random position
#    on its own chromosome (same anchor width, same anchor-to-anchor distance,
#    same bin grid), so both anchors keep their separation.
# ══════════════════════════════════════════════════════════════════════════════
# Pre-computed overlap tables: every grid bin (fully inside its chromosome) for each
# anchor width present in the loop sets
bin_rows <- expand.grid(width=sort(unique(loops_all$width)), chr=STD_CHROMS, stringsAsFactors=FALSE)
bin_rows$n_bins <- floor(unname(seqlens[bin_rows$chr]) / bin_rows$width)
bin_rows$offset <- c(0, cumsum(bin_rows$n_bins))[seq_len(nrow(bin_rows))]   # 0-based offset
bin_rows$key    <- paste(bin_rows$width, bin_rows$chr)
bin_chr   <- rep(bin_rows$chr, bin_rows$n_bins)
bin_width <- rep(as.numeric(bin_rows$width), bin_rows$n_bins)
bin_k     <- sequence(bin_rows$n_bins) - 1
bin_start <- bin_k * bin_width + 1
bin_end   <- (bin_k + 1) * bin_width
bin_exon  <- overlap_exon(bin_chr, bin_start, bin_end)
bin_gene  <- overlap_gene(bin_chr, bin_start, bin_end)
bin_ccre  <- overlap_ccre(bin_chr, bin_start, bin_end)
message("\nOverlap tables: ", length(bin_exon), " grid bins across ",
        length(unique(bin_rows$width)), " anchor widths")
rm(bin_chr, bin_width, bin_k, bin_start, bin_end)

# Table row (offset, number of bins) for each anchor of a set
anchor_bin_rows <- function(av) {
  width  <- av$end - av$start + 1L
  row_id <- match(paste(width, av$chr), bin_rows$key)
  stopifnot(!anyNA(row_id))
  list(base=bin_rows$offset[row_id], nb=bin_rows$n_bins[row_id], width=width)
}

# Live checks: observed anchors sit on the grid, and the table lookup reproduces the
# direct overlap calls for every observed anchor
for (key in names(obs_sets)) {
  s  <- obs_sets[[key]]
  br <- anchor_bin_rows(s$av)
  stopifnot(all((s$av$start - 1L) %% br$width == 0))                 # anchors on the bin grid
  k_obs <- (s$av$start - 1L) %/% br$width
  stopifnot(all(k_obs < br$nb))
  stopifnot(identical(bin_exon[br$base + k_obs + 1], s$hit$exon))
  stopifnot(identical(bin_gene[br$base + k_obs + 1], s$hit$gene))
  stopifnot(identical(bin_ccre[br$base + k_obs + 1], s$hit$ccre))
}
message("Grid and overlap-table checks passed for ", length(obs_sets), " condition x width sets")

simulate_null <- function(av, n_perm) {
  br     <- anchor_bin_rows(av)
  n_loop <- av$n_loops
  first  <- seq_len(n_loop)
  second <- n_loop + first
  k_obs  <- (av$start - 1L) %/% br$width
  span   <- k_obs[second] - k_obs[first]                              # anchor separation in bins
  stopifnot(all(span >= 0), all(br$base[first] == br$base[second]),
            all(br$nb[first] == br$nb[second]), all(br$nb[first] > span))
  n_pos  <- br$nb[first] - span                                       # positions keeping both anchors on the chromosome
  draws  <- matrix(NA_real_, nrow=n_perm, ncol=length(STAT_IDS), dimnames=list(NULL, STAT_IDS))
  for (i in seq_len(n_perm)) {
    i1  <- br$base[first] + floor(stats::runif(n_loop) * n_pos) + 1
    idx <- c(i1, i1 + span)
    if (i == 1L) stopifnot(all(idx >= 1), all(idx <= length(bin_exon)),
                           all(idx[second] - idx[first] == span))
    draws[i, ] <- anchor_stats(bin_exon[idx], bin_gene[idx], bin_ccre[idx], n_loop)
  }
  draws
}

set.seed(SEED)
null_draws <- list()
null_rows  <- list()
for (key in names(obs_sets)) {
  s <- obs_sets[[key]]
  message("Permuting ", key, " (", N_PERM, " permutations)...")
  d <- simulate_null(s$av, N_PERM)
  null_draws[[key]] <- d
  # Calibration: one null draw treated as "observed" should not be extreme
  pseudo_p <- (1 + sum(d[-1, "E"] >= d[1, "E"])) / nrow(d)
  message("  calibration p (null draw vs rest, statistic E): ", signif(pseudo_p, 3))
  for (st in STAT_IDS) {
    obs <- unname(s$stats[st])
    nd  <- d[, st]
    null_rows[[length(null_rows) + 1]] <- data.frame(
      cond=s$cond, width_class=s$width_class, statistic=st, N_loops=s$av$n_loops,
      observed_pct=100 * obs,
      null_mean_pct=100 * mean(nd, na.rm=TRUE),
      null_lo_pct=100 * unname(quantile(nd, 0.025, na.rm=TRUE)),
      null_hi_pct=100 * unname(quantile(nd, 0.975, na.rm=TRUE)),
      fold_enrichment=obs / mean(nd, na.rm=TRUE),
      p_enriched=(1 + sum(nd >= obs - 1e-12, na.rm=TRUE)) / (N_PERM + 1),
      p_depleted=(1 + sum(nd <= obs + 1e-12, na.rm=TRUE)) / (N_PERM + 1),
      stringsAsFactors=FALSE)
  }
}
null_summary <- dplyr::bind_rows(null_rows)
null_summary$statistic_label <- STAT_LABELS[null_summary$statistic]
null_summary$p_direction     <- ifelse(null_summary$fold_enrichment >= 1,
                                       null_summary$p_enriched, null_summary$p_depleted)
message("\n── Item 2: observed vs background (", N_PERM, " permutations; minimum p = ",
        signif(1 / (N_PERM + 1), 3), ") ──────────────")
print(null_summary %>%
        dplyr::mutate(dplyr::across(c(observed_pct, null_mean_pct, null_lo_pct, null_hi_pct), ~round(.x, 2)),
                      fold_enrichment=round(fold_enrichment, 2),
                      dplyr::across(c(p_enriched, p_depleted, p_direction), ~signif(.x, 3))) %>%
        dplyr::select(-statistic_label))

# Direct answer: do anchors primarily (>50%) fall on functional elements?
message("\n── Functional elements (exon or cCRE), statistic E ─────────────────")
print(null_summary %>%
        dplyr::filter(statistic == "E") %>%
        dplyr::transmute(cond, width_class, N_loops, observed_pct=round(observed_pct, 1),
                         random_placement_pct=round(null_mean_pct, 1),
                         fold=round(fold_enrichment, 2), p_enriched=signif(p_enriched, 3),
                         majority_functional=observed_pct > 50))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Plots
# ══════════════════════════════════════════════════════════════════════════════
CLASS_COLORS <- c("Exon+cCRE"="#a50f15", "Exon"="#d6604d", "cCRE"="#f4a582",
                  "Intron"="#2166ac", "Intergenic"="#bdbdbd")
ANCHOR_TYPE_COLORS <- c("Exon-Exon"="#d6604d", "Exon-cCRE"="#984ea3", "Exon-Intron"="#4DAF4A",
                        "cCRE-cCRE"="#ff7f00", "Intron-Intron"="#2166ac")   # as in the Hi-C loop-strength figure
OBS_CI_COLOR <- "black"       # bootstrap CI of the observed value
NULL_COLOR   <- "#8c8c8c"     # random-placement mean and 95% range

# Plots are printed to the live device, then written as SVG to OUT_DIR.
# svglite writes plain SVG (no compositing groups, opens cleanly in Illustrator);
# base svg() is the fallback if svglite is not installed.
save_svg <- function(plot_obj, path, width, height) {
  use_svglite <- requireNamespace("svglite", quietly=TRUE)
  if (!use_svglite) message("svglite not installed; saving with grDevices::svg()")
  dev <- if (use_svglite) svglite::svglite else grDevices::svg
  ggplot2::ggsave(path, plot_obj, device=dev, width=width, height=height, units="in")
  stopifnot(file.exists(path), file.size(path) > 0)
  stopifnot(any(grepl("<svg", readLines(path, n=5, warn=FALSE))))
  message("Saved: ", path)
}

# Bootstrap 95% CI over loops for the share (%) of each class. cls is an n_loops x m integer
# matrix of class codes 1..nbins (m = 1: one class per loop; m = 2: the two anchors of a loop).
# Loops are resampled with replacement, so the two anchors of a loop stay together. Also
# returns the CI of the cumulative share (classes 1..k), the CI of each boundary in a stacked bar.
boot_class_ci <- function(cls, nbins, n_boot=N_BOOT) {
  cls <- as.matrix(cls)
  n   <- nrow(cls)
  stopifnot(n > 1, nbins > 1, all(cls >= 1L & cls <= nbins))
  pct  <- 100 * tabulate(cls, nbins=nbins) / length(cls)
  draw <- vapply(seq_len(n_boot), function(b)
    tabulate(cls[sample.int(n, n, replace=TRUE), ], nbins=nbins), numeric(nbins)) * 100 / length(cls)
  cdraw <- apply(draw, 2, cumsum)
  q     <- function(m) t(apply(m, 1, stats::quantile, probs=c(0.025, 0.975), names=FALSE))
  qd    <- q(draw); qc <- q(cdraw)
  data.frame(class=seq_len(nbins), pct=pct, lo=qd[, 1], hi=qd[, 2],
             cum=cumsum(pct), cum_lo=qc[, 1], cum_hi=qc[, 2])
}
# Wilson 95% interval (%) for k of n; used as an independent check on the bootstrap
wilson_pct <- function(k, n) {
  z <- stats::qnorm(0.975); p <- k / n
  ctr <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  100 * cbind(lo=ctr - hw, hi=ctr + hw)
}

# 5a. Anchor classes
plot1_df <-dplyr::bind_rows(lapply(CLASS_ORDER, function(k)
  anchor_summary %>% dplyr::transmute(cond, width_class, anchor_class=k, pct=.data[[paste0("pct_", k)]]))) %>%
  dplyr::mutate(cond=factor(cond, levels=COND_LEVELS),
                width_class=factor(width_class, levels=WIDTH_LEVELS),
                anchor_class=factor(anchor_class, levels=rev(CLASS_ORDER)))
n_df <- anchor_summary %>%
  dplyr::transmute(cond=factor(cond, levels=COND_LEVELS),
                   width_class=factor(width_class, levels=WIDTH_LEVELS),
                   label=paste0("N=", format(n_anchors, big.mark=",")))

# Bootstrap 95% CI of every boundary between stacked classes (classes stack bottom to top in
# CLASS_ORDER, so the boundary at the top of class k is the cumulative share of classes 1..k)
set.seed(SEED)
class_ci <- dplyr::bind_rows(lapply(obs_sets, function(s) {
  f  <- seq_len(s$av$n_loops)
  ci <- boot_class_ci(cbind(s$code[f], s$code[s$av$n_loops + f]), length(CLASS_ORDER))
  stopifnot(all(abs(ci$pct - 100 * tabulate(s$code, nbins=length(CLASS_ORDER)) / length(s$code)) < 1e-9))
  stopifnot(all(ci$cum_lo <= ci$cum + 1e-9), all(ci$cum_hi >= ci$cum - 1e-9))   # CI brackets the observed value
  data.frame(cond=s$cond, width_class=s$width_class, anchor_class=CLASS_ORDER[ci$class],
             pct=ci$pct, ci_lo=ci$lo, ci_hi=ci$hi, cum=ci$cum, cum_lo=ci$cum_lo, cum_hi=ci$cum_hi,
             stringsAsFactors=FALSE)
}))
rownames(class_ci) <- NULL
stopifnot(nrow(class_ci) == nrow(plot1_df))
chk <- plot1_df %>% dplyr::transmute(cond, width_class, anchor_class=as.character(anchor_class), pct_plot=pct) %>%
  dplyr::inner_join(class_ci %>% dplyr::transmute(cond=factor(cond, levels=COND_LEVELS),
                                                  width_class=factor(width_class, levels=WIDTH_LEVELS),
                                                  anchor_class, pct), by=c("cond", "width_class", "anchor_class"))
stopifnot(nrow(chk) == nrow(plot1_df), all(abs(chk$pct_plot - chk$pct) < 1e-9))   # same percentages as the bars
message("\nObserved anchor classes (%) with bootstrap 95% CI over loops:")
print(class_ci %>% dplyr::transmute(cond, width_class, anchor_class, pct=round(pct, 1),
                                    ci_lo=round(ci_lo, 1), ci_hi=round(ci_hi, 1)))
bnd_df <- class_ci %>%
  dplyr::filter(cum < 100 - 1e-9) %>%                                  # the top boundary is 100% by definition
  dplyr::mutate(cond=factor(cond, levels=COND_LEVELS), width_class=factor(width_class, levels=WIDTH_LEVELS))

fig_anchor_classes <- ggplot(plot1_df, aes(x=cond, y=pct, fill=anchor_class)) +
  geom_col(width=0.7, color="white", linewidth=0.3) +
  geom_errorbar(data=bnd_df, aes(x=cond, ymin=cum_lo, ymax=cum_hi), inherit.aes=FALSE,
                width=0.25, linewidth=0.5, color="black") +
  geom_text(data=n_df, aes(x=cond, y=102, label=label), inherit.aes=FALSE, size=3.5, vjust=0) +
  facet_wrap(~width_class) +
  scale_fill_manual(values=CLASS_COLORS, name=NULL, breaks=CLASS_ORDER) +
  scale_y_continuous(name="Loop anchors (%)", limits=c(0, 110), breaks=seq(0, 100, 25),
                     expand=expansion(mult=c(0, 0))) +
  labs(x=NULL, title="Micro-C loop anchors by genomic class",
       subtitle="Whiskers: bootstrap 95% CI (over loops) of each boundary between classes") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=14), plot.subtitle=element_text(size=9),
        axis.text.x=element_text(angle=30, hjust=1, vjust=1),
        strip.background=element_blank(), strip.text=element_text(face="bold"))
print(fig_anchor_classes)
save_svg(fig_anchor_classes, file.path(OUT_DIR, "microc_loop_anchor_classes.svg"), width=7.5, height=5)

# 5b. Permutation test: loops by anchor-type pair, observed vs random placement
#     Observed: bar with a bootstrap 95% CI over loops (black). Random placement: mean (grey
#     line across the bar) and 95% range of the permutations (grey whisker).
set.seed(SEED)
type_ci <- dplyr::bind_rows(lapply(obs_sets, function(s) {
  ci <- boot_class_ci(s$loop_type, length(LOOP_TYPES))
  n  <- s$av$n_loops
  stopifnot(all(abs(ci$pct - 100 * tabulate(s$loop_type, nbins=length(LOOP_TYPES)) / n) < 1e-9))
  stopifnot(all(ci$lo <= ci$pct + 1e-9), all(ci$hi >= ci$pct - 1e-9))   # CI brackets the observed value
  # independent check: for loop-level classes the bootstrap CI width matches the Wilson interval
  k  <- tabulate(s$loop_type, nbins=length(LOOP_TYPES))
  ok <- k >= 30 & n - k >= 30
  w  <- wilson_pct(k, n)
  stopifnot(all(abs(log((ci$hi - ci$lo)[ok] / (w[, "hi"] - w[, "lo"])[ok])) < log(1.3)))
  data.frame(cond=s$cond, width_class=s$width_class, statistic=LOOP_TYPES[ci$class],
             ci_lo_pct=ci$lo, ci_hi_pct=ci$hi, stringsAsFactors=FALSE)
}))
fig_df <- null_summary %>%
  dplyr::filter(statistic %in% names(ANCHOR_TYPE_COLORS)) %>%
  dplyr::inner_join(type_ci, by=c("cond", "width_class", "statistic")) %>%
  dplyr::mutate(cond=factor(cond, levels=COND_LEVELS),
                width_class=factor(width_class, levels=WIDTH_LEVELS),
                loop_type=factor(statistic, levels=names(ANCHOR_TYPE_COLORS)),
                fold_lo=ci_lo_pct / null_mean_pct, fold_hi=ci_hi_pct / null_mean_pct,
                p_label=ifelse(p_direction <= 1 / (N_PERM + 1) + 1e-12,
                               paste0("p<", signif(1 / (N_PERM + 1), 1)),
                               paste0("p=", signif(p_direction, 2))),
                top=pmax(ci_hi_pct, null_hi_pct))
stopifnot(nrow(fig_df) == length(unique(null_summary$cond)) * length(unique(null_summary$width_class)) *
            length(ANCHOR_TYPE_COLORS))
message("\nLoops by anchor-type pair: observed % with bootstrap 95% CI (over loops), random-placement mean and 95% range:")
print(fig_df %>% dplyr::transmute(cond, width_class, loop_type, observed_pct=round(observed_pct, 2),
                                  ci_lo=round(ci_lo_pct, 2), ci_hi=round(ci_hi_pct, 2),
                                  null_mean=round(null_mean_pct, 2), null_lo=round(null_lo_pct, 2),
                                  null_hi=round(null_hi_pct, 2), fold=round(fold_enrichment, 2),
                                  fold_ci_lo=round(fold_lo, 2), fold_ci_hi=round(fold_hi, 2)))

fig_loop_types <- ggplot(fig_df, aes(x=loop_type)) +
  geom_col(aes(y=observed_pct, fill=loop_type), width=0.6) +
  geom_errorbar(aes(ymin=null_mean_pct, ymax=null_mean_pct), width=0.6, linewidth=0.8, color=NULL_COLOR) +
  geom_errorbar(aes(ymin=null_lo_pct, ymax=null_hi_pct), width=0.32, linewidth=0.5, color=NULL_COLOR) +
  geom_errorbar(aes(ymin=ci_lo_pct, ymax=ci_hi_pct), width=0.16, linewidth=0.5, color=OBS_CI_COLOR) +
  geom_text(aes(y=top, label=paste0(round(fold_enrichment, 2), "x\n", p_label)),
            size=2.7, vjust=-0.2, lineheight=0.9) +
  facet_grid(width_class ~ cond, scales="free_y") +
  scale_fill_manual(values=ANCHOR_TYPE_COLORS, guide="none") +
  scale_y_continuous(name="Loops (%)", expand=expansion(mult=c(0, 0.2))) +
  labs(x=NULL, title="Loops by anchor-type pair vs random placement",
       subtitle="Bars with black whisker: observed and its bootstrap 95% CI (over loops).\nGrey line and whisker: random-placement mean and 95% range.\nLabels: fold (observed / random mean) and one-sided p, in the direction of the difference.\nOther (any intergenic anchor, or a cCRE-Intron pair) not shown") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=14), plot.subtitle=element_text(size=8),
        axis.text.x=element_text(angle=30, hjust=1, vjust=1),
        strip.background=element_blank(), strip.text=element_text(face="bold"))
print(fig_loop_types)
save_svg(fig_loop_types, file.path(OUT_DIR, "microc_loop_types_vs_random.svg"), width=9, height=7)

# 5c. Permutation distributions, small-anchor loops: exon / cCRE / functional anchors
SUMMARY_STATS <- c("A", "D", "E")
small_keys <- names(obs_sets)[vapply(obs_sets, function(s) s$width_class == WIDTH_LEVELS[1], logical(1))]
plot2_draws <- dplyr::bind_rows(lapply(small_keys, function(key) {
  s <- obs_sets[[key]]
  dplyr::bind_rows(lapply(SUMMARY_STATS, function(st)
    data.frame(cond=s$cond, statistic=st, value_pct=100 * null_draws[[key]][, st],
               stringsAsFactors=FALSE)))
})) %>% dplyr::mutate(cond=factor(cond, levels=COND_LEVELS),
                      statistic_label=factor(STAT_LABELS[statistic], levels=STAT_LABELS[SUMMARY_STATS]))
plot2_obs <- null_summary %>%
  dplyr::filter(width_class == WIDTH_LEVELS[1], statistic %in% SUMMARY_STATS) %>%
  dplyr::mutate(cond=factor(cond, levels=COND_LEVELS),
                statistic_label=factor(STAT_LABELS[statistic], levels=STAT_LABELS[SUMMARY_STATS]))

fig_background <- ggplot(plot2_draws, aes(x=value_pct)) +
  geom_histogram(bins=50, fill="gray70", color="white", linewidth=0.15) +
  geom_vline(data=plot2_obs, aes(xintercept=observed_pct, color=cond), linewidth=0.9) +
  facet_grid(statistic_label ~ cond, scales="free", labeller=labeller(statistic_label=label_wrap_gen(24))) +
  scale_color_manual(values=COND_COLORS, guide="none") +
  labs(x="Anchors overlapping element (%)", y="Permutations",
       title=paste0("Small-anchor loops (<=", SMALL_ANCHOR_MAX / 1000, " kb): observed vs random placement")) +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=14),
        strip.background=element_blank(), strip.text=element_text(face="bold"))
print(fig_background)
save_svg(fig_background, file.path(OUT_DIR, "microc_small_anchor_background.svg"), width=7.5, height=8)

# 5d. Wildtype loop length by anchor-type pair: observed vs random placement
#     Random placement moves each loop as a unit, so its length never changes; the null
#     length distribution of a category is the lengths of loops weighted by how often they
#     fall in that category across permutations.
N_LEN_PERM <- 1000L
LEN_TYPES  <- names(ANCHOR_TYPE_COLORS)          # Exon-Exon (top) ... Intron-Intron (bottom)
wt_keys    <- names(obs_sets)[vapply(obs_sets, function(s) s$cond == COND_LEVELS[1], logical(1))]

# n_loops x loop-type matrix: fraction of permutations in which each loop falls in each type
len_null_weights <- function(av, n_perm) {
  br     <- anchor_bin_rows(av)
  n_loop <- av$n_loops
  first  <- seq_len(n_loop)
  second <- n_loop + first
  k_obs  <- (av$start - 1L) %/% br$width
  span   <- k_obs[second] - k_obs[first]
  stopifnot(all(span >= 0), all(br$base[first] == br$base[second]),
            all(br$nb[first] == br$nb[second]), all(br$nb[first] > span))
  n_pos  <- br$nb[first] - span
  counts <- matrix(0L, nrow=n_loop, ncol=length(LOOP_TYPES))
  for (i in seq_len(n_perm)) {
    i1   <- br$base[first] + floor(stats::runif(n_loop) * n_pos) + 1
    idx  <- c(i1, i1 + span)
    code <- anchor_class_code(bin_exon[idx], bin_gene[idx], bin_ccre[idx])
    tp   <- loop_type_code(code[first], code[second])
    counts[cbind(first, tp)] <- counts[cbind(first, tp)] + 1L
  }
  counts / n_perm
}

set.seed(SEED)
len_list <- lapply(wt_keys, function(key) {
  s   <- obs_sets[[key]]
  n   <- s$av$n_loops
  len <- s$av$start[n + seq_len(n)] - s$av$start[seq_len(n)]          # anchor 1 start to anchor 2 start (bp)
  stopifnot(all(len == (s$loops$y1 - s$loops$x1)), all(len > 0))
  message("Length null for ", key, " (", N_LEN_PERM, " permutations)...")
  W <- len_null_weights(s$av, N_LEN_PERM)
  stopifnot(all(abs(rowSums(W) - 1) < 1e-9))
  # type fractions of this null agree with the main permutation test
  ref  <- null_summary %>% dplyr::filter(cond == s$cond, width_class == s$width_class, statistic %in% LOOP_TYPES)
  dmax <- max(abs(100 * colMeans(W) - ref$null_mean_pct[match(LOOP_TYPES, ref$statistic)]))
  message("  max |null % - main permutation null %|: ", round(dmax, 2))
  stopifnot(dmax < 1)
  list(len=len, type=s$loop_type, W=W)
})
len_all  <- unlist(lapply(len_list, `[[`, "len"))
type_all <- unlist(lapply(len_list, `[[`, "type"))
W_all    <- do.call(rbind, lapply(len_list, `[[`, "W"))
loglen   <- log10(len_all)
len_bw   <- stats::bw.nrd0(loglen)                                      # one bandwidth for every curve
len_grid <- seq(min(loglen) - 3 * len_bw, max(loglen) + 3 * len_bw, length.out=512)
message("Wildtype loops: ", length(len_all), "; length range ", min(len_all), "-", max(len_all),
        " bp; log10 bandwidth ", round(len_bw, 3))

ridge_density <- function(x, w) {
  stats::density(x, weights=w / sum(w), bw=len_bw, from=min(len_grid), to=max(len_grid),
                 n=length(len_grid))$y
}
w_quantile <- function(x, w, p) {
  o <- order(x); cw <- cumsum(w[o]) / sum(w)
  vapply(p, function(pp) x[o][which(cw >= pp)[1]], numeric(1))
}

ridge_rows <- lapply(seq_along(LEN_TYPES), function(r) {
  k   <- match(LEN_TYPES[r], LOOP_TYPES)
  sel <- type_all == k
  stopifnot(sum(sel) > 10)
  data.frame(loop_type=LEN_TYPES[r], rank=r, x=len_grid,
             obs=ridge_density(loglen[sel], rep(1, sum(sel))),
             null=ridge_density(loglen, W_all[, k]),
             n_obs=sum(sel), stringsAsFactors=FALSE)
})
ridge_df <- dplyr::bind_rows(ridge_rows)
ridge_scale <- 1.6 / max(c(ridge_df$obs, ridge_df$null))               # shared height scale, ridges overlap
ridge_df <- ridge_df %>%
  dplyr::mutate(base=length(LEN_TYPES) - rank,
                obs_top=base + obs * ridge_scale, null_top=base + null * ridge_scale,
                loop_type=factor(loop_type, levels=LEN_TYPES))
ridge_labels <- ridge_df %>% dplyr::distinct(loop_type, base, n_obs) %>% dplyr::arrange(dplyr::desc(base))

len_summary <- dplyr::bind_rows(lapply(LEN_TYPES, function(tp) {
  k   <- match(tp, LOOP_TYPES); sel <- type_all == k
  qo  <- w_quantile(len_all[sel], rep(1, sum(sel)), c(0.5, 0.9))
  qn  <- w_quantile(len_all, W_all[, k], c(0.5, 0.9))
  data.frame(loop_type=tp, N=sum(sel), median_kb_obs=qo[1] / 1e3, median_kb_null=qn[1] / 1e3,
             p90_kb_obs=qo[2] / 1e3, p90_kb_null=qn[2] / 1e3, stringsAsFactors=FALSE)
}))
message("\nWildtype loop length (kb): observed vs random placement")
print(len_summary %>% dplyr::mutate(dplyr::across(-c(loop_type, N), ~round(.x, 1))))

len_breaks <- c(1e4, 3e4, 1e5, 3e5, 1e6, 3e6)
fig_length <- ggplot(ridge_df, aes(x=x)) +
  geom_ribbon(aes(ymin=base, ymax=obs_top, fill=loop_type, group=loop_type),
              alpha=0.75, color="gray20", linewidth=0.3) +
  geom_line(aes(y=null_top, group=loop_type), color="black", linewidth=0.5, linetype="dashed") +
  scale_fill_manual(values=ANCHOR_TYPE_COLORS, guide="none") +
  scale_y_continuous(name=NULL, breaks=ridge_labels$base,
                     labels=paste0(ridge_labels$loop_type, "\n(N=", format(ridge_labels$n_obs, big.mark=",", trim=TRUE), ")"),
                     expand=expansion(mult=c(0.02, 0.05))) +
  scale_x_continuous(name="Loop length (anchor-to-anchor distance)", breaks=log10(len_breaks),
                     labels=c("10 kb", "30 kb", "100 kb", "300 kb", "1 Mb", "3 Mb")) +
  labs(title="Wildtype loop length by anchor-type pair",
       subtitle="Filled: observed. Dashed: random placement (loops moved as a unit).") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=14), plot.subtitle=element_text(size=9),
        axis.line.y=element_blank(), axis.ticks.y=element_blank())
print(fig_length)
save_svg(fig_length, file.path(OUT_DIR, "microc_wt_loop_length_by_anchor_type.svg"), width=7.5, height=6)

# 5e. Exon-Exon and Exon-cCRE loops: same gene, adjacent genes, or further apart (Wildtype)
#     Genes = merged gene-body loci (gene_gr), numbered in genomic order. An anchor belongs
#     to every locus it overlaps. Same gene: the two anchors share a locus. Adjacent genes:
#     no shared locus, but the loci are consecutive in genomic order. Background: the same
#     loop-as-a-unit permutations as above; in each permutation the relation is tabulated
#     among the loops that fall in the category.
N_GENE_PERM    <- 1000L
REL_LOOP_TYPES <- 1:2                                   # loop-type codes: Exon-Exon, Exon-cCRE
REL_LEVELS     <- c("Same gene", "Adjacent genes", "Non-adjacent genes", "cCRE anchor outside gene bodies",
                    "Same or adjacent gene")            # code 5 = codes 1 + 2
REL_ORDER      <- REL_LEVELS[c(1, 2, 5, 3, 4)]
wt_keys        <- names(obs_sets)[vapply(obs_sets, function(s) s$cond == COND_LEVELS[1], logical(1))]

# First and last gene locus overlapped by each interval (lo > hi: no locus)
make_locus_fn <- function(gr, chroms, chrom_lengths) {
  gr      <- GenomicRanges::reduce(gr, ignore.strand=TRUE)
  offsets <- setNames(c(0, cumsum(as.numeric(chrom_lengths)))[seq_along(chroms)], chroms)
  chr_gr  <- as.character(seqnames(gr))
  r_start <- unname(offsets[chr_gr]) + start(gr)
  r_end   <- unname(offsets[chr_gr]) + end(gr)
  o       <- order(r_start)
  r_start <- r_start[o]; r_end <- r_end[o]
  stopifnot(all(r_start[-1] > r_end[-length(r_end)]))
  function(chr, start, end) {
    g_start <- unname(offsets[chr]) + start
    g_end   <- unname(offsets[chr]) + end
    list(lo=findInterval(g_start - 1, r_end) + 1L, hi=findInterval(g_end, r_start))
  }
}
gene_loci <- make_locus_fn(gene_gr, STD_CHROMS, seqlens)

# 1 same gene, 2 adjacent genes, 3 non-adjacent genes, 4 an anchor has no gene body
gene_relation <- function(lo1, hi1, lo2, hi2) {
  v1   <- hi1 >= lo1
  v2   <- hi2 >= lo2
  same <- v1 & v2 & pmax(lo1, lo2) <= pmin(hi1, hi2)
  adj  <- v1 & v2 & !same & (hi1 == lo2 - 1L | hi2 == lo1 - 1L)
  out  <- rep(3L, length(lo1))
  out[adj]         <- 2L
  out[same]        <- 1L
  out[!(v1 & v2)]  <- 4L
  out
}

# Locus range of every grid bin (same grid as the overlap tables)
g_chr   <- rep(bin_rows$chr, bin_rows$n_bins)
g_width <- rep(as.numeric(bin_rows$width), bin_rows$n_bins)
g_k     <- sequence(bin_rows$n_bins) - 1
g_loc   <- gene_loci(g_chr, g_k * g_width + 1, (g_k + 1) * g_width)
bin_lo  <- g_loc$lo
bin_hi  <- g_loc$hi
stopifnot(length(bin_lo) == length(bin_gene), identical(bin_hi >= bin_lo, bin_gene))
rm(g_chr, g_width, g_k, g_loc)

gene_rel_null <- function(av, n_perm) {
  br     <- anchor_bin_rows(av)
  n_loop <- av$n_loops
  first  <- seq_len(n_loop)
  second <- n_loop + first
  k_obs  <- (av$start - 1L) %/% br$width
  span   <- k_obs[second] - k_obs[first]
  stopifnot(all(span >= 0), all(br$base[first] == br$base[second]),
            all(br$nb[first] == br$nb[second]), all(br$nb[first] > span))
  n_pos  <- br$nb[first] - span
  cnt    <- array(0L, dim=c(n_perm, length(REL_LOOP_TYPES), 4L))
  frac   <- matrix(0, nrow=n_perm, ncol=length(LOOP_TYPES))
  for (i in seq_len(n_perm)) {
    i1   <- br$base[first] + floor(stats::runif(n_loop) * n_pos) + 1
    i2   <- i1 + span
    idx  <- c(i1, i2)
    code <- anchor_class_code(bin_exon[idx], bin_gene[idx], bin_ccre[idx])
    tp   <- loop_type_code(code[first], code[second])
    rel  <- gene_relation(bin_lo[i1], bin_hi[i1], bin_lo[i2], bin_hi[i2])
    for (j in seq_along(REL_LOOP_TYPES)) cnt[i, j, ] <- tabulate(rel[tp == REL_LOOP_TYPES[j]], nbins=4L)
    frac[i, ] <- tabulate(tp, nbins=length(LOOP_TYPES)) / n_loop
  }
  list(cnt=cnt, frac=frac)
}

# rel_by_type: list (one element per REL_LOOP_TYPES) of the per-loop relation codes 1..4 of the
# observed loops; gives the bootstrap CI over loops of each relation's share. Relation 5
# (same or adjacent gene) takes the CI of the cumulative share of relations 1 and 2.
rel_summary <- function(cond, width_class, obs_cnt, cnt, rel_by_type) {
  dplyr::bind_rows(lapply(seq_along(REL_LOOP_TYPES), function(j) {
    stopifnot(all(tabulate(rel_by_type[[j]], nbins=4L) == obs_cnt[j, ]))
    ci  <- boot_class_ci(rel_by_type[[j]], 4L)
    oc  <- c(obs_cnt[j, ], obs_cnt[j, 1] + obs_cnt[j, 2])
    nc  <- cbind(cnt[, j, ], cnt[, j, 1] + cnt[, j, 2])
    opc <- oc / sum(obs_cnt[j, ])
    npc <- nc / rowSums(cnt[, j, ])
    ci_lo  <- c(ci$lo, ci$cum_lo[2]);  ci_hi  <- c(ci$hi, ci$cum_hi[2])
    cum    <- c(ci$cum, NA);           cum_lo <- c(ci$cum_lo[1:4], NA); cum_hi <- c(ci$cum_hi[1:4], NA)
    dplyr::bind_rows(lapply(1:5, function(r) {
      nd <- npc[, r]; ob <- opc[r]
      data.frame(cond=cond, width_class=width_class, loop_type=LOOP_TYPES[REL_LOOP_TYPES[j]],
                 relation=REL_LEVELS[r], N_loops=sum(obs_cnt[j, ]), observed_n=oc[r],
                 observed_pct=100 * ob, ci_lo_pct=ci_lo[r], ci_hi_pct=ci_hi[r],
                 cum_pct=cum[r], cum_lo_pct=cum_lo[r], cum_hi_pct=cum_hi[r],
                 null_mean_pct=100 * mean(nd, na.rm=TRUE),
                 null_lo_pct=100 * unname(quantile(nd, 0.025, na.rm=TRUE)),
                 null_hi_pct=100 * unname(quantile(nd, 0.975, na.rm=TRUE)),
                 fold=ob / mean(nd, na.rm=TRUE),
                 p_enriched=(1 + sum(nd >= ob - 1e-12, na.rm=TRUE)) / (N_GENE_PERM + 1),
                 p_depleted=(1 + sum(nd <= ob + 1e-12, na.rm=TRUE)) / (N_GENE_PERM + 1),
                 stringsAsFactors=FALSE)
    }))
  }))
}

set.seed(SEED)
rel_res <- lapply(wt_keys, function(key) {
  s  <- obs_sets[[key]]
  n  <- s$av$n_loops
  f  <- seq_len(n)
  # observed relation from a direct locus lookup; agrees with the grid-table lookup
  loc <- gene_loci(s$av$chr, s$av$start, s$av$end)
  br  <- anchor_bin_rows(s$av)
  k_o <- (s$av$start - 1L) %/% br$width
  stopifnot(identical(loc$lo, bin_lo[br$base + k_o + 1]), identical(loc$hi, bin_hi[br$base + k_o + 1]))
  gr  <- GRanges(s$av$chr, IRanges(s$av$start, s$av$end))
  stopifnot(identical(pmax(loc$hi - loc$lo + 1L, 0L),
                      as.integer(suppressWarnings(countOverlaps(gr, gene_gr, ignore.strand=TRUE)))))
  rel <- gene_relation(loc$lo[f], loc$hi[f], loc$lo[n + f], loc$hi[n + f])
  obs_cnt <- t(vapply(REL_LOOP_TYPES, function(k) tabulate(rel[s$loop_type == k], nbins=4L), integer(4)))
  stopifnot(all(rowSums(obs_cnt) == tabulate(s$loop_type, nbins=length(LOOP_TYPES))[REL_LOOP_TYPES]))
  stopifnot(obs_cnt[1, 4] == 0)                                  # exon anchors always lie in a gene body
  message("Gene-relation null for ", key, " (", N_GENE_PERM, " permutations)...")
  nl <- gene_rel_null(s$av, N_GENE_PERM)
  stopifnot(all(nl$cnt[, 1, 4] == 0))
  # type fractions of this null agree with the main permutation test
  ref  <- null_summary %>% dplyr::filter(cond == s$cond, width_class == s$width_class, statistic %in% LOOP_TYPES)
  dmax <- max(abs(100 * colMeans(nl$frac) - ref$null_mean_pct[match(LOOP_TYPES, ref$statistic)]))
  message("  max |null % - main permutation null %|: ", round(dmax, 2))
  stopifnot(dmax < 1)
  list(cond=s$cond, width_class=s$width_class, obs_cnt=obs_cnt, cnt=nl$cnt,
       rel_by_type=lapply(REL_LOOP_TYPES, function(k) rel[s$loop_type == k]))
})

pool_obs <- Reduce(`+`, lapply(rel_res, `[[`, "obs_cnt"))
pool_cnt <- Reduce(`+`, lapply(rel_res, `[[`, "cnt"))
pool_rel <- lapply(seq_along(REL_LOOP_TYPES), function(j) unlist(lapply(rel_res, function(r) r$rel_by_type[[j]])))
set.seed(SEED)
gene_rel_summary <- dplyr::bind_rows(c(
  lapply(rel_res, function(r) rel_summary(r$cond, r$width_class, r$obs_cnt, r$cnt, r$rel_by_type)),
  list(rel_summary(COND_LEVELS[1], "All anchor widths", pool_obs, pool_cnt, pool_rel)))) %>%
  dplyr::filter(!(loop_type == "Exon-Exon" & relation == REL_LEVELS[4])) %>%
  dplyr::mutate(relation=factor(relation, levels=REL_ORDER),
                p_direction=ifelse(fold >= 1, p_enriched, p_depleted)) %>%
  dplyr::arrange(width_class != "All anchor widths", width_class, loop_type, relation)
message("\n── Wildtype Exon-Exon and Exon-cCRE loops: same gene / adjacent genes (",
        N_GENE_PERM, " permutations; minimum p = ", signif(1 / (N_GENE_PERM + 1), 3), ") ──")
print(as.data.frame(gene_rel_summary %>% dplyr::select(-cum_pct, -cum_lo_pct, -cum_hi_pct) %>%
                      dplyr::mutate(dplyr::across(c(observed_pct, ci_lo_pct, ci_hi_pct, null_mean_pct, null_lo_pct, null_hi_pct), ~round(.x, 1)),
                                    fold=round(fold, 2),
                                    dplyr::across(c(p_enriched, p_depleted, p_direction), ~signif(.x, 3)))))

# Plot 1: observed relation of the two anchors' genes (Wildtype)
REL_COLORS <- c("Same gene"="#08306b", "Adjacent genes"="#4292c6",
                "Non-adjacent genes"="#c6dbef", "cCRE anchor outside gene bodies"="#bdbdbd")
stk_df <- gene_rel_summary %>%
  dplyr::filter(relation != REL_LEVELS[5], width_class != "All anchor widths") %>%
  dplyr::mutate(relation=factor(as.character(relation), levels=rev(REL_LEVELS[1:4])),   # first level stacks on top
                width_class=factor(width_class, levels=WIDTH_LEVELS),
                loop_type=factor(loop_type, levels=LOOP_TYPES[REL_LOOP_TYPES]),
                x_label=paste0(loop_type, "\n(N=", format(N_loops, big.mark=",", trim=TRUE), ")"),
                text_col=ifelse(relation %in% REL_LEVELS[1:2], "white", "black")) %>%
  dplyr::arrange(width_class, loop_type)
stk_df$x_label <- factor(stk_df$x_label, levels=unique(stk_df$x_label))
stopifnot(!anyNA(stk_df$width_class), !anyNA(stk_df$loop_type))
chk <- stk_df %>% dplyr::group_by(width_class, loop_type) %>% dplyr::summarise(obs=sum(observed_pct), .groups="drop")
stopifnot(all(abs(chk$obs - 100) < 1e-6))                         # classes of each loop type sum to 100%
# Whiskers: bootstrap 95% CI of each boundary between stacked classes (classes stack bottom to
# top as relation codes 1..4, so the boundary at the top of class r is the cumulative share of 1..r)
chk <- stk_df %>% dplyr::mutate(code=match(as.character(relation), REL_LEVELS[1:4])) %>%
  dplyr::arrange(width_class, loop_type, code) %>% dplyr::group_by(width_class, loop_type) %>%
  dplyr::mutate(cum_obs=cumsum(observed_pct)) %>% dplyr::ungroup()
stopifnot(all(abs(chk$cum_obs - chk$cum_pct) < 1e-6))              # boundary CIs are for the plotted boundaries
stopifnot(all(chk$cum_lo_pct <= chk$cum_pct + 1e-9), all(chk$cum_hi_pct >= chk$cum_pct - 1e-9))
bnd_gene_rel <- stk_df %>% dplyr::filter(cum_pct < 100 - 1e-9)    # the top boundary is 100% by definition

fig_gene_rel <- ggplot(stk_df, aes(x=x_label, y=observed_pct, fill=relation)) +
  geom_col(width=0.7, color="white", linewidth=0.3) +
  geom_errorbar(data=bnd_gene_rel, aes(x=x_label, ymin=cum_lo_pct, ymax=cum_hi_pct), inherit.aes=FALSE,
                width=0.25, linewidth=0.5, color="black") +
  geom_text(aes(label=ifelse(observed_pct >= 3, paste0(round(observed_pct), "%"), ""), color=text_col),
            position=position_stack(vjust=0.5), size=3.2) +
  scale_color_identity() +
  facet_wrap(~width_class, nrow=1, scales="free_x") +
  scale_fill_manual(values=REL_COLORS, name=NULL, guide=guide_legend(reverse=TRUE)) +
  scale_y_continuous(name="Loops (%)", expand=expansion(mult=c(0, 0.02))) +
  labs(x=NULL, title="Wildtype loops: gene relation of the two anchors",
       subtitle="Genes = merged gene-body loci; adjacent = consecutive loci in genomic order.\nAn Exon-cCRE loop has its cCRE anchor outside gene bodies when that anchor overlaps no gene body.\nWhiskers: bootstrap 95% CI (over loops) of each boundary between classes.") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13), plot.subtitle=element_text(size=8),
        strip.background=element_blank(), strip.text=element_text(face="bold"),
        legend.position="right")
print(fig_gene_rel)
save_svg(fig_gene_rel, file.path(OUT_DIR, "microc_wt_exon_loops_gene_relation.svg"), width=8.5, height=5)

# Plot 2: all Wildtype loops (every loop type): same / adjacent / non-adjacent gene vs random
#         placement of each loop as a unit (same length), no conditioning on loop type.
#         Loops with an anchor outside gene bodies are split by whether that anchor is on a cCRE,
#         so the "no cCRE" class holds no loop with a regulatory-feature anchor.
ALL_REL_LABELS <- c(REL_LEVELS[1:3], "An anchor outside gene bodies, none on a cCRE",
                    "An anchor on a cCRE outside gene bodies")
N_ALL_REL      <- length(ALL_REL_LABELS)
# relation 4 (an anchor has no gene body) -> 4 if no outside anchor is on a cCRE, 5 if one is
split_outside_ccre <- function(rel, in_gene1, in_gene2, ccre1, ccre2) {
  rel[rel == 4L & ((!in_gene1 & ccre1) | (!in_gene2 & ccre2))] <- 5L
  rel
}
gene_rel_null_all <- function(av, n_perm) {
  br     <- anchor_bin_rows(av)
  n_loop <- av$n_loops
  first  <- seq_len(n_loop)
  second <- n_loop + first
  k_obs  <- (av$start - 1L) %/% br$width
  span   <- k_obs[second] - k_obs[first]
  stopifnot(all(span >= 0), all(br$base[first] == br$base[second]),
            all(br$nb[first] == br$nb[second]), all(br$nb[first] > span))
  n_pos  <- br$nb[first] - span
  cnt    <- matrix(0L, nrow=n_perm, ncol=N_ALL_REL)
  for (i in seq_len(n_perm)) {
    i1  <- br$base[first] + floor(stats::runif(n_loop) * n_pos) + 1
    i2  <- i1 + span
    rel <- gene_relation(bin_lo[i1], bin_hi[i1], bin_lo[i2], bin_hi[i2])
    rel <- split_outside_ccre(rel, bin_gene[i1], bin_gene[i2], bin_ccre[i1], bin_ccre[i2])
    cnt[i, ] <- tabulate(rel, nbins=N_ALL_REL)
  }
  cnt
}
all_rel_summary <- function(cond, width_class, obs_cnt, cnt, rel) {
  n <- sum(obs_cnt); opc <- obs_cnt / n; npc <- cnt / n
  stopifnot(all(rowSums(cnt) == n), all(tabulate(rel, nbins=N_ALL_REL) == obs_cnt))
  ci <- boot_class_ci(rel, N_ALL_REL)                  # bootstrap 95% CI over loops
  dplyr::bind_rows(lapply(seq_len(N_ALL_REL), function(r) {
    nd <- npc[, r]; ob <- opc[r]
    data.frame(cond=cond, width_class=width_class, relation=ALL_REL_LABELS[r], N_loops=n,
               observed_n=obs_cnt[r], observed_pct=100 * ob, ci_lo_pct=ci$lo[r], ci_hi_pct=ci$hi[r],
               null_mean_pct=100 * mean(nd),
               null_lo_pct=100 * unname(quantile(nd, 0.025)), null_hi_pct=100 * unname(quantile(nd, 0.975)),
               fold=ob / mean(nd),
               p_enriched=(1 + sum(nd >= ob - 1e-12)) / (N_GENE_PERM + 1),
               p_depleted=(1 + sum(nd <= ob + 1e-12)) / (N_GENE_PERM + 1),
               stringsAsFactors=FALSE)
  }))
}

set.seed(SEED)
all_res <- lapply(seq_along(wt_keys), function(i) {
  key <- wt_keys[i]
  s   <- obs_sets[[key]]
  n   <- s$av$n_loops
  f   <- seq_len(n)
  loc <- gene_loci(s$av$chr, s$av$start, s$av$end)
  rel <- gene_relation(loc$lo[f], loc$hi[f], loc$lo[n + f], loc$hi[n + f])
  rel <- split_outside_ccre(rel, loc$hi[f] >= loc$lo[f], loc$hi[n + f] >= loc$lo[n + f],
                            s$hit$ccre[f], s$hit$ccre[n + f])
  obs_cnt <- tabulate(rel, nbins=N_ALL_REL)
  stopifnot(sum(obs_cnt) == n)
  # agrees with the Exon-Exon / Exon-cCRE counts above (classes 4 + 5 together = relation 4 there)
  for (k in REL_LOOP_TYPES) {
    tb <- tabulate(rel[s$loop_type == k], nbins=N_ALL_REL)
    stopifnot(all(c(tb[1:3], tb[4] + tb[5]) == rel_res[[i]]$obs_cnt[k, ]))
  }
  # the cCRE anchor of an Exon-cCRE loop outside gene bodies is on a cCRE by definition
  stopifnot(tabulate(rel[s$loop_type == 2], nbins=N_ALL_REL)[4] == 0)
  message("All-loop gene-relation null for ", key, " (", N_GENE_PERM, " permutations)...")
  list(cond=s$cond, width_class=s$width_class, obs_cnt=obs_cnt, rel=rel,
       cnt=gene_rel_null_all(s$av, N_GENE_PERM))
})
set.seed(SEED)
all_gene_rel_summary <- dplyr::bind_rows(c(
  lapply(all_res, function(r) all_rel_summary(r$cond, r$width_class, r$obs_cnt, r$cnt, r$rel)),
  list(all_rel_summary(COND_LEVELS[1], "All anchor widths",
                       Reduce(`+`, lapply(all_res, `[[`, "obs_cnt")),
                       Reduce(`+`, lapply(all_res, `[[`, "cnt")),
                       unlist(lapply(all_res, `[[`, "rel"))))))
all_gene_rel_summary$p_direction <- ifelse(all_gene_rel_summary$fold >= 1,
                                           all_gene_rel_summary$p_enriched, all_gene_rel_summary$p_depleted)
message("\n── All Wildtype loops: same gene / adjacent genes vs random placement (",
        N_GENE_PERM, " permutations; minimum p = ", signif(1 / (N_GENE_PERM + 1), 3), ") ──")
print(all_gene_rel_summary %>%
        dplyr::mutate(dplyr::across(c(observed_pct, ci_lo_pct, ci_hi_pct, null_mean_pct, null_lo_pct, null_hi_pct), ~round(.x, 1)),
                      fold=round(fold, 2),
                      dplyr::across(c(p_enriched, p_depleted, p_direction), ~signif(.x, 3))))

all_plot_df <- all_gene_rel_summary %>%
  dplyr::filter(width_class != "All anchor widths") %>%
  dplyr::mutate(relation=factor(relation, levels=ALL_REL_LABELS),
                width_class=factor(width_class, levels=WIDTH_LEVELS),
                p_label=ifelse(p_direction <= 1 / (N_GENE_PERM + 1) + 1e-12,
                               paste0("p<", signif(1 / (N_GENE_PERM + 1), 1)),
                               paste0("p=", signif(p_direction, 2))),
                fold_label=ifelse(is.finite(fold), paste0(round(fold, 2), "x"), "NA"),
                top=pmax(ci_hi_pct, null_hi_pct))
stopifnot(!anyNA(all_plot_df$width_class), !anyNA(all_plot_df$relation))
stopifnot(all(all_plot_df$ci_lo_pct <= all_plot_df$observed_pct + 1e-9),
          all(all_plot_df$ci_hi_pct >= all_plot_df$observed_pct - 1e-9))   # CI brackets the observed value
chk <- all_plot_df %>% dplyr::group_by(width_class) %>% dplyr::summarise(obs=sum(observed_pct), nul=sum(null_mean_pct))
stopifnot(all(abs(chk$obs - 100) < 1e-6), all(abs(chk$nul - 100) < 1e-6))   # classes sum to 100%

fig_all_rel <- ggplot(all_plot_df, aes(x=relation)) +
  geom_col(aes(y=observed_pct), fill="#08519c", width=0.7) +
  geom_errorbar(aes(ymin=null_mean_pct, ymax=null_mean_pct), width=0.7, linewidth=0.8, color=NULL_COLOR) +
  geom_errorbar(aes(ymin=null_lo_pct, ymax=null_hi_pct), width=0.32, linewidth=0.5, color=NULL_COLOR) +
  geom_errorbar(aes(ymin=ci_lo_pct, ymax=ci_hi_pct), width=0.16, linewidth=0.5, color=OBS_CI_COLOR) +
  geom_text(aes(y=top, label=paste0(fold_label, "\n", p_label)), size=2.7, vjust=-0.2, lineheight=0.9) +
  facet_wrap(~width_class, nrow=1) +
  scale_x_discrete(labels=c("Same gene"="Same\ngene", "Adjacent genes"="Adjacent\ngenes",
                            "Non-adjacent genes"="Non-adjacent\ngenes",
                            "An anchor outside gene bodies, none on a cCRE"="Outside genes,\nno cCRE",
                            "An anchor on a cCRE outside gene bodies"="cCRE outside\ngenes")) +
  scale_y_continuous(name="Loops (%)", expand=expansion(mult=c(0, 0.25))) +
  labs(x=NULL, title="All Wildtype loops: gene relation of the two anchors\nvs random placement",
       subtitle="Bars with black whisker: observed and its bootstrap 95% CI (over loops).\nGrey line and whisker: mean and 95% range of random placement\n(each loop moved to a random position on its chromosome, same length).\nLabels: fold (observed / random mean) and one-sided p, in the direction of the difference.") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13), plot.subtitle=element_text(size=8),
        axis.text.x=element_text(size=8.5),
        strip.background=element_blank(), strip.text=element_text(face="bold"))
print(fig_all_rel)
save_svg(fig_all_rel, file.path(OUT_DIR, "microc_wt_all_loops_gene_relation_vs_random.svg"), width=11, height=5.5)

# 5f. Wildtype loop length by gene relation (same / adjacent / non-adjacent genes), Exon-Exon and
#     Exon-cCRE loops, all anchor widths pooled: observed vs random placement. Each loop is
#     moved as a unit, so its length never changes; the null length distribution of a class is
#     the lengths of all loops weighted by how often each falls in that class across permutations.
REL_RIDGE <- 1:3                                         # same gene (top), adjacent, non-adjacent (bottom)

# n_loops x (loop type x gene relation) matrix: fraction of permutations in which each loop falls in
# each class; column = (loop type - 1) * length(REL_RIDGE) + relation
gene_rel_len_weights <- function(av, n_perm) {
  br     <- anchor_bin_rows(av)
  n_loop <- av$n_loops
  first  <- seq_len(n_loop)
  second <- n_loop + first
  k_obs  <- (av$start - 1L) %/% br$width
  span   <- k_obs[second] - k_obs[first]
  stopifnot(all(span >= 0), all(br$base[first] == br$base[second]),
            all(br$nb[first] == br$nb[second]), all(br$nb[first] > span))
  n_pos  <- br$nb[first] - span
  W      <- matrix(0L, nrow=n_loop, ncol=length(REL_LOOP_TYPES) * length(REL_RIDGE))
  for (i in seq_len(n_perm)) {
    i1   <- br$base[first] + floor(stats::runif(n_loop) * n_pos) + 1
    i2   <- i1 + span
    idx  <- c(i1, i2)
    code <- anchor_class_code(bin_exon[idx], bin_gene[idx], bin_ccre[idx])
    tp   <- loop_type_code(code[first], code[second])
    rel  <- gene_relation(bin_lo[i1], bin_hi[i1], bin_lo[i2], bin_hi[i2])
    ok   <- tp %in% REL_LOOP_TYPES & rel %in% REL_RIDGE
    cm   <- cbind(first[ok], (tp[ok] - 1L) * length(REL_RIDGE) + rel[ok])
    W[cm] <- W[cm] + 1L
  }
  W / n_perm
}

set.seed(SEED)
rel_len_list <- lapply(seq_along(wt_keys), function(i) {
  key <- wt_keys[i]
  s   <- obs_sets[[key]]
  n   <- s$av$n_loops
  f   <- seq_len(n)
  len <- s$av$start[n + f] - s$av$start[f]                             # anchor 1 start to anchor 2 start (bp)
  stopifnot(all(len == (s$loops$y1 - s$loops$x1)), all(len > 0))
  loc <- gene_loci(s$av$chr, s$av$start, s$av$end)
  rel <- gene_relation(loc$lo[f], loc$hi[f], loc$lo[n + f], loc$hi[n + f])
  message("Length null for ", key, " (", N_LEN_PERM, " permutations)...")
  W   <- gene_rel_len_weights(s$av, N_LEN_PERM)
  # expected loops per class agree with the gene-relation permutation test above
  exp_cnt <- matrix(colSums(W), nrow=length(REL_LOOP_TYPES), byrow=TRUE)
  ref_cnt <- apply(rel_res[[i]]$cnt[, , REL_RIDGE, drop=FALSE], c(2, 3), mean)
  dmax    <- max(abs(exp_cnt - ref_cnt)) / n * 100
  message("  max |expected loops - gene-relation null| (% of loops): ", round(dmax, 2))
  stopifnot(dmax < 1)
  list(len=len, type=s$loop_type, rel=rel, W=W)
})
rl_len  <- unlist(lapply(rel_len_list, `[[`, "len"))
rl_type <- unlist(lapply(rel_len_list, `[[`, "type"))
rl_rel  <- unlist(lapply(rel_len_list, `[[`, "rel"))
rl_W    <- do.call(rbind, lapply(rel_len_list, `[[`, "W"))
rl_log  <- log10(rl_len)
rl_bw   <- stats::bw.nrd0(rl_log)                                      # one bandwidth for every curve
rl_grid <- seq(min(rl_log) - 3 * rl_bw, max(rl_log) + 3 * rl_bw, length.out=512)
message("Wildtype loops: ", length(rl_len), "; log10 bandwidth ", round(rl_bw, 3))

rl_density <- function(x, w) {
  stats::density(x, weights=w / sum(w), bw=rl_bw, from=min(rl_grid), to=max(rl_grid),
                 n=length(rl_grid))$y
}
rl_quantile <- function(x, w, p) {
  o <- order(x); cw <- cumsum(w[o]) / sum(w)
  vapply(p, function(pp) x[o][which(cw >= pp)[1]], numeric(1))
}

rel_len_df <- dplyr::bind_rows(lapply(REL_LOOP_TYPES, function(j) {
  dplyr::bind_rows(lapply(REL_RIDGE, function(r) {
    sel <- rl_type == j & rl_rel == r
    stopifnot(sum(sel) > 1)
    col <- (j - 1L) * length(REL_RIDGE) + r
    data.frame(loop_type=LOOP_TYPES[j], relation=REL_LEVELS[r], rank=r, x=rl_grid,
               obs=rl_density(rl_log[sel], rep(1, sum(sel))),
               null=rl_density(rl_log, rl_W[, col]),
               n_obs=sum(sel), stringsAsFactors=FALSE)
  }))
}))
rl_scale <- 1.6 / max(c(rel_len_df$obs, rel_len_df$null))              # shared height scale, ridges overlap
rel_len_df <- rel_len_df %>%
  dplyr::mutate(base=length(REL_RIDGE) - rank,
                obs_top=base + obs * rl_scale, null_top=base + null * rl_scale,
                loop_type=factor(loop_type, levels=LOOP_TYPES[REL_LOOP_TYPES]))
rel_len_n <- rel_len_df %>% dplyr::distinct(loop_type, relation, base, n_obs)

rel_len_summary <- dplyr::bind_rows(lapply(REL_LOOP_TYPES, function(j) {
  dplyr::bind_rows(lapply(REL_RIDGE, function(r) {
    sel <- rl_type == j & rl_rel == r
    col <- (j - 1L) * length(REL_RIDGE) + r
    qo  <- rl_quantile(rl_len[sel], rep(1, sum(sel)), c(0.5, 0.9))
    qn  <- rl_quantile(rl_len, rl_W[, col], c(0.5, 0.9))
    data.frame(loop_type=LOOP_TYPES[j], relation=REL_LEVELS[r], N=sum(sel),
               median_kb_obs=qo[1] / 1e3, median_kb_null=qn[1] / 1e3,
               p90_kb_obs=qo[2] / 1e3, p90_kb_null=qn[2] / 1e3, stringsAsFactors=FALSE)
  }))
}))
message("\nWildtype loop length (kb) by gene relation: observed vs random placement")
print(as.data.frame(rel_len_summary %>% dplyr::mutate(dplyr::across(-c(loop_type, relation, N), ~round(.x, 1)))))

rl_breaks <- c(1e4, 1e5, 1e6)
fig_gene_rel_len <- ggplot(rel_len_df, aes(x=x)) +
  geom_ribbon(aes(ymin=base, ymax=obs_top, fill=loop_type, group=relation),
              alpha=0.75, color="gray20", linewidth=0.3) +
  geom_line(aes(y=null_top, group=relation), color="black", linewidth=0.5, linetype="dashed") +
  geom_text(data=rel_len_n, aes(x=max(rl_grid), y=base + 0.15,
                                label=paste0("N=", format(n_obs, big.mark=",", trim=TRUE))),
            hjust=1, size=3, color="gray20", inherit.aes=FALSE) +
  facet_wrap(~loop_type, nrow=1) +
  scale_fill_manual(values=ANCHOR_TYPE_COLORS, guide="none") +
  scale_y_continuous(name=NULL, breaks=length(REL_RIDGE) - REL_RIDGE, labels=REL_LEVELS[REL_RIDGE],
                     expand=expansion(mult=c(0.02, 0.05))) +
  scale_x_continuous(name="Loop length (anchor-to-anchor distance)", breaks=log10(rl_breaks),
                     labels=c("10 kb", "100 kb", "1 Mb")) +
  labs(title="Wildtype loop length by gene relation of the two anchors",
       subtitle="Filled: observed. Dashed: random placement (loops moved as a unit). All anchor widths pooled.") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13), plot.subtitle=element_text(size=8),
        axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background=element_blank(), strip.text=element_text(face="bold"))
print(fig_gene_rel_len)
save_svg(fig_gene_rel_len, file.path(OUT_DIR, "microc_wt_exon_loops_gene_relation_length.svg"), width=10, height=5)