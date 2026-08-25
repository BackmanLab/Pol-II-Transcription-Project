##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# POLR2A ChIA-PET loop strength by anchor type
# (Exon-Exon / Exon-cCRE / Exon-Intron / cCRE-cCRE / Intron-Intron).
#
# Data: public POLR2A ChIA-PET, HCT116 (wildtype, no degron), untreated,
# ENCODE accession ENCSR035PVZ (see Table 5).
#
# Classification/plotting logic (classify_loop_anchors,
# plot_loop_strength_raincloud) is reused from
# Hic_loop_analysis/hic_general_analysis.R, where it already existed
# correctly titled/labeled for ChIA-PET but was never wired up to
# real ChIA-PET input in that file. The ChIA-PET loop loader below
# (load_chiapet_bed12) is new -- it was not present anywhere in the
# original codebase and has not been tested against a real downloaded
# ENCSR035PVZ file. The parsing logic should be checked against the
# actual file before the output is trusted; see the loader's own
# comments for the expected format.
##################################################################

library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(ggdist)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ══════════════════════════════════════════════════════════════════════════════
# 0. Paths -- fill in the actual downloaded file location
# ══════════════════════════════════════════════════════════════════════════════

CHIAPET_BED12 <- "/path/to/ENCSR035PVZ.bed"   # <- CHECK: fill in actual path
ENCODE_CCRE   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/encode/wt_hct116.bed"

# ── 1. Load ChIA-PET loops -- new, unverified against a real file ──────────
#
# ENCODE ChIA-PET interaction calls are distributed as BED12, where the two
# interacting anchors are the two "blocks," and the `name` field packs both
# anchor coordinates plus the PET count as one string, e.g.:
#   chr1  875400  877590  chr1:875400..877590-chr8:126185955..126188744,2  200  .  ...
# i.e. name = "<anchor1_coords>-<anchor2_coords>,<PET_count>"
#
# The name-field format above is the assumed layout and has not been
# confirmed against the actual downloaded ENCSR035PVZ file.

load_chiapet_bed12 <- function(path) {
  raw <- read.table(path, header=FALSE, sep="\t", stringsAsFactors=FALSE)

  # Parse the name field: chrA:startA..endA-chrB:startB..endB,PETcount
  name_pattern <- "^([^:]+):(\\d+)\\.\\.(\\d+)-([^:]+):(\\d+)\\.\\.(\\d+),(\\d+)$"
  m <- regmatches(raw$V4, regexec(name_pattern, raw$V4))

  parsed <- do.call(rbind, lapply(m, function(x) {
    if (length(x) != 8) return(rep(NA, 7))
    x[2:8]
  }))
  parsed <- as.data.frame(parsed, stringsAsFactors=FALSE)
  colnames(parsed) <- c("chr1","start1","end1","chr2","start2","end2","pet_count")

  out <- parsed %>%
    dplyr::mutate(
      start1=as.integer(start1), end1=as.integer(end1),
      start2=as.integer(start2), end2=as.integer(end2),
      pet_count=as.integer(pet_count),
      centroid1=as.integer((start1+end1)/2), centroid2=as.integer((start2+end2)/2)
    ) %>%
    filter(!is.na(start1), !is.na(start2))

  message("Loaded ", nrow(out), " ChIA-PET loops from ", basename(path),
          " (", sum(is.na(parsed$chr1)), " rows failed name-field parsing -- ",
          "if this is large relative to the total, the parser likely does NOT ",
          "match your actual file's name-field format and needs fixing)")
  out
}

chiapet_loops <- load_chiapet_bed12(CHIAPET_BED12)

# ══════════════════════════════════════════════════════════════════════════════
# 2. Anchor classification: exon / intron / cCRE
#    (reused as-is from hic_general_analysis.R's classify_loop_anchors())
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

classify_loop_anchors <- function(chiapet_loops, exon_gr, intron_gr, ccre_gr=NULL,
                                  anchor_window=5000L) {
  chiapet_loops <- chiapet_loops %>% dplyr::mutate(loop_idx=row_number())
  fix_chr <- function(gr, reference_gr) {
    ref_chroms <- as.character(seqnames(reference_gr)); has_chr_ref <- any(grepl("^chr", ref_chroms))
    gr_chroms  <- as.character(seqnames(gr)); has_chr_gr <- any(grepl("^chr", gr_chroms))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr", seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
    gr
  }
  # Tight, resolution-sized window centered on each anchor's own midpoint --
  # NOT the raw start1:end1/start2:end2 span, which can be much wider than
  # the loop's true anchor and would act as an unintended tolerance against
  # small features like exons (same rationale as the Hi-C HiCCUPS classifier).
  half_w <- anchor_window / 2
  anchor1_gr <- GRanges(seqnames=chiapet_loops$chr1,
                        ranges=IRanges(chiapet_loops$centroid1-half_w, chiapet_loops$centroid1+half_w))
  anchor2_gr <- GRanges(seqnames=chiapet_loops$chr2,
                        ranges=IRanges(chiapet_loops$centroid2-half_w, chiapet_loops$centroid2+half_w))
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
  chiapet_loops %>% dplyr::mutate(
    anchor1_class=anchor1_class, anchor2_class=anchor2_class,
    loop_anchor_type=factor(loop_type, levels=c("Exon-Exon","Exon-cCRE","Exon-Intron",
                                                 "cCRE-cCRE","Intron-Intron","Other")))
}

chiapet_classified <- classify_loop_anchors(chiapet_loops, exon_gr, intron_gr, ccre_gr=ccre_gr)

message("\nLoop anchor type breakdown:")
print(table(chiapet_classified$loop_anchor_type))

# ══════════════════════════════════════════════════════════════════════════════
# 3. Plot (reused as-is from hic_general_analysis.R's plot_loop_strength_raincloud())
# ══════════════════════════════════════════════════════════════════════════════

ANCHOR_TYPE_COLORS <- c("Exon-Exon"="#d6604d", "Exon-cCRE"="#984ea3", "Exon-Intron"="#4DAF4A",
                        "cCRE-cCRE"="#ff7f00", "Intron-Intron"="#2166ac")
ANCHOR_TYPE_LEVELS <- names(ANCHOR_TYPE_COLORS)

plot_loop_strength_raincloud <- function(classified_loops) {
  plot_df <- classified_loops %>%
    filter(loop_anchor_type %in% ANCHOR_TYPE_LEVELS, !is.na(pet_count), pet_count > 0) %>%
    dplyr::mutate(log10_score_capped=pmin(log10(pet_count), 2.15),
                  loop_anchor_type=factor(loop_anchor_type, levels=ANCHOR_TYPE_LEVELS))
  counts <- plot_df %>% dplyr::count(loop_anchor_type) %>%
    dplyr::mutate(loop_anchor_type=factor(loop_anchor_type, levels=ANCHOR_TYPE_LEVELS),
                  label=paste0(loop_anchor_type, "\n(n=", n, ")")) %>%
    dplyr::arrange(loop_anchor_type)
  plot_df <- plot_df %>% left_join(counts %>% dplyr::select(loop_anchor_type, label), by="loop_anchor_type") %>%
    dplyr::mutate(label=factor(label, levels=counts$label))

  ggplot(plot_df, aes(x=label, y=log10_score_capped, fill=loop_anchor_type, color=loop_anchor_type)) +
    ggdist::stat_halfeye(aes(thickness=after_stat(pdf * n)), side="right", scale=0.4, width=0.6,
                         .width=0, point_colour=NA, alpha=0.75, normalize="groups",
                         position=position_nudge(x=0.05)) +
    geom_boxplot(width=0.12, outlier.shape=NA, color="gray20", alpha=0.9, linewidth=0.6,
                position=position_nudge(x=-0.07)) +
    scale_fill_manual(values=ANCHOR_TYPE_COLORS, guide="none") +
    scale_color_manual(values=ANCHOR_TYPE_COLORS, guide="none") +
    scale_y_continuous(name="Loop strength (PET count)", limits=c(NA, 2.15),
                       breaks=c(log10(2), log10(5), 1, log10(20), log10(50), 2),
                       labels=c("2","5","10","20","50","100"), expand=expansion(mult=c(0.02, 0))) +
    scale_x_discrete(expand=expansion(add=c(0.3, 0.5))) +
    labs(title="POLR2A ChIA-PET loop strength by anchor type", x=NULL) +
    theme_classic(base_size=13) +
    theme(plot.title=element_text(face="bold", size=14),
          axis.text.x=element_text(size=10, angle=30, hjust=1, vjust=1),
          plot.margin=margin(8, 16, 8, 8), aspect.ratio=1.4)
}

fig5f <- plot_loop_strength_raincloud(chiapet_classified)
print(fig5f)
