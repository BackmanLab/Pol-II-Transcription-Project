##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# 1D insulation score pipeline. 4-way TAD analysis: WT/POLR2A-degraded
# x CTCF-anchored/non-CTCF. For each subset, TAD positions fixed,
# insulation scored in all conditions. Expression bins: Low
# (incl. unexpressed) / Medium / High.
#
# Figure output: Fig. 4A-B (probable); SI Fig. 7A-F (probable).
##################################################################
library(strawr)
library(data.table)
library(parallel)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(GenomicRanges)
library(zoo)

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════

DS_GROUP <- "ds_wt_pol2"   # ← "ds_wt_pol2" | "ds_all" | "unsampled"

BASE_NEW <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
BASE_OLD <- "/projects/p32171/juicer/work/112123_HiC/juicer_analysis"

if (DS_GROUP == "ds_wt_pol2") {
  wt_hic_path   <- file.path(BASE_NEW, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic_path <- file.path(BASE_NEW, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  actd_hic_path <- NULL
  wt_tads       <- file.path(BASE_NEW, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")
  pol2_tads     <- file.path(BASE_NEW, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")
  actd_tads     <- NULL
  CONDITIONS_USE <- c("wt","pol2")
  
} else if (DS_GROUP == "ds_all") {
  wt_hic_path   <- file.path(BASE_NEW, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic_path <- file.path(BASE_NEW, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  actd_hic_path <- file.path(BASE_NEW, "mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic")
  wt_tads       <- file.path(BASE_NEW, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")
  pol2_tads     <- file.path(BASE_NEW, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")
  actd_tads     <- file.path(BASE_NEW, "mega_ds_all/1hr_ActD/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")
  CONDITIONS_USE <- c("wt","pol2","actd")
  
} else if (DS_GROUP == "unsampled") {
  wt_hic_path   <- file.path(BASE_OLD, "WT_HCT116_CTRL/mega/aligned/inter_30.hic")
  pol2_hic_path <- file.path(BASE_OLD, "Pol2_6hrs_Aux/mega/aligned/inter_30.hic")
  actd_hic_path <- file.path(BASE_OLD, "1hr_ActD/mega/aligned/inter_30.hic")
  wt_tads       <- file.path(BASE_OLD, "WT_HCT116_CTRL/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")
  pol2_tads     <- file.path(BASE_OLD, "Pol2_6hrs_Aux/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")
  actd_tads     <- file.path(BASE_OLD, "1hr_ActD/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")
  CONDITIONS_USE <- c("wt","pol2","actd")
  
} else {
  stop("Unknown DS_GROUP: '", DS_GROUP, "'")
}
message("DS_GROUP: ", DS_GROUP)

EXPR_FILE   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
DESEQ_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
ENCODE_DIR  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/encode"
RDS_DIR     <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/insulation"

RESOLUTION    <- 25000
INS_WINDOW    <- 10
INS_FLANK_KB  <- 500
ROLL_WINDOW   <- 5
TAD_MIN_KB    <- 50
TAD_MAX_KB    <- 5000
GENE_TOL_BP   <- 5000
CTCF_TOL_BP   <- 10000
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1
STRAW_NORM    <- "SCALE"

get_slurm_cores <- function(default=2) {
  cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (nchar(cpus) > 0) return(as.integer(cpus))
  ntasks <- Sys.getenv("SLURM_NTASKS")
  if (nchar(ntasks) > 0) return(as.integer(ntasks))
  default
}
N_CORES <- max(1L, get_slurm_cores(default=2) - 1L)
message("Using ", N_CORES, " cores")

ALL_CONDITIONS <- list(
  wt   = list(label="Wildtype",   hic=wt_hic_path,   color="#88CCEE"),
  pol2 = list(label="POLR2A deg", hic=pol2_hic_path, color="#CC6677"),
  actd = list(label="ActD 1hr",   hic=actd_hic_path, color="#4DAF4A")
)
CONDITIONS  <- ALL_CONDITIONS[CONDITIONS_USE]
COND_LEVELS <- sapply(CONDITIONS, `[[`, "label")
COND_COLORS <- setNames(sapply(CONDITIONS, `[[`, "color"), COND_LEVELS)

# Expression bins — Low includes unexpressed (TPM=0)
BINS       <- c("Low","Medium","High")
BIN_COLORS <- c("Low"="#fdae61", "Medium"="#f46d43", "High"="#d73027")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Expression data
# ══════════════════════════════════════════════════════════════════════════════
expr <- read_tsv(EXPR_FILE) %>% dplyr::select(gene_id, TPM)

deseq_res <- read_csv(DESEQ_FILE) %>%
  dplyr::mutate(
    gene_id = sub("\\..*", "", gene_id),
    symbol  = dplyr::case_when(!is.na(symbol) & symbol != "" ~ symbol, TRUE ~ gene_id)
  ) %>%
  left_join(expr, by="gene_id")

expressed  <- deseq_res %>% filter(!is.na(TPM), TPM > 0) %>% pull(TPM)
tpm_medium <- quantile(expressed, 1/3, na.rm=TRUE)
tpm_high   <- quantile(expressed, 2/3, na.rm=TRUE)

message("Expression bin thresholds (WT TPM tertiles):")
message("  Low (incl. unexpressed): TPM <= ", round(tpm_medium, 2))
message("  Medium: ", round(tpm_medium, 2), " < TPM <= ", round(tpm_high, 2))
message("  High:   TPM > ", round(tpm_high, 2))

deseq_res <- deseq_res %>%
  dplyr::mutate(
    expr_bin = dplyr::case_when(
      is.na(TPM) | TPM == 0 ~ "Low",
      TPM <= tpm_medium     ~ "Low",
      TPM <= tpm_high       ~ "Medium",
      TRUE                  ~ "High"
    ),
    expr_bin = factor(expr_bin, levels=BINS)
  )
message("Expression bin counts:"); print(table(deseq_res$expr_bin))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Gene GRanges
# ══════════════════════════════════════════════════════════════════════════════
gene_coords <- deseq_res %>%
  filter(!is.na(gene_seq_start), !is.na(gene_seq_end), !is.na(seq_name)) %>%
  distinct(gene_id, .keep_all=TRUE)

gene_gr <- GRanges(
  seqnames = as.character(gene_coords$seq_name),
  ranges   = IRanges(gene_coords$gene_seq_start, gene_coords$gene_seq_end),
  symbol   = gene_coords$symbol,
  gene_id  = gene_coords$gene_id,
  expr_bin = as.character(gene_coords$expr_bin),
  TPM      = gene_coords$TPM
)
message("gene_gr: ", length(gene_gr), " genes")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Load TADs
# ══════════════════════════════════════════════════════════════════════════════
load_arrowhead <- function(path, min_size_kb=TAD_MIN_KB, max_size_kb=TAD_MAX_KB) {
  all_lines  <- readLines(path)
  data_lines <- all_lines[!grepl("^#", all_lines)]
  df <- read.table(text=paste(data_lines, collapse="\n"),
                   header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
  colnames(df)[1:6] <- c("chr1","x1","x2","chr2","y1","y2")
  df %>%
    dplyr::mutate(tad_size_kb=(x2-x1)/1000) %>%
    filter(!is.na(x1), !is.na(x2), x2>x1,
           tad_size_kb>=min_size_kb, tad_size_kb<=max_size_kb)
}

tads_raw <- list(wt=load_arrowhead(wt_tads), pol2=load_arrowhead(pol2_tads))
if (!is.null(actd_tads)) tads_raw$actd <- load_arrowhead(actd_tads)
invisible(lapply(names(tads_raw), function(n)
  message(CONDITIONS[[n]]$label, " TADs loaded: ", nrow(tads_raw[[n]]))))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Annotate TADs by expression bin (boundary ±GENE_TOL_BP)
# ══════════════════════════════════════════════════════════════════════════════
link_tads_anchor <- function(tads_df, gene_gr, tol_bp=GENE_TOL_BP) {
  tads_df <- tads_df %>% dplyr::mutate(tad_idx=row_number())
  fix_chr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr",seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr","",seqlevels(gr))
    gr
  }
  gene_gr_exp <- trim(resize(gene_gr, width(gene_gr)+2*tol_bp, fix="center"))
  left_gr  <- fix_chr(GRanges(tads_df$chr1,
                              IRanges(pmax(0,tads_df$x1-tol_bp), tads_df$x1+tol_bp)), gene_gr)
  right_gr <- fix_chr(GRanges(tads_df$chr1,
                              IRanges(pmax(0,tads_df$x2-tol_bp), tads_df$x2+tol_bp)), gene_gr)
  hits_left  <- findOverlaps(left_gr,  gene_gr_exp)
  hits_right <- findOverlaps(right_gr, gene_gr_exp)
  all_pairs <- bind_rows(
    data.frame(tad_idx=queryHits(hits_left),  gene_idx=subjectHits(hits_left)),
    data.frame(tad_idx=queryHits(hits_right), gene_idx=subjectHits(hits_right))
  ) %>%
    distinct(tad_idx, gene_idx) %>%
    dplyr::mutate(
      expr_bin = gene_gr$expr_bin[gene_idx],
      TPM      = gene_gr$TPM[gene_idx],
      symbol   = gene_gr$symbol[gene_idx]
    )
  tad_annotations <- all_pairs %>%
    dplyr::group_by(tad_idx) %>%
    dplyr::slice_max(TPM, n=1, with_ties=FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(tad_idx, expr_bin, TPM, symbol)
  tads_df %>%
    left_join(tad_annotations, by="tad_idx") %>%
    dplyr::mutate(
      expr_bin = replace_na(expr_bin, "Low"),   # unannotated → Low
      symbol   = replace_na(symbol, "Unknown")
    )
}

message("Annotating TADs by expression bin...")
tads_ann <- lapply(tads_raw, link_tads_anchor, gene_gr=gene_gr)
invisible(lapply(names(tads_ann), function(n) {
  message(CONDITIONS[[n]]$label, " TAD expression breakdown:")
  print(table(tads_ann[[n]]$expr_bin))
}))

# ══════════════════════════════════════════════════════════════════════════════
# 5. CTCF filtering — split each condition's TADs into CTCF / non-CTCF
# ══════════════════════════════════════════════════════════════════════════════
ctcf_sites <- read.table(file.path(ENCODE_DIR, "wt_hct116.bed"),
                         header=FALSE, sep="\t", stringsAsFactors=FALSE) %>%
  dplyr::select(chr=V1, start=V2, end=V3, id=V4, ccre_type=V10) %>%
  filter(ccre_type=="CA-CTCF")
ctcf_gr <- GRanges(seqnames=ctcf_sites$chr,
                   ranges=IRanges(ctcf_sites$start, ctcf_sites$end))
message("CA-CTCF sites: ", length(ctcf_gr))

filter_ctcf_tads <- function(tads_df, ctcf_gr, tol_bp=CTCF_TOL_BP) {
  fix_chr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr",seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr","",seqlevels(gr))
    gr
  }
  left_gr  <- fix_chr(GRanges(tads_df$chr1,
                              IRanges(pmax(0,tads_df$x1-tol_bp), tads_df$x1+tol_bp)), ctcf_gr)
  right_gr <- fix_chr(GRanges(tads_df$chr1,
                              IRanges(pmax(0,tads_df$x2-tol_bp), tads_df$x2+tol_bp)), ctcf_gr)
  left_hits  <- unique(queryHits(findOverlaps(left_gr,  ctcf_gr)))
  right_hits <- unique(queryHits(findOverlaps(right_gr, ctcf_gr)))
  both_ctcf  <- intersect(left_hits, right_hits)
  list(
    ctcf     = tads_df[both_ctcf, ],
    non_ctcf = tads_df[setdiff(seq_len(nrow(tads_df)), both_ctcf), ]
  )
}

message("\nFiltering TADs by CTCF anchor status...")
message("WT TADs:")
wt_split   <- filter_ctcf_tads(tads_ann$wt,   ctcf_gr)
message("  CTCF-anchored:     ", nrow(wt_split$ctcf))
message("  Non-CTCF-anchored: ", nrow(wt_split$non_ctcf))

message("pol2 TADs:")
pol2_split <- filter_ctcf_tads(tads_ann$pol2, ctcf_gr)
message("  CTCF-anchored:     ", nrow(pol2_split$ctcf))
message("  Non-CTCF-anchored: ", nrow(pol2_split$non_ctcf))

# The 4 TAD subsets to analyse
TAD_SUBSETS <- list(
  wt_non_ctcf   = list(tads=wt_split$non_ctcf,   label="WT TADs — non-CTCF",   source="wt"),
  pol2_non_ctcf = list(tads=pol2_split$non_ctcf, label="pol2 TADs — non-CTCF", source="pol2"),
  wt_ctcf       = list(tads=wt_split$ctcf,       label="WT TADs — CTCF",       source="wt"),
  pol2_ctcf     = list(tads=pol2_split$ctcf,     label="pol2 TADs — CTCF",     source="pol2")
)

# ══════════════════════════════════════════════════════════════════════════════
# 6. Compute insulation scores
# ══════════════════════════════════════════════════════════════════════════════
compute_genome_insulation <- function(hic_file, sample_name,
                                      norm=STRAW_NORM, resolution=RESOLUTION,
                                      window=INS_WINDOW, n_cores=N_CORES) {
  .compute_chr <- function(hic_file, chrom, resolution, window, norm) {
    library(strawr); library(data.table)
    raw <- tryCatch(
      strawr::straw(norm, hic_file, as.character(chrom), as.character(chrom), "BP", resolution),
      error=function(e) NULL
    )
    if (is.null(raw) || nrow(raw)==0) return(NULL)
    dt <- data.table::as.data.table(raw)
    dt[, bin1 := as.integer(x/resolution)]
    dt[, bin2 := as.integer(y/resolution)]
    all_bins <- sort(unique(c(dt$bin1, dt$bin2)))
    if (length(all_bins) < 2*window+1) return(NULL)
    dt[, b_lo := pmax(bin1, bin2-window)+1L]
    dt[, b_hi := pmin(bin1+window-1L, bin2)]
    dt_contrib <- dt[b_lo <= b_hi]
    if (nrow(dt_contrib)==0) return(NULL)
    dt_exp   <- dt_contrib[, .(b=seq.int(b_lo,b_hi), counts=counts), by=seq_len(nrow(dt_contrib))]
    win_sums <- dt_exp[, .(window_sum=sum(counts, na.rm=TRUE)), by=b]
    result   <- merge(data.table(b=all_bins), win_sums, by="b", all.x=TRUE)
    data.frame(chr=as.character(chrom), start=result$b*resolution,
               end=result$b*resolution+resolution, bin_idx=result$b,
               window_sum=result$window_sum, stringsAsFactors=FALSE)
  }
  chroms <- strawr::readHicChroms(hic_file) %>%
    filter(name!="ALL", name!="M", name!="Y", !grepl("_",name)) %>%
    filter(grepl("^[0-9]+$",name) | name=="X") %>% pull(name)
  message("Computing insulation scores — ", sample_name,
          "  norm=", norm, "  res=", resolution/1000, "kb",
          "  window=", window, " bins  cores=", n_cores)
  results_list <- parallel::mclapply(chroms,
                                     function(chrom) .compute_chr(hic_file, chrom, resolution, window, norm),
                                     mc.cores=n_cores, mc.preschedule=FALSE)
  message("  ", sum(!sapply(results_list,is.null)), "/", length(chroms), " chromosomes OK")
  results <- bind_rows(Filter(Negate(is.null), results_list))
  if (nrow(results)==0) return(NULL)
  results %>% dplyr::mutate(ins_score=window_sum, sample=sample_name)
}

dir.create(RDS_DIR, showWarnings=FALSE, recursive=TRUE)
ins_rds_raw <- file.path(RDS_DIR, paste0("insulation_raw_",     DS_GROUP, "_", STRAW_NORM, ".rds"))
ins_rds     <- file.path(RDS_DIR, paste0("insulation_zscored_", DS_GROUP, "_", STRAW_NORM, ".rds"))

if (file.exists(ins_rds)) {
  message("Loading cached z-scored insulation scores: ", basename(ins_rds))
  ins_all <- readRDS(ins_rds)
} else {
  if (file.exists(ins_rds_raw)) {
    message("Loading cached raw scores...")
    ins_all <- readRDS(ins_rds_raw)
  } else {
    ins_list <- lapply(CONDITIONS_USE, function(n)
      compute_genome_insulation(CONDITIONS[[n]]$hic, CONDITIONS[[n]]$label,
                                norm=STRAW_NORM, n_cores=N_CORES))
    ins_all <- bind_rows(ins_list) %>%
      dplyr::mutate(sample=factor(sample, levels=COND_LEVELS))
    saveRDS(ins_all, ins_rds_raw)
    message("Raw scores saved: ", ins_rds_raw)
  }
  message("Z-score normalizing per chromosome per condition...")
  ins_all <- ins_all %>%
    dplyr::group_by(sample, chr) %>%
    dplyr::mutate(ins_zscore=(ins_score-mean(ins_score,na.rm=TRUE))/sd(ins_score,na.rm=TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample=factor(sample, levels=COND_LEVELS))
  saveRDS(ins_all, ins_rds)
  message("Z-scored scores saved: ", ins_rds)
}
message("\nInsulation score summary:")
print(ins_all %>% dplyr::group_by(sample) %>%
        dplyr::summarise(n=dplyr::n(),
                         mean=round(mean(ins_zscore,na.rm=TRUE),3),
                         sd  =round(sd(ins_zscore,na.rm=TRUE),3),
                         .groups="drop"))

# ══════════════════════════════════════════════════════════════════════════════
# 7. Boundary profile extraction
# ══════════════════════════════════════════════════════════════════════════════
get_boundary_profiles <- function(ins_all, tads_df,
                                  flank_kb=INS_FLANK_KB,
                                  resolution=RESOLUTION,
                                  roll_window=ROLL_WINDOW) {
  flank_bp   <- flank_kb * 1000
  boundaries <- bind_rows(
    tads_df %>% dplyr::mutate(boundary=x1, side="left"),
    tads_df %>% dplyr::mutate(boundary=x2, side="right")
  ) %>% dplyr::mutate(chr=as.character(chr1))
  
  profiles <- bind_rows(lapply(COND_LEVELS, function(samp) {
    ins_cond <- ins_all %>% filter(sample==samp) %>% dplyr::select(chr, start, ins_zscore)
    bp_data  <- bind_rows(lapply(seq_len(nrow(boundaries)), function(i) {
      ins_cond %>%
        filter(chr==boundaries$chr[i],
               start>=boundaries$boundary[i]-flank_bp,
               start<=boundaries$boundary[i]+flank_bp) %>%
        dplyr::mutate(offset_kb=(start-boundaries$boundary[i])/1000) %>%
        dplyr::select(offset_kb, ins_zscore)
    }))
    if (is.null(bp_data) || nrow(bp_data)==0) return(NULL)
    bp_data %>%
      dplyr::group_by(offset_kb) %>%
      dplyr::summarise(
        mean_score=mean(ins_zscore, na.rm=TRUE),
        se_score  =sd(ins_zscore,   na.rm=TRUE)/sqrt(sum(!is.na(ins_zscore))),
        n         =sum(!is.na(ins_zscore)),
        .groups   ="drop"
      ) %>%
      dplyr::mutate(sample=samp) %>%
      dplyr::arrange(offset_kb) %>%
      dplyr::mutate(
        mean_score=zoo::rollmean(mean_score, k=roll_window, fill=NA, align="center", na.rm=TRUE),
        se_score  =zoo::rollmean(se_score,   k=roll_window, fill=NA, align="center", na.rm=TRUE)
      )
  }))
  profiles %>% dplyr::mutate(sample=factor(sample, levels=COND_LEVELS))
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Plot functions
# ══════════════════════════════════════════════════════════════════════════════
make_insulation_plot <- function(profiles, plot_title, plot_subtitle,
                                 bin_color="gray30") {
  if (is.null(profiles) || nrow(profiles)==0) {
    message("No data: ", plot_title); return(NULL)
  }
  ggplot(profiles, aes(x=offset_kb, y=mean_score, color=sample, fill=sample)) +
    geom_ribbon(aes(ymin=mean_score-se_score, ymax=mean_score+se_score),
                alpha=0.15, color=NA) +
    geom_line(linewidth=0.9) +
    geom_vline(xintercept=0, linetype="dashed", color="gray40", linewidth=0.5) +
    scale_color_manual(values=COND_COLORS, name=NULL) +
    scale_fill_manual(values=COND_COLORS,  name=NULL) +
    scale_x_continuous(name="Distance from TAD boundary (kb)",
                       breaks=seq(-INS_FLANK_KB, INS_FLANK_KB, by=200)) +
    scale_y_continuous(name="Mean insulation score\n(z-score, per chromosome)") +
    labs(title=plot_title, subtitle=plot_subtitle) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold", size=10, color=bin_color, hjust=0.5),
          plot.subtitle=element_text(size=8, color="gray50", hjust=0.5),
          legend.position="none", plot.margin=margin(4,8,4,8))
}

shared_legend <- cowplot::get_legend(
  ggplot(data.frame(x=1, y=1, sample=factor(COND_LEVELS, levels=COND_LEVELS)),
         aes(x=x, y=y, color=sample)) +
    geom_line(linewidth=1) +
    scale_color_manual(values=COND_COLORS, name=NULL) +
    theme_void() +
    theme(legend.position="bottom", legend.text=element_text(size=10),
          legend.key.size=unit(1,"cm"))
)

assemble_insulation <- function(plots, title, subtitle) {
  plots <- Filter(Negate(is.null), plots)
  if (length(plots)==0) return(invisible(NULL))
  row  <- wrap_plots(plots, nrow=1) +
    plot_layout(guides="collect") & theme(legend.position="none")
  comp <- cowplot::plot_grid(row, shared_legend, ncol=1, rel_heights=c(1,0.07))
  ttl  <- cowplot::ggdraw() +
    cowplot::draw_label(title,    fontface="bold", size=13, x=0.01, hjust=0) +
    cowplot::draw_label(subtitle, size=8, color="gray40", x=0.01, hjust=0, y=0.2)
  final <- cowplot::plot_grid(ttl, comp, ncol=1, rel_heights=c(0.07,1))
  print(final); invisible(final)
}

subtitle_base <- paste0(
  "\u00b1", INS_FLANK_KB, "kb around TAD boundaries  |  ",
  "window=", INS_WINDOW, " bins (", INS_WINDOW*RESOLUTION/1000, "kb)  |  ",
  "z-score per chr  |  rolling mean k=", ROLL_WINDOW,
  "  |  ribbon=\u00b11 SE  |  norm=", STRAW_NORM, "  |  ", DS_GROUP
)

# ══════════════════════════════════════════════════════════════════════════════
# 9. Core analysis runner
#    For each TAD subset: score ALL conditions at FIXED TAD positions,
#    plot globally and by expression bin.
# ══════════════════════════════════════════════════════════════════════════════
run_subset_analysis <- function(subset_name, tads_df, subset_label) {
  if (nrow(tads_df) < 10) {
    message("Skipping ", subset_label, " — too few TADs (n=", nrow(tads_df), ")")
    return(invisible(NULL))
  }
  message("\n══ ", subset_label, " (n=", nrow(tads_df), " TADs) ════════════════")
  sub <- paste0(subset_label, "  |  ", nrow(tads_df), " TADs  |  ", subtitle_base)
  
  # Global — all TADs in subset
  assemble_insulation(
    list(make_insulation_plot(get_boundary_profiles(ins_all, tads_df),
                              paste0(subset_label, " (global)"), sub)),
    title    = paste0("Insulation — ", subset_label),
    subtitle = sub
  )
  
  # By expression bin
  assemble_insulation(
    lapply(BINS, function(bin) {
      bin_tads <- tads_df %>% filter(expr_bin == bin)
      if (nrow(bin_tads) < 10) return(NULL)
      message("  Bin ", bin, ": ", nrow(bin_tads), " TADs")
      make_insulation_plot(
        get_boundary_profiles(ins_all, bin_tads),
        plot_title    = bin,
        plot_subtitle = paste0("n=", nrow(bin_tads), " TADs"),
        bin_color     = BIN_COLORS[bin]
      )
    }),
    title    = paste0("Insulation — ", subset_label, " — by expression"),
    subtitle = sub
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. Run all 4 subsets
# ══════════════════════════════════════════════════════════════════════════════
for (nm in names(TAD_SUBSETS)) {
  s <- TAD_SUBSETS[[nm]]
  run_subset_analysis(nm, s$tads, s$label)
}

message("\nAll insulation analyses complete.")

# ══════════════════════════════════════════════════════════════════════════════
# Pure CTCF TADs — CTCF at both boundaries, no Promoter or Enhancer cCRE,
# but must have a gene nearby (already ensured by link_tads_anchor annotation)
# ══════════════════════════════════════════════════════════════════════════════

# Load promoter and enhancer cCRE sites from WT ENCODE file
ccre_all <- read.table(file.path(ENCODE_DIR, "wt_hct116.bed"),
                       header=FALSE, sep="\t", stringsAsFactors=FALSE) %>%
  dplyr::select(chr=V1, start=V2, end=V3, ccre_type=V10)

promoter_gr <- GRanges(
  seqnames = ccre_all %>% filter(ccre_type=="PLS") %>% pull(chr),
  ranges   = IRanges(
    ccre_all %>% filter(ccre_type=="PLS") %>% pull(start),
    ccre_all %>% filter(ccre_type=="PLS") %>% pull(end)
  )
)
enhancer_gr <- GRanges(
  seqnames = ccre_all %>% filter(ccre_type %in% c("pELS","dELS")) %>% pull(chr),
  ranges   = IRanges(
    ccre_all %>% filter(ccre_type %in% c("pELS","dELS")) %>% pull(start),
    ccre_all %>% filter(ccre_type %in% c("pELS","dELS")) %>% pull(end)
  )
)
message("Promoter cCREs (PLS): ",      length(promoter_gr))
message("Enhancer cCREs (pELS/dELS): ", length(enhancer_gr))

filter_pure_ctcf_tads <- function(ctcf_tads, promoter_gr, enhancer_gr,
                                  tol_bp=CTCF_TOL_BP) {
  fix_chr <- function(gr, ref_gr) {
    has_chr_ref <- any(grepl("^chr", as.character(seqnames(ref_gr))))
    has_chr_gr  <- any(grepl("^chr", as.character(seqnames(gr))))
    if (has_chr_ref && !has_chr_gr) seqlevels(gr) <- paste0("chr",seqlevels(gr))
    else if (!has_chr_ref && has_chr_gr) seqlevels(gr) <- gsub("^chr","",seqlevels(gr))
    gr
  }
  left_gr  <- fix_chr(GRanges(ctcf_tads$chr1,
                              IRanges(pmax(0,ctcf_tads$x1-tol_bp), ctcf_tads$x1+tol_bp)),
                      promoter_gr)
  right_gr <- fix_chr(GRanges(ctcf_tads$chr1,
                              IRanges(pmax(0,ctcf_tads$x2-tol_bp), ctcf_tads$x2+tol_bp)),
                      promoter_gr)
  
  # TADs where either boundary hits a promoter or enhancer
  prom_hits <- unique(c(
    queryHits(findOverlaps(left_gr,  promoter_gr)),
    queryHits(findOverlaps(right_gr, promoter_gr))
  ))
  enh_hits  <- unique(c(
    queryHits(findOverlaps(left_gr,  enhancer_gr)),
    queryHits(findOverlaps(right_gr, enhancer_gr))
  ))
  exclude   <- unique(c(prom_hits, enh_hits))
  keep      <- setdiff(seq_len(nrow(ctcf_tads)), exclude)
  
  message("  CTCF TADs input:              ", nrow(ctcf_tads))
  message("  Excluded (promoter boundary): ", length(prom_hits))
  message("  Excluded (enhancer boundary): ", length(enh_hits))
  message("  Pure CTCF remaining:          ", length(keep))
  
  ctcf_tads[keep, ]
}

message("\nFiltering pure CTCF TADs (no promoter or enhancer at boundaries)...")
message("WT:")
wt_pure_ctcf   <- filter_pure_ctcf_tads(wt_split$ctcf,   promoter_gr, enhancer_gr)
message("pol2:")
pol2_pure_ctcf <- filter_pure_ctcf_tads(pol2_split$ctcf, promoter_gr, enhancer_gr)

# Require gene overlap — drop TADs where expr_bin is Low from no annotation
# (keep Low-expression genes but not completely unannotated boundaries)
wt_pure_ctcf_gene   <- wt_pure_ctcf   %>% filter(!is.na(symbol), symbol != "Unknown")
pol2_pure_ctcf_gene <- pol2_pure_ctcf %>% filter(!is.na(symbol), symbol != "Unknown")

message("WT pure CTCF with gene overlap:   ", nrow(wt_pure_ctcf_gene))
message("pol2 pure CTCF with gene overlap: ", nrow(pol2_pure_ctcf_gene))

# Add to TAD_SUBSETS
TAD_SUBSETS$wt_pure_ctcf   <- list(tads=wt_pure_ctcf_gene,
                                   label="WT TADs — pure CTCF (no promoter/enhancer)",
                                   source="wt")
TAD_SUBSETS$pol2_pure_ctcf <- list(tads=pol2_pure_ctcf_gene,
                                   label="pol2 TADs — pure CTCF (no promoter/enhancer)",
                                   source="pol2")

run_subset_analysis("wt_pure_ctcf",
                    wt_pure_ctcf_gene,
                    "WT TADs — pure CTCF (no promoter/enhancer)")

run_subset_analysis("pol2_pure_ctcf",
                    pol2_pure_ctcf_gene,
                    "pol2 TADs — pure CTCF (no promoter/enhancer)")
