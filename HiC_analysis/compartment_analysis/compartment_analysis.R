##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Coding/noncoding compartment saddle analysis, pipeline 2.
# DS_GROUP switch for downsampled/unsampled Hi-C files. Eigenvectors
# loaded from bedgraph files.
#
# Three series:
#   1. E/I boxplot marginals (per coding bin, separate figures)
#   2. PC1 bar marginals (per coding bin, separate figures)
#   3. Expression-stratified saddle with PC1 marginals
#
# Figure output: SI Fig. 6C-E (probable).
##################################################################
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(GenomicRanges)
library(strawr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(ggside)
library(ggnewscale)
library(scales)
library(parallel)

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════
DS_GROUP <- "ds_wt_pol2"   # ← "ds_wt_pol2" | "ds_all" | "unsampled"

BASE_NEW  <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
BASE_OLD  <- "/projects/p32171/juicer/work/112123_HiC/juicer_analysis"
COMP_BASE <- file.path(BASE_NEW, "compartments")

if (DS_GROUP == "ds_wt_pol2") {
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE",
                hic     =file.path(BASE_NEW,"mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                comp_bed=file.path(COMP_BASE,"ds_wt_pol2_WT/ds_wt_pol2_WT_compartments.bedgraph")),
    pol2 = list(label="POLR2A deg", color="#CC6677",
                hic     =file.path(BASE_NEW,"mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                comp_bed=file.path(COMP_BASE,"ds_wt_pol2_Pol2/ds_wt_pol2_Pol2_compartments.bedgraph"))
  )
} else if (DS_GROUP == "ds_all") {
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE",
                hic     =file.path(BASE_NEW,"mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                comp_bed=file.path(COMP_BASE,"ds_all_WT/ds_all_WT_compartments.bedgraph")),
    pol2 = list(label="POLR2A deg", color="#CC6677",
                hic     =file.path(BASE_NEW,"mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                comp_bed=file.path(COMP_BASE,"ds_all_Pol2/ds_all_Pol2_compartments.bedgraph")),
    actd = list(label="ActD 1hr",   color="#44AA99",
                hic     =file.path(BASE_NEW,"mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic"),
                comp_bed=file.path(COMP_BASE,"ds_all_ActD/ds_all_ActD_compartments.bedgraph"))
  )
} else if (DS_GROUP == "unsampled") {
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE",
                hic     =file.path(BASE_OLD,"WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                comp_bed=NULL),
    pol2 = list(label="POLR2A deg", color="#CC6677",
                hic     =file.path(BASE_OLD,"Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                comp_bed=NULL),
    actd = list(label="ActD 1hr",   color="#44AA99",
                hic     =file.path(BASE_OLD,"1hr_ActD/mega/aligned/inter_30.hic"),
                comp_bed=file.path(COMP_BASE,"unsampled_ActD/unsampled_ActD_compartments.bedgraph"))
  )
} else { stop("Unknown DS_GROUP: '", DS_GROUP, "'") }
message("DS_GROUP: ", DS_GROUP)

COND_LEVELS <- sapply(CONDITIONS, `[[`, "label")
COND_COLORS <- setNames(sapply(CONDITIONS, `[[`, "color"), COND_LEVELS)

EXPR_FILE      <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
DESEQ_FILE     <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
RDS_DIR        <- file.path(BASE_NEW, "compartments/rds", DS_GROUP)
dir.create(RDS_DIR, showWarnings=FALSE, recursive=TRUE)

RESOLUTION     <- 500000
SADDLE_BINS    <- 20L
STRAW_NORM     <- "SCALE"
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1

# Colorblind-safe divergent pair
COL_A    <- "#E69F00"   # amber  — A compartment / exon-enriched
COL_B    <- "#56B4E9"   # sky blue — B compartment / intron-enriched

get_slurm_cores <- function(default=2) {
  cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (nchar(cpus)>0) { message("SLURM_CPUS_PER_TASK: ",as.integer(cpus)); return(as.integer(cpus)) }
  ntasks <- Sys.getenv("SLURM_NTASKS")
  if (nchar(ntasks)>0) { message("SLURM_NTASKS: ",as.integer(ntasks)); return(as.integer(ntasks)) }
  message("No SLURM — using ",default," cores"); return(default)
}
N_CORES <- get_slurm_cores(default=2)
BINS    <- c("Low","Medium","High")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Expression data
# ══════════════════════════════════════════════════════════════════════════════
message("Loading expression data...")
expr <- read_tsv(EXPR_FILE, show_col_types=FALSE) %>%
  dplyr::select(gene_id, TPM)
deseq_res <- read_csv(DESEQ_FILE, show_col_types=FALSE) %>%
  dplyr::mutate(
    gene_id = sub("\\..*","",gene_id),
    symbol  = dplyr::case_when(!is.na(symbol)&symbol!=""~symbol, TRUE~gene_id)
  ) %>% left_join(expr, by="gene_id")
expressed  <- deseq_res %>% filter(!is.na(TPM),TPM>0) %>% pull(TPM)
tpm_medium <- quantile(expressed,1/3,na.rm=TRUE)
tpm_high   <- quantile(expressed,2/3,na.rm=TRUE)
deseq_res  <- deseq_res %>%
  dplyr::mutate(
    expr_bin = dplyr::case_when(
      is.na(TPM)|TPM==0~"Low", TPM<=tpm_medium~"Low",
      TPM<=tpm_high~"Medium", TRUE~"High"))
message("Expression bins:"); print(table(deseq_res$expr_bin))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Load eigenvectors from bedgraph
# ══════════════════════════════════════════════════════════════════════════════
load_eigenvectors <- function(comp_bed_path, sample_label) {
  if (is.null(comp_bed_path)||!file.exists(comp_bed_path)) {
    message("  No bedgraph for ",sample_label," — skipping"); return(NULL)
  }
  df <- read.table(comp_bed_path, header=FALSE, sep="\t",
                   stringsAsFactors=FALSE,
                   col.names=c("chr","start","end","eigen")) %>%
    dplyr::mutate(chr=as.character(chr)) %>%
    filter(!is.na(eigen))
  message("  ",sample_label,": ",nrow(df)," bins across ",
          length(unique(df$chr))," chromosomes")
  df
}

check_eigen_sign <- function(eigen_df, deseq_res, sample_name) {
  gene_coords <- deseq_res %>%
    filter(!is.na(gene_seq_start),!is.na(gene_seq_end),
           !is.na(seq_name),expr_bin=="High") %>%
    distinct(gene_id,.keep_all=TRUE)
  gene_gr <- GRanges(seqnames=paste0("chr",gene_coords$seq_name),
                     ranges=IRanges(gene_coords$gene_seq_start,
                                    gene_coords$gene_seq_end))
  bin_gr  <- GRanges(seqnames=paste0("chr",eigen_df$chr),
                     ranges=IRanges(eigen_df$start+1, eigen_df$end))
  hits      <- findOverlaps(bin_gr, gene_gr)
  high_bins <- unique(queryHits(hits))
  mean_high <- mean(eigen_df$eigen[high_bins], na.rm=TRUE)
  pct_pos   <- mean(eigen_df$eigen[high_bins]>0, na.rm=TRUE)*100
  message("  ",sample_name,": mean eigen (High-expr bins) = ",round(mean_high,4),
          "  (",round(pct_pos,1),"% positive)")
  if (mean_high<0) {
    message("  *** Sign flipped — correcting ***")
    eigen_df <- eigen_df %>% dplyr::mutate(eigen=-eigen)
  }
  eigen_df
}

message("\nLoading eigenvectors...")
eigen_list <- lapply(names(CONDITIONS), function(n) {
  cond <- CONDITIONS[[n]]
  df   <- load_eigenvectors(cond$comp_bed, cond$label)
  if (!is.null(df)) df %>% dplyr::mutate(cond_label=cond$label) else NULL
})
names(eigen_list) <- names(CONDITIONS)

message("\nRaw eigenvector sign diagnostic:")
invisible(lapply(names(eigen_list), function(n) {
  if (is.null(eigen_list[[n]])) return(NULL)
  ef <- eigen_list[[n]]
  gene_coords <- deseq_res %>%
    filter(!is.na(gene_seq_start),!is.na(gene_seq_end),
           !is.na(seq_name),expr_bin=="High") %>%
    distinct(gene_id,.keep_all=TRUE)
  gene_gr <- GRanges(seqnames=paste0("chr",gene_coords$seq_name),
                     ranges=IRanges(gene_coords$gene_seq_start,gene_coords$gene_seq_end))
  bin_gr  <- GRanges(seqnames=paste0("chr",ef$chr),
                     ranges=IRanges(ef$start+1,ef$end))
  hits      <- findOverlaps(bin_gr,gene_gr)
  high_bins <- unique(queryHits(hits))
  message(CONDITIONS[[n]]$label,":")
  message("  overall pct_positive: ",round(mean(ef$eigen>0,na.rm=TRUE)*100,1),"%")
  if (length(high_bins)>0) {
    message("  mean eigen in High-expr bins: ",
            round(mean(ef$eigen[high_bins],na.rm=TRUE),4),
            "  pct_positive=",
            round(mean(ef$eigen[high_bins]>0,na.rm=TRUE)*100,1),"%")
    message("  → sign ",
            ifelse(mean(ef$eigen[high_bins],na.rm=TRUE)<0,"NEEDS FLIPPING","is correct"))
  } else {
    message("  *** WARNING: no overlaps found ***")
  }
}))

message("\nApplying sign correction...")
eigen_list_signed <- lapply(names(eigen_list), function(n) {
  if (is.null(eigen_list[[n]])) return(NULL)
  check_eigen_sign(eigen_list[[n]], deseq_res, CONDITIONS[[n]]$label)
})
names(eigen_list_signed) <- names(eigen_list)
eigen_list <- Filter(Negate(is.null), eigen_list_signed)
message("Conditions with eigenvectors: ", paste(names(eigen_list), collapse=", "))

# ══════════════════════════════════════════════════════════════════════════════
# 3. Coding/noncoding ratio per 500kb bin
# ══════════════════════════════════════════════════════════════════════════════
compute_bin_coding_ratio <- function(resolution=RESOLUTION) {
  message("Loading genome annotation from TxDb...")
  txdb       <- TxDb.Hsapiens.UCSC.hg38.knownGene
  exons_gr   <- exons(txdb)  %>% GenomicRanges::reduce()
  genes_gr   <- genes(txdb)  %>% GenomicRanges::reduce()
  introns_gr <- GenomicRanges::setdiff(genes_gr, exons_gr)
  message("  Exon intervals: ",length(exons_gr),
          "  Intron intervals: ",length(introns_gr))
  chroms  <- paste0("chr",c(as.character(1:22),"X"))
  seqlens <- seqlengths(txdb)[chroms]
  seqlens <- seqlens[!is.na(seqlens)]
  message("  Tiling genome at ",resolution/1e3,"kb...")
  bin_gr <- tileGenome(seqlengths=seqlens, tilewidth=resolution,
                       cut.last.tile.in.chrom=TRUE)
  compute_overlap_bp <- function(bin_gr, feature_gr) {
    hits       <- findOverlaps(bin_gr, feature_gr)
    overlaps   <- GenomicRanges::pintersect(bin_gr[queryHits(hits)],
                                            feature_gr[subjectHits(hits)])
    bp_per_bin <- tapply(width(overlaps), queryHits(hits), sum)
    result     <- rep(0L, length(bin_gr))
    result[as.integer(names(bp_per_bin))] <- as.integer(bp_per_bin)
    result
  }
  message("  Computing exon coverage...")
  exon_bp   <- compute_overlap_bp(bin_gr, exons_gr)
  message("  Computing intron coverage...")
  intron_bp <- compute_overlap_bp(bin_gr, introns_gr)
  data.frame(
    chr       = gsub("chr","",as.character(seqnames(bin_gr))),
    start     = start(bin_gr)-1,
    end       = end(bin_gr),
    bin_width = width(bin_gr),
    exon_bp   = exon_bp,
    intron_bp = intron_bp,
    stringsAsFactors=FALSE
  ) %>%
    dplyr::mutate(
      genic_bp      = exon_bp+intron_bp,
      intergenic_bp = pmax(bin_width-genic_bp,0),
      noncoding_bp  = intron_bp+intergenic_bp,
      coding_ratio  = exon_bp/bin_width,
      ei_ratio      = ifelse(genic_bp>0, exon_bp/genic_bp, NA_real_)
    )
}

coding_rds <- file.path(RDS_DIR, "coding_bins_500kb.rds")
if (file.exists(coding_rds)) {
  message("Loading saved coding bins...")
  coding_bins <- readRDS(coding_rds)
} else {
  message("\n══ Computing coding/noncoding ratios ════════════════════════════")
  coding_bins <- compute_bin_coding_ratio(RESOLUTION)
  saveRDS(coding_bins, coding_rds)
  message("Saved: ",coding_rds)
}

coding_expressed   <- coding_bins %>% filter(coding_ratio>0) %>% pull(coding_ratio)
coding_low_thresh  <- quantile(coding_expressed,1/3,na.rm=TRUE)
coding_high_thresh <- quantile(coding_expressed,2/3,na.rm=TRUE)
coding_bins <- coding_bins %>%
  dplyr::mutate(
    coding_bin = dplyr::case_when(
      coding_ratio==0                  ~ "Low",
      coding_ratio<=coding_low_thresh  ~ "Low",
      coding_ratio<=coding_high_thresh ~ "Medium",
      TRUE                             ~ "High"
    ),
    coding_bin = factor(coding_bin, levels=c("Low","Medium","High"))
  )
message("Coding bin counts:"); print(table(coding_bins$coding_bin))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Merge coding annotation into eigenvectors
# ══════════════════════════════════════════════════════════════════════════════
merge_coding_to_eigen <- function(eigen_df, coding_bins, sample_label) {
  coding_sub <- coding_bins %>%
    dplyr::select(chr, start, coding_ratio, ei_ratio,
                  exon_bp, intron_bp, noncoding_bp, coding_bin)
  result <- eigen_df %>%
    left_join(coding_sub, by=c("chr","start")) %>%
    dplyr::mutate(
      coding_ratio = replace_na(coding_ratio, 0),
      exon_bp      = replace_na(exon_bp,      0),
      intron_bp    = replace_na(intron_bp,     0),
      noncoding_bp = replace_na(noncoding_bp,  0),
      coding_bin   = replace_na(as.character(coding_bin),"Low"),
      coding_bin   = factor(coding_bin, levels=c("Low","Medium","High"))
    )
  message("  ",sample_label," coding bins:"); print(table(result$coding_bin))
  result
}

message("\nMerging coding annotation...")
eigen_coding_list <- lapply(names(eigen_list), function(n)
  merge_coding_to_eigen(eigen_list[[n]], coding_bins, CONDITIONS[[n]]$label))
names(eigen_coding_list) <- names(eigen_list)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Assign expression bins to 500kb genomic bins
# ══════════════════════════════════════════════════════════════════════════════
assign_expr_bins_to_genome <- function(eigen_df, deseq_res) {
  gene_coords <- deseq_res %>%
    filter(!is.na(gene_seq_start),!is.na(gene_seq_end),!is.na(seq_name)) %>%
    distinct(gene_id,.keep_all=TRUE)
  gene_gr <- GRanges(
    seqnames = paste0("chr",gene_coords$seq_name),
    ranges   = IRanges(gene_coords$gene_seq_start,gene_coords$gene_seq_end),
    TPM      = gene_coords$TPM
  )
  bin_gr <- GRanges(
    seqnames = paste0("chr",eigen_df$chr),
    ranges   = IRanges(eigen_df$start+1, eigen_df$end)
  )
  hits <- findOverlaps(bin_gr, gene_gr)
  bin_tpm <- data.frame(
    bin_idx = queryHits(hits),
    TPM     = gene_gr$TPM[subjectHits(hits)]
  ) %>%
    dplyr::group_by(bin_idx) %>%
    dplyr::summarise(mean_TPM=mean(TPM,na.rm=TRUE), .groups="drop")
  bin_expressed  <- bin_tpm %>% filter(!is.na(mean_TPM),mean_TPM>0) %>% pull(mean_TPM)
  bin_tpm_medium <- quantile(bin_expressed,1/3,na.rm=TRUE)
  bin_tpm_high   <- quantile(bin_expressed,2/3,na.rm=TRUE)
  message("  Bin-level TPM thresholds: Low<=",round(bin_tpm_medium,2),
          "  Med<=",round(bin_tpm_high,2),
          "  High>",round(bin_tpm_high,2))
  eigen_df %>%
    dplyr::mutate(bin_idx=row_number()) %>%
    left_join(bin_tpm, by="bin_idx") %>%
    dplyr::mutate(
      expr_bin = dplyr::case_when(
        is.na(mean_TPM)|mean_TPM==0 ~ "Low",
        mean_TPM<=bin_tpm_medium    ~ "Low",
        mean_TPM<=bin_tpm_high      ~ "Medium",
        TRUE                        ~ "High"
      ),
      expr_bin = factor(expr_bin, levels=BINS)
    ) %>%
    dplyr::select(-bin_idx)
}

message("\nAssigning expression bins to 500kb bins...")
eigen_expr_list <- lapply(names(eigen_list), function(n) {
  message(CONDITIONS[[n]]$label,":")
  assign_expr_bins_to_genome(eigen_list[[n]], deseq_res)
})
names(eigen_expr_list) <- names(eigen_list)
invisible(lapply(names(eigen_expr_list), function(n) {
  message(CONDITIONS[[n]]$label," expr bin counts:")
  print(table(eigen_expr_list[[n]]$expr_bin))
}))

# ══════════════════════════════════════════════════════════════════════════════
# 6. Saddle computation
# ══════════════════════════════════════════════════════════════════════════════
compute_saddle_chr <- function(hic_file, chrom, eigen_chr,
                               n_bins=SADDLE_BINS, resolution=RESOLUTION,
                               norm=STRAW_NORM) {
  library(strawr); library(dplyr)
  eigen_valid <- eigen_chr %>%
    filter(!is.na(eigen)) %>%
    dplyr::arrange(eigen) %>%
    dplyr::mutate(quantile_bin=ntile(eigen, n_bins))
  if (nrow(eigen_valid)<n_bins) return(NULL)
  raw <- tryCatch(
    strawr::straw(norm, hic_file, as.character(chrom),
                  as.character(chrom), "BP", resolution),
    error=function(e) NULL)
  if (is.null(raw)||nrow(raw)==0) return(NULL)
  raw <- raw %>% dplyr::mutate(bin1=x/resolution, bin2=y/resolution)
  bin_to_quantile <- setNames(eigen_valid$quantile_bin,
                              as.character(eigen_valid$start/resolution))
  raw <- raw %>%
    dplyr::mutate(
      q1        = bin_to_quantile[as.character(bin1)],
      q2        = bin_to_quantile[as.character(bin2)],
      dist_bins = abs(bin2-bin1)
    ) %>% filter(!is.na(q1),!is.na(q2))
  if (nrow(raw)==0) return(NULL)
  expected_by_dist <- raw %>%
    dplyr::group_by(dist_bins) %>%
    dplyr::summarise(expected=mean(counts,na.rm=TRUE),.groups="drop")
  raw <- raw %>%
    left_join(expected_by_dist,by="dist_bins") %>%
    dplyr::mutate(oe=counts/expected) %>%
    filter(!is.na(oe),is.finite(oe))
  if (nrow(raw)==0) return(NULL)
  saddle_upper <- raw %>%
    dplyr::group_by(q1,q2) %>%
    dplyr::summarise(mean_oe=mean(oe,na.rm=TRUE),n=dplyr::n(),.groups="drop")
  bind_rows(saddle_upper,
            saddle_upper %>% dplyr::rename(q1=q2,q2=q1)) %>%
    distinct(q1,q2,.keep_all=TRUE) %>%
    dplyr::mutate(chr=as.character(chrom))
}

compute_stratified_saddle <- function(hic_file, eigen_annotated,
                                      sample_name, filter_col, filter_val,
                                      n_bins=SADDLE_BINS,
                                      resolution=RESOLUTION,
                                      norm=STRAW_NORM,
                                      n_cores=N_CORES) {
  eigen_subset <- eigen_annotated %>%
    filter(.data[[filter_col]]==filter_val)
  message("  ",filter_val,": ",nrow(eigen_subset)," bins")
  if (nrow(eigen_subset)<n_bins*2) { message("  Too few — skipping"); return(NULL) }
  chroms <- unique(eigen_subset$chr)
  results_list <- parallel::mclapply(chroms, function(chrom) {
    eigen_chr <- eigen_subset %>% filter(chr==chrom)
    if (nrow(eigen_chr)<n_bins) return(NULL)
    compute_saddle_chr(hic_file, chrom, eigen_chr, n_bins, resolution, norm)
  }, mc.cores=n_cores, mc.preschedule=FALSE)
  result <- bind_rows(Filter(Negate(is.null),results_list))
  if (nrow(result)==0) return(NULL)
  result %>%
    dplyr::group_by(q1,q2) %>%
    dplyr::summarise(mean_oe=mean(mean_oe,na.rm=TRUE),
                     n=sum(n),.groups="drop") %>%
    dplyr::mutate(
      cond_label = sample_name,
      log2_oe    = log2(mean_oe),
      filter_col = filter_col,
      filter_val = filter_val
    )
}

# ── Coding-stratified saddle ──────────────────────────────────────────────────
coding_saddle_rds <- file.path(RDS_DIR, "saddle_coding_ratio_500kb_SCALE.rds")
if (file.exists(coding_saddle_rds)) {
  message("Loading saved coding saddle matrices...")
  saddle_coding <- readRDS(coding_saddle_rds)
} else {
  message("\n══ Computing coding-stratified saddle matrices ══════════════════")
  coding_results <- list()
  for (n in names(eigen_coding_list)) {
    hic_file <- CONDITIONS[[n]]$hic; samp <- CONDITIONS[[n]]$label
    message("\n--- ",samp," ---")
    for (cb in c("Low","Medium","High")) {
      res <- compute_stratified_saddle(hic_file, eigen_coding_list[[n]],
                                       samp, "coding_bin", cb, norm=STRAW_NORM)
      if (!is.null(res)) coding_results[[paste(n,cb,sep="_")]] <- res
    }
  }
  saddle_coding <- bind_rows(coding_results) %>%
    dplyr::rename(sample=cond_label) %>%
    dplyr::mutate(sample=factor(.data$sample,levels=COND_LEVELS))
  saveRDS(saddle_coding, coding_saddle_rds)
  message("Saved: ",coding_saddle_rds)
}

# ── Expression-stratified saddle ──────────────────────────────────────────────
expr_saddle_rds <- file.path(RDS_DIR, "saddle_expr_bin_500kb_SCALE.rds")
if (file.exists(expr_saddle_rds)) {
  message("Loading saved expression saddle matrices...")
  saddle_expr <- readRDS(expr_saddle_rds)
} else {
  message("\n══ Computing expression-stratified saddle matrices ══════════════")
  expr_results <- list()
  for (n in names(eigen_expr_list)) {
    hic_file <- CONDITIONS[[n]]$hic; samp <- CONDITIONS[[n]]$label
    message("\n--- ",samp," ---")
    for (eb in c("Low","Medium","High")) {
      res <- compute_stratified_saddle(hic_file, eigen_expr_list[[n]],
                                       samp, "expr_bin", eb, norm=STRAW_NORM)
      if (!is.null(res)) expr_results[[paste(n,eb,sep="_")]] <- res
    }
  }
  saddle_expr <- bind_rows(expr_results) %>%
    dplyr::rename(sample=cond_label) %>%
    dplyr::mutate(sample=factor(.data$sample,levels=COND_LEVELS))
  saveRDS(saddle_expr, expr_saddle_rds)
  message("Saved: ",expr_saddle_rds)
}

message("\nSaddle coding summary:")
print(table(as.character(saddle_coding$sample), saddle_coding$filter_val))
message("\nSaddle expression summary:")
print(table(as.character(saddle_expr$sample), saddle_expr$filter_val))

# ══════════════════════════════════════════════════════════════════════════════
# 7. Core helpers
# ══════════════════════════════════════════════════════════════════════════════
compartment_strength <- function(saddle_df, n_bins=SADDLE_BINS, corner_size=5) {
  lo <- 1:corner_size; hi <- (n_bins-corner_size+1):n_bins
  get_mean <- function(q1r,q2r)
    saddle_df %>% filter(q1%in%q1r,q2%in%q2r) %>%
    summarise(m=mean(mean_oe,na.rm=TRUE)) %>% pull(m)
  aa <- get_mean(hi,hi); bb <- get_mean(lo,lo)
  ab <- get_mean(lo,hi); ba <- get_mean(hi,lo)
  list(AA=aa,BB=bb,AB=ab,BA=ba,CS=(aa+bb)/(ab+ba))
}

make_corner_scores <- function(cs_vals, n_bins) {
  data.frame(
    q1   =c(1.5,       1.5,        n_bins-0.5, n_bins-0.5),
    q2   =c(n_bins-0.5,1.5,        n_bins-0.5, 1.5),
    label=c(paste0("BA\n",round(log2(cs_vals$BA),2)),
            paste0("BB\n",round(log2(cs_vals$BB),2)),
            paste0("AA\n",round(log2(cs_vals$AA),2)),
            paste0("AB\n",round(log2(cs_vals$AB),2))),
    hjust=c(0,0,1,1), vjust=c(1,0,1,0)
  )
}

# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# 8. Shared axis helpers
# ══════════════════════════════════════════════════════════════════════════════
get_shared_ei_lim <- function(eigen_coding_df, coding_bins=c("Low","Medium","High")) {
  all_vals <- unlist(lapply(coding_bins, function(cb) {
    eigen_coding_df %>%
      dplyr::filter(as.character(coding_bin)==cb, !is.na(eigen),
                    exon_bp>0, intron_bp>0) %>%
      dplyr::mutate(log2_ei=log2(exon_bp/intron_bp)) %>%
      pull(log2_ei)
  }))
  all_vals <- all_vals[is.finite(all_vals)]
  quantile(abs(all_vals), 0.98, na.rm=TRUE)*1.1
}

get_shared_eigen_lim <- function(eigen_df, filter_col,
                                 bin_vals=c("Low","Medium","High"),
                                 n_bins=SADDLE_BINS) {
  all_means <- unlist(lapply(bin_vals, function(b) {
    sub <- eigen_df %>% dplyr::filter(.data[[filter_col]]==b, !is.na(eigen))
    if (nrow(sub)==0) return(NULL)
    sub %>%
      dplyr::mutate(quantile_bin=ntile(eigen, n_bins)) %>%
      dplyr::group_by(quantile_bin) %>%
      dplyr::summarise(mean_eigen=mean(eigen,na.rm=TRUE),.groups="drop") %>%
      pull(mean_eigen)
  }))
  max(abs(all_means), na.rm=TRUE)*1.1
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Series 1 — E/I boxplot marginals
# ══════════════════════════════════════════════════════════════════════════════
plot_saddle_ei <- function(saddle_df, eigen_coding_df,
                           sample_name, coding_bin_val,
                           n_bins=SADDLE_BINS,
                           title=NULL, color_limit=NULL,
                           ei_lim=NULL) {
  mat_df <- saddle_df %>%
    dplyr::filter(as.character(.data$sample)==sample_name,
                  as.character(filter_val)==coding_bin_val) %>%
    tidyr::complete(q1=1:n_bins, q2=1:n_bins, fill=list(log2_oe=NA))
  if (nrow(mat_df)==0||all(is.na(mat_df$log2_oe))) return(NULL)
  if (is.null(color_limit))
    color_limit <- quantile(abs(mat_df$log2_oe),0.95,na.rm=TRUE)
  saddle_sub    <- saddle_df %>%
    dplyr::filter(as.character(.data$sample)==sample_name,
                  as.character(filter_val)==coding_bin_val)
  cs_vals       <- compartment_strength(saddle_sub)
  corner_scores <- make_corner_scores(cs_vals, n_bins)
  eigen_sub <- eigen_coding_df %>%
    dplyr::filter(as.character(coding_bin)==coding_bin_val, !is.na(eigen)) %>%
    dplyr::mutate(quantile_bin=as.integer(ntile(eigen, n_bins)))
  if (nrow(eigen_sub)==0) return(NULL)
  # Per-quantile aggregate E/I — determines box color
  ei_per_bin <- eigen_sub %>%
    dplyr::group_by(quantile_bin) %>%
    dplyr::summarise(
      total_exon_bp   = sum(exon_bp,   na.rm=TRUE),
      total_intron_bp = sum(intron_bp, na.rm=TRUE),
      .groups="drop") %>%
    dplyr::mutate(
      log2_ei_agg = log2(total_exon_bp/total_intron_bp),
      fill_color  = ifelse(log2_ei_agg>=0, COL_A, COL_B)
    ) %>%
    filter(is.finite(log2_ei_agg))
  # Per-bin values for boxplot — color via I() to bypass scale system
  ei_boxplot <- eigen_sub %>%
    dplyr::filter(exon_bp>0, intron_bp>0) %>%
    dplyr::mutate(log2_ei=log2(exon_bp/intron_bp)) %>%
    filter(is.finite(log2_ei)) %>%
    left_join(ei_per_bin %>% dplyr::select(quantile_bin, fill_color),
              by="quantile_bin")
  if (is.null(ei_lim))
    ei_lim <- quantile(abs(ei_boxplot$log2_ei),0.98,na.rm=TRUE)*1.1
  ei_breaks      <- round(c(-ei_lim*0.8,0,ei_lim*0.8),2)
  subtitle_label <- paste0(nrow(eigen_sub)," bins  |  coding = ",coding_bin_val)
  ggplot(mat_df,aes(x=q1,y=q2))+
    geom_tile(aes(fill=pmin(pmax(log2_oe,-color_limit),color_limit)))+
    scale_fill_distiller(palette="RdBu",direction=-1,
                         limits=c(-color_limit,color_limit),
                         name="log2 O/E",na.value="gray80")+
    geom_text(data=corner_scores,
              aes(x=q1,y=q2,label=label,hjust=hjust,vjust=vjust),
              color="black",size=2.8,fontface="bold",
              lineheight=0.85,inherit.aes=FALSE)+
    annotate("text",x=n_bins/2,y=n_bins/2,
             label=paste0("CS=",round(cs_vals$CS,2)),
             size=3.2,fontface="bold",color="gray20")+
    ggnewscale::new_scale_fill()+
    ggside::geom_xsideboxplot(
      data=ei_boxplot,
      aes(x=quantile_bin, y=log2_ei, group=quantile_bin, fill=I(fill_color)),
      inherit.aes=FALSE, outlier.size=0.3, outlier.alpha=0.3, linewidth=0.3)+
    ggside::geom_xsidehline(yintercept=0,color="gray30",linewidth=0.5,linetype="dashed")+
    scale_xsidey_continuous(name="log2(E/I)",limits=c(-ei_lim,ei_lim),
                            breaks=ei_breaks,
                            labels=scales::number(ei_breaks,accuracy=0.01))+
    ggside::geom_ysideboxplot(
      data=ei_boxplot,
      aes(y=quantile_bin, x=log2_ei, group=quantile_bin, fill=I(fill_color)),
      inherit.aes=FALSE, outlier.size=0.3, outlier.alpha=0.3, linewidth=0.3)+
    ggside::geom_ysidehline(yintercept=0,color="gray30",linewidth=0.5,linetype="dashed")+
    scale_ysidex_continuous(name="log2(E/I)",limits=c(-ei_lim,ei_lim),
                            breaks=ei_breaks,
                            labels=scales::number(ei_breaks,accuracy=0.01))+
    scale_x_continuous(name="Eigenvector quantile",
                       breaks=c(1,n_bins/2,n_bins),
                       labels=c("B","0","A"),expand=c(0,0))+
    scale_y_continuous(name="Eigenvector quantile",
                       breaks=c(1,n_bins/2,n_bins),
                       labels=c("B","0","A"),expand=c(0,0))+
    coord_fixed()+labs(title=title,subtitle=subtitle_label)+
    theme_classic(base_size=11)+ggside::theme_ggside_void()+
    theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
          plot.subtitle=element_text(size=7,color="gray40",hjust=0.5),
          legend.position="none",
          axis.text=element_text(size=9,color="black"),
          axis.title=element_text(size=10,color="black",face="bold"),
          axis.ticks=element_line(color="black",linewidth=0.4),
          plot.margin=margin(4,4,4,4),
          ggside.panel.scale.x=0.3,ggside.panel.scale.y=0.3,
          ggside.axis.text.x=element_text(size=7,color="black"),
          ggside.axis.text.y=element_text(size=7,color="black"),
          ggside.axis.ticks.x=element_line(color="black",linewidth=0.3),
          ggside.axis.ticks.y=element_line(color="black",linewidth=0.3))
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. Series 2 — PC1 bar marginals (coding bins)
# ══════════════════════════════════════════════════════════════════════════════
plot_saddle_pc1 <- function(saddle_df, eigen_coding_df,
                            sample_name, coding_bin_val,
                            n_bins=SADDLE_BINS,
                            title=NULL, color_limit=NULL,
                            eigen_lim=NULL) {
  mat_df <- saddle_df %>%
    dplyr::filter(as.character(.data$sample)==sample_name,
                  as.character(filter_val)==coding_bin_val) %>%
    tidyr::complete(q1=1:n_bins, q2=1:n_bins, fill=list(log2_oe=NA))
  if (nrow(mat_df)==0||all(is.na(mat_df$log2_oe))) return(NULL)
  if (is.null(color_limit))
    color_limit <- quantile(abs(mat_df$log2_oe),0.95,na.rm=TRUE)
  saddle_sub    <- saddle_df %>%
    dplyr::filter(as.character(.data$sample)==sample_name,
                  as.character(filter_val)==coding_bin_val)
  cs_vals       <- compartment_strength(saddle_sub)
  corner_scores <- make_corner_scores(cs_vals, n_bins)
  eigen_sub <- eigen_coding_df %>%
    dplyr::filter(as.character(coding_bin)==coding_bin_val, !is.na(eigen)) %>%
    dplyr::mutate(quantile_bin=as.numeric(ntile(eigen, n_bins)))
  if (nrow(eigen_sub)==0) return(NULL)
  pc1_marginals <- eigen_sub %>%
    dplyr::group_by(quantile_bin) %>%
    dplyr::summarise(
      mean_eigen = mean(eigen,na.rm=TRUE),
      .groups="drop") %>%
    dplyr::mutate(
      fill_color = ifelse(mean_eigen>=0, COL_A, COL_B)
    )
  if (is.null(eigen_lim))
    eigen_lim <- max(abs(pc1_marginals$mean_eigen),na.rm=TRUE)*1.1
  eigen_breaks   <- round(c(-eigen_lim*0.8,0,eigen_lim*0.8),2)
  subtitle_label <- paste0(nrow(eigen_sub)," bins  |  coding = ",coding_bin_val)
  ggplot(mat_df,aes(x=q1,y=q2))+
    geom_tile(aes(fill=pmin(pmax(log2_oe,-color_limit),color_limit)))+
    scale_fill_distiller(palette="RdBu",direction=-1,
                         limits=c(-color_limit,color_limit),
                         name="log2 O/E",na.value="gray80")+
    geom_text(data=corner_scores,
              aes(x=q1,y=q2,label=label,hjust=hjust,vjust=vjust),
              color="black",size=2.8,fontface="bold",
              lineheight=0.85,inherit.aes=FALSE)+
    annotate("text",x=n_bins/2,y=n_bins/2,
             label=paste0("CS=",round(cs_vals$CS,2)),
             size=3.2,fontface="bold",color="gray20")+
    ggnewscale::new_scale_fill()+
    ggside::geom_xsidecol(data=pc1_marginals,
                          aes(x=quantile_bin,y=mean_eigen,fill=I(fill_color)),
                          inherit.aes=FALSE,width=0.9)+
    ggside::geom_xsidehline(yintercept=0,color="gray30",linewidth=0.4)+
    scale_xsidey_continuous(name="PC1",limits=c(-eigen_lim,eigen_lim),
                            breaks=eigen_breaks,
                            labels=scales::number(eigen_breaks,accuracy=0.01))+
    ggside::geom_ysidecol(data=pc1_marginals,
                          aes(y=quantile_bin,x=mean_eigen,fill=I(fill_color)),
                          inherit.aes=FALSE,width=0.9)+
    ggside::geom_ysidehline(yintercept=0,color="gray30",linewidth=0.4)+
    scale_ysidex_continuous(name="PC1",limits=c(-eigen_lim,eigen_lim),
                            breaks=eigen_breaks,
                            labels=scales::number(eigen_breaks,accuracy=0.01))+
    scale_x_continuous(name="Eigenvector quantile",
                       breaks=c(1,n_bins/2,n_bins),
                       labels=c("B","0","A"),expand=c(0,0))+
    scale_y_continuous(name="Eigenvector quantile",
                       breaks=c(1,n_bins/2,n_bins),
                       labels=c("B","0","A"),expand=c(0,0))+
    coord_fixed()+labs(title=title,subtitle=subtitle_label)+
    theme_classic(base_size=11)+ggside::theme_ggside_void()+
    theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
          plot.subtitle=element_text(size=7,color="gray40",hjust=0.5),
          legend.position="none",
          axis.text=element_text(size=9,color="black"),
          axis.title=element_text(size=10,color="black",face="bold"),
          axis.ticks=element_line(color="black",linewidth=0.4),
          plot.margin=margin(4,4,4,4),
          ggside.panel.scale.x=0.25,ggside.panel.scale.y=0.25,
          ggside.axis.text.x=element_text(size=7,color="black"),
          ggside.axis.text.y=element_text(size=7,color="black"),
          ggside.axis.ticks.x=element_line(color="black",linewidth=0.3),
          ggside.axis.ticks.y=element_line(color="black",linewidth=0.3))
}

# ══════════════════════════════════════════════════════════════════════════════
# 11. Series 3 — Expression-stratified saddle with PC1 marginals
# ══════════════════════════════════════════════════════════════════════════════
plot_saddle_expr_pc1 <- function(saddle_df, eigen_expr_df,
                                 sample_name, expr_bin_val,
                                 n_bins=SADDLE_BINS,
                                 title=NULL, color_limit=NULL,
                                 eigen_lim=NULL) {
  mat_df <- saddle_df %>%
    dplyr::filter(as.character(.data$sample)==sample_name,
                  as.character(filter_val)==expr_bin_val) %>%
    tidyr::complete(q1=1:n_bins, q2=1:n_bins, fill=list(log2_oe=NA))
  if (nrow(mat_df)==0||all(is.na(mat_df$log2_oe))) return(NULL)
  if (is.null(color_limit))
    color_limit <- quantile(abs(mat_df$log2_oe),0.95,na.rm=TRUE)
  saddle_sub    <- saddle_df %>%
    dplyr::filter(as.character(.data$sample)==sample_name,
                  as.character(filter_val)==expr_bin_val)
  cs_vals       <- compartment_strength(saddle_sub)
  corner_scores <- make_corner_scores(cs_vals, n_bins)
  eigen_sub <- eigen_expr_df %>%
    dplyr::filter(as.character(expr_bin)==expr_bin_val, !is.na(eigen)) %>%
    dplyr::mutate(quantile_bin=as.numeric(ntile(eigen, n_bins)))
  if (nrow(eigen_sub)==0) return(NULL)
  pc1_marginals <- eigen_sub %>%
    dplyr::group_by(quantile_bin) %>%
    dplyr::summarise(
      mean_eigen = mean(eigen,na.rm=TRUE),
      .groups="drop") %>%
    dplyr::mutate(
      fill_color = ifelse(mean_eigen>=0, COL_A, COL_B)
    )
  if (is.null(eigen_lim))
    eigen_lim <- max(abs(pc1_marginals$mean_eigen),na.rm=TRUE)*1.1
  eigen_breaks   <- round(c(-eigen_lim*0.8,0,eigen_lim*0.8),2)
  subtitle_label <- paste0(nrow(eigen_sub)," bins  |  expression = ",expr_bin_val)
  ggplot(mat_df,aes(x=q1,y=q2))+
    geom_tile(aes(fill=pmin(pmax(log2_oe,-color_limit),color_limit)))+
    scale_fill_distiller(palette="RdBu",direction=-1,
                         limits=c(-color_limit,color_limit),
                         name="log2 O/E",na.value="gray80")+
    geom_text(data=corner_scores,
              aes(x=q1,y=q2,label=label,hjust=hjust,vjust=vjust),
              color="black",size=2.8,fontface="bold",
              lineheight=0.85,inherit.aes=FALSE)+
    annotate("text",x=n_bins/2,y=n_bins/2,
             label=paste0("CS=",round(cs_vals$CS,2)),
             size=3.2,fontface="bold",color="gray20")+
    ggnewscale::new_scale_fill()+
    ggside::geom_xsidecol(data=pc1_marginals,
                          aes(x=quantile_bin,y=mean_eigen,fill=I(fill_color)),
                          inherit.aes=FALSE,width=0.9)+
    ggside::geom_xsidehline(yintercept=0,color="gray30",linewidth=0.4)+
    scale_xsidey_continuous(name="PC1",limits=c(-eigen_lim,eigen_lim),
                            breaks=eigen_breaks,
                            labels=scales::number(eigen_breaks,accuracy=0.01))+
    ggside::geom_ysidecol(data=pc1_marginals,
                          aes(y=quantile_bin,x=mean_eigen,fill=I(fill_color)),
                          inherit.aes=FALSE,width=0.9)+
    ggside::geom_ysidehline(yintercept=0,color="gray30",linewidth=0.4)+
    scale_ysidex_continuous(name="PC1",limits=c(-eigen_lim,eigen_lim),
                            breaks=eigen_breaks,
                            labels=scales::number(eigen_breaks,accuracy=0.01))+
    scale_x_continuous(name="Eigenvector quantile",
                       breaks=c(1,n_bins/2,n_bins),
                       labels=c("B","0","A"),expand=c(0,0))+
    scale_y_continuous(name="Eigenvector quantile",
                       breaks=c(1,n_bins/2,n_bins),
                       labels=c("B","0","A"),expand=c(0,0))+
    coord_fixed()+labs(title=title,subtitle=subtitle_label)+
    theme_classic(base_size=11)+ggside::theme_ggside_void()+
    theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
          plot.subtitle=element_text(size=7,color="gray40",hjust=0.5),
          legend.position="none",
          axis.text=element_text(size=9,color="black"),
          axis.title=element_text(size=10,color="black",face="bold"),
          axis.ticks=element_line(color="black",linewidth=0.4),
          plot.margin=margin(4,4,4,4),
          ggside.panel.scale.x=0.25,ggside.panel.scale.y=0.25,
          ggside.axis.text.x=element_text(size=7,color="black"),
          ggside.axis.text.y=element_text(size=7,color="black"),
          ggside.axis.ticks.x=element_line(color="black",linewidth=0.3),
          ggside.axis.ticks.y=element_line(color="black",linewidth=0.3))
}

# ══════════════════════════════════════════════════════════════════════════════
# 12. Legend helpers
# ══════════════════════════════════════════════════════════════════════════════
make_heatmap_legend <- function(color_limit) {
  cowplot::get_legend(
    ggplot(data.frame(x=1:100,y=1,fill=seq(-color_limit,color_limit,length.out=100)),
           aes(x=x,y=y,fill=fill))+geom_tile()+
      scale_fill_distiller(palette="RdBu",direction=-1,
                           limits=c(-color_limit,color_limit),name="log2 O/E")+
      theme_void()+
      theme(legend.position="right",legend.key.height=unit(1.5,"cm"),
            legend.key.width=unit(0.4,"cm"),legend.text=element_text(size=9),
            legend.title=element_text(size=10)))
}

make_ei_legend <- function() {
  cowplot::get_legend(
    ggplot(data.frame(
      direction=factor(c("Exon-enriched","Intron-enriched"),
                       levels=c("Intron-enriched","Exon-enriched")),
      x=c(1,1),y=c(1,2)),aes(x=x,y=y,fill=direction))+
      geom_col()+
      scale_fill_manual(values=c("Exon-enriched"=COL_A,"Intron-enriched"=COL_B),
                        name="log2(E/I)\nmarginal")+
      theme_void()+
      theme(legend.position="right",legend.key.size=unit(0.5,"cm"),
            legend.text=element_text(size=9),
            legend.title=element_text(size=10,face="bold")))
}

make_pc1_legend <- function() {
  cowplot::get_legend(
    ggplot(data.frame(compartment=factor(c("A","B"),levels=c("B","A")),
                      x=c(1,1),y=c(1,2)),aes(x=x,y=y,fill=compartment))+
      geom_col()+
      scale_fill_manual(values=c(A=COL_A,B=COL_B),name="PC1\nmarginal")+
      theme_void()+
      theme(legend.position="right",legend.key.size=unit(0.5,"cm"),
            legend.text=element_text(size=9),
            legend.title=element_text(size=10,face="bold")))
}

assemble_panel <- function(panels, legends, cond_color, title, subtitle) {
  panels <- Filter(Negate(is.null), panels)
  if (length(panels)==0) return(invisible(NULL))
  main <- cowplot::plot_grid(plotlist=panels, nrow=1,
                             labels=c("a","b","c"), label_size=12)
  composite <- cowplot::plot_grid(main, legends, ncol=2, rel_widths=c(1,0.1))
  title_grob <- cowplot::ggdraw()+
    cowplot::draw_label(title,fontface="bold",size=13,x=0.01,hjust=0,
                        color=cond_color)+
    cowplot::draw_label(subtitle,size=8,color="gray40",x=0.01,hjust=0,y=0.2)
  final <- cowplot::plot_grid(title_grob,composite,ncol=1,rel_heights=c(0.06,1))
  print(final); invisible(final)
}

# ══════════════════════════════════════════════════════════════════════════════
# 13. Run — Series 1: E/I boxplot marginals (Low/Med/High in one panel)
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Series 1: E/I boxplot marginals ═════════════════════════════════")
for (n in names(eigen_coding_list)) {
  samp        <- CONDITIONS[[n]]$label
  cond_color  <- CONDITIONS[[n]]$color
  saddle_sub  <- saddle_coding %>% dplyr::filter(as.character(.data$sample)==samp)
  color_limit <- quantile(abs(saddle_sub$log2_oe),0.95,na.rm=TRUE)
  ei_lim      <- get_shared_ei_lim(eigen_coding_list[[n]])
  message(samp," — shared ei_lim = ",round(ei_lim,3))
  panels <- lapply(c("Low","Medium","High"), function(cb) {
    plot_saddle_ei(saddle_coding, eigen_coding_list[[n]], samp, cb,
                   color_limit=color_limit, ei_lim=ei_lim,
                   title=paste0(cb," coding density"))
  })
  legs <- cowplot::plot_grid(make_heatmap_legend(color_limit),
                             make_ei_legend(), ncol=1)
  assemble_panel(panels, legs, cond_color,
                 title=paste0(samp,"  |  log2(E/I) marginals  |  ",DS_GROUP),
                 subtitle=paste0("500kb  |  ",STRAW_NORM,"  |  ",SADDLE_BINS,
                                 " quantile bins  |  boxplots per eigenvector quantile  |  ",
                                 "amber = exon-enriched  |  blue = intron-enriched"))
}

# ══════════════════════════════════════════════════════════════════════════════
# 14. Run — Series 2: PC1 bar marginals (coding bins)
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Series 2: PC1 bar marginals (coding bins) ═══════════════════════")
for (n in names(eigen_coding_list)) {
  samp        <- CONDITIONS[[n]]$label
  cond_color  <- CONDITIONS[[n]]$color
  saddle_sub  <- saddle_coding %>% dplyr::filter(as.character(.data$sample)==samp)
  color_limit <- quantile(abs(saddle_sub$log2_oe),0.95,na.rm=TRUE)
  eigen_lim   <- get_shared_eigen_lim(eigen_coding_list[[n]], "coding_bin")
  message(samp," — shared eigen_lim = ",round(eigen_lim,3))
  panels <- lapply(c("Low","Medium","High"), function(cb) {
    plot_saddle_pc1(saddle_coding, eigen_coding_list[[n]], samp, cb,
                    color_limit=color_limit, eigen_lim=eigen_lim,
                    title=paste0(cb," coding density"))
  })
  legs <- cowplot::plot_grid(make_heatmap_legend(color_limit),
                             make_pc1_legend(), ncol=1)
  assemble_panel(panels, legs, cond_color,
                 title=paste0(samp,"  |  PC1 marginals  |  coding bins  |  ",DS_GROUP),
                 subtitle=paste0("500kb  |  ",STRAW_NORM,"  |  ",SADDLE_BINS,
                                 " quantile bins  |  amber = A compartment  |  blue = B compartment"))
}

# ══════════════════════════════════════════════════════════════════════════════
# 15. Run — Series 3: Expression-stratified saddle with PC1 marginals
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Series 3: Expression-stratified saddle — PC1 marginals ══════════")
for (n in names(eigen_expr_list)) {
  samp        <- CONDITIONS[[n]]$label
  cond_color  <- CONDITIONS[[n]]$color
  saddle_sub  <- saddle_expr %>% dplyr::filter(as.character(.data$sample)==samp)
  color_limit <- quantile(abs(saddle_sub$log2_oe),0.95,na.rm=TRUE)
  eigen_lim   <- get_shared_eigen_lim(eigen_expr_list[[n]], "expr_bin")
  message(samp," — shared eigen_lim = ",round(eigen_lim,3))
  panels <- lapply(c("Low","Medium","High"), function(eb) {
    plot_saddle_expr_pc1(saddle_expr, eigen_expr_list[[n]], samp, eb,
                         color_limit=color_limit, eigen_lim=eigen_lim,
                         title=paste0(eb," expression"))
  })
  legs <- cowplot::plot_grid(make_heatmap_legend(color_limit),
                             make_pc1_legend(), ncol=1)
  assemble_panel(panels, legs, cond_color,
                 title=paste0(samp,"  |  PC1 marginals  |  expression bins  |  ",DS_GROUP),
                 subtitle=paste0("500kb  |  ",STRAW_NORM,"  |  ",SADDLE_BINS,
                                 " quantile bins  |  bin-level mean TPM tertiles  |  ",
                                 "amber = A compartment  |  blue = B compartment"))
}


# ══════════════════════════════════════════════════════════════════════════════
# 15. CS helper functions and summary plots
# ══════════════════════════════════════════════════════════════════════════════
extract_cs_table <- function(saddle_df, label,
                             n_bins=SADDLE_BINS, corner_size=5) {
  bind_rows(lapply(COND_LEVELS, function(samp) {
    cs <- compartment_strength(saddle_df %>% filter(.data$sample==samp),
                               n_bins=n_bins, corner_size=corner_size)
    data.frame(sample=samp, category=label,
               AA=cs$AA, BB=cs$BB, AB=cs$AB, BA=cs$BA, CS=cs$CS,
               log2_AA=log2(cs$AA), log2_BB=log2(cs$BB),
               log2_AB=log2(cs$AB), log2_BA=log2(cs$BA),
               stringsAsFactors=FALSE)
  }))
}

plot_cs_faceted <- function(cs_df, categories,
                            title="Compartment strength by category",
                            ncol_wrap=3) {
  plot_df <- cs_df %>%
    filter(category %in% categories) %>%
    dplyr::mutate(category=factor(category,levels=categories),
                  sample=factor(.data$sample,levels=COND_LEVELS))
  ggplot(plot_df,aes(x=sample,y=CS,fill=sample,color=sample))+
    geom_col(alpha=0.8,width=0.6)+
    geom_text(aes(label=round(CS,2)),vjust=-0.4,size=3,fontface="bold")+
    scale_fill_manual(values=COND_COLORS,  guide="none")+
    scale_color_manual(values=COND_COLORS, guide="none")+
    scale_y_continuous(name="Compartment strength\n(AA+BB)/(AB+BA)",
                       expand=expansion(mult=c(0,0.2)))+
    facet_wrap(~category,ncol=ncol_wrap,scales="free_y")+
    labs(title=title,x=NULL)+
    theme_classic(base_size=11)+
    theme(plot.title=element_text(face="bold",size=13,hjust=0.5),
          strip.text=element_text(face="bold",size=10),
          strip.background=element_blank(),
          axis.text.x=element_text(size=8,angle=30,hjust=1),
          plot.margin=margin(6,6,6,6))
}

plot_corner_scores <- function(cs_df, categories,
                               title="Per-corner log2 O/E scores") {
  plot_df <- cs_df %>%
    filter(category %in% categories) %>%
    dplyr::select(sample,category,log2_AA,log2_BB,log2_AB,log2_BA) %>%
    tidyr::pivot_longer(cols=c(log2_AA,log2_BB,log2_AB,log2_BA),
                        names_to="corner",values_to="log2_oe") %>%
    dplyr::mutate(
      corner   = gsub("log2_","",corner),
      corner   = factor(corner,levels=c("AA","BB","AB","BA")),
      category = factor(category,levels=categories),
      sample   = factor(.data$sample,levels=COND_LEVELS)
    )
  ggplot(plot_df,aes(x=corner,y=log2_oe,color=sample,group=sample))+
    geom_hline(yintercept=0,linetype="dashed",color="gray50",linewidth=0.4)+
    geom_line(linewidth=0.8)+geom_point(size=3)+
    geom_text(aes(label=round(log2_oe,2)),vjust=-0.8,size=2.5,show.legend=FALSE)+
    scale_color_manual(values=COND_COLORS,name=NULL)+
    scale_y_continuous(name="log2 O/E")+
    facet_wrap(~category,ncol=3,scales="free_y")+
    labs(title=title,x="Compartment interaction")+
    theme_classic(base_size=11)+
    theme(plot.title=element_text(face="bold",size=13,hjust=0.5),
          strip.text=element_text(face="bold",size=10),
          strip.background=element_blank(),
          legend.position="bottom",
          axis.text.x=element_text(size=9,face="bold"),
          plot.margin=margin(6,6,6,6))
}

cs_coding <- bind_rows(lapply(c("Low","Medium","High"), function(cb) {
  sub <- saddle_coding %>% filter(filter_val==cb)
  if (nrow(sub)==0) return(NULL)
  extract_cs_table(sub, cb)
})) %>%
  dplyr::mutate(
    sample   = factor(.data$sample,  levels=COND_LEVELS),
    category = factor(category, levels=c("Low","Medium","High"))
  )

cs_expr <- bind_rows(lapply(c("Low","Medium","High"), function(eb) {
  sub <- saddle_expr %>% filter(filter_val==eb)
  if (nrow(sub)==0) return(NULL)
  extract_cs_table(sub, eb)
})) %>%
  dplyr::mutate(
    sample   = factor(.data$sample,  levels=COND_LEVELS),
    category = factor(category, levels=c("Low","Medium","High"))
  )

message("\n══ CS by coding bin ═════════════════════════════════════════════════")
print(plot_cs_faceted(cs_coding, c("Low","Medium","High"),
                      title="Compartment strength by coding density"))
print(plot_corner_scores(cs_coding, c("Low","Medium","High"),
                         title="Per-corner log2 O/E by coding density"))

message("\n══ CS by expression bin ═════════════════════════════════════════════")
print(plot_cs_faceted(cs_expr, c("Low","Medium","High"),
                      title="Compartment strength by expression bin"))
print(plot_corner_scores(cs_expr, c("Low","Medium","High"),
                         title="Per-corner log2 O/E by expression bin"))

message("\nCoding/noncoding compartment analysis complete.")