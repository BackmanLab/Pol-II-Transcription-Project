##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Coding/noncoding compartment saddle analysis. DS_GROUP switch for
# downsampled/unsampled Hi-C files. Eigenvectors loaded from
# bedgraph files (SCALE-normalized, 500kb, from
# HiC_HPC_processing/juicer/scripts/run_compartments.sh).
#
# Figure output: Fig. 4C-E (probable); SI Fig. 7G (confirmed --
# literal title "log2(coding/noncoding) ratio per compartment" and
# "Rows = coding density bin | columns = condition"); SI Fig. 7H
# (probable, CS=(AA+BB)/(AB+BA) formula and log2 O/E labels).
##################################################################
library(strawr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(parallel)
library(ggside)
library(ggnewscale)
library(scales)

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

RDS_DIR <- file.path(BASE_NEW, "compartments/rds", DS_GROUP)
dir.create(RDS_DIR, showWarnings=FALSE, recursive=TRUE)

EXPR_FILE      <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
DESEQ_FILE     <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
RESOLUTION     <- 500000
SADDLE_BINS    <- 20L
STRAW_NORM     <- "SCALE"
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1

get_slurm_cores <- function(default=2) {
  cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (nchar(cpus)>0) return(as.integer(cpus))
  ntasks <- Sys.getenv("SLURM_NTASKS")
  if (nchar(ntasks)>0) return(as.integer(ntasks))
  default
}
N_CORES <- get_slurm_cores(default=2)
message("Using ", N_CORES, " cores")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Expression data
# ══════════════════════════════════════════════════════════════════════════════
message("Loading expression data...")
expr <- read_tsv(EXPR_FILE, show_col_types=FALSE) %>%
  dplyr::select(gene_id, TPM)
deseq_res <- read_csv(DESEQ_FILE, show_col_types=FALSE) %>%
  dplyr::mutate(
    gene_id = sub("\\..*","",gene_id),
    symbol  = dplyr::case_when(!is.na(symbol)&symbol!=""~symbol,TRUE~gene_id)
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
                     ranges=IRanges(eigen_df$start+1,eigen_df$end))
  hits      <- findOverlaps(bin_gr,gene_gr)
  high_bins <- unique(queryHits(hits))
  mean_high <- mean(eigen_df$eigen[high_bins],na.rm=TRUE)
  pct_pos   <- mean(eigen_df$eigen[high_bins]>0,na.rm=TRUE)*100
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

message("\nChecking eigenvector signs...")
eigen_list <- lapply(names(eigen_list), function(n) {
  if (is.null(eigen_list[[n]])) return(NULL)
  check_eigen_sign(eigen_list[[n]], deseq_res, CONDITIONS[[n]]$label)
})
names(eigen_list) <- names(CONDITIONS)

eigen_list <- Filter(Negate(is.null), eigen_list)
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
  bin_gr  <- tileGenome(seqlengths=seqlens, tilewidth=resolution,
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
                  exon_bp, noncoding_bp, coding_bin)
  result <- eigen_df %>%
    left_join(coding_sub, by=c("chr","start")) %>%
    dplyr::mutate(
      coding_ratio = replace_na(coding_ratio, 0),
      exon_bp      = replace_na(exon_bp,      0),
      noncoding_bp = replace_na(noncoding_bp, 0),
      coding_bin   = replace_na(as.character(coding_bin),"Low"),
      coding_bin   = factor(coding_bin, levels=c("Low","Medium","High"))
    )
  message("  ",sample_label," coding bins:"); print(table(result$coding_bin))
  result
}

message("\nMerging coding annotation...")
eigen_coding_list <- lapply(names(eigen_list), function(n) {
  merge_coding_to_eigen(eigen_list[[n]], coding_bins, CONDITIONS[[n]]$label)
})
names(eigen_coding_list) <- names(eigen_list)

message("\nCoding ratio vs eigenvector correlation:")
invisible(lapply(names(eigen_coding_list), function(n) {
  cor_df <- eigen_coding_list[[n]] %>% filter(!is.na(eigen), coding_ratio>0)
  message("  ",CONDITIONS[[n]]$label,": r = ",
          round(cor(cor_df$coding_ratio, cor_df$eigen, use="complete.obs"),3))
}))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Saddle computation
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
    ) %>% filter(!is.na(q1), !is.na(q2))
  if (nrow(raw)==0) return(NULL)
  expected_by_dist <- raw %>%
    dplyr::group_by(dist_bins) %>%
    dplyr::summarise(expected=mean(counts, na.rm=TRUE), .groups="drop")
  raw <- raw %>%
    left_join(expected_by_dist, by="dist_bins") %>%
    dplyr::mutate(oe=counts/expected) %>%
    filter(!is.na(oe), is.finite(oe))
  if (nrow(raw)==0) return(NULL)
  saddle_upper <- raw %>%
    dplyr::group_by(q1, q2) %>%
    dplyr::summarise(mean_oe=mean(oe, na.rm=TRUE), n=dplyr::n(), .groups="drop")
  bind_rows(saddle_upper,
            saddle_upper %>% dplyr::rename(q1=q2, q2=q1)) %>%
    distinct(q1, q2, .keep_all=TRUE) %>%
    dplyr::mutate(chr=as.character(chrom))
}

compute_stratified_saddle_coding <- function(hic_file, eigen_coding,
                                             sample_name, coding_bin_val,
                                             n_bins=SADDLE_BINS,
                                             resolution=RESOLUTION,
                                             norm=STRAW_NORM,
                                             n_cores=N_CORES) {
  eigen_subset <- eigen_coding %>%
    filter(as.character(coding_bin)==coding_bin_val)
  message("  ",coding_bin_val,": ",nrow(eigen_subset)," bins")
  if (nrow(eigen_subset)<n_bins*2) { message("  Too few — skipping"); return(NULL) }
  chroms <- unique(eigen_subset$chr)
  results_list <- parallel::mclapply(chroms, function(chrom) {
    eigen_chr <- eigen_subset %>% filter(chr==chrom)
    if (nrow(eigen_chr)<n_bins) return(NULL)
    compute_saddle_chr(hic_file, chrom, eigen_chr,
                       n_bins, resolution, norm)
  }, mc.cores=n_cores, mc.preschedule=FALSE)
  result <- bind_rows(Filter(Negate(is.null), results_list))
  if (nrow(result)==0) return(NULL)
  result %>%
    dplyr::group_by(q1, q2) %>%
    dplyr::summarise(mean_oe=mean(mean_oe, na.rm=TRUE),
                     n=sum(n), .groups="drop") %>%
    dplyr::mutate(
      cond_label = sample_name,
      log2_oe    = log2(mean_oe),
      filter_col = "coding_bin",
      filter_val = coding_bin_val
    )
}

coding_saddle_rds <- file.path(RDS_DIR, "saddle_coding_ratio_500kb_SCALE.rds")
if (file.exists(coding_saddle_rds)) {
  message("Loading saved coding saddle matrices...")
  saddle_coding <- readRDS(coding_saddle_rds)
} else {
  message("\n══ Computing coding-stratified saddle matrices ══════════════════")
  coding_results <- list()
  for (n in names(eigen_coding_list)) {
    hic_file    <- CONDITIONS[[n]]$hic
    samp        <- CONDITIONS[[n]]$label
    eigen_annot <- eigen_coding_list[[n]]
    message("\n--- ",samp," ---")
    for (cb in c("Low","Medium","High")) {
      res <- compute_stratified_saddle_coding(
        hic_file=hic_file, eigen_coding=eigen_annot,
        sample_name=samp, coding_bin_val=cb,
        norm=STRAW_NORM)
      if (!is.null(res)) coding_results[[paste(n,cb,sep="_")]] <- res
    }
  }
  saddle_coding <- bind_rows(coding_results) %>%
    dplyr::rename(sample=cond_label) %>%
    dplyr::mutate(sample=factor(.data$sample, levels=COND_LEVELS))
  saveRDS(saddle_coding, coding_saddle_rds)
  message("Saved: ",coding_saddle_rds)
}
message("\nSaddle coding summary:")
print(table(as.character(saddle_coding$sample), saddle_coding$filter_val))

# ══════════════════════════════════════════════════════════════════════════════
# 6. Core helpers
# ══════════════════════════════════════════════════════════════════════════════
compartment_strength <- function(saddle_df, n_bins=SADDLE_BINS, corner_size=5) {
  lo <- 1:corner_size; hi <- (n_bins-corner_size+1):n_bins
  get_mean <- function(q1r, q2r)
    saddle_df %>% filter(q1%in%q1r, q2%in%q2r) %>%
    summarise(m=mean(mean_oe, na.rm=TRUE)) %>% pull(m)
  aa <- get_mean(hi,hi); bb <- get_mean(lo,lo)
  ab <- get_mean(lo,hi); ba <- get_mean(hi,lo)
  list(AA=aa, BB=bb, AB=ab, BA=ba, CS=(aa+bb)/(ab+ba))
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. plot_saddle_coding_single
# ══════════════════════════════════════════════════════════════════════════════
plot_saddle_coding_single <- function(saddle_df, eigen_coding_df,
                                      sample_name, coding_bin_val,
                                      n_bins=SADDLE_BINS,
                                      title=NULL, color_limit=NULL) {
  mat_df <- saddle_df %>%
    dplyr::filter(as.character(.data$sample)==sample_name,
                  as.character(filter_val)==coding_bin_val) %>%
    tidyr::complete(q1=1:n_bins, q2=1:n_bins, fill=list(log2_oe=NA))
  if (nrow(mat_df)==0||all(is.na(mat_df$log2_oe))) return(NULL)
  if (is.null(color_limit))
    color_limit <- quantile(abs(mat_df$log2_oe),0.95,na.rm=TRUE)
  saddle_sub <- saddle_df %>%
    dplyr::filter(as.character(.data$sample)==sample_name,
                  as.character(filter_val)==coding_bin_val)
  cs_vals  <- compartment_strength(saddle_sub)
  cs_label <- paste0("CS=",round(cs_vals$CS,2))
  corner_scores <- data.frame(
    q1   =c(1.5,       1.5,        n_bins-0.5, n_bins-0.5),
    q2   =c(n_bins-0.5,1.5,        n_bins-0.5, 1.5),
    label=c(paste0("BA\n",round(log2(cs_vals$BA),2)),
            paste0("BB\n",round(log2(cs_vals$BB),2)),
            paste0("AA\n",round(log2(cs_vals$AA),2)),
            paste0("AB\n",round(log2(cs_vals$AB),2))),
    hjust=c(0,0,1,1), vjust=c(1,0,1,0)
  )
  eigen_sub <- eigen_coding_df %>%
    dplyr::filter(as.character(coding_bin)==coding_bin_val, !is.na(eigen)) %>%
    dplyr::mutate(quantile_bin=ntile(eigen, n_bins))
  if (nrow(eigen_sub)==0) return(NULL)
  coding_marginals <- eigen_sub %>%
    dplyr::group_by(quantile_bin) %>%
    dplyr::summarise(
      mean_exon_bp     =mean(exon_bp,     na.rm=TRUE),
      mean_noncoding_bp=mean(noncoding_bp,na.rm=TRUE),
      .groups="drop") %>%
    dplyr::mutate(
      quantile_bin      =as.numeric(quantile_bin),
      log2_coding_ratio =log2((mean_exon_bp+1)/(mean_noncoding_bp+1)),
      direction=factor(
        ifelse(log2_coding_ratio>=0,"Coding-enriched","Noncoding-enriched"),
        levels=c("Noncoding-enriched","Coding-enriched"))
    )
  ratio_lim    <- max(abs(coding_marginals$log2_coding_ratio),na.rm=TRUE)*1.1
  ratio_breaks <- round(c(-ratio_lim*0.8,0,ratio_lim*0.8),2)
  subtitle_label <- paste0(nrow(eigen_sub)," bins  |  mean coding ratio = ",
                           round(mean(eigen_sub$coding_ratio,na.rm=TRUE),4))
  ggplot(mat_df,aes(x=q1,y=q2))+
    geom_tile(aes(fill=pmin(pmax(log2_oe,-color_limit),color_limit)))+
    scale_fill_distiller(palette="RdBu",direction=-1,
                         limits=c(-color_limit,color_limit),
                         name="log2 O/E",na.value="gray80")+
    geom_text(data=corner_scores,
              aes(x=q1,y=q2,label=label,hjust=hjust,vjust=vjust),
              color="black",size=2.8,fontface="bold",
              lineheight=0.85,inherit.aes=FALSE)+
    annotate("text",x=n_bins/2,y=n_bins/2,label=cs_label,
             size=3.2,fontface="bold",color="gray20")+
    ggnewscale::new_scale_fill()+
    ggside::geom_xsidecol(data=coding_marginals,
                          aes(x=quantile_bin,y=log2_coding_ratio,fill=direction),
                          inherit.aes=FALSE,width=0.9)+
    ggside::geom_xsidehline(yintercept=0,color="gray30",linewidth=0.4)+
    ggside::scale_xfill_manual(
      values=c("Coding-enriched"="#e6550d","Noncoding-enriched"="#bdbdbd"),
      guide="none")+
    scale_xsidey_continuous(name="log2(E/I)",limits=c(-ratio_lim,ratio_lim),
                            breaks=ratio_breaks,
                            labels=scales::number(ratio_breaks,accuracy=0.01))+
    ggside::geom_ysidecol(data=coding_marginals,
                          aes(y=quantile_bin,x=log2_coding_ratio,fill=direction),
                          inherit.aes=FALSE,width=0.9)+
    ggside::geom_ysidehline(yintercept=0,color="gray30",linewidth=0.4)+
    ggside::scale_yfill_manual(
      values=c("Coding-enriched"="#e6550d","Noncoding-enriched"="#bdbdbd"),
      guide="none")+
    scale_ysidex_continuous(name="log2(E/I)",limits=c(-ratio_lim,ratio_lim),
                            breaks=ratio_breaks,
                            labels=scales::number(ratio_breaks,accuracy=0.01))+
    scale_x_continuous(name="Eigenvector quantile",
                       breaks=c(1,n_bins/2,n_bins),
                       labels=c("B","0","A"),expand=c(0,0))+
    scale_y_continuous(name="Eigenvector quantile",
                       breaks=c(1,n_bins/2,n_bins),
                       labels=c("B","0","A"),expand=c(0,0))+
    coord_fixed()+
    labs(title=title,subtitle=subtitle_label)+
    theme_classic(base_size=11)+
    ggside::theme_ggside_void()+
    theme(
      plot.title          =element_text(face="bold",size=11,hjust=0.5),
      plot.subtitle       =element_text(size=7,color="gray40",hjust=0.5),
      legend.position     ="none",
      axis.text           =element_text(size=9,color="black"),
      axis.title          =element_text(size=10,color="black",face="bold"),
      axis.ticks          =element_line(color="black",linewidth=0.4),
      plot.margin         =margin(4,4,4,4),
      ggside.panel.scale.x=0.25,
      ggside.panel.scale.y=0.25,
      ggside.axis.text.x  =element_text(size=7,color="black"),
      ggside.axis.text.y  =element_text(size=7,color="black"),
      ggside.axis.ticks.x =element_line(color="black",linewidth=0.3),
      ggside.axis.ticks.y =element_line(color="black",linewidth=0.3)
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Assemble per-condition figure
# ══════════════════════════════════════════════════════════════════════════════
assemble_saddle_coding_per_condition <- function(saddle_coding, eigen_coding_list,
                                                 condition_name, sample_label) {
  saddle_sub <- saddle_coding %>%
    dplyr::filter(as.character(.data$sample)==sample_label)
  message("  Rows for ",sample_label,": ",nrow(saddle_sub))
  if (nrow(saddle_sub)==0) return(NULL)
  color_limit <- quantile(abs(saddle_sub$log2_oe),0.95,na.rm=TRUE)
  eigen_cd    <- eigen_coding_list[[condition_name]]
  panels <- lapply(c("Low","Medium","High"), function(cb) {
    message("    Building panel: ",cb)
    plot_saddle_coding_single(
      saddle_df=saddle_sub, eigen_coding_df=eigen_cd,
      sample_name=sample_label, coding_bin_val=cb,
      n_bins=SADDLE_BINS, color_limit=color_limit,
      title=paste0(cb," coding density"))
  })
  panels <- Filter(Negate(is.null), panels)
  if (length(panels)==0) return(NULL)
  heatmap_legend <- cowplot::get_legend(
    ggplot(data.frame(x=1:100,y=1,
                      fill=seq(-color_limit,color_limit,length.out=100)),
           aes(x=x,y=y,fill=fill))+geom_tile()+
      scale_fill_distiller(palette="RdBu",direction=-1,
                           limits=c(-color_limit,color_limit),name="log2 O/E")+
      theme_void()+
      theme(legend.position="right",legend.key.height=unit(1.5,"cm"),
            legend.key.width=unit(0.4,"cm"),legend.text=element_text(size=9),
            legend.title=element_text(size=10)))
  marginal_legend <- cowplot::get_legend(
    ggplot(data.frame(
      direction=factor(c("Coding-enriched","Noncoding-enriched"),
                       levels=c("Noncoding-enriched","Coding-enriched")),
      x=c(1,1),y=c(1,2)),aes(x=x,y=y,fill=direction))+
      geom_col()+
      scale_fill_manual(values=c("Coding-enriched"="#e6550d",
                                 "Noncoding-enriched"="#bdbdbd"),
                        name="log2(E/I)\nmarginal")+
      theme_void()+
      theme(legend.position="right",legend.key.size=unit(0.5,"cm"),
            legend.text=element_text(size=9),
            legend.title=element_text(size=10,face="bold")))
  legends   <- cowplot::plot_grid(heatmap_legend,marginal_legend,
                                  ncol=1,rel_heights=c(1,0.5))
  panel_row <- cowplot::plot_grid(plotlist=panels,nrow=1,
                                  labels=c("a","b","c"),label_size=12)
  composite <- cowplot::plot_grid(panel_row,legends,ncol=2,rel_widths=c(1,0.08))
  cond_color <- COND_COLORS[sample_label]
  title_grob <- cowplot::ggdraw()+
    cowplot::draw_label(
      paste0(sample_label," — compartment contacts by coding/noncoding ratio"),
      fontface="bold",size=13,x=0.01,hjust=0,color=cond_color)+
    cowplot::draw_label(
      paste0("500kb  |  ",STRAW_NORM,"  |  ",SADDLE_BINS,
             " quantile bins  |  marginals = log2(exon/noncoding)  |  ",DS_GROUP),
      size=8,color="gray40",x=0.01,hjust=0,y=0.2)
  final <- cowplot::plot_grid(title_grob,composite,ncol=1,rel_heights=c(0.05,1))
  print(final); invisible(final)
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Run saddle plots
# ══════════════════════════════════════════════════════════════════════════════
for (n in names(eigen_coding_list)) {
  samp <- CONDITIONS[[n]]$label
  message("\n══ ",samp," — coding/noncoding saddle ═══════════════════════════")
  assemble_saddle_coding_per_condition(
    saddle_coding=saddle_coding, eigen_coding_list=eigen_coding_list,
    condition_name=n, sample_label=samp)
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. Supplementary plots
# ══════════════════════════════════════════════════════════════════════════════
eigen_coding_all <- bind_rows(lapply(names(eigen_coding_list), function(n)
  eigen_coding_list[[n]] %>% dplyr::mutate(condition=CONDITIONS[[n]]$label)
)) %>%
  dplyr::rename(sample=condition) %>%
  dplyr::mutate(
    sample        = factor(sample,     levels=COND_LEVELS),
    coding_bin    = factor(coding_bin, levels=c("Low","Medium","High")),
    log2_cn_ratio = log2((exon_bp+1)/(noncoding_bp+1)),
    compartment   = dplyr::case_when(
      is.na(eigen) ~ "Unassigned",
      eigen >  0   ~ "A",
      eigen <  0   ~ "B",
      TRUE         ~ "Boundary"),
    compartment   = factor(compartment, levels=c("A","B","Boundary","Unassigned"))
  )
message("Combined eigen_coding_all: ",nrow(eigen_coding_all)," rows")

compartment_colors <- c(A="#d6604d", B="#4393c3")

plot_coding_ratio_histogram <- function(eigen_coding_all,
                                        title="log2(coding/noncoding) ratio per compartment") {
  plot_df <- eigen_coding_all %>%
    filter(!is.na(log2_cn_ratio),!is.na(eigen),compartment%in%c("A","B"))
  medians <- plot_df %>%
    dplyr::group_by(sample,compartment) %>%
    dplyr::summarise(median_ratio=median(log2_cn_ratio,na.rm=TRUE),
                     n_bins=dplyr::n(),.groups="drop")
  p <- ggplot(plot_df,aes(x=log2_cn_ratio,fill=compartment,color=compartment))+
    geom_histogram(bins=40,alpha=0.6,position="identity",linewidth=0.2)+
    geom_vline(data=medians,aes(xintercept=median_ratio,color=compartment),
               linetype="dashed",linewidth=0.7)+
    geom_vline(xintercept=0,color="gray30",linewidth=0.4)+
    scale_fill_manual(values=compartment_colors,name="Compartment")+
    scale_color_manual(values=compartment_colors,name="Compartment")+
    scale_x_continuous(name="log2(exon bp / noncoding bp)",
                       breaks=scales::pretty_breaks(n=6))+
    scale_y_continuous(name="Number of 500kb bins",
                       expand=expansion(mult=c(0,0.1)))+
    facet_wrap(~sample,nrow=1)+
    labs(title=title,
         subtitle="dashed = median per compartment  |  vertical line = equal coding/noncoding")+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=13,hjust=0.5),
          plot.subtitle=element_text(size=8,color="gray40",hjust=0.5),
          strip.text=element_text(face="bold",size=11),
          strip.background=element_blank(),
          legend.position="bottom")
  print(p); invisible(p)
}

plot_coding_ratio_histogram_by_bin <- function(eigen_coding_all,
                                               title="log2(coding/noncoding) by coding density bin") {
  plot_df <- eigen_coding_all %>%
    filter(!is.na(log2_cn_ratio),!is.na(eigen),
           compartment%in%c("A","B"),!is.na(coding_bin))
  p <- ggplot(plot_df,aes(x=log2_cn_ratio,fill=compartment,color=compartment))+
    geom_histogram(bins=30,alpha=0.6,position="identity",linewidth=0.2)+
    geom_vline(xintercept=0,color="gray30",linewidth=0.4)+
    scale_fill_manual(values=compartment_colors,name="Compartment")+
    scale_color_manual(values=compartment_colors,name="Compartment")+
    scale_x_continuous(name="log2(exon bp / noncoding bp)",
                       breaks=scales::pretty_breaks(n=4))+
    scale_y_continuous(name="Number of 500kb bins",
                       expand=expansion(mult=c(0,0.1)))+
    facet_grid(coding_bin~sample)+
    labs(title=title,
         subtitle="Rows = coding density bin  |  columns = condition")+
    theme_classic(base_size=11)+
    theme(plot.title=element_text(face="bold",size=13,hjust=0.5),
          strip.text=element_text(face="bold",size=10),
          strip.background=element_blank(),
          legend.position="bottom")
  print(p); invisible(p)
}

plot_eigen_vs_coding_ratio <- function(eigen_coding_all,
                                       title="Eigenvector vs log2(coding/noncoding) ratio") {
  plot_df <- eigen_coding_all %>%
    filter(!is.na(log2_cn_ratio),!is.na(eigen),compartment%in%c("A","B"))
  cors <- plot_df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(r=round(cor(eigen,log2_cn_ratio,use="complete.obs"),3),
                     n=dplyr::n(),.groups="drop") %>%
    dplyr::mutate(label=paste0("r = ",r,"\nn = ",n))
  p <- ggplot(plot_df,aes(x=log2_cn_ratio,y=eigen,color=compartment))+
    geom_hline(yintercept=0,linetype="dashed",color="gray50",linewidth=0.4)+
    geom_vline(xintercept=0,linetype="dashed",color="gray50",linewidth=0.4)+
    geom_point(alpha=0.5,size=1.2)+
    geom_smooth(aes(group=1),method="loess",formula=y~x,
                color="gray20",fill="gray80",linewidth=0.8,alpha=0.3,se=TRUE)+
    geom_text(data=cors,aes(x=-Inf,y=Inf,label=label),
              hjust=-0.1,vjust=1.3,size=3.2,color="black",
              fontface="italic",inherit.aes=FALSE)+
    scale_color_manual(values=compartment_colors,name="Compartment")+
    scale_x_continuous(name="log2(exon bp / noncoding bp)",
                       breaks=scales::pretty_breaks(n=5))+
    scale_y_continuous(name="PC1 (eigenvector)",
                       breaks=scales::pretty_breaks(n=5))+
    facet_wrap(~sample,nrow=1)+
    labs(title=title,
         subtitle="Each point = one 500kb bin  |  loess smooth across all bins")+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=13,hjust=0.5),
          plot.subtitle=element_text(size=8,color="gray40",hjust=0.5),
          strip.text=element_text(face="bold",size=11),
          strip.background=element_blank(),
          legend.position="bottom")
  print(p); invisible(p)
}

plot_eigen_vs_coding_ratio_by_bin <- function(eigen_coding_all,
                                              title="Eigenvector vs log2(C/N) by coding density bin") {
  plot_df <- eigen_coding_all %>%
    filter(!is.na(log2_cn_ratio),!is.na(eigen),!is.na(coding_bin))
  bin_colors <- c(Low="#4393c3",Medium="#f7f7f7",High="#d6604d")
  cors <- plot_df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(r=round(cor(eigen,log2_cn_ratio,use="complete.obs"),3),
                     .groups="drop") %>%
    dplyr::mutate(label=paste0("r = ",r))
  p <- ggplot(plot_df,aes(x=log2_cn_ratio,y=eigen,color=coding_bin))+
    geom_hline(yintercept=0,linetype="dashed",color="gray50",linewidth=0.4)+
    geom_vline(xintercept=0,linetype="dashed",color="gray50",linewidth=0.4)+
    geom_point(alpha=0.5,size=1.2)+
    geom_smooth(aes(group=coding_bin,color=coding_bin),
                method="lm",formula=y~x,linewidth=0.8,se=FALSE)+
    geom_text(data=cors,aes(x=-Inf,y=Inf,label=label),
              hjust=-0.1,vjust=1.3,size=3.2,color="black",
              fontface="italic",inherit.aes=FALSE)+
    scale_color_manual(values=bin_colors,name="Coding density")+
    scale_x_continuous(name="log2(exon bp / noncoding bp)",
                       breaks=scales::pretty_breaks(n=5))+
    scale_y_continuous(name="PC1 (eigenvector)",
                       breaks=scales::pretty_breaks(n=5))+
    facet_wrap(~sample,nrow=1)+
    labs(title=title,
         subtitle="Color = coding density tertile  |  lines = per-tertile linear fit")+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=13,hjust=0.5),
          strip.text=element_text(face="bold",size=11),
          strip.background=element_blank(),
          legend.position="bottom")
  print(p); invisible(p)
}

message("\n══ Supplementary plots ══════════════════════════════════════════════")
plot_coding_ratio_histogram(eigen_coding_all,
                            title="log2(coding/noncoding) ratio distribution by compartment")
plot_coding_ratio_histogram_by_bin(eigen_coding_all,
                                   title="log2(coding/noncoding) distribution by coding density and condition")
plot_eigen_vs_coding_ratio(eigen_coding_all,
                           title="PC1 eigenvector vs log2(coding/noncoding) ratio")
plot_eigen_vs_coding_ratio_by_bin(eigen_coding_all,
                                  title="PC1 eigenvector vs log2(coding/noncoding) by coding density bin")

message("\nCoding/noncoding compartment analysis complete.")