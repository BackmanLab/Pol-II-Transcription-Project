##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Aggregate TAD Analysis (ATA) pipeline. Uses GENOVA::ATA for
# canonical ATA computation. Chromosome-selective strawr loader
# replaces load_contacts.
#
# Figure output: SI Fig. 6A-B (probable).
##################################################################
library(GENOVA)
library(strawr)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(cowplot)

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════
DS_GROUP <- "ds_wt_pol2"   # ← "ds_wt_pol2" | "ds_all" | "unsampled"

BASE_NEW <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
BASE_OLD <- "/projects/p32171/juicer/work/112123_HiC/juicer_analysis"

if (DS_GROUP == "ds_wt_pol2") {
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE", scale=NULL,
                hic =file.path(BASE_NEW,"mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                tads=file.path(BASE_NEW,"mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")),
    pol2 = list(label="POLR2A deg", color="#CC6677", scale=NULL,
                hic =file.path(BASE_NEW,"mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                tads=file.path(BASE_NEW,"mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe"))
  )
} else if (DS_GROUP == "ds_all") {
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE", scale=NULL,
                hic =file.path(BASE_NEW,"mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                tads=file.path(BASE_NEW,"mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")),
    pol2 = list(label="POLR2A deg", color="#CC6677", scale=NULL,
                hic =file.path(BASE_NEW,"mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                tads=file.path(BASE_NEW,"mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")),
    actd = list(label="ActD 1hr",   color="#44AA99", scale=NULL,
                hic =file.path(BASE_NEW,"mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic"),
                tads=file.path(BASE_NEW,"mega_ds_all/1hr_ActD/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe"))
  )
} else if (DS_GROUP == "unsampled") {
  depth_wt <- 639409109 + 377220354
  CONDITIONS <- list(
    wt   = list(label="Wildtype",   color="#88CCEE", scale=depth_wt,
                hic =file.path(BASE_OLD,"WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
                tads=file.path(BASE_OLD,"WT_HCT116_CTRL/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")),
    pol2 = list(label="POLR2A deg", color="#CC6677", scale=depth_wt,
                hic =file.path(BASE_OLD,"Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
                tads=file.path(BASE_OLD,"Pol2_6hrs_Aux/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe")),
    actd = list(label="ActD 1hr",   color="#44AA99", scale=depth_wt,
                hic =file.path(BASE_OLD,"1hr_ActD/mega/aligned/inter_30.hic"),
                tads=file.path(BASE_OLD,"1hr_ActD/mega/aligned/inter_30_contact_domains/5000_blocks.bedpe"))
  )
} else { stop("Unknown DS_GROUP: '", DS_GROUP, "'") }
message("DS_GROUP: ", DS_GROUP)

EXPR_FILE      <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
DESEQ_FILE     <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
ENCODE_DIR     <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/encode"
TAD_MIN_KB     <- 50
TAD_MAX_KB     <- 5000
GENE_TOL_BP    <- 5000
CTCF_TOL_BP    <- 10000
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1
ATA_SIZE       <- 100L
RESOLUTION     <- 25000
STRAW_NORM     <- "SCALE"
BINS           <- c("Low","Medium","High")
BIN_COLORS     <- c("Low"="#fdae61","Medium"="#f46d43","High"="#d73027")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Expression data
# ══════════════════════════════════════════════════════════════════════════════
expr <- read_tsv(EXPR_FILE) %>% dplyr::select(gene_id, TPM)
deseq_res <- read_csv(DESEQ_FILE) %>%
  dplyr::mutate(
    gene_id   = sub("\\..*","",gene_id),
    de_status = dplyr::case_when(
      is.na(padj)                                              ~ "Filtered",
      padj >= PADJ_THRESHOLD                                   ~ "Not DE",
      padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD ~ "Down",
      TRUE                                                     ~ "Not DE"),
    de_status = factor(de_status,levels=c("Up","Not DE","Down","Filtered")),
    symbol    = dplyr::case_when(!is.na(symbol)&symbol!=""~symbol,TRUE~gene_id)
  ) %>% left_join(expr,by="gene_id")
expressed  <- deseq_res %>% filter(!is.na(TPM),TPM>0) %>% pull(TPM)
tpm_medium <- quantile(expressed,1/3,na.rm=TRUE)
tpm_high   <- quantile(expressed,2/3,na.rm=TRUE)
deseq_res  <- deseq_res %>%
  dplyr::mutate(
    expr_bin = dplyr::case_when(
      is.na(TPM)|TPM==0~"Low",TPM<=tpm_medium~"Low",
      TPM<=tpm_high~"Medium",TRUE~"High"),
    expr_bin = factor(expr_bin,levels=BINS))
message("DE breakdown:"); print(table(deseq_res$de_status))

gene_coords <- deseq_res %>%
  filter(!is.na(gene_seq_start),!is.na(gene_seq_end),!is.na(seq_name)) %>%
  distinct(gene_id,.keep_all=TRUE)
gene_gr <- GRanges(
  seqnames=as.character(gene_coords$seq_name),
  ranges  =IRanges(gene_coords$gene_seq_start,gene_coords$gene_seq_end),
  symbol  =gene_coords$symbol, gene_id=gene_coords$gene_id,
  de_status=as.character(gene_coords$de_status),
  expr_bin =as.character(gene_coords$expr_bin), TPM=gene_coords$TPM)
message("gene_gr: ",length(gene_gr)," genes")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Load TADs
# ══════════════════════════════════════════════════════════════════════════════
load_arrowhead <- function(path,min_size_kb=TAD_MIN_KB,max_size_kb=TAD_MAX_KB) {
  lines <- readLines(path)
  df    <- read.table(text=paste(lines[!grepl("^#",lines)],collapse="\n"),
                      header=FALSE,sep="\t",stringsAsFactors=FALSE,fill=TRUE)
  colnames(df)[1:6] <- c("chr1","x1","x2","chr2","y1","y2")
  df %>% dplyr::mutate(tad_size_kb=(x2-x1)/1000) %>%
    filter(!is.na(x1),!is.na(x2),x2>x1,
           tad_size_kb>=min_size_kb,tad_size_kb<=max_size_kb)
}
tads_raw <- lapply(names(CONDITIONS),function(n) load_arrowhead(CONDITIONS[[n]]$tads))
names(tads_raw) <- names(CONDITIONS)
invisible(lapply(names(tads_raw),function(n)
  message(CONDITIONS[[n]]$label," TADs: ",nrow(tads_raw[[n]]))))

# ══════════════════════════════════════════════════════════════════════════════
# 3. Link TADs to genes
# ══════════════════════════════════════════════════════════════════════════════
link_tads_anchor <- function(tads_df,gene_gr,tol_bp=GENE_TOL_BP) {
  tads_df <- tads_df %>% dplyr::mutate(tad_idx=row_number())
  fix_chr <- function(gr,ref_gr) {
    has_chr_ref<-any(grepl("^chr",as.character(seqnames(ref_gr))))
    has_chr_gr <-any(grepl("^chr",as.character(seqnames(gr))))
    if (has_chr_ref&&!has_chr_gr) seqlevels(gr)<-paste0("chr",seqlevels(gr))
    else if (!has_chr_ref&&has_chr_gr) seqlevels(gr)<-gsub("^chr","",seqlevels(gr))
    gr
  }
  gene_gr_exp <- trim(resize(gene_gr,width(gene_gr)+2*tol_bp,fix="center"))
  left_gr  <- fix_chr(GRanges(tads_df$chr1,IRanges(pmax(0,tads_df$x1-tol_bp),tads_df$x1+tol_bp)),gene_gr)
  right_gr <- fix_chr(GRanges(tads_df$chr1,IRanges(pmax(0,tads_df$x2-tol_bp),tads_df$x2+tol_bp)),gene_gr)
  all_pairs <- bind_rows(
    data.frame(tad_idx=queryHits(findOverlaps(left_gr, gene_gr_exp)),
               gene_idx=subjectHits(findOverlaps(left_gr, gene_gr_exp))),
    data.frame(tad_idx=queryHits(findOverlaps(right_gr,gene_gr_exp)),
               gene_idx=subjectHits(findOverlaps(right_gr,gene_gr_exp)))
  ) %>% distinct(tad_idx,gene_idx) %>%
    dplyr::mutate(
      expr_bin =gene_gr$expr_bin[gene_idx],
      de_status=gene_gr$de_status[gene_idx],
      TPM      =gene_gr$TPM[gene_idx],
      symbol   =gene_gr$symbol[gene_idx])
  tad_annotations <- all_pairs %>%
    dplyr::group_by(tad_idx) %>%
    dplyr::slice_max(TPM,n=1,with_ties=FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(tad_idx,expr_bin,de_status,TPM,symbol)
  tads_df %>%
    left_join(tad_annotations,by="tad_idx") %>%
    dplyr::mutate(
      expr_bin =replace_na(expr_bin, "Low"),
      de_status=replace_na(de_status,"Not DE"),
      symbol   =replace_na(symbol,  "Unknown"))
}
message("Linking TADs to genes...")
tads_anchor <- lapply(tads_raw,link_tads_anchor,gene_gr=gene_gr)
invisible(lapply(names(tads_anchor),function(n) {
  message(CONDITIONS[[n]]$label," bin breakdown:")
  print(table(tads_anchor[[n]]$expr_bin))
}))

# ══════════════════════════════════════════════════════════════════════════════
# 4. CTCF filtering
# ══════════════════════════════════════════════════════════════════════════════
ccre_all <- read.table(file.path(ENCODE_DIR,"wt_hct116.bed"),
                       header=FALSE,sep="\t",stringsAsFactors=FALSE) %>%
  dplyr::select(chr=V1,start=V2,end=V3,id=V4,ccre_type=V10)
ctcf_gr <- GRanges(seqnames=ccre_all$chr[ccre_all$ccre_type=="CA-CTCF"],
                   ranges=IRanges(ccre_all$start[ccre_all$ccre_type=="CA-CTCF"],
                                  ccre_all$end[ccre_all$ccre_type=="CA-CTCF"]))
promoter_gr <- GRanges(seqnames=ccre_all$chr[ccre_all$ccre_type=="PLS"],
                       ranges=IRanges(ccre_all$start[ccre_all$ccre_type=="PLS"],
                                      ccre_all$end[ccre_all$ccre_type=="PLS"]))
enhancer_gr <- GRanges(seqnames=ccre_all$chr[ccre_all$ccre_type %in% c("pELS","dELS")],
                       ranges=IRanges(ccre_all$start[ccre_all$ccre_type %in% c("pELS","dELS")],
                                      ccre_all$end[ccre_all$ccre_type %in% c("pELS","dELS")]))

filter_ctcf_tads <- function(tads_df,ctcf_gr,tol_bp=CTCF_TOL_BP) {
  fix_chr <- function(gr,ref_gr) {
    has_chr_ref<-any(grepl("^chr",as.character(seqnames(ref_gr))))
    has_chr_gr <-any(grepl("^chr",as.character(seqnames(gr))))
    if (has_chr_ref&&!has_chr_gr) seqlevels(gr)<-paste0("chr",seqlevels(gr))
    else if (!has_chr_ref&&has_chr_gr) seqlevels(gr)<-gsub("^chr","",seqlevels(gr))
    gr
  }
  left_gr  <- fix_chr(GRanges(tads_df$chr1,IRanges(pmax(0,tads_df$x1-tol_bp),tads_df$x1+tol_bp)),ctcf_gr)
  right_gr <- fix_chr(GRanges(tads_df$chr1,IRanges(pmax(0,tads_df$x2-tol_bp),tads_df$x2+tol_bp)),ctcf_gr)
  both_ctcf <- intersect(unique(queryHits(findOverlaps(left_gr, ctcf_gr))),
                         unique(queryHits(findOverlaps(right_gr,ctcf_gr))))
  list(ctcf    =tads_df[both_ctcf,],
       non_ctcf=tads_df[setdiff(seq_len(nrow(tads_df)),both_ctcf),])
}
filter_pure_ctcf <- function(ctcf_tads,promoter_gr,enhancer_gr,tol_bp=CTCF_TOL_BP) {
  fix_chr <- function(gr,ref_gr) {
    has_chr_ref<-any(grepl("^chr",as.character(seqnames(ref_gr))))
    has_chr_gr <-any(grepl("^chr",as.character(seqnames(gr))))
    if (has_chr_ref&&!has_chr_gr) seqlevels(gr)<-paste0("chr",seqlevels(gr))
    else if (!has_chr_ref&&has_chr_gr) seqlevels(gr)<-gsub("^chr","",seqlevels(gr))
    gr
  }
  left_gr  <- fix_chr(GRanges(ctcf_tads$chr1,IRanges(pmax(0,ctcf_tads$x1-tol_bp),ctcf_tads$x1+tol_bp)),promoter_gr)
  right_gr <- fix_chr(GRanges(ctcf_tads$chr1,IRanges(pmax(0,ctcf_tads$x2-tol_bp),ctcf_tads$x2+tol_bp)),promoter_gr)
  excl <- unique(c(queryHits(findOverlaps(left_gr, promoter_gr)),
                   queryHits(findOverlaps(right_gr,promoter_gr)),
                   queryHits(findOverlaps(left_gr, enhancer_gr)),
                   queryHits(findOverlaps(right_gr,enhancer_gr))))
  keep <- setdiff(seq_len(nrow(ctcf_tads)),excl)
  message("  CTCF input:",nrow(ctcf_tads)," excluded:",length(excl)," pure CTCF:",length(keep))
  ctcf_tads[keep,]
}

message("\nSplitting TADs by CTCF status...")
message("WT:");   wt_split       <- filter_ctcf_tads(tads_anchor$wt,  ctcf_gr)
message("  CTCF:",nrow(wt_split$ctcf),"  Non-CTCF:",nrow(wt_split$non_ctcf))
wt_pure_ctcf   <- filter_pure_ctcf(wt_split$ctcf,promoter_gr,enhancer_gr)
message("pol2:"); pol2_split     <- filter_ctcf_tads(tads_anchor$pol2,ctcf_gr)
message("  CTCF:",nrow(pol2_split$ctcf),"  Non-CTCF:",nrow(pol2_split$non_ctcf))
pol2_pure_ctcf <- filter_pure_ctcf(pol2_split$ctcf,promoter_gr,enhancer_gr)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Chromosome-selective GENOVA experiment loader
#
#    All TAD subsets are collected first to determine which chromosomes are
#    needed, then only those chromosomes are queried from each .hic file.
#    A single shared IDX is built from the WT .hic chromosome sizes so all
#    experiments have identical bin coordinate systems — no sync_indices needed.
# ══════════════════════════════════════════════════════════════════════════════

# Collect all chromosomes needed across every TAD subset
all_tad_chroms <- unique(gsub("^chr","", c(
  tads_anchor$wt$chr1,
  tads_anchor$pol2$chr1)))
message("\nChromosomes needed: ",paste(sort(all_tad_chroms),collapse=", "))

# Build shared IDX from WT .hic — used for all conditions
message("Building shared IDX from WT .hic...")
wt_hic_chroms <- strawr::readHicChroms(CONDITIONS$wt$hic) %>%
  filter(name != "ALL", !grepl("_",name)) %>%
  filter(grepl("^[0-9]+$",name) | name == "X") %>%
  filter(name %in% all_tad_chroms) %>%
  arrange(suppressWarnings(as.integer(ifelse(name=="X",23,name))))

idx_list   <- vector("list",nrow(wt_hic_chroms))
bin_offset <- 0L
for (i in seq_len(nrow(wt_hic_chroms))) {
  chrom     <- wt_hic_chroms$name[i]
  chrom_len <- wt_hic_chroms$length[i]
  n_bins    <- ceiling(chrom_len / RESOLUTION)
  starts    <- seq(0L, by=RESOLUTION, length.out=n_bins)
  idx_list[[i]] <- data.table::data.table(
    V1 = chrom,
    V2 = as.integer(starts),
    V3 = as.integer(pmin(starts + RESOLUTION, chrom_len)),
    V4 = seq_len(n_bins) + bin_offset
  )
  bin_offset <- bin_offset + n_bins
}
SHARED_IDX <- data.table::rbindlist(idx_list)
data.table::setkey(SHARED_IDX, V1, V2)
message("Shared IDX: ",nrow(SHARED_IDX)," bins  max_bin=",max(SHARED_IDX$V4))
message("IDX chr names: ",paste(unique(SHARED_IDX$V1),collapse=", "))

# Verify TAD chr names match IDX
tad_chroms_check <- unique(gsub("^chr","",tads_anchor$wt$chr1))
missing <- setdiff(tad_chroms_check, unique(SHARED_IDX$V1))
if (length(missing) > 0)
  warning("TAD chromosomes not in IDX: ", paste(missing, collapse=", "))

# Build GENOVA contacts object from strawr — per condition

load_contacts_selective <- function(hic_file, shared_idx,
                                    resolution  = RESOLUTION,
                                    norm        = STRAW_NORM,
                                    sample_name = "sample",
                                    colour      = "#888888",
                                    scale_bp    = NULL) {
  tad_chroms <- unique(shared_idx$V1)
  message("  Chromosomes: ", paste(sort(tad_chroms), collapse=", "))
  mat_list       <- vector("list", length(tad_chroms))
  total_contacts <- 0L
  for (i in seq_along(tad_chroms)) {
    chrom <- tad_chroms[i]
    raw   <- tryCatch(
      strawr::straw(norm, hic_file, chrom, chrom, "BP", resolution),
      error=function(e) { message("    Warning chr",chrom,": ",e$message); NULL })
    if (is.null(raw)||nrow(raw)==0) { message("    chr",chrom,": no data"); next }
    chr_idx    <- shared_idx[V1==chrom]
    bin_offset <- min(chr_idx$V4) - 1L
    bin1  <- as.integer(raw$x / resolution) + 1L + bin_offset
    bin2  <- as.integer(raw$y / resolution) + 1L + bin_offset
    valid <- bin1>=1L & bin2>=1L &
      bin1<=max(shared_idx$V4) & bin2<=max(shared_idx$V4)
    if (!any(valid)) { message("    chr",chrom,": no valid bins"); next }
    mat_list[[i]] <- data.table::data.table(
      V1 = as.integer(bin1[valid]),
      V2 = as.integer(bin2[valid]),
      V3 = as.numeric(raw$counts[valid])
    )
    total_contacts <- total_contacts + sum(valid)
    message("    chr",chrom,": ",sum(valid)," contacts")
  }
  MAT <- data.table::rbindlist(Filter(Negate(is.null), mat_list))
  MAT <- MAT[V1 <= V2]
  data.table::setkey(MAT, V1, V2)
  if (!is.null(scale_bp) && total_contacts > 0)
    MAT[, V3 := V3 * (scale_bp / total_contacts)]
  message("  Total: ", nrow(MAT), " contact pairs")
  exp <- list(
    MAT  = MAT,
    IDX  = data.table::copy(shared_idx),
    RES  = as.integer(resolution),
    ATTR = list(
      samplename = sample_name,
      colour     = colour,
      resolution = as.integer(resolution),
      scale_cis  = 1L,
      balanced   = FALSE
    )
  )
  class(exp) <- c("contacts", "list")
  attr(exp, "res")        <- as.integer(resolution)
  attr(exp, "samplename") <- sample_name
  exp
}

# Load all conditions — shared IDX means no sync_indices needed
message("\nLoading contacts (TAD chromosomes only)...")
exp_list <- lapply(names(CONDITIONS), function(n) {
  cond <- CONDITIONS[[n]]
  message(cond$label,":")
  load_contacts_selective(
    hic_file    = cond$hic,
    shared_idx  = SHARED_IDX,
    resolution  = RESOLUTION,
    norm        = STRAW_NORM,
    sample_name = cond$label,
    colour      = cond$color,
    scale_bp    = cond$scale
  )
})
names(exp_list) <- names(CONDITIONS)

# Sanity check
invisible(lapply(names(exp_list),function(n) {
  e <- exp_list[[n]]
  message(CONDITIONS[[n]]$label,
          "  bins:",nrow(e$IDX),
          "  contacts:",nrow(e$MAT),
          "  res:",e$RES)
}))

# ══════════════════════════════════════════════════════════════════════════════
# 6. BED converter for GENOVA::ATA
#    Chromosome names in BED must match IDX — strip chr prefix to match .hic
# ══════════════════════════════════════════════════════════════════════════════
tads_to_bed <- function(tads_df, bin=NULL, de=NULL) {
  df <- tads_df
  if (!is.null(bin)) df <- df %>% filter(expr_bin == bin)
  if (!is.null(de))  df <- df %>% filter(de_status == de)
  df %>%
    dplyr::select(chr=chr1, start=x1, end=x2) %>%
    dplyr::mutate(chr = gsub("^chr","", as.character(chr)))
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. ATA runner
# ══════════════════════════════════════════════════════════════════════════════
run_ata <- function(tads_df, exp_list, bins=BINS, min_tads=10, size=ATA_SIZE) {
  setNames(lapply(bins, function(bin) {
    bed <- tads_to_bed(tads_df, bin=bin)
    if (nrow(bed) < min_tads) {
      message("Skipping ",bin," — only ",nrow(bed)," TADs"); return(NULL)
    }
    message("ATA — ",bin," (",nrow(bed)," TADs)...")
    tryCatch({
      res <- GENOVA::ATA(explist    = unname(exp_list),
                         bed        = bed,
                         dist_thres = c(TAD_MIN_KB*1000, TAD_MAX_KB*1000),
                         size       = size)
      message("  Done"); res
    }, error=function(e) { message("  Error: ",e$message); NULL })
  }), bins)
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Plot grid
# ══════════════════════════════════════════════════════════════════════════════
plot_ata_grid <- function(ata_results, bins=BINS, title="ATA") {
  get_mat <- function(res, idx) {
    if (is.null(res)) return(NULL)
    tryCatch(res$signal[,,idx], error=function(e) NULL)
  }
  ata_score_mat <- function(mat, cs=10) {
    if (is.null(mat)) return(NA_real_)
    n  <- nrow(mat); i1 <- floor(n/3); i2 <- ceiling(2*n/3)
    inner <- mat[i1:i2,i1:i2]; sz <- nrow(inner)
    for (d in -2:2) for (i in seq_len(sz)) { j<-i-d; if(j>=1&&j<=sz) inner[i,j]<-NA }
    fg <- mean(inner,na.rm=TRUE)
    bg <- mean(c(mat[1:cs,1:cs],mat[1:cs,(n-cs+1):n],
                 mat[(n-cs+1):n,1:cs],mat[(n-cs+1):n,(n-cs+1):n]),na.rm=TRUE)
    if (is.na(bg)||bg==0) return(NA_real_); fg/bg
  }
  mask_diag <- function(mat,d=2) {
    m <- matrix(TRUE,nrow(mat),ncol(mat))
    for (i in -d:d) m[row(m)==col(m)-i] <- FALSE
    mat[!m] <- NA; mat
  }
  valid_bins <- bins[!sapply(ata_results[bins], is.null)]
  all_ind <- unlist(lapply(valid_bins,function(b) {
    mat <- ata_results[[b]]$signal; n <- nrow(mat[,,1])
    m   <- matrix(TRUE,n,n); for (d in -2:2) m[row(m)==col(m)-d] <- FALSE
    as.vector(mat[rep(m,dim(mat)[3])])
  }))
  all_ind <- all_ind[!is.na(all_ind)&is.finite(all_ind)]
  g_max   <- quantile(all_ind,0.98,na.rm=TRUE)
  g_min   <- quantile(all_ind,0.02,na.rm=TRUE)
  all_diff <- unlist(lapply(valid_bins,function(b) {
    if (dim(ata_results[[b]]$signal)[3]<2) return(NULL)
    wt <- ata_results[[b]]$signal[,,1]
    lapply(2:dim(ata_results[[b]]$signal)[3],function(i)
      as.vector(mask_diag(ata_results[[b]]$signal[,,i]-wt)))
  }))
  all_diff  <- unlist(all_diff); all_diff <- all_diff[!is.na(all_diff)&is.finite(all_diff)]
  d_abs_max <- quantile(abs(all_diff),0.98,na.rm=TRUE)
  
  n_per_bin <- setNames(sapply(bins,function(b) {
    if (is.null(ata_results[[b]])) return(NA_integer_)
    tryCatch(dim(ata_results[[b]]$signal_raw[[1]])[1],error=function(e) NA_integer_)
  }),bins)
  
  make_ind <- function(mat,score=NULL,n_tads=NULL) {
    if (is.null(mat)) return(ggplot()+theme_void())
    n  <- nrow(mat)
    df <- expand.grid(row=1:n,col=1:n) %>%
      dplyr::mutate(value=as.vector(mat),x=col,y=-row)
    p <- ggplot(df,aes(x=x,y=y,fill=pmin(pmax(value,g_min),g_max)))+
      geom_tile()+
      scale_fill_viridis_c(option="viridis",limits=c(g_min,g_max),
                           name="Contacts",oob=scales::squish)+
      coord_fixed()+theme_void()+
      theme(legend.position="none",plot.margin=margin(2,2,2,2))
    if (!is.null(score)&&!is.na(score))
      p <- p+annotate("text",x=3,y=-3,label=paste0("ATA=",round(score,2)),
                      size=2.5,hjust=0,fontface="bold",color="white")
    if (!is.null(n_tads)&&!is.na(n_tads))
      p <- p+annotate("text",x=n-2,y=-n+2,label=paste0("n=",n_tads),
                      size=2.3,hjust=1,color="gray80")
    p
  }
  make_diff <- function(diff_mat,score_diff=NULL,n_tads=NULL) {
    if (is.null(diff_mat)) return(ggplot()+theme_void())
    n  <- nrow(diff_mat)
    df <- expand.grid(row=1:n,col=1:n) %>%
      dplyr::mutate(value=as.vector(diff_mat),x=col,y=-row)
    p <- ggplot(df,aes(x=x,y=y,fill=pmin(pmax(value,-d_abs_max),d_abs_max)))+
      geom_tile()+
      scale_fill_gradientn(colors=RColorBrewer::brewer.pal(11,"PRGn"),
                           limits=c(-d_abs_max,d_abs_max),name="\u0394")+
      coord_fixed()+theme_void()+
      theme(legend.position="none",plot.margin=margin(2,2,2,2))
    if (!is.null(score_diff)&&!is.na(score_diff))
      p <- p+annotate("text",x=3,y=-3,label=paste0("\u0394ATA=",round(score_diff,2)),
                      size=2.5,hjust=0,fontface="bold",color="black")
    if (!is.null(n_tads)&&!is.na(n_tads))
      p <- p+annotate("text",x=n-2,y=-n+2,label=paste0("n=",n_tads),
                      size=2.3,hjust=1,color="gray40")
    p
  }
  
  col_headers <- lapply(bins,function(b)
    ggplot()+annotate("text",x=0.5,y=0.65,label=b,
                      color=BIN_COLORS[b],size=3.5,fontface="bold")+theme_void())
  all_rows <- list(header=wrap_plots(col_headers,nrow=1))
  for (i in seq_along(CONDITIONS)) {
    n   <- names(CONDITIONS)[i]
    row <- lapply(bins,function(b) {
      mat <- get_mat(ata_results[[b]],i)
      make_ind(mat,score=ata_score_mat(mat),n_tads=n_per_bin[[b]])
    })
    lbl <- ggplot()+
      annotate("text",x=0.5,y=0.5,label=CONDITIONS[[n]]$label,
               color=CONDITIONS[[n]]$color,size=3.5,fontface="bold",angle=90)+theme_void()
    all_rows[[paste0("cond_",n)]] <- lbl+wrap_plots(row,nrow=1)+plot_layout(widths=c(0.08,1))
  }
  for (n in names(CONDITIONS)[names(CONDITIONS)!="wt"]) {
    ci  <- which(names(CONDITIONS)==n)
    row <- lapply(bins,function(b) {
      wt_mat <- get_mat(ata_results[[b]],1)
      tr_mat <- get_mat(ata_results[[b]],ci)
      if (is.null(wt_mat)||is.null(tr_mat)) return(make_diff(NULL))
      diff   <- tr_mat-wt_mat
      s_diff <- ata_score_mat(tr_mat)-ata_score_mat(wt_mat)
      make_diff(diff,score_diff=s_diff,n_tads=n_per_bin[[b]])
    })
    lbl <- ggplot()+
      annotate("text",x=0.5,y=0.5,label=paste0(CONDITIONS[[n]]$label,"\n\u2212 WT"),
               color=CONDITIONS[[n]]$color,size=3,fontface="bold",angle=90)+theme_void()
    all_rows[[paste0("diff_",n)]] <- lbl+wrap_plots(row,nrow=1)+plot_layout(widths=c(0.08,1))
  }
  n_cond <- length(CONDITIONS); n_diff <- n_cond-1
  main   <- wrap_plots(all_rows,ncol=1,heights=c(0.15,rep(1,n_cond),rep(1,n_diff)))
  make_leg <- function(limits,palette,direction,ttl) {
    legend_df <- data.frame(x=1:100,y=1,fill=seq(limits[1],limits[2],length.out=100))
    cowplot::get_legend(
      ggplot(legend_df,aes(x=x,y=y,fill=fill))+geom_tile()+
        scale_fill_distiller(palette=palette,direction=direction,limits=limits,name=ttl)+
        theme_void()+theme(legend.position="right",legend.key.height=unit(1.2,"cm"),
                           legend.key.width=unit(0.3,"cm"),legend.text=element_text(size=8),
                           legend.title=element_text(size=9)))
  }
  make_viridis_leg <- function(limits, ttl) {
    legend_df <- data.frame(x=1:100,y=1,fill=seq(limits[1],limits[2],length.out=100))
    cowplot::get_legend(
      ggplot(legend_df,aes(x=x,y=y,fill=fill))+geom_tile()+
        scale_fill_viridis_c(option="viridis",limits=limits,
                             name=ttl,oob=scales::squish)+
        theme_void()+theme(legend.position="right",legend.key.height=unit(1.2,"cm"),
                           legend.key.width=unit(0.3,"cm"),legend.text=element_text(size=8),
                           legend.title=element_text(size=9)))
  }
  legends <- cowplot::plot_grid(
    make_viridis_leg(c(g_min,g_max),"Contacts"),
    make_leg(c(-d_abs_max,d_abs_max),"PRGn",-1,"\u0394 Contacts\n(vs WT)"),ncol=1)
  composite <- cowplot::plot_grid(main,legends,ncol=2,rel_widths=c(1,0.08))
  title_grob <- cowplot::ggdraw()+
    cowplot::draw_label(title,fontface="bold",size=13,x=0,hjust=0)+
    cowplot::draw_label(
      paste0("GENOVA ATA  |  selective chr load  |  ",ATA_SIZE,"\u00d7",ATA_SIZE,
             " bins  |  ",STRAW_NORM,"  |  ",DS_GROUP,"  |  \u0394 = PRGn  |  n = TADs"),
      size=8,color="gray40",x=0,hjust=0,y=0.2)
  final <- cowplot::plot_grid(title_grob,composite,ncol=1,rel_heights=c(0.06,1))
  print(final); invisible(final)
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Run analyses
# ══════════════════════════════════════════════════════════════════════════════
# expr.all.tads.ata
message("\n══ ATA — all WT TADs ═══════════════════════════════════════════════") 
ata_all <- run_ata(tads_anchor$wt, exp_list)
plot_ata_grid(ata_all, title=paste0("ATA — all WT TADs  |  ",DS_GROUP))

#expr.wt.tads.noctcf.ata
message("\n══ ATA — WT non-CTCF ══════════════════════════════════════════════")
ata_wt_non_ctcf <- run_ata(wt_split$non_ctcf, exp_list)
plot_ata_grid(ata_wt_non_ctcf, title=paste0("ATA — WT non-CTCF TADs  |  ",DS_GROUP))

message("\n══ ATA — pol2 non-CTCF ════════════════════════════════════════════")
ata_pol2_non_ctcf <- run_ata(pol2_split$non_ctcf, exp_list)
plot_ata_grid(ata_pol2_non_ctcf, title=paste0("ATA — pol2 non-CTCF TADs  |  ",DS_GROUP))

message("\n══ ATA — WT CTCF ══════════════════════════════════════════════════")
ata_wt_ctcf <- run_ata(wt_split$ctcf, exp_list)
plot_ata_grid(ata_wt_ctcf, title=paste0("ATA — WT CTCF TADs  |  ",DS_GROUP))

message("\n══ ATA — pol2 CTCF ════════════════════════════════════════════════")
ata_pol2_ctcf <- run_ata(pol2_split$ctcf, exp_list)
plot_ata_grid(ata_pol2_ctcf, title=paste0("ATA — pol2 CTCF TADs  |  ",DS_GROUP))

message("\n══ ATA — WT pure CTCF ═════════════════════════════════════════════")
ata_wt_pure <- run_ata(wt_pure_ctcf, exp_list)
plot_ata_grid(ata_wt_pure, title=paste0("ATA — WT pure CTCF TADs  |  ",DS_GROUP))

message("\n══ ATA — pol2 pure CTCF ═══════════════════════════════════════════")
ata_pol2_pure <- run_ata(pol2_pure_ctcf, exp_list)
plot_ata_grid(ata_pol2_pure, title=paste0("ATA — pol2 pure CTCF TADs  |  ",DS_GROUP))

# ══════════════════════════════════════════════════════════════════════════════
# Final summary plots
# ══════════════════════════════════════════════════════════════════════════════
ata_score_mat <- function(mat, cs=10) {
  if (is.null(mat)) return(NA_real_)
  n  <- nrow(mat); i1 <- floor(n/3); i2 <- ceiling(2*n/3)
  inner <- mat[i1:i2,i1:i2]; sz <- nrow(inner)
  for (d in -2:2) for (i in seq_len(sz)) { j<-i-d; if(j>=1&&j<=sz) inner[i,j]<-NA }
  fg <- mean(inner,na.rm=TRUE)
  bg <- mean(c(mat[1:cs,1:cs],mat[1:cs,(n-cs+1):n],
               mat[(n-cs+1):n,1:cs],mat[(n-cs+1):n,(n-cs+1):n]),na.rm=TRUE)
  if (is.na(bg)||bg==0) return(NA_real_); fg/bg
}

make_ind_panel <- function(mat, g_min, g_max, label, n_tads=NA) {
  if (is.null(mat)) return(ggplot()+theme_void())
  n  <- nrow(mat); sc <- ata_score_mat(mat)
  df <- expand.grid(row=1:n,col=1:n) %>%
    dplyr::mutate(value=as.vector(mat),x=col,y=-row)
  p <- ggplot(df,aes(x=x,y=y,fill=pmin(pmax(value,g_min),g_max)))+
    geom_tile()+
    scale_fill_viridis_c(option="viridis",limits=c(g_min,g_max),
                         name="Contacts",oob=scales::squish)+
    coord_fixed()+labs(title=label)+theme_void()+
    theme(legend.position="none",
          plot.title=element_text(size=9,hjust=0.5,face="bold",color="gray20"),
          plot.margin=margin(4,4,4,4))
  if (!is.na(sc))
    p <- p+annotate("text",x=3,y=-3,label=paste0("ATA=",round(sc,2)),
                    size=2.5,hjust=0,fontface="bold",color="white")
  if (!is.na(n_tads))
    p <- p+annotate("text",x=n-2,y=-n+2,label=paste0("n=",n_tads),
                    size=2.3,hjust=1,color="gray80")
  p
}

make_diff_panel_final <- function(wt_mat, tr_mat, d_max, label, n_tads=NA) {
  if (is.null(wt_mat)||is.null(tr_mat)) return(ggplot()+theme_void())
  diff <- tr_mat - wt_mat
  n    <- nrow(diff); sc <- ata_score_mat(diff)
  df   <- expand.grid(row=1:n,col=1:n) %>%
    dplyr::mutate(value=as.vector(diff),x=col,y=-row)
  p <- ggplot(df,aes(x=x,y=y,fill=pmin(pmax(value,-d_max),d_max)))+
    geom_tile()+
    scale_fill_gradientn(colors=RColorBrewer::brewer.pal(11,"PRGn"),
                         limits=c(-d_max,d_max),name="\u0394 Contacts",
                         oob=scales::squish)+
    coord_fixed()+labs(title=label)+theme_void()+
    theme(legend.position="none",
          plot.title=element_text(size=9,hjust=0.5,face="bold",color="gray20"),
          plot.margin=margin(4,4,4,4))
  if (!is.na(sc))
    p <- p+annotate("text",x=3,y=-3,label=paste0("\u0394ATA=",round(sc,2)),
                    size=2.5,hjust=0,fontface="bold",color="black")
  if (!is.na(n_tads))
    p <- p+annotate("text",x=n-2,y=-n+2,label=paste0("n=",n_tads),
                    size=2.3,hjust=1,color="gray40")
  p
}

plot_ata_summary <- function(ata_res, title,
                             label_cond1="WT",
                             label_cond2="pol2") {
  wt_hi  <- tryCatch(ata_res$High$signal[,,1], error=function(e) NULL)
  tr_hi  <- tryCatch(ata_res$High$signal[,,2], error=function(e) NULL)
  wt_lo  <- tryCatch(ata_res$Low$signal[,,1],  error=function(e) NULL)
  tr_lo  <- tryCatch(ata_res$Low$signal[,,2],  error=function(e) NULL)
  n_hi   <- tryCatch(dim(ata_res$High$signal_raw[[1]])[1], error=function(e) NA)
  n_lo   <- tryCatch(dim(ata_res$Low$signal_raw[[1]])[1],  error=function(e) NA)
  
  all_ind  <- c(as.vector(wt_hi),as.vector(tr_hi),
                as.vector(wt_lo),as.vector(tr_lo))
  all_ind  <- all_ind[!is.na(all_ind)&is.finite(all_ind)]
  g_min    <- quantile(all_ind,0.02,na.rm=TRUE)
  g_max    <- quantile(all_ind,0.98,na.rm=TRUE)
  all_diff <- c(as.vector(tr_hi-wt_hi),as.vector(tr_lo-wt_lo))
  all_diff <- all_diff[!is.na(all_diff)&is.finite(all_diff)]
  d_max    <- quantile(abs(all_diff),0.98,na.rm=TRUE)
  
  p1 <- make_ind_panel(wt_hi,g_min,g_max,paste0(label_cond1,"  |  High"),n_tads=n_hi)
  p2 <- make_ind_panel(tr_hi,g_min,g_max,paste0(label_cond2,"  |  High"),n_tads=n_hi)
  p3 <- make_diff_panel_final(wt_hi,tr_hi,d_max,paste0(label_cond2,"\u2212",label_cond1,"  |  High"),n_tads=n_hi)
  p4 <- make_ind_panel(wt_lo,g_min,g_max,paste0(label_cond1,"  |  Low"), n_tads=n_lo)
  p5 <- make_ind_panel(tr_lo,g_min,g_max,paste0(label_cond2,"  |  Low"), n_tads=n_lo)
  p6 <- make_diff_panel_final(wt_lo,tr_lo,d_max,paste0(label_cond2,"\u2212",label_cond1,"  |  Low"), n_tads=n_lo)
  
  main <- patchwork::wrap_plots(p1,p2,p3,p4,p5,p6,ncol=3)
  
  leg_ind <- cowplot::get_legend(
    ggplot(data.frame(x=1:100,y=1,fill=seq(g_min,g_max,length.out=100)),
           aes(x=x,y=y,fill=fill))+geom_tile()+
      scale_fill_viridis_c(option="viridis",limits=c(g_min,g_max),
                           name="Contacts",oob=scales::squish)+
      theme_void()+
      theme(legend.position="right",legend.key.height=unit(1.2,"cm"),
            legend.key.width=unit(0.3,"cm"),legend.text=element_text(size=8),
            legend.title=element_text(size=9)))
  leg_diff <- cowplot::get_legend(
    ggplot(data.frame(x=1:100,y=1,fill=seq(-d_max,d_max,length.out=100)),
           aes(x=x,y=y,fill=fill))+geom_tile()+
      scale_fill_gradientn(colors=RColorBrewer::brewer.pal(11,"PRGn"),
                           limits=c(-d_max,d_max),name="\u0394 Contacts")+
      theme_void()+
      theme(legend.position="right",legend.key.height=unit(1.2,"cm"),
            legend.key.width=unit(0.3,"cm"),legend.text=element_text(size=8),
            legend.title=element_text(size=9)))
  
  legs      <- cowplot::plot_grid(leg_ind,leg_diff,ncol=1)
  composite <- cowplot::plot_grid(main,legs,ncol=2,rel_widths=c(1,0.08))
  title_grob <- cowplot::ggdraw()+
    cowplot::draw_label(title,fontface="bold",size=13,x=0,hjust=0)
  final <- cowplot::plot_grid(title_grob,composite,ncol=1,rel_heights=c(0.06,1))
  print(final); invisible(final)
}

# ── Plot 1: WT CTCF positions scored in both conditions ─────────────────────── wt.ctcf.ata.finalplot
plot_ata_summary(ata_wt_ctcf,
                 title=paste0("ATA — CTCF TADs (WT positions)  |  ",DS_GROUP),
                 label_cond1="WT",
                 label_cond2="pol2 at WT pos")

# ── Plot 2: pol2 CTCF positions scored in both conditions ─────────────────────
plot_ata_summary(ata_pol2_ctcf,
                 title=paste0("ATA — CTCF TADs (pol2 positions)  |  ",DS_GROUP),
                 label_cond1="WT at pol2 pos",
                 label_cond2="pol2")

# ── Plot 3: WT non-CTCF positions scored in both conditions ───────────────────
plot_ata_summary(ata_wt_non_ctcf,
                 title=paste0("ATA — non-CTCF TADs (WT positions)  |  ",DS_GROUP),
                 label_cond1="WT",
                 label_cond2="pol2 at WT pos")

# ── Plot 4: pol2 non-CTCF positions scored in both conditions ─────────────────
plot_ata_summary(ata_pol2_non_ctcf,
                 title=paste0("ATA — non-CTCF TADs (pol2 positions)  |  ",DS_GROUP),
                 label_cond1="WT at pol2 pos",
                 label_cond2="pol2")

library(ggpubr)

get_pertad_scores <- function(ata_res, bins=BINS) {
  bind_rows(lapply(bins, function(b) {
    if (is.null(ata_res[[b]])) return(NULL)
    quantify(ata_res[[b]])$per_TAD %>% dplyr::mutate(expr_bin=b)
  })) %>%
    dplyr::mutate(
      expr_bin = factor(expr_bin, levels=bins),
      sample   = factor(sample, levels=sapply(CONDITIONS, `[[`, "label"))
    )
}

cond_colors <- setNames(sapply(CONDITIONS, `[[`, "color"),
                        sapply(CONDITIONS, `[[`, "label"))

make_score_box <- function(scores, subtitle) {
  # Compute Wilcoxon p-values per expression bin
  lvls <- levels(scores$sample)
  pvals <- bind_rows(lapply(levels(scores$expr_bin), function(b) {
    g1 <- scores %>% filter(expr_bin==b, sample==lvls[1]) %>% pull(foldchange)
    g2 <- scores %>% filter(expr_bin==b, sample==lvls[2]) %>% pull(foldchange)
    p  <- wilcox.test(g1, g2)$p.value
    sig <- dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns")
    data.frame(expr_bin=b, p=p, sig=sig,
               y=max(c(g1,g2), na.rm=TRUE)*1.08)
  })) %>% dplyr::mutate(expr_bin=factor(expr_bin, levels=levels(scores$expr_bin)))
  
  ggplot(scores, aes(x=expr_bin, y=foldchange, fill=sample)) +
    geom_boxplot(outlier.size=0.5, outlier.alpha=0.3,
                 position=position_dodge(0.8), width=0.6,
                 na.rm=TRUE) +
    geom_text(data=pvals, aes(x=expr_bin, y=y, label=sig),
              inherit.aes=FALSE, size=4, vjust=0) +
    scale_fill_manual(values=cond_colors, name=NULL) +
    scale_y_continuous(name="ATA score (foreground/background)",
                       expand=expansion(mult=c(0.05,0.2))) +
    labs(subtitle=subtitle, x=NULL) +
    theme_classic(base_size=11) +
    theme(plot.subtitle=element_text(face="bold",size=10,color="gray20"),
          legend.position="none")
}

p1 <- make_score_box(get_pertad_scores(ata_wt_ctcf),     "WT TAD positions  |  CTCF")
p2 <- make_score_box(get_pertad_scores(ata_pol2_ctcf),   "pol2 TAD positions  |  CTCF")
p3 <- make_score_box(get_pertad_scores(ata_wt_non_ctcf), "WT TAD positions  |  non-CTCF")
p4 <- make_score_box(get_pertad_scores(ata_pol2_non_ctcf),"pol2 TAD positions  |  non-CTCF")

shared_legend <- cowplot::get_legend(
  ggplot(get_pertad_scores(ata_wt_ctcf),
         aes(x=expr_bin, y=foldchange, fill=sample)) +
    geom_boxplot() +
    scale_fill_manual(values=cond_colors, name=NULL) +
    theme_classic() +
    theme(legend.position="bottom",
          legend.text=element_text(size=11),
          legend.key.size=unit(0.8,"cm")))

main <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=c("A","B","C","D"))
title_grob <- cowplot::ggdraw() +
  cowplot::draw_label("ATA scores by expression bin and CTCF status",
                      fontface="bold", size=13, x=0, hjust=0)
final <- cowplot::plot_grid(title_grob,
                            main,
                            shared_legend,
                            ncol=1, rel_heights=c(0.05, 1, 0.06))
print(final)

# Mirror of plot_ata_summary heatmaps — per-TAD score boxplots
scores_wt_ctcf   <- get_pertad_scores(ata_wt_ctcf)   %>% dplyr::mutate(tad_set="WT positions")
scores_pol2_ctcf <- get_pertad_scores(ata_pol2_ctcf) %>% dplyr::mutate(tad_set="pol2 positions")

p_box_ctcf_wt   <- make_score_box(scores_wt_ctcf,   "WT TAD positions  |  CTCF")
p_box_ctcf_pol2 <- make_score_box(scores_pol2_ctcf, "pol2 TAD positions  |  CTCF")

shared_legend <- cowplot::get_legend(
  ggplot(scores_wt_ctcf, aes(x=expr_bin, y=foldchange, fill=sample)) +
    geom_boxplot() +
    scale_fill_manual(values=cond_colors, name=NULL) +
    theme_classic() +
    theme(legend.position="bottom",
          legend.text=element_text(size=11),
          legend.key.size=unit(0.8,"cm")))

main <- cowplot::plot_grid(p_box_ctcf_wt, p_box_ctcf_pol2, ncol=2, labels=c("A","B"))
title_grob <- cowplot::ggdraw() +
  cowplot::draw_label("ATA scores — CTCF TADs  |  WT vs pol2 TAD positions",
                      fontface="bold", size=13, x=0, hjust=0)
final_ctcf_box <- cowplot::plot_grid(title_grob, main, shared_legend,
                                     ncol=1, rel_heights=c(0.05,1,0.06))
print(final_ctcf_box)
