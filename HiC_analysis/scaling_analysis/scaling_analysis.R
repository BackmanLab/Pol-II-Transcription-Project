##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Hi-C contact scaling by E/I ratio. DS_GROUP switch for
# downsampled/unsampled Hi-C files. Contacts extracted via strawr
# (SCALE norm, 5kb resolution).
#
# Faceted analyses: E/I bins | Expression | Pol2 occupancy | Intergenic E/I
#
# Genes (>=100kb) stratified by exon-to-intron (E/I) ratio into 100
# quantile bins, or by expression into tertiles; per-bin contact
# probability P(s) pooled across genes, S = slope of
# log10(P(s)) vs. log10(distance). Sliding-window S(d) by TPM
# quartile for the distance-resolved panel. Short genes (<100kb)
# use a separate tandem-cluster branch (exon/intergenic ratio
# rather than exon/intron) -- see tandem_gene_analysis.R for the
# diagnostic behind the 50kb cluster-gap threshold.
#
# Figure output: Fig. 4F-H (confirmed).
##################################################################
library(strawr)
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggdist)
library(cowplot)
library(GenomicRanges)
library(GenomicFeatures)
library(plyranges)
library(parallel)
library(rstatix)
library(ggpubr)

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════
DS_GROUP <- "ds_wt_pol2"   # ← "ds_wt_pol2" | "ds_all" | "unsampled"

BASE_NEW  <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
BASE_OLD  <- "/projects/p32171/juicer/work/112123_HiC/juicer_analysis"

if (DS_GROUP == "ds_wt_pol2") {
  CONDITIONS <- list(
    wt   = list(label="WT_HCT116_CTRL", display="Wildtype",   color="#88CCEE",
                hic=file.path(BASE_NEW,"mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic")),
    pol2 = list(label="Pol2_6hrs_Aux",  display="POLR2A deg", color="#CC6677",
                hic=file.path(BASE_NEW,"mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"))
  )
} else if (DS_GROUP == "ds_all") {
  CONDITIONS <- list(
    wt   = list(label="WT_HCT116_CTRL", display="Wildtype",   color="#88CCEE",
                hic=file.path(BASE_NEW,"mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic")),
    pol2 = list(label="Pol2_6hrs_Aux",  display="POLR2A deg", color="#CC6677",
                hic=file.path(BASE_NEW,"mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"))
  )
} else if (DS_GROUP == "unsampled") {
  CONDITIONS <- list(
    wt   = list(label="WT_HCT116_CTRL", display="Wildtype",   color="#88CCEE",
                hic=file.path(BASE_OLD,"WT_HCT116_CTRL/mega/aligned/inter_30.hic")),
    pol2 = list(label="Pol2_6hrs_Aux",  display="POLR2A deg", color="#CC6677",
                hic=file.path(BASE_OLD,"Pol2_6hrs_Aux/mega/aligned/inter_30.hic"))
  )
} else { stop("Unknown DS_GROUP: '", DS_GROUP, "'") }
message("DS_GROUP: ", DS_GROUP)

COND_LEVELS  <- sapply(CONDITIONS, `[[`, "display")
COND_COLORS  <- setNames(sapply(CONDITIONS, `[[`, "color"), COND_LEVELS)

GTF_FILE   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
EXPR_FILE  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/WT_averaged_TPM.tsv"
DESEQ_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/Conditions_Pol2_6hrs_vs_CTRL_0hrs.genes.csv"
POL2_PEAK  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/peaks/pol2.narrowPeak"

RDS_DIR    <- file.path(BASE_NEW, "ei_scaling/rds", DS_GROUP)
dir.create(RDS_DIR, showWarnings=FALSE, recursive=TRUE)

RESOLUTION <- 5000L
STRAW_NORM <- "SCALE"
EI_BINS    <- 100L
MIN_GENE_SIZE <- 100000L
MAX_GAP    <- 50000L    # max inter-gene gap within a tandem cluster
MIN_SPAN   <- 100000L   # min cluster span to be kept
MIN_GAP    <- 250L
STD_CHROMS <- paste0("chr", c(1:22, "X"))
N_CORES    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset="2"))
message("Using ", N_CORES, " cores")

COND_DISPLAY <- setNames(sapply(CONDITIONS, `[[`, "display"),
                         sapply(CONDITIONS, `[[`, "label"))

# ══════════════════════════════════════════════════════════════════════════════
# 1. Build TxDb and compute EI
# ══════════════════════════════════════════════════════════════════════════════
ei_rds <- file.path(RDS_DIR, "ei_table.rds")
if (file.exists(ei_rds)) {
  message("Loading cached EI table...")
  ei_df <- readRDS(ei_rds)
} else {
  message("Building TxDb from GTF...")
  txdb     <- makeTxDbFromGFF(GTF_FILE, format="gtf")
  exons_gr <- exons(txdb) %>% GenomicRanges::reduce()
  genes_gr <- genes(txdb)
  seqlevelsStyle(genes_gr) <- "UCSC"
  seqlevelsStyle(exons_gr) <- "UCSC"
  genes_gr <- keepSeqlevels(genes_gr,
                            intersect(seqlevels(genes_gr), STD_CHROMS),
                            pruning.mode="coarse")
  exons_gr <- keepSeqlevels(exons_gr,
                            intersect(seqlevels(exons_gr), STD_CHROMS),
                            pruning.mode="coarse")
  
  message("Computing exon bp per gene...")
  genes_gr$gene_id <- sub("\\..*","", names(genes_gr))
  
  # Per-gene exon coverage
  hits     <- findOverlaps(genes_gr, exons_gr)
  overlaps <- GenomicRanges::pintersect(genes_gr[queryHits(hits)],
                                        exons_gr[subjectHits(hits)])
  exon_bp_per_gene <- tapply(width(overlaps), queryHits(hits), sum)
  exon_bp_vec <- rep(0L, length(genes_gr))
  exon_bp_vec[as.integer(names(exon_bp_per_gene))] <- as.integer(exon_bp_per_gene)
  
  ei_df <- data.frame(
    gene_id  = genes_gr$gene_id,
    chr      = as.character(seqnames(genes_gr)),
    start    = start(genes_gr),
    end      = end(genes_gr),
    gene_bp  = width(genes_gr),
    exon_bp  = exon_bp_vec,
    stringsAsFactors=FALSE
  ) %>%
    dplyr::mutate(
      intron_bp = gene_bp - exon_bp,
      EI        = exon_bp / pmax(intron_bp, 1),
      size      = end - start
    ) %>%
    dplyr::filter(
      exon_bp > 0,
      !is.infinite(EI),
      EI > 0,
      EI != 1
    )
  
  saveRDS(ei_df, ei_rds)
  message("Saved EI table: ", nrow(ei_df), " genes")
}
message("EI table: ", nrow(ei_df), " genes")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Expression data
# ══════════════════════════════════════════════════════════════════════════════
message("Loading expression data...")
expr <- read_tsv(EXPR_FILE, show_col_types=FALSE) %>%
  dplyr::select(gene_id, TPM) %>%
  dplyr::mutate(gene_id=sub("\\..*","",gene_id))

expressed    <- expr %>% filter(TPM>0) %>% pull(TPM)
tpm_medium   <- quantile(expressed, 1/3, na.rm=TRUE)
tpm_high     <- quantile(expressed, 2/3, na.rm=TRUE)
message("  TPM thresholds: Low<=",round(tpm_medium,2),
        "  Med<=",round(tpm_high,2),
        "  High>",round(tpm_high,2))

expr <- expr %>%
  dplyr::mutate(
    expr_bin = dplyr::case_when(
      is.na(TPM)|TPM==0       ~ "Low",
      TPM <= tpm_medium       ~ "Low",
      TPM <= tpm_high         ~ "Medium",
      TRUE                    ~ "High"
    ),
    expr_bin = factor(expr_bin, levels=c("Low","Medium","High"))
  )

# Join expression onto EI table
ei_df <- ei_df %>%
  left_join(expr, by="gene_id") %>%
  dplyr::mutate(
    TPM      = replace_na(TPM, 0),
    expr_bin = replace_na(as.character(expr_bin), "Low"),
    expr_bin = factor(expr_bin, levels=c("Low","Medium","High"))
  )

# ══════════════════════════════════════════════════════════════════════════════
# 3. Pol2 ChIP occupancy
# ══════════════════════════════════════════════════════════════════════════════
pol2_rds <- file.path(RDS_DIR, "pol2_occupancy.rds")
if (file.exists(pol2_rds)) {
  message("Loading cached Pol2 occupancy...")
  pol2_occ <- readRDS(pol2_rds)
} else {
  message("Computing Pol2 occupancy from ChIP peaks...")
  pol2_peaks <- read.table(POL2_PEAK, header=FALSE, sep="\t",
                           colClasses=c("character","integer","integer",
                                        "character","integer","character",
                                        "numeric","numeric","numeric","integer")) %>%
    dplyr::rename(seqnames=V1,start=V2,end=V3,name=V4,score=V5,
                  strand=V6,signalValue=V7,pValue=V8,qValue=V9,peak=V10) %>%
    dplyr::mutate(start=start+1L,
                  strand=dplyr::if_else(strand==".","*",strand)) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  seqlevelsStyle(pol2_peaks) <- "UCSC"
  pol2_peaks <- keepSeqlevels(pol2_peaks,
                              intersect(seqlevels(pol2_peaks),STD_CHROMS),
                              pruning.mode="coarse")
  
  genes_for_pol2 <- makeGRangesFromDataFrame(ei_df, keep.extra.columns=TRUE)
  hits      <- findOverlaps(genes_for_pol2, pol2_peaks)
  pol2_signal <- tapply(pol2_peaks$signalValue[subjectHits(hits)],
                        queryHits(hits), sum)
  pol2_vec    <- rep(0, nrow(ei_df))
  pol2_vec[as.integer(names(pol2_signal))] <- as.numeric(pol2_signal)
  
  pol2_occ <- data.frame(gene_id=ei_df$gene_id, pol2_signal=pol2_vec,
                         stringsAsFactors=FALSE)
  saveRDS(pol2_occ, pol2_rds)
  message("Saved Pol2 occupancy")
}

# Join Pol2 occupancy and assign bins
ei_df <- ei_df %>%
  left_join(pol2_occ, by="gene_id") %>%
  dplyr::mutate(pol2_signal=replace_na(pol2_signal, 0))

pol2_expressed    <- ei_df %>% filter(pol2_signal>0) %>% pull(pol2_signal)
pol2_medium_thresh <- quantile(pol2_expressed, 1/3, na.rm=TRUE)
pol2_high_thresh   <- quantile(pol2_expressed, 2/3, na.rm=TRUE)
ei_df <- ei_df %>%
  dplyr::mutate(
    pol2_bin = dplyr::case_when(
      pol2_signal==0                  ~ "Low",
      pol2_signal<=pol2_medium_thresh ~ "Low",
      pol2_signal<=pol2_high_thresh   ~ "Medium",
      TRUE                            ~ "High"
    ),
    pol2_bin = factor(pol2_bin, levels=c("Low","Medium","High"))
  )
message("Pol2 bin counts:"); print(table(ei_df$pol2_bin))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Tandem gene cluster analysis
#    - Input: genes < 100kb (single genes >= 100kb handled separately)
#    - Cluster: consecutive genes within MAX_GAP on same chromosome
#    - Full: cluster span >= MIN_SPAN; last gene stays in that cluster
#    - EI_cluster = total exon bp / total intergenic bp within cluster
#    - Minimum 2 genes per cluster
# ══════════════════════════════════════════════════════════════════════════════
cluster_rds <- file.path(RDS_DIR, "tandem_clusters.rds")
if (file.exists(cluster_rds)) {
  message("Loading cached tandem clusters...")
  cluster_summary <- readRDS(cluster_rds)
} else {
  message("Building tandem gene clusters...")
  
  # Pool: genes < 100kb with exon annotation
  small_genes <- ei_df %>%
    dplyr::filter(size < MIN_GENE_SIZE) %>%
    dplyr::arrange(chr, start)
  
  message("  Input genes < 100kb: ", nrow(small_genes))
  
  build_clusters <- function(genes_chr, max_gap=MAX_GAP, min_span=MIN_SPAN) {
    if (nrow(genes_chr) < 2) return(NULL)
    clusters   <- list()
    cluster_id <- 1L
    current    <- list(genes_chr[1, , drop=FALSE])
    
    for (i in 2:nrow(genes_chr)) {
      gene <- genes_chr[i, , drop=FALSE]
      
      if (length(current) == 0) {
        current <- list(gene)
        next
      }
      
      prev          <- current[[length(current)]]
      gap           <- gene$start - prev$end - 1L
      span_if_added <- gene$end - current[[1]]$start + 1L
      
      if (gap <= max_gap) {
        current <- c(current, list(gene))
        if (span_if_added >= min_span) {
          cluster_df <- bind_rows(current) %>%
            dplyr::mutate(cluster_id=cluster_id)
          clusters[[cluster_id]] <- cluster_df
          cluster_id <- cluster_id + 1L
          current    <- list()
        }
      } else {
        current <- list(gene)
      }
    }
    if (length(clusters)==0) return(NULL)
    bind_rows(clusters)
  }
  
  cluster_genes <- bind_rows(Filter(Negate(is.null),
                                    lapply(split(small_genes, small_genes$chr),
                                           build_clusters)))
  
  message("  Gene-cluster assignments: ", nrow(cluster_genes))
  message("  Total clusters: ", length(unique(cluster_genes$cluster_id)))
  
  # Cluster-level summary
  cluster_summary <- cluster_genes %>%
    dplyr::group_by(cluster_id, chr) %>%
    dplyr::summarise(
      n_genes             = dplyr::n(),
      cluster_start       = min(start),
      cluster_end         = max(end),
      cluster_span        = max(end) - min(start) + 1L,
      total_exon_bp       = sum(exon_bp),
      total_gene_bp       = sum(size),
      gene_ids            = paste(gene_id, collapse=";"),
      .groups="drop"
    ) %>%
    dplyr::mutate(
      total_intergenic_bp = cluster_span - total_gene_bp,
      EI_cluster          = dplyr::if_else(
        total_intergenic_bp > 0,
        total_exon_bp / total_intergenic_bp,
        NA_real_),
      coding_density      = total_exon_bp / cluster_span
    ) %>%
    dplyr::filter(n_genes >= 2, total_intergenic_bp > 0)
  
  saveRDS(cluster_summary, cluster_rds)
  message("Saved cluster summary: ", nrow(cluster_summary), " clusters")
}

message("\nCluster summary:")
message("  Total clusters: ", nrow(cluster_summary))
message("  Span range: ",
        round(min(cluster_summary$cluster_span)/1000,1),"kb - ",
        round(max(cluster_summary$cluster_span)/1000,1),"kb")
message("  Genes per cluster:"); print(table(cluster_summary$n_genes))
message("  EI_cluster quantiles:")
print(quantile(cluster_summary$EI_cluster,
               probs=c(0,0.1,0.25,0.5,0.75,0.9,1), na.rm=TRUE))

# ── Cluster diagnostic plots ──────────────────────────────────────────────────
message("\n── Cluster diagnostics ──────────────────────────────────────────────")

# Summary table by n_genes
cluster_tbl <- cluster_summary %>%
  dplyr::group_by(n_genes) %>%
  dplyr::summarise(
    n_clusters        = dplyr::n(),
    median_span_kb    = round(median(cluster_span)/1000, 1),
    median_EI         = round(median(EI_cluster, na.rm=TRUE), 3),
    median_coding_pct = round(median(coding_density)*100, 1),
    .groups="drop"
  )
message("Cluster summary by n_genes:")
print(cluster_tbl)

p_span <- ggplot(cluster_summary, aes(x=cluster_span/1000))+
  geom_histogram(bins=80, fill="#0072B2", color="white",
                 linewidth=0.15, alpha=0.85)+
  geom_vline(xintercept=100, linetype="dashed",
             color="#D55E00", linewidth=0.8)+
  annotate("text", x=105, y=Inf,
           label="100kb", hjust=0, vjust=1.5,
           size=3.5, color="#D55E00", fontface="bold")+
  scale_x_continuous(labels=scales::comma)+
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.05)))+
  labs(x="Cluster span (kb)", y="Number of clusters",
       title="Tandem gene cluster span distribution",
       subtitle=paste0("n = ",scales::comma(nrow(cluster_summary)),
                       " clusters  |  max gap = ",MAX_GAP/1000,"kb  |  ",
                       "genes < 100kb only"))+
  theme_classic(base_size=13)+
  theme(plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=45,hjust=1))
print(p_span)

p_ngenes <- ggplot(cluster_summary, aes(x=n_genes))+
  geom_bar(fill="#0072B2", color="white", alpha=0.85)+
  geom_text(stat="count",
            aes(label=scales::comma(after_stat(count))),
            vjust=-0.4, size=3.5)+
  scale_x_continuous(breaks=seq(2, max(cluster_summary$n_genes), by=1))+
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.12)))+
  labs(x="Number of genes in cluster", y="Number of clusters",
       title="Genes per tandem cluster",
       subtitle=paste0("max gap = ",MAX_GAP/1000,
                       "kb  |  min span = ",MIN_SPAN/1000,"kb"))+
  theme_classic(base_size=13)+
  theme(plot.title=element_text(face="bold"))
print(p_ngenes)

p_ei_clust <- ggplot(cluster_summary %>% filter(!is.na(EI_cluster)),
                     aes(x=EI_cluster))+
  geom_histogram(bins=60, fill="#009E73", color="white",
                 linewidth=0.15, alpha=0.85)+
  geom_vline(xintercept=1, linetype="dashed",
             color="#D55E00", linewidth=0.8)+
  annotate("text", x=1.05, y=Inf,
           label="E/I = 1", hjust=0, vjust=1.5,
           size=3.5, color="#D55E00", fontface="bold")+
  scale_x_log10()+
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.05)))+
  labs(x="Cluster E/I ratio (log scale)", y="Number of clusters",
       title="E/I distribution for tandem clusters",
       subtitle="E = total exon bp  |  I = total intergenic bp within cluster")+
  theme_classic(base_size=13)+
  theme(plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=45,hjust=1))
print(p_ei_clust)

p_span_ngenes <- ggplot(cluster_summary,
                        aes(x=factor(n_genes), y=cluster_span/1000))+
  geom_boxplot(fill="#0072B2", alpha=0.7, outlier.size=0.5,
               outlier.alpha=0.3)+
  geom_hline(yintercept=100, linetype="dashed",
             color="#D55E00", linewidth=0.7)+
  scale_y_continuous(labels=scales::comma)+
  labs(x="Number of genes in cluster", y="Cluster span (kb)",
       title="Cluster span by number of genes",
       subtitle="Dashed line = 100kb minimum span threshold")+
  theme_classic(base_size=13)+
  theme(plot.title=element_text(face="bold"))
print(p_span_ngenes)

# Assign EI bins, coordinates, and expression for contact extraction
cluster_expr <- cluster_summary %>%
  dplyr::mutate(gene_id_list=strsplit(gene_ids, ";")) %>%
  tidyr::unnest(gene_id_list) %>%
  dplyr::left_join(ei_df %>% dplyr::select(gene_id, TPM),
                   by=c("gene_id_list"="gene_id")) %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(mean_TPM=mean(TPM, na.rm=TRUE), .groups="drop")

cluster_expressed  <- cluster_expr %>% filter(mean_TPM>0) %>% pull(mean_TPM)
cluster_tpm_medium <- quantile(cluster_expressed, 1/3, na.rm=TRUE)
cluster_tpm_high   <- quantile(cluster_expressed, 2/3, na.rm=TRUE)

cluster_summary <- cluster_summary %>%
  dplyr::left_join(cluster_expr, by="cluster_id") %>%
  dplyr::mutate(
    EI_bin   = as.integer(ntile(EI_cluster, EI_BINS)),
    start    = cluster_start,
    end      = cluster_end,
    size     = cluster_span,
    mean_TPM = replace_na(mean_TPM, 0),
    expr_bin = dplyr::case_when(
      is.na(mean_TPM)|mean_TPM==0    ~ "Low",
      mean_TPM <= cluster_tpm_medium ~ "Low",
      mean_TPM <= cluster_tpm_high   ~ "Medium",
      TRUE                           ~ "High"
    ),
    expr_bin = factor(expr_bin, levels=c("Low","Medium","High"))
  ) %>%
  dplyr::rename(gene_id = cluster_id)
message("Cluster expression bins:"); print(table(cluster_summary$expr_bin))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Assign EI bins — single genes >= 100kb
# ══════════════════════════════════════════════════════════════════════════════
ei_binned <- ei_df %>%
  dplyr::filter(size >= MIN_GENE_SIZE) %>%
  dplyr::mutate(
    EI_bin   = as.integer(ntile(EI, EI_BINS)),
    EI_bin_3 = as.integer(ntile(EI, 3L))
  )
message("EI binned single genes (size>=100kb): ", nrow(ei_binned), " genes")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Extract contacts from .hic via strawr
#    - Extract once per condition for ALL genes/clusters
#    - Carry all bin columns (EI_bin, expr_bin, pol2_bin) in contacts table
#    - Split by whichever column needed at slope computation time
# ══════════════════════════════════════════════════════════════════════════════
extract_contacts_for_genes <- function(hic_file, genes_df, resolution=RESOLUTION,
                                       norm=STRAW_NORM, n_cores=N_CORES) {
  bin_cols <- intersect(names(genes_df),
                        c("EI_bin","expr_bin","pol2_bin","gene_id"))
  chroms <- intersect(
    unique(genes_df$chr),
    strawr::readHicChroms(hic_file) %>%
      dplyr::filter(name!="ALL", !grepl("_",name)) %>%
      dplyr::pull(name) %>% paste0("chr",.))
  message("  Extracting contacts for ", length(chroms), " chromosomes...")
  
  results <- parallel::mclapply(chroms, function(chrom) {
    chrom_bare <- gsub("^chr","",chrom)
    raw <- tryCatch(
      strawr::straw(norm, hic_file, chrom_bare, chrom_bare, "BP", resolution),
      error=function(e) { message("  Warning ",chrom,": ",e$message); NULL })
    if (is.null(raw)||nrow(raw)==0) return(NULL)
    
    genes_chr <- genes_df %>% dplyr::filter(chr==chrom)
    if (nrow(genes_chr)==0) return(NULL)
    
    raw_dt <- data.table::as.data.table(raw)
    raw_dt <- raw_dt[!is.nan(counts) & !is.na(counts)]
    if (nrow(raw_dt)==0) return(NULL)
    
    gene_contacts <- lapply(seq_len(nrow(genes_chr)), function(i) {
      g   <- genes_chr[i,]
      sub <- raw_dt[x >= g$start & x <= g$end]
      if (nrow(sub)==0) return(NULL)
      sub[, chr := chrom]
      for (col in bin_cols) sub[, (col) := g[[col]]]
      sub
    })
    data.table::rbindlist(Filter(Negate(is.null), gene_contacts))
  }, mc.cores=n_cores, mc.preschedule=FALSE)
  
  data.table::rbindlist(Filter(Negate(is.null), results))
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. Slope and P(s) computation
#    col_sum = total contact probability across ALL bins for the experiment
# ══════════════════════════════════════════════════════════════════════════════
compute_slopes <- function(contacts_dt, filter_col="EI_bin",
                           resolution=RESOLUTION) {
  contacts_dt <- data.table::copy(contacts_dt)
  contacts_dt[, index := abs(x - y)]
  contacts_dt <- contacts_dt[!is.nan(counts) & !is.na(counts)]
  
  bin_list <- split(contacts_dt, contacts_dt[[filter_col]])
  bin_list <- Filter(function(x) nrow(x) > 0, bin_list)
  
  dist_list <- lapply(names(bin_list), function(b) {
    df <- bin_list[[b]]
    # Weighted possible.contacts: mean of max observed contact distance per gene
    gene_max_dist  <- df[, .(max_dist=max(index)), by=gene_id]
    n_gene         <- floor(gene_max_dist$max_dist / resolution)
    total_possible <- mean(n_gene, na.rm=TRUE)
    
    dist <- df[, .(counts=sum(counts)), by=index]
    data.table::setorder(dist, index)
    n <- nrow(dist)
    dist[, possible.contacts := seq(total_possible*resolution, resolution,
                                    by=-resolution)[seq_len(n)]]
    dist[, contact.prob      := counts / possible.contacts]
    dist[, bin_val           := as.integer(b)]
    dist
  })
  
  all_dist <- data.table::rbindlist(dist_list)
  col_sum  <- sum(all_dist$contact.prob, na.rm=TRUE)
  message("  col_sum (experiment-level): ", round(col_sum, 6))
  if (col_sum == 0) {
    warning("  col_sum is 0 — no valid contacts found for this subset")
    return(NULL)
  }
  
  bind_rows(lapply(dist_list, function(dist) {
    b <- dist$bin_val[1]
    dist[, contact.prob := contact.prob / col_sum]
    dist[index==0, index := 1L]
    range <- dist[, .(distance=log10(index), prob=log10(contact.prob))]
    range <- range[is.finite(prob) & distance > 0 & distance < 6]
    if (nrow(range) < 3) return(NULL)
    model <- lm(prob~distance, data=range, na.action="na.exclude")
    r2    <- summary(model)$r.squared
    data.frame(bin_val   =as.integer(b),
               slope     =coef(model)[2],
               intercept =coef(model)[1],
               r2        =r2,
               filter_col=filter_col)
  }))
}

compute_ps_curves <- function(contacts_dt, filter_col="EI_bin",
                              resolution=RESOLUTION,
                              max_distance=6) {
  contacts_dt <- data.table::copy(contacts_dt)
  contacts_dt[, index := abs(x - y)]
  contacts_dt <- contacts_dt[!is.nan(counts) & !is.na(counts)]
  
  bin_list  <- split(contacts_dt, contacts_dt[[filter_col]])
  bin_list  <- Filter(function(x) nrow(x) > 0, bin_list)
  
  dist_list <- lapply(names(bin_list), function(b) {
    df <- bin_list[[b]]
    # Weighted possible.contacts: mean of max observed contact distance per gene
    gene_max_dist  <- df[, .(max_dist=max(index)), by=gene_id]
    n_gene         <- floor(gene_max_dist$max_dist / resolution)
    total_possible <- mean(n_gene, na.rm=TRUE)
    
    dist <- df[, .(counts=sum(counts)), by=index]
    data.table::setorder(dist, index)
    n <- nrow(dist)
    dist[, possible.contacts := seq(total_possible*resolution, resolution,
                                    by=-resolution)[seq_len(n)]]
    dist[, contact.prob      := counts / possible.contacts]
    dist[, bin_val           := as.integer(b)]
  })
  
  all_dist <- data.table::rbindlist(dist_list)
  col_sum  <- sum(all_dist$contact.prob, na.rm=TRUE)
  if (col_sum == 0) return(NULL)
  
  all_dist[, contact.prob := contact.prob / col_sum]
  all_dist[index==0, index := 1L]
  all_dist[, `:=`(distance=log10(index), prob=log10(contact.prob))]
  all_dist[is.infinite(prob), prob := NA]
  all_dist[!is.na(prob) & distance > 0 & distance < max_distance]
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Extract contacts once per condition, compute slopes for each analysis
# ══════════════════════════════════════════════════════════════════════════════
# To clear all caches and start fresh, run:
# file.remove(list.files(RDS_DIR, pattern="\\.rds$", full.names=TRUE))

# ── 8a. Extract single-gene contacts once per condition ───────────────────────
message("\n══ Extracting single-gene contacts ══════════════════════════════════")
gene_contacts_list <- lapply(names(CONDITIONS), function(n) {
  cond         <- CONDITIONS[[n]]
  contacts_rds <- file.path(RDS_DIR, paste0("contacts_",n,"_genes.rds"))
  if (file.exists(contacts_rds)) {
    message(cond$display,": loading cached contacts...")
    readRDS(contacts_rds)
  } else {
    message(cond$display,": extracting contacts...")
    contacts <- extract_contacts_for_genes(cond$hic, ei_binned)
    contacts[, cond := cond$display]
    saveRDS(contacts, contacts_rds)
    message("  Saved: ",basename(contacts_rds))
    contacts
  }
})
names(gene_contacts_list) <- names(CONDITIONS)

# ── 8b. Extract cluster contacts once per condition ───────────────────────────
message("\n══ Extracting cluster contacts ══════════════════════════════════════")
cluster_contacts_list <- lapply(names(CONDITIONS), function(n) {
  cond         <- CONDITIONS[[n]]
  contacts_rds <- file.path(RDS_DIR, paste0("contacts_",n,"_clusters.rds"))
  if (file.exists(contacts_rds)) {
    message(cond$display,": loading cached cluster contacts...")
    readRDS(contacts_rds)
  } else {
    message(cond$display,": extracting cluster contacts...")
    contacts <- extract_contacts_for_genes(cond$hic, cluster_summary)
    contacts[, cond := cond$display]
    saveRDS(contacts, contacts_rds)
    message("  Saved: ",basename(contacts_rds))
    contacts
  }
})
names(cluster_contacts_list) <- names(CONDITIONS)

# ── Helper: slopes from pre-loaded contacts ───────────────────────────────────
get_slopes_from_contacts <- function(contacts, filter_col,
                                     cond_display, cache_path) {
  if (file.exists(cache_path)) {
    message("  Loading cached slopes: ",basename(cache_path))
    return(readRDS(cache_path))
  }
  message("  Computing slopes (filter_col=",filter_col,")...")
  slopes <- compute_slopes(contacts, filter_col=filter_col)
  if (!is.null(slopes)) slopes$cond <- cond_display
  saveRDS(slopes, cache_path)
  message("  Saved: ",basename(cache_path))
  slopes
}

# ── 8c. Main EI slopes ────────────────────────────────────────────────────────
message("\n══ Main EI slopes ═══════════════════════════════════════════════════")
slopes_ei <- bind_rows(lapply(names(CONDITIONS), function(n) {
  get_slopes_from_contacts(
    gene_contacts_list[[n]], "EI_bin",
    CONDITIONS[[n]]$display,
    file.path(RDS_DIR, paste0("slopes_",n,"_ei.rds")))
})) %>% dplyr::mutate(cond=factor(cond, levels=COND_LEVELS))

# ── 8d. Expression-stratified slopes ─────────────────────────────────────────
message("\n══ Expression-stratified slopes ═════════════════════════════════════")
slopes_expr <- bind_rows(lapply(names(CONDITIONS), function(n) {
  contacts <- data.table::copy(gene_contacts_list[[n]])
  bind_rows(lapply(c("Low","Medium","High"), function(eb) {
    sub <- contacts[as.character(expr_bin)==eb]
    if (nrow(sub)==0) return(NULL)
    s <- get_slopes_from_contacts(
      sub, "EI_bin",
      CONDITIONS[[n]]$display,
      file.path(RDS_DIR, paste0("slopes_",n,"_expr_",eb,".rds")))
    if (!is.null(s)) {
      s$expr_bin <- eb
      s$EI_label <- dplyr::case_when(
        s$bin_val <= 33 ~ "Low EI",
        s$bin_val <= 67 ~ "Medium EI",
        TRUE            ~ "High EI")
    }
    s
  }))
})) %>%
  dplyr::mutate(
    cond     = factor(cond,     levels=COND_LEVELS),
    expr_bin = factor(expr_bin, levels=c("Low","Medium","High")),
    EI_label = factor(EI_label, levels=c("Low EI","Medium EI","High EI"))
  )

# ── 8e. Pol2-stratified slopes ────────────────────────────────────────────────
message("\n══ Pol2-stratified slopes ═══════════════════════════════════════════")
slopes_pol2 <- bind_rows(lapply(names(CONDITIONS), function(n) {
  contacts <- data.table::copy(gene_contacts_list[[n]])
  bind_rows(lapply(c("Low","Medium","High"), function(pb) {
    sub <- contacts[as.character(pol2_bin)==pb]
    if (nrow(sub)==0) return(NULL)
    s <- get_slopes_from_contacts(
      sub, "EI_bin",
      CONDITIONS[[n]]$display,
      file.path(RDS_DIR, paste0("slopes_",n,"_pol2_",pb,".rds")))
    if (!is.null(s)) {
      s$pol2_bin <- pb
      s$EI_label <- dplyr::case_when(
        s$bin_val <= 33 ~ "Low EI",
        s$bin_val <= 67 ~ "Medium EI",
        TRUE            ~ "High EI")
    }
    s
  }))
})) %>%
  dplyr::mutate(
    cond     = factor(cond,     levels=COND_LEVELS),
    pol2_bin = factor(pol2_bin, levels=c("Low","Medium","High")),
    EI_label = factor(EI_label, levels=c("Low EI","Medium EI","High EI"))
  )

# ── 8f. Tandem cluster slopes ─────────────────────────────────────────────────
message("\n══ Tandem cluster slopes ════════════════════════════════════════════")
slopes_cluster <- bind_rows(lapply(names(CONDITIONS), function(n) {
  get_slopes_from_contacts(
    cluster_contacts_list[[n]], "EI_bin",
    CONDITIONS[[n]]$display,
    file.path(RDS_DIR, paste0("slopes_",n,"_cluster.rds")))
})) %>% dplyr::mutate(cond=factor(cond, levels=COND_LEVELS))

# ── 8g. Cluster expression-stratified slopes ──────────────────────────────────
message("\n══ Cluster expression-stratified slopes ═════════════════════════════")
slopes_cluster_expr <- bind_rows(lapply(names(CONDITIONS), function(n) {
  contacts <- data.table::copy(cluster_contacts_list[[n]])
  bind_rows(lapply(c("Low","Medium","High"), function(eb) {
    sub <- contacts[as.character(expr_bin)==eb]
    if (nrow(sub)==0) return(NULL)
    s <- get_slopes_from_contacts(
      sub, "EI_bin",
      CONDITIONS[[n]]$display,
      file.path(RDS_DIR, paste0("slopes_",n,"_cluster_expr_",eb,".rds")))
    if (!is.null(s)) {
      s$expr_bin <- eb
      s$EI_label <- dplyr::case_when(
        s$bin_val <= 33 ~ "Low EI",
        s$bin_val <= 67 ~ "Medium EI",
        TRUE            ~ "High EI")
    }
    s
  }))
})) %>%
  dplyr::mutate(
    cond     = factor(cond,     levels=COND_LEVELS),
    expr_bin = factor(expr_bin, levels=c("Low","Medium","High")),
    EI_label = factor(EI_label, levels=c("Low EI","Medium EI","High EI"))
  )

# ══════════════════════════════════════════════════════════════════════════════
# 8h. Diagnostics — genes per EI bin and P(S) curves
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Diagnostics ══════════════════════════════════════════════════════")

ei_bin_category <- function(df, bin_col="EI_bin", n_bins=EI_BINS) {
  df %>%
    dplyr::mutate(
      EI_category = dplyr::case_when(
        .data[[bin_col]] <= n_bins/3     ~ "Low EI",
        .data[[bin_col]] <= 2*(n_bins/3) ~ "Medium EI",
        TRUE                             ~ "High EI"
      ),
      EI_category = factor(EI_category,
                           levels=c("Low EI","Medium EI","High EI"))
    )
}

genes_per_bin <- ei_binned %>%
  dplyr::count(EI_bin) %>%
  ei_bin_category()
message("Genes per EI bin — single genes:")
message("  Min=",min(genes_per_bin$n),
        "  Median=",median(genes_per_bin$n),
        "  Max=",max(genes_per_bin$n),
        "  Bins<10: ",sum(genes_per_bin$n<10))

p_genes_per_bin <- ggplot(genes_per_bin,
                          aes(x=EI_bin, y=n, fill=EI_category))+
  geom_col(width=0.9)+
  geom_hline(yintercept=10,linetype="dashed",color="gray30",linewidth=0.6)+
  scale_fill_manual(values=c("Low EI"="#56B4E9",
                             "Medium EI"="#E69F00",
                             "High EI"="#D55E00"),name=NULL)+
  scale_x_continuous(breaks=c(1,EI_BINS/3,2*EI_BINS/3,EI_BINS),
                     labels=c("1","Low|Med","Med|High",as.character(EI_BINS)))+
  scale_y_continuous(expand=expansion(mult=c(0,0.05)))+
  labs(x="EI bin",y="Number of genes",
       title="Genes per EI bin — single genes ≥ 100kb",
       subtitle=paste0(EI_BINS," bins  |  dashed = 10 gene minimum"))+
  theme_classic(base_size=12)+
  theme(plot.title=element_text(face="bold"),legend.position="bottom")
print(p_genes_per_bin)

cluster_per_bin <- cluster_summary %>%
  dplyr::count(EI_bin) %>%
  ei_bin_category()
message("Clusters per EI bin:")
message("  Min=",min(cluster_per_bin$n),
        "  Median=",median(cluster_per_bin$n),
        "  Max=",max(cluster_per_bin$n),
        "  Bins<10: ",sum(cluster_per_bin$n<10))

p_cluster_per_bin <- ggplot(cluster_per_bin,
                            aes(x=EI_bin, y=n, fill=EI_category))+
  geom_col(width=0.9)+
  geom_hline(yintercept=10,linetype="dashed",color="gray30",linewidth=0.6)+
  scale_fill_manual(values=c("Low EI"="#56B4E9",
                             "Medium EI"="#E69F00",
                             "High EI"="#D55E00"),name=NULL)+
  scale_x_continuous(breaks=c(1,EI_BINS/3,2*EI_BINS/3,EI_BINS),
                     labels=c("1","Low|Med","Med|High",as.character(EI_BINS)))+
  scale_y_continuous(expand=expansion(mult=c(0,0.05)))+
  labs(x="EI bin",y="Number of clusters",
       title="Clusters per EI bin — tandem clusters",
       subtitle=paste0(EI_BINS," bins  |  dashed = 10 cluster minimum"))+
  theme_classic(base_size=12)+
  theme(plot.title=element_text(face="bold"),legend.position="bottom")
print(p_cluster_per_bin)

message("Computing P(S) curves (WT)...")
ps_curves <- compute_ps_curves(gene_contacts_list$wt,
                               filter_col="EI_bin") %>%
  as.data.frame() %>%
  ei_bin_category(bin_col="bin_val")

ps_summary <- ps_curves %>%
  dplyr::group_by(bin_val, EI_category,
                  distance=round(distance,2)) %>%
  dplyr::summarise(prob=mean(prob,na.rm=TRUE),.groups="drop")

p_ps <- ggplot(ps_summary,
               aes(x=distance,y=prob,
                   group=bin_val,color=EI_category))+
  geom_line(alpha=0.6,linewidth=0.5)+
  scale_color_manual(values=c("Low EI"="#56B4E9",
                              "Medium EI"="#E69F00",
                              "High EI"="#D55E00"),name=NULL)+
  scale_x_continuous(breaks=c(4,4.5,5,5.5,6),
                     labels=c("10kb","30kb","100kb","300kb","1Mb"))+
  facet_wrap(~EI_category,nrow=1)+
  labs(x="Genomic distance (log scale)",y="log10 P(s)",
       title=paste0("P(s) curves by EI bin  |  Wildtype  |  ",DS_GROUP),
       subtitle=paste0("Each line = one EI bin  |  ",EI_BINS,
                       " bins  |  colored by Low/Med/High EI"))+
  theme_classic(base_size=12)+
  theme(plot.title=element_text(face="bold"),
        strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        legend.position="none")
print(p_ps)

# Extended P(s) — full range to 10^8
ps_curves_wide <- compute_ps_curves(gene_contacts_list$wt,
                                    filter_col="EI_bin",
                                    max_distance=8) %>%
  as.data.frame() %>%
  ei_bin_category(bin_col="bin_val")

ps_summary_wide <- ps_curves_wide %>%
  dplyr::group_by(bin_val, EI_category,
                  distance=round(distance,2)) %>%
  dplyr::summarise(prob=mean(prob,na.rm=TRUE),.groups="drop")

p_ps_wide <- ggplot(ps_summary_wide,
                    aes(x=distance,y=prob,
                        group=bin_val,color=EI_category))+
  geom_line(alpha=0.6,linewidth=0.5)+
  scale_color_manual(values=c("Low EI"="#56B4E9",
                              "Medium EI"="#E69F00",
                              "High EI"="#D55E00"),name=NULL)+
  scale_x_continuous(breaks=c(4,5,6,7,8),
                     labels=c("10kb","100kb","1Mb","10Mb","100Mb"))+
  facet_wrap(~EI_category,nrow=1)+
  labs(x="Genomic distance (log scale)",y="log10 P(s)",
       title=paste0("P(s) curves — full range  |  Wildtype  |  ",DS_GROUP),
       subtitle=paste0("Each line = one EI bin  |  ",EI_BINS,
                       " bins  |  extended to 10^8 bp"))+
  theme_classic(base_size=12)+
  theme(plot.title=element_text(face="bold"),
        strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        legend.position="none")
print(p_ps_wide)


# ══════════════════════════════════════════════════════════════════════════════
# 9. Plot helpers
# ══════════════════════════════════════════════════════════════════════════════
ei_range_labels <- function(slopes_df, ei_df_src, n_bins=EI_BINS) {
  # Map EI bin index to EI range for axis labels
  breaks <- quantile(ei_df_src$EI, probs=seq(0,1,length.out=n_bins+1), na.rm=TRUE)
  yax    <- format(breaks[seq(1,n_bins+1,length.out=5)], scientific=FALSE, digits=3)
  list(breaks=seq(1,n_bins,length.out=5), labels=as.numeric(yax))
}

theme_scaling <- function() {
  theme_classic(base_size=13) +
    theme(
      plot.title     = element_text(face="bold",size=14,hjust=0.5),
      plot.subtitle  = element_text(color="grey40",size=11,hjust=0.5),
      legend.position= "bottom",
      text           = element_text(family="Arial"),
      panel.border   = element_rect(color="black",fill=NA,linewidth=0.5)
    )
}


# ══════════════════════════════════════════════════════════════════════════════
# Helper: make_halfeye — vertical raincloud with median annotation
# ══════════════════════════════════════════════════════════════════════════════
make_halfeye <- function(slopes_df, x_col, fill_col, color_col,
                         x_labels=NULL, title="", subtitle="",
                         facet_vars=NULL) {
  group_cols <- unique(c(x_col, facet_vars))
  med_df <- slopes_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(med=median(slope, na.rm=TRUE), .groups="drop")
  
  p <- ggplot(slopes_df,
              aes(x=factor(.data[[x_col]]), y=slope,
                  fill=.data[[fill_col]], color=.data[[color_col]]))+
    ggdist::stat_dots(
      side="left", justification=1.12,
      binwidth=0.008, dotsize=0.8, alpha=0.4)+
    geom_boxplot(
      width=0.15, outlier.shape=NA, alpha=0.8, linewidth=0.5)+
    ggdist::stat_halfeye(
      adjust=0.6, width=0.5, justification=-0.25,
      .width=0, point_colour=NA, slab_alpha=0.75)+
    geom_text(
      data=med_df,
      aes(x=factor(.data[[x_col]]), y=med,
          label=round(med, 3)),
      inherit.aes=FALSE,
      hjust=-0.2, vjust=0.4,
      size=2.8, fontface="bold", color="gray20")+
    scale_fill_manual(values=COND_COLORS, guide="none")+
    scale_color_manual(values=COND_COLORS, guide="none")+
    labs(title=title, subtitle=subtitle, y="|S|", x="")+
    { if (!is.null(x_labels))
      scale_x_discrete(labels=x_labels)
      else
        scale_x_discrete() }+
    theme_classic(base_size=13)+
    theme(
      plot.title         = element_text(face="bold", size=14, hjust=0.5),
      plot.subtitle      = element_text(color="grey40", size=10, hjust=0.5),
      axis.text          = element_text(size=11, color="black"),
      axis.title.y       = element_text(size=13, face="bold"),
      axis.title.x       = element_blank(),
      axis.ticks         = element_line(color="black", linewidth=0.4),
      axis.line          = element_line(color="black", linewidth=0.4),
      panel.grid.major.y = element_line(color="grey92", linewidth=0.3),
      panel.grid.major.x = element_blank(),
      strip.text         = element_text(face="bold", size=11),
      strip.background   = element_blank(),
      legend.position    = "none",
      plot.margin        = margin(8, 16, 8, 8)
    )
  p
}

# Helper: delta line panel (WT - pol2 per bin_val)
make_delta_panel <- function(slopes_df, bin_col="bin_val",
                             x_medians=NULL,
                             title="ΔS (WT − POLR2A deg)") {
  wt_label   <- COND_LEVELS[1]
  pol2_label <- COND_LEVELS[2]
  
  delta_df <- slopes_df %>%
    dplyr::select(bin_val=dplyr::all_of(bin_col), slope, cond) %>%
    tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
    dplyr::mutate(delta = .data[[wt_label]] - .data[[pol2_label]]) %>%
    dplyr::filter(!is.na(delta))
  
  # 5 evenly spaced x-axis labels from medians
  x_breaks <- x_labels <- NULL
  if (!is.null(x_medians)) {
    idx       <- round(seq(1, nrow(x_medians), length.out=5))
    x_breaks  <- x_medians$EI_bin[idx]
    x_labels  <- sprintf("%.2g", x_medians$median_EI[idx])
  }
  
  p <- ggplot(delta_df, aes(x=bin_val, y=delta))+
    geom_hline(yintercept=0, linetype="dashed",
               color="gray40", linewidth=0.5)+
    geom_line(color="#D55E00", linewidth=0.8)+
    geom_area(alpha=0.15, fill="#D55E00")+
    labs(x="E/I Ratio", y="ΔS", title=title)+
    theme_classic(base_size=11)+
    theme(
      plot.title         = element_text(face="bold", size=11, hjust=0.5),
      axis.text          = element_text(size=9, color="black"),
      axis.title         = element_text(size=10, color="black"),
      panel.grid.major.y = element_line(color="grey92", linewidth=0.3),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(4,16,8,8)
    )
  if (!is.null(x_breaks))
    p <- p + scale_x_continuous(breaks=x_breaks, labels=x_labels)
  p
}

# Helper: assemble halfeye + delta with cowplot
assemble_with_delta <- function(p_half, p_delta, title,
                                rel_heights=c(2,1)) {
  title_grob <- cowplot::ggdraw()+
    cowplot::draw_label(title, fontface="bold", size=13, x=0.01, hjust=0)
  cowplot::plot_grid(
    title_grob,
    cowplot::plot_grid(p_half, p_delta, ncol=1, rel_heights=rel_heights),
    ncol=1, rel_heights=c(0.05,1))
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. Plots — Main EI analysis
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Plotting: Main EI analysis ═══════════════════════════════════════")

ei_bin_medians <- ei_binned %>%
  dplyr::group_by(EI_bin) %>%
  dplyr::summarise(median_EI=median(EI, na.rm=TRUE), .groups="drop")

cluster_bin_medians <- cluster_summary %>%
  dplyr::group_by(EI_bin) %>%
  dplyr::summarise(median_EI=median(EI_cluster, na.rm=TRUE), .groups="drop")

slopes_ei_plot <- slopes_ei %>%
  dplyr::left_join(ei_bin_medians, by=c("bin_val"="EI_bin"))

p_ei_scatter <- ggplot(slopes_ei_plot,
                       aes(x=slope, y=median_EI, color=cond))+
  geom_point(size=4, alpha=0.7)+
  scale_color_manual(values=COND_COLORS, name=NULL)+
  scale_y_log10(
    breaks=c(0.01,0.1,1,10,100,1000,10000),
    labels=c("0.01","0.1","1","10","100","1,000","10,000")
  )+
  labs(title=paste0("Contact Scaling by E/I Ratio  |  ",DS_GROUP),
       subtitle="Genes > 100kb  |  100 EI bins",
       x="|S|", y="E/I Ratio")+
  theme_scaling()
print(p_ei_scatter)

# Add EI_label to slopes_ei before filtering
slopes_ei_4 <- slopes_ei %>%
  dplyr::mutate(
    bin_val3 = as.integer(ntile(bin_val, 3)),
    EI_label = dplyr::case_when(
      bin_val3==1 ~ "Low EI",
      bin_val3==2 ~ "Medium EI",
      TRUE        ~ "High EI"),
    EI_label = factor(EI_label, levels=c("Low EI","Medium EI","High EI"))
  ) %>%
  dplyr::filter(EI_label != "Medium EI") %>%
  dplyr::mutate(EI_label=factor(EI_label, levels=c("Low EI","High EI")))

stat_ei <- slopes_ei_4 %>%
  group_by(EI_label) %>%
  dplyr::filter(all(table(cond) >= 2)) %>%
  rstatix::t_test(slope~cond) %>%
  rstatix::add_significance() %>%
  dplyr::mutate(y.position=max(slopes_ei_4$slope,na.rm=TRUE)*1.05)

p_ei_half <- make_halfeye(slopes_ei_4,"cond","cond","cond",
                          x_labels=c("WT","POLR2A deg"),
                          facet_vars="EI_label",
                          title=paste0("Contact Scaling  |  ",DS_GROUP),
                          subtitle="Genes > 100kb  |  Low & High EI")+
  facet_wrap(~EI_label, nrow=1)+
  stat_pvalue_manual(stat_ei, inherit.aes=FALSE, label="p.signif")

p_ei_delta <- make_delta_panel(slopes_ei, x_medians=ei_bin_medians,
                               title="ΔS (WT − POLR2A deg)  |  Main EI")

print(assemble_with_delta(p_ei_half, p_ei_delta,
                          title=paste0("Contact Scaling by E/I  |  ",DS_GROUP)))

# ══════════════════════════════════════════════════════════════════════════════
# 11. Plots — Expression-stratified
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Plotting: Expression-stratified ═════════════════════════════════")

# EI_label already on slopes_expr from computation
slopes_expr_filt <- slopes_expr %>%
  dplyr::filter(EI_label != "Medium EI") %>%
  dplyr::mutate(EI_label=factor(EI_label, levels=c("Low EI","High EI")))

for (eb in c("High","Low")) {
  sub <- slopes_expr_filt %>% dplyr::filter(as.character(expr_bin)==eb)
  if (nrow(sub)==0) next
  stat_e <- sub %>%
    group_by(EI_label) %>%
    dplyr::filter(all(table(cond) >= 2)) %>%
    rstatix::t_test(slope~cond) %>%
    rstatix::add_significance() %>%
    dplyr::mutate(y.position=max(sub$slope,na.rm=TRUE)*1.05)
  
  p_half <- make_halfeye(sub,"cond","cond","cond",
                         x_labels=c("WT","POLR2A deg"),
                         facet_vars="EI_label",
                         title=paste0(eb," Expression  |  ",DS_GROUP),
                         subtitle="100 EI bins  |  Low & High EI")+
    facet_wrap(~EI_label, nrow=1)+
    stat_pvalue_manual(stat_e, inherit.aes=FALSE, label="p.signif")
  
  p_delta <- make_delta_panel(
    slopes_expr %>% dplyr::filter(as.character(expr_bin)==eb),
    x_medians=ei_bin_medians,
    title=paste0("ΔS (WT − POLR2A deg)  |  ",eb," expression"))
  
  print(assemble_with_delta(p_half, p_delta,
                            title=paste0("Scaling by Expression — ",eb,"  |  ",DS_GROUP)))
}

p_expr_scatter <- ggplot(
  slopes_expr %>% dplyr::left_join(ei_bin_medians, by=c("bin_val"="EI_bin")),
  aes(x=slope, y=median_EI, color=cond))+
  geom_point(size=3, alpha=0.7)+
  scale_color_manual(values=COND_COLORS, name=NULL)+
  scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000),
                labels=c("0.01","0.1","1","10","100","1,000","10,000"))+
  labs(title=paste0("Scaling vs EI by Expression  |  ",DS_GROUP),
       x="|S|", y="E/I Ratio")+
  facet_wrap(~expr_bin, nrow=1)+
  theme_scaling()
print(p_expr_scatter)

# ══════════════════════════════════════════════════════════════════════════════
# 12. Plots — Pol2 occupancy-stratified
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Plotting: Pol2 occupancy-stratified ══════════════════════════════")

# EI_label already on slopes_pol2 from computation
slopes_pol2_filt <- slopes_pol2 %>%
  dplyr::filter(EI_label != "Medium EI") %>%
  dplyr::mutate(EI_label=factor(EI_label, levels=c("Low EI","High EI")))

for (pb in c("High","Low")) {
  sub <- slopes_pol2_filt %>% dplyr::filter(as.character(pol2_bin)==pb)
  if (nrow(sub)==0) next
  stat_p <- sub %>%
    group_by(EI_label) %>%
    dplyr::filter(all(table(cond) >= 2)) %>%
    rstatix::t_test(slope~cond) %>%
    rstatix::add_significance() %>%
    dplyr::mutate(y.position=max(sub$slope,na.rm=TRUE)*1.05)
  
  p_half <- make_halfeye(sub,"cond","cond","cond",
                         x_labels=c("WT","POLR2A deg"),
                         facet_vars="EI_label",
                         title=paste0(pb," Pol2 Occupancy  |  ",DS_GROUP),
                         subtitle="100 EI bins  |  Low & High EI")+
    facet_wrap(~EI_label, nrow=1)+
    stat_pvalue_manual(stat_p, inherit.aes=FALSE, label="p.signif")
  
  p_delta <- make_delta_panel(
    slopes_pol2 %>% dplyr::filter(as.character(pol2_bin)==pb),
    x_medians=ei_bin_medians,
    title=paste0("ΔS (WT − POLR2A deg)  |  ",pb," Pol2"))
  
  print(assemble_with_delta(p_half, p_delta,
                            title=paste0("Scaling by Pol2 — ",pb,"  |  ",DS_GROUP)))
}

p_pol2_scatter <- ggplot(
  slopes_pol2 %>% dplyr::left_join(ei_bin_medians, by=c("bin_val"="EI_bin")),
  aes(x=slope, y=median_EI, color=cond))+
  geom_point(size=3, alpha=0.7)+
  scale_color_manual(values=COND_COLORS, name=NULL)+
  scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000),
                labels=c("0.01","0.1","1","10","100","1,000","10,000"))+
  labs(title=paste0("Scaling vs EI by Pol2 Occupancy  |  ",DS_GROUP),
       x="|S|", y="E/I Ratio")+
  facet_wrap(~pol2_bin, nrow=1)+
  theme_scaling()
print(p_pol2_scatter)

# ══════════════════════════════════════════════════════════════════════════════
# 13. Plots — Tandem cluster EI
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Plotting: Tandem cluster EI ══════════════════════════════════════")

slopes_cluster_plot <- slopes_cluster %>%
  dplyr::left_join(cluster_bin_medians, by=c("bin_val"="EI_bin"))

p_cluster_scatter <- ggplot(slopes_cluster_plot,
                            aes(x=slope, y=median_EI, color=cond))+
  geom_point(size=4, alpha=0.7)+
  scale_color_manual(values=COND_COLORS, name=NULL)+
  scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000),
                labels=c("0.01","0.1","1","10","100","1,000","10,000"))+
  labs(title=paste0("Scaling by Tandem Cluster E/I  |  ",DS_GROUP),
       subtitle=paste0("Genes < 100kb  |  max gap=",MAX_GAP/1000,"kb  |  100 EI bins"),
       x="|S|", y="Cluster E/I Ratio")+
  theme_scaling()
print(p_cluster_scatter)

# Add EI_label to slopes_cluster before filtering
slopes_cluster_4 <- slopes_cluster %>%
  dplyr::mutate(
    bin_val3 = as.integer(ntile(bin_val, 3)),
    EI_label = dplyr::case_when(
      bin_val3==1 ~ "Low E/I",
      bin_val3==2 ~ "Medium E/I",
      TRUE        ~ "High E/I"),
    EI_label = factor(EI_label, levels=c("Low E/I","Medium E/I","High E/I"))
  ) %>%
  dplyr::filter(EI_label != "Medium E/I") %>%
  dplyr::mutate(EI_label=factor(EI_label, levels=c("Low E/I","High E/I")))

stat_cluster <- slopes_cluster_4 %>%
  group_by(EI_label) %>%
  dplyr::filter(all(table(cond) >= 2)) %>%
  rstatix::t_test(slope~cond) %>%
  rstatix::add_significance() %>%
  dplyr::mutate(y.position=max(slopes_cluster_4$slope,na.rm=TRUE)*1.05)

p_cluster_half <- make_halfeye(slopes_cluster_4,"cond","cond","cond",
                               x_labels=c("WT","POLR2A deg"),
                               facet_vars="EI_label",
                               title=paste0("Cluster Scaling  |  ",DS_GROUP),
                               subtitle="E = cluster exon bp  |  I = cluster intergenic bp")+
  facet_wrap(~EI_label, nrow=1)+
  stat_pvalue_manual(stat_cluster, inherit.aes=FALSE, label="p.signif")

p_cluster_delta <- make_delta_panel(slopes_cluster,
                                    x_medians=cluster_bin_medians,
                                    title="ΔS (WT − POLR2A deg)  |  Tandem clusters")

print(assemble_with_delta(p_cluster_half, p_cluster_delta,
                          title=paste0("Cluster Scaling by E/I  |  ",DS_GROUP)))

# ══════════════════════════════════════════════════════════════════════════════
# 14. Plots — Tandem cluster expression-stratified
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Plotting: Tandem cluster expression-stratified ═══════════════════")

# EI_label already on slopes_cluster_expr from computation
slopes_cluster_expr_filt <- slopes_cluster_expr %>%
  dplyr::filter(EI_label != "Medium EI") %>%
  dplyr::mutate(EI_label=factor(EI_label, levels=c("Low EI","High EI")))

for (eb in c("High","Low")) {
  sub <- slopes_cluster_expr_filt %>% dplyr::filter(as.character(expr_bin)==eb)
  if (nrow(sub)==0) next
  stat_ce <- sub %>%
    group_by(EI_label) %>%
    dplyr::filter(all(table(cond) >= 2)) %>%
    rstatix::t_test(slope~cond) %>%
    rstatix::add_significance() %>%
    dplyr::mutate(y.position=max(sub$slope,na.rm=TRUE)*1.05)
  
  p_half <- make_halfeye(sub,"cond","cond","cond",
                         x_labels=c("WT","POLR2A deg"),
                         facet_vars="EI_label",
                         title=paste0("Cluster — ",eb," Expression  |  ",DS_GROUP),
                         subtitle="Mean cluster TPM  |  Low & High EI")+
    facet_wrap(~EI_label, nrow=1)+
    stat_pvalue_manual(stat_ce, inherit.aes=FALSE, label="p.signif")
  
  p_delta <- make_delta_panel(
    slopes_cluster_expr %>% dplyr::filter(as.character(expr_bin)==eb),
    x_medians=cluster_bin_medians,
    title=paste0("ΔS (WT − POLR2A deg)  |  Cluster ",eb," expression"))
  
  print(assemble_with_delta(p_half, p_delta,
                            title=paste0("Cluster Scaling — ",eb," Expression  |  ",DS_GROUP)))
}

p_cluster_expr_scatter <- ggplot(
  slopes_cluster_expr %>%
    dplyr::left_join(cluster_bin_medians, by=c("bin_val"="EI_bin")),
  aes(x=slope, y=median_EI, color=cond))+
  geom_point(size=3, alpha=0.7)+
  scale_color_manual(values=COND_COLORS, name=NULL)+
  scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000),
                labels=c("0.01","0.1","1","10","100","1,000","10,000"))+
  labs(title=paste0("Cluster Scaling vs EI by Expression  |  ",DS_GROUP),
       x="|S|", y="Cluster E/I Ratio")+
  facet_wrap(~expr_bin, nrow=1)+
  theme_scaling()
print(p_cluster_expr_scatter)

# ══════════════════════════════════════════════════════════════════════════════
# 15. Expression delta analysis
#     - 100 expression bins computed from TPM
#     - Joined onto contacts in memory (no re-extraction)
#     - S computed per expression bin
#     - Delta S (WT - pol2) plotted vs mean TPM per bin
#     - Also: halfeye for Low/High EI × Low/High expression with shared y-axis
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Section 15: Expression delta analysis ════════════════════════════")

# ── 15a. Add 100-bin expression to ei_binned ──────────────────────────────────
ei_binned <- ei_binned %>%
  dplyr::mutate(expr_bin_100 = as.integer(ntile(TPM, 100L)))

# ── 15b. Join expr_bin_100 onto loaded contacts in memory ─────────────────────
message("Joining expr_bin_100 onto gene contacts...")
gene_contacts_list <- lapply(gene_contacts_list, function(ct) {
  data.table::as.data.table(ct) %>%
    dplyr::left_join(
      ei_binned %>% dplyr::select(gene_id, expr_bin_100),
      by="gene_id"
    ) %>%
    data.table::as.data.table()
})

# ── 15c. Compute slopes by expr_bin_100 ───────────────────────────────────────
message("Computing slopes by expression bin (100 bins)...")
slopes_expr_100 <- bind_rows(lapply(names(CONDITIONS), function(n) {
  get_slopes_from_contacts(
    gene_contacts_list[[n]], "expr_bin_100",
    CONDITIONS[[n]]$display,
    file.path(RDS_DIR, paste0("slopes_",n,"_expr100.rds")))
})) %>% dplyr::mutate(cond=factor(cond, levels=COND_LEVELS))

# ── 15d. Delta S vs mean TPM ──────────────────────────────────────────────────
expr_bin_tpm <- ei_binned %>%
  dplyr::group_by(expr_bin_100) %>%
  dplyr::summarise(mean_TPM=mean(TPM, na.rm=TRUE), .groups="drop")

delta_expr <- slopes_expr_100 %>%
  dplyr::select(bin_val, slope, cond) %>%
  tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
  dplyr::mutate(delta = .data[[COND_LEVELS[1]]] - .data[[COND_LEVELS[2]]]) %>%
  dplyr::filter(!is.na(delta)) %>%
  dplyr::left_join(expr_bin_tpm, by=c("bin_val"="expr_bin_100"))

p_delta_expr <- ggplot(delta_expr, aes(x=mean_TPM, y=delta))+
  geom_hline(yintercept=0, linetype="dashed",
             color="gray40", linewidth=0.5)+
  geom_line(color="#D55E00", linewidth=0.8)+
  geom_area(alpha=0.15, fill="#D55E00")+
  scale_x_log10(labels=scales::comma)+
  labs(x="Mean TPM (log scale)", y="ΔS",
       title=paste0("ΔS (WT − POLR2A deg) vs Expression  |  ",DS_GROUP),
       subtitle="100 expression bins  |  genes ≥ 100kb")+
  theme_classic(base_size=13)+
  theme(plot.title=element_text(face="bold",size=14,hjust=0.5),
        plot.subtitle=element_text(color="grey40",size=10,hjust=0.5),
        axis.text=element_text(size=11,color="black"),
        axis.title=element_text(size=12,color="black"),
        panel.grid.major.y=element_line(color="grey92",linewidth=0.3),
        panel.grid.major.x=element_blank())
print(p_delta_expr)

# ── 15e. Scatter: S vs TPM per bin ───────────────────────────────────────────
slopes_expr_100_plot <- slopes_expr_100 %>%
  dplyr::left_join(expr_bin_tpm, by=c("bin_val"="expr_bin_100"))

p_expr100_scatter <- ggplot(slopes_expr_100_plot,
                            aes(x=slope, y=mean_TPM, color=cond))+
  geom_point(size=3, alpha=0.7)+
  scale_color_manual(values=COND_COLORS, name=NULL)+
  scale_y_log10(labels=scales::comma)+
  labs(title=paste0("Contact Scaling vs Expression  |  ",DS_GROUP),
       subtitle="100 expression bins  |  genes ≥ 100kb",
       x="|S|", y="Mean TPM")+
  theme_scaling()
print(p_expr100_scatter)

# ── 15f. Halfeye: Low vs High expression, shared y-axis shown once ────────────
message("Halfeye plots: Low vs High expression...")

slopes_expr_lh <- slopes_expr %>%
  dplyr::filter(as.character(expr_bin) != "Medium") %>%
  dplyr::mutate(expr_bin=factor(expr_bin, levels=c("Low","High")))

shared_ylim <- range(slopes_expr_lh$slope, na.rm=TRUE)
shared_ylim <- c(floor(shared_ylim[1]*10)/10 - 0.05,
                 ceiling(shared_ylim[2]*10)/10 + 0.05)

# Low expression halfeye
stat_low <- slopes_expr_lh %>%
  dplyr::filter(as.character(expr_bin)=="Low") %>%
  dplyr::filter(all(table(cond) >= 2)) %>%
  rstatix::t_test(slope~cond) %>%
  rstatix::add_significance() %>%
  dplyr::mutate(y.position=shared_ylim[2]*0.95)

p_low <- make_halfeye(slopes_expr_lh %>% dplyr::filter(as.character(expr_bin)=="Low"),
                      "cond","cond","cond",
                      x_labels=c("WT","POLR2A deg"),
                      title="Low Expression")+
  coord_cartesian(ylim=shared_ylim)+
  stat_pvalue_manual(stat_low, inherit.aes=FALSE, label="p.signif")

# High expression halfeye — remove y-axis so it's only shown once on the left
stat_high <- slopes_expr_lh %>%
  dplyr::filter(as.character(expr_bin)=="High") %>%
  dplyr::filter(all(table(cond) >= 2)) %>%
  rstatix::t_test(slope~cond) %>%
  rstatix::add_significance() %>%
  dplyr::mutate(y.position=shared_ylim[2]*0.95)

p_high <- make_halfeye(slopes_expr_lh %>% dplyr::filter(as.character(expr_bin)=="High"),
                       "cond","cond","cond",
                       x_labels=c("WT","POLR2A deg"),
                       title="High Expression")+
  coord_cartesian(ylim=shared_ylim)+
  stat_pvalue_manual(stat_high, inherit.aes=FALSE, label="p.signif")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())

# Delta line plot
p_delta <- ggplot(delta_expr, aes(x=mean_TPM, y=delta))+
  geom_hline(yintercept=0, linetype="dashed",
             color="gray40", linewidth=0.5)+
  geom_line(color="#D55E00", linewidth=0.8)+
  geom_area(alpha=0.15, fill="#D55E00")+
  scale_x_log10(labels=scales::comma)+
  labs(x="Mean TPM (log scale)", y="ΔS",
       title="ΔS (WT − POLR2A deg) vs Expression")+
  theme_classic(base_size=11)+
  theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
        axis.text=element_text(size=9,color="black"),
        axis.title=element_text(size=10,color="black"),
        panel.grid.major.y=element_line(color="grey92",linewidth=0.3),
        panel.grid.major.x=element_blank(),
        plot.margin=margin(4,16,8,8))

top_row    <- cowplot::plot_grid(p_low, p_high, nrow=1, rel_widths=c(1,1))
title_grob <- cowplot::ggdraw()+
  cowplot::draw_label(paste0("Scaling by Expression  |  ",DS_GROUP),
                      fontface="bold", size=13, x=0.01, hjust=0)
final <- cowplot::plot_grid(
  title_grob,
  cowplot::plot_grid(top_row, p_delta, ncol=1, rel_heights=c(2,1)),
  ncol=1, rel_heights=c(0.05,1))
print(final)


# ══════════════════════════════════════════════════════════════════════════════
# 15b. Expression delta analysis — Tandem clusters
#     - 100 expression bins computed from mean cluster TPM
#     - Joined onto cluster contacts in memory
#     - S computed per expression bin
#     - Delta S (WT - pol2) plotted vs mean TPM per bin
#     - Halfeye for Low/High expression with shared y-axis
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Section 15b: Cluster expression delta analysis ═══════════════════")

# ── 15b-i. Add 100-bin expression to cluster_summary ─────────────────────────
cluster_summary <- cluster_summary %>%
  dplyr::mutate(
    expr_bin_100 = dplyr::if_else(
      !is.na(mean_TPM) & mean_TPM > 0,
      as.integer(ntile(mean_TPM, 100L)),
      NA_integer_)
  )

message("expr_bin_100 range: ",
        paste(range(cluster_summary$expr_bin_100, na.rm=TRUE), collapse=" - "))
message("NAs: ", sum(is.na(cluster_summary$expr_bin_100)))

# ── 15b-ii. Join expr_bin_100 onto cluster contacts in memory ─────────────────
message("Joining expr_bin_100 onto cluster contacts...")
cluster_contacts_list <- lapply(cluster_contacts_list, function(ct) {
  data.table::as.data.table(ct) %>%
    dplyr::left_join(
      cluster_summary %>% dplyr::select(gene_id, expr_bin_100),
      by="gene_id"
    ) %>%
    data.table::as.data.table()
})

# ── 15b-iii. Compute slopes by expr_bin_100 ───────────────────────────────────
message("Computing cluster slopes by expression bin (100 bins)...")
slopes_cluster_expr_100 <- bind_rows(lapply(names(CONDITIONS), function(n) {
  get_slopes_from_contacts(
    cluster_contacts_list[[n]], "expr_bin_100",
    CONDITIONS[[n]]$display,
    file.path(RDS_DIR, paste0("slopes_",n,"_cluster_expr100.rds")))
})) %>% dplyr::mutate(cond=factor(cond, levels=COND_LEVELS))

# ── 15b-iv. Delta S vs mean cluster TPM ───────────────────────────────────────
cluster_expr_bin_tpm <- cluster_summary %>%
  dplyr::filter(!is.na(expr_bin_100)) %>%
  dplyr::group_by(expr_bin_100) %>%
  dplyr::summarise(mean_TPM=mean(mean_TPM, na.rm=TRUE), .groups="drop")

delta_cluster_expr <- slopes_cluster_expr_100 %>%
  dplyr::select(bin_val, slope, cond) %>%
  tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
  dplyr::mutate(delta = .data[[COND_LEVELS[1]]] - .data[[COND_LEVELS[2]]]) %>%
  dplyr::filter(!is.na(delta)) %>%
  dplyr::left_join(cluster_expr_bin_tpm, by=c("bin_val"="expr_bin_100"))

p_delta_cluster_expr <- ggplot(delta_cluster_expr, aes(x=mean_TPM, y=delta))+
  geom_hline(yintercept=0, linetype="dashed",
             color="gray40", linewidth=0.5)+
  geom_line(color="#D55E00", linewidth=0.8)+
  geom_area(alpha=0.15, fill="#D55E00")+
  scale_x_log10(labels=scales::comma)+
  labs(x="Mean cluster TPM (log scale)", y="ΔS",
       title=paste0("ΔS (WT − POLR2A deg) vs Expression  |  ",DS_GROUP),
       subtitle="100 expression bins  |  tandem clusters")+
  theme_classic(base_size=13)+
  theme(plot.title=element_text(face="bold",size=14,hjust=0.5),
        plot.subtitle=element_text(color="grey40",size=10,hjust=0.5),
        axis.text=element_text(size=11,color="black"),
        axis.title=element_text(size=12,color="black"),
        panel.grid.major.y=element_line(color="grey92",linewidth=0.3),
        panel.grid.major.x=element_blank())
print(p_delta_cluster_expr)

# ── 15b-v. Scatter: S vs mean cluster TPM per bin ────────────────────────────
slopes_cluster_expr_100_plot <- slopes_cluster_expr_100 %>%
  dplyr::left_join(cluster_expr_bin_tpm, by=c("bin_val"="expr_bin_100"))

p_cluster_expr100_scatter <- ggplot(slopes_cluster_expr_100_plot,
                                    aes(x=slope, y=mean_TPM, color=cond))+
  geom_point(size=3, alpha=0.7)+
  scale_color_manual(values=COND_COLORS, name=NULL)+
  scale_y_log10(labels=scales::comma)+
  labs(title=paste0("Cluster Scaling vs Expression  |  ",DS_GROUP),
       subtitle="100 expression bins  |  tandem clusters",
       x="|S|", y="Mean cluster TPM")+
  theme_scaling()
print(p_cluster_expr100_scatter)

# ── 15b-vi. Halfeye: Low vs High expression, shared y-axis ───────────────────
message("Cluster halfeye: Low vs High expression...")

slopes_cluster_expr_lh <- slopes_cluster_expr %>%
  dplyr::filter(as.character(expr_bin) != "Medium") %>%
  dplyr::mutate(expr_bin=factor(expr_bin, levels=c("Low","High")))

shared_ylim_cl <- range(slopes_cluster_expr_lh$slope, na.rm=TRUE)
shared_ylim_cl <- c(floor(shared_ylim_cl[1]*10)/10 - 0.05,
                    ceiling(shared_ylim_cl[2]*10)/10 + 0.05)

stat_cl_low <- slopes_cluster_expr_lh %>%
  dplyr::filter(as.character(expr_bin)=="Low") %>%
  dplyr::filter(all(table(cond) >= 2)) %>%
  rstatix::t_test(slope~cond) %>%
  rstatix::add_significance() %>%
  dplyr::mutate(y.position=shared_ylim_cl[2]*0.95)

p_cl_low <- make_halfeye(
  slopes_cluster_expr_lh %>% dplyr::filter(as.character(expr_bin)=="Low"),
  "cond","cond","cond",
  x_labels=c("WT","POLR2A deg"),
  title="Low Expression")+
  coord_cartesian(ylim=shared_ylim_cl)+
  stat_pvalue_manual(stat_cl_low, inherit.aes=FALSE, label="p.signif")

stat_cl_high <- slopes_cluster_expr_lh %>%
  dplyr::filter(as.character(expr_bin)=="High") %>%
  dplyr::filter(all(table(cond) >= 2)) %>%
  rstatix::t_test(slope~cond) %>%
  rstatix::add_significance() %>%
  dplyr::mutate(y.position=shared_ylim_cl[2]*0.95)

p_cl_high <- make_halfeye(
  slopes_cluster_expr_lh %>% dplyr::filter(as.character(expr_bin)=="High"),
  "cond","cond","cond",
  x_labels=c("WT","POLR2A deg"),
  title="High Expression")+
  coord_cartesian(ylim=shared_ylim_cl)+
  stat_pvalue_manual(stat_cl_high, inherit.aes=FALSE, label="p.signif")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())

p_cl_delta <- ggplot(delta_cluster_expr, aes(x=mean_TPM, y=delta))+
  geom_hline(yintercept=0, linetype="dashed",
             color="gray40", linewidth=0.5)+
  geom_line(color="#D55E00", linewidth=0.8)+
  geom_area(alpha=0.15, fill="#D55E00")+
  scale_x_log10(labels=scales::comma)+
  labs(x="Mean cluster TPM (log scale)", y="ΔS",
       title="ΔS (WT − POLR2A deg) vs Expression")+
  theme_classic(base_size=11)+
  theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
        axis.text=element_text(size=9,color="black"),
        axis.title=element_text(size=10,color="black"),
        panel.grid.major.y=element_line(color="grey92",linewidth=0.3),
        panel.grid.major.x=element_blank(),
        plot.margin=margin(4,16,8,8))

top_row    <- cowplot::plot_grid(p_cl_low, p_cl_high, nrow=1, rel_widths=c(1,1))
title_grob <- cowplot::ggdraw()+
  cowplot::draw_label(paste0("Cluster Scaling by Expression  |  ",DS_GROUP),
                      fontface="bold", size=13, x=0.01, hjust=0)
final <- cowplot::plot_grid(
  title_grob,
  cowplot::plot_grid(top_row, p_cl_delta, ncol=1, rel_heights=c(2,1)),
  ncol=1, rel_heights=c(0.05,1))
print(final)

# ══════════════════════════════════════════════════════════════════════════════
# 16. Rolling window S(d) and loess derivative — expression stratified
#     - 4 TPM quartile bins
#     - Sliding window: 250kb width, 50kb step, linear distance space
#     - Loess derivative: smooth P(s) then take numerical derivative
#     - Distance range: 5kb to 5Mb
#     - Plots: WT/pol2 lines + delta, plus diagnostics
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Section 16: Rolling window S(d) ═════════════════════════════════")

# ── 16a. Add 4-bin TPM quartiles to ei_binned ─────────────────────────────────
expressed_idx        <- which(ei_binned$TPM > 0)
ei_binned$tpm_bin4   <- NA_integer_
ei_binned$tpm_bin4[expressed_idx] <- as.integer(
  ntile(ei_binned$TPM[expressed_idx], 4L))

tpm_bin_labels <- c("Q1 (Low)","Q2","Q3","Q4 (High)")
tpm_bin_colors <- c("#4575B4","#74ADD1","#FDAE61","#D73027")

tpm_bin_medians <- ei_binned %>%
  dplyr::filter(!is.na(tpm_bin4)) %>%
  dplyr::group_by(tpm_bin4) %>%
  dplyr::summarise(median_TPM=median(TPM, na.rm=TRUE),
                   label=tpm_bin_labels[unique(tpm_bin4)],
                   .groups="drop")
message("TPM bin medians:")
print(tpm_bin_medians)

# ── 16b. Compute P(s) curves per condition × TPM bin ─────────────────────────
# Join tpm_bin4 onto contacts in memory
gene_contacts_list <- lapply(gene_contacts_list, function(ct) {
  ct %>%
    dplyr::left_join(
      ei_binned %>% dplyr::select(gene_id, tpm_bin4),
      by="gene_id"
    ) %>%
    data.table::as.data.table()
})

# Compute P(s) curves — full range to 5Mb (log10(5e6) ≈ 6.7)
message("Computing P(s) curves per TPM bin...")
ps_tpm_list <- lapply(names(CONDITIONS), function(n) {
  ct <- gene_contacts_list[[n]]
  ct_expr <- ct[!is.na(tpm_bin4)]
  ps <- compute_ps_curves(ct_expr, filter_col="tpm_bin4", max_distance=6.7)
  if (is.null(ps)) return(NULL)
  as.data.frame(ps) %>%
    dplyr::mutate(cond=CONDITIONS[[n]]$display)
})
ps_tpm <- bind_rows(Filter(Negate(is.null), ps_tpm_list)) %>%
  dplyr::mutate(
    cond    = factor(cond, levels=COND_LEVELS),
    tpm_bin4= as.integer(bin_val),
    tpm_label = factor(tpm_bin_labels[tpm_bin4],
                       levels=tpm_bin_labels)
  )

message("P(s) rows: ", nrow(ps_tpm))

# ── 16c. Sliding window slope ─────────────────────────────────────────────────
WINDOW_BP <- 250000L   # 250kb window
STEP_BP   <- 50000L    # 50kb step
MAX_DIST_BP <- 5000000L  # 5Mb

compute_sliding_slopes <- function(ps_df) {
  # ps_df: distance (log10), prob (log10), cond, tpm_bin4
  # Work in linear distance space
  ps_df <- ps_df %>%
    dplyr::mutate(dist_bp = 10^distance) %>%
    dplyr::filter(dist_bp >= 5000, dist_bp <= MAX_DIST_BP)
  
  window_starts <- seq(5000, MAX_DIST_BP - WINDOW_BP, by=STEP_BP)
  
  bind_rows(lapply(window_starts, function(ws) {
    we  <- ws + WINDOW_BP
    mid <- ws + WINDOW_BP/2
    sub <- ps_df %>% dplyr::filter(dist_bp >= ws & dist_bp <= we)
    if (nrow(sub) < 3) return(NULL)
    model <- tryCatch(
      lm(prob ~ distance, data=sub, na.action="na.exclude"),
      error=function(e) NULL)
    if (is.null(model)) return(NULL)
    data.frame(
      window_mid_bp = mid,
      window_mid_log= log10(mid),
      slope         = coef(model)[2],
      r2            = summary(model)$r.squared,
      n_points      = nrow(sub)
    )
  }))
}

message("Computing sliding window slopes...")
sliding_slopes <- ps_tpm %>%
  dplyr::group_by(cond, tpm_bin4, tpm_label) %>%
  dplyr::group_modify(~compute_sliding_slopes(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(slope), slope < 0)  # biologically valid slopes only

# ── 16d. Loess derivative ─────────────────────────────────────────────────────
compute_loess_deriv <- function(ps_df, span=0.3) {
  ps_df <- ps_df %>%
    dplyr::filter(distance >= log10(5000),
                  distance <= log10(MAX_DIST_BP)) %>%
    dplyr::arrange(distance)
  if (nrow(ps_df) < 10) return(NULL)
  fit <- tryCatch(
    loess(prob ~ distance, data=ps_df, span=span),
    error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  x_seq  <- seq(min(ps_df$distance), max(ps_df$distance), length.out=200)
  y_pred <- predict(fit, newdata=data.frame(distance=x_seq))
  # Numerical derivative: d(prob)/d(distance) — this IS the scaling exponent
  dy <- diff(y_pred) / diff(x_seq)
  dx <- (x_seq[-1] + x_seq[-length(x_seq)]) / 2
  data.frame(
    distance    = dx,
    dist_bp     = 10^dx,
    prob_smooth = predict(fit, newdata=data.frame(distance=dx)),
    slope       = dy
  )
}

message("Computing loess derivatives...")
loess_deriv <- ps_tpm %>%
  dplyr::group_by(cond, tpm_bin4, tpm_label) %>%
  dplyr::group_modify(~compute_loess_deriv(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(slope))

# ── 16e. Main plots — sliding window ─────────────────────────────────────────
wt_label   <- COND_LEVELS[1]
pol2_label <- COND_LEVELS[2]

# Top panel: S(d) curves per TPM bin, WT and pol2
p_slide_main <- ggplot(sliding_slopes,
                       aes(x=window_mid_bp/1000, y=slope,
                           color=tpm_label,
                           linetype=cond))+
  geom_hline(yintercept=0, linetype="dotted",
             color="gray50", linewidth=0.4)+
  geom_line(linewidth=0.7, alpha=0.85)+
  scale_color_manual(values=setNames(tpm_bin_colors, tpm_bin_labels),
                     name="TPM quartile")+
  scale_linetype_manual(values=c("solid","dashed"), name=NULL)+
  scale_x_continuous(labels=function(x) paste0(x,"kb"),
                     breaks=c(0,500,1000,2000,3000,4000,5000))+
  labs(x="Genomic distance", y="S(d)",
       title=paste0("Rolling window S(d)  |  ",DS_GROUP),
       subtitle=paste0("Window=250kb  |  Step=50kb  |  by expression quartile"))+
  theme_classic(base_size=13)+
  theme(plot.title=element_text(face="bold",size=14,hjust=0.5),
        plot.subtitle=element_text(color="grey40",size=10,hjust=0.5),
        legend.position="right",
        panel.grid.major.y=element_line(color="grey92",linewidth=0.3))

# Bottom panel: delta S(d) per TPM bin
delta_slide <- sliding_slopes %>%
  dplyr::select(window_mid_bp, slope, cond, tpm_bin4, tpm_label) %>%
  tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
  dplyr::mutate(delta = .data[[wt_label]] - .data[[pol2_label]]) %>%
  dplyr::filter(!is.na(delta))

p_slide_delta <- ggplot(delta_slide,
                        aes(x=window_mid_bp/1000, y=delta,
                            color=tpm_label))+
  geom_hline(yintercept=0, linetype="dashed",
             color="gray40", linewidth=0.5)+
  geom_line(linewidth=0.7)+
  scale_color_manual(values=setNames(tpm_bin_colors, tpm_bin_labels),
                     name="TPM quartile")+
  scale_x_continuous(labels=function(x) paste0(x,"kb"),
                     breaks=c(0,500,1000,2000,3000,4000,5000))+
  labs(x="Genomic distance", y="ΔS(d)",
       title="ΔS(d) = WT − POLR2A deg")+
  theme_classic(base_size=12)+
  theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
        legend.position="none",
        panel.grid.major.y=element_line(color="grey92",linewidth=0.3))

print(cowplot::plot_grid(p_slide_main, p_slide_delta,
                         ncol=1, rel_heights=c(1.5,1)))

# ── 16f. Main plots — loess derivative ───────────────────────────────────────

# Smart x-axis: kb below 1Mb, Mb at 1Mb and above
dist_axis <- function() {
  scale_x_continuous(
    breaks = c(100, 250, 500, 1000, 2000, 3000, 5000),
    labels = c("100kb","250kb","500kb","1Mb","2Mb","3Mb","5Mb")
  )
}

make_loess_main <- function(deriv_df, title_suffix="") {
  ggplot(deriv_df,
         aes(x=dist_bp/1000, y=slope,
             color=tpm_label, linetype=cond))+
    geom_hline(yintercept=0, linetype="dotted",
               color="gray50", linewidth=0.4)+
    geom_line(linewidth=0.8, alpha=0.85)+
    scale_color_manual(values=setNames(tpm_bin_colors, tpm_bin_labels),
                       name="TPM quartile")+
    scale_linetype_manual(values=c("solid","dashed"), name=NULL)+
    dist_axis()+
    labs(x="Genomic distance", y="S(d)",
         title=paste0("Loess S(d)  |  ",DS_GROUP,title_suffix),
         subtitle="span=0.3  |  by expression quartile")+
    theme_classic(base_size=13)+
    theme(plot.title=element_text(face="bold",size=14,hjust=0.5),
          plot.subtitle=element_text(color="grey40",size=10,hjust=0.5),
          legend.position="right",
          axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
}

make_loess_delta <- function(delta_df, title_suffix="") {
  ggplot(delta_df,
         aes(x=dist_bp/1000, y=delta, color=tpm_label))+
    geom_hline(yintercept=0, linetype="dashed",
               color="gray40", linewidth=0.5)+
    geom_line(linewidth=0.8)+
    scale_color_manual(values=setNames(tpm_bin_colors, tpm_bin_labels),
                       name="TPM quartile")+
    dist_axis()+
    labs(x="Genomic distance", y="ΔS(d)",
         title=paste0("ΔS(d) = WT − POLR2A deg",title_suffix))+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          legend.position="none",
          axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
}

delta_loess <- loess_deriv %>%
  dplyr::select(dist_bp, slope, cond, tpm_bin4, tpm_label) %>%
  tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
  dplyr::mutate(delta = .data[[wt_label]] - .data[[pol2_label]]) %>%
  dplyr::filter(!is.na(delta))

# All 4 quartiles
print(cowplot::plot_grid(
  make_loess_main(loess_deriv),
  make_loess_delta(delta_loess),
  ncol=1, rel_heights=c(1.5,1)))

# Q1 and Q4 only
loess_deriv_q14  <- loess_deriv  %>% dplyr::filter(tpm_bin4 %in% c(1L,4L))
delta_loess_q14  <- delta_loess  %>% dplyr::filter(tpm_bin4 %in% c(1L,4L))

print(cowplot::plot_grid(
  make_loess_main(loess_deriv_q14,  " — Q1 vs Q4"),
  make_loess_delta(delta_loess_q14, " — Q1 vs Q4"),
  ncol=1, rel_heights=c(1.5,1)))



# ── 16g. Diagnostic 1: P(s) curve with sliding window slope segments ──────────
message("Diagnostic plots...")

# One diagnostic per condition per TPM bin — too many to do all; do WT only Q1 and Q4
for (qbin in c(1L, 4L)) {
  ps_sub <- ps_tpm %>%
    dplyr::filter(cond==wt_label, tpm_bin4==qbin,
                  distance >= log10(5000),
                  distance <= log10(MAX_DIST_BP))
  
  slopes_sub <- sliding_slopes %>%
    dplyr::filter(cond==wt_label, tpm_bin4==qbin)
  
  # Build slope segment lines for each window
  seg_df <- bind_rows(lapply(seq_len(nrow(slopes_sub)), function(i) {
    row   <- slopes_sub[i,]
    x_lo  <- log10(row$window_mid_bp - WINDOW_BP/2)
    x_hi  <- log10(row$window_mid_bp + WINDOW_BP/2)
    x_lo  <- max(x_lo, min(ps_sub$distance))
    x_hi  <- min(x_hi, max(ps_sub$distance))
    # Predict y from slope + intercept — refit to get intercept
    sub   <- ps_sub %>% dplyr::filter(distance >= x_lo, distance <= x_hi)
    if (nrow(sub) < 3) return(NULL)
    model <- tryCatch(lm(prob~distance, data=sub), error=function(e) NULL)
    if (is.null(model)) return(NULL)
    data.frame(
      x=c(x_lo, x_hi),
      y=predict(model, newdata=data.frame(distance=c(x_lo,x_hi))),
      window_mid_bp=row$window_mid_bp,
      slope=row$slope
    )
  }))
  
  p_diag1 <- ggplot()+
    geom_point(data=ps_sub, aes(x=distance, y=prob),
               size=0.6, alpha=0.4, color="gray50")+
    geom_line(data=ps_sub, aes(x=distance, y=prob),
              color="gray50", linewidth=0.5)+
    geom_line(data=seg_df,
              aes(x=x, y=y, group=window_mid_bp,
                  color=slope),
              linewidth=1.2)+
    scale_color_distiller(palette="RdBu", direction=1,
                          name="S(d)")+
    scale_x_continuous(breaks=log10(c(5000,50000,100000,250000,500000,1e6,2e6,3e6,5e6)),
                       labels=c("5kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","5Mb"))+
    labs(x="Genomic distance (log scale)", y="log10 P(s)",
         title=paste0("Diagnostic: sliding window slopes  |  WT  |  ",
                      tpm_bin_labels[qbin]),
         subtitle="Each colored segment = one 250kb window slope")+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle=element_text(color="grey40",size=9,hjust=0.5),
          panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
  print(p_diag1)
}

# ── 16h. Diagnostic 2: P(s) with loess curve and derivative points ────────────
for (qbin in c(1L, 4L)) {
  ps_sub <- ps_tpm %>%
    dplyr::filter(cond==wt_label, tpm_bin4==qbin,
                  distance >= log10(5000),
                  distance <= log10(MAX_DIST_BP))
  deriv_sub <- loess_deriv %>%
    dplyr::filter(cond==wt_label, tpm_bin4==qbin)
  
  p_diag2_top <- ggplot(ps_sub, aes(x=distance, y=prob))+
    geom_point(size=0.6, alpha=0.4, color="gray50")+
    geom_line(color="gray50", linewidth=0.5)+
    geom_line(data=deriv_sub,
              aes(x=distance, y=prob_smooth),
              color=tpm_bin_colors[qbin], linewidth=1.2)+
    scale_x_continuous(breaks=log10(c(5000,50000,100000,250000,500000,1e6,2e6,3e6,5e6)),
                       labels=c("5kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","5Mb"))+
    labs(x=NULL, y="log10 P(s)",
         title=paste0("Diagnostic: loess P(s)  |  WT  |  ",tpm_bin_labels[qbin]),
         subtitle="Gray = raw  |  Colored = loess smooth")+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle=element_text(color="grey40",size=9,hjust=0.5),
          axis.text.x=element_blank(),
          panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
  
  p_diag2_bot <- ggplot(deriv_sub, aes(x=distance, y=slope))+
    geom_hline(yintercept=0, linetype="dashed",
               color="gray40", linewidth=0.5)+
    geom_point(aes(color=slope), size=1.5)+
    geom_line(color=tpm_bin_colors[qbin], linewidth=0.8)+
    scale_color_distiller(palette="RdBu", direction=1, guide="none")+
    scale_x_continuous(breaks=log10(c(5000,50000,100000,250000,500000,1e6,2e6,3e6,5e6)),
                       labels=c("5kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","5Mb"))+
    labs(x="Genomic distance (log scale)", y="S(d) derivative")+
    theme_classic(base_size=12)+
    theme(panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
  
  print(cowplot::plot_grid(p_diag2_top, p_diag2_bot,
                           ncol=1, rel_heights=c(1.5,1),
                           align="v", axis="lr"))
}


# ══════════════════════════════════════════════════════════════════════════════
# 17. Loess derivative S(d) — tandem gene clusters × expression
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Section 17: Cluster loess S(d) ══════════════════════════════════")

tpm_bin_labels <- c("Q1 (Low)","Q2","Q3","Q4 (High)")
tpm_bin_colors <- c("#4575B4","#74ADD1","#FDAE61","#D73027")

# ── 17a. Add tpm_bin4 to cluster_summary ─────────────────────────────────────
cluster_expressed_idx <- which(cluster_summary$mean_TPM > 0)
cluster_summary$tpm_bin4 <- NA_integer_
cluster_summary$tpm_bin4[cluster_expressed_idx] <- as.integer(
  ntile(cluster_summary$mean_TPM[cluster_expressed_idx], 4L))

message("Cluster TPM bin counts:")
print(table(cluster_summary$tpm_bin4, useNA="always"))

# ── 17b. Join tpm_bin4 onto cluster contacts in memory ───────────────────────
message("Joining tpm_bin4 onto cluster contacts...")
cluster_contacts_list <- lapply(cluster_contacts_list, function(ct) {
  data.table::as.data.table(ct) %>%
    dplyr::left_join(
      cluster_summary %>% dplyr::select(gene_id, tpm_bin4),
      by="gene_id"
    ) %>%
    data.table::as.data.table()
})

# ── 17c. Compute P(s) curves per condition × TPM bin ─────────────────────────
message("Computing cluster P(s) curves per TPM bin...")
ps_cluster_tpm_list <- lapply(names(CONDITIONS), function(n) {
  ct      <- cluster_contacts_list[[n]]
  ct_expr <- ct[!is.na(tpm_bin4)]
  ps      <- compute_ps_curves(ct_expr, filter_col="tpm_bin4", max_distance=6.7)
  if (is.null(ps)) return(NULL)
  as.data.frame(ps) %>%
    dplyr::mutate(cond=CONDITIONS[[n]]$display)
})
ps_cluster_tpm <- bind_rows(Filter(Negate(is.null), ps_cluster_tpm_list)) %>%
  dplyr::mutate(
    cond      = factor(cond, levels=COND_LEVELS),
    tpm_bin4  = as.integer(bin_val),
    tpm_label = factor(tpm_bin_labels[tpm_bin4], levels=tpm_bin_labels)
  )
message("Cluster P(s) rows: ", nrow(ps_cluster_tpm))

# ── 17d. Loess derivative ─────────────────────────────────────────────────────
message("Computing cluster loess derivatives...")
loess_deriv_cluster <- ps_cluster_tpm %>%
  dplyr::group_by(cond, tpm_bin4, tpm_label) %>%
  dplyr::group_modify(~compute_loess_deriv(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(slope))

delta_loess_cluster <- loess_deriv_cluster %>%
  dplyr::select(dist_bp, slope, cond, tpm_bin4, tpm_label) %>%
  tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
  dplyr::mutate(delta = .data[[wt_label]] - .data[[pol2_label]]) %>%
  dplyr::filter(!is.na(delta))

# ── 17e. Plots — all 4 quartiles ─────────────────────────────────────────────
print(cowplot::plot_grid(
  make_loess_main(loess_deriv_cluster,
                  title_suffix=" — Tandem clusters"),
  make_loess_delta(delta_loess_cluster,
                   title_suffix=" — Tandem clusters"),
  ncol=1, rel_heights=c(1.5,1)))

# ── 17f. Plots — Q1 and Q4 only ──────────────────────────────────────────────
loess_cluster_q14 <- loess_deriv_cluster %>% dplyr::filter(tpm_bin4 %in% c(1L,4L))
delta_cluster_q14 <- delta_loess_cluster  %>% dplyr::filter(tpm_bin4 %in% c(1L,4L))

print(cowplot::plot_grid(
  make_loess_main(loess_cluster_q14,
                  title_suffix=" — Tandem clusters  |  Q1 vs Q4"),
  make_loess_delta(delta_cluster_q14,
                   title_suffix=" — Tandem clusters  |  Q1 vs Q4"),
  ncol=1, rel_heights=c(1.5,1)))

# ── 17g. Diagnostic: P(s) + loess + derivative for Q1 and Q4 ─────────────────
for (qbin in c(1L, 4L)) {
  ps_sub    <- ps_cluster_tpm %>%
    dplyr::filter(cond==wt_label, tpm_bin4==qbin,
                  distance >= log10(5000),
                  distance <= log10(MAX_DIST_BP))
  deriv_sub <- loess_deriv_cluster %>%
    dplyr::filter(cond==wt_label, tpm_bin4==qbin)
  
  p_top <- ggplot(ps_sub, aes(x=distance, y=prob))+
    geom_point(size=0.6, alpha=0.4, color="gray50")+
    geom_line(color="gray50", linewidth=0.5)+
    geom_line(data=deriv_sub,
              aes(x=distance, y=prob_smooth),
              color=tpm_bin_colors[qbin], linewidth=1.2)+
    scale_x_continuous(
      breaks=log10(c(5000,50000,100000,250000,500000,1e6,2e6,3e6,5e6)),
      labels=c("5kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","5Mb"))+
    labs(x=NULL, y="log10 P(s)",
         title=paste0("Cluster P(s)  |  WT  |  ",tpm_bin_labels[qbin]),
         subtitle="Gray = raw  |  Colored = loess smooth")+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle=element_text(color="grey40",size=9,hjust=0.5),
          axis.text.x=element_blank(),
          panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
  
  p_bot <- ggplot(deriv_sub, aes(x=distance, y=slope))+
    geom_hline(yintercept=0, linetype="dashed",
               color="gray40", linewidth=0.5)+
    geom_point(aes(color=slope), size=1.5)+
    geom_line(color=tpm_bin_colors[qbin], linewidth=0.8)+
    scale_color_distiller(palette="RdBu", direction=1, guide="none")+
    scale_x_continuous(
      breaks=log10(c(5000,50000,100000,250000,500000,1e6,2e6,3e6,5e6)),
      labels=c("5kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","5Mb"))+
    labs(x="Genomic distance", y="S(d) derivative")+
    theme_classic(base_size=12)+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
  
  print(cowplot::plot_grid(p_top, p_bot, ncol=1, rel_heights=c(1.5,1),
                           align="v", axis="lr"))
}


# ══════════════════════════════════════════════════════════════════════════════
# 17. Loess derivative S(d) — tandem gene clusters, expression stratified
# ══════════════════════════════════════════════════════════════════════════════
message("\n══ Section 17: Cluster loess S(d) ═══════════════════════════════════")

# ── 17a. Add tpm_bin4 to cluster_summary ──────────────────────────────────────
cluster_expressed_idx          <- which(cluster_summary$mean_TPM > 0)
cluster_summary$tpm_bin4       <- NA_integer_
cluster_summary$tpm_bin4[cluster_expressed_idx] <- as.integer(
  ntile(cluster_summary$mean_TPM[cluster_expressed_idx], 4L))

message("Cluster TPM bin counts:")
print(table(cluster_summary$tpm_bin4, useNA="always"))

# ── 17b. Join tpm_bin4 onto cluster contacts in memory ────────────────────────
message("Joining tpm_bin4 onto cluster contacts...")
cluster_contacts_list <- lapply(cluster_contacts_list, function(ct) {
  ct %>%
    dplyr::left_join(
      cluster_summary %>% dplyr::select(gene_id, chr, tpm_bin4),
      by=c("gene_id","chr")
    ) %>%
    data.table::as.data.table()
})

# ── 17c. Compute P(s) curves per condition × tpm_bin4 ─────────────────────────
message("Computing cluster P(s) curves per TPM bin...")
ps_cluster_tpm_list <- lapply(names(CONDITIONS), function(n) {
  ct      <- cluster_contacts_list[[n]]
  ct_expr <- ct[!is.na(tpm_bin4)]
  ps      <- compute_ps_curves(ct_expr, filter_col="tpm_bin4",
                               max_distance=6.7)
  if (is.null(ps)) return(NULL)
  as.data.frame(ps) %>%
    dplyr::mutate(cond=CONDITIONS[[n]]$display)
})
ps_cluster_tpm <- bind_rows(Filter(Negate(is.null), ps_cluster_tpm_list)) %>%
  dplyr::mutate(
    cond      = factor(cond, levels=COND_LEVELS),
    tpm_bin4  = as.integer(bin_val),
    tpm_label = factor(tpm_bin_labels[tpm_bin4], levels=tpm_bin_labels)
  )
message("Cluster P(s) rows: ", nrow(ps_cluster_tpm))

# ── 17d. Loess derivative ─────────────────────────────────────────────────────
message("Computing cluster loess derivatives...")
loess_deriv_cluster <- ps_cluster_tpm %>%
  dplyr::group_by(cond, tpm_bin4, tpm_label) %>%
  dplyr::group_modify(~compute_loess_deriv(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(slope))

delta_loess_cluster <- loess_deriv_cluster %>%
  dplyr::select(dist_bp, slope, cond, tpm_bin4, tpm_label) %>%
  tidyr::pivot_wider(names_from=cond, values_from=slope, values_fn=mean) %>%
  dplyr::mutate(delta = .data[[wt_label]] - .data[[pol2_label]]) %>%
  dplyr::filter(!is.na(delta))

# ── 17e. Plots — all 4 quartiles ─────────────────────────────────────────────
print(cowplot::plot_grid(
  make_loess_main(loess_deriv_cluster,
                  " — Tandem clusters"),
  make_loess_delta(delta_loess_cluster,
                   " — Tandem clusters"),
  ncol=1, rel_heights=c(1.5,1)))

# ── 17f. Plots — Q1 and Q4 only ──────────────────────────────────────────────
loess_cluster_q14 <- loess_deriv_cluster %>%
  dplyr::filter(tpm_bin4 %in% c(1L,4L))
delta_cluster_q14 <- delta_loess_cluster %>%
  dplyr::filter(tpm_bin4 %in% c(1L,4L))

print(cowplot::plot_grid(
  make_loess_main(loess_cluster_q14,
                  " — Tandem clusters  |  Q1 vs Q4"),
  make_loess_delta(delta_cluster_q14,
                   " — Tandem clusters  |  Q1 vs Q4"),
  ncol=1, rel_heights=c(1.5,1)))

# ── 17g. Diagnostic: loess P(s) + derivative for Q1 and Q4 ──────────────────
for (qbin in c(1L, 4L)) {
  ps_sub    <- ps_cluster_tpm %>%
    dplyr::filter(cond==wt_label, tpm_bin4==qbin,
                  distance >= log10(5000),
                  distance <= log10(MAX_DIST_BP))
  deriv_sub <- loess_deriv_cluster %>%
    dplyr::filter(cond==wt_label, tpm_bin4==qbin)
  
  p_top <- ggplot(ps_sub, aes(x=distance, y=prob))+
    geom_point(size=0.6, alpha=0.4, color="gray50")+
    geom_line(color="gray50", linewidth=0.5)+
    geom_line(data=deriv_sub,
              aes(x=distance, y=prob_smooth),
              color=tpm_bin_colors[qbin], linewidth=1.2)+
    scale_x_continuous(
      breaks=log10(c(5000,50000,100000,250000,500000,1e6,2e6,3e6,5e6)),
      labels=c("5kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","5Mb"))+
    labs(x=NULL, y="log10 P(s)",
         title=paste0("Cluster loess P(s)  |  WT  |  ",tpm_bin_labels[qbin]),
         subtitle="Gray = raw  |  Colored = loess smooth")+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle=element_text(color="grey40",size=9,hjust=0.5),
          axis.text.x=element_blank(),
          panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
  
  p_bot <- ggplot(deriv_sub, aes(x=distance, y=slope))+
    geom_hline(yintercept=0, linetype="dashed",
               color="gray40", linewidth=0.5)+
    geom_point(aes(color=slope), size=1.5)+
    geom_line(color=tpm_bin_colors[qbin], linewidth=0.8)+
    scale_color_distiller(palette="RdBu", direction=1, guide="none")+
    scale_x_continuous(
      breaks=log10(c(5000,50000,100000,250000,500000,1e6,2e6,3e6,5e6)),
      labels=c("5kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","5Mb"))+
    labs(x="Genomic distance", y="S(d) derivative")+
    theme_classic(base_size=12)+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major.y=element_line(color="grey92",linewidth=0.3))
  
  print(cowplot::plot_grid(p_top, p_bot, ncol=1,
                           rel_heights=c(1.5,1),
                           align="v", axis="lr"))
}

message("\nEI scaling analysis complete.")