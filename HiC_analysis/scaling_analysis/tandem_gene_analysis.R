##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Diagnostic only -- genome-wide inter-gene gap distribution, used
# to justify the 50kb MAX_GAP tandem-cluster threshold used by
# scaling_analysis.R's short-gene (<100kb) cluster branch
# (exon/intergenic-ratio analysis, separate from the long-gene
# exon/intron branch used for Fig. 4F-H).
#
# Figure output: none -- not a figure-generating script.
##################################################################
library(GenomicFeatures)
library(dplyr)
library(ggplot2)

GTF_FILE   <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
STD_CHROMS <- paste0("chr", c(1:22, "X"))

message("Building TxDb...")
txdb     <- makeTxDbFromGFF(GTF_FILE, format="gtf")
genes_gr <- genes(txdb) %>% GenomicRanges::reduce()
seqlevelsStyle(genes_gr) <- "UCSC"
genes_gr <- keepSeqlevels(genes_gr,
                          intersect(seqlevels(genes_gr), STD_CHROMS),
                          pruning.mode="coarse")

# Compute inter-gene gaps
genes_df <- as.data.frame(genes_gr) %>%
  dplyr::rename(chr=seqnames) %>%
  dplyr::arrange(chr, start) %>%
  dplyr::group_by(chr) %>%
  dplyr::mutate(
    gap = start - dplyr::lag(end) - 1L
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(gap), gap >= 0)

message("Total inter-gene gaps: ", nrow(genes_df))
message("Gap summary:")
print(quantile(genes_df$gap, probs=c(0,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1), na.rm=TRUE))

# Count tandem arrays at different thresholds
thresholds <- c(5000, 10000, 25000, 50000, 100000, 200000, 500000)
array_counts <- sapply(thresholds, function(thresh) {
  genes_df %>%
    dplyr::mutate(in_array = gap <= thresh) %>%
    dplyr::summarise(n_pairs=sum(in_array, na.rm=TRUE)) %>%
    pull(n_pairs)
})
threshold_df <- data.frame(threshold_kb=thresholds/1000, n_pairs=array_counts)
message("\nPairs within threshold:")
print(threshold_df)

# Plot gap distribution
p1 <- ggplot(genes_df %>% filter(gap <= 1e6),
             aes(x=gap)) +
  geom_histogram(bins=100, fill="#0072B2", color="white",
                 linewidth=0.15, alpha=0.85) +
  geom_vline(xintercept=c(5000,10000,25000,50000,100000),
             linetype="dashed", color=c("#E69F00","#56B4E9","#009E73","#D55E00","#CC79A7"),
             linewidth=0.7) +
  annotate("text",
           x=c(5000,10000,25000,50000,100000)*1.05,
           y=Inf,
           label=c("5kb","10kb","25kb","50kb","100kb"),
           hjust=0, vjust=1.5, size=3.2,
           color=c("#E69F00","#56B4E9","#009E73","#D55E00","#CC79A7"),
           fontface="bold") +
  scale_x_log10(
    breaks=c(100,500,1000,5000,10000,50000,100000,500000,1e6),
    labels=c("100bp","500bp","1kb","5kb","10kb","50kb","100kb","500kb","1Mb")
  ) +
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.05))) +
  labs(x="Inter-gene gap (log scale)",
       y="Number of gene pairs",
       title="Distribution of inter-gene gaps",
       subtitle=paste0("n = ", scales::comma(nrow(genes_df)),
                       " consecutive gene pairs  |  autosomes + chrX  |  gaps ≤ 1Mb shown")) +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(face="bold"))
print(p1)

# Plot number of tandem pairs by threshold
p2 <- ggplot(threshold_df, aes(x=threshold_kb, y=n_pairs)) +
  geom_line(color="#0072B2", linewidth=0.8) +
  geom_point(color="#0072B2", size=3) +
  geom_text(aes(label=scales::comma(n_pairs)),
            vjust=-0.8, size=3.2) +
  scale_x_log10(breaks=thresholds/1000,
                labels=paste0(thresholds/1000,"kb")) +
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.15))) +
  labs(x="Max gap threshold (log scale)",
       y="Number of tandem gene pairs",
       title="Tandem gene pairs by gap threshold",
       subtitle="Number of consecutive gene pairs within each threshold") +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(face="bold"))
print(p2)


# Summary of what the data shows
message("=== Gap Distribution Summary ===")
message("Total gene pairs: 43,174")
message("")
message("Key percentiles:")
message("  10th: ~891 bp   — very tight clusters")
message("  25th: ~4.4 kb   — compact gene arrays")
message("  50th: ~14.8 kb  — median gap")
message("  75th: ~41.2 kb  — most pairs within TAD range")
message("  90th: ~93.5 kb  — upper range of loop domains")
message("  95th: ~145 kb")
message("  99th: ~306 kb")
message("")
message("=== Tandem pairs captured at each threshold ===")
message("  5kb:   11,740 pairs (27% of all pairs)")
message("  10kb:  17,620 pairs (41%)")
message("  25kb:  27,255 pairs (63%)")
message("  50kb:  34,190 pairs (79%)")
message("  100kb: 39,223 pairs (91%)")
message("  200kb: 41,971 pairs (97%)")
message("  500kb: 43,052 pairs (>99%)")
message("")
message("=== Inflection points ===")
pct_gains <- diff(c(0, 11740, 17620, 27255, 34190, 39223, 41971, 43052)) / 43174 * 100
thresholds <- c(5, 10, 25, 50, 100, 200, 500)
for (i in seq_along(thresholds))
  message("  ",thresholds[i],"kb → ",thresholds[i+1],"kb: +",
          round(pct_gains[i+1],1),"% additional pairs")
