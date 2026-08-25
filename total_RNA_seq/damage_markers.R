##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Splicing / Apoptosis / UPR pathway gene panels. Fold-change
# heatmaps and MA plot (ActD, POLR2A vs. CTRL) for three curated
# gene panels: Core Splicing, Core Apoptosis, Core UPR. Sources
# functions.R for shared plotting utilities.
#
# Figure output: SI Fig. 2A (confirmed -- literal title match
# "Core Pathway Gene Expression Is Stable Upon Pol2 Depletion";
# exact color scheme Core Splicing=blue #0072B2, Core
# Apoptosis=orange #E69F00, Core UPR=pink #CC79A7).
##################################################################

# upregulated color: "#D46780"
# downregulated color: "#798234"

library(here)
root <- here::here()

dir.rsem <- here::here("counts/RSEM")

## get custom functions
source(here::here("scripts/R/functions.R"))

# import packages
library("ggplot2")
library("dplyr")
library("magrittr")
library("tximport")
library(forcats)
library(stringr)
library(DESeq2)

#########################################################################################
# Set Directory, pass in counts, filter counts for DGE analysis (RSEM counts)
#########################################################################################

sample <- read.table(file.path(dir.rsem, "metadata.csv"), sep= ",", header= TRUE)

sample <- sample[-c(6),] ## remove pol 1 r3
sample <- sample[-c(7),] ## remove pol 2 r1

newcolnames <-sample$Alias
sample

## Import files
files <- file.path(dir.rsem, paste0(sample$ID))
files

## Imports RSEM TPM files into a list genes.results
cts <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

## Move col names from sample info to counts matrix
colnames(cts$counts)<-newcolnames

#Replace any transcripts with length 0 with a 1
cts$length[cts$length == 0] <- 1

# logical test to see if ^ is true now
identical(sample$Alias, colnames(cts$counts))

# Store the file names in the sample info matrix as row names
row.names(sample) <- sample$Alias
head(sample, n=10)

dput(as.character(unique(sample$Conditions)))

#factor out dependencies
sample$Conditions <- factor(sample$Conditions, levels = c("CTRL_0hrs","ActD_6hrs","Pol1_6hrs","Pol2_6hrs"))
sample<- mutate(sample, Conditions = fct_relevel(Conditions, c( "CTRL_0hrs","ActD_6hrs","Pol1_6hrs","Pol2_6hrs")))

sample$Label <- factor(sample$Label, levels = dput(as.character(unique(sample$Label))))

#########################################################################################
# Create DESeq data object and get normalized counts
#########################################################################################

# Create DDS object using RSEM counts
dds <- DESeqDataSetFromTximport(cts, colData = sample, design = ~ Conditions)

# Find size factors to normalize differences in sequencing depth across samples
dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- data.frame(counts(dds, normalized=TRUE))

#########################################################################################
# Differential expression analysis and 'results' matrices
#########################################################################################

## Generate DGE DESeqDataSet data frame
dds.res <- DESeq(dds)

## Which contrasts to use?
resultsNames(dds.res)

## Generate a results matrix for Treatment vs. Control with an alpha of 0.05
res.actd <- results(dds.res, contrast = c("Conditions","ActD_6hrs", "CTRL_0hrs"), alpha = 0.05)
res.pol2 <- results(dds.res, contrast = c("Conditions", "Pol2_6hrs", "CTRL_0hrs"), alpha = 0.05)

#########################################################################################
# Annotating the results using Annotation Hub
#########################################################################################

ts = F
res.actd <- Annotate_results(res.actd, transcript=ts)
res.pol2 <- Annotate_results(res.pol2, transcript=ts)

#########################################################################################
# Core gene lists for Apoptosis, Splicing, and UPR panels (fallback if full lists missing)
#########################################################################################

gene.ls.1 <- c("U1", "U2", "U4", "U6", "U5", "RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", ## Key Splicing Genes
               "SNRNP48", "RNPC3", "PHF5A", "SF3B1", "SF3B2", "SF3B3", "SF3B4", "SF3B5", "SF3B6",
               "ZRSR2", "U2AF1")
gene.ls.2 <- c("APAF1", "ATM", "BAD", "BAX", "TP53", "CASP3", "CASP8", "FAS", "FADD", "CASP3", "CASP6", "CASP7") ## Key Apoptosis Genes
gene.ls.3 <- c("ERN1", "ATF6", "EIF2AK3", "XBP1", "ATF4", "DDIT3", "HSPA5", "HSP90B1", "PPP1R15A") ## Key UPR Genes

#########################################################################################
# Standard DE gene call, used everywhere in this script
#########################################################################################

## padj < 0.05 & |log2FC| > 1 -- applied for the heatmap, the summary counts, and any
## other place genes get called "differentially expressed". Explicitly NA-safe: any NA
## input returns FALSE rather than propagating NA (which would otherwise poison sum()
## downstream if even one gene has an NA padj/log2FC).
is_sig <- function(padj, log2fc) {
  ifelse(is.na(padj) | is.na(log2fc), FALSE, padj < 0.05 & abs(log2fc) > 1)
}

#########################################################################################
# Fold-change heatmap: Apoptosis / Splicing / UPR core genes, ActD & Pol2 vs CTRL
#########################################################################################

build_fc_df <- function(gene.ls, names, res.list, res.names){
  out <- list()
  for (i in seq_along(gene.ls)){
    y <- data.frame(genes = gene.ls[[i]])
    for (j in seq_along(res.list)){
      x <- res.list[[j]] %>% dplyr::filter(symbol %in% y$genes) %>%
        dplyr::select(symbol, baseMean, log2FoldChange, padj) %>%
        dplyr::mutate(conds = names[i], comparison = res.names[j])
      out[[length(out)+1]] <- x
    }
  }
  df <- do.call(rbind, out)
  ## Some symbols (esp. snRNA families like U1/U2/U4/U6/U5, RNU4ATAC, RNU6ATAC, RNU11,
  ## RNU12) map to many genomic loci/pseudogene copies sharing the same symbol. Keep
  ## only the most-expressed locus per symbol so each gene contributes one row instead
  ## of one row per pseudogene copy.
  df <- df %>%
    dplyr::group_by(symbol, conds, comparison) %>%
    dplyr::slice_max(baseMean, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  df
}

fc.df <- build_fc_df(
  gene.ls = list(gene.ls.1, gene.ls.2, gene.ls.3),
  names = c("Core Splicing", "Core Apoptosis", "Core UPR"),
  res.list = list(res.actd, res.pol2),
  res.names = c("ActD", "Pol2")
)

fc.df$comparison <- factor(fc.df$comparison, levels = c("ActD", "Pol2"))
fc.df$sig <- ifelse(is_sig(fc.df$padj, fc.df$log2FoldChange), "sig", "ns")

p <- ggplot(fc.df, aes(x = symbol, y = comparison, fill = log2FoldChange)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(data = subset(fc.df, sig == "sig"), aes(x = symbol, y = comparison),
             shape = 8, size = 1.2, color = "black", inherit.aes = FALSE) +
  scale_fill_gradient2(low = "#798234", mid = "white", high = "#D46780",
                       midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                       name = "Log2FC\nvs CTRL") +
  facet_grid(. ~ conds, scales = "free_x", space = "free_x") +
  theme_classic() + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(p)

#########################################################################################
# MA plot: Pol2 vs CTRL, all genes, colored by pathway category
#########################################################################################

## Use the broader annotated gene lists (same source/paths as the original GO-list
## panels) when available; fall back to the small "core" lists per category otherwise.
category_files <- list(
  "Core Splicing"   = here::here("Annotations/genelists/KEGG_SPLICEOSOME.txt"),
  "Core Apoptosis"  = here::here("Annotations/genelists/KEGG_APOPTOSIS.txt"),
  "Core UPR"        = here::here("Annotations/genelists/MERGED_UPR.txt")
)

category_core <- list(
  "Core Splicing"  = gene.ls.1,
  "Core Apoptosis" = gene.ls.2,
  "Core UPR"       = gene.ls.3
)

panel_genes <- list()
for (cat in names(category_files)){
  f <- category_files[[cat]]
  if (file.exists(f)){
    y <- read.table(f, sep = "\t")
    ## DIAGNOSTIC -- how many rows/cols did we actually read, and how many of those
    ## entries match a real gene symbol in res.pol2? A file that's one gene per line
    ## should show ncol = 1 and a match count close to nrow. If match count is much
    ## smaller than nrow, or ncol > 1, the file format isn't what this loop assumes.
    message(paste0(cat, ": read ", nrow(y), " row(s), ", ncol(y), " column(s) from ", basename(f)))
    colnames(y)[1] <- "genes"
    n.match <- sum(y$genes %in% data.frame(res.pol2)$symbol)
    message(paste0(cat, ": ", n.match, " of ", nrow(y), " entries match a symbol in res.pol2"))
    panel_genes[[cat]] <- data.frame(symbol = y$genes, conds = cat)
  } else {
    message(paste0(cat, ": annotation file not found, falling back to core gene list (", length(category_core[[cat]]), " genes)"))
    panel_genes[[cat]] <- data.frame(symbol = category_core[[cat]], conds = cat)
  }
}
panel_genes <- do.call(rbind, panel_genes) %>% dplyr::distinct(symbol, .keep_all = TRUE)

ma.df <- data.frame(res.pol2) %>%
  dplyr::left_join(panel_genes, by = "symbol") %>%
  dplyr::mutate(conds = ifelse(is.na(conds), "Background", conds))

## Trim the y-axis (and the ridge-plot x-axis below) to the bulk of the distribution --
## a handful of extreme outliers were otherwise creating a lot of dead space. Computed
## from the data rather than hardcoded, so it adapts if the outliers change.
ma.range <- quantile(ma.df$log2FoldChange, probs = c(0.005, 0.995), na.rm = TRUE)
ma.pad <- diff(ma.range) * 0.15
ma.ylim <- c(ma.range[1] - ma.pad, ma.range[2] + ma.pad)
ma.ylim

p.ma <- ggplot() +
  geom_point(data = subset(ma.df, conds == "Background"), aes(x = log10(baseMean + 1), y = log2FoldChange),
             color = "grey80", size = 0.6, alpha = 0.5) +
  geom_point(data = subset(ma.df, conds != "Background"), aes(x = log10(baseMean + 1), y = log2FoldChange, color = conds),
             size = 1.6) +
  scale_color_manual(values = c("Core Splicing" = "#0072B2", "Core Apoptosis" = "#E69F00", "Core UPR" = "#CC79A7"), name = "") +
  theme_classic() + xlab("Log10 Mean Expression") + ylab("Log2FC (Pol2 vs CTRL)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(ylim = ma.ylim) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))

#########################################################################################
# Ridge plot: LFC distribution per category, paired below the MA plot
#########################################################################################

library(ggridges)
library(patchwork)

ridge.counts <- ma.df %>%
  dplyr::filter(conds != "Background") %>%
  dplyr::count(conds, name = "n")
ridge.counts

p.ridge <- ggplot(subset(ma.df, conds != "Background"), aes(x = log2FoldChange, y = conds, fill = conds)) +
  geom_density_ridges(alpha = 0.85, color = "white", scale = 1.1,
                      quantile_lines = TRUE, quantiles = 4, vline_size = 0.7, vline_color = "black") +
  scale_fill_manual(values = c("Core Splicing" = "#0072B2", "Core Apoptosis" = "#E69F00", "Core UPR" = "#CC79A7"), guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(data = ridge.counts, aes(x = ma.ylim[2] - diff(ma.ylim) * 0.05, y = conds, label = paste0("N = ", n)),
            inherit.aes = FALSE, hjust = 0, vjust = -1.2, size = 3.5) +
  coord_cartesian(xlim = ma.ylim) +
  scale_x_continuous(breaks = seq(floor(ma.ylim[1]), ceiling(ma.ylim[2]), by = 1)) +
  theme_classic() + xlab("Log2FC (Pol2 vs CTRL)") + ylab("") +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))

print(p.ma / p.ridge +
        patchwork::plot_layout(heights = c(2, 1)) +
        patchwork::plot_annotation(title = "Core Pathway Gene Expression Is Stable Upon Pol2 Depletion") &
        theme(plot.title = element_text(size = 18, face = "bold")))

#########################################################################################
# Summary counts: significant genes per panel / comparison
#########################################################################################

fc.summary <- fc.df %>%
  dplyr::group_by(conds, comparison) %>%
  dplyr::summarise(n_sig = sum(sig == "sig", na.rm = TRUE), n_total = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0(n_sig, "/", n_total, " padj<0.05 & |log2FC|>1"))

print(fc.summary)

#########################################################################################
# RNA log2FC vs quantified western band intensity (needs your densitometry values)
#########################################################################################

## Fill in with your own quantified band intensities (e.g. normalized to loading
## control, expressed as log2FC vs Auxin- for each blot). One row per gene x
## comparison. This is a scaffold only -- band_log2fc is NA until you supply real values.
western.intensity <- data.frame(
  symbol = c("TP53", "XBP1"),          ## e.g. p53 blot, spliced XBP1 blot
  comparison = c("Pol2", "Pol2"),      ## which RNA-seq contrast this blot matches
  band_log2fc = c(NA, NA)              ## <-- plug in your quantified values here
)

western.plot.df <- fc.df %>%
  dplyr::filter(comparison == "Pol2") %>%
  dplyr::inner_join(western.intensity, by = c("symbol", "comparison"))

p.west <- ggplot(western.plot.df, aes(x = log2FoldChange, y = band_log2fc, label = symbol)) +
  geom_point(size = 2) + geom_text(vjust = -1) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() + xlab("RNA Log2FC (vs CTRL)") + ylab("Western Log2FC (band intensity)")
print(p.west)