##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Shared plotting/analysis utilities. Sourced by other scripts in
# this folder (e.g. damage_markers.R), not run standalone. Includes
# DE_Vol_Plot() (volcano plot generator).
#
# Figure output: none directly -- shared utility functions.
##################################################################

DE_Vol_Plot <- function(results, plt.title, lfc.up=0.58, lfc.dwn=-0.58,pvadj=0.05,
                        y.limits, y.inc, x.limits=c(-7,7),dncol="#D46780", upcol="#798234") {
  require(ggrepel)
  require(ggplot2)
  
  ## Generate data frame
  genes.rna <- data.frame(symbol = results$symbol,lfc = results$log2FoldChange, padj=results$padj, base=results$baseMean)
  
  resultsDN <- as.data.frame(genes.rna[which(genes.rna$lfc < lfc.dwn),])
  resultsDN <- as.data.frame(resultsDN[which(resultsDN$padj < pvadj),])
  length(rownames(resultsDN)) ## Gives count of down regulated
  
  resultsUP <- as.data.frame(genes.rna[which(genes.rna$lfc > lfc.up),])
  resultsUP <- as.data.frame(resultsUP[which(resultsUP$padj < pvadj),])
  length(rownames(resultsUP)) ## Gives count of up regulated  
  
  ## add a column of NAs
  genes.rna$diffexpressed <- "NO"
  ## If log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  genes.rna$diffexpressed[genes.rna$lfc > lfc.up & genes.rna$padj < pvadj] <- "Upregulated"
  ## If log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  genes.rna$diffexpressed[genes.rna$lfc < lfc.dwn & genes.rna$padj < pvadj] <- "Downregulated"
  
  ## Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  genes.rna$delabel <- NA
  genes.rna$delabel[genes.rna$diffexpressed != "NO"] <- genes.rna$symbol[genes.rna$diffexpressed != "NO"]
  
  ## Jitter labeling
  pos <- position_jitter(width = 0.2, seed = 1)
  
  vplot <- genes.rna %>% mutate(mean.intensity = base/max(base)) %>% 
    ggplot(aes(x=lfc, y=-log10(padj), col=diffexpressed, size = mean.intensity , label=delabel
    )) +
    geom_point(position = pos, alpha=0.7)+ theme_classic()+
    scale_color_manual(name="Upregulated", values=c(dncol, "lightgrey", upcol))+
    geom_vline(xintercept=c(lfc.dwn, lfc.up), col="red", lty = "dashed") +
    geom_hline(yintercept=-log10(pvadj), col="red", lty = "dashed")+ 
    geom_text_repel(segment.color = 'transparent', position = pos, size = 3, xlim  = c(-6,6)) +
    labs(title=plt.title, subtitle=paste(paste("Downregulated:",length(rownames(resultsDN))), paste("| Upregulated:",length(rownames(resultsUP))))) +
    scale_y_continuous(name="-Log 10 Padj", limits=c(0,y.limits), breaks=seq(0,y.limits, by =y.inc))+
    scale_x_continuous(name="Log 2 Fold Change", limits=x.limits, breaks=c(-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7))+
    theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) 
  
  return(vplot)
  
}

DE_Vol_Plot_Facet <- function(results, plt.title, lfc.up=0.58, lfc.dwn=-0.58,pvadj=0.05,
                        y.limits, y.inc, x.limits=c(-7,7),dncol="#D46780", upcol="#798234") {
  require(ggrepel)
  require(ggplot2)
  
  ## Generate data frame
  genes.rna <- data.frame(symbol = results$symbol,lfc = results$log2FoldChange, padj=results$padj, base=results$baseMean, conds = results$conds)
  
  resultsDN <- as.data.frame(genes.rna[which(genes.rna$lfc < lfc.dwn),])
  resultsDN <- as.data.frame(resultsDN[which(resultsDN$padj < pvadj),])
  length(rownames(resultsDN)) ## Gives count of down regulated
  
  resultsUP <- as.data.frame(genes.rna[which(genes.rna$lfc > lfc.up),])
  resultsUP <- as.data.frame(resultsUP[which(resultsUP$padj < pvadj),])
  length(rownames(resultsUP)) ## Gives count of up regulated  
  
  ## add a column of NAs
  genes.rna$diffexpressed <- "NO"
  ## If log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  genes.rna$diffexpressed[genes.rna$lfc > lfc.up & genes.rna$padj < pvadj] <- "Upregulated"
  ## If log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  genes.rna$diffexpressed[genes.rna$lfc < lfc.dwn & genes.rna$padj < pvadj] <- "Downregulated"
  
  ## Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  genes.rna$delabel <- NA
  genes.rna$delabel[genes.rna$diffexpressed != "NO"] <- genes.rna$symbol[genes.rna$diffexpressed != "NO"]
  
  ## Jitter labeling
  pos <- position_jitter(width = 0.2, seed = 1)
  
  vplot <- genes.rna %>% mutate(mean.intensity = base/max(base)) %>% 
    ggplot(aes(x=lfc, y=-log10(padj), col=diffexpressed, size = mean.intensity , label=delabel
    )) +
    geom_point(position = pos, alpha=0.7)+ theme_classic()+
    scale_color_manual(name="Upregulated", values=c(dncol, "lightgrey", upcol))+
    geom_vline(xintercept=c(lfc.dwn, lfc.up), col="red", lty = "dashed") +
    geom_hline(yintercept=-log10(pvadj), col="red", lty = "dashed")+ 
    geom_text_repel(segment.color = 'transparent', position = pos, size = 3, xlim  = c(-6,6)) +
    scale_y_continuous(name="-Log 10 Padj", limits=c(0,y.limits), breaks=seq(0,y.limits, by =y.inc))+
    scale_x_continuous(name="Log 2 Fold Change", limits=x.limits, breaks=c(-6,-4,-2,0,2,4,6))+
    theme( legend.position="none",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) 
  
  return(vplot)
  
}

Annotate_results <- function(results, transcript){
  
  require("AnnotationHub")
  require(AnnotationDbi)
  
  ah <- AnnotationHub() ## Connect to Annotation Hub
  
  db <- query(ah, c("EnsDb","Homo sapiens","111"))[[1]] ## Make a query to your organism/database
  results <- data.frame(results)
  
  if (transcript == F) {
    
    annotations <- genes(db, return.type = "data.frame")
    annot <- annotations %>%
      dplyr::select(gene_id, canonical_transcript, symbol,gene_biotype, seq_name, gene_seq_start, gene_seq_end) %>%
      dplyr::filter(gene_id %in% rownames(results)) ## genes.results
    
    ## Order results to match them
    order <- match(rownames(results), annot$gene_id) ## genes.results
    annot <- annot[order,]
    
    ## Annotations
    results$transcript <- annot$canonical_transcript
    results$symbol <- annot$symbol
    results$gene_biotype <- annot$gene_biotype
    results$seqnames <- annot$seq_name
    results$start <- annot$gene_seq_start
    results$end <- annot$gene_seq_end
    
  } else {
    
    keys <- keys(db, c("TXID"))
    annotations <- AnnotationDbi::select(db, keys, c("GENEID", "TXNAME", "SYMBOL", "GENENAME","GENEBIOTYPE","TXBIOTYPE", "SEQNAME", "TXSEQSTART","TXSEQEND","EXONID","EXONSEQSTART", "EXONSEQEND"), "TXID")
    
    annot <- annotations %>%
      dplyr::select("GENEID", "TXID", "TXNAME", "SYMBOL", "GENENAME","GENEBIOTYPE","TXBIOTYPE", "SEQNAME", "TXSEQSTART","TXSEQEND","EXONID","EXONSEQSTART", "EXONSEQEND") %>%
      dplyr::filter(TXID %in% rownames(results)) ## genes.results
    
    ## Order results to match them
    order <- match(rownames(results), annot$TXID) ## genes.results
    annot <- annot[order,]
    
    ## Annotations
    results$transcript <- annot$TXID
    results$symbol <- annot$GENENAME
    results$gene_biotype <- annot$GENEBIOTYPE
    results$tx_biotype <- annot$TXBIOTYPE
    results$seqnames <- annot$SEQNAME
    results$start <- annot$TXSEQSTART
    results$end <- annot$TXSEQEND
    
  }
  return(results)
}

Filter_Results <- function(results, lfc.up=0.58, lfc.dwn=-0.58,pvadj=0.05) {
  
  resultsDN <- as.data.frame(results[which(results$log2FoldChange < lfc.dwn),])
  resultsDN <- as.data.frame(resultsDN[which(resultsDN$padj <= pvadj),])
  length(rownames(resultsDN)) ## Gives count of down regulated
  
  resultsUP <- as.data.frame(results[which(results$log2FoldChange  > lfc.up),])
  resultsUP <- as.data.frame(resultsUP[which(resultsUP$padj <= pvadj),])
  length(rownames(resultsUP)) ## Gives count of up regulated  
  
  ## add a column of NAs
  results$diffexpressed <- "NO"
  ## If log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results$diffexpressed[results$log2FoldChange > lfc.up & results$padj <= pvadj] <- "Upregulated"
  ## If log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results$diffexpressed[results$log2FoldChange  < lfc.dwn & results$padj <= pvadj] <- "Downregulated"
  
  ## Make counts vectors
  results$up <- ifelse(results$diffexpressed == "Upregulated", 1, 0)
  results$down <- ifelse(results$diffexpressed == "Downregulated", 1, 0)
  
  ## Make lfc vectors
  results$up.lfc <- ifelse(results$diffexpressed == "Upregulated", results$log2FoldChange, 0)
  results$down.lfc <- ifelse(results$diffexpressed == "Downregulated", results$log2FoldChange, 0)
  
  ## Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  results$delabel <- NA
  results$delabel[results$diffexpressed != "NO"] <- results$symbol[results$diffexpressed != "NO"]
  
  return(results)
  
}