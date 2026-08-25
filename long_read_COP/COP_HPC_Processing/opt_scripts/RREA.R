#!/usr/bin/env Rscript
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Read-Relative Exon Annotations (RREA), data.table-optimised
# implementation. Maps aligned long reads to exon coordinates from
# a GTF annotation, producing a read-relative exon coordinate table
# for downstream analysis (per-chromosome, per-read exon overlap
# and intron gap extraction).
#
# Key optimisations:
#   1. BAM streamed in chunks (yieldSize) -- never fully loaded into memory
#   2. findOverlaps() called once per chunk (vectorised)
#   3. CIGAR parsing only for reads with exon overlaps
#   4. All position lookups done with data.table keyed binary search (no R loops)
#   5. Per-chromosome parallelism via BiocParallel SnowParam
#
# Figure output: none directly -- upstream data-prep for
# long_read_COP/RREA_analysis.R (SI Fig. 13G-H).
##################################################################

# ── Usage ─────────────────────────────────────────────────────────────────────
usage <- "
RREA — Read-Relative Exon Annotations
======================================
Maps aligned long reads to exon coordinates from a GTF annotation,
producing a read-relative exon coordinate table for downstream analysis.

Usage:
  Rscript RREA.R --bam <file> --gtf <file> --out <file> [options]

Required arguments:
  --bam   STR   Input indexed BAM file (.bam, must have .bai index)
  --gtf   STR   Genome annotation file in GTF format
  --out   STR   Output path for tab-delimited results file

Optional arguments:
  --workers INT  Number of parallel workers for chromosome-level
                 processing (default: 4, requires Linux/macOS)
  --yield   INT  Number of reads to stream per chunk (default: 50000,
                 reduce if memory is limited)

Output columns:
  read_name        Read identifier
  chromosome       Chromosome name
  feature          Feature type (exon)
  ref_start        Exon start in reference coordinates
  ref_end          Exon end in reference coordinates
  read_start       Exon start in read coordinates
  read_end         Exon end in read coordinates
  strand           Strand of the annotated exon (+/-)
  gene_id          Ensembl gene ID
  gene_type        Gene biotype (e.g. protein_coding, lncRNA)
  gene_name        Gene symbol
  transcript_id    Ensembl transcript ID
  transcript_type  Transcript biotype
  transcript_name  Transcript name
  exon_number      Exon number within transcript
  exon_id          Ensembl exon ID
  soft_clip_start  Soft-clipped bases at the start of the read
  soft_clip_end    Soft-clipped bases at the end of the read
  fastq_length     Total read length

Notes:
  - BAM file must be coordinate-sorted and indexed (samtools index)
  - Chromosome naming (chr-prefixed vs bare) is reconciled automatically
  - GTF columns gene_biotype/transcript_biotype are normalised to
    gene_type/transcript_type automatically
  - Parallel execution uses forking (MulticoreParam) and requires
    Linux or macOS. On Windows execution falls back to serial mode.
  - Unplaced/decoy contigs (chrUn_*, *_random, EBV) are skipped

Examples:
  # Basic usage
  Rscript RREA.R --bam sample.bam --gtf gencode.v44.gtf --out rrea_out.tsv

  # With 8 cores and smaller chunks for lower memory usage
  Rscript RREA.R --bam sample.bam --gtf gencode.v44.gtf --out rrea_out.tsv \\
    --workers 8 --yield 10000

Dependencies:
  Rsamtools, GenomicAlignments, rtracklayer, GenomicRanges,
  BiocParallel, data.table, optparse
"

cat(usage)

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(rtracklayer)
  library(GenomicRanges)
  library(BiocParallel)
  library(data.table)
  library(optparse)
  library(dplyr)
})

# ── 1. Argument parsing ───────────────────────────────────────────────────────
option_list <- list(
  make_option("--bam",     type = "character", help = "Input indexed BAM file"),
  make_option("--gtf",     type = "character", help = "Annotation GTF file"),
  make_option("--out",     type = "character", help = "Output tab-delimited file"),
  make_option("--workers", type = "integer",   default = 4L,
              help = "Number of parallel workers [default: 4]"),
  make_option("--yield",   type = "integer",   default = 50000L,
              help = "Reads per streaming chunk [default: 50000]")
)
parser <- OptionParser(
  usage       = usage,
  option_list = option_list
)

opt <- parse_args(parser)

# Debug: print all parsed arguments
cat("BAM:", opt$bam,  "\n")
cat("GTF:", opt$gtf,  "\n")
cat("OUT:", opt$out,  "\n")


# ── 2. Load GTF exons once ────────────────────────────────────────────────────
cat("Loading GTF …\n")


# ── Helper function 1: GTF column name normalisation ─────────────────────────────────────────────

# Some GTF sources (e.g. Ensembl) use 'gene_biotype' / 'transcript_biotype'
# while others (e.g. GENCODE) use 'gene_type' / 'transcript_type'.
# This normalises to 'gene_type' / 'transcript_type' for downstream consistency.
normalise_gtf_columns <- function(exons) {
  col_map <- c(
    gene_type       = "gene_biotype",
    transcript_type = "transcript_biotype"
  )
  current_cols <- colnames(mcols(exons))
  for (old_name in names(col_map)) {
    new_name <- col_map[[old_name]]
    if (old_name %in% current_cols && !new_name %in% current_cols) {
      colnames(mcols(exons))[current_cols == old_name] <- new_name
      current_cols <- colnames(mcols(exons))  # refresh after rename
      cat(sprintf("  Renamed GTF column '%s' -> '%s'\n", old_name, new_name))
    }
  }
  exons
}

# ── Helper function 2: Chromosome naming reconciliation ─────────────────────────────────────────
# Detects whether BAM and exon GRanges use 'chr'-prefixed or bare chromosome
# names and strips/adds 'chr' on the exons to match the BAM if needed.
reconcile_chr_names <- function(exons, bam_path) {
  bam_chrs  <- names(scanBamHeader(bam_path)[[1]]$targets)
  exon_chrs <- unique(as.character(seqnames(exons)))
  
  bam_has_chr  <- any(grepl("^chr", bam_chrs))
  exon_has_chr <- any(grepl("^chr", exon_chrs))
  
  if (bam_has_chr == exon_has_chr) {
    cat("Chromosome names already match — no reconciliation needed\n")
    return(exons)
  }
  
  if (bam_has_chr && !exon_has_chr) {
    cat("Adding 'chr' prefix to exon chromosome names to match BAM\n")
    new_levels <- paste0("chr", seqlevels(exons))
  } else {
    cat("Stripping 'chr' prefix from exon chromosome names to match BAM\n")
    new_levels <- sub("^chr", "", seqlevels(exons))
  }
  
  # Rebuild GRanges cleanly rather than modifying in place
  names(new_levels) <- seqlevels(exons)
  exons <- renameSeqlevels(exons, new_levels)
  exons <- dropSeqlevels(exons, 
                         setdiff(seqlevels(exons), 
                                 unique(as.character(seqnames(exons)))),
                         pruning.mode = "coarse")
  
  # Verify
  n_matched <- length(intersect(bam_chrs, unique(as.character(seqnames(exons)))))
  cat(sprintf("  %d chromosomes matched between BAM and GTF\n", n_matched))
  if (n_matched == 0)
    warning("No chromosomes matched after reconciliation — check your GTF and BAM")
  
  exons
}

# load GTF
gtf   <- rtracklayer::import(opt$gtf)
exons <- gtf[gtf$type == "exon"]
exons <- normalise_gtf_columns(exons)
exons <- reconcile_chr_names(exons, opt$bam)

required_cols <- c("gene_id", "gene_biotype", "gene_name",
                   "transcript_id", "transcript_biotype", "transcript_name",
                   "exon_number", "exon_id")
missing <- setdiff(required_cols, colnames(mcols(exons)))
if (length(missing))
  stop("GTF is missing expected columns: ", paste(missing, collapse = ", "))

# ── 3. CIGAR parsing → data.table of aligned pairs ───────────────────────────
# Returns a data.table(read_idx, ref_pos, read_pos) for a batch of reads.
# Only emits rows where both ref_pos and read_pos are defined (matches only).
parse_cigars_dt <- function(cigars, ref_starts, read_indices) {
  ops_list  <- GenomicAlignments::explodeCigarOps(cigars)
  lens_list <- GenomicAlignments::explodeCigarOpLengths(cigars)
  
  results <- vector("list", length(cigars))
  
  for (i in seq_along(cigars)) {
    ops  <- ops_list[[i]]
    lens <- lens_list[[i]]
    r    <- 1L
    g    <- ref_starts[i]
    
    ref_out  <- vector("list", length(ops))
    read_out <- vector("list", length(ops))
    
    for (j in seq_along(ops)) {
      op <- ops[j]; n <- lens[j]
      if (op %in% c("M", "=", "X")) {
        ref_out[[j]]  <- g:(g + n - 1L)
        read_out[[j]] <- r:(r + n - 1L)
        r <- r + n; g <- g + n
      } else if (op == "I") {
        r <- r + n          # consumes read only
      } else if (op %in% c("D", "N")) {
        g <- g + n          # consumes ref only
      } else if (op == "S") {
        r <- r + n          # soft clip: advance read cursor only
      }
      # H / P: neither cursor advances
    }
    
    results[[i]] <- data.table(
      read_idx = read_indices[i],
      ref_pos  = unlist(ref_out,  use.names = FALSE),
      read_pos = unlist(read_out, use.names = FALSE)
    )
  }
  
  rbindlist(results)
}

# Soft-clip lengths for a vector of cigar strings
parse_soft_clips <- function(cigars, read_indices, qwidths, rnames) {
  ops_list  <- GenomicAlignments::explodeCigarOps(cigars)
  lens_list <- GenomicAlignments::explodeCigarOpLengths(cigars)
  sc_start  <- integer(length(cigars))
  sc_end    <- integer(length(cigars))
  for (i in seq_along(cigars)) {
    ops  <- ops_list[[i]]
    lens <- lens_list[[i]]
    if (length(ops) >= 1 && ops[1]           == "S") sc_start[i] <- lens[1]
    if (length(ops) >= 1 && ops[length(ops)] == "S") sc_end[i]   <- lens[length(lens)]
  }
  data.table(
    read_idx     = read_indices,
    read_name    = rnames,
    fastq_length = qwidths,
    sc_start     = sc_start,
    sc_end       = sc_end
  )
}

# ── 4. Position lookup via data.table rolling joins ───────────────────────────
# pairs_dt : data.table(read_idx, ref_pos, read_pos)
# hits_dt  : data.table(exon_idx, read_idx, exon_start, exon_end)
# Adds read_start and read_end columns to hits_dt in-place
map_positions <- function(pairs_dt, hits_dt) {
  setkeyv(pairs_dt, c("read_idx", "ref_pos"))
  
  # Start: roll forward — nearest ref_pos >= exon_start
  start_q <- data.table(read_idx = hits_dt$read_idx, ref_pos = hits_dt$exon_start)
  start_m <- pairs_dt[start_q, roll = -Inf, on = c("read_idx", "ref_pos")]
  hits_dt[, read_start := start_m$read_pos]
  
  # End: roll backward — nearest ref_pos <= exon_end
  end_q <- data.table(read_idx = hits_dt$read_idx, ref_pos = hits_dt$exon_end)
  end_m <- pairs_dt[end_q, roll = Inf, on = c("read_idx", "ref_pos")]
  hits_dt[, read_end := end_m$read_pos]
  
  hits_dt
}

# ── 5. Per-chromosome worker ──────────────────────────────────────────────────
process_chromosome <- function(chr, bam_path, exons, yield_size) {
  
  chr_len <- tryCatch(
    seqlengths(BamFile(bam_path))[[chr]],
    error = function(e) 536870912L
  )
  param <- ScanBamParam(
    which = GRanges(chr, IRanges(1L, chr_len)),
    what  = c("qname", "cigar", "qwidth")
  )
  
  chr_exons <- exons[as.character(seqnames(exons)) == chr]
  if (length(chr_exons) == 0) return(NULL)
  
  # Exon metadata as a keyed data.table for fast joining
  exon_meta <- data.table(
    exon_idx        = seq_along(chr_exons),
    chromosome      = as.character(seqnames(chr_exons)),
    feature         = as.character(chr_exons$type),
    exon_start      = start(chr_exons),
    exon_end        = end(chr_exons),
    strand          = as.character(strand(chr_exons)),
    gene_id         = chr_exons$gene_id,
    gene_biotype       = chr_exons$gene_biotype,
    gene_name       = chr_exons$gene_name,
    transcript_id   = chr_exons$transcript_id,
    transcript_biotype = chr_exons$transcript_biotype,
    transcript_name = chr_exons$transcript_name,
    exon_number     = as.character(chr_exons$exon_number),
    exon_id         = chr_exons$exon_id
  )
  setkey(exon_meta, exon_idx)
  
  bam_file <- BamFile(bam_path, yieldSize = yield_size)
  open(bam_file)
  on.exit(close(bam_file), add = TRUE)
  
  output_chunks <- list()
  chunk_i       <- 0L
  read_count    <- 0L
  
  repeat {
    chunk <- readGAlignments(bam_file, use.names = TRUE, param = param)
    if (length(chunk) == 0) break
    chunk_i <- chunk_i + 1L
    
    cat(sprintf("  [%s] chunk %d: %d reads\n", chr, chunk_i , length(chunk)))
    
    # ── 5a. Vectorised overlap across all reads in chunk ─────────────────────
    chunk_blocks_list <- grglist(chunk)
    chunk_blocks      <- unlist(chunk_blocks_list)
    block_to_read     <- rep(seq_along(chunk),
                             times = elementNROWS(chunk_blocks_list))
    
    hits <- findOverlaps(chr_exons, chunk_blocks, ignore.strand = TRUE)
    if (length(hits) == 0) next
    
    # Deduplicate to one row per unique (read, exon) pair
    hits_dt <- unique(data.table(
      exon_idx = queryHits(hits),
      read_idx = block_to_read[subjectHits(hits)]
    ))
    
    # ── 5b. Parse CIGARs only for reads that overlapped ──────────────────────
    unique_reads <- unique(hits_dt$read_idx)
    read_count   <- read_count + length(unique_reads)
    
    cigars     <- cigar(chunk)[unique_reads]
    ref_starts <- start(chunk)[unique_reads]
    qwidths    <- qwidth(chunk)[unique_reads]
    rnames     <- names(chunk)[unique_reads]
    
    pairs_dt <- parse_cigars_dt(cigars, ref_starts, unique_reads)
    
    # ── 5c. Soft clips ────────────────────────────────────────────────────────
    sc_dt <- parse_soft_clips(cigars, unique_reads, qwidths, rnames)
    
    # ── 5d. Map exon boundaries → read positions (rolling joins) ─────────────
    hits_dt <- merge(hits_dt,
                     exon_meta[, .(exon_idx, exon_start, exon_end)],
                     by = "exon_idx")
    hits_dt <- map_positions(pairs_dt, hits_dt)
    
    # ── 5e. Attach exon metadata and read metadata ────────────────────────────
    result <- merge(hits_dt, exon_meta, by = c("exon_idx", "exon_start", "exon_end"))
    result <- merge(result,  sc_dt,     by = "read_idx")
    
    # Clean up and reorder columns to match original output schema
    result[, c("read_idx", "exon_idx") := NULL]
    setcolorder(result, c(
      "read_name", "chromosome", "feature",
      "exon_start", "exon_end", "read_start", "read_end",
      "strand", "gene_id", "gene_biotype", "gene_name",
      "transcript_id", "transcript_biotype", "transcript_name",
      "exon_number", "exon_id", "sc_start", "sc_end", "fastq_length"
    ))
    setnames(result,
             c("exon_start", "exon_end", "sc_start", "sc_end"),
             c("ref_start",  "ref_end",  "soft_clip_start", "soft_clip_end"))
    
    output_chunks[[chunk_i]] <- result
  }
  
  cat(sprintf("  [%s] %d reads with exon overlaps\n", chr, read_count))
  if (length(output_chunks) == 0) return(NULL)
  rbindlist(output_chunks)
}

# ── 6. Chromosome list ────────────────────────────────────────────────────────

cat("Scanning BAM header …\n")
bam_header <- scanBamHeader(opt$bam)
all_chrs   <- names(bam_header[[1]]$targets)
gtf_chrs   <- unique(as.character(seqnames(exons)))

# Intersect BAM and GTF chromosomes; skip unplaced/decoy contigs
chrs <- intersect(all_chrs, gtf_chrs)
chrs <- chrs[grepl("^(chr)?([0-9]+|[XYM]|MT)$", chrs)]

cat(sprintf("Processing %d chromosomes with %d workers …\n",
            length(chrs), opt$workers))

# ── 7. Parallel dispatch ──────────────────────────────────────────────────────
# MulticoreParam inherits the full session — no export needed
# Only works on Linux/Mac, not Windows
if (.Platform$OS.type != "windows" && opt$workers > 1L) {
  register(MulticoreParam(workers = opt$workers, progressbar = TRUE))
} else {
  register(SnowParam(
    workers      = opt$workers,
    type         = "SOCK",
    exportglobals = TRUE
  ))
}

results <- bplapply(
  chrs,
  process_chromosome,
  bam_path   = opt$bam,
  exons      = exons,
  yield_size = opt$yield
)

# ── 8. Write output ───────────────────────────────────────────────────────────
cat("Writing output …\n")
output_dt <- rbindlist(Filter(Negate(is.null), results))
fwrite(output_dt, file = opt$out, sep = "\t")

cat(sprintf("Done. %d rows written to %s\n", nrow(output_dt), opt$out))

fwrite(output_dt, file = "your_output.tsv", sep = "\t")
