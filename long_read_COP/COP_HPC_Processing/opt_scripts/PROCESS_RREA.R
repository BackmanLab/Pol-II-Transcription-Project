#!/usr/bin/env Rscript
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Takes the output of RREA.R and computes per-read intron lengths
# and exon/intron coverage fractions.
#
# Usage:
#   Rscript PROCESS_RREA.R <input_file> <output_prefix>
#
# Outputs:
#   <prefix>_introns.rds           -- per-exon rows with intron gap annotations
#   <prefix>_introns_collapsed.rds -- one row per read (largest intron)
#   <prefix>_fractions.rds         -- per-read exon/intron coverage fractions
#
# Explanation of output:
#   read_exon_fraction = read_exon_length / read_length
#     Fraction of the aligned read that maps to exonic sequence. The
#     denominator is the read length minus soft clips -- how many
#     bases of the read actually aligned to the genome. The
#     numerator is how many of those aligned bases fell within
#     annotated exon boundaries.
#   ref_exon_fraction = ref_exon_length / ref_length
#     Fraction of the genomic span covered by the read that is
#     annotated as exonic. The denominator is the total reference
#     span from the first to last exon coordinate in the read
#     (max(ref_end) - min(ref_start)). The numerator is how many
#     bases within that span are covered by annotated exons.
#   read_intron_fraction = read_intron_length / read_length
#     Fraction of the aligned read that maps to intronic sequence --
#     bases in the read that fall in the gaps between annotated
#     exons. This is the most direct measure of intron retention. A
#     high read_intron_fraction means the read contains a lot of
#     sequence that should have been spliced out.
#   ref_intron_fraction = ref_intron_length / ref_length
#     Fraction of the genomic span covered by the read that is
#     annotated as intronic. This is essentially a property of the
#     reference annotation across the region the read spans, not a
#     property of the read itself -- it tells you how intron-dense
#     the genomic locus is.
#
# Figure output: none directly -- upstream data-prep for
# long_read_COP/RREA_analysis.R (SI Fig. 13G-H).
##################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(IRanges)
  library(optparse)
})

# ── Usage ─────────────────────────────────────────────────────────────────────
usage <- "
process_RREA.R — Intron Retention Metrics from RREA Output
===========================================================
Computes per-read intron gap lengths and exon/intron coverage fractions
from the tab-delimited output of RREA.R.

Usage:
  Rscript process_RREA.R --input <file> --output <prefix>

Required arguments:
  --input   STR   Tab-delimited output file from RREA.R
  --output  STR   Prefix for output .rds files

Output files:
  <prefix>_introns.rds
    Full RREA table with intron gap and transcript span columns appended.
    One row per read-exon overlap.

  <prefix>_introns_collapsed.rds
    One row per read, retaining the exon preceding the largest intron.

  <prefix>_fractions.rds
    One row per read with exon/intron length and fraction summaries.
    Key columns:
      read_intron_fraction  fraction of aligned read that is intronic
      ref_intron_fraction   fraction of reference span that is intronic
      read_exon_fraction    fraction of aligned read that is exonic
      ref_exon_fraction     fraction of reference span that is exonic
      longest_intron        length of the largest intron overlapped by the read

Notes:
  - read_intron_fraction is the most direct measure of intron retention
  - Single-base exon overlaps (read_start == read_end) can be filtered
    out if uninformative: dt <- dt[read_start != read_end]

Examples:
  Rscript process_RREA.R --input rrea_output.tsv --output sample_processed

Dependencies:
  data.table, IRanges, optparse
"

cat(usage)

# ── 1. Arguments ──────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input",  type = "character", help = "RREA output file (.tsv)"),
  make_option("--output", type = "character", help = "Output file prefix")
)

parser <- OptionParser(
  description = "Compute intron lengths and coverage fractions from RREA output",
  usage       = usage,
  option_list = option_list
)

opt <- parse_args(parser)


if (is.null(opt$input) || is.null(opt$output))
  stop("Both --input and --output are required. ",
       "Usage: Rscript process_rrea.R --input <file> --output <prefix>")

# ── 2. Load ───────────────────────────────────────────────────────────────────
cat("Loading RREA output …\n")
rrea <- fread(opt$input)

# ── 3. Intron gap annotation ──────────────────────────────────────────────────
# For each (read, transcript) pair, compute the gap between consecutive exons
# in both read-space and reference-space, and derive transcript-level metrics.
calc_introns <- function(dt) {
  setorder(dt, read_name, transcript_id, read_start)
  
  dt[, `:=`(
    next_read_start = shift(read_start, type = "lead"),
    next_ref_start  = shift(ref_start,  type = "lead")
  ), by = .(read_name, transcript_id)]
  
  dt[, `:=`(
    next_read_intron = as.numeric(next_read_start - read_end),
    next_ref_intron  = as.numeric(next_ref_start  - ref_end)
  )]
  
  dt[, delta_intron := abs(next_ref_intron - next_read_intron)]
  
  dt[, `:=`(
    read_tx_start = min(read_start),
    read_tx_end   = max(read_end),
    ref_tx_start  = min(ref_start),
    ref_tx_end    = max(ref_end)
  ), by = .(read_name, transcript_id)]
  
  dt[, `:=`(
    read_tx_length = as.numeric(read_tx_end - read_tx_start),
    ref_tx_length  = as.numeric(ref_tx_end  - ref_tx_start)
  )]
  
  dt[, delta_length := abs(ref_tx_length - read_tx_length)]
  
  dt[]
}

# ── 4. Coverage fractions ─────────────────────────────────────────────────────
# Summarise each read into total exon/intron lengths and their fractions of the
# aligned read length. IRanges::reduce() merges overlapping exon blocks before
# summing, avoiding double-counting.
calc_fractions <- function(dt) {
  dt[, aligned_length := fastq_length - soft_clip_start - soft_clip_end]
  
  # ── Pre-compute reduced read-space ranges across all reads at once ──────────
  read_ir  <- IRanges(dt$read_start, dt$read_end)
  ref_ir   <- IRanges(dt$ref_start,  dt$ref_end)
  read_grp <- dt$read_name
  
  # Split into per-read IRangesList, reduce, then compute widths in one pass
  read_reduced <- reduce(split(read_ir, read_grp))
  ref_reduced  <- reduce(split(ref_ir,  read_grp))
  
  read_exon_lengths   <- sum(width(read_reduced))
  ref_exon_lengths    <- sum(width(ref_reduced))
  read_intron_lengths <- sum(width(gaps(read_reduced)))
  ref_intron_lengths  <- sum(width(gaps(ref_reduced)))
  
  # Collect into a data.table keyed on read_name
  ranges_dt <- data.table(
    read_name          = names(read_exon_lengths),
    read_exon_length   = as.numeric(read_exon_lengths),
    ref_exon_length    = as.numeric(ref_exon_lengths),
    read_intron_length = as.numeric(read_intron_lengths),
    ref_intron_length  = as.numeric(ref_intron_lengths)
  )
  
  # ── Per-read scalar summaries via data.table ────────────────────────────────
  scalars_dt <- dt[, .(
    read_length    = aligned_length[1L],
    ref_length     = as.numeric(max(ref_end) - min(ref_start)),
    longest_intron = if (all(is.na(next_ref_intron))) NA_real_
    else max(next_ref_intron, na.rm = TRUE)
  ), by = read_name]
  
  # ── Join and compute fractions ───────────────────────────────────────────────
  result <- ranges_dt[scalars_dt, on = "read_name"]
  result[, `:=`(
    read_exon_fraction   = read_exon_length   / read_length,
    ref_exon_fraction    = ref_exon_length    / ref_length,
    read_intron_fraction = read_intron_length / read_length,
    ref_intron_fraction  = ref_intron_length  / ref_length
  )][]
}

# ── 5. Run ────────────────────────────────────────────────────────────────────
cat("Annotating intron gaps …\n")
introns <- calc_introns(rrea)

cat("Computing coverage fractions …\n")
fractions <- calc_fractions(introns)

## explanation
## ref_intron_fraction = ref_intron_length / ref_length
## This is the fraction of the **genomic span** covered by the read that is annotated as intronic. It's essentially a property of the reference annotation across the region the read spans, not a property of the read itself. It tells you how intron-dense the genomic locus is.


# One row per read: the exon preceding the largest reference intron
cat("Collapsing to largest intron per read …\n")
collapsed <- introns[
  !is.na(next_ref_intron),
][
  introns[!is.na(next_ref_intron), .I[which.max(next_ref_intron)], by = read_name]$V1
]

# ── 6. Write outputs ──────────────────────────────────────────────────────────
cat("Writing outputs …\n")
saveRDS(introns,   file = paste0(opt$output, "_introns.rds"))
saveRDS(collapsed, file = paste0(opt$output, "_introns_collapsed.rds"))
saveRDS(fractions, file = paste0(opt$output, "_fractions.rds"))

cat(sprintf(
  "Done.\n  Introns:   %d rows\n  Collapsed: %d rows\n  Fractions: %d rows\n",
  nrow(introns), nrow(collapsed), nrow(fractions)
))