#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Calculates Hi-C resolution (Rao 2014 definition) for
# all mega .hic files. Resolution = finest bin size where >=80% of
# genome-wide intra-chromosomal bins have >=1000 contacts (Rao et
# al. 2014, Cell).
#
# Uses strawr R package to stream contact data directly from .hic
# files. Runs inside chiapet-hic-env conda environment.
#
# Output: /projects/b1042/BackmanLab/Lucas/060825_hic-samp/qc_plots/hic_resolution.tsv
# Usage:  sbatch calc_resolution.sh
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --mem=150G
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH -o /projects/b1042/BackmanLab/Lucas/060825_hic-samp/logs/hic_resolution-%j.out
#SBATCH -e /projects/b1042/BackmanLab/Lucas/060825_hic-samp/logs/hic_resolution-%j.err
#SBATCH -J hic_resolution

source activate chiapet-hic-env

BASE=/projects/b1042/BackmanLab/Lucas/060825_hic-samp
mkdir -p ${BASE}/logs

Rscript - << 'REOF'
suppressPackageStartupMessages({
  library(strawr)
  library(dplyr)
  library(readr)
})

BASE    <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
OUTDIR  <- file.path(BASE, "qc_plots")
OUTFILE <- file.path(OUTDIR, "hic_resolution.tsv")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ── .hic files to evaluate ─────────────────────────────────────────────────
hic_files <- list(
  "ds_wt_pol2 / WT"       = file.path(BASE, "mega_ds_wt_pol2/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
  "ds_wt_pol2 / Pol2 6hr" = file.path(BASE, "mega_ds_wt_pol2/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
  "ds_all / WT"           = file.path(BASE, "mega_ds_all/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
  "ds_all / Pol2 6hr"     = file.path(BASE, "mega_ds_all/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
  "ds_all / ActD 1hr"     = file.path(BASE, "mega_ds_all/1hr_ActD/mega/aligned/inter_30.hic"),
  "unsampled / WT"        = file.path(BASE, "mega_unsampled/WT_HCT116_CTRL/mega/aligned/inter_30.hic"),
  "unsampled / Pol2 6hr"  = file.path(BASE, "mega_unsampled/Pol2_6hrs_Aux/mega/aligned/inter_30.hic"),
  "unsampled / ActD 1hr"  = file.path(BASE, "mega_unsampled/1hr_ActD/mega/aligned/inter_30.hic")
)

# Chromosomes (bare number format matching .hic files)
CHROMS      <- c(as.character(1:22), "X")
MIN_CONTACTS <- 1000L
MIN_COVERAGE <- 0.80

# ── Resolution calculation ──────────────────────────────────────────────────
calc_resolution <- function(hic_path, label) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Processing:", label, "\n")

  if (!file.exists(hic_path)) {
    cat("  SKIP: file not found\n")
    return(NULL)
  }

  chroms_in_file <- tryCatch(readHicChroms(hic_path),   error = function(e) NULL)
  resolutions    <- tryCatch(readHicBpResolutions(hic_path), error = function(e) NULL)

  if (is.null(chroms_in_file) || is.null(resolutions)) {
    cat("  ERROR: could not read .hic metadata\n")
    return(NULL)
  }

  # Keep only standard chromosomes present in the file
  valid_chroms <- chroms_in_file[chroms_in_file$name %in% CHROMS, ]
  cat("  Chromosomes found:", paste(valid_chroms$name, collapse = ", "), "\n")
  cat("  Resolutions:", paste(resolutions, collapse = ", "), "\n")

  results     <- list()
  called_res  <- NA_integer_

  for (res in sort(resolutions)) {
    total_bins   <- 0L
    covered_bins <- 0L

    for (i in seq_len(nrow(valid_chroms))) {
      chrom     <- valid_chroms$name[i]
      chrom_len <- valid_chroms$length[i]
      n_bins    <- ceiling(chrom_len / res)
      total_bins <- total_bins + n_bins

      contacts <- tryCatch(
        straw("NONE", hic_path, chrom, chrom, "BP", res),
        error = function(e) NULL
      )

      if (!is.null(contacts) && nrow(contacts) > 0) {
        # Sum contacts per bin (both x and y positions)
        bin_totals <- contacts %>%
          tidyr::pivot_longer(c(x, y), values_to = "bin") %>%
          group_by(bin) %>%
          summarise(total = sum(counts), .groups = "drop")
        covered_bins <- covered_bins + sum(bin_totals$total >= MIN_CONTACTS)
      }
    }

    pct <- if (total_bins > 0) covered_bins / total_bins * 100 else 0

    cat(sprintf("  %10s bp  covered=%7d  total=%7d  (%.1f%%)\n",
                format(res, big.mark = ","), covered_bins, total_bins, pct))

    results[[length(results) + 1]] <- list(
      resolution_bp   = res,
      bins_covered    = covered_bins,
      bins_total      = total_bins,
      pct_covered     = round(pct, 2)
    )

    if (pct >= MIN_COVERAGE * 100 && is.na(called_res)) {
      called_res <- res
    }
  }

  if (!is.na(called_res)) {
    cat(sprintf("  ✓ Called resolution: %s bp\n",
                format(called_res, big.mark = ",")))
  } else {
    cat("  ✗ Did not reach 80% coverage at any resolution\n")
  }

  parts <- strsplit(label, " / ")[[1]]
  bind_rows(lapply(results, function(r) {
    tibble(
      dataset          = trimws(parts[1]),
      condition        = trimws(parts[2]),
      resolution_bp    = r$resolution_bp,
      bins_covered     = r$bins_covered,
      bins_total       = r$bins_total,
      pct_covered      = r$pct_covered,
      called_resolution = ifelse(is.na(called_res), "not reached",
                                 as.character(called_res))
    )
  }))
}

# ── Run all files ───────────────────────────────────────────────────────────
all_results <- bind_rows(mapply(
  calc_resolution,
  hic_path = unlist(hic_files),
  label    = names(hic_files),
  SIMPLIFY = FALSE
))

# ── Write output ────────────────────────────────────────────────────────────
if (nrow(all_results) > 0) {
  write_tsv(all_results, OUTFILE)
  cat("\nResults written to:", OUTFILE, "\n")
} else {
  cat("\nNo results to write.\n")
}

# ── Summary ─────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n", sep = "")
cat("CALLED RESOLUTIONS SUMMARY\n")
cat(strrep("=", 60), "\n", sep = "")
all_results %>%
  distinct(dataset, condition, called_resolution) %>%
  arrange(dataset, condition) %>%
  mutate(label = sprintf("  %-20s %-15s → %s bp",
                         dataset, condition, called_resolution)) %>%
  pull(label) %>%
  cat(sep = "\n")
cat("\n")
REOF
