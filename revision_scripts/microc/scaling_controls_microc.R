##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Micro-C fork of the contact-scaling control analysis (ENCODE HCT116
# POLR2A-degron Micro-C). Runs, in order:
#   long  : genes >= 100 kb, res 5 kb, fit 5 kb - 1 Mb
#   short : genes 1 kb - 100 kb, res 500 bp, fit 1 kb - 100 kb
# Controls are length-matched non-expressed genes and non-genic windows,
# both drawn only from unexpressed neighborhoods; shared functions are in
# scaling_controls_main.R (same directory). Use >= 24 cores for the 500 bp
# extraction (it prints a time projection after the first tiles).
##################################################################
# ── Load the shared functions from this script's own directory ────────────────
.find_script_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {                       # source("this_file.R")
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(dirname(normalizePath(f)))
  }
  a <- grep("^--file=", commandArgs(FALSE), value=TRUE)          # Rscript this_file.R
  if (length(a) > 0) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error=function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))             # RStudio / Quest OnDemand
  }
  getwd()
}
.main_file <- file.path(.find_script_dir(), "scaling_controls_main.R")
if (!file.exists(.main_file))
  stop("scaling_controls_main.R not found next to this script (looked in ",
       dirname(.main_file), "). Keep the four scaling_controls_*.R files in one directory.")
source(.main_file)
message("Loaded shared functions from ", .main_file)

# ══════════════════════════════════════════════════════════════════════════════
# Micro-C settings
# ══════════════════════════════════════════════════════════════════════════════
MICROC_DIR <- "/projects/b1042/BackmanLab/Lucas/microc_encode"

ASSAY <- "Micro-C"
CONDITIONS <- list(
  wt   = list(display="Wildtype",   color="#88CCEE",     # ENCSR533VUS, auxin untreated
              hic=file.path(MICROC_DIR, "ENCFF534TYD.hic")),
  pol2 = list(display="POLR2A deg", color="#CC6677",     # ENCSR218URU, 6 h 5-Ph-IAA
              hic=file.path(MICROC_DIR, "ENCFF302CUI.hic"))
)
ANALYSES    <- c("long", "short")
COMPARTMENT <- "all"
RDS_BASE    <- "/projects/b1042/BackmanLab/Lucas/scaling_controls_rds/microc"

# every plot is saved here as an SVG (folder must already exist)
PLOT_DIR    <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/micro-c/plots/scaling_microc_plots"

results_microc <- run_scaling_analysis()