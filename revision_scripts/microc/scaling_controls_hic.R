##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Hi-C fork of the contact-scaling control analysis (genes >= 100 kb).
# Downsampled Hi-C: WT vs POLR2A degron (6 h auxin). Controls are
# length-matched non-expressed genes and non-genic windows, both drawn only
# from unexpressed neighborhoods; shared functions are in
# scaling_controls_main.R (same directory).
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
# Hi-C settings
# ══════════════════════════════════════════════════════════════════════════════
BASE_NEW <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
HIC_DIR  <- file.path(BASE_NEW, "mega_ds_wt_pol2")

ASSAY <- "Hi-C"
CONDITIONS <- list(
  wt   = list(display="Wildtype",   color="#88CCEE",
              hic=file.path(HIC_DIR, "WT_HCT116_CTRL/mega/aligned/inter_30.hic")),
  pol2 = list(display="POLR2A deg", color="#CC6677",
              hic=file.path(HIC_DIR, "Pol2_6hrs_Aux/mega/aligned/inter_30.hic"))
)
ANALYSES    <- "long"      # genes >= 100 kb, res 5 kb, fit 5 kb - 1 Mb
COMPARTMENT <- "all"
RDS_BASE    <- "/projects/b1042/BackmanLab/Lucas/scaling_controls_rds/hic"

# every plot is saved here as an SVG (folder must already exist)
PLOT_DIR    <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/micro-c/plots/scaling_hic_plots"

results_hic <- run_scaling_analysis()