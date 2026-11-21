##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Micro-C A/B compartment fork: the Micro-C analysis (long and short) run
# separately inside A and inside B compartments. Compartments come from the
# WT Hi-C eigenvector (500 kb bedgraph); expressed genes are restricted to the
# compartment. The controls are genome-wide and identical in the A, B and
# all-compartment runs. Shared functions are in scaling_controls_main.R (same
# directory).
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
# Micro-C A/B settings
# ══════════════════════════════════════════════════════════════════════════════
MICROC_DIR <- "/projects/b1042/BackmanLab/Lucas/microc_encode"
BASE_NEW   <- "/projects/b1042/BackmanLab/Lucas/060825_hic-samp"
COMP_BED   <- file.path(BASE_NEW, "compartments/ds_wt_pol2_WT/ds_wt_pol2_WT_compartments.bedgraph")

ASSAY <- "Micro-C"
CONDITIONS <- list(
  wt   = list(display="Wildtype",   color="#88CCEE",
              hic=file.path(MICROC_DIR, "ENCFF534TYD.hic")),
  pol2 = list(display="POLR2A deg", color="#CC6677",
              hic=file.path(MICROC_DIR, "ENCFF302CUI.hic"))
)
ANALYSES <- c("long", "short")
RDS_BASE <- "/projects/b1042/BackmanLab/Lucas/scaling_controls_rds/microc"

# every plot is saved here as an SVG (folder must already exist)
PLOT_DIR <- "/home/lmc0633/Transcription Paper/101426_revision_scripts/micro-c/plots/scaling_microc_AB_plots"

results_microc_AB <- list()
for (COMPARTMENT in c("A", "B")) {
  message("\n████████ Compartment ", COMPARTMENT, " ████████")
  results_microc_AB[[COMPARTMENT]] <- run_scaling_analysis()
}