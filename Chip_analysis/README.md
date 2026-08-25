# ChIP-seq / CUT&Tag Coverage Analysis

Downstream R analysis relating ChIP-seq and CUT&Tag bigWig/peak signal to
gene structure (exons, introns, intergenic regions) and expression. This
folder does **not** perform alignment or peak calling — it consumes
pre-generated bigWigs and narrowPeak files.

## Data sources

- **ChIP-seq**: public ENCODE data (H3K27ac, H3K9me3, POL2, ± POLR2A-KO).
  See Table 5 for accessions.
- **CUT&Tag**: generated in-house (H3K4me3, POL2). See `../CutTag/` for the
  processing pipeline that produced these bigWigs.

## Scripts

### `ChIP_seq_compute_coverage.R` / `CT_compute_coverage.R`
Both compute normalized exon/intron gene-body coverage profiles (TSS→TES,
20 bins) over differentially expressed genes. Originally one script
(`MAIN_chip_bw_compute_coverage.R`) that toggled between ChIP-seq and
CUT&Tag data via which of two data-source blocks was left active — split
into two single-purpose scripts for clarity. Same pipeline structure,
deliberately duplicated rather than factored into shared code.

- `CT_compute_coverage.R` → **Fig. 3D**; **SI Fig. 3G** (full 4-row expression facet)
- `ChIP_seq_compute_coverage.R` → exploratory code; no figure match yet

### `intergenic_coverage.R`
Builds intergenic intervals from gene-to-gene gaps and relates ChIP-seq
signal/distance in those intervals to upstream gene expression tertile.

- `p_density` → **Fig. 7G**

### `peak_scatter.R`
Relates narrowPeak signal (H3K27ac, H3K9me3, POL2 ± POLR2A-KO) to intron
length, stratified by DEG status and expression.

- **Fig. 7F** (`"Cumulative fraction of total signal"` /
  `"ChIP mark"` — exact match)

## Color convention

All scripts share a Tol "Muted" palette (rose/sand/indigo/cyan/teal) for
consistency across ChIP and CUT&Tag panels — marks common to both data
sources (e.g., H3K27ac/H3K4me3, POL2) intentionally share colors:

```r
sample_colors <- c(
  H3K27ac    = "#CC6677",  # rose
  H3K4me3    = "#CC6677",  # rose (same as H3K27ac, by design)
  H3K27ac_ko = "#DDCC77",  # sand
  H3K9me3    = "#332288",  # indigo
  H3K9me3_ko = "#88CCEE",  # cyan
  POL2       = "#44AA99"   # teal
)
```

Some scripts contain additional, earlier color-palette definitions that
were never used in a final figure (exploratory only) — the palette above is
the one that matches published panels.
