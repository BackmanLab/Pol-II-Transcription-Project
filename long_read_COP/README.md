# Long-Read Nascent RNA-seq (RREA)

Read-Relative Exon Annotations (RREA) analysis of public long-read nascent
RNA-seq data (Drexler et al., *Nat. Protoc.* 2021 -- K562 DMSO/PlaB, GM19102;
see Table 5). For each aligned read, computes per-exon overlap and intron
length in both read-space and reference-space, to assess splicing
completeness and exon/intron ratio.


## Pipeline (`COP_HPC_Processing/`)

1. `minimap.sh` -- bare single-sample minimap2 alignment (splice-aware, ONT
   preset). For the full alignment-to-primary-only-BAM pipeline (unique-read
   filtering + primary-alignment filtering), use `align_and_filter.sh`
   instead -- `opt_scripts/RREA.sh` expects `*_primary_only.bam` as input.
2. `align_and_filter.sh` -- minimap2 alignment -> unique-read filtering ->
   primary-only BAM (indexed). Usage: `sbatch align_and_filter.sh <NAME> <fastq> <reference.fa>`.
3. `opt_scripts/RREA.sh` -> `opt_scripts/RREA.R` -- per-chromosome, per-read
   exon overlap and intron gap extraction (data.table-optimized: chunked BAM
   streaming, vectorized `findOverlaps()`, keyed binary search). Kept
   alongside its own wrapper script, matching how they're actually
   referenced -- `RREA.sh`'s own `SCRIPT` variable expects `RREA.R` in this
   same directory.
4. `opt_scripts/PROCESS_RREA.sh` -> `opt_scripts/PROCESS_RREA.R` -- per-read
   intron lengths and exon/intron coverage fractions from `RREA.R` output.
   Same co-location rationale as above.
5. `RREA_analysis.R` (top level) -- splicing category classification
   (unspliced/partially retained/retained), per-read and per-gene E/I
   ratio, intron-length tertile stratification, termination position
   analysis. Includes a deduplication step for a known duplicate-row issue
   in `RREA.R`'s `process_chromosome()` (multiple transcript hits per exon
   from `findOverlaps()`) -- see the script's own header. Stays at the top
   level (not paired with a wrapper script) -- same convention as other
   downstream analysis scripts in this repo (e.g. `RIPseq_analysis/`,
   `atac_analysis/`).

Output: **SI Fig. 13G-H** (see `../FIGURES.md`).
