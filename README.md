# Gene transcription and chromatin packing domains form a self-organizing system

Code for Carter, Almassalha, et al. — *"Gene transcription and chromatin
packing domains form a self-organizing system."* bioRxiv (2026).
doi: [10.64898/2026.03.15.711889](https://doi.org/10.64898/2026.03.15.711889)
· PMID: [41889979](https://pubmed.ncbi.nlm.nih.gov/41889979/)

Backman Lab, Northwestern University.

## Overview

This repository contains the processing pipelines and downstream analysis
scripts for a study of RNA Polymerase II's role in chromatin packing domain
formation and maintenance, using an auxin-inducible POLR2A degron system in
HCT116 cells paired with Hi-C, RNA-seq/EU-seq/RIP-seq, ChIA-PET, ATAC-seq,
CUT&Tag, and multimodal imaging (ChromSTEM, STORM, PWS).

## Repository structure

| Folder | Contents |
|---|---|
| `HiC_analysis/` | Hi-C processing (Juicer), loop calling, compartments, TAD/insulation, contact scaling |
| `chia_loop_analysis/` | ChIA-PET locus plotting and loop anchor analysis |
| `RIPseq_analysis/` | RIP-seq (POLR2A, EZH2) processing and feature-level quantification |
| `EU_seq_analysis/` | Nascent transcription (EU-seq) processing and analysis |
| `total_RNA_seq/` | Total RNA-seq differential expression |
| `atac_analysis/` | ATAC-seq downstream analysis (public ENCODE data) |
| `Chip_analysis/` | ChIP-seq / CUT&Tag coverage analysis |
| `CutTag/` | CUT&Tag processing pipeline |
| `long_read_COP/` | Long-read nascent RNA-seq (RREA) splicing analysis |
| `imaging_analysis_scripts/` | PWS, EU-STORM, EdU/K9-STORM, ChromSTEM analysis |

Each folder has its own `README.md` with a script-by-script breakdown and,
where confirmed, the specific figure/panel each script produces. See
`FIGURES.md` for the full figure-to-script index across the whole paper.

## Requirements

See `ENVIRONMENT.md` for the full list of software versions and R/Python
packages used across this repository. Most processing pipelines were run on
Northwestern's Quest HPC cluster (SLURM); paths in submission scripts will
need to be updated for other environments.

## Public data

All publicly available datasets used in this study (ENCODE ChIA-PET/ATAC-seq/
ChIP-seq, GEO nanopore direct RNA-seq) are listed with accessions in Table 5
of the manuscript.

## Third-party code

`HiC_analysis/HiC_HPC_processing/juicer_src/` is excluded from this
repository — it is a modified copy of the
[Juicer 2.0 toolkit](https://github.com/aidenlab/juicer) (Aiden Lab, MIT
License). Install Juicer separately. The Quest-specific adaptation of the
Juicer pipeline for this study is documented in `HiC_analysis/README.md`.

## License

See `LICENSE`.

## Citation

See `CITATION.cff`, or cite:

> Carter LM, Almassalha LM, et al. Gene transcription and chromatin packing
> domains form a self-organizing system. bioRxiv. 2026.
> doi:10.64898/2026.03.15.711889
