# Hi-C Analysis

Hi-C for POLR2A-degron HCT116 cells (WT, POLR2A 6h auxin, 1h ActD), processed
with three downsampling groups (`ds_wt_pol2`, `ds_all`, `unsampled` — see
`HiC_HPC_processing/juicer/scripts/` for how each is built) and analyzed for
loops, compartments, TAD/insulation, and contact scaling.

## Structure

| Folder | Contents |
|---|---|
| `HiC_HPC_processing/` | Downsampling, Juicer alignment/loop-calling, mega merging |
| `Hic_loop_analysis/` | APA, differential loops (Gained/Stable/Lost), loop anchor classification, locus tracks |
| `compartment_analysis/` | A/B compartment eigenvectors, saddle plots, compartment strength |
| `tad_analysis/` | Aggregate TAD analysis (ATA), insulation score |
| `scaling_analysis/` | Contact scaling (S) by E/I ratio and expression, tandem-gene-cluster distance analysis |
| `HiC_stats_analysis/` | HiCRep stratum-adjusted correlation (SCC), general QC plots |

## A note on `HiC_HPC_processing/juicer_src/`

**This folder is excluded from the repository.** It is a modified copy of
the [Juicer 2.0 toolkit](https://github.com/aidenlab/juicer) (Aiden Lab, MIT
License). The `juicer_src` scripts were adapted for use on Northwestern's Quest HPC.

## Pipeline order (`HiC_HPC_processing/juicer/scripts/`)

Three downsampling groups run through this pipeline in parallel:
**`ds_wt_pol2`** (WT + Pol2 only, depth-matched to Pol2, higher depth),
**`ds_all`** (WT + Pol2 + ActD, depth-matched to ActD -- this is the group
that includes ActD), and **`unsampled`** (full depth, all conditions).

1. `00_count_reads.sh` -> `00b_make_samplesheet.sh` -- read counts and
   downsampling fractions per sample
2. `01_downsample_fastqs.sh` / `01_downsampled_wt_pol2.sh` -- depth-matched
   FASTQ downsampling (see script headers for which downsampling group each
   produces)
3. `02_run_juicer_ds_all.sh` / `02_run_juicer_ds_all_genhm.sh` (genomics-himem
   partition variant) / `02_run_juicer_ds_wt_pol2.sh` / `02_run_juicer_unsampled.sh`
   -- per-replicate Juicer alignment + loop calling, one script per
   downsampling group. Run from the login node -- `juicer.sh` submits its
   own SLURM jobs.
4. `03_merge_mega.sh`, and the per-group variants `run_mega_all.sh`
   (ds_all), `run_mega_ds_wt_pol2.sh` (ds_wt_pol2), `run_mega_unsampled.sh`
   (unsampled) -- merge per-replicate outputs into a mega `.hic` per
   condition. How `03_merge_mega.sh` relates to the three group-specific
   `run_mega_*.sh` scripts (generic dispatcher vs. an earlier/superseded
   version) isn't documented anywhere in the scripts themselves.
5. `run_hiccups.sh` / `run_hiccups_kr.sh` (KR-normalized variant, requires
   `juicer_tools addNorm -k KR` first) -- GPU-accelerated loop calling on
   merged and per-replicate `.hic` files
6. `run_compartments.sh` -- A/B compartment eigenvectors (Juicer eigenvector
   tool, SCALE-normalized, 500kb; sign-oriented by GC content -- see
   `compartment_analysis/` for the downstream saddle/strength analysis)
7. `calc_resolution.sh` -- Hi-C resolution per Rao et al. 2014's definition
   (finest bin size where >=80% of intra-chromosomal bins have >=1000 contacts)

`verify_downsampled.sh` and `test_hiccups.sh` are standalone diagnostics, not
part of the main sequence. `mega.sh`, `check.sh`, `cleanup.sh`,
`countligations.sh`, `relaunch_prep.sh`, `juicer_arrowhead.sh`,
`juicer_hiccups.sh`, `juicer_postprocessing.sh`, `juicer_genhm.sh`, and
`juicer.sh` itself are Juicer's own common scripts, present here because
Juicer's install expects them at this path -- not original code.

See each analysis subfolder for its own scripts and figure mapping.
