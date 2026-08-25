# CUT&Tag Data Processing

## Background

Cleavage Under Targets & Tagmentation (CUT&Tag) is a ChIP-like chromatin
profiling method in which primary antibodies are bound to chromatin proteins
in situ after nuclei are permeabilized. A secondary antibody conjugated to
Tn5 transposase then binds the primary antibody; activating the transposase
simultaneously cuts DNA and adds sequencing adapters ("tagmentation") for
paired-end sequencing. Compared to conventional ChIP, CUT&Tag libraries have
higher signal-to-noise and require far fewer input cells.

This pipeline follows the general approach described in the Henikoff Lab's
[CUT&Tag data processing tutorial](https://yezhengstat.github.io/CUTTag_tutorial/#I_Introduction)
and [protocol paper](https://www.nature.com/articles/s41596-020-0373-x)
([primary paper](https://www.nature.com/articles/s41467-019-09982-5)), run
as a standalone wrapper (`CutTag.sh` + `stats.R`) rather than through the
existing [nf-core/cutandrun](https://nf-co.re/cutandrun/3.2.2/) pipeline.

Further reading:
- [A 2023 review on CUT&Tag](https://www.tandfonline.com/doi/full/10.1080/15592294.2023.2293411#abstract)
- [Integrative analysis of CUT&Tag and RNA-seq data](https://link.springer.com/protocol/10.1007/978-1-0716-4071-5_13)
- [CUT&RUNTools 2.0](https://academic.oup.com/bioinformatics/article/38/1/252/6318389)
- [Spatial-CUT&Tag](https://www.science.org/doi/full/10.1126/science.abg7216)
- [Active Motif's spike-in strategy](https://www.activemotif.com/documents/AACR-2024-Spike-In-Poster.pdf)
- [GoPeaks](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02707-w)
- [CUT&Tag recovers up to half of ENCODE ChIP-seq peaks (H3K27 marks)](https://www.biorxiv.org/content/10.1101/2022.03.30.486382v2.full)
- [ChIPseqSpikeInFree](https://github.com/stjude/ChIPseqSpikeInFree)
- [SEACR](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4)
- [HBC ChIP-seq/CUT&RUN workshop](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/tree/main)

## Samples

Three biological replicate submissions. The second and third submissions
include technical replicates that should be merged before downstream analysis.

**Submission 1:**

| Cells | Antibody | Index    |
|-------|----------|----------|
| WT    | Pol1AB   | i71-i51  |
| WT    | K4Me3AB  | i71-i54  |
| WT    | Pol2AB   | i72-i51  |
| WT    | K9Me3AB  | i72-i53  |

**Submission 2:**

| Cells (Rep1) | Antibody | Index   |
|--------------|----------|---------|
| WT           | Pol1AB   | i71-i51 |
| WT           | Pol2AB   | i71-i54 |
| WT           | Pol1AB   | i73-i51 |
| WT           | Pol2AB   | i73-i54 |

**Submission 3:**

| Cells (Rep1) | Antibody | Index   |
|--------------|----------|---------|
| WT           | Pol1AB   | i71-i51 |
| WT           | K9AB     | i73-i51 |
| WT           | K9AB     | i73-i54 |

## Pipeline

### 1. Alignment

Alignment uses [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

Build a Bowtie2 index once per reference genome:
```shell
bowtie2-build path_to_reference_genome.fa prefix_to_name_indexes
```

For combined rDNA + genome mapping, use
[this rDNA-annotated hg38 build](https://github.com/vikramparalkar/rDNA-Mapping-Genomes/blob/main/Human_hg38-rDNA_genome_v1.0_annotation.tar.gz).
For general mapping, either Ensembl:

```shell
VERSION=108
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz
```

or UCSC (preferred; four GTF references available there):

```shell
wget -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
wget -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
```

Example SLURM submission for building the index on Quest:

```shell
#!/bin/bash
#SBATCH -A <allocation>
#SBATCH -p genhimem
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=10G
#SBATCH -t 24:00:00
#SBATCH --job-name=build_genome
#SBATCH --output=outlog
#SBATCH --error=errlog

module load bowtie2/2.4.5
bowtie2-build /path/to/genome/hg38.fa hg38
```

### 2. Running CutTag.sh

```shell
for file in *R1_001.fastq.gz; do
  bash CutTag.sh -f "${file}" -w <fastq dir> -g <bowtie2 index prefix> \
    -r <annotation GTF> -c <chrom sizes file> -t 16 -m 2 -d
done
```

See `CutTag.sh -h` (or the script's `helpFunction`) for the full flag list —
this covers QC, trimming, alignment, optional Picard duplicate removal,
sorting, BED/BEDGRAPH conversion, and peak calling.

### 3. Peak calling

Two peak callers are supported:
- **[MACS2](https://github.com/macs3-project/MACS)** — more lenient/traditional; higher sensitivity, higher false-positive rate.
- **[SEACR](https://github.com/FredHutch/SEACR/tree/master)** — purpose-built for high-signal-to-noise CUT&Run/Tag data.

([GoPeaks](https://github.com/maxsonBraunLab/gopeaks) is a newer CUT&Tag-specific
alternative not used here; see
[this comparison](https://help.pluto.bio/en/articles/choosing-between-macs2-and-seacr-for-peak-calling-with-dna-sequencing-datasets)
of MACS2 vs. SEACR.)

To use SEACR, fetch both of its scripts into a local directory first:
```shell
mkdir SEACR
wget https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.sh
wget https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.R
```
```shell
bash SEACR_1.3.sh <experimental bedgraph> <control bedgraph | numeric threshold> \
  ["norm" | "non"] ["relaxed" | "stringent"] <output prefix>
```

Both callers' invocations for this study are handled inside `CutTag.sh`.

### 4. Replicate correlation QC

Reproducibility between replicates was assessed by binned-fragment
correlation:

```R
library(corrplot)
library(dplyr)

hists <- list.files(bed_dir, pattern = ".bin.bed")
fragCount <- NULL
for (hist in hists) {
  tmp <- read.table(file.path(bed_dir, hist), header = FALSE)
  colnames(tmp) <- c("chrom", "bin", hist)
  fragCount <- if (is.null(fragCount)) tmp else full_join(fragCount, tmp, by = c("chrom", "bin"))
}

M <- cor(fragCount %>% select(-chrom, -bin) %>% log2(), use = "complete.obs")
corrplot(M, method = "color", outline = TRUE, addgrid.col = "darkgray",
         order = "hclust", addrect = 3, rect.col = "black", rect.lwd = 3,
         cl.pos = "b", tl.col = "indianred4", tl.cex = 0.5, cl.cex = 0.5,
         addCoef.col = "black", number.digits = 2, number.cex = 0.5,
         col = colorRampPalette(c("midnightblue", "white", "darkred"))(100))
```

### 5. Duplication rate and downsampling

Duplication rates were high across samples (consistent with PCR
overamplification combined with oversequencing by the sequencing core);
see [this preprint](https://www.biorxiv.org/content/10.1101/2022.03.30.486382v1.full.pdf)
for a discussion of similar duplication issues in CUT&Tag data.

Duplication was summarized from Picard's output:
```R
dupResult <- data.frame()
for (hist in list.files(sam_dir, pattern = ".dupMark.txt")) {
  dupRes <- read.table(file.path(sam_dir, hist), header = TRUE, fill = TRUE)
  info <- strsplit(hist, "_")[[1]]
  dupResult <- rbind(dupResult, data.frame(
    Histone = info[1], Replicate = info[2],
    MappedFragNum_hg38 = as.numeric(dupRes$READ_PAIRS_EXAMINED[1]),
    DuplicationRate = as.numeric(dupRes$PERCENT_DUPLICATION[1]) * 100,
    EstimatedLibrarySize = as.numeric(dupRes$ESTIMATED_LIBRARY_SIZE[1])
  ))
}
```

Read counts in deduplicated, indexed BAMs:
```shell
for file in *.ind.bam; do
  READS=$(samtools view -c "${file}")
  echo -e "${file}\t${READS}" >> reads.txt
done
```

Samples were then downsampled to a common read depth per condition (example
for one condition/read target):
```shell
module load samtools
outdir=<path/to/sampled/BAM>
seed=100
reads=3000000
array=(sample1.ind.bam sample2.ind.bam)

for bam in "${array[@]}"; do
  bamS="${bam/.ind.bam/.samp.ind.bam}"
  fraction=$(samtools idxstats "${bam}" | cut -f3 | \
    awk -v ct=${reads} 'BEGIN{t=0}{t+=$1}END{print ct/t}')
  samtools view -b -s ${seed}$(echo ${fraction} | sed 's/^0*//') "${bam}" > "${outdir}/${bamS}"
done
```

### 6. Merging downsampled replicates

Downsampled replicate BAMs were merged per condition to maximize signal for
peak calling:
```shell
module load samtools
samtools merge merged.bam rep1.samp.ind.bam rep2.samp.ind.bam rep3.samp.ind.bam -@ 24
```

### 7. Gene-body coverage profiles

Genes were intersected against nucleolus-associated domain (NAD) coordinates
to compute enrichment profiles (K4, K9, Pol1, Pol2 bigWigs) in vs. outside
NADs at gene TSS:

```shell
sed 's/^chr\|%$//g' nads.bed > nads_clean.bed
bedtools intersect -wa -a genes.gtf -b nads_clean.bed    > gene.nad.ov.gtf
bedtools intersect -wa -a genes.gtf -b nads_clean.bed -v > gene.nad.nonov.gtf

module load deeptools
computeMatrix scale-regions -S wt.bigWig actd.bigWig pol1.bigWig pol2.bigWig \
  -R gene.nad.ov.gtf gene.nad.nonov.gtf --skipZeros \
  -o matrix.gz -p 10 -a 2000 -b 2000

plotHeatmap -m matrix.gz -out heatmap.png --colorMap viridis --sortUsing sum \
  --legendLocation 'lower-center' --regionsLabel "NADs" "non-NADs" \
  --startLabel "Start" --endLabel "End" --samplesLabel "WT" "ActD" "Pol1" "Pol2"
```

The same approach was used with Pol1 MACS2 peaks (rather than NAD
coordinates) to compare K4/K9/Pol1/Pol2 signal at Pol1 peaks in vs. outside
NADs, using `computeMatrix reference-point` centered on peak midpoints.

## Scripts in this folder

| Script | Purpose |
|---|---|
| `CutTag.sh` | Per-sample QC, trimming, alignment, sorting, and peak calling |
| `stats.R` | Aggregates alignment/duplication QC across all processed samples |
