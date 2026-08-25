#!/usr/bin/env Rscript
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Splicing category classification from RREA.R/PROCESS_RREA.R
# output: per-read and per-gene exon/intron (E/I) ratio,
# intron-length tertile stratification, and read-termination
# position analysis.
#
# Includes deduplication of RREA output rows (a known source of
# duplicate read-exon pairs when findOverlaps() returns multiple
# transcript hits for the same exon in process_chromosome()) --
# see lines ~106-110.
#
# Figure output: SI Fig. 13G-H (confirmed).
##################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

######################################################################
## INTRON COVERAGE
######################################################################

# ── 1. Load and label each condition ─────────────────────────────────────────
path <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/out/"
list.files(path)

files <- list(
  K562_DMSO = paste0(path, "K562_DMSO_rrea_introns.rds"),
  K562_DMSO = paste0(path, "K562_DMSO_rep2_rrea_introns.rds"),
  K562_PlaB = paste0(path, "K562_PlaB_rep1_rrea_introns.rds"),
  K562_PlaB = paste0(path, "K562_PlaB_rep2_rrea_introns.rds")
  #K562_dRNA = paste0(path, "K562_dRNA_rep1_rrea_introns.rds"),
  #K562_dRNA = paste0(path, "K562_dRNA_rep2_rrea_introns.rds")
)

dt <- rbindlist(
  lapply(names(files), function(cond) {
    rds <- readRDS(files[[cond]])
    rds[, condition := cond]
  }),
  use.names = TRUE,
  fill      = TRUE
)

# ── 2. Identify transcripts with next_read_intron > 0 in ANY condition ────────
keep_transcripts <- dt[
  !is.na(next_read_intron) & next_read_intron > 0,
  unique(transcript_id)
]

# ── Transcripts with at least one non-NA next_read_intron ────────────────────
has_non_na <- dt[
  !is.na(next_read_intron),
  unique(transcript_id)
]

# ── 3. Filter — keep all rows for those transcripts across all conditions ─────
dt_filtered <- dt[
  transcript_id %in% keep_transcripts &
    transcript_id %in% has_non_na
]

cat(sprintf(
  "Kept %d / %d transcripts\n",
  uniqueN(dt_filtered$transcript_id),
  uniqueN(dt$transcript_id)
))

# ── 4. Normalise read abundance per condition ────────────────────────────────────
# Count unique reads per condition, then compute each condition's fraction
# of its total so all conditions sum to 1

dt_filtered[, 
            read_fraction := 1 / uniqueN(read_name), 
            by = condition
]

## verify that values sum to 1
dt_filtered[, .(total = sum(unique(read_fraction) * uniqueN(read_name))), by = condition]
#    condition total
# 1:  GM19102     1
# 2: K562_4sU     1

## deduplicate
cat("Rows before deduplication:", nrow(dt_filtered), "\n")
dt_filtered <- unique(dt_filtered)
cat("Rows after deduplication:", nrow(dt_filtered), "\n")
cat("Duplicate rows removed:", nrow(dt_filtered) - nrow(unique(dt_filtered)), "\n")

# ── 5. Plot histogram of intron lengths vs fraction ────────────────────────────────────

dt_plot <- dt_filtered[!is.na(next_read_intron) & next_read_intron > 0]

ggplot(dt_plot, aes(x = next_read_intron, weight = read_fraction, colour = condition)) +
  geom_freqpoly(bins = 50, linewidth = 0.8) +
  scale_x_log10(labels = scales::comma) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Intron length in read (log10)",
    y      = "Fraction of reads",
    colour = "Condition",
    title  = "Distribution of intron lengths across conditions"
  ) +
  theme_bw()

# Check the distribution of values
dt_filtered[, .(
  n_total    = .N,
  n_na       = sum(is.na(next_read_intron)),
  n_zero     = sum(next_read_intron == 0, na.rm = TRUE),
  n_negative = sum(next_read_intron < 0,  na.rm = TRUE),
  n_valid    = sum(next_read_intron > 0,  na.rm = TRUE)
), by = condition]


# ── 6. Compute statisicts for single barplot of RI, unspliced ────────────────────────────────────────

## filter conditions of interest (if needed)
conditions_of_interest <- c("DMSO", "PlaB", "4sU", "GM19102")
dt_plot <- dt_plot[condition %in% conditions_of_interest]

## get intron fraction
dt_plot[, intron_coverage_fraction := next_read_intron / next_ref_intron]

# Summarise to one row per read ───────────────────────────────────────────────
read_intron_summary <- dt_plot[!is.na(next_ref_intron) & next_ref_intron > 0, .(
  mean_intron_fraction   = mean(intron_coverage_fraction,   na.rm = TRUE),
  median_intron_fraction = median(intron_coverage_fraction, na.rm = TRUE),
  max_intron_fraction    = max(intron_coverage_fraction,    na.rm = TRUE),
  n_introns              = .N
), by = .(read_name, condition)]

read_intron_summary[, splicing_category := fcase(
  mean_intron_fraction <  0.5,                           "unspliced",
  mean_intron_fraction >= 0.5 & mean_intron_fraction < 0.75, "partially_retained",
  mean_intron_fraction >= 0.75,                          "retained"
)]

read_intron_summary[, .(
  n_reads              = .N,
  mean_frac            = mean(mean_intron_fraction),
  n_unspliced            = sum(splicing_category == "unspliced"),
  n_partially_retained = sum(splicing_category == "partially_retained"),
  n_retained           = sum(splicing_category == "retained")
), by = condition]

# ── Compute directly in long format ──────────────────────────────────────────
category_long <- read_intron_summary[, .(
  n_reads = .N
), by = .(condition, splicing_category)]

# Add total per condition for proportion calculation
category_long[, total := sum(n_reads), by = condition]
category_long[, proportion := n_reads / total]

# Set factor order
category_long[, splicing_category := factor(
  splicing_category,
  levels = c("unspliced", "partially_retained", "retained")
)]

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(category_long, aes(x = condition, y = proportion, fill = splicing_category)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(
    "unspliced"            = "#2166ac",
    "partially_retained" = "#f4a582",
    "retained"           = "#d6604d"
  )) +
  labs(
    x     = "Condition",
    y     = "Fraction of reads",
    fill  = "Splicing category",
    title = "Intron retention categories by condition"
  ) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# ── 7. Compute statisicts for barplot of RI, unspliced, by length ────────────────────────────────────────

dt_plot <- dt_plot %>% 
  dplyr::as_tibble() %>% 
  mutate(
    length = dplyr::case_when(
      next_ref_intron >= quantile(next_ref_intron, 2/3) & next_ref_intron <= max(next_ref_intron) ~ "long",
      next_ref_intron >= quantile(next_ref_intron, 1/3) & next_ref_intron < quantile(next_ref_intron, 2/3) ~ "average",
      next_ref_intron >= min(next_ref_intron) & next_ref_intron< quantile(next_ref_intron, 1/3) ~ "short"
    ),
    length = factor(length, levels = c("short", "average", "long"))
  )

# ── Recompute intron coverage fraction and splicing category ───────────────
dt_plot <- data.table(dt_plot)
read_intron_summary_by_length <- dt_plot[
  !is.na(next_ref_intron) & next_ref_intron > 0,
  .(
    mean_intron_fraction = mean(next_read_intron / next_ref_intron, na.rm = TRUE)
  ),
  by = .(read_name, condition, length)
][, splicing_category := fcase(
  mean_intron_fraction <  0.5,                                "unspliced",
  mean_intron_fraction >= 0.5 & mean_intron_fraction < 0.75,  "partially_retained",
  mean_intron_fraction >= 0.75,                               "retained"
)]

# ──  Compute proportions ────────────────────────────────────────────────────
category_by_length <- read_intron_summary_by_length[, .(
  n_reads = .N
), by = .(condition, length, splicing_category)]

category_by_length[, total := sum(n_reads), by = .(condition, length)]
category_by_length[, proportion := n_reads / total]

# Set factor orders for logical plotting order
category_by_length[, splicing_category := factor(
  splicing_category,
  levels = c("unspliced", "partially_retained", "retained")
)]

category_by_length[, length := factor(
  length,
  levels = c("short", "average", "long")
)]

# ── Compute label data using unique introns ───────────────────────────────────
length_labels <- dt_plot[
  !is.na(next_ref_intron) & next_ref_intron > 0,
  .SD[1],   # take one row per unique intron
  by = .(chromosome, transcript_id, exon_number, length)
][, .(
  n_introns   = .N,
  mean_length = round(mean(next_ref_intron, na.rm = TRUE))
),
by = length
][, label := sprintf("n = %s\nmean = %s bp",
                     scales::comma(n_introns),
                     scales::comma(mean_length))]

# Set factor order to match plot
length_labels[, length := factor(
  length,
  levels = c("short", "average", "long")
)]

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(category_by_length,
       aes(x = condition, y = proportion, fill = splicing_category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_label(
    data        = length_labels,
    aes(label   = label),
    x           = 1.5,        # center horizontally (middle condition)
    y           = 0.98,      # center vertically
    hjust       = 0.5,
    vjust       = 0.5,
    inherit.aes = FALSE,
    size        = 3,
    colour      = "grey30",
    fill        = "white",
    label.size  = 0.3,      # border thickness
    label.padding = unit(0.4, "lines")   # padding around text
  ) +
  facet_wrap(~ length,
             labeller = labeller(length = c(
               short   = "Short (< 1kb)",
               average = "Average (1–10kb)",
               long    = "Long (> 10kb)"
             ))) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(
    "unspliced"            = "#2166ac",
    "partially_retained" = "#f4a582",
    "retained"           = "#d6604d"
  )) +
  labs(
    x     = "Condition",
    y     = "Fraction of reads",
    fill  = "Splicing category",
    title = "Intron retention by condition and intron length"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

# ── 8. Compute intron coverage fraction and read fraction ────────────────────────

# Summarise to one row per read per intron length category
read_summary <- dt_plot[, .(
  intron_coverage_fraction = mean(intron_coverage_fraction, na.rm = TRUE),
  read_fraction            = sum(read_fraction)
), by = .(read_name, condition, length)]

# Set factor order
read_summary[, length:= factor(
  length,
  levels = c("short", "average", "long")
)]

# ── Plot hist ──────────────────────────────────────────────────────────────────────

# Recompute length_labels if needed
length_labels <- dt_plot[
  !is.na(next_ref_intron) & next_ref_intron > 0,
  .SD[1],
  by = .(chromosome, transcript_id, exon_number, length)
][, .(
  n_introns   = .N,
  mean_length = round(mean(next_ref_intron, na.rm = TRUE))
),
by = length
][, label := sprintf("n = %s\nmean = %s bp",
                     scales::comma(n_introns),
                     scales::comma(mean_length))]

length_labels[,length := factor(
  length,
  levels = c("short", "average", "long")
)]

ggplot(read_summary,
       aes(x = intron_coverage_fraction, colour = condition)) +
  geom_freqpoly(aes(weight = read_fraction), bins = 50, linewidth = 0.8) +
  geom_label(
    data          = length_labels,
    aes(x = 0.5, y = Inf, label = label),
    hjust         = 0.5,
    vjust         = 1.2,
    inherit.aes   = FALSE,
    size          = 3,
    colour        = "grey30",
    fill          = "white",
    label.size    = 0,
    label.padding = unit(0.4, "lines")
  ) +
  facet_wrap(~ length,
             scales   = "free_y",
             labeller = labeller(length = c(
               short   = "Short (< 1kb)",
               average = "Average (1–10kb)",
               long    = "Long (> 10kb)"
             ))) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Intron coverage fraction",
    y      = "Fraction of reads",
    colour = "Condition",
    title  = "Intron coverage fraction by condition and intron length"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

# ── Plot density ──────────────────────────────────────────────────────────────────────
ggplot(read_summary,
       aes(x = intron_coverage_fraction, colour = condition)) +
  geom_density(aes(weight = read_fraction), linewidth = 0.8, adjust = 1.5) +
  facet_wrap(~ length,
             scales   = "free_y",
             labeller = labeller(length = c(
               short   = "Short (< 1kb)",
               average = "Average (1–10kb)",
               long    = "Long (> 10kb)"
             ))) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Intron coverage fraction",
    y      = "Density",
    colour = "Condition",
    title  = "Intron coverage fraction by condition and intron length"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

######################################################################
## EXON COVERAGE
######################################################################

# ── 1. Compute exon coverage fraction ─────────────────────────────────────────
exon_dt <- dt_filtered[
  condition != "dRNA" & 
    !is.na(read_start) & !is.na(read_end) &
    !is.na(ref_start)  & !is.na(ref_end)
][, `:=`(
  read_exon_length = as.numeric(read_end - read_start),
  ref_exon_length  = as.numeric(ref_end  - ref_start)
)][ref_exon_length > 0 & read_start != read_end]

exon_dt[, exon_coverage_fraction := read_exon_length / ref_exon_length]

# ── 2. Define exon length categories ──────────────────────────────────────────
# Check quantiles first
exon_dt <- exon_dt %>% 
  dplyr::as_tibble() %>% 
  mutate(
    length = dplyr::case_when(
      ref_exon_length >= quantile(ref_exon_length, 2/3) & ref_exon_length <= max(ref_exon_length) ~ "long",
      ref_exon_length >= quantile(ref_exon_length, 1/3) & ref_exon_length < quantile(ref_exon_length, 2/3) ~ "average",
      ref_exon_length >= min(ref_exon_length) & ref_exon_length< quantile(ref_exon_length, 1/3) ~ "short"
    ),
    length = factor(length, levels = c("short", "average", "long"))
  )


# ── 3. Normalise read fraction across all conditions ──────────────────────────
exon_dt <- data.table(exon_dt)
total_reads <- exon_dt[, uniqueN(read_name), by = condition][, sum(V1)]
exon_dt[, read_fraction := 1 / total_reads]

# ── 4. Summarise to one row per read per exon length category ─────────────────
exon_summary <- exon_dt[, .(
  exon_coverage_fraction = mean(exon_coverage_fraction, na.rm = TRUE),
  read_fraction          = sum(read_fraction)
), by = .(read_name, condition, length)]

exon_summary[, length := factor(
  length,
  levels = c("short", "average", "long")
)]

# ── 5. Compute labels ─────────────────────────────────────────────────────────
# Exon label — unique exons seen by at least one read
exon_labels <- exon_dt[, .SD[1],
                       by = .(chromosome, transcript_id, exon_number, length)
][, .(
  n_exons     = .N,
  mean_length = round(mean(ref_exon_length, na.rm = TRUE))
),
by = length
][, label := sprintf("n = %s\nmean = %s bp",
                     scales::comma(n_exons),
                     scales::comma(mean_length))]


exon_labels[, length := factor(
  length,
  levels = c("short", "average", "long")
)]

# ── 6. Plot ───────────────────────────────────────────────────────────────────
ggplot(exon_summary,
       aes(x = exon_coverage_fraction, colour = condition)) +
  geom_freqpoly(aes(weight = read_fraction), bins = 50, linewidth = 0.8) +
  geom_label(
    data          = exon_labels,
    aes(x = 0.5, y = Inf, label = label),
    hjust         = 0.5,
    vjust         = 1.2,
    inherit.aes   = FALSE,
    size          = 3,
    colour        = "grey30",
    fill          = "white",
    label.size    = 0,
    label.padding = unit(0.4, "lines")
  ) +
  facet_wrap(~ length,
             scales   = "free_y",
             labeller = labeller(length = c(
               short   = "Short Exons",
               average = "Average Exons",
               long    = "Long Exons"
             ))) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Exon coverage fraction",
    y      = "Fraction of reads",
    colour = "Condition",
    title  = "Exon coverage fraction by condition and exon length"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

print(exon_labels)
print(length_labels)

######################################################################
## TERMINATION PLOT
######################################################################

# ── 1. Filter to reads that terminate inside an intron ────────────────────────
# These are reads where next_read_intron > 0 but < next_ref_intron
# i.e. the read entered the intron but didn't make it to the next exon
termination_dt <- dt_plot[
  !is.na(next_read_intron)  & next_read_intron  > 0 &
    !is.na(next_ref_intron)   & next_ref_intron   > 0 &
    next_read_intron < next_ref_intron 
]

# ── 2. Normalise read fraction ────────────────────────────────────────────────
total_reads <- termination_dt[, uniqueN(read_name), by = condition][, sum(V1)]
termination_dt[, read_fraction := 1 / total_reads]

# ── 3. Set factor order ───────────────────────────────────────────────────────
termination_dt[,length := factor(
  length,
  levels = c("short", "average", "long")
)]

termination_dt <- termination_dt[next_read_intron > 1]
# ── 4. Plot pileup by nucleotide position ─────────────────────────────────────
ggplot(termination_dt,
       aes(x = next_read_intron, colour = condition, linetype = condition)) +
  geom_freqpoly(aes(weight = read_fraction),bins = 100, linewidth = 0.8) +
  scale_linetype_manual(values = c(
    "K562_DMSO"    = "11",
    "K562_PlaB"    = "solid"
  )) +
  geom_label(
    data          = length_labels,
    aes(x = Inf, y = Inf, label = label),
    hjust         = 1.1,
    vjust         = 1.2,
    inherit.aes   = FALSE,
    size          = 3,
    colour        = "grey30",
    fill          = "white",
    label.size    = 0,
    label.padding = unit(0.4, "lines")
  ) +
  facet_wrap(~ length,
             scales   = "free",
             labeller = labeller(length= c(
               short   = "Short Introns",
               average = "Average Introns",
               long    = "Long Introns"
             ))) +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Nucleotide position inside intron (bp)",
    y      = "Fraction of reads",
    colour = "Condition",
    title  = "Read termination pileup within introns by intron length"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

######################################################################
## EI PLOT
######################################################################

library(plyranges)

# ── Genes with any read coverage ──────────────────────────────────────────────
covered_genes <- dt_filtered[, unique(gene_id)]

# ── Helper: assign tertile-based length categories ────────────────────────────
assign_length_category <- function(dt, col) {
  q33 <- quantile(dt[[col]], 1/3, na.rm = TRUE)
  q66 <- quantile(dt[[col]], 2/3, na.rm = TRUE)
  
  dt[, intron_length_category := fcase(
    get(col) <  q33, "short",
    get(col) <  q66, "average",
    get(col) >= q66, "long"
  )]
  
  dt[, intron_length_category := factor(
    intron_length_category,
    levels = c("short", "average", "long")
  )]
  
  dt[]
}

##-------------------------------------------------------## Reference E/I ratio per gene

# bring in GTF of interest
gen.root <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/"
gen.ls <- list.files(paste0(gen.root))
gen.ls

gtf_gr <- plyranges::read_gff(paste0(gen.root,"/", gen.ls[1]))
seqlevelsStyle(gtf_gr) <- "NCBI"

## Get lengths of features
length(na.omit(unique(gtf_gr$transcript_id)))
length(na.omit(unique(gtf_gr$exon_id)))
length(na.omit(unique(gtf_gr$gene_id)))
length(na.omit(unique(gtf_gr$gene_name)))

## get exons in read data only
covered_exons <- dt_filtered[,unique(exon_id)]

## get exon total lengths for each gene
exon_length <- as.data.frame(gtf_gr) %>%
  filter(type == "exon" & exon_id %in% covered_exons) %>%
  distinct(exon_id, .keep_all = TRUE) %>%          # one row per unique exon
  group_by(gene_id) %>%
  reframe(
    score = sum(width(IRanges::reduce(IRanges(start, end))))  # merge overlapping exons
  )

## get gene total lengths for each gene
gene_length <-  as.data.frame(gtf_gr) %>% filter(type == "exon" & exon_id %in% covered_exons)  %>%  group_by(gene_id) %>% reframe(
  score = max(end)-min(start), chr = seqnames, start = start, end = end
) %>% distinct(gene_id, score, .keep_all = TRUE)

## join table and get ratio
n <-merge(exon_length, gene_length , by = c("gene_id")) %>% mutate(total_exon_length = score.x, total_gene_length = score.y) %>% select(-score.x, -score.y)
ref_ei <- n %>% mutate(ref_ei_ratio  = (total_exon_length/(total_gene_length-total_exon_length) ) ) ## get EI ratio

ref_ei <- data.table(ref_ei)
ref_ei <- ref_ei[,total_intron_length := total_gene_length - total_exon_length] ## get intron lengths
               
##-------------------------------------------------------## Read E/I ratio per gene

# ── 1. Total exon length per gene per read from read coordinates ──────────────
read_exon_lengths <- dt_filtered[
  gene_id %in% covered_genes,
  .SD[1],
  by = .(gene_id, condition, read_name, exon_id)
][, .(
  total_read_exon_length = sum(read_end - read_start)
), by = .(gene_id, condition, read_name)]


# ── 2. Total intron length per gene per read from next_read_intron ────────────

# ── Read intron length per gene per read ──────────────────────────────────────
read_intron_lengths <- dt_plot[
 gene_id %in% covered_genes &
    !is.na(next_read_intron) & next_read_intron > 0,
  .(
    total_read_intron_length = sum(next_read_intron)
  ),
  by = .(gene_id,condition, read_name)
]


# ── 3. Join and calculate read E/I ratio ──────────────────────────────────────
read_ei <- merge(
  read_exon_lengths,
  read_intron_lengths,
  by = c("gene_id", "condition", "read_name")
)[total_read_intron_length > 0][,
                                read_ei_ratio := total_read_exon_length / total_read_intron_length
]


# ── Join and add intron length category ───────────────────────────────────────

ei_combined <- merge(
  read_ei,
  ref_ei[, .(gene_id, ref_ei_ratio, 
             total_intron_length)],
  by = c("gene_id")
)

# Intron length category based on mean intron length per gene
# Check quantiles first
ei_combined <- ei_combined %>% 
  dplyr::as_tibble() %>% 
  mutate(
    length = dplyr::case_when(
      total_intron_length >= quantile(total_intron_length, 2/3) & total_intron_length <= max(total_intron_length) ~ "long",
      total_intron_length >= quantile(total_intron_length, 1/3) & total_intron_length < quantile(total_intron_length, 2/3) ~ "average",
      total_intron_length >= min(total_intron_length) & total_intron_length< quantile(total_intron_length, 1/3) ~ "short"
    ),
    length = factor(length, levels = c("short", "average", "long"))
  )

ei_combined <- data.table(ei_combined)

# ── Melt for plotting ─────────────────────────────────────────────────────────

ei_long <- melt(
  ei_combined,
  id.vars       = c("gene_id", "condition","read_name",
                    "total_intron_length", "length"),
  measure.vars  = c("read_ei_ratio", "ref_ei_ratio"),
  variable.name = "source",
  value.name    = "ei_ratio"
)

## move Ref over to condition
ei_long[source == "ref_ei_ratio", condition := "Reference"]

## rename source
ei_long[, source := fcase(
  source == "read_ei_ratio", "Read",
  source == "ref_ei_ratio",  "Reference"
)]

## factor out lengths
ei_long[, length := factor(
  length,
  levels = c("short", "average", "long")
)]

## ── Plot Density First  ──────────────────────────────────────────────────────────────────────

# ── Compute stats  ────────────────────────

# One row per unique intron from the GTF
ref_introns_dt <- as.data.frame(gtf_gr) %>%
  filter(type == "exon" & exon_id %in% covered_exons) %>%
  distinct(exon_id, .keep_all = TRUE) %>%
  arrange(gene_id, start) %>%
  group_by(gene_id) %>%
  mutate(
    next_start      = lead(start),
    intron_length   = next_start - end
  ) %>%
  filter(!is.na(intron_length) & intron_length > 0) %>%
  ungroup() %>%
  data.table()

# Add length category using same breakpoints
ref_introns_dt  <- assign_length_category(ref_introns_dt,  "intron_length")

# ── N introns per category from GTF ──────────────────────────────────────────
n_introns_by_category <- ref_introns_dt[, .(
  n_introns = .N
), by = intron_length_category]

# ── Mean intron length per category from ei_long ─────────────────────────────
setnames(ei_long, "length", "intron_length_category")

mean_length_by_category <- ei_long[
  source == "Reference",
  .(mean_intron_bp = round(mean(total_intron_length, na.rm = TRUE))),
  by = intron_length_category
]

# ── Join and build label ──────────────────────────────────────────────────────
intron_stats <- merge(
  n_introns_by_category,
  mean_length_by_category,
  by = "intron_length_category"
)[, label := sprintf("N = %s introns", scales::comma(n_introns))]

intron_stats[, intron_length_category := factor(
  intron_length_category,
  levels = c("short", "average", "long")
)]

# ── Facet labels  ────────────────────────

facet_labels <- intron_stats[, setNames(
  sprintf("%s\n(mean intron: %s bp)", 
          tools::toTitleCase(as.character(intron_length_category)),
          scales::comma(mean_intron_bp)),
  as.character(intron_length_category)
)]

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(ei_long[!is.na(ei_ratio) & !is.na(total_intron_length)],
       aes(x = ei_ratio, colour = condition, linetype = source)) +
  geom_density(linewidth = 0.8, adjust = 1.5) +
  geom_label(
    data          = intron_stats,
    aes(label     = label),
    x             = Inf,
    y             = Inf,
    hjust         = 1.1,
    vjust         = 1.2,
    inherit.aes   = FALSE,
    size          = 3,
    colour        = "grey30",
    fill          = "white",
    label.size    = 0,
    label.padding = unit(0.4, "lines")
  ) +scale_y_continuous(breaks = function(x) round(seq(min(x), max(x), length.out = 4), 2))+
  facet_wrap(~ intron_length_category,
             labeller = labeller(intron_length_category = facet_labels)) +
  scale_y_continuous(
    limits = c(0, 1.5),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5),
    labels = c("0", "0.25", "0.50", "0.75", "1.0", "1.25", "1.50")
  )+
  
  scale_linetype_manual(values = c(
    "Read"      = "21",
    "Reference" = "solid"
  )) +
  scale_x_log10(
    limits = c(1e-02, 1e+04),
    breaks = c(1e-02, 1e-01, 1e+00, 1e+01, 1e+02, 1e+03, 1e+04),
    labels = scales::scientific
  ) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x        = "Exon / Intron ratio (log10)",
    y        = "Density",
    colour   = "Condition",
    linetype = "Source",
    title    = "Exon/Intron ratio per gene by intron length"
  ) +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial")
  )



# ── Plot dotplot ──────────────────────────────────────────────────────────────────────


ggplot(ei_long[!is.na(ei_ratio) & !is.na(intron_length_category)],
       aes(x = total_intron_length, y = ei_ratio, colour = condition)) +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_wrap(~ intron_length_category,
             scales   = "free",
             labeller = labeller(intron_length_category = c(
               short   = "Short (< 1kb)",
               average = "Average (1–10kb)",
               long    = "Long (> 10kb)"
             ))) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Total reference intron length (bp, log10)",
    y      = "E/I ratio (log10)",
    colour = "Condition",
    title  = "Exon/Intron ratio vs reference intron length by condition"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )



library(ggdist)
library(gghalves)
library(ggpubr)

# ── Plot violin plot ────────────────────────────────────────────────────────

comparisons <- list(
  c("Reference", "K562_DMSO"),   # adjust to match your exact names
  c("Reference", "K562_PlaB")
)
ggplot(ei_long[!is.na(ei_ratio) & !is.na(intron_length_category)],
       aes(x = condition, y = ei_ratio, fill = condition, colour = condition)) +
  
  # Half violin on the left — no outline
  geom_half_violin(
    side      = "l",
    alpha     = 0.4,
    linewidth = 0,
    trim      = FALSE,
    scale     = "width"
  ) +
  
  # Boxplot nudged to the right
  geom_boxplot(
    width         = 0.2,
    outlier.shape = NA,
    alpha         = 0.6,
    linewidth     = 0.5,
    position      = position_nudge(x = 0.15)
  ) +
  
  # N introns label in upper right
  geom_label(
    data          = intron_stats,
    aes(label     = label),
    x             = -Inf,
    y             = Inf,
    hjust         = -0.1,
    vjust         = 1.2,
    inherit.aes   = FALSE,
    size          = 3,
    colour        = "grey30",
    fill          = "white",
    label.size    = 0,
    label.padding = unit(0.4, "lines")
  ) +
  
  # Statistical comparisons
  geom_signif(
    comparisons      = comparisons,
    map_signif_level = TRUE,
    test             = "wilcox.test",
    step_increase    = 0.15,
    tip_length       = 0.01,
    size             = 0.5,
    textsize         = 3,
    vjust            = 0.5,
    y_position       = c(4, 2), 
    color = "black"
  ) +
  
  facet_wrap(~ intron_length_category,
             labeller = labeller(intron_length_category = facet_labels)) +
  scale_y_log10(
    limits = c(10^-2, 10^5),
    breaks = 10^(-1:5),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
  ) +
  scale_x_discrete(labels = c(
    "K562_DMSO" = "DMSO",
    "K562_PlaB" = "PlaB",
    "Reference" = "Reference"
  )) +
  scale_fill_brewer(palette   = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Condition",
    y      = "E/I ratio (log10)",
    fill   = "Condition",
    colour = "Condition",
    title  = "Exon/Intron ratio by condition and intron length"
  ) +
  theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12)
  )

######################################################################
## EI PLOT CONTROLS
######################################################################

# ── 1. Classify genes by intron length composition ────────────────────────────
# For each gene compute the proportion of short vs long introns

# ── 1. Compute intron composition per gene from GTF ───────────────────────────

# ref_introns_dt already has one row per unique intron from the GTF
q33 <- quantile(ref_introns_dt$intron_length, 1/3, na.rm = TRUE)
q66 <- quantile(ref_introns_dt$intron_length, 2/3, na.rm = TRUE)

cat(sprintf("Quantile breakpoints:\n  Short  < %s bp\n  Average %s-%s bp\n  Long   > %s bp\n",
            scales::comma(q33), scales::comma(q33),
            scales::comma(q66), scales::comma(q66)))

ref_gene_composition <- ref_introns_dt[, .(
  n_short   = sum(intron_length <  q33),
  n_average = sum(intron_length >= q33 & intron_length < q66),
  n_long    = sum(intron_length >= q66),
  n_total   = .N
), by = gene_id][, `:=`(
  pct_short = n_short / n_total,
  pct_long  = n_long  / n_total
)]

# ── 2. Classify genes into categories ────────────────────────────────────────
ref_gene_composition[, gene_category := fcase(
  pct_short >= 0.8,                        "short_only",
  pct_long  >= 0.8,                        "long_only",
  pct_short >= 0.33 & pct_long >= 0.33,      "mixed",
  default = NA_character_
)]

cat("Genes per category in GTF:\n")
print(ref_gene_composition[!is.na(gene_category), .N, by = gene_category])

# ── 3. Sample 100 genes per category ─────────────────────────────────────────
set.seed(42)

ref_gene_sets <- ref_gene_composition[
  !is.na(gene_category),
  .SD[sample(.N, min(.N, 1000))],
  by = gene_category
]

# ── 4. Subset read data using GTF-derived gene classifications ────────────────
read_selected <- dt_plot[
  gene_id %in% ref_gene_sets$gene_id &
    !is.na(next_ref_intron) & next_ref_intron > 0
]

# Join gene category from ref_gene_sets
read_selected <- merge(
  read_selected,
  ref_gene_sets[, .(gene_id, gene_category)],
  by = "gene_id"
)

cat("Read data rows per gene category:\n")
print(read_selected[, .N, by = gene_category])


# ── 1. Read E/I ratio for selected genes ──────────────────────────────────────

# Total exon length per gene per read per condition
selected_exon_lengths <- dt_filtered[
  gene_id %in% ref_gene_sets$gene_id & condition != "dRNA",
  .SD[1],
  by = .(gene_id, condition, read_name, exon_id)
][, .(
  total_read_exon_length = sum(read_end - read_start)
), by = .(gene_id, condition, read_name)]

# Total intron length per gene per read per condition
selected_intron_lengths <- dt_plot[
  gene_id %in% ref_gene_sets$gene_id & condition != "dRNA" &
    !is.na(next_read_intron) & next_read_intron > 0,
  .(
    total_read_intron_length = sum(next_read_intron)
  ),
  by = .(gene_id, condition, read_name)
]

# Join and compute read E/I ratio
selected_read_ei <- merge(
  selected_exon_lengths,
  selected_intron_lengths,
  by = c("gene_id", "condition", "read_name")
)[total_read_intron_length > 0][,
                                read_ei_ratio := total_read_exon_length / total_read_intron_length
]

# Join gene category
selected_read_ei <- merge(
  selected_read_ei,
  ref_gene_sets[, .(gene_id, gene_category)],
  by = "gene_id"
)

# ── 2. Reference E/I ratio for selected genes ─────────────────────────────────
selected_ref_ei <- ref_ei[
  gene_id %in% ref_gene_sets$gene_id
]

# Join gene category
selected_ref_ei <- merge(
  selected_ref_ei,
  ref_gene_sets[, .(gene_id, gene_category)],
  by = "gene_id"
)[, condition := "Reference"]

# ── 3. Melt and combine ───────────────────────────────────────────────────────
selected_read_long <- melt(
  selected_read_ei,
  id.vars      = c("gene_id", "condition", "read_name", "gene_category"),
  measure.vars = "read_ei_ratio",
  variable.name = "source",
  value.name    = "ei_ratio"
)[, source := "Read"]

selected_ref_long <- melt(
  selected_ref_ei,
  id.vars      = c("gene_id", "condition", "gene_category"),
  measure.vars = "ref_ei_ratio",
  variable.name = "source",
  value.name    = "ei_ratio"
)[, source := "Reference"]

selected_ei_long <- rbindlist(
  list(selected_read_long, selected_ref_long),
  use.names = TRUE,
  fill      = TRUE
)

# Set factor orders
selected_ei_long[, gene_category := factor(
  gene_category,
  levels = c("short_only", "mixed", "long_only")
)]

selected_ei_long[, condition := factor(
  condition,
  levels = c("K562_DMSO", "K562_PlaB", "Reference")
)]

# ── 4. Sanity check ───────────────────────────────────────────────────────────
selected_ei_long[, .(
  n_reads   = .N,
  median_ei = median(ei_ratio, na.rm = TRUE)
), by = .(gene_category, condition)]

selected_read_ei[, uniqueN(gene_id), by = .(gene_category, condition)]

# ── Join total intron length from ref_ei onto selected_ei_long ────────────────
selected_ei_long <- merge(
  selected_ei_long,
  ref_ei[, .(gene_id, total_intron_length)],
  by = "gene_id"
)

# ── density plot ────────────────
ggplot(selected_ei_long[!is.na(ei_ratio)],
       aes(x = ei_ratio, colour = condition, linetype = source)) +
  geom_density(linewidth = 0.8, adjust = 1.5) +
  facet_wrap(~ gene_category,
             scales   = "free_y",
             labeller = labeller(gene_category = c(
               short_only = "Short Only",
               mixed      = "Mixed",
               long_only  = "Long Only"
             ))) +
  scale_x_log10(
    limits = c(1e-02, 1e+04),
    breaks = 10^(-2:4),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_linetype_manual(values = c("Read" = "21", "Reference" = "solid")) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x        = "E/I ratio (log10)",
    y        = "Density",
    colour   = "Condition",
    linetype = "Source",
    title    = "E/I ratio by gene category"
  ) +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12)
  )


# Filter to mixed genes with intron termination events
mixed_termination <- dt_plot[
  gene_id %in% ref_gene_sets[gene_category == "mixed", gene_id] &
    !is.na(next_read_intron)  & next_read_intron  > 1  &
    !is.na(next_ref_intron)   & next_ref_intron   > 0  &
    next_read_intron < next_ref_intron &
    condition != "dRNA"
]

# Normalise read fraction
total_reads_mixed <- mixed_termination[, uniqueN(read_name), 
                                       by = condition][, sum(V1)]
mixed_termination[, read_fraction := 1 / total_reads_mixed]

# Add intron length category using same quantile breakpoints
mixed_termination[, intron_length_category := fcase(
  next_ref_intron <  q33,                           "short",
  next_ref_intron >= q33 & next_ref_intron <  q66,  "average",
  next_ref_intron >= q66,                           "long"
)]

mixed_termination[, intron_length_category := factor(
  intron_length_category,
  levels = c("short", "average", "long")
)]

ggplot(mixed_termination[!is.na(intron_length_category)],
       aes(x = next_read_intron, colour = condition, linetype = condition)) +
  geom_freqpoly(aes(weight = read_fraction), bins = 100, linewidth = 0.8) +
  facet_wrap(~ intron_length_category,
             scales   = "free",
             labeller = labeller(intron_length_category = c(
               short   = "Short Introns",
               average = "Average Introns",
               long    = "Long Introns"
             ))) +
  scale_x_continuous(labels = scales::comma) +
  scale_linetype_manual(values = c(
    "K562_DMSO" = "21",
    "K562_PlaB" = "21"
  )) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Nucleotide position inside intron (bp)",
    y      = "Fraction of reads",
    colour = "Condition",
    linetype = "Condition",
    title  = "Read termination pileup in mixed-intron genes"
  ) +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12)
  )


## ── Termination point vs reference intron length ──────────────────────────────
termination_scatter <- unique(dt_plot[
  !is.na(next_read_intron)  & next_read_intron  > 1  &
    !is.na(next_ref_intron)   & next_ref_intron   > 0  &
    next_read_intron < next_ref_intron &
    condition != "dRNA"
])

termination_scatter <- unique(
  termination_scatter,
  by = c("read_name", "gene_id", "condition", 
         "ref_start", "ref_end", "next_ref_start",
         "next_read_intron", "next_ref_intron")
)

# Verify
nrow(termination_scatter) - nrow(unique(termination_scatter))
termination_scatter[, uniqueN(gene_id), by = condition]

# Re-add intron length category
termination_scatter[, intron_length_category := fcase(
  next_ref_intron <  q33,                          "short",
  next_ref_intron >= q33 & next_ref_intron <  q66, "average",
  next_ref_intron >= q66,                          "long"
)]

termination_scatter[, intron_length_category := factor(
  intron_length_category,
  levels = c("short", "average", "long")
)]

# Re-add intron coverage percentile
termination_scatter[, intron_coverage_pct := (next_read_intron / next_ref_intron) * 100]

# Re-add read fraction
total_reads <- termination_scatter[, uniqueN(read_name), 
                                   by = condition][, sum(V1)]
termination_scatter[, read_fraction := 1 / total_reads]

# Verify
nrow(termination_scatter) - nrow(unique(termination_scatter))
termination_scatter[, uniqueN(gene_id), by = condition]

ggplot(termination_scatter,
       aes(x = next_ref_intron, y = next_read_intron, colour = condition)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", colour = "black", linewidth = 0.5) +
  scale_x_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Reference intron length (bp, log10)",
    y      = "Read termination point (bp, log10)",
    colour = "Condition",
    title  = "Read termination point vs reference intron length"
  ) +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12)
  )


# ── Compute intron coverage percentile ────────────────────────────────────────
termination_scatter[, intron_coverage_pct := (next_read_intron / next_ref_intron) * 100]

# Add intron length category
termination_scatter[, intron_length_category := fcase(
  next_ref_intron <  q33,                          "short",
  next_ref_intron >= q33 & next_ref_intron <  q66, "average",
  next_ref_intron >= q66,                          "long"
)]

termination_scatter[, intron_length_category := factor(
  intron_length_category,
  levels = c("short", "average", "long")
)]

# Normalise read fraction across all conditions
total_reads <- termination_scatter[, uniqueN(read_name), 
                                   by = condition][, sum(V1)]
termination_scatter[, read_fraction := 1 / total_reads]

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(termination_scatter[!is.na(intron_length_category)],
       aes(x = intron_coverage_pct, colour = condition)) +
  geom_density(aes(weight = read_fraction), linewidth = 0.8, adjust = 1.5) +
  facet_wrap(~ intron_length_category,
             labeller = labeller(intron_length_category = c(
               short   = "Short Introns",
               average = "Average Introns",
               long    = "Long Introns"
             ))) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Intron coverage (% of reference intron length)",
    y      = "Density",
    colour = "Condition",
    title  = "Read intron coverage as percentage of reference intron length"
  ) +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12)
  )


ggplot(termination_scatter[!is.na(intron_length_category)],
       aes(x = intron_coverage_pct, fill = condition, colour = condition)) +
  geom_histogram(
    aes(weight = read_fraction),
    bins     = 20,
    alpha    = 0.4,
    position = "identity"
  ) +
  facet_wrap(~ intron_length_category,
             labeller = labeller(intron_length_category = c(
               short   = "Short Introns",
               average = "Average Introns",
               long    = "Long Introns"
             ))) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  scale_fill_brewer(palette   = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x      = "Intron coverage (% of reference intron length)",
    y      = "Fraction of reads",
    fill   = "Condition",
    colour = "Condition",
    title  = "Read intron coverage as percentage of reference intron length"
  ) +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(size = 14),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", family = "Arial", size = 12)
  )

termination_scatter[, uniqueN(gene_id), by = condition]
nrow(termination_scatter) - nrow(unique(termination_scatter))
termination_scatter <- unique(dt_plot[
  !is.na(next_read_intron)  & next_read_intron  > 1  &
    !is.na(next_ref_intron)   & next_ref_intron   > 0  &
    next_read_intron < next_ref_intron &
    condition != "dRNA"
])

termination_scatter <- unique(
  termination_scatter,
  by = c("read_name", "gene_id", "condition", "ref_end", "next_ref_start")
)
