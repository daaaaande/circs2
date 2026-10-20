# circRNA Pipeline — Output Files & Downstream R Analysis

This document describes the output files produced by the pipeline and walks through a complete downstream analysis in R: data loading, QC, cross-tool comparison, and differential expression of circRNAs.

---

## 1. Output Files

### 1.1 Per-sample consensus files

```
results/consensus/{sample}/{sample}_consensus.tsv
```

One file per sample. Produced by `consensus_filter.py` after the 2-of-2 vote between CIRCexplorer2 and CIRI2. Only circRNAs detected by **both** tools are retained.

| Column | Description |
|---|---|
| `sample` | Sample identifier |
| `coord` | Genomic coordinate key: `chr:start-end:strand` (0-based start) |
| `chrom` | Chromosome |
| `start` | BSJ start, 0-based |
| `end` | BSJ end |
| `strand` | `+` or `-` |
| `geneName` | Host gene (from CIRCexplorer2 annotation) |
| `circType` | `circRNA` or `ciRNA` |
| `cx2_BSJ` | Back-splice junction reads (CIRCexplorer2) |
| `ciri2_BSJ` | Back-splice junction reads (CIRI2/CIRI-full) |
| `total_mapped` | Uniquely mapped reads from STAR (denominator for RPM) |
| `cx2_RPM` | RPM normalised from CIRCexplorer2 BSJ count |
| `ciri2_RPM` | RPM normalised from CIRI2 BSJ count |
| `mean_RPM` | Mean of `cx2_RPM` and `ciri2_RPM` |

### 1.2 Combined featureCounts file

```
results/featurecounts/all_samples_counts.txt
```

Standard featureCounts output for **linear** gene expression, all samples in one matrix. The first few metadata columns (Chr, Start, End, Strand, Length) are followed by one count column per sample. This is the linear RNA counterpart used to contextualise circRNA abundance.

---

## 2. R Analysis

### Setup and dependencies

```r
# Install if needed:
# BiocManager::install(c("DESeq2", "ComplexHeatmap", "EnhancedVolcano"))
# install.packages(c("tidyverse", "circlize", "ggrepel", "patchwork", "pheatmap"))

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(ggrepel)
library(patchwork)
library(circlize)
```

---

## 3. Part 1 — Data Loading & Preparation

This section loads all data, builds count matrices, and defines sample groups. Doing this cleanly in one place means every downstream step uses the same objects and any metadata change only needs to happen once.

### 3.1 Sample metadata

Define your groups here. The `condition` column drives all differential expression analyses below.

```r
# Edit this to match your experiment.
# Each row is one sample; names must match the sample IDs in the consensus files
# and the column headers in the featureCounts output.

sample_meta <- tribble(
  ~sample,         ~condition,  ~batch,
  "1412-10_S10",   "control",   "B1",
  "1412-20_S20",   "control",   "B1",
  "1412-30_S30",   "treated",   "B1",
  "1412-40_S40",   "treated",   "B2",
  "1412-50_S50",   "treated",   "B2",
  # ... add all samples
) %>%
  mutate(condition = factor(condition, levels = c("control", "treated")))
```

### 3.2 Load consensus circRNA files

Each per-sample TSV is read and stacked into a single long-format table. We keep only the RPM columns for the main analysis, but retain BSJ counts separately for DESeq2, which requires integer counts.

```r
consensus_dir <- "results/consensus"

# Read all per-sample files and combine
circ_long <- sample_meta$sample %>%
  set_names() %>%
  map_dfr(function(s) {
    path <- file.path(consensus_dir, s, paste0(s, "_consensus.tsv"))
    read_tsv(path, col_types = cols(
      sample        = col_character(),
      coord         = col_character(),
      chrom         = col_character(),
      start         = col_integer(),
      end           = col_integer(),
      strand        = col_character(),
      geneName      = col_character(),
      circType      = col_character(),
      cx2_BSJ       = col_integer(),
      ciri2_BSJ     = col_integer(),
      total_mapped  = col_integer(),
      cx2_RPM       = col_double(),
      ciri2_RPM     = col_double(),
      mean_RPM      = col_double()
    ))
  })

# Coordinate annotation table (one row per unique circRNA)
circ_annot <- circ_long %>%
  distinct(coord, chrom, start, end, strand, geneName, circType)
```

### 3.3 Build RPM matrices (one per tool + mean)

Wide matrices with circRNAs as rows and samples as columns. circRNAs not detected in a sample get RPM = 0.

```r
make_rpm_matrix <- function(df, rpm_col) {
  df %>%
    select(coord, sample, rpm = {{ rpm_col }}) %>%
    pivot_wider(names_from = sample, values_from = rpm, values_fill = 0) %>%
    column_to_rownames("coord") %>%
    as.matrix()
}

rpm_cx2   <- make_rpm_matrix(circ_long, cx2_RPM)
rpm_ciri2 <- make_rpm_matrix(circ_long, ciri2_RPM)
rpm_mean  <- make_rpm_matrix(circ_long, mean_RPM)

# Align columns to sample_meta order
rpm_cx2   <- rpm_cx2[,   sample_meta$sample]
rpm_ciri2 <- rpm_ciri2[, sample_meta$sample]
rpm_mean  <- rpm_mean[,  sample_meta$sample]
```

### 3.4 Build BSJ count matrices (for DESeq2)

DESeq2 needs raw integer counts — never RPM. We use CIRCexplorer2 and CIRI2 BSJ counts separately, with the same zero-fill approach.

```r
make_count_matrix <- function(df, count_col) {
  df %>%
    select(coord, sample, count = {{ count_col }}) %>%
    pivot_wider(names_from = sample, values_from = count, values_fill = 0L) %>%
    column_to_rownames("coord") %>%
    as.matrix()
}

counts_cx2   <- make_count_matrix(circ_long, cx2_BSJ)
counts_ciri2 <- make_count_matrix(circ_long, ciri2_BSJ)

counts_cx2   <- counts_cx2[,   sample_meta$sample]
counts_ciri2 <- counts_ciri2[, sample_meta$sample]
```

### 3.5 Load linear featureCounts

The featureCounts file has metadata columns before the count columns. We strip those and align the sample names.

```r
fc_raw <- read_tsv(
  "results/featurecounts/all_samples_counts.txt",
  comment = "#"
)

# featureCounts prepends the full BAM path as column names — strip to sample ID
colnames(fc_raw) <- colnames(fc_raw) %>%
  basename() %>%
  str_remove("\\.Aligned\\..*\\.bam$")   # adjust regex to your BAM naming

linear_counts <- fc_raw %>%
  select(Geneid, all_of(sample_meta$sample)) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# RPM normalisation for linear counts (for comparison plots)
lib_sizes    <- colSums(linear_counts)
linear_rpm   <- sweep(linear_counts, 2, lib_sizes / 1e6, "/")
```

---

## 4. Part 2 — Simple QC

QC before any statistical analysis catches library size outliers, tool disagreements, and batch effects early. The plots below give a quick read on data quality.

### 4.1 circRNA detection counts per sample and tool

How many circRNAs were detected per sample? Large differences may flag failed samples or extreme library size differences.

```r
detection_counts <- circ_long %>%
  group_by(sample) %>%
  summarise(
    n_cx2   = sum(cx2_BSJ   > 0),
    n_ciri2 = sum(ciri2_BSJ > 0),
    n_consensus = n()
  ) %>%
  left_join(sample_meta, by = "sample") %>%
  pivot_longer(cols = starts_with("n_"), names_to = "tool", values_to = "n_detected") %>%
  mutate(tool = recode(tool,
    n_cx2       = "CIRCexplorer2",
    n_ciri2     = "CIRI2",
    n_consensus = "Consensus (2/2)"
  ))

ggplot(detection_counts, aes(x = sample, y = n_detected, fill = condition)) +
  geom_col() +
  facet_wrap(~tool, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "circRNA detection counts per sample",
    x = NULL, y = "Number of circRNAs detected"
  )
```

### 4.2 BSJ RPM distribution — boxplots across samples

Log-transformed RPM distributions should be roughly comparable across samples after normalisation. Heavy outliers suggest a problem with a specific sample or library.

```r
# Reshape all three RPM matrices for plotting
rpm_long <- bind_rows(
  as.data.frame(rpm_cx2)   %>% rownames_to_column("coord") %>%
    pivot_longer(-coord, names_to = "sample", values_to = "RPM") %>%
    mutate(tool = "CIRCexplorer2"),
  as.data.frame(rpm_ciri2) %>% rownames_to_column("coord") %>%
    pivot_longer(-coord, names_to = "sample", values_to = "RPM") %>%
    mutate(tool = "CIRI2"),
  as.data.frame(rpm_mean)  %>% rownames_to_column("coord") %>%
    pivot_longer(-coord, names_to = "sample", values_to = "RPM") %>%
    mutate(tool = "Mean (consensus)")
) %>%
  filter(RPM > 0) %>%
  left_join(sample_meta, by = "sample")

ggplot(rpm_long, aes(x = sample, y = log2(RPM + 0.01), fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4) +
  facet_wrap(~tool, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "circRNA BSJ RPM distributions per sample",
    x = NULL, y = "log2(RPM + 0.01)"
  )
```

### 4.3 circRNA vs linear RNA expression — global comparison

Do samples with more linear reads also have more circRNA signal? A strong deviation from the trend might indicate a biological effect or a technical artefact.

```r
# Mean RPM across all circRNAs per sample
circ_mean_per_sample <- circ_long %>%
  group_by(sample) %>%
  summarise(mean_circ_RPM = mean(mean_RPM[mean_RPM > 0]))

# Mean RPM across all linear genes per sample
linear_mean_per_sample <- tibble(
  sample          = colnames(linear_rpm),
  mean_linear_RPM = colMeans(linear_rpm[rowSums(linear_rpm) > 0, ])
)

circ_vs_linear <- left_join(circ_mean_per_sample, linear_mean_per_sample, by = "sample") %>%
  left_join(sample_meta, by = "sample")

ggplot(circ_vs_linear, aes(x = log2(mean_linear_RPM), y = log2(mean_circ_RPM),
                            colour = condition, label = sample)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey50", linetype = "dashed") +
  geom_text_repel(size = 3) +
  theme_bw() +
  labs(
    title  = "Mean circRNA RPM vs mean linear RNA RPM per sample",
    x      = "log2(mean linear RPM)",
    y      = "log2(mean circRNA RPM)"
  )
```

### 4.4 CIRCexplorer2 vs CIRI2 RPM correlation

For the 2-of-2 consensus circRNAs, how well do the two tools agree on abundance? Systematic divergence could indicate a tool-specific bias.

```r
# Per-circRNA mean RPM across samples for each tool
tool_comparison <- tibble(
  coord     = rownames(rpm_cx2),
  mean_cx2  = rowMeans(rpm_cx2),
  mean_ciri2 = rowMeans(rpm_ciri2)
) %>%
  filter(mean_cx2 > 0, mean_ciri2 > 0)

ggplot(tool_comparison, aes(x = log2(mean_cx2 + 0.01), y = log2(mean_ciri2 + 0.01))) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  geom_smooth(method = "lm", colour = "steelblue") +
  theme_bw() +
  labs(
    title = "CIRCexplorer2 vs CIRI2 mean RPM (consensus circRNAs)",
    x     = "log2(mean CIRCexplorer2 RPM)",
    y     = "log2(mean CIRI2 RPM)"
  )

cor_val <- cor(log2(tool_comparison$mean_cx2 + 0.01),
               log2(tool_comparison$mean_ciri2 + 0.01),
               method = "pearson")
message(sprintf("Pearson r = %.3f", cor_val))
```

### 4.5 Heatmap — top variable circRNAs (mean RPM matrix)

A quick unsupervised look at the data structure. Using the top 100 most variable circRNAs. If samples cluster by condition without supervision, the signal is likely real.

```r
# Select top 100 most variable circRNAs by inter-sample variance of log2 RPM
log2_rpm_mean <- log2(rpm_mean + 0.01)
top100_var    <- order(apply(log2_rpm_mean, 1, var), decreasing = TRUE)[1:100]
mat_top100    <- log2_rpm_mean[top100_var, ]

# Annotation for columns (samples)
col_annot <- sample_meta %>%
  select(sample, condition) %>%
  column_to_rownames("sample")

pheatmap(
  mat_top100,
  annotation_col   = col_annot,
  scale            = "row",       # z-score per circRNA for visual contrast
  show_rownames    = FALSE,
  clustering_method = "ward.D2",
  fontsize_col     = 8,
  main             = "Top 100 variable circRNAs — mean RPM (z-scored)"
)
```

---

## 5. Part 3 — Cross-dataset Comparison

Here we compare circRNA and linear RNA expression in more detail: are the most abundant circRNAs in the same genes as the most expressed linear transcripts? This can reveal host gene effects.

### 5.1 Host gene linear expression vs circRNA abundance

For each consensus circRNA, look up the linear expression of its host gene and plot both.

```r
# Mean linear RPM per gene across all samples
gene_linear_mean <- tibble(
  geneName       = rownames(linear_rpm),
  linear_mean_RPM = rowMeans(linear_rpm)
)

# Mean circRNA RPM per host gene (sum of all circRNAs from that gene)
circ_by_gene <- circ_long %>%
  group_by(geneName, sample) %>%
  summarise(total_circ_RPM = sum(mean_RPM), .groups = "drop") %>%
  group_by(geneName) %>%
  summarise(circ_mean_RPM = mean(total_circ_RPM))

host_comparison <- inner_join(circ_by_gene, gene_linear_mean, by = "geneName") %>%
  filter(circ_mean_RPM > 0, linear_mean_RPM > 0)

ggplot(host_comparison, aes(x = log2(linear_mean_RPM + 0.01),
                             y = log2(circ_mean_RPM  + 0.01))) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", colour = "tomato") +
  geom_text_repel(
    data = slice_max(host_comparison, circ_mean_RPM, n = 15),
    aes(label = geneName), size = 2.8
  ) +
  theme_bw() +
  labs(
    title = "Host gene linear expression vs total circRNA output",
    x     = "log2(linear RPM)",
    y     = "log2(circRNA RPM, all isoforms summed)"
  )
```

### 5.2 Clustering comparison — circRNA tools and linear RNA side by side

Run PCA on each dataset separately and plot on a common layout. This lets you compare how samples separate in circRNA vs linear space.

```r
run_pca <- function(mat, meta, label) {
  # Filter to expressed features, log-transform, transpose for PCA
  keep <- rowSums(mat > 0) >= ncol(mat) * 0.5  # present in ≥50% of samples
  pca  <- prcomp(t(log2(mat[keep, ] + 0.01)), scale. = TRUE)
  as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(meta, by = "sample") %>%
    mutate(
      dataset = label,
      var1    = round(summary(pca)$importance[2, 1] * 100, 1),
      var2    = round(summary(pca)$importance[2, 2] * 100, 1)
    )
}

pca_cx2    <- run_pca(rpm_cx2,     sample_meta, "CIRCexplorer2")
pca_ciri2  <- run_pca(rpm_ciri2,   sample_meta, "CIRI2")
pca_mean   <- run_pca(rpm_mean,    sample_meta, "Consensus mean")
pca_linear <- run_pca(linear_rpm,  sample_meta, "Linear RNA")

pca_all <- bind_rows(pca_cx2, pca_ciri2, pca_mean, pca_linear)

ggplot(pca_all, aes(x = PC1, y = PC2, colour = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5) +
  facet_wrap(~dataset, scales = "free") +
  theme_bw() +
  labs(title = "PCA — circRNA tools vs linear RNA")
```

---

## 6. Part 4 — Differential Expression Analysis

We use DESeq2 for all DE analyses. DESeq2 was designed for RNA-seq count data and handles the overdispersion typical of BSJ read counts well, provided counts are not too sparse. We run it separately for CIRCexplorer2 counts, CIRI2 counts, and linear RNA, then compare results.

**Why DESeq2 for circRNA?** BSJ counts are count data with overdispersion, the same statistical model applies. The key requirement — integer counts and a sensible size factor — is met by using raw BSJ reads and the STAR-mapped read counts as a normalisation reference. Filter out very lowly detected circRNAs before running to reduce noise.

### 6.1 Helper function — run DESeq2

```r
run_deseq2 <- function(count_mat, meta, design_formula = ~ condition) {
  # Filter: keep circRNAs with ≥5 reads in ≥2 samples
  keep <- rowSums(count_mat >= 5) >= 2
  message(sprintf("Retaining %d / %d features after count filter", sum(keep), nrow(count_mat)))
  count_mat <- count_mat[keep, ]

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = meta,
    design    = design_formula
  )
  dds <- DESeq(dds)
  dds
}
```

### 6.2 Run DESeq2 on each dataset

```r
# Ensure row order of meta matches column order of matrices
meta_df <- as.data.frame(sample_meta) %>% column_to_rownames("sample")

dds_cx2    <- run_deseq2(counts_cx2[,   rownames(meta_df)], meta_df)
dds_ciri2  <- run_deseq2(counts_ciri2[, rownames(meta_df)], meta_df)
dds_linear <- run_deseq2(linear_counts[, rownames(meta_df)], meta_df)
```

### 6.3 Extract results

```r
get_results_df <- function(dds, contrast = c("condition", "treated", "control"), label) {
  res <- results(dds, contrast = contrast, alpha = 0.05)
  as.data.frame(res) %>%
    rownames_to_column("feature") %>%
    mutate(dataset = label) %>%
    arrange(padj)
}

res_cx2    <- get_results_df(dds_cx2,    label = "CIRCexplorer2")
res_ciri2  <- get_results_df(dds_ciri2,  label = "CIRI2")
res_linear <- get_results_df(dds_linear, label = "Linear RNA")

# How many significant at FDR 5%?
map_dfr(list(res_cx2, res_ciri2, res_linear), function(r) {
  tibble(
    dataset = r$dataset[1],
    n_sig   = sum(r$padj < 0.05 & abs(r$log2FoldChange) >= 1, na.rm = TRUE)
  )
})
```

### 6.4 Volcano plots

```r
make_volcano <- function(res_df, title) {
  EnhancedVolcano(
    res_df,
    lab        = res_df$feature,
    x          = "log2FoldChange",
    y          = "padj",
    title      = title,
    pCutoff    = 0.05,
    FCcutoff   = 1,
    pointSize  = 1.5,
    labSize    = 3,
    drawConnectors = TRUE,
    max.overlaps   = 15
  )
}

make_volcano(res_cx2,    "DESeq2 — CIRCexplorer2 BSJ counts")
make_volcano(res_ciri2,  "DESeq2 — CIRI2 BSJ counts")
make_volcano(res_linear, "DESeq2 — Linear RNA")
```

### 6.5 Heatmap of top 500 DE circRNAs per tool

Visualise expression patterns of the most significant circRNAs. Using variance-stabilised counts (`vst`) from DESeq2 gives better visual contrast than raw RPM for heatmaps.

```r
plot_de_heatmap <- function(dds, res_df, meta, n_top = 500, title) {
  sig <- res_df %>%
    filter(!is.na(padj)) %>%
    slice_min(padj, n = n_top) %>%
    pull(feature)

  if (length(sig) < 2) {
    message("Too few significant features for heatmap: ", title)
    return(invisible(NULL))
  }

  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)[intersect(sig, rownames(assay(vsd))), ]
  mat <- mat - rowMeans(mat)  # centre rows

  col_annot <- meta %>%
    select(condition) %>%
    as.data.frame()

  pheatmap(
    mat,
    annotation_col    = col_annot,
    show_rownames     = FALSE,
    clustering_method = "ward.D2",
    fontsize_col      = 8,
    main              = title,
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(50)
  )
}

plot_de_heatmap(dds_cx2,    res_cx2,    meta_df, n_top = 500,
                title = "Top 500 DE circRNAs — CIRCexplorer2")
plot_de_heatmap(dds_ciri2,  res_ciri2,  meta_df, n_top = 500,
                title = "Top 500 DE circRNAs — CIRI2")
plot_de_heatmap(dds_linear, res_linear, meta_df, n_top = 500,
                title = "Top 500 DE linear genes")
```

### 6.6 Overlap between tools — which DE circRNAs are called by both?

A circRNA called DE by both CIRCexplorer2 and CIRI2 is a high-confidence hit. Concordant direction of change increases confidence further.

```r
sig_cx2   <- res_cx2   %>% filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% pull(feature)
sig_ciri2 <- res_ciri2 %>% filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% pull(feature)

both_sig <- intersect(sig_cx2, sig_ciri2)
message(sprintf(
  "Significant DE circRNAs: CX2=%d, CIRI2=%d, both=%d",
  length(sig_cx2), length(sig_ciri2), length(both_sig)
))

# Compare log2FC direction for shared hits
concordance <- inner_join(
  res_cx2   %>% select(feature, lfc_cx2   = log2FoldChange, padj_cx2   = padj),
  res_ciri2 %>% select(feature, lfc_ciri2 = log2FoldChange, padj_ciri2 = padj),
  by = "feature"
) %>%
  filter(padj_cx2 < 0.05 | padj_ciri2 < 0.05) %>%
  mutate(
    sig_cx2   = padj_cx2   < 0.05 & abs(lfc_cx2)   >= 1,
    sig_ciri2 = padj_ciri2 < 0.05 & abs(lfc_ciri2) >= 1,
    concordant = sign(lfc_cx2) == sign(lfc_ciri2)
  )

ggplot(concordance, aes(x = lfc_cx2, y = lfc_ciri2,
                         colour = concordant, alpha = sig_cx2 | sig_ciri2)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  scale_alpha_manual(values = c(0.15, 0.8)) +
  theme_bw() +
  labs(
    title   = "log2FC concordance between tools (significant in either)",
    x       = "log2FC — CIRCexplorer2",
    y       = "log2FC — CIRI2",
    colour  = "Concordant direction"
  )
```

### 6.7 High-confidence DE circRNA table

Produce a final ranked table of circRNAs significant in both tools, annotated with host gene and linear DE status of that gene.

```r
# Linear DE results indexed by gene
linear_sig <- res_linear %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
  select(geneName = feature, linear_lfc = log2FoldChange, linear_padj = padj)

high_conf <- inner_join(
  res_cx2   %>% select(coord = feature, lfc_cx2   = log2FoldChange, padj_cx2   = padj),
  res_ciri2 %>% select(coord = feature, lfc_ciri2 = log2FoldChange, padj_ciri2 = padj),
  by = "coord"
) %>%
  filter(padj_cx2 < 0.05, padj_ciri2 < 0.05,
         abs(lfc_cx2) >= 1, abs(lfc_ciri2) >= 1,
         sign(lfc_cx2) == sign(lfc_ciri2)) %>%
  left_join(circ_annot, by = "coord") %>%
  left_join(linear_sig, by = "geneName") %>%
  mutate(
    mean_lfc          = (lfc_cx2 + lfc_ciri2) / 2,
    host_gene_also_DE = !is.na(linear_padj)
  ) %>%
  arrange(padj_cx2)

write_tsv(high_conf, "results/high_confidence_DE_circRNAs.tsv")
head(high_conf, 20)
```

---

## 7. Additional Visualisations

### 7.1 Circos-style genomic distribution of DE circRNAs

Gives an overview of whether DE circRNAs are scattered across the genome or concentrated in specific loci.

```r
library(circlize)

# Use high-confidence DE circRNAs
de_circ_bed <- high_conf %>%
  filter(!is.na(chrom)) %>%
  transmute(
    chr   = chrom,
    start = start,
    end   = end,
    value = mean_lfc
  )

circos.clear()
circos.par("track.height" = 0.1, start.degree = 90)
circos.initializeWithIdeogram(species = "hg38", plotType = c("ideogram", "labels"))

circos.genomicTrack(
  de_circ_bed,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, pch = 16, cex = 0.5,
                         col = ifelse(value[[1]] > 0, "firebrick", "steelblue"), ...)
  }
)
title("Genomic distribution of high-confidence DE circRNAs\n(red = up, blue = down)")
```

### 7.2 Top DE circRNAs — expression across samples (dot plot)

```r
top_hits <- high_conf %>% slice_head(n = 20) %>% pull(coord)

circ_long %>%
  filter(coord %in% top_hits) %>%
  left_join(sample_meta, by = "sample") %>%
  mutate(coord = factor(coord, levels = rev(top_hits))) %>%
  ggplot(aes(x = sample, y = coord, size = log2(mean_RPM + 0.01),
             colour = condition)) +
  geom_point() +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Top 20 high-confidence DE circRNAs — mean RPM",
    x = NULL, y = NULL, size = "log2(RPM)"
  )
```

### 7.3 MA plot — mean expression vs fold change

Useful for checking whether DE calls are driven by highly expressed or lowly expressed circRNAs. Ideally significant hits should be spread across the expression range.

```r
make_ma <- function(res_df, title) {
  res_df %>%
    filter(!is.na(padj), !is.na(baseMean)) %>%
    mutate(significant = padj < 0.05 & abs(log2FoldChange) >= 1) %>%
    ggplot(aes(x = log2(baseMean + 1), y = log2FoldChange,
               colour = significant, alpha = significant)) +
    geom_point(size = 0.8) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("grey70", "firebrick")) +
    scale_alpha_manual(values  = c(0.3, 0.9)) +
    theme_bw() +
    labs(title = title, x = "log2(mean normalised count)", y = "log2 fold change")
}

make_ma(res_cx2,    "MA plot — CIRCexplorer2") /
make_ma(res_ciri2,  "MA plot — CIRI2") /
make_ma(res_linear, "MA plot — Linear RNA")
```

---

## 8. Session Info

Always record your session info for reproducibility.

```r
sessionInfo()
```