# RNA-seq Differential Expression Analysis of NRDE2 Knockdown

## Overview
This project analzyes RNA-seq data comparing NRDE2 knockdown samples against control samples to identify differentailly expressed genes. This workflow uses transcript quantification data, imports counts into R, performs differential expression analysis, and generates visualizations to summarize the results.

## Objective
The goal of this prroject is to determine how NRDE2 knockdown affects gene expression samples compared with control samples

### Tools used
R
Salmon
DESeq2
ggplot2
pheatmap

## Workflow
1. Import transcript-level quantification files from Salmon.
2. Summarrize transcript-level abundance estimates to gene-level counts.
3. Run differential expression analysis with DESeq2.
4. Generate PCA, MA, and volcano plots.
5. Export significant differentially expressed genes forr downstream analysis.

## Skills Demonstrated
RNA-seq differential expression analysis, transcript-to-gene count summarization, statistical analysis of sequenecing data, data visualization in R, and reproducible bioinformatics workflow development.
