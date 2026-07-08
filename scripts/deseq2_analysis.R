############################################################
# RNA-seq Differential Expression Analysis Using DESeq2
#
# Author: Emil Mendoza
#
# Description:
# This script imports Salmon transcript quantification files,
# summarizes transcript abundance to gene-level counts using tximport,
# performs differential expression analysis with DESeq2,
# applies log2 fold change shrinkage with apeglm,
# and exports figures and result tables.
############################################################

#############################
# Load Required Packages
#############################

library(tximport)
library(readr)
library(DESeq2)
library(apeglm)

#############################
# Define File Paths
#############################

sample_file <- "pipeline/samplesheet.csv"
salmon_dir <- "results_rnaseq/star_salmon"
tx2gene_file <- file.path(salmon_dir, "tx2gene.tsv")

figures_dir <- "results/figures"
tables_dir <- "results/tables"

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

#############################
# Import Sample Metadata
#############################

samples <- read.csv(sample_file)
samples$sample <- as.character(samples$sample)

#############################
# Import Salmon Quantifications
#############################

files <- file.path(salmon_dir, samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read_tsv(
  tx2gene_file,
  col_names = c("TXNAME", "GENEID")
)

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene
)

#############################
# Create DESeq2 Dataset
#############################

dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~ condition
)

dds$condition <- relevel(dds$condition, ref = "control")

#############################
# Run Differential Expression
#############################

dds <- DESeq(dds)

res <- results(dds)

resLFC <- lfcShrink(
  dds,
  coef = 2,
  type = "apeglm"
)

#############################
# Summarize Results
#############################

total_genes <- nrow(resLFC)

sig_genes <- sum(resLFC$padj < 0.1, na.rm = TRUE)

upregulated <- sum(
  resLFC$padj < 0.1 &
    resLFC$log2FoldChange > 0,
  na.rm = TRUE
)

downregulated <- sum(
  resLFC$padj < 0.1 &
    resLFC$log2FoldChange < 0,
  na.rm = TRUE
)

cat("\nAnalysis Summary\n")
cat("----------------\n")
cat("Total genes analyzed:", total_genes, "\n")
cat("Significant genes (FDR < 0.1):", sig_genes, "\n")
cat("Upregulated genes:", upregulated, "\n")
cat("Downregulated genes:", downregulated, "\n\n")

#############################
# Generate Figures
#############################

vsd <- vst(dds, blind = FALSE)

png(file.path(figures_dir, "pca_plot.png"), width = 1800, height = 1400, res = 300)
plotPCA(vsd, intgroup = "condition")
dev.off()

png(file.path(figures_dir, "ma_plot.png"), width = 1800, height = 1400, res = 300)
plotMA(resLFC, main = "MA Plot with Log2 Fold Change Shrinkage", ylim = c(-3, 3))
dev.off()

png(file.path(figures_dir, "dispersion_plot.png"), width = 1800, height = 1400, res = 300)
plotDispEsts(dds)
dev.off()

png(file.path(figures_dir, "pvalue_histogram.png"), width = 1800, height = 1400, res = 300)
hist(
  res$pvalue,
  breaks = 50,
  main = "Raw P-value Distribution",
  xlab = "p-value"
)
dev.off()

#############################
# Export Result Tables
#############################

resLFC_df <- as.data.frame(resLFC)

write.csv(
  resLFC_df,
  file = file.path(tables_dir, "deseq2_results.csv")
)

significant_genes <- resLFC_df[
  !is.na(resLFC_df$padj) & resLFC_df$padj < 0.1,
]

upregulated_genes <- resLFC_df[
  !is.na(resLFC_df$padj) &
    resLFC_df$padj < 0.1 &
    resLFC_df$log2FoldChange > 0,
]

downregulated_genes <- resLFC_df[
  !is.na(resLFC_df$padj) &
    resLFC_df$padj < 0.1 &
    resLFC_df$log2FoldChange < 0,
]

write.csv(
  significant_genes,
  file = file.path(tables_dir, "significant_genes.csv")
)

write.csv(
  upregulated_genes,
  file = file.path(tables_dir, "upregulated_genes.csv")
)

write.csv(
  downregulated_genes,
  file = file.path(tables_dir, "downregulated_genes.csv")
)

#############################
# Completion Message
#############################

cat("Analysis complete.\n")
cat("Figures saved to:", figures_dir, "\n")
cat("Tables saved to:", tables_dir, "\n")
