library(tximport)
library(readr)
library(DESeq2)
library(apeglm)

samples <- read.csv("samplesheet.csv")
samples$sample <- as.character(samples$sample)
files <- file.path("results_rnaseq/star_salmon", samples$sample, "quant.sf")
names(files) <- samples$sample
tx2gene <- read_tsv("results_rnaseq/star_salmon/tx2gene.tsv", col_names = c("TXNAME", "GENEID"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)
res <- results(dds)
resLFC <- lfcShrink(dds, coef=2, type="apeglm")

summary(resLFC)
sig_genes <- sum(resLFC$padj < 0.1, na.rm = TRUE)
upregulated <- sum(resLFC$padj < 0.1 & resLFC$log2FoldChange > 0, na.rm = TRUE)
downregulated <- sum(resLFC$padj < 0.1 & resLFC$log2FoldChange < 0, na.rm = TRUE)

head(resLFC[order(resLFC$padj), ], 10)

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

plotMA(resLFC, main="MA Plot (LFC Shrinkage)", ylim=c(-3,3))

plotDispEsts(dds)

hist(res$pvalue, breaks=50, main="Raw P-value Histogram", xlab="p-value")

write.csv(as.data.frame(resLFC), file = "deseq2_results.csv")
