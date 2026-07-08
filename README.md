# RNA-seq Differential Expression Analysis Using DESeq2
![R](https://img.shields.io/badge/R-4.x-blue)
![DESeq2](https://img.shields.io/badge/DESeq2-Bioconductor-green)
![Nextflow](https://img.shields.io/badge/Nextflow-Workflow-orange)
![License](https://img.shields.io/badge/License-MIT-yellow)

Analysis of the publicly available **PRJNA490376** RNA-seq dataset comparing **NRDE2 knockdown** and **control** samples using the **nf-core/rnaseq** workflow, **Salmon**, **tximport**, and **DESeq2**.

This project demonstrates an end-to-end RNA-seq analysis workflow, from transcript quantification to differential gene expression analysis and downstream visualization.
## Overview

This repository presents a complete RNA-seq differential expression workflow for identifying genes affected by NRDE2 knockdown. The analysis begins with Salmon transcript quantification, summarizes transcript abundance to gene-level counts using tximport, and identifies differentially expressed genes with DESeq2.
## Objective

The goal of this project is to identify genes that are differentially expressed between NRDE2 knockdown and control samples using an RNA-seq differential expression workflow.

## Tools Used

**Programming**
- R

**RNA-seq Analysis**
- Salmon
- tximport
- DESeq2
- apeglm

**Visualization**
- ggplot2

**Data Import**
- readr

**Workflow**
- nf-core/rnaseq
- Nextflow
- SLURM

## Repository Structure

```text
scripts/      R analysis scripts
pipeline/     Nextflow and SLURM execution files
docs/         Analysis report and workflow reports
results/      Figures and output tables
data/         Dataset information
environment/  Session information and package versions
```
# Data

This repository documents the analysis performed on the publicly available **PRJNA490376** RNA-seq dataset.

Due to storage limitations and because some intermediate files are no longer available, the repository does **not** include:

- Raw FASTQ files
- Salmon transcript quantification files (`quant.sf`)
- Transcript-to-gene mapping (`tx2gene.tsv`)

The repository focuses on the downstream differential expression analysis and includes:

- DESeq2 analysis script
- Analysis report
- Workflow documentation
- SLURM submission script
- Sample metadata
  
## Input

The analysis was performed using:

- Salmon transcript quantification files (`quant.sf`)
- Transcript-to-gene mapping (`tx2gene.tsv`)
- Sample metadata (`samplesheet.csv`)

Only the sample metadata is included in this repository. The original quantification and mapping files are not included.
Note: The original Salmon quantification files and transcript-to-gene mapping file are not included, so the script documents the analysis workflow but cannot be rerun from this repository alone without regenerating those files.
## Workflow

1. Import Salmon transcript quantifications.
2. Summarize transcript-level abundance to gene-level counts using tximport.
3. Perform differential expression analysis with DESeq2.
4. Apply log2 fold change shrinkage using apeglm.
5. Generate quality-control visualizations.
6. Export differential expression results.

## Results

Running the workflow generates:

- Differential expression results (`deseq2_results.csv`)
- Significant genes table (`significant_genes.csv`)
- Upregulated genes table (`upregulated_genes.csv`)
- Downregulated genes table (`downregulated_genes.csv`)
- PCA plot
- MA plot
- Dispersion estimate plot
- P-value histogram

## Skills Demonstrated

- RNA-seq differential expression analysis
- Transcript quantification with Salmon
- Gene-level count summarization using tximport
- Differential expression analysis with DESeq2
- Log2 fold change shrinkage with apeglm
- Statistical analysis of sequencing data
- RNA-seq quality assessment
- Data visualization in R
- Reproducible bioinformatics workflow development

## Acknowledgments

This project uses the publicly available RNA-seq dataset **PRJNA490376** from the NCBI Sequence Read Archive (SRA).
