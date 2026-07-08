#!/bin/bash

#SBATCH --job-name=rnaseq

#SBATCH --account=pr_284_general

#SBATCH --output=rnaseq_%j.out

#SBATCH --error=rnaseq_%j.err

#SBATCH --time=24:00:00

#SBATCH --cpus-per-task=12

#SBATCH --mem=96G

module load nextflow/24.10.4

nextflow run nf-core/rnaseq -r 3.14.0 -profile singularity --input samplesheet.csv --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf Homo_sapiens.GRCh38.111.gtf --outdir results_rnaseq --saveReference --salmon_gc_bias
