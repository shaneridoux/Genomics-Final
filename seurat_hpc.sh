#!/bin/bash
#SBATCH --job-name=seurat_hpc
#SBATCH --output=seurat_output.log
#SBATCH --error=seurat_error.log
#SBATCH --time=24:00:00
#SBATCH --mem=64G 
#SBATCH --cpus-per-task=8 
#SBATCH --ntasks=1
#SBATCH --partition=compute

# Load R module
module load R/4.2.0

# Run the R script
Rscript t-SNE-HPC.R
