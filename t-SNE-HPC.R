#!/usr/bin/env Rscript

# Shane Ridoux
# 250224
# Genomics tSNE pipeline for HPC

# Load required libraries
library(Seurat)
library(ggplot2)
library(Matrix)
library(future.apply)

# Set up parallel processing
plan("multicore", workers = 8)

# Increase memory limit
options(future.globals.maxSize = 10000 * 1024^2)

# Define data directory
data_dir <- "/home/ridoux/bio/project/GSE141784_RAW/"

# Load each time point separately
# 4-week dataset
data_4w <- ReadMtx(mtx = paste0(data_dir, "GSM4213199_NOD_4w_2849_matrix.mtx.gz"),
                   features = paste0(data_dir, "GSM4213199_NOD_4w_2849_genes.tsv.gz"),
                   cells = paste0(data_dir, "GSM4213199_NOD_4w_2849_barcodes.tsv.gz"))
seurat_4w <- CreateSeuratObject(counts = data_4w, project = "NOD_4w")
seurat_4w$time_point <- "4w"

# 8-week dataset
data_8w <- ReadMtx(mtx = paste0(data_dir, "GSM4213200_NOD_8w_2849_matrix.mtx.gz"),
                   features = paste0(data_dir, "GSM4213200_NOD_8w_2849_genes.tsv.gz"),
                   cells = paste0(data_dir, "GSM4213200_NOD_8w_2849_barcodes.tsv.gz"))
seurat_8w <- CreateSeuratObject(counts = data_8w, project = "NOD_8w")
seurat_8w$time_point <- "8w"

# 15-week dataset
data_15w <- ReadMtx(mtx = paste0(data_dir, "GSM4213198_NOD_15w_2734_matrix.mtx.gz"),
                    features = paste0(data_dir, "GSM4213198_NOD_15w_2734_genes.tsv.gz"),
                    cells = paste0(data_dir, "GSM4213198_NOD_15w_2734_barcodes.tsv.gz"))
seurat_15w <- CreateSeuratObject(counts = data_15w, project = "NOD_15w")
seurat_15w$time_point <- "15w"

# Merge datasets
seurat_combined <- merge(seurat_4w, y = c(seurat_8w, seurat_15w), add.cell.ids = c("4w", "8w", "15w"))

# Remove intermediate objects to free memory
rm(list = c("seurat_4w", "seurat_8w", "seurat_15w", "data_4w", "data_8w", "data_15w"))
gc()  # Run garbage collection to free RAM

########## QC & Normalization ##########
seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)

# Calculate ribosomal gene percentage
seurat_combined[["percent.ribo"]] <- PercentageFeatureSet(seurat_combined, pattern = "^Rp[sl]")

# Define variables to regress out
unwanted <- c("nCount_RNA", "percent.ribo")

# Remove cells with zero total counts
seurat_combined <- subset(seurat_combined, subset = unwanted[1] > 0)
seurat_combined <- subset(seurat_combined, subset = unwanted[2] > 0)

# Identify variable genes
# seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)

# Batch processing for ScaleData using parallel execution
batch_size <- 500
num_batches <- ceiling(nrow(seurat_combined) / batch_size)

# Function to process a batch
process_batch <- function(i) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, nrow(seurat_combined))
  batch_features <- rownames(seurat_combined)[start_idx:end_idx]
  
  message("Processing batch ", i, " of ", num_batches)
  
  # Scale only the batch of genes
  scaled_batch <- ScaleData(seurat_combined, vars.to.regress = unwanted, 
                            features = batch_features, do.return = TRUE)
  
  return(list(features = batch_features, scaled_data = scaled_batch@assays$RNA@scale.data[batch_features, ]))
}

# Run parallel processing
batch_results <- future_lapply(seq_len(num_batches), process_batch)

# Combine results
scaled_matrix <- do.call(rbind, lapply(batch_results, function(x) x$scaled_data))

# Assign scaled data back to Seurat object
seurat_combined@assays$RNA@scale.data <- scaled_matrix

# Save processed data
saveRDS(seurat_combined, "seurat_combined_scaled.rds")

# Shut down parallel workers & free memory
plan("sequential")
gc()

message("Seurat scaling complete on HPC!")