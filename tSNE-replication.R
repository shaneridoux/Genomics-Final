rm(list=ls())
cat("\014")
library(Seurat)
library(ggplot2)
library(Seurat)
library(Matrix)

# Load the dataset
data_dir <- "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/nodrag"
sc_data <- Read10X(data.dir = data_dir)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "scRNAseq")

# View metadata summary
print(seurat_obj)
head(seurat_obj@meta.data)

# Make counts_matrix
counts_matrix <- GetAssayData(seurat_obj, layer = "counts")
dim(counts_matrix)

# Load each time point separately
# 4-week dataset
data_4w <- ReadMtx(mtx = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213199_NOD_4w_2849_matrix.mtx.gz",
                   features = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213199_NOD_4w_2849_genes.tsv.gz",
                   cells = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213199_NOD_4w_2849_barcodes.tsv.gz")
seurat_4w <- CreateSeuratObject(counts = data_4w, project = "NOD_4w")
seurat_4w$time_point <- "4w"

# 8-week dataset
data_8w <- ReadMtx(mtx = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213200_NOD_8w_2849_matrix.mtx.gz",
                   features = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213200_NOD_8w_2849_genes.tsv.gz",
                   cells = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213200_NOD_8w_2849_barcodes.tsv.gz")
seurat_8w <- CreateSeuratObject(counts = data_8w, project = "NOD_8w")
seurat_8w$time_point <- "8w"

# 15-week dataset
data_15w <- ReadMtx(mtx = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213198_NOD_15w_2734_matrix.mtx.gz",
                    features = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213198_NOD_15w_2734_genes.tsv.gz",
                    cells = "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/GSM4213198_NOD_15w_2734_barcodes.tsv.gz")
seurat_15w <- CreateSeuratObject(counts = data_15w, project = "NOD_15w")
seurat_15w$time_point <- "15w"

# Merge datasets with unique cell identifiers for each time point
seurat_combined <- merge(seurat_4w, y = c(seurat_8w, seurat_15w), 
                         add.cell.ids = c("4w", "8w", "15w"))
rm(list = c("data_15w", "data_8w", "data_4w"))
rm(list = c("seurat_15w", "seurat_8w", "seurat_4w"))

# Check time points
table(seurat_combined$time_point)

# Add mitochondrial gene percentage
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^mt-")

# Visualize quality control metrics
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter
seurat_combined <- subset(seurat_combined, 
                          subset = nFeature_RNA > 500 & 
                            nFeature_RNA < 2500 & 
                            nCount_RNA < 20000 & 
                            percent.mt < 10)

# Check the number of cells remaining after filtering
table(seurat_combined$time_point)

# Re-plot QC metrics to confirm filtering
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalize the data and find variable features
seurat_combined <- NormalizeData(seurat_combined)
seurat_combined <- FindVariableFeatures(seurat_combined)

# Scale the data and run PCA
seurat_combined <- ScaleData(seurat_combined)
seurat_combined <- RunPCA(seurat_combined, npcs = 20)

# Visualize PCA results
ElbowPlot(seurat_combined)

# Find neighbors and clusters
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:7)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)  # Adjust resolution as needed

# Check cluster levels again
levels(seurat_combined$seurat_clusters)

# Assign new cell type labels based on identified markers
new_cluster_ids <- c("B cells",         # Cluster 0
                     "T cells",         # Cluster 1
                     "Mac",             # Cluster 2 (Macrophages)
                     "cDC",             # Cluster 3 (Conventional Dendritic Cells)
                     "pDC",             # Cluster 4 (Plasmacytoid Dendritic Cells)
                     "Plasma cells",    # Cluster 5
                     "Other Lymphoid",  # Cluster 6 (NK cells, etc.)
                     "Unknown",         # Remaining unknown clusters
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown")

# Apply labels to clusters
names(new_cluster_ids) <- levels(seurat_combined$seurat_clusters)
seurat_combined <- RenameIdents(seurat_combined, new_cluster_ids)

# Confirm updated cluster assignments
table(Idents(seurat_combined))

# Run t-SNE using the top principal components (adjust dims if needed)
seurat_combined <- RunTSNE(seurat_combined, dims = 1:7)  # Use 7 PCs based on your ElbowPlot

# Confirm t-SNE has run
DimPlot(seurat_combined, reduction = "tsne", label = TRUE) +
  ggtitle("t-SNE Clusters")

#Define custom colors to match the paper's scheme
custom_colors <- c("B cells" = "#E41A1C",         # Red
                   "T cells" = "#984EA3",         # Purple
                   "Mac" = "#4DAF4A",             # Green
                   "cDC" = "#FF7F00",             # Orange
                   "pDC" = "#377EB8",             # Blue
                   "Plasma cells" = "#6A3D9A",    # Violet
                   "Other Lymphoid" = "#00CED1",  # Cyan
                   "Unknown" = "gray")            # Gray for unknown

# Plot final t-SNE with updated labels and colors
DimPlot(seurat_combined, reduction = "tsne", split.by = "time_point", group.by = "ident", 
        cols = custom_colors, ncol = 3, label = FALSE) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = "dashed")) +
  ggtitle("t-SNE of Immune Cell Heterogeneity Across Time Points")
# 
# # Add dashed outlines around clusters
# DimPlot(seurat_combined, reduction = "tsne", split.by = "time_point", group.by = "ident",
#         cols = custom_colors, ncol = 3, label = TRUE) +
#   theme_minimal() +
#   theme(panel.border = element_rect(color = "black", fill = NA, linetype = "dashed")) +
#   ggtitle("t-SNE of Immune Cell Heterogeneity Across Time Points")
# 
# DimPlot(seurat_combined, reduction = "tsne", split.by = "time_point", group.by = "ident",
#         cols = custom_colors, ncol = 3, label = TRUE, label.size = 6, repel = TRUE) +
#   theme_minimal(base_size = 14) +  # Increase overall font size
#   ggtitle("t-SNE of Immune Cell Heterogeneity Across Time Points")

# Set time_point as a factor to control the plot order
seurat_combined$time_point <- factor(seurat_combined$time_point, levels = c("4w", "8w", "15w"))

# Re-plot with the correct order
DimPlot(seurat_combined, reduction = "tsne", split.by = "time_point", group.by = "ident",
        cols = custom_colors, ncol = 3, label = TRUE, repel = TRUE) +
  theme_minimal() +
  ggtitle("t-SNE of Immune Cell Heterogeneity Across Time Points")

# Plasma Cell Markers
FeaturePlot(seurat_combined, features = c("Prdm1", "Xbp1", "Sdc1"), min.cutoff = "q9")

# Plasmacytoid Dendritic Cell (pDC) Markers
FeaturePlot(seurat_combined, features = c("Siglech", "Bst2", "Irf7"), min.cutoff = "q9")

# Other Lymphoid Markers (if relevant)
FeaturePlot(seurat_combined, features = c("Nkg7", "Gzmb"), min.cutoff = "q9")

# Re-run t-SNE with optimized perplexity
seurat_combined <- RunTSNE(seurat_combined, dims = 1:7, perplexity = 30, max_iter = 1000)

# Re-plot with the optimized layout
DimPlot(seurat_combined, reduction = "tsne", split.by = "time_point", group.by = "ident", 
        cols = custom_colors, ncol = 3, label = FALSE) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = "dashed")) +
  ggtitle("t-SNE of Immune Cell Heterogeneity Across Time Points")
