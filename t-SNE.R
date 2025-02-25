# Shane Ridoux
# 250224
# Genomics tSNE pipeline

rm(list=ls())
cat("\014")
library(Seurat)
library(ggplot2)
library(Seurat)
library(Matrix)
library(future)
library(RCurl)

############# SINGLE-CELL CLUSTERING USING SEURAT ###################
######### Load data
data_dir <- "/Users/shane/School/CU-Denver/24-Fall/Genomics-Bioinformatics/project/GSE141784_RAW/"

data_4w <- ReadMtx(mtx = paste0(data_dir, "GSM4213199_NOD_4w_2849_matrix.mtx.gz"),
                   features = paste0(data_dir, "GSM4213199_NOD_4w_2849_genes.tsv.gz"),
                   cells = paste0(data_dir, "GSM4213199_NOD_4w_2849_barcodes.tsv.gz"))
data_8w <- ReadMtx(mtx = paste0(data_dir, "GSM4213200_NOD_8w_2849_matrix.mtx.gz"),
                   features = paste0(data_dir, "GSM4213200_NOD_8w_2849_genes.tsv.gz"),
                   cells = paste0(data_dir, "GSM4213200_NOD_8w_2849_barcodes.tsv.gz"))
data_15w <- ReadMtx(mtx = paste0(data_dir, "GSM4213198_NOD_15w_2734_matrix.mtx.gz"),
                    features = paste0(data_dir, "GSM4213198_NOD_15w_2734_genes.tsv.gz"),
                    cells = paste0(data_dir, "GSM4213198_NOD_15w_2734_barcodes.tsv.gz"))

seurat_list <- list(
  "4w" = CreateSeuratObject(counts = data_4w, project = "NOD_4w"),
  "8w" = CreateSeuratObject(counts = data_8w, project = "NOD_8w"),
  "15w" = CreateSeuratObject(counts = data_15w, project = "NOD_15w")
)
rm(list = c("data_15w","data_8w","data_4w"))

# LogNormalize() data
seurat_list <- lapply(seurat_list, NormalizeData)

# Add ribosomal percentage
seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  return(obj)
})

# Find variable features
seurat_list <- lapply(seurat_list, FindVariableFeatures, selection.method = "vst", nfeatures = 2000)

# Filter based on violin plots
seurat_list <- lapply(seurat_list, function(obj) {
  subset(obj, subset = nFeature_RNA > 500 & 
           nFeature_RNA < 2500 & 
           nCount_RNA < 10000 & 
           percent.ribo < 30)
})

# Subset to just leukocytes like in the paper
leukocyte_list <- lapply(seurat_list, function(obj) {
  subset(obj, subset = Ptprc > 0)
})

######## RunCCA --> FindIntegrationAnchors
anchors <- FindIntegrationAnchors(object.list = leukocyte_list, anchor.features = 2000)
seurat_combined <- IntegrateData(anchorset = anchors)

# Default assay should be the integrated one now
DefaultAssay(seurat_combined) <- "integrated"

########## Scale PCA and tSNE
seurat_combined <- ScaleData(seurat_combined)
seurat_combined <- RunPCA(seurat_combined, npcs = 30)
ElbowPlot(seurat_combined) # I picked 15 cuz it seemed to flatten out there
seurat_combined <- RunTSNE(seurat_combined, dims = 1:15)

# Find clusters
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:15)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.8)

# Factor time points for ordering
seurat_combined$time_point <- factor(seurat_combined$orig.ident,
                                     levels = c("NOD_4w", "NOD_8w", "NOD_15w"),
                                     labels = c("4w","8w","15w"))

FeaturePlot(seurat_combined, 
            features = c("Cd79a", "Cd3e", "C1qa", "Emr1", 
                         "H2-Oa", "Flt3", "H2-Ab1", "Siglech", "Igj"), 
            cols = c("lightgray", "red")) & NoAxes()

DotPlot(seurat_combined, features = c("Cd79a", "Cd3e", "C1qa", "Emr1", 
                                      "H2-Oa", "Flt3", "H2-Ab1", "Siglech", "Igj")) +
  RotatedAxis()

# Make labels based on DotPlot
new_cluster_ids <- c(
  "Mac",    # Cluster 0 
  "T cell",    # Cluster 1 
  "B cell",       # Cluster 2 
  "cDC",       # Cluster 3 
  "Mac",       # Cluster 4 
  "Mac", # Cluster 5 
  "T cell", # Cluster 6 
  "cDC",       # Cluster 7 
  "T cell",    # Cluster 8 
  "B cell",    # Cluster 9 
  "Other Lymphoid", # Cluster 10 
  "T cell", # Cluster 11 
  "cDC",       # Cluster 12 
  "Mac",       # Cluster 13 
  "pDC",       # Cluster 14
  "cDC",    # Cluster 15 
  "Mac",    # Cluster 16
  "Mac",       # Cluster 17 
  "Other Lymphoid",       # Cluster 18 
  "Other Lymphoid",       # Cluster 19 
  "Plasma cell"  # Cluster 20 
)

names(new_cluster_ids) <- levels(seurat_combined)
new_cluster_ids <- factor(new_cluster_ids,
                          levels = c("B cell", "T cell", "Mac", "cDC", "pDC", "Plasma cell", "Other Lymphoid"))

# Rename clusters
seurat_combined_renamed <- RenameIdents(seurat_combined, new_cluster_ids)

# Set colors similar to paper
custom_colors <- c("B cell" = "#E41A1C", 
                   "T cell" = "pink", 
                   "Mac" = "#4DAF4A", 
                   "cDC" = "#FF7F00",
                   "pDC" = "#377EB8", 
                   "Plasma cell" = "#6A3D9A",
                   "Other Lymphoid" = "lightgreen")
                   


DimPlot(seurat_combined_renamed, reduction = "tsne", label = FALSE, pt.size = 0.5) +
  scale_color_manual(values = custom_colors) +
  ggtitle("Figure 1C") & NoAxes()

DotPlot(seurat_combined_renamed, features = c("Cd79a", "Cd3e", "C1qa", "Emr1", 
                                              "H2-Oa", "Flt3", "H2-Ab1", "Siglech", "Igj")) +
  RotatedAxis()

DimPlot(seurat_combined_renamed, reduction = "tsne", label = FALSE, pt.size = 0.5, split.by = "time_point") +
  scale_color_manual(values = custom_colors) +
  ggtitle("Figure 1E") & NoAxes()

plot <- FeaturePlot(seurat_combined_renamed,  
                    features = c("Cd79a", "Cd3e", "C1qa", "Emr1", 
                                 "H2-Oa", "Flt3", "H2-Ab1", "Siglech", "Igj"),  
                    cols = c("lightgray", "red")) & NoAxes()
plot + patchwork::plot_annotation(title = "Figure 1D")
