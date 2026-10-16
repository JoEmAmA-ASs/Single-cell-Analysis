library(dplyr)
library(Seurat)
library(patchwork)


file_path <- "D:/GITHUB_MAXXING/Single Cell Studies (FUCK YOU VEENA)/GSE77288_reads-raw-single-per-sample.txt.gz"

df <- read.delim(gzfile(file_path), header = TRUE, sep = "\t")

# METADATA EXTRACTION
meta_data <- df[, c("individual", "replicate", "well"), drop = FALSE]
meta_data$cell_id <- paste0(meta_data$individual, "_", meta_data$replicate, "_", meta_data$well)
rownames(meta_data) <- make.unique(meta_data$cell_id)

# EXPRESSION MATRIX
expr_matrix <- as.matrix(df[, 4:ncol(df), drop = FALSE])
# Rows are cells/wells; columns are genes. Seurat expects genes x cells.
expr_matrix_t <- t(expr_matrix)
colnames(expr_matrix_t) <- rownames(meta_data)
rownames(expr_matrix_t) <- colnames(df)[4:ncol(df)]

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = expr_matrix_t,
  meta.data = meta_data,
  min.cells = 3,       # Keep genes expressed in at least 3 wells
  min.features = 500   # Keep wells with at least 500 detected genes
)

print(seurat_obj)

#QUALITY CONTROL

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) #QC METRIC PLOT
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 800 & nCount_RNA > 1500) #FILTER LOW QUALITY CELLS

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000) #NORMALIZATION
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000) #HVG SELECTION
seurat_obj <- ScaleData(seurat_obj) #DATA SCALING

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj)) #PCA
ElbowPlot(seurat_obj, ndims  = 5) #APPROPRIATE NUMBER OF PC SELECTION

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) #CLUSTERING USING LOUVAIN ALGORITHM

seurat_obj <- RunUMAP(seurat_obj, dims = 1:20) #UMAP DIMENSIONALITY REDUCTION

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "individual")
print(p1 + p2)

cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #FINDING MARKER GENES
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) #TOP 5 MARKERS PER CLUSTER
print(top_markers)