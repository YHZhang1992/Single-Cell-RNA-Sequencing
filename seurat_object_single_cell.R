library(Seurat)
library(dplyr)
library(patchwork)

data_dir <- "path/to/filtered_feature_bc_matrix"
sc_data <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "10X_ScRNA", min.cells = 3, min.features = 200)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)  # Unsupervised clustering
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)  # Visualize clusters

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = top_markers$gene) + NoLegend()

