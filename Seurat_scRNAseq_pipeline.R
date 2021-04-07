#Initialize and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(patchwork)
library(ggplot2)

#Load GBM data downloaded from 10X genomics
pbmc.data <-Read10X ('GBM/GBM_10X_Genomics/Parent_SC3v3_Human_Glioblastoma/filtered_feature_bc_matrix/')

#Examine the memory savings between regular and sparse matrices
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size

# Initialize the Seurat object with the raw (non-normalized data)
pbmc <- new("seurat", raw.data = pbmc.data)

#Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Calculate the proportion of transcripts mapping to mitochondrial genes as a QC metric
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2)) #CombinePlots is deprecated, hence usign patchwork library below

#Filter cells that have unique feature counts over 2,500 or less than 200 and have >5% mitochondrial counts

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge=0, ynudge=0)
plot1 + plot2 + plot_layout(ncol = 1)

#scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Perform linear dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results in a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Visualize the single-cell clusters either using UMAP or tSNE
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne")

saveRDS(pbmc, file = "GBM/GBM_test.rds")

# Find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# Find the markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>% group_by(cluster) %>% top_n(2)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#violin plots shows expression probability distribution across clusters
VlnPlot(pbmc, features = c("CCL3", "APOE"))

# Violin Plot using raw UMI counts
VlnPlot(pbmc, features = c("ANXA1", "VEGFA"), slot = "counts", log = TRUE)

#Visualize individul tSNE plots per feature using FeaturePlots
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD14", "A2M", "FCGR3A", "LYZ", "CCR7", "APOE", "MBP"))

#Adding labels to the tSNE/UMAP clusters
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "GBM/GBM_final.rds")

pbmc[["ClusterNames_0.6"]] <- Idents(object = pbmc)
pbmc <- FindClusters(pbmc,resolution = 0.8)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = FALSE)

#Find T-cell markers and viualize using feature plot
tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)
FeaturePlot(object = pbmc, features = c("S100A4", "CCR7"), cols= c("green", "red"))

