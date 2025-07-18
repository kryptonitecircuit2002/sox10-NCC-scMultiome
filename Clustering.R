library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(future)
library(patchwork)
library(hdf5r)
library(readr)
library(pheatmap)
library(ggrepel)
library(LSD)
library(MASS)
library(ensembldb)
library(BSgenome.Drerio.UCSC.danRer11)

set.seed(1234)

#Read R object
data <- readRDS(".../data_linked.rds")

# PCA, clustering, and UMAP visualization
DefaultAssay(data) <- "SCT"
data <- RunPCA(data, verbose = FALSE, reduction.name = "pca")
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.7)
data <- RunUMAP(data, reduction = "pca", reduction.name = "umap.rna", dims = 1:30, verbose = FALSE, spread = 0.32, min.dist = 0.35)
DimPlot(data, reduction = "umap.rna", label = TRUE)

#ATAC Assay
DefaultAssay(data) <- "ATAC"
data <- RunTFIDF(data)
data <- FindTopFeatures(data)
data <- RunSVD(data, reduction.name = "lsi.atac")
data <- FindNeighbors(data, reduction = "lsi.atac", dims = 2:10, assay = "ATAC")
data <- FindClusters(data, resolution = 0.7)
data <- RunUMAP(data, reduction = 'lsi.atac', dims = 2:10, assay = 'ATAC',
                        reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(data, reduction = "umap.atac", label = T)

#Peaks Assay
DefaultAssay(data) <- "peaks"
data <- RunTFIDF(data)
data <- FindTopFeatures(data)
data <- RunSVD(data, reduction.name = "lsi.peaks")
data <- FindNeighbors(data, reduction = "lsi.peaks", dims = 2:10, assay = "ATAC")
data <- FindClusters(data, resolution = 0.5)
data <- RunUMAP(data, reduction = 'lsi.peaks', dims = 2:10, assay = 'peaks',
                        reduction.name = "umap.peaks", reduction.key = "peakUMAP_", spread = 0.3)
DimPlot(data, reduction = "umap.peaks", label = T)

saveRDS(".../data_clustered.rds")
