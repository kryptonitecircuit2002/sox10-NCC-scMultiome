lapply(required_packages, library, character.only = TRUE)
set.seed(1234)

H3.data <- readRDS("..../sox10ncc_h33_kd.rds")
Control <- readRDS("..../sox10_ncc_linked.rds")
Control <- Control[, sample(colnames(Control), size =1406, replace=F)]

# PCA, clustering, and UMAP visualization
DefaultAssay(H3.data) <- "SCT"
H3.data <- RunPCA(H3.data, verbose = FALSE, reduction.name = "pca.h3")
H3.data <- FindNeighbors(H3.data, reduction = "pca.h3", dims = 1:30)
H3.data <- FindClusters(H3.data, resolution = 0.7)
H3.data <- RunUMAP(H3.data, reduction = "pca.h3", reduction.name = "umap.rna.h3", dims = 1:30, verbose = FALSE, spread = 0.32, min.dist = 0.35)
DimPlot(H3.data, reduction = "umap.rna.h3", label = TRUE)

##Merge using RNA anchors
Merge_H3.3 <- merge(Control, y = H3.data, add.cell.ids = c("Control", "H3.3 KD"))
DefaultAssay(Merge_H3.3) <- "peaks"
Merge_H3.3 <- NucleosomeSignal(Merge_H3.3)
Merge_H3.3  <- TSSEnrichment(Merge_H3.3, fast = FALSE)

DensityScatter(Merge_H3.3, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
VlnPlot(
  object = Merge_H3.3,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  ncol = 4,
  pt.size = 0)

#Filter out low quality cells
Merge_H3.3  <- subset(
  x = Merge_H3.3,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 40000 &
    nucleosome_signal < 1.75 &
    TSS.enrichment < 10
) 

# Perform SCTransform normalization on RNA assay
DefaultAssay(Merge_H3.3) <- "RNA"
Merge_H3.3 <- SCTransform(Merge_H3.3, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

# Integration step
H3.list <- SplitObject(Merge_H3.3, split.by = "orig.ident")
Control_H3 <- H3.list[["Control"]]
H3.3 <- H3.list[["H3.3 KD"]]

features.h3 <- SelectIntegrationFeatures(object.list = H3.list, nfeatures = 2000)
H3.list <- PrepSCTIntegration(object.list = H3.list, anchor.features = features.h3)
H3.anchors <- FindIntegrationAnchors(object.list = H3.list, normalization.method = "SCT", anchor.features = features.h3)
H3.combined.sct <- IntegrateData(anchorset = H3.anchors, normalization.method = "SCT")

# PCA, clustering, and UMAP visualization
H3.combined.sct <- RunPCA(H3.combined.sct, verbose = FALSE, reduction.name = "pca.h3")
H3.combined.sct <- FindNeighbors(H3.combined.sct, reduction = "pca.h3", dims = 1:30)
H3.combined.sct <- FindClusters(H3.combined.sct, resolution = 0.5)
H3.combined.sct <- RunUMAP(H3.combined.sct, reduction = "pca.h3", reduction.name = "umap.rna.h3", dims = 1:30, verbose = FALSE,  spread = 0.15, min.dist = 0.25)
DimPlot(H3.combined.sct, reduction = "umap.rna.h3", split.by = "orig.ident", label = TRUE)

saveRDS(H3.combined.sct, file = "..../Integrated_H3.rds")
