lapply(required_packages, library, character.only = TRUE)
set.seed(1234)

H2.data <- readRDS(".../sox10ncc_h2az_kd.rds")
Control <- readRDS(".../sox10_ncc_linked.rds")
Control <- Control[, sample(colnames(Control), size =667, replace=F)]

##Merging using RNA anchors
DefaultAssay(Control) <- "SCT"
DefaultAssay(H2.data) <- "SCT"
Merge_H2 <- merge(Control, y = H2.data, add.cell.ids = c("Control", "H2A.Z KD"))
DefaultAssay(Merge_H2) <- "ATAC"
Merge_H2 <- NucleosomeSignal(Merge_H2)
Merge_H2  <- TSSEnrichment(Merge_H2)

DensityScatter(Merge_H2, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
VlnPlot(
  object = Merge_H2,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  ncol = 4,
  pt.size = 0)

#Filtering
Merge_H2  <- subset(
  x = Merge_H2,
  subset = nCount_ATAC < 60000 &
    nCount_RNA < 40000 &
    nucleosome_signal < 1.75 &
    TSS.enrichment < 10
)
# Perform SCTransform normalization 
DefaultAssay(Merge_H2) <- "RNA"
Merge_H2 <- SCTransform(Merge_H2, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

# Integratio
H2.list <- SplitObject(Merge_H2, split.by = "orig.ident")
Control_H2 <- H2.list[["Control"]]
H2 <- H2.list[["H2A.Z KD"]]

features.h2 <- SelectIntegrationFeatures(object.list = H2.list, nfeatures = 1500)
H2.list <- PrepSCTIntegration(object.list = H2.list, anchor.features = features.h2)
H2.anchors <- FindIntegrationAnchors(object.list = H2.list, normalization.method = "SCT", anchor.features = features.h2)
H2.combined <- IntegrateData(anchorset = H2.anchors, normalization.method = "SCT")

# PCA, clustering, and UMAP visualization
H2.combined <- RunPCA(H2.combined, verbose = FALSE, reduction.name = "pca.h2")
H2.combined <- FindNeighbors(H2.combined, reduction = "pca.h2", dims = 1:30)
H2.combined <- FindClusters(H2.combined, resolution = 0.5)
H2.combined <- RunUMAP(H2.combined, reduction = "pca.h2", reduction.name = "umap.rna.h2", dims = 1:30, verbose = FALSE,  spread = 0.25, min.dist = 0.35)
DimPlot(H2.combined, reduction = "umap.rna.h2", split.by = "orig.ident", label = TRUE)

saveRDS(H2.combined, file = "Integrated_H2.rds")
