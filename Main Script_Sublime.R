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

setwd("....../Control+H3.3")

#Read the combined count matrix
inputdata.10X.Control <- Read10X_h5(".../filtered_feature_bc_matrix_control.h5")
inputdata.10X.H33 <- Read10X_h5("..../filtered_feature_bc_matrix_H33.h5")

#path to atac fragment file
fragpath.Control <- "...../atac_fragments_control.tsv.gz"
fragpath.H33 <- "..../atac_fragments_H33.tsv.gz"

#create counts file for Control and H3.3
rna.control <- inputdata.10X.Control$`Gene Expression`
atac.control <- inputdata.10X.Control$Peaks
rna_counts.H33 <- inputdata.10X.H33$`Gene Expression`
atac_counts.H33 <- inputdata.10X.H33$Peaks

#feed metadata and create seurat object
meta_data.control <- read.csv(
  file = ".../per_barcode_metrics_control.csv",
  header = TRUE,
  row.names = 1
)

meta_data.H33 <- read.csv(
  file = "..../per_barcode_metrics_H33.csv",
  header = TRUE,
  row.names = 1
)

#Create gene annotaion from GTF
gref.path = "D:/Transit/Single_Cell+ATAC/Control/ATAC/Data+Analysis/Analysis_Control/Danio_rerio.GRCz11.105.gtf"
gtf_zf <- rtracklayer::import(gref.path)
gene.coords.zf <- gtf_zf
gene.coords.zf <- gene.coords.zf[! is.na(gene.coords.zf$gene_name),]
gene.coords.zf <- keepStandardChromosomes(gene.coords.zf, pruning.mode = 'coarse')
genome(gene.coords.zf) <- 'GRCz11'
# copy the "gene_id" for the "tx_id" and "transcript_id" 
gene.coords.zf$tx_id <- gene.coords.zf$gene_id
gene.coords.zf$transcript_id <- gene.coords.zf$gene_id

#Create Seurat objects for both
Control.data <- CreateSeuratObject(counts = rna_counts.control,  meta.data = meta_data.control, project = "Control")
chrom_assay_control <- CreateChromatinAssay(
  counts = atac_counts.control,
  sep = c(":", "-"),
  genome = 'GRCz11',
  fragments = fragpath.Control,
  annotation = gene.coords.zf,
  min.cells = 10,
  min.features = 1
)
Control.data[["ATAC"]]<- chrom_assay_control

Seurat.H33 <- CreateSeuratObject(counts = rna_counts.H33, project = "H3.3", meta.data = meta_data.H33)
chrom_assay_H33 <- CreateChromatinAssay(
  counts = atac_counts.H33,
  sep = c(":", "-"),
  genome = 'GRCz11',
  fragments = fragpath.H33,
  annotation = gene.coords.zf,
  min.cells = 10,
  min.features = 1
)
Seurat.H33[["ATAC"]] <- chrom_assay_H33

#since the number of cells in sample was less than control, we downsampled at the very beginning and started with equal number of cells across control and condition
#perform downsampling
Control.H3 <- Control.data[, sample(colnames(Control.data), size =1346, replace=F)]

#Merge the two seurat objects
Merge_H3.3 <- merge(x = Control.H3, y = Seurat.H33, add.cell.ids = c("Control", "H3.3"))

#QC and preprocessing
DefaultAssay(Merge_H3.3) <- "ATAC"
Merge_H3.3  <- NucleosomeSignal(Merge_H3.3)
Merge_H3.3  <- TSSEnrichment(Merge_H3.3)
DensityScatter(Merge_H3.3, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
VlnPlot(
  object = Merge_H3.3,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

#Filter out low quality cells
Merge_H3.3  <- subset(
  x = Merge_H3.3,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 12000 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 1 &
    TSS.enrichment < 10
)
 
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#For integration, we first perform normalisation across individual RNA layers, find integration anchors, integrate the datasets and then 
#Perform pre-analysis
#FOR RNA
#SCT
DefaultAssay(Merge_H3.3) <- "RNA"
H3.3_list <- SplitObject(Merge_H3.3, split.by = "orig.ident")
Control <- H3.3_list[["Control"]]
H3.3 <- H3.3_list[["H3.3"]]

Control <- SCTransform(Control, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

H3.3 <- SCTransform(H3.3, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)


#Create integration lists 
H3.3_list <- list(Control = Control, H3.3 = H3.3)
features <- SelectIntegrationFeatures(object.list = H3.3_list, nfeatures = 3000)
H3.3_list <- PrepSCTIntegration(object.list = H3.3_list, anchor.features = features)

#Find anchors
H3.3_rna.anchors <- FindIntegrationAnchors(object.list = H3.3_list, normalization.method = "SCT",
                                     anchor.features = features)
H3.3.combined <- IntegrateData(anchorset = H3.3_rna.anchors, normalization.method = "SCT")

#perform Analysis
H3.3.combined <- RunPCA(H3.3.combined, dims = 1:30, verbose = FALSE)
H3.3.combined <- FindNeighbors(H3.3.combined, reduction = "pca", dims = 1:30)
H3.3.combined <- FindClusters(H3.3.combined, resolution = 0.3)
H3.3.combined <- RunUMAP(H3.3.combined, reduction = "pca", reduction.name = "umap.rna", dims = 1:30, verbose = FALSE,  spread = 0.1, min.dist = 0.15)
DimPlot(H3.3.combined, reduction = "umap.rna", split.by = "orig.ident", label = T, pt.size = 1.1)

DotPlot(H3.3.combined, features = c("twist1a", "grem2b", "col11a1b", "lamc3", "mitfa", "aox5", "tfec", "paics", "gfap", "pou3f1", "crestin", "sox10", "oc90", "epcam", "elavl3", "elavl4", "foxd3", "zeb2a", "sox9b", "dct", "kita", "tpma", "musk", "ltk", "alx4a")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="viridis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() + 
  coord_flip()

RNA_cluster.ids <- c("0" = "mesenchymal", "1" = "twist1a+ progenitor", "2" = "neural/glial", "3" = "MIX+", "4" = "glial/neural",
                     "5" = "foxd3 progenitor2", "6" = "otic", "7" = "foxd3 progenitor 1", "8" = "neuronal", "9" = "melanoblast",
                     "10" = "skeletal muscle progenitor", "11" = "twist1a+ progenitor", "12" = "iridoblast", "13" = "skeletal muscle progenitor",
                     "14" = "neuronal", "15" = "otic", "16" = "neural/glial")

New_cluster.ids <- c("0" = "twist1a+ progenitor", "1" = "twist1a+ progenitor", "2" = "mesenchymal",
                     "3" = "neural/glial", "4" = "otic", "5" = "foxd3 progenitor 1", "6" = "MIX+", "7" = "neuronal", "8" = "foxd3 progenitor2",
                    "9" = "skeletal muscle progenitor", "10" = "melanoblast",  "11" = "twist1a+ progenitor", "12" = "iridoblast")
H3.3.combined <- RenameIdents(H3.3.combined, New_cluster.ids)
#for ATAC
#first I run lsi reduction for combined object
DefaultAssay(H3.3.combined) <- "ATAC"
H3.3.combined<- RunTFIDF(H3.3.combined)
H3.3.combined <- FindTopFeatures(H3.3.combined)
H3.3.combined <- RunSVD(H3.3.combined)

H3.3.combined <- FindNeighbors(H3.3.combined, reduction = "lsi", assay = "ATAC", dims = 2:30)
H3.3.combined  <- FindClusters(H3.3.combined, resolution = 0.5)
H3.3.combined <- RunUMAP(H3.3.combined, dims = 2:30, reduction = "lsi", reduction.name = "umap.atac", spread=0.5) 
DimPlot(H3.3.combined, reduction = "umap.atac", label = T, label.size = 4, repel = TRUE, split.by = "orig.ident")  + ggtitle("ATAC UMAP")

DefaultAssay(Control) <- "ATAC"
Control<- RunTFIDF(Control)
Control <- FindTopFeatures(Control)
Control <- RunSVD(Control)

Control <- FindNeighbors(Control, reduction = "lsi", assay = "ATAC", dims = 2:30)
Control  <- FindClusters(Control, resolution = 0.5)
Control <- RunUMAP(Control, dims = 2:30, reduction = "lsi") 
DimPlot(Control , reduction = "umap", label = T, label.size = 4, repel = TRUE, split.by = "orig.ident")  + ggtitle("ATAC UMAP")

DefaultAssay(H3.3) <- "ATAC"
H3.3 <- RunTFIDF(H3.3)
H3.3 <- FindTopFeatures(H3.3)
H3.3 <- RunSVD(H3.3)

H3.3 <- FindNeighbors(H3.3, reduction = "lsi", assay = "ATAC", dims = 2:30)
H3.3  <- FindClusters(H3.3, resolution = 0.5)
H3.3 <- RunUMAP(H3.3, dims = 2:30, reduction = "lsi") 
DimPlot(H3.3, reduction = "umap", label = T, label.size = 4, repel = TRUE, split.by = "orig.ident")  + ggtitle("ATAC UMAP")

#Perform integration
integration.anchors.H33.atac <- FindIntegrationAnchors(
  object.list = list(Control, H3.3),
  anchor.features = rownames(Control),
  reduction = "rlsi",
  dims = 2:30
)

H3.3.combined_integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors.H33.atac,
  reductions = H3.3.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)
H3.3.combined[['integrated_lsi_atac']] <- CreateDimReducObject(
  Embeddings(H3.3.combined_integrated, "integrated_lsi")[colnames(H3.3.combined),], key="INTEGRATEDLSIATAC_", assay="ATAC"
)

H3.3.combined <- FindNeighbors(H3.3.combined, reduction = "integrated_lsi_atac", assay = "ATAC", dims = 2:30)
H3.3.combined  <- FindClusters(H3.3.combined, resolution = 0.3)
H3.3.combined <- RunUMAP(H3.3.combined,
                           reduction = "integrated_lsi_atac",
                           dims = 2:30,
                           spread = 1,
                         min.dist = 0.8,
                           reduction.name = "umap_atac.int",
                           reduction.key = "UMAPSEURATATAC_")
DimPlot(H3.3.combined, reduction = "umap_atac.int", label = T, label.size = 4, repel = TRUE, split.by = "orig.ident")  + ggtitle("ATAC UMAP")

# Coverage Plots
DefaultAssay(H3.3.combined) <- "ATAC"
library(BSgenome.Drerio.UCSC.danRer11)
H3.3.combined <- RegionStats(H3.3.combined, genome = BSgenome.Drerio.UCSC.danRer11)
CoveragePlot(
  object = H3.3.combined,
  region = "mitfa",
  features = "mitfa",
  expression.assay = "SCT",
  extend.downstream = 1000,
  split.by = "orig.ident"
)

#make stacked bar plot
ncells_H3 <- table(Idents(H3.3.combined), H3.3.combined$orig.ident)
ncells_H3 <- as.data.frame(ncells_H3)
ncells_H3$Var1 <- as.character(ncells_H3$Var1)
ncells_H3 <- ncells_H3 %>%
group_by(Var2) %>%
mutate(Proportion = Freq / sum(Freq))

ggplot(ncells_H3, aes(x = Var2, y = Proportion, fill = Var1)) +
  +     theme_bw(base_size = 15) +
  +     geom_col(position = "fill", width = 0.5) +
  +     geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
                  +               position = position_fill(vjust = 0.5),  # Center labels in each segment
                  +               size = 4, color = "red") +  # Adjust size and color of labels
  +     xlab("Sample") +
  +     ylab("Proportion") +
  +     scale_fill_viridis_d(option = "mako")  # Use viridis palette for colors