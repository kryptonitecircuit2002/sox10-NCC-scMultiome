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

data <- readRDS( ".../data_linked.rds")
DotPlot(data, features = c("twist1a", "grem2b", "col11a1b", "lamc3", "foxd3", "crestin", "sox10", "zeb2a", "aox5", "paics", "dct", "mitfa", "tyr", "ltk", "alx4a", "gfap", "pou3f1", "her12", "elavl3", "elavl4", "oc90", "epcam", "neurod1", "musk", "tpma")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="viridis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() + 
  coord_flip() + ggtitle("Markers for Cluster Annotation")

#RNA Derived
Idents(data) <- data$SCT_snn_res.5
celltype <- rep(NA, length = ncol(data))

celltype[which(Idents(data) %in% c(0,10))] <- 'twist 1a+ NCC'
celltype[which(Idents(data) %in% c(1))] <- 'mesenchymal'
celltype[which(Idents(data) %in% c(4))] <- 'unknown'
celltype[which(Idents(data) %in% c(2))] <- 'MX+'
celltype[which(Idents(data) %in% c(5))] <- 'pigment progenitor'
celltype[which(Idents(data) %in% c(7))] <- 'foxd3+ NCC'
celltype[which(Idents(data) %in% c(9))] <- 'melanoblast'
celltype[which(Idents(data) %in% c(12))] <- 'iridoblast'
celltype[which(Idents(data) %in% c(6))] <- 'neuronal'
celltype[which(Idents(data) %in% c(8))] <- 'otic'
celltype[which(Idents(data) %in% c(3))] <- 'glial/neural'
celltype[which(Idents(data) %in% c(13))] <- 'immature neurons'
celltype[which(Idents(data) %in% c(11))] <- 'skeletal muscle progenitor'
celltype <- factor(celltype, 
                   levels = c('twist 1a+ NCC', 'mesenchymal', 'foxd3+ NCC', 'pigment progenitor', 'MX+', 
                              'melanoblast', 'iridoblast', 'unknown', 'glial/neural', 
                              'neuronal', 'otic', 'immature neurons', 'skeletal muscle progenitor'), 
                   ordered = T)

data$celltype.rna <- celltype
Idents(data) <- data$celltype.rna
DimPlot(data, reduction = "umap.rna.control", label = T, label.size = 5, repel = T, stroke.size = 1.1) + ggtitle("RNA clusters in sox10+ NCCs")

##ATAC derived
DefaultAssay(data) <- "peaks"
gene.activities <- GeneActivity(data)
data[["Gene Activity"]] <- gene.activites
DefaultAssay(data) <- "Gene Activity"
DotPlot(data, features = c("twist1a", "grem2b", "col11a1b", "lamc3", "foxd3", "crestin", "sox10", "zeb2a", "aox5", "paics", "dct", "mitfa", "tyr", "ltk", "alx4a", "gfap", "pou3f1", "her12", "elavl3", "elavl4", "oc90", "epcam", "neurod1", "musk", "tpma")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="viridis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() + 
  coord_flip() + ggtitle("Markers for Cluster Annotation")

DefaultAssay(data) <- "peaks"
Idents(data) <- data$ATAC_snn_res.5
celltype <- rep(NA, length = ncol(data))

celltype[which(Idents(data) %in% c(0,7,12))] <- 'twist 1a+ NCC'
celltype[which(Idents(data) %in% c(1,14,15))] <- 'mesenchymal'
celltype[which(Idents(data) %in% c(2))] <- 'MX+'
celltype[which(Idents(data) %in% c(11))] <- 'pigment progenitor'
celltype[which(Idents(data) %in% c(4))] <- 'foxd3+ NCC'
celltype[which(Idents(data) %in% c(6))] <- 'melanoblast'
celltype[which(Idents(data) %in% c(16,17))] <- 'neuronal'
celltype[which(Idents(data) %in% c(9,10))] <- 'otic'
celltype[which(Idents(data) %in% c(5,8))] <- 'glial/neural'
celltype[which(Idents(data) %in% c(13))] <- 'skeletal muscle progenitor'
celltype <- factor(celltype, 
                   levels = c('twist 1a+ NCC', 'mesenchymal', 'foxd3+ NCC', 'pigment progenitor', 'MX+', 
                              'melanoblast', 'glial/neural', 
                              'neuronal', 'otic', 'skeletal muscle progenitor'), 
                   ordered = T)
data$celltype.atac <- celltype
DimPlot(data, reduction = "umap.peaks", label = T, label.size = 5, repel = T, stroke.size = 1.1) + ggtitle("ATAC clusters in sox10+ NCCs")

saveRDS("..../data_annotated.rds")
