library(Seurat)
library(SingleCellExperiment)
library(scater)
library(TSCAN)

set.seed(1234)

# Making sure Seurat clusters exist
clusters <- Idents(data)

# Converting to SCE
sce <- as.SingleCellExperiment(data)

# Assigning celltype.rna clusters as colLabels
colLabels(sce) <- clusters
by.cluster <- aggregateAcrossCells(sce, ids = colLabels(sce))
centroids <- reducedDim(by.cluster, "PCA")

mst <- createClusterMST(centroids, clusters=NULL)

line.data <- reportEdges(
  by.cluster,
  mst = mst,
  clusters = NULL,
  use.dimred = "UMAP.RNA"
)

#Projecting mst on UMAP
p <- plotReducedDim(sce, dimred="UMAP.RNA", colour_by="label") + 
  geom_line(data=line.data, mapping=aes(x=umaprna_1, y=umaprna_2, group=edge))
p + scale_color_manual(values = color_layers) +
  geom_line(
    data = line.data,
    mapping = aes(x = umaprnacontrol_1, y = umaprnacontrol_2, group = edge)
  )

#Pseudotime coloring
map.tscan <- mapCellsToEdges(sce, mst=mst, use.dimred="PCA")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)

common.pseudo <- averagePseudotime(tscan.pseudo) 
plotReducedDim(sce, dimred="UMAP.RNA", colour_by=I(common.pseudo)) +
  geom_line(data=line.data, mapping=aes(x=umaprna_1, y=umaprna_2, group=edge))

saveRDS(sce, file = "..../sox10_ncc_sce.rds")
