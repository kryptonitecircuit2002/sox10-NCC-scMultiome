library(tidyverse)
library(RColorBrewer)
library(Seurat)
library(SeuratData)
library(future)

set.seed(1234)
data <- readRDS("......../data_motifenriched.rds")

module_scores <- read.csv(".../Module_scoring.csv")
melanocytic_genes <- na.omit(module_scores$Melanocytic.lineage.score)
glial_genes       <- na.omit(module_scores$Glial.Neuronal.lineage.score)
muscle_genes      <- na.omit(module_scores$Muscle.lineage.score)

# Add module scores
data <- AddModuleScore(data, features = list(melanocytic_genes), name = "Melanocytic_score")
data <- AddModuleScore(data, features = list(glial_genes), name = "Neural_score")
data <- AddModuleScore(data, features = list(muscle_genes), name = "Muscle_score")

#Represent using FeaturePlot
FeaturePlot(data, features = "Melanocytic_score1", reduction = "umap.rna")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(data, features = "Neural_score1", reduction = "umap.rna")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(data, features = "Muscle_score1", reduction = "umap.rna")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

saveRDS(data, file = "...../data_module.rds")
