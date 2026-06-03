lapply(required_packages, library, character.only = TRUE)
library(tidyverse)
library(RColorBrewer)
library(Seurat)
library(SeuratData)
library(future)
library(ggplot2)

set.seed(1234)

H2.combined <- readRDS("E:/Transit/Single_Cell+ATAC/New_Integrative_Analysis/H2AZ_sample/Data/07.08.2025/Integrated_H2_annotated.rds")
module_scores <- read.csv("E:/Transit/Single_Cell+ATAC/Module_scoring.csv")
pigment_genes <- na.omit(module_scores$Pigment.lineage.score)
glial_genes <- na.omit(module_scores$Glial.lineage.score)
neuronal_genes <- na.omit(module_scores$Neuronal.lineage.score)
muscle_genes <- na.omit(module_scores$Muscle.lineage.score)

# Add module scores
H2.combined <- AddModuleScore(H2.combined, features = list(pigment_genes), name = "Pigment_score")
H2.combined <- AddModuleScore(H2.combined, features = list(neuronal_genes), name = "Neuronal_score")

#Module scoring heatmap
celltypes_to_show <- c("mesenchymal", "twist1a+ NCC", "unknown", "glial/neural", "neuronal", "otic", "foxd3+ NCC", "pigment progenitor", "MX+", "MI+", "melanoblast")
meta_df <- H2.combined@meta.data %>%
  dplyr::select(celltype.rna_n, orig.ident, Pigment_score1, Neuronal_score1) %>%
  dplyr::filter(celltype.rna_n %in% celltypes_to_show)

#Mean score
agg_scores <- meta_df %>%
  group_by(celltype.rna_n, orig.ident) %>%
  summarise(across(c(Pigment_score1, Neuronal_score1), mean, na.rm = TRUE)) %>%
  ungroup()

#heatmap matrix
heatmap_mat <- agg_scores %>%
  pivot_wider(names_from = orig.ident, values_from = c(Pigment_score1, Neuronal_score1))
heatmap_mat <- as.data.frame(heatmap_mat)
rownames(heatmap_mat) <- heatmap_mat$celltype.rna_n
heatmap_mat <- heatmap_mat %>% dplyr::select(-celltype.rna_n) %>% as.matrix()

#Scaling
scaled_mat <- t(scale(t(heatmap_mat), center = TRUE, scale = TRUE))
scaled_mat[is.na(scaled_mat)] <- 0
scaled_mat_pos <- scaled_mat
scaled_mat_pos[scaled_mat_pos < 0] <- 0

pheatmap(
  scaled_mat_pos,
  color = viridis::viridis(100, option = "magma"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  scale = "none",
  main = "Module scores (row z-score)",
  display_numbers = FALSE,     
  number_format = "%.2f"
)

write.csv(scaled_mat_pos, file = ".../module_scores.csv")
