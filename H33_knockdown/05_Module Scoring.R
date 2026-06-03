lapply(required_packages, library, character.only = TRUE)
set.seed(1234)
library(RColorBrewer)
library(future)
library(ggplot2)
library(viridis)

module_scores <- read.csv(".../Module_scoring.csv")
pigment_genes <- na.omit(module_scores$Pigment.lineage.score)
muscle_genes <- na.omit(module_scores$Muscle.lineage.score)

# Add module scores
H3.combined <- AddModuleScore(H3.combined, features = list(pigment_genes), name = "Pigment_score")
H3.combined <- AddModuleScore(H3.combined, features = list(muscle_genes), name = "Muscle_score")

##Module Scoring Heatmap
# Metadata
celltypes_to_show <- c("mesenchymal", "twist1a+ NCC", "unknown", "muscle progenitor", 
                       "foxd3+ NCC", "pigment progenitor", "MX+", "MI+", "melanoblast")
meta_df <- H3.combined@meta.data %>%
  dplyr::select(celltype.rna_n, orig.ident, Pigment_score1, Muscle_score1) %>%
  dplyr::filter(celltype.rna_n %in% celltypes_to_show)

#Aggregate mean scores
agg_scores <- meta_df %>%
  group_by(celltype.rna_n, orig.ident) %>%
  summarise(across(c(Pigment_score1, Muscle_score1), mean, na.rm = TRUE)) %>%
  ungroup()

#Heatmap matrix
heatmap_mat <- agg_scores %>%
  pivot_wider(names_from = orig.ident, values_from = c(Pigment_score1, Muscle_score1))
heatmap_mat <- as.data.frame(heatmap_mat)
rownames(heatmap_mat) <- heatmap_mat$celltype.rna_n
heatmap_mat <- heatmap_mat %>% dplyr::select(-celltype.rna_n) %>% as.matrix()

#Row-scale
scaled_mat <- t(scale(t(heatmap_mat), center = TRUE, scale = TRUE))
scaled_mat[is.na(scaled_mat)] <- 0
scaled_mat_pos <- scaled_mat
scaled_mat_pos[scaled_mat_pos < 0] <- 0

pheatmap(
  scaled_mat_pos,
  color = viridis::viridis(200, option = "magma"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  scale = "none",
  main = "Module scores (row z-score)"
)
saveRDS(H3.combined, file = ".../h3_modulescore.rds")
write.csv(scaled_mat_pos, file = ".../Module_scores.csv")
