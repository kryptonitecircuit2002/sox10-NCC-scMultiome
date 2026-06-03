library(devtools)
library(sva)
library(reticulate) 
library(ggplot2)

set.seed(1234)
H2.data <- readRDS(".../Integrated_H2_annotated.rds")
data_control_h2 <- subset(H2.data, idents = "Control")
data_h2 <- subset(H2.data, idents = "H2A.Z KD")

umapembeddings_h2 <- Embeddings(H2.data, reduction = "umap.rna.h2")
controlraw_h2 <- GetAssayData(data_control_h2, layer = "counts")
h2raw <- GetAssayData(data_h2, layer = "counts")
controlraw_matrix_h2 <- as.matrix(controlraw_h2)
h2raw_matrix <- as.matrix(h2raw)

#read in data matrix nuclei
datasets_h2 <- list(controlraw_matrix_h2, h2raw_matrix)
results_h2 <- iCytoTRACE(datasets_h2)
plotCytoTRACE(results_h2, emb = umapembeddings_h2)

# Extract CytoTRACE scores
cyto_scores_h2 <- results_h2$CytoTRACE

#Plotting
df_h2 <- data.frame(
  UMAP1 = umapembeddings_h2[,1],
  UMAP2 = umapembeddings_h2[,2],
  CytoTRACE = cyto_scores_h2,
  Condition = H2.data$orig.ident,
  Identity = H2.data$celltype.rna
)


# Plot with faceting
ggplot(df_h2, aes(x = UMAP1, y = UMAP2, color = CytoTRACE)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "H") +
  facet_wrap(~ Condition) +
  theme_classic() +
  labs(title = "CytoTRACE scores by condition")

#Statistical Plot
ggplot(df_h2, aes(x = Condition, y = CytoTRACE, fill = Condition)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) +
  facet_wrap(~Identity, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    title = "CytoTRACE Score by Condition per Cell Type",
    x = "Condition",
    y = "CytoTRACE Score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

