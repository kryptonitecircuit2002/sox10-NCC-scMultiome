library(devtools)
library(sva)
library(reticulate)
library(CytoTRACE)

use_condaenv("cytotrace39")
py_install("scanoramaCT",pip=TRUE)
numpy <- import("numpy")
scanoramaCT <- import("scanoramaCT")


H3.data <- readRDS(".../h3_integrated_linked.rds")
data_control <- subset(H3.data, idents = "Control")
data_h3 <- subset(H3.data, idents = "H3.3 KD")

umapembeddings <- Embeddings(H3.data, reduction = "umap.rna.h3")
controlraw <- GetAssayData(data_control, layer = "counts.Control.1")
h3raw <- GetAssayData(data_h3, layer = "counts.H3.3 KD.2")
controlraw_matrix <- as.matrix(controlraw)
h3raw_matrix <- as.matrix(h3raw)

#read in data matrix nuclei
datasets <- list(controlraw_matrix, h3raw_matrix)
results <- iCytoTRACE(datasets)
plotCytoTRACE(results, emb = umapembeddings, phenotype = "orig.ident")

# Extract CytoTRACE scores
cyto_scores <- results$CytoTRACE

#Plotting
df <- data.frame(
  UMAP1 = umapembeddings[,1],
  UMAP2 = umapembeddings[,2],
  CytoTRACE = cyto_scores,
  Condition = H3.data$orig.ident,
  Identity = H3.data$celltype.rna
)

ggplot(df, aes(x = UMAP1, y = UMAP2, color = CytoTRACE)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "H") +
  facet_wrap(~ Condition) +
  theme_classic() +
  labs(title = "CytoTRACE scores by condition")

ggplot(df, aes(x = Identity, y = CytoTRACE, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "white") +
  theme_minimal(base_size = 14) +
  labs(
    title = "CytoTRACE Score per Cell Type and Condition",
    x = "Cell Type",
    y = "CytoTRACE Score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Statistical Plot
ggplot(df, aes(x = Condition, y = CytoTRACE, fill = Condition)) +
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
