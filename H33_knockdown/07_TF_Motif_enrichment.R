lapply(required_packages, library, character.only = TRUE)
set.seed(1234)
library(dplyr)
library(tidyr)

tfs_of_interest <- c("MA0139.1", "MA0499.2", "MA0500.2", "MA1641.1", "MA1123.2", "MA0003.4", "MA0814.2",
                     "MA1569.1", "MA0811.1", "MA0442.2", "MA0620.3")

#TF Activity Score
DefaultAssay(data) <- "chromvar"
chromvar_data <- GetAssayData(data, assay = "chromvar", layer = "data")

tf_summary <- data.frame()
for (tf in tfs_of_interest) {
  tf_activity <- chromvar_data[tf, ]
  temp_df <- data.frame(
    celltype = data$celltype.rna,
    condition = data$orig.ident,
    activity = tf_activity,
    tf = tf
  )
  tf_summary <- rbind(tf_summary, temp_df)
}

#Summarise 
tf_stats <- tf_summary %>%
  group_by(celltype, condition, tf) %>%
  summarise(
    mean_activity = mean(activity, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

#Statistical Testing
p_values <- tf_summary %>%
  group_by(celltype, tf) %>%
  summarise(
    p_value = if (n_distinct(condition) == 2) {
      tryCatch({
        wilcox.test(activity ~ condition)$p.value
      }, error = function(e) NA)
    } else {
      NA
    },
    .groups = "drop"
  )

tf_stats <- tf_stats %>%
  left_join(p_values, by = c("celltype", "tf"))
tf_stats$p_adj <- p.adjust(tf_stats$p_value, method = "BH")
tf_stats <- tf_stats %>%
  mutate(mean_activity = ifelse(mean_activity < 0, 0, mean_activity))
tf_stats$tf <- factor(tf_stats$tf, levels = tfs_of_interest)
tf_stats$condition <- factor(tf_stats$condition, levels = c("Control", "H3.3 KD"))

#Labels
tf_stats <- tf_stats %>%
  arrange(tf, condition) %>%
  mutate(x_label = paste0(tf, "\n", condition))

tf_cond_order <- tf_stats %>%
  arrange(tf, condition) %>%
  distinct(tf, condition, x_label) %>%
  pull(x_label)

tf_stats$x_label <- factor(tf_stats$x_label, levels = tf_cond_order)

#Slope Plot
celltype_levels <- unique(tf_stats$celltype)
n_celltypes <- length(celltype_levels)

celltype_cols   <- c(
  '#e24041', '#efd129', '#b06500', '#808000', '#7ac143', '#40733e', '#ff9900',
  '#3ac4e7', '#277ea6', '#182953', '#ab93c6', '#b22987', '#f64a8a'
)
celltype_colors <- setNames(
  celltype_cols[seq_len(n_celltypes)],
  celltype_levels
)

p_slope <- ggplot(tf_stats,
                  aes(x = condition,
                      y = mean_activity,
                      group = celltype,
                      color = celltype)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             color      = "grey60",
             linewidth  = 0.4) +
  scale_color_manual(values = celltype_colors, name = "Cell type") +
  facet_wrap(~ tf, ncol = 3, scales = "free_y") +
  labs(
    title = "TF chromatin activity: Control vs H3.3 KD",
    x = NULL,
    y = "Mean chromvar activity"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold",
                        size = 14, colour = "black"),
    strip.text = element_text(size = 14, face = "bold", family = "Times New Roman"),
    strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.4),
    axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.key.size = unit(0.5, "cm"),
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

p_slope

write.csv(tf_stats, file = ".../Motif_activity.csv")
