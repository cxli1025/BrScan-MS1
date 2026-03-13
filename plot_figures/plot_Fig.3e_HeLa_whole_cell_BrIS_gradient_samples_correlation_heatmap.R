library(ggplot2)
library(reshape2)

file <- 'D:/project/04BrIS/HeLa_whole_cell_BrIS_gradient_samples_peaks/human/lfq.proteins.csv'
df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
cols_to_keep <- grep("Accession|Sample.*Area", colnames(df), value = TRUE)
df_subset <- df[, cols_to_keep, drop = FALSE]

sample_cols_indices <- which(colnames(df_subset) != "Accession")
current_names <- colnames(df_subset)[sample_cols_indices]

new_sample_names <- paste0("S", seq(1, length(current_names)))
colnames(df_subset)[sample_cols_indices] <- new_sample_names

rownames(df_subset) <- df_subset$Accession
df_subset$Accession <- NULL

expression_matrix <- as.matrix(df_subset)
expression_matrix[expression_matrix == 0] <- NA
expression_matrix_log2 <- log2(expression_matrix)

correlation_matrix <- cor(expression_matrix_log2, use = "pairwise.complete.obs", method = "pearson")
my_breaks <- seq(0.5, 1, length.out = 101)
melted_cor <- melt(correlation_matrix)
colnames(melted_cor) <- c("Var1", "Var2", "Correlation")

sample_levels <- paste0("S", 1:ncol(correlation_matrix))
melted_cor$Var1 <- factor(melted_cor$Var1, levels = rev(sample_levels)) # 纵轴逆序显示，使 S1 在最上方
melted_cor$Var2 <- factor(melted_cor$Var2, levels = sample_levels)      # 横轴正序

my_colors <- colorRampPalette(c("navy", "white", "#B2182B"))(100)

fontsize <- 6.5

p <- ggplot(melted_cor, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = my_colors, 
    limits = c(0.5, 1), 
    breaks = c(0.5, 0.75, 1.0), 
    labels = c("0.5", "0.75", "1.0"),
    na.value = "grey90",
    guide = guide_colorbar(
      title = NULL,
      frame.colour = NA,
      ticks.colour = NA,
      barwidth = unit(20, "mm"),    
      barheight = unit(2.0, "mm"),
    )
  ) +
  scale_x_discrete(position = "top") + 
  scale_y_discrete(position = "left") +
  labs(x = NULL, y = NULL) +
  theme(
    plot.margin = margin(1, 1, 1, 1, "pt"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = NA),
    axis.text.x = element_text(size = fontsize, colour = "black", margin = margin(t = 0, b = 0)),
    axis.text.y = element_text(size = fontsize, colour = "black", margin = margin(r = 0, l = 0)),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = fontsize, colour = "black"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt"),
    aspect.ratio = 1 
  )

print(p)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure3")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
ggsave(
  file.path(output_path, "Fig.3e_HeLa_whole_cell_BrIS_gradient_samples_correlation_heatmap.pdf"),
  plot = p,
  width = 1.2,
  height = 1.6,
  family = 'Arial',
  device = cairo_pdf
  )