library(ggplot2)
library(dplyr)

file <- "D:/project/04BrIS/HeLa_BrIS_gradient_samples_human_proteins_peptides_count.csv"
data <- read.csv(file)

sample_order <- c("0ng", "0.1ng", "0.5ng", "1ng", "5ng", "10ng", "50ng", "100ng")
level_order <- c("protein", "peptide")
data$Gradient <- factor(data$Gradient, levels = sample_order)
data$Level <- factor(data$Level, levels = level_order)

colors <- c("#31859B", "#EA700D")

create_figure <- function(df, colors, fontsize) {
  p <- ggplot(df, aes(x = Gradient, y = Count, fill = Level, color = Level)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.6,
      linewidth = 0.5,
      alpha = 0.7
    ) +
    geom_text(
      aes(
        label = Count,
        group = Level,
        angle = ifelse(Level == "protein", 90, 0),
        hjust = ifelse(Level == "protein", -0.1, 0.5),
        vjust = ifelse(Level == "protein", 0.5, -0.5)
      ),
      position = position_dodge(width = 0.8),
      size = fontsize / .pt,
      color = "black",
      show.legend = FALSE
    ) +
    scale_fill_manual(values = colors, name = NULL) +
    scale_color_manual(values = colors, name = NULL) +
    scale_x_discrete(labels = function(x) gsub("ng", "", x)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, 40000), breaks = seq(0, 40000, 10000), labels = c(0, 1, 2, 3, 4)) +
    labs(x = NULL, y = expression(paste("Count (", "\u00D7", 10^4, ")"))) +
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(fill = NA, colour = NA),
      axis.line = element_line(colour = "black", linewidth = 0.25),
      axis.text.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = "black", margin = margin(t = 3, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = "black", margin = margin(r = 1, l = 0)),
      axis.ticks = element_line(colour = "black", linewidth = 0.25),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = fontsize, colour = "black"),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt")
    )
  return(p)
}

p <- create_figure(data, colors, fontsize = 6.5)
print(p)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure3")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
ggsave(
  file.path(output_path, "Fig.3d_HeLa_whole_cell_BrIS_gradient_samples_human_proteins_peptides_count.pdf"),
  plot = p,
  width = 3.3,
  height = 1.45,
  family = "Arial",
  device = cairo_pdf
)
