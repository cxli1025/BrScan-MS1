library(ggplot2)
library(dplyr)

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')
file <- 'Fig.3c_HeLa_whole_cell_BrIS_gradient_samples_BrIS_intensity_ratio_BrScan-MS1_PEAKS-Online.csv'
data <- read.csv(file)

sample_order <- c('0.1ng', '0.5ng', '1ng', '5ng', '10ng', '50ng', '100ng')
data$sample <- factor(data$sample, levels = sample_order)
theory_ratio <- unique(data[, c("sample", "theory_ratio")])

colors <- c("#31859B", "#EA700D") 

create_figure <- function(df, colors, fontsize, theory_ratio) {
  p <- ggplot(df, aes(x = sample, y = ratio)) + 
    geom_boxplot(
      aes(fill = method, color = method), 
      position = position_dodge(width = 0.8), 
      width = 0.6, alpha = 0.7, linewidth = 0.5,
      outlier.shape = NA
      ) + 
    geom_point(
      aes(fill = method, color = method),
      position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.8),
      size = 0.5, alpha = 1.0
    ) +
    geom_errorbar(
      data = theory_ratio, 
      aes(x = sample, ymin = theory_ratio, ymax = theory_ratio),
      color = 'darkgrey', 
      width = 1.0, linetype = "dotted", linewidth = 0.5, alpha = 1.0,
      inherit.aes = FALSE 
    ) +
    labs(
      x = NULL,
      y = expression(log[2]~(Intensity~ratio))
    ) +
    scale_x_discrete(labels = function(x) gsub("ng", "", x)) +
    scale_y_continuous(limits = c(-9, 9), breaks = seq(-9, 9, 3)) +
    scale_fill_manual(values = colors, name = NULL) +
    scale_color_manual(values = colors, name = NULL) +
    theme(
      plot.margin  = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = NA),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.ticks = element_line(colour = 'black', linewidth = 0.5),
      axis.line = element_line(colour = 'black', linewidth = 0.5),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = fontsize, colour = 'black'),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt")
    )
    return(p)
}

p <- create_figure(data, colors, 6.5, theory_ratio)
print(p)