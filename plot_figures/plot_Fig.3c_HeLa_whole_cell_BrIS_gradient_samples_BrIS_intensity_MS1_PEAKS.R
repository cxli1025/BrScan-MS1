library(ggplot2)
library(dplyr)

df_ms1 <- read.csv('D:/project/04BrIS/HeLa_whole_cell_BrIS_gradient_samples/areas.csv')
df_peaks <- read.csv('D:/project/04BrIS/HeLa_whole_cell_BrIS_gradient_samples_peaks/BrIS/areas.csv')

df_ms1 <- df_ms1 %>% mutate(method = 'BrScan-MS1')
df_peaks <- df_peaks %>% mutate(method = 'PEAKS Online')

df <- rbind(df_ms1, df_peaks)

sample_order <- c('0.1ng', '0.5ng', '1ng', '5ng', '10ng', '50ng', '100ng')
df$sample <- factor(df$sample, levels = sample_order)
theory_ratio <- unique(df[, c("sample", "theory_ratio")])

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
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
      axis.line = element_line(colour = 'black', linewidth = 0.25),
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

p <- create_figure(df, colors, 6.5, theory_ratio)
print(p)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure3")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
ggsave(
  file.path(output_path, "Fig.3c_HeLa_whole_cell_BrIS_gradient_samples_BrIS_intensity_MS1_PEAKS.pdf"),
  plot = p,
  width = 3.1,
  height = 1.7,
  family = 'Arial',
  device = cairo_pdf)