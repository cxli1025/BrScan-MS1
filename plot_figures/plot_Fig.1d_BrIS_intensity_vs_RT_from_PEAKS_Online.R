library(ggplot2)
library(ggrepel)

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')
file <- 'Fig.1d_wgl_dda_30min_20240705_10_Slot1-36_1_27757_results_PEAKS_Online.csv'
data <- read.csv(file)
data <- data[, c('Accession', 'RT', 'Intensity.Sample.1')]
data$log10_Intensity <- log10(data$Intensity.Sample.1)

colors <- c(
  "#31859B", "#51A3B8", "#1D5667", 
  "#EA700D", "#F4A460", "#A34D09", 
  "#7E649E", "#A391BD", "#5F4B75"  
)

create_figure <- function(data, colors, fontsize) {
  p <- ggplot(data, aes(x = RT, y = log10_Intensity, color = Accession)) +
    geom_segment(aes(x = RT, xend = RT, y = 0, yend = log10_Intensity), linewidth = 0.5, alpha = 1.0) +
    labs(x = "RT (min)", y = expression(log[10]~(Intensity))) +
    scale_x_continuous(breaks = c(10, 15, 20, 25, 30), limits = c(10, 30)) +
    scale_y_continuous(breaks = c(0, 3, 6, 9), limits = c(0, 9.5), expand = c(0, 0)) +
    scale_color_manual(values = colors, name = NULL) +
    theme(
      plot.margin = margin(1,1,1,1, "pt"),  
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = NA),
      legend.position = "right",
      legend.text = element_text(size = fontsize, colour = 'black'),
      legend.title = element_text(size = fontsize, colour = 'black'),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.ticks = element_line(colour = 'black', linewidth = 0.5),
      axis.line = element_line(colour = 'black', linewidth = 0.5),
    )
  return(p)
}

p <- create_figure(data, colors, fontsize = unit(6.5, "pt"))
print(p)