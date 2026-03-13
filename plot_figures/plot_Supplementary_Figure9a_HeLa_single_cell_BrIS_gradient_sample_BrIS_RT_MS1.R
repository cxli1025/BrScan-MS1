library(ggplot2)

file <- 'D:/project/04BrIS/hela_single_cell_gradient_samples/BrIS/new_ms1/RTs.csv'
data <- read.csv(file)

sample_order <- c('0.1pg', '0.5pg', '1pg', '5pg', '10pg', '50pg', '100pg')
data$sample <- factor(data$sample, levels = sample_order)

colors <- c(
  "#31859B", "#51A3B8", "#1D5667", 
  "#EA700D", "#F4A460", "#A34D09", 
  "#7E649E", "#A391BD", "#5F4B75"  
)

create_figure <- function(data, colors, fontsize){
  p <- ggplot(data, aes(x = sample, y = RT)) + 
    geom_line(aes(color = peptide,group = peptide), linetype = 'solid', linewidth = 0.5, alpha = 1.0) + 
    geom_point(aes(color = peptide), size = 0.5, alpha = 1.0) +
    labs(x = NULL, y = 'RT (min)', title = "BrScan-MS1") + 
    scale_x_discrete(drop = FALSE, labels = function(x) gsub("pg", "", x)) +  
    scale_y_continuous(limits = c(10, 30), breaks = seq(10, 30, 5)) +
    scale_color_manual(values = colors, name = "BrIS", labels = function(x) gsub("^BrIS-", "", x)) +
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      plot.title = element_text(size = fontsize, colour = 'black', margin = margin(t= 2, b = 2), hjust = 0.5),
      axis.ticks = element_line(colour = 'black', linewidth = 0.5),
      legend.position = "right",
      legend.text = element_text(size = fontsize, colour = 'black'),
      legend.title = element_text(size = fontsize, colour = 'black'),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt")
    )
  return(p)
}

p <- create_figure(data, colors, fontsize = unit(6.5, "pt"))
print(p)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Supplementary_Figures", "Figure9")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
ggsave(
  file.path(output_path, "Supplementary_Figure9a_HeLa_single_cell_BrIS_gradient_sample_BrIS_RT_MS1.pdf"),
  plot = p,
  width = 2.0,
  height = 1.9,
  family = 'Arial',
  device = cairo_pdf
  )