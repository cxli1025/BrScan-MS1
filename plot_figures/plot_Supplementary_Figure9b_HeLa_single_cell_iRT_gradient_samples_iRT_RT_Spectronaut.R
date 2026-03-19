library(ggplot2)

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')
file <- 'Supplementary_Figure9b_HeLa_single_cell_BrIS_gradient_samples_iRT_RT_Spectronaut.csv'
data <- read.csv(file)

sample_order <- c('0.1pg', '0.5pg', '1pg', '5pg', '10pg', '50pg', '100pg')
data$sample <- factor(data$sample, levels = sample_order)

colors <- c(
  "#31859B", "#51A3B8", "#1D5667",  
  "#EA700D", "#F4A460", "#A34D09",  
  "#7E649E", "#A391BD", "#5F4B75",  
  "#76933C", "#9BBB59"
)

create_figure <- function(data, colors, fontsize){
  p <- ggplot(data, aes(x = sample, y = RT)) + 
    geom_line(aes(color = peptide,group = peptide), linetype = 'solid', linewidth = 0.5, alpha = 1.0) + 
    geom_point(aes(color = peptide), size = 0.5, alpha = 1.0) +
    labs(x = NULL, y = 'RT (min)', title = "Spectronaut") + 
    scale_x_discrete(drop = FALSE, labels = function(x) gsub("pg", "", x)) +  
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    scale_color_manual(values = colors, name = 'iRT', labels = function(x) gsub("^RT-pep ", "", x)) +
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      plot.title = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 2), hjust = 0.5),
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