library(ggplot2)

file <- 'D:/project/04BrIS/three_cancer_single_cell/iRT/n460/RTs.csv'
data <- read.csv(file)

colors <- c(
  "#86b637", "#7bb02b", "#a6d854", 
  "#de77ae", "#bf5a95", "#ae447a", "#e78ac3", 
  "#6d8ab8", "#4a66a0", "#5e77ad", "#8da0cb"
)

create_figure <- function(data, colors, fontsize){
  p <- ggplot(data, aes(x = sample, y = RT)) + 
    geom_line(aes(color = peptide,group = peptide), linetype = 'solid', linewidth = 0.25, alpha = 0.8) + 
    geom_point(aes(color = peptide), size = 1) +
    labs(x = 'sample', y = 'RT (min)', title = 'N460') + 
    scale_x_continuous(breaks = seq(1,8,1)) +  
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    scale_color_manual(values = colors, name = 'iRT', labels = function(x) gsub("^RT-", "", x)) +
    facet_wrap(~ method, ncol = 2, scales = 'free_x') +
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.25),
      plot.title = element_text(size = 10, colour = 'black', face = 'bold'),
      strip.text = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 2)),
      strip.background = element_rect(fill = '#F0F0F0', colour = 'black', linewidth = 0.25),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 3, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
      legend.position = "right",
      legend.text = element_text(size = fontsize, colour = 'black'),
      legend.title = element_text(size = fontsize, colour = 'black'),
      legend.key.size = unit(0.2, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(3, "pt")
    )
  return(p)
}

p <- create_figure(data, colors, fontsize = unit(6.5, "pt"))
print(p)