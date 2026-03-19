library(ggplot2)

create_figure <- function(data, colors, fontsize) {
  p <- ggplot(data, aes(x = sample, y = RT)) + 
    geom_line(aes(color = peptide, group = peptide), linetype = 'solid', linewidth = 0.25, alpha = 1.0) + 
    geom_point(aes(color = peptide), size = 0.5, alpha = 1.0) +
    labs(x = 'Sample', y = "RT (min)") +
    scale_x_discrete(labels = function(x) gsub("^S", "", x)) +
    scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
    scale_color_manual(values = colors, name = 'iRT', labels = function(x) gsub("^RT-pep ", "", x)) +
    facet_wrap(~ method) +
    theme(
      plot.margin = margin(1,1,1,1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.25),
      strip.text = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 2)),
      strip.background = element_blank(),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0), 
                                 angle = 90, vjust = 0.5),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
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

colors <- c(
  "#31859B", "#51A3B8", "#1D5667",  
  "#EA700D", "#F4A460", "#A34D09",  
  "#7E649E", "#A391BD", "#5F4B75",  
  "#76933C", "#9BBB59"
)

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')

data1 <- read.csv('Supplementary_Figure3e_plasma_BrIS_iRT_samples_on_TIMS2_iRT_RT_Spectronaut.csv')
p1 <- create_figure(data1, colors, unit(6.5, "pt"))
print(p1)

data2 <- read.csv('Supplementary_Figure3g_plasma_BrIS_iRT_samples_on_E480_iRT_RT_Spectronaut.csv')
p2 <- create_figure(data2, colors, unit(6.5, "pt"))
print(p2)