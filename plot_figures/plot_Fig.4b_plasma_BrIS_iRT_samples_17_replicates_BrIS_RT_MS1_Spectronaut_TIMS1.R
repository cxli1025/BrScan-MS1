library(ggplot2)
library(dplyr)

setwd('D:/project/04BrIS/across_platform/BrIS/tims1')
df1 <- read.csv('ms1/samples/RTs.csv')
df2 <- read.csv('Spectronaut/RTs.csv')
df1$method <- 'BrScan-MS1'
df2$method <- 'Spectronaut'
df <- rbind(df1, df2)

colors <- c(
  "#31859B", "#51A3B8", "#1D5667", 
  "#EA700D", "#F4A460", "#A34D09", 
  "#7E649E", "#A391BD", "#5F4B75"  
)

create_figure <- function(df, colors, fontsize) {
  p <- ggplot(df, aes(x = sample, y = RT)) +
    geom_line(aes(color = peptide, group = peptide), linetype = 'solid', linewidth = 0.25, alpha = 1.0) + 
    geom_point(aes(color = peptide), size = 0.5, alpha = 1.0) +
    labs(x = 'Sample', y = "RT (min)") +
    scale_x_discrete(labels = function(x) gsub("^S", "", x)) +
    scale_y_continuous(limits = c(10, 60), breaks = seq(10, 60, 10)) +
    scale_color_manual(values = colors, name = 'BrIS', labels = function(x) gsub("^BrIS-", "", x)) +
    facet_wrap(~ method, ncol = 2, scales = 'free_x') +
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt"),
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

p1 <- create_figure(df, colors, fontsize = unit(6.5, "pt"))
print(p1)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure4")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

ggsave(
  file.path(output_path, "Fig.4b_plasma_BrIS_iRT_samples_17_replicates_BrIS_RT_MS1_Spectronaut_TIMS1.pdf"),
  plot = p,
  width = 4.2,
  height = 1.9,
  family = 'Arial',
  device = cairo_pdf
)