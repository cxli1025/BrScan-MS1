library(ggplot2)
library(dplyr)

load_data <- function(workdir) {
  setwd(workdir)
  df1 <- read.csv('ms1/samples/RTs.csv')
  df2 <- read.csv('Spectronaut/RTs.csv')
  df1 <- df1 %>% mutate(method = 'BrScan-MS1')
  df2 <- df2 %>% mutate(method = 'Spectronaut')
  df <- rbind(df1, df2)
  return(df)
}

create_figure <- function(df, colors, fontsize) {
  p <- ggplot(df, aes(x = sample, y = RT)) +
    geom_line(aes(color = peptide, group = peptide), linetype = 'solid', linewidth = 0.25, alpha = 1.0) + 
    geom_point(aes(color = peptide), size = 0.5, alpha = 1.0) +
    labs(x = 'Sample', y = "RT (min)") +
    scale_x_discrete(labels = function(x) gsub("^S", "", x)) +
    scale_y_continuous(limits = c(20, 60), breaks = c(20, 30, 40, 50, 60)) +
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

colors <- c(
  "#31859B", "#51A3B8", "#1D5667", 
  "#EA700D", "#F4A460", "#A34D09", 
  "#7E649E", "#A391BD", "#5F4B75"  
)

df1 <- load_data('D:/project/04BrIS/across_platform/BrIS/tims2')
p1 <- create_figure(df1, colors, fontsize = unit(6.5, "pt"))
print(p1)

df2 <- load_data('D:/project/04BrIS/across_platform/BrIS/480')
p2 <- create_figure(df2, colors, fontsize = unit(6.5, "pt"))
print(p2)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Supplementary_Figures", "Figure3")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

ggsave(
  file.path(output_path, "Supplementary_Figure3a_plasma_BrIS_iRT_samples_BrIS_RT_MS1_Spectronaut_TIMS2.pdf"),
  plot = p1,
  width = 3.5,
  height = 1.9,
  family = 'Arial',
  device = cairo_pdf
)

ggsave(
  file.path(output_path, "Supplementary_Figure3c_plasma_BrIS_iRT_samples_BrIS_RT_MS1_Spectronaut_E480.pdf"),
  plot = p2,
  width = 3.5,
  height = 1.9,
  family = 'Arial',
  device = cairo_pdf
)