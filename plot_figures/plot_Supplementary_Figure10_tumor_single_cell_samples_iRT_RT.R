library(tidyverse)
library(ggplot2)
library(patchwork)

load_data <- function(file, sample_type) {
  df <- read.csv(file)
  df_clean <- df %>% 
    select(Peptide, contains("ApexRT"))
  colnames(df_clean) <- c(
    "peptide", 
    paste0(sample_type, "B", 1:5), 
    paste0(sample_type, "I", 1:5), 
    paste0(sample_type, "C", 1:5)
  )
  df_long <- df_clean %>%
    mutate(across(-peptide, as.character)) %>%
    pivot_longer(
      cols = -peptide, 
      names_to = "sample", 
      values_to = "RT"
    ) %>%
    mutate(RT = ifelse(RT == "Filtered", NA, RT)) %>%
    mutate(RT = as.numeric(RT))
  sample_order <- c(
    paste0(sample_type, "C", 1:5), 
    paste0(sample_type, "B", 1:5), 
    paste0(sample_type, "I", 1:5)
  )
  df_long$sample <- factor(df_long$sample, levels = sample_order)
  return(df_long)
}

create_figure <- function(df, colors, fontsize) {
  p <- ggplot(df, aes(x = sample, y = RT)) +
    geom_line(
      data = df %>% filter(!is.na(RT)),
      aes(color = peptide, group = peptide),
      linetype = 'solid', linewidth = 0.5, alpha = 1.0
    ) + 
    geom_point(aes(color = peptide), size = 0.5, alpha = 1.0) +   
    labs(x = "Sample", y = "RT (min)") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    scale_color_manual(values = colors, name = 'iRT') + 
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 2)),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 10)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
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

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')
hela_data <- load_data("Supplementary_Figure10_HeLa_single_cell_samples_iRT_RT.csv", "H")
n460_data <- load_data("Supplementary_Figure10_NCI-H460_single_cell_samples_iRT_RT.csv", "N")
m231_data <- load_data("Supplementary_Figure10_MDA-MB-231_single_cell_samples_iRT_RT.csv", "M")

colors <- c(
  "#31859B", "#51A3B8", "#1D5667",  
  "#EA700D", "#F4A460", "#A34D09",  
  "#7E649E", "#A391BD", "#5F4B75",  
  "#76933C", "#9BBB59"
)
p_hela <- create_figure(hela_data, colors, fontsize = unit(6.5, "pt"))
p_n460 <- create_figure(n460_data, colors, fontsize = unit(6.5, "pt"))
p_m231 <- create_figure(m231_data, colors, fontsize = unit(6.5, "pt"))

p <- p_hela / p_n460 / p_m231
print(p)