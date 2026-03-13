library(ggplot2)
library(dplyr)

load_data <- function(workdir) {
  setwd(workdir)
  df1 <- read.csv('ms1/samples/areas.csv')
  df2 <- read.csv('Spectronaut/areas.csv')
  df1 <- df1 %>% mutate(method = 'BrScan-MS1')
  df2 <- df2 %>% mutate(method = 'Spectronaut')
  df <- rbind(df1, df2)
  return(df)
}

create_intensity_plot <- function(df, colors, fontsize) {
  summary_df <- df %>%
    group_by(method, peptide) %>%
    summarise(
      mean_area = mean(areas, na.rm = TRUE),
      sd_area   = sd(areas, na.rm = TRUE),
      n         = n(),
      se_area   = sd_area / sqrt(n),
      .groups   = 'drop'
    )
  p <- ggplot() +
    geom_bar(
      data = summary_df,
      aes(x = peptide, y = mean_area, color = peptide),
      stat = "identity",
      width = 0.6,
      fill = "white",
      linewidth = 0.25,
      alpha = 1.0
    ) +
    geom_point(
      data = df,
      aes(x = peptide, y = areas, color = peptide),
      size = 0.6, alpha = 1.0,
      position = position_jitter(width = 0.1, height = 0)
    ) +
    geom_errorbar(
      data = summary_df,
      aes(x = peptide,
          ymin = mean_area - se_area,
          ymax = mean_area + se_area),
      width = 0.3,
      linewidth = 0.25,
      color = "black",
      alpha = 1.0
    ) +
    scale_color_manual(values = colors) +
    scale_x_discrete(labels = function(x) gsub("^BrIS-", "", x)) +
    labs(x = "BrIS", y = expression(log[10]~(Intensity))) +
    facet_wrap(~ method, ncol = 2, scales = "free_x") +
    theme(
      plot.margin  = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.25),
      strip.text = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 2)),
      strip.background = element_blank(),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y  = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.ticks   = element_line(colour = 'black', linewidth = 0.25),
      legend.position = "none"
    )
  return(p)
}

colors <- c(
  "#31859B", "#51A3B8", "#1D5667", 
  "#EA700D", "#F4A460", "#A34D09", 
  "#7E649E", "#A391BD", "#5F4B75"  
)

df1 <- load_data('D:/project/04BrIS/across_platform/BrIS/tims2')
p1 <- create_intensity_plot(df1, colors, fontsize = unit(6.5, "pt")) +
  scale_y_continuous(breaks = seq(0, 8, 2), limits = c(0, 8))
print(p1)

df2 <- load_data('D:/project/04BrIS/across_platform/BrIS/480')
p2 <- create_intensity_plot(df2, colors, fontsize = unit(6.5, "pt")) +
  scale_y_continuous(breaks = seq(0, 12, 3), limits = c(0, 12))
print(p2)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Supplementary_Figures", "Figure3")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

ggsave(
  file.path(output_path, "Supplementary_Figure3b_plasma_BrIS_iRT_samples_BrIS_intensity_MS1_Spectronaut_TIMS2.pdf"),
  plot = p1,
  width = 3.15,
  height = 1.9,
  family = 'Arial',
  device = cairo_pdf
)

ggsave(
  file.path(output_path, "Supplementary_Figure3d_plasma_BrIS_iRT_samples_BrIS_intensity_MS1_Spectronaut_E480.pdf"),
  plot = p2,
  width = 3.15,
  height = 1.9,
  family = 'Arial',
  device = cairo_pdf
)