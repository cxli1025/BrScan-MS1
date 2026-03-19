library(ggplot2)

create_figure <-  function(df, colors, fontsize) {
  summary_df <- df %>%
    group_by(peptide) %>%
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
      linewidth = 0.25
    ) +
    geom_point(
      data = df,
      aes(x = peptide, y = areas, color = peptide),
      size = 0.5, alpha = 1.0,
      position = position_jitter(width = 0.1, height = 0)
    ) +
    geom_errorbar(
      data = summary_df,
      aes(x = peptide,
          ymin = mean_area - se_area,
          ymax = mean_area + se_area),
      width = 0.3,
      linewidth = 0.25,
      color = "black"
    ) +
    scale_color_manual(values = colors) +
    scale_x_discrete(labels = function(x) gsub("^RT-pep ", "", x)) +
    labs(x = "iRT peptide", y = expression(log[10](Intensity))) +
    facet_wrap(~method) +
    theme(
      plot.margin  = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.25),
      strip.text = element_text(size = fontsize, colour = 'black', margin = margin(t = 0, b = 2)),
      strip.background = element_blank(),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
      legend.position = "none"
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

df1 <- read.csv('Supplementary_Figure3f_plasma_BrIS_iRT_samples_on_TIMS2_iRT_intensity_Spectronaut.csv')
p1 <- create_figure(df1, colors, unit(6.5, "pt")) +
  scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0, 5))
print(p1)

df2 <- read.csv('Supplementary_Figure3h_plasma_BrIS_iRT_samples_on_E480_iRT_intensity_Spectronaut.csv')
p2 <- create_figure(df2, colors, unit(6.5, "pt")) +
  scale_y_continuous(breaks = seq(0, 8, 2), limits = c(0, 8))
print(p2)