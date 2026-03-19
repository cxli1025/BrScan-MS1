library(ggplot2)
library(dplyr)
library(ggh4x)
library(cowplot)

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')
file <- 'Fig.4e_plasma_BrIS_iRT_samples_across_platforms_BrIS_iRT_RT_Intensity_CV.csv'
df <- read.csv(file)

instrument_levels <- c("timsTOF Pro 1", "timsTOF Pro 2", "Exploris 480")
instrument_symbols <- c(
  "timsTOF Pro 1" = "A", 
  "timsTOF Pro 2" = "B",
  "Exploris 480"  = "C"
)

method_levels <- c("BrIS-MS1", "BrIS-Spectronaut", "iRT-Spectronaut")
method_symbols <- c("BrIS-MS1" = "BrIS (BrScan-MS1)", "BrIS-Spectronaut" = "BrIS (Spectronaut)", "iRT-Spectronaut" = "iRT (Spectronaut)")
method_colors <- c("#31859B", "#EA700D", "#7E649E") 
names(method_colors) <- method_levels

df <- df %>%
  mutate(
    instrument = factor(instrument, levels = instrument_levels),
    value_name = factor(value_name, levels = c("RT", "areas"), labels = c("RT", "Intensity")),
    method = factor(method, levels = method_levels)
  )

create_figure <- function(df, colors, fontsize) {
  p <- ggplot(df, aes(x = instrument, y = cv)) +
    geom_boxplot(
      aes(fill = method, color = method),
      width = 0.6,
      linewidth = 0.5,
      alpha = 0.6,
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    geom_point(
      aes(fill = method, color = method),
      position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0),
      size = 0.5,
      alpha = 1.0
    ) +
    scale_fill_manual(values = colors, name = NULL, labels = method_symbols) +
    scale_color_manual(values = colors, name = NULL, labels = method_symbols) +
    scale_x_discrete(labels = instrument_symbols) +
    labs(
      x = "A: timsTOF Pro 1  B: timsTOF Pro 2  C: Exploris 480", 
      y = NULL
    ) +
    facet_wrap2(
      ~ value_name, 
      ncol = 2, 
      scales = 'free',
      strip.position = "left",
      labeller = as_labeller(
        c("RT" = "CV of RT (%)", "Intensity" = "CV of Intensity (%)")
        )
      ) +
    theme(
      plot.margin = margin(1,1,1,1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(fill = NA, colour = NA),
      strip.background = element_blank(), 
      strip.text = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      strip.placement = "outside",
      axis.line = element_line(colour = 'black', linewidth = 0.25),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = fontsize, colour = 'black'),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt")
    ) +
  facetted_pos_scales(
    y = list(
      value_name == "RT" ~ scale_y_continuous(breaks = seq(0, 4, 1),limits = c(0, 4)),
      value_name == "Intensity" ~ scale_y_continuous(breaks = seq(0, 12, 3), limits = c(0, 12))
    )
  )
  return(p)
}

p <- create_figure(df, method_colors, fontsize = unit(6.5, "pt"))
print(p)