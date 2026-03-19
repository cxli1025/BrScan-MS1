library(ggplot2)
library(dplyr)
library(ggh4x)

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')
file <- 'Fig.4d_plasma_BrIS_iRT_samples_on_TIMS1_BrIS_iRT_intensity.csv'
df <- read.csv(file)
df$method <- factor(df$method, levels = c('BrIS (BrScan-MS1)', 'BrIS (Spectronaut)', 'iRT (Spectronaut)'))

method_labels <- c(
  'BrIS (BrScan-MS1)' = 'BrScan-MS1',
  'BrIS (Spectronaut)' = 'Spectronaut',
  'iRT (Spectronaut)'  = 'Spectronaut'
)

# BrIS
colors1 <- c(
  "#31859B", "#51A3B8", "#1D5667", 
  "#EA700D", "#F4A460", "#A34D09", 
  "#7E649E", "#A391BD", "#5F4B75"  
)

# iRT 
colors2 <- c(
  "#31859B", "#51A3B8", "#1D5667",  
  "#EA700D", "#F4A460", "#A34D09",  
  "#7E649E", "#A391BD", "#5F4B75",  
  "#76933C", "#9BBB59"
)

names(colors1) <- c("BrIS-1", "BrIS-2", "BrIS-3", "BrIS-4", "BrIS-5", "BrIS-6", "BrIS-7", "BrIS-8", "BrIS-9")
names(colors2) <- c("RT-pep a", "RT-pep b", "RT-pep c", "RT-pep d", "RT-pep e", "RT-pep f", 
                    "RT-pep g", "RT-pep h", "RT-pep i", "RT-pep k", "RT-pep l")
all_colors <- c(colors1, colors2)

summary_df <- df %>%
  group_by(method, peptide) %>%
  summarise(
    mean_area = mean(areas, na.rm = TRUE),
    sd_area   = sd(areas, na.rm = TRUE),
    n         = n(),
    se_area   = sd_area / sqrt(n),
    .groups   = 'drop'
  )

fontsize <- 6.5

p <- ggplot() +
  geom_bar(
    data = summary_df,
    aes(x = peptide, y = mean_area, color = peptide),
    stat = "identity", width = 0.6, linewidth = 0.25, fill = "white"
  ) +
  geom_point(
    data = df,
    aes(x = peptide, y = areas, color = peptide),
    size = 0.5, alpha = 1.0,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = peptide, ymin = mean_area - se_area, ymax = mean_area + se_area),
    width = 0.3, linewidth = 0.5, color = "black", alpha = 1.0
  ) +
  facet_wrap(~ method, ncol = 3, scales = "free_x",labeller = labeller(method = method_labels)) +
  scale_color_manual(values = all_colors) +
  scale_y_continuous(breaks = seq(0, 8, 2), limits = c(0, 8.5)) +
  labs(x = 'title', y = expression(log[10]~(Intensity))) +
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
  )+
  facetted_pos_scales(
    x = list(
      method %in% c("BrIS (BrScan-MS1)", "BrIS (Spectronaut)") ~ scale_x_discrete(
        labels = function(x) gsub("^BrIS-", "", x)),
      method == "iRT (Spectronaut)" ~ scale_x_discrete(
        labels = function(x) gsub("^RT-pep ", "", x))
      )
    )
print(p)