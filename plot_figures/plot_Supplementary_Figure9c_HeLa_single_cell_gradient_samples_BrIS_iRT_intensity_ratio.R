library(ggplot2)
library(stringr)
library(tidyverse)

f1 <- 'D:/project/04BrIS/hela_single_cell_gradient_samples/BrIS/new_ms1/areas_ratios.csv'
f2 <- 'D:/project/04BrIS/hela_single_cell_gradient_samples/iRT/areas_ratios.csv'
df1 <- read.csv(f1)
df1 <- df1 %>% mutate(peptide = 'BrIS peptides')
df2 <- read.csv(f2)
df2 <- df2 %>% mutate(peptide = 'iRT peptides')
df <- rbind(df1, df2)

sample_order <- c('0.1pg', '0.5pg', '1pg', '5pg', '10pg', '50pg', '100pg')
df$sample <- factor(df$sample, levels = sample_order)
df$peptide <- factor(df$peptide, levels = c("BrIS peptides", "iRT peptides"))
theory_df <- unique(df[, c("sample", "theory_ratio")])

colors <- c("BrIS peptides" = "#31859B", "iRT peptides" = "#EA700D")
fontsize <- 6.5
p <- ggplot(df, aes(x = sample, y = ratio)) + 
  geom_boxplot(
    aes(fill = peptide, color = peptide), 
    width = 0.6, alpha = 0.7, lwd = 0.5, 
    position = position_dodge(width = 0.8, preserve = "single"),
    outlier.shape = NA
    ) + 
  geom_point(
    aes(fill = peptide, color = peptide),
    position = position_dodge(width = 0.8, preserve = "single"),
    size = 0.5,
    alpha = 1
  ) +
  geom_errorbar(
    data = theory_df, 
    aes(x = sample, ymin = theory_ratio, ymax = theory_ratio),
    color = 'darkgrey', 
    width = 1.0,        
    linetype = "dotted", 
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x = NULL, y = expression(log[2](Intensity~ratio))) + 
  scale_x_discrete(labels = function(x) gsub("pg", "", x)) +
  scale_y_continuous(limits = c(-9, 9), breaks = seq(-9, 9, 3)) +
  scale_fill_manual(values = colors, name = NULL) +
  scale_color_manual(values = colors, name = NULL) +
  theme(
    plot.margin  = margin(1, 1, 1, 1, "pt"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
    axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
    axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = 'black', linewidth = 0.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = fontsize, colour = "black"),
    legend.title = element_text(size = fontsize, colour = "black"),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt")
  )
print(p)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Supplementary_Figures", "Figure9")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
ggsave(
  file.path(output_path, "Supplementary_Figure9c_HeLa_single_cell_gradient_samples_BrIS_iRT_intensity_ratio.pdf"),
  plot = p,
  width = 4.35,
  height = 1.9,
  family = 'Arial',
  device = cairo_pdf
)