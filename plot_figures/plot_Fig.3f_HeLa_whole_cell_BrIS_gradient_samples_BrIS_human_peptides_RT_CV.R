library(ggplot2)
library(dplyr)
library(stringr)
library(ggbreak)

data1 <- read.csv("D:/project/04BrIS/HeLa_BrIS_gradient_samples_BrIS_peptides_RT_CV.csv")
data2 <- read.csv("D:/project/04BrIS/HeLa_BrIS_gradient_samples_human_peptides_RT_CV.csv")
data2 <- data2 %>%
  select(Peptide, CV, Method, Category)
data <- rbind(data1, data2)
colnames(data) <- tolower(colnames(data))
data$category <- factor(data$category, levels = c("BrIS_1", "BrIS_2", "Human"))
data$method <- factor(data$method, levels = c("BrScan-MS1", "PEAKS Online"))

low_cv_data <- data %>% filter(cv <= 10)

# Calculate summary stats
summary_stats <- low_cv_data %>%
    group_by(category, method) %>%
    summarise(max_cv = max(cv, na.rm = TRUE), .groups = "drop")

colors <- c("#31859B", "#EA700D")

fontsize <- unit(6.5, "pt")
p <- ggplot() +
    geom_boxplot(
      data = low_cv_data, 
      aes(x = category, y = cv, fill = method, colour = method),
      width = 0.3, linewidth = 0.5, alpha = 0.7, outlier.shape = NA
      ) +
    geom_point(
        data = summary_stats, 
        aes(x = category, y = max_cv, fill = method, color = method), 
        size = 0.5, alpha = 1.0
    ) +
    scale_fill_manual(values = colors, name = NULL) +
    scale_color_manual(values = colors, name = NULL) +
    scale_x_discrete(labels = c("BrIS", "BrIS", "Human")) +
    scale_y_break(c(0.7, 0.8), scales = c(1, 1), space = 0.3, expand = expansion(mult = c(0, 0)))+
    scale_y_continuous(breaks = c(0, 0.7, 0.8, 10.0), limits = c(0, 10.5)) +
    labs(x = NULL, y = "CV of RT (%)") +
    theme(
        plot.margin = margin(0,0,0,0, "pt"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        axis.line.x = element_line(colour = "black", linewidth = 0.25),
        axis.line.y.left = element_line(colour = "black", linewidth = 0.25),
        axis.line.y.right = element_blank(),
        axis.text.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
        axis.text.y.left = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
        axis.text.y.right = element_blank(),
        axis.title.y = element_text(size = fontsize, colour = "black", margin = margin(r = 0, l = 0), angle = 90),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.25),
        axis.ticks.y.left = element_line(colour = "black", linewidth = 0.25),
        axis.ticks.y.right = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.text = element_text(size = fontsize, colour = "black"),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
    )

print(p)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure3")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
ggsave(
  file.path(output_path, "Fig.3f_HeLa_whole_cell_BrIS_gradient_samples_BrIS_human_peptides_RT_CV.pdf"),
  plot = last_plot(),
  width = 2.3,
  height = 1.8,
  family = 'Arial',
  device = cairo_pdf
  )