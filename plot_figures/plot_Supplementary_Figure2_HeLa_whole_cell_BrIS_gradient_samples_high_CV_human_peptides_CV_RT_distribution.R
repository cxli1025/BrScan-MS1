library(ggplot2)
library(stringr)
library(dplyr)
library(patchwork)

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')
data <- read.csv('Fig.3f_HeLa_whole_cell_BrIS_gradient_samples_human_peptides_RT_CV.csv')
colnames(data) <- tolower(colnames(data))
high_cv_data <- data %>% filter(cv > 10)
total_count <- nrow(high_cv_data)
fontsize = 6.5
p_hist <- ggplot(high_cv_data, aes(x = cv)) +
  geom_histogram(binwidth = 10, boundary = 0, fill = "skyblue", color = "white") +
  stat_bin(binwidth = 10, boundary = 0, geom = "text", 
           aes(label = after_stat(count)), 
           vjust = -0.5, size = fontsize/.pt, color = "black") +
  scale_x_continuous(
    limits = c(0, 240),
    breaks = seq(0, 240, 40),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = paste0("CV Distribution of High-CV Peptides (CV > 10%, n = ", total_count, ")"),
    x = 'CV of Retention Time (%)',
    y = "Count"
  ) +
  theme(
    plot.margin = margin(1, 2, 5, 2, "pt"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.text.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
    axis.text.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
    axis.title.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
    axis.title.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
    plot.title = element_text(size = fontsize, colour = "black", margin = margin(t = 0, b = 10), hjust = 0.5)
  )

# ---
target_cols <- paste0("x", 2:9)
high_cv_data$median_rt <- apply(high_cv_data[, target_cols], 1, median, na.rm = TRUE)

count_0_6  <- sum(high_cv_data$median_rt <= 6, na.rm = TRUE)
count_0_10 <- sum(high_cv_data$median_rt <= 10, na.rm = TRUE)
pct_0_6  <- round((count_0_6 / total_count) * 100, 1)
pct_0_10 <- round((count_0_10 / total_count) * 100, 1)

p_median <- ggplot(high_cv_data, aes(x = median_rt)) +
  geom_histogram(binwidth = 1, boundary = 0, fill = "skyblue", color = "white") +
  stat_bin(binwidth = 1, boundary = 0, geom = "text", 
           aes(label = after_stat(count)), 
           vjust = -0.5, size = fontsize/.pt, color = "black") +
  geom_vline(xintercept = 6, linetype = "dashed", color = "#E41A1C", linewidth = 0.5) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "#377EB8", linewidth = 0.5) +
  annotate("text", x = 3.5, y = Inf, label = paste0("0-6min: ", pct_0_6, "%"), 
           vjust = 1.0, size = fontsize / .pt, color = "#E41A1C") +
  annotate("text", x = 13, y = Inf, label = paste0("0-10min: ", pct_0_10, "%"), 
           vjust = 1.0, size = fontsize / .pt, color = "#377EB8") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = paste0("RT Distribution of High-CV Peptides (CV > 10%, n = ", total_count, ")"),
    x = "Median Retention Time (min)",
    y = "Count"
  ) +
  theme(
    plot.margin = margin(5,2,1,2, "pt"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = NA),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.text.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
    axis.text.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
    axis.title.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
    axis.title.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
    plot.title = element_text(size = fontsize, colour = "black", margin = margin(t = 0, b = 10), hjust = 0.5)
  )

p <- p_hist / p_median
print(p)