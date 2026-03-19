library(tidyverse)
library(stringr)
library(patchwork)

count_proteins <- function(df, selected_list, gradient_map, gradient_order, category_order) {
  proteins <- df %>%
    select(ProteinGroups = PG.ProteinGroups, contains("PG.Quantity")) %>%
    rename_with(~ str_split_i(.x, "_", 3), contains("PG.Quantity")) %>%
    select(ProteinGroups, all_of(intersect(selected_list, names(.)))) %>%
    mutate(ProteinGroups = str_split_i(ProteinGroups, ";", 1)) %>%
    distinct(ProteinGroups, .keep_all = TRUE)

  sample_cols <- setdiff(names(proteins), "ProteinGroups")
  threshold <- length(sample_cols) * 0.5

  proteins_filtered <- proteins %>%
    rowwise() %>%
    mutate(na_count = sum(is.na(c_across(all_of(sample_cols))))) %>%
    ungroup() %>%
    filter(na_count <= threshold) %>%
    select(-na_count)

  stats <- proteins_filtered %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "Quantity") %>%
    filter(!is.na(Quantity)) %>%
    group_by(sample) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      gradient = gradient_map[sample],
      category = case_when(
        sample %in% c("A2") ~ "No spike-in",
        sample %in% c("2", "3", "C1", "D1", "6", "7", "8") ~ "BrIS peptides spike-in",
        sample %in% c("18", "19", "20", "B1", "22", "23", "24") ~ "iRT peptides spike-in",
      )
    )

  stats$gradient <- factor(stats$gradient, levels = gradient_order)
  stats$category <- factor(stats$category, levels = category_order)

  return(stats)
}

count_peptides <- function(df, selected_list, gradient_map, gradient_order, category_order) {
  peptides <- df %>%
    select(PeptideSequence = PEP.StrippedSequence, contains("PEP.Quantity")) %>%
    rename_with(~ str_split_i(.x, "_", 3), contains("PEP.Quantity")) %>%
    select(PeptideSequence, all_of(intersect(selected_list, names(.)))) %>%
    distinct(PeptideSequence, .keep_all = TRUE)

  sample_cols <- setdiff(names(peptides), "PeptideSequence")
  threshold <- length(sample_cols) * 0.5

  peptides_filtered <- peptides %>%
    rowwise() %>%
    mutate(na_count = sum(is.na(c_across(all_of(sample_cols))))) %>%
    ungroup() %>%
    filter(na_count <= threshold) %>%
    select(-na_count)

  stats <- peptides_filtered %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "Quantity") %>%
    filter(!is.na(Quantity)) %>%
    group_by(sample) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      gradient = gradient_map[sample],
      category = case_when(
        sample %in% c("A2") ~ "No spike-in",
        sample %in% c("2", "3", "C1", "D1", "6", "7", "8") ~ "BrIS peptides spike-in",
        sample %in% c("18", "19", "20", "B1", "22", "23", "24") ~ "iRT peptides spike-in",
      )
    )

  stats$gradient <- factor(stats$gradient, levels = gradient_order)
  stats$category <- factor(stats$category, levels = category_order)

  return(stats)
}

create_figure <- function(stats_data, colors, fontsize, level) {
  p <- ggplot(stats_data, aes(x = gradient, y = count, fill = category)) +
    geom_bar(
      stat = "identity",
      width = 0.6,
      position = position_dodge(width = 0.8),
      color = "black",
      linewidth = 0.5
    ) +
    scale_fill_manual(values = colors, name = NULL) +
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = NA),
      axis.text.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.line = element_line(colour = "black", linewidth = 0.5),
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
  return(p)
}

setwd('C:/Users/chenxin/Desktop/BrIS_project/source_data')
file <- 'Supplementary_Figure9d_HeLa_single_cell_BrIS_iRT_gradient_samples_results.csv'
df <- read.csv(file)

selected_samples <- c("A2", "2", "3", "C1", "D1", "6", "7", "8", "18", "19", "20", "B1", "22", "23", "24")
dilution_gradient <- c(
  "A2" = "0",
  "2" = "0.1", "3" = "0.5", "C1" = "1", "D1" = "5", "6" = "10", "7" = "50", "8" = "100",
  "18" = "0.1", "19" = "0.5", "20" = "1", "B1" = "5", "22" = "10", "23" = "50", "24" = "100"
)

gradient_order <- c("0", "0.1", "0.5", "1", "5", "10", "50", "100")
category_order <- c("No spike-in", "BrIS peptides spike-in", "iRT peptides spike-in")

protein_stats_df <- count_proteins(df, selected_samples, dilution_gradient, gradient_order, category_order)
peptide_stats_df <- count_peptides(df, selected_samples, dilution_gradient, gradient_order, category_order)

fontsize <- 6.5
colors <- c("No spike-in" = "white", "BrIS peptides spike-in" = "#8DA0CB", "iRT peptides spike-in" = "#66C2A5")

p_protein <- create_figure(protein_stats_df, colors, fontsize) +
  scale_y_continuous(
    limits = c(0, 3500),
    breaks = seq(0, 3000, 1000),
    labels = c("0", "1", "2", "3"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = expression(paste("Protein Count (", "\u00D7", 10^3, ")")))

p_peptide <- create_figure(peptide_stats_df, colors, fontsize) +
  scale_y_continuous(
    limits = c(0, 31000),
    breaks = seq(0, 30000, 10000),
    labels = c("0", "1", "2", "3"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = expression(paste("Peptide Count (", "\u00D7", 10^4, ")")))

p <- p_protein / p_peptide

print(p)