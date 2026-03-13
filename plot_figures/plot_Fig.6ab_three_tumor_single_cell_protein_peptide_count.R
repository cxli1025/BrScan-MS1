library(tidyverse)
library(ggplot2)
library(ggprism)
library(rstatix)

# count proteins
count_proteins <- function(df, selected_list) {
  df_clean <- df %>%
    select(PG.ProteinGroups, contains("PG.Quantity")) %>%
    mutate(PG.ProteinGroups = str_split_i(PG.ProteinGroups, ";", 1)) %>%
    distinct(PG.ProteinGroups, .keep_all = TRUE)
  
  long_df <- df_clean %>%
    pivot_longer(cols = contains("PG.Quantity"), names_to = "raw_name", values_to = "Quantity") %>%
    mutate(sample_id = paste(str_split_i(raw_name, "_", 3), str_split_i(raw_name, "_", 4), sep = "_")) %>%
    filter(sample_id %in% selected_list) %>%
    mutate(condition = case_when(
      grepl("^.B_[0-9]", sample_id) ~ "BrIS",
      grepl("^.I_[0-9]", sample_id) ~ "iRT",
      grepl("^._[0-9]", sample_id) ~ "Ctrl"
    ))
  
  valid_proteins <- long_df %>%
    group_by(PG.ProteinGroups, condition) %>%
    summarise(missing_rate = mean(is.na(Quantity)), .groups = "drop") %>%
    group_by(PG.ProteinGroups) %>%
    filter(any(missing_rate <= 0.5)) %>%
    pull(PG.ProteinGroups) %>%
    unique()
  
  stats <- long_df %>%
    filter(PG.ProteinGroups %in% valid_proteins) %>%
    filter(!is.na(Quantity)) %>%
    group_by(sample_id, condition) %>%
    summarise(count = n_distinct(PG.ProteinGroups), .groups = "drop") %>%
    mutate(group = str_sub(sample_id, 1, 1), level = "Protein")
  return(stats)
}

# count peptides
count_peptides <- function(df, selected_list) {
  df_clean <- df %>%
    select(PEP.StrippedSequence, contains("PEP.Quantity")) %>%
    distinct(PEP.StrippedSequence, .keep_all = TRUE)
  
  long_df <- df_clean %>%
    pivot_longer(cols = contains("PEP.Quantity"), names_to = "raw_name", values_to = "Quantity") %>%
    mutate(sample_id = paste(str_split_i(raw_name, "_", 3), str_split_i(raw_name, "_", 4), sep = "_")) %>%
    filter(sample_id %in% selected_list) %>%
    mutate(condition = case_when(
      grepl("^.B_[0-9]", sample_id) ~ "BrIS",
      grepl("^.I_[0-9]", sample_id) ~ "iRT",
      grepl("^._[0-9]", sample_id) ~ "Ctrl"
    ))
  
  valid_peptides <- long_df %>%
    group_by(PEP.StrippedSequence, condition) %>%
    summarise(missing_rate = mean(is.na(Quantity)), .groups = "drop") %>%
    group_by(PEP.StrippedSequence) %>%
    filter(any(missing_rate <= 0.5)) %>%
    pull(PEP.StrippedSequence) %>%
    unique()
  
  stats <- long_df %>%
    filter(PEP.StrippedSequence %in% valid_peptides) %>%
    filter(!is.na(Quantity)) %>%
    group_by(sample_id, condition) %>%
    summarise(count = n_distinct(PEP.StrippedSequence), .groups = "drop") %>%
    mutate(group = str_sub(sample_id, 1, 1), level = "Peptide")
  return(stats)
}

# stats_analysis
stats_analysis <- function(hela_df, hela_sel, n460_df, n460_sel, m231_df, m231_sel, level) {
  if (level == "peptide") {
    count_df <- rbind(count_peptides(hela_df, hela_sel), count_peptides(n460_df, n460_sel), count_peptides(m231_df, m231_sel))
  } else if (level == "protein") {
    count_df <- rbind(count_proteins(hela_df, hela_sel), count_proteins(n460_df, n460_sel), count_proteins(m231_df, m231_sel))
  }
  count_df <- count_df %>%
    mutate(
      group = case_when(group == "H" ~ "HeLa", group == "N" ~ "NCI-H460", group == "M" ~ "MDA-MB-231", TRUE ~ group),
      condition = factor(condition, levels = c("Ctrl", "BrIS", "iRT")),
      group = factor(group, levels = c("HeLa", "NCI-H460", "MDA-MB-231"))
    )
  summary_df <- count_df %>%
    group_by(condition, group) %>%
    summarise(mean_val = mean(count), se_val = sd(count) / sqrt(n()), .groups = "drop")
  return(list(count_df = count_df, summary_df = summary_df))
}

create_figure <- function(count_df, summary_df, comparison_subject, colors, fontsize) {
  p <- ggplot(data = count_df, aes(x = condition, y = count)) +
    geom_bar(
      data = summary_df,
      aes(x = condition, y = mean_val, fill = condition, color = condition),
      stat = "identity", width = 0.6, linewidth = 0.5, alpha = 0.6
    ) +
    geom_point(
      aes(fill = condition, color = condition),
      position = position_jitter(width = 0.2), size = 1.0, alpha = 1.0
    ) +
    geom_errorbar(
      data = summary_df,
      aes(x = condition, ymin = mean_val - se_val, ymax = mean_val + se_val),
      width = 0.3, linewidth = 0.5, color = "black", inherit.aes = FALSE
    ) +
    facet_wrap(~group, scales = "free_x") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt"), 
      panel.grid = element_blank(), 
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = NA), 
      strip.background = element_rect(fill = "#E0E0E0", colour = "#E0E0E0", linewidth = 0.5),
      strip.text = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 2)),
      axis.line = element_line(colour = "black", linewidth = 0.5), 
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.text.x = element_text(size = fontsize, colour = "black", margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
      axis.title.y = element_text(size = fontsize, colour = "black", margin = margin(r = 2, l = 0)),
      panel.spacing = unit(3, "pt"),
      legend.position = "none"
    )
  return(p)
}

# load data
setwd("D:/project/04BrIS/three_cancer_single_cell")
hela_df <- read.csv("20260129_hela_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")
n460_df <- read.csv("20260129_n460_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")
m231_df <- read.csv("20260129_m231_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")

hela_sel <- c("HB_1", "HB_2", "HB_4", "HB_6", "HB_8", "HI_1", "HI_2", "HI_5", "HI_7", "HI_8", "H_4", "H_5", "H_6", "H_7", "H_8")
n460_sel <- c("NB_2", "NB_3", "NB_4", "NB_7", "NB_8", "NI_3", "NI_4", "NI_6", "NI_7", "NI_8", "N_4", "N_5", "N_6", "N_7", "N_8")
m231_sel <- c("MB_2", "MB_3", "MB_5", "MB_6", "MB_8", "MI_1", "MI_2", "MI_4", "MI_5", "MI_7", "M_3", "M_4", "M_6", "M_7", "M_8")

protein_stats <- stats_analysis(hela_df, hela_sel, n460_df, n460_sel, m231_df, m231_sel, "protein")
peptide_stats <- stats_analysis(hela_df, hela_sel, n460_df, n460_sel, m231_df, m231_sel, "peptide")

colors <- c("Ctrl" = "#31859B", "BrIS" = "#EA700D", "iRT" = "#7E649E")
fontsize <- 6.5
comparison_subject <- list(c("Ctrl", "BrIS"), c("BrIS", "iRT"), c("Ctrl", "iRT"))

# --- Protein ---
df_p_protein <- protein_stats$count_df %>%
  group_by(group) %>%
  t_test(count ~ condition, comparisons = comparison_subject) %>%
  add_significance() %>%
  mutate(y.position = rep(c(4200, 4200, 4600), 3)) 

p_protein <- create_figure(protein_stats$count_df, protein_stats$summary_df, comparison_subject, colors, fontsize) +
  add_pvalue(
    df_p_protein, 
    label = "p.adj.signif",
    label.size = fontsize / .pt,
    tip.length = 0.05,
    bracket.shorten = 0.1, 
    vjust = 0.3
  ) +
  scale_y_continuous(
    limits = c(0, 5000),
    breaks = seq(0, 4000, 1000),
    labels = c("0", "1", "2", "3", "4"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = expression(paste("Protein Count (", "\u00D7", 10^3, ")")))

# --- Peptide ---
df_p_peptide <- peptide_stats$count_df %>%
  group_by(group) %>%
  t_test(count ~ condition, comparisons = comparison_subject) %>%
  add_significance() %>%
  mutate(y.position = rep(c(33000, 33000, 36500), 3))

p_peptide <- create_figure(peptide_stats$count_df, peptide_stats$summary_df, comparison_subject, colors, fontsize) +
  add_pvalue(
    df_p_peptide, 
    label = "p.adj.signif",
    label.size = fontsize / .pt,
    tip.length = 0.05,
    bracket.shorten = 0.1, 
    vjust = 0.3
  ) +
  scale_y_continuous(
    limits = c(0, 39500),
    breaks = seq(0, 30000, 10000),
    labels = c("0", "1", "2", "3"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = expression(paste("Peptide Count (", "\u00D7", 10^4, ")")))

print(p_protein)
print(p_peptide)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure6")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

ggsave(
  file.path(output_path, "Fig.6a_three_tumor_single_cell_protein_count.pdf"),
  plot = p_protein,
  width = 3.3,
  height = 1.5,
  family = 'Arial',
  device = cairo_pdf
)

ggsave(
  file.path(output_path, "Fig.6b_three_tumor_single_cell_peptide_count.pdf"),
  plot = p_peptide,
  width = 3.3,
  height = 1.5,
  family = 'Arial',
  device = cairo_pdf
)