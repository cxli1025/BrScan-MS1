library(tidyverse)
library(FactoMineR)
library(factoextra)

process_data <- function(df, selected_list) {
  df_cleaned <- df %>% 
    select(PG.ProteinGroups, contains("PG.Quantity")) %>%
    mutate(PG.ProteinGroups = str_split_i(PG.ProteinGroups, ";", 1)) %>%
    distinct(PG.ProteinGroups, .keep_all = TRUE)
  
  new_colnames <- c(
    "PG.ProteinGroups", 
    sapply(colnames(df_cleaned)[-1], function(x) {paste(str_split_i(x, "_", 3), str_split_i(x, "_", 4), sep = "_")})
  )
  colnames(df_cleaned) <- new_colnames
  
  df_cleaned <- df_cleaned %>% select(PG.ProteinGroups, all_of(selected_list))
  
  expr_matrix <- as.matrix(df_cleaned[, -1])
  rownames(expr_matrix) <- df_cleaned$PG.ProteinGroups
  
  is_complete <- rowSums(is.na(expr_matrix)) == 0
  expr_filtered <- expr_matrix[is_complete, ]
  expr_final <- log2(expr_filtered + 1)
  
  res_df <- as.data.frame(expr_final)
  res_df$PG.ProteinGroups <- rownames(expr_final)
  return(res_df)
}

# load data
setwd("D:/project/04BrIS/three_cancer_single_cell")

hela_df <- read.csv("20260129_hela_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")
n460_df <- read.csv("20260129_n460_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")
m231_df <- read.csv("20260129_m231_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")

hela_sel <- c("HB_1", "HB_2", "HB_4", "HB_6", "HB_8", "HI_1", "HI_2", "HI_5", "HI_7", "HI_8", "H_4", "H_5", "H_6", "H_7", "H_8")
n460_sel <- c("NB_2", "NB_3", "NB_4", "NB_7", "NB_8", "NI_3", "NI_4", "NI_6", "NI_7", "NI_8", "N_4", "N_5", "N_6", "N_7", "N_8")
m231_sel <- c("MB_2", "MB_3", "MB_5", "MB_6", "MB_8", "MI_1", "MI_2", "MI_4", "MI_5", "MI_7", "M_3", "M_4", "M_6", "M_7", "M_8")

hela <- process_data(hela_df, hela_sel)
n460 <- process_data(n460_df, n460_sel)
m231 <- process_data(m231_df, m231_sel)
df <- hela %>% inner_join(m231, by = "PG.ProteinGroups") %>% inner_join(n460, by = "PG.ProteinGroups")

pca_input <- t(df %>% select(-PG.ProteinGroups))
res.pca <- PCA(pca_input, scale.unit = TRUE, ncp = 5, graph = FALSE)

condition_groups <- case_when(
  grepl("^._[0-9]", rownames(pca_input)) ~ "C",
  grepl("^.B_[0-9]", rownames(pca_input)) ~ "B",
  grepl("^.I_[0-9]", rownames(pca_input)) ~ "I"
  
)

cell_line_groups <- case_when(
  grepl("^H", rownames(pca_input)) ~ "H",
  grepl("^N", rownames(pca_input)) ~ "N",
  grepl("^M", rownames(pca_input)) ~ "M"
)

colors <- c(
  "HC" = "#000080", "HB" = "#D95F02", "HI" = "#B2182B", 
  "NC" = "#1B9E77", "NB" = "#7570B3", "NI" = "#E6AB02",
  "MC" = "#E7298A", "MB" = "#66A61E", "MI" = "#A6761D"
)

groups <- factor(
  paste0(cell_line_groups, condition_groups),
  levels = c("HC", "HB", "HI", "NC", "NB", "NI","MC", "MB", "MI")
)

p <- fviz_pca_ind(
  res.pca,
  geom.ind = "point",      
  pointsize = 0.5,
  pointshape = 19,
  col.ind = groups,
  palette = colors,
  addEllipses = TRUE,
  ellipse.level = 0.95,
  ellipse.type = "t",
  repel = TRUE
) +
  coord_cartesian(xlim = c(-55, 55)) + 
  scale_x_continuous(breaks = c(-50, -25, 0, 25, 50))+
  theme(
    plot.margin = margin(1,1,1,1, "pt"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.25),
    axis.ticks = element_line(colour = 'black', linewidth = 0.25),
    axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = -5, b = 0)),
    axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = -5, l = 0)),
    axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
    axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
    plot.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 6.5, colour = 'black'),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt")
  ) +
  guides(color = guide_legend(override.aes = list(size = 0.5, shape = 19)))
print(p)

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure6")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

ggsave(
  file.path(output_path, "Fig.6c_three_tumor_single_cell_pca.pdf"),
  plot = p,
  width = 1.9,    
  height = 1.5,    
  device = cairo_pdf,
  family = "Arial"
)