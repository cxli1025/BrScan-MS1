library(tidyverse)
library(ComplexHeatmap)
library(circlize)

ht_opt$message = FALSE

process_single_dataset <- function(df, selected_list) {
  df_cleaned <- df %>% 
    select(PG.ProteinGroups, contains("PG.Quantity")) %>%
    mutate(PG.ProteinGroups = str_split_i(PG.ProteinGroups, ";", 1)) %>%
    distinct(PG.ProteinGroups, .keep_all = TRUE)
  
  raw_colnames <- colnames(df_cleaned)
  new_colnames <- c("PG.ProteinGroups", 
                    sapply(raw_colnames[-1], function(x) paste(str_split_i(x, "_", 3), str_split_i(x, "_", 4), sep = "_")))
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

calc_corr <- function(df_res) {
  mat <- df_res %>% select(-PG.ProteinGroups) %>% as.matrix()
  return(cor(mat, method = "pearson"))
}

setwd("D:/project/04BrIS/three_cancer_single_cell")

hela_df <- read.csv("20260129_hela_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")
n460_df <- read.csv("20260129_n460_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")
m231_df <- read.csv("20260129_m231_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")

hela_sel <- c("H_4", "H_5", "H_6", "H_7", "H_8", "HB_1", "HB_2", "HB_4", "HB_6", "HB_8", "HI_1", "HI_2", "HI_5", "HI_7", "HI_8")
n460_sel <- c("N_4", "N_5", "N_6", "N_7", "N_8", "NB_2", "NB_3", "NB_4", "NB_7", "NB_8", "NI_3", "NI_4", "NI_6", "NI_7", "NI_8")
m231_sel <- c("M_3", "M_4", "M_6", "M_7", "M_8", "MB_2", "MB_3", "MB_5", "MB_6", "MB_8", "MI_1", "MI_2", "MI_4", "MI_5", "MI_7")

mat_h <- calc_corr(process_single_dataset(hela_df, hela_sel))
mat_n <- calc_corr(process_single_dataset(n460_df, n460_sel))
mat_m <- calc_corr(process_single_dataset(m231_df, m231_sel))

fontsize <- 6.5
grid_h <- unit(2, "mm") 
grid_w <- unit(3, "mm")
anno_height <- unit(2, "mm") 

col_fun = colorRamp2(c(0.5, 0.75, 1.0), c("navy", "#F7F7F7", "#B2182B"))

cell_line_cols = c(
  "HeLa"       = "#31859B", 
  "NCI-H460"   = "#EA700D", 
  "MDA-MB-231" = "#7E649E"  
)

cond_cols = c(
  "Ctrl" = "#66C2A5",
  "BrIS" = "#8DA0CB",
  "iRT"  = "#E78AC3"
)

lgd_cell = Legend(title = "Cell Line", labels = names(cell_line_cols), 
                  legend_gp = gpar(fill = cell_line_cols), 
                  ncol = 1,
                  title_gp = gpar(fontsize = fontsize), 
                  labels_gp = gpar(fontsize = fontsize),
                  grid_height = grid_h, grid_width = grid_w)

lgd_cond = Legend(title = "Condition", labels = names(cond_cols), 
                  legend_gp = gpar(fill = cond_cols), 
                  ncol = 1, 
                  title_gp = gpar(fontsize = fontsize), 
                  labels_gp = gpar(fontsize = fontsize),
                  grid_height = grid_h, grid_width = grid_w)

lgd_corr = Legend(title = "", col_fun = col_fun, 
                  direction = "horizontal", 
                  at = seq(0.5, 1.0, 0.25), 
                  labels_gp = gpar(fontsize = fontsize),
                  grid_height = grid_h, legend_width = unit(1.7, "cm"), border = NA)

ha_h = HeatmapAnnotation(CellLine = rep("HeLa", 15), 
                         Condition = factor(rep(c("Ctrl", "BrIS", "iRT"), each = 5), levels = names(cond_cols)),
                         col = list(CellLine = cell_line_cols, Condition = cond_cols), 
                         show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = anno_height)

ha_n = HeatmapAnnotation(CellLine = rep("NCI-H460", 15), 
                         Condition = factor(rep(c("Ctrl", "BrIS", "iRT"), each = 5), levels = names(cond_cols)),
                         col = list(CellLine = cell_line_cols, Condition = cond_cols), 
                         show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = anno_height)

ha_m = HeatmapAnnotation(CellLine = rep("MDA-MB-231", 15), 
                         Condition = factor(rep(c("Ctrl", "BrIS", "iRT"), each = 5), levels = names(cond_cols)),
                         col = list(CellLine = cell_line_cols, Condition = cond_cols), 
                         show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = anno_height)

h_labels <- c(paste0("HC", 1:5), paste0("HB", 1:5), paste0("HI", 1:5))
n_labels <- c(paste0("NC", 1:5), paste0("NB", 1:5), paste0("NI", 1:5))
m_labels <- c(paste0("MC", 1:5), paste0("MB", 1:5), paste0("MI", 1:5))

ht_list = 
  Heatmap(mat_h, name = "H", col = col_fun, top_annotation = ha_h, 
          cluster_rows = FALSE, cluster_columns = FALSE, 
          show_column_names = TRUE, column_names_side = "bottom",
          column_labels = h_labels, 
          column_names_gp = gpar(fontsize = fontsize),
          show_row_names = FALSE, show_heatmap_legend = FALSE) + 
  Heatmap(mat_n, name = "N", col = col_fun, top_annotation = ha_n, 
          cluster_rows = FALSE, cluster_columns = FALSE, 
          show_column_names = TRUE, column_names_side = "bottom",
          column_labels = n_labels, 
          column_names_gp = gpar(fontsize = fontsize),
          show_row_names = FALSE, show_heatmap_legend = FALSE) + 
  Heatmap(mat_m, name = "M", col = col_fun, top_annotation = ha_m, 
          cluster_rows = FALSE, cluster_columns = FALSE, 
          show_column_names = TRUE, column_names_side = "bottom",
          column_labels = m_labels, 
          column_names_gp = gpar(fontsize = fontsize),
          show_row_names = FALSE, show_heatmap_legend = FALSE)

pd = packLegend(lgd_cell, lgd_cond, lgd_corr, 
                direction = "vertical", 
                gap = unit(2, "mm"))

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure6")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

cairo_pdf(
  file.path(output_path, "Fig.6d_three_tumor_single_cell_heapmap.pdf"),
  width = 4.9, 
  height = 1.5, 
  family = "Arial")

draw(ht_list, 
     ht_gap = unit(1, "mm"),
     annotation_legend_list = pd,
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     merge_legends = TRUE,
     padding = unit(c(1,1,1,1), "pt") 
)

dev.off()

draw(ht_list, 
     ht_gap = unit(1, "mm"),
     annotation_legend_list = pd,
     merge_legends = TRUE,
     padding = unit(c(1,1,1,1), "pt"))