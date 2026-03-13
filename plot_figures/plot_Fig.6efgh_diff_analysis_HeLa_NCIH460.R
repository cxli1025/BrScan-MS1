library(tidyverse)
library(limma)
library(VennDiagram)
library(ggrepel)
library(grid)

process_data <- function(df, selected_list) {
  df_clean <- df %>% 
    select(PG.ProteinGroups, contains("PG.Quantity")) %>%
    mutate(PG.ProteinGroups = str_split_i(PG.ProteinGroups, ";", 1)) %>%
    distinct(PG.ProteinGroups, .keep_all = TRUE)
  
  raw_cols <- colnames(df_clean)
  colnames(df_clean) <- c(
    "PG.ProteinGroups", 
    sapply(raw_cols[-1], function(x) {paste(str_split_i(x, "_", 3), str_split_i(x, "_", 4), sep = "_")})
  )
  
  df_clean <- df_clean %>% select(PG.ProteinGroups, all_of(selected_list))
  
  expr_mat <- as.matrix(df_clean[, -1])
  rownames(expr_mat) <- df_clean$PG.ProteinGroups
  
  gps <- case_when(
    grepl("^._[0-9]", colnames(expr_mat)) ~ "Ctrl",
    grepl("^.B_[0-9]", colnames(expr_mat)) ~ "BrIS",
    grepl("^.I_[0-9]", colnames(expr_mat)) ~ "iRT"
    )
  
  # 过滤：在每个组内缺失率需 <= 50%
  keep <- apply(expr_mat, 1, function(row) {
    group_na_rates <- tapply(row, gps, function(x) mean(is.na(x)))
    all(group_na_rates <= 0.5, na.rm = TRUE) 
  })
  expr_mat <- expr_mat[keep, ]
  
  # 组内中位数填充
  expr_filled <- t(apply(expr_mat, 1, function(row) {
    for (g in unique(gps)) {
      idx <- which(gps == g)
      if (any(is.na(row[idx]))) {
        m <- median(row[idx], na.rm = TRUE)
        row[idx][is.na(row[idx])] <- ifelse(is.na(m), 0, m)
      }
    }
    return(row)
  }))
  
  expr_log <- log2(expr_filled + 1)
  res_df <- as.data.frame(expr_log)
  res_df$PG.ProteinGroups <- rownames(expr_log)
  return(res_df)
}

# read  data
setwd("D:/project/04BrIS/three_cancer_single_cell")

hela_df <- read.csv("20260129_hela_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")
n460_df <- read.csv("20260129_n460_human_iRT_fasta_spectronaut_19_exclude_iRT.csv")

# select samples
hela_sel <- c("HB_1", "HB_2", "HB_4", "HB_6", "HB_8", "HI_1", "HI_2", "HI_5", "HI_7", "HI_8", "H_4", "H_5", "H_6", "H_7", "H_8")
n460_sel <- c("NB_2", "NB_3", "NB_4", "NB_7", "NB_8", "NI_3", "NI_4", "NI_6", "NI_7", "NI_8", "N_4", "N_5", "N_6", "N_7", "N_8")

hela <- process_data(hela_df, hela_sel)
n460 <- process_data(n460_df, n460_sel)

# 合并交集蛋白
combined <- inner_join(hela, n460, by = "PG.ProteinGroups")
data_mat <- as.matrix(combined %>% select(-PG.ProteinGroups))
rownames(data_mat) <- combined$PG.ProteinGroups

# 4. 构建 limma 差异分析模型
group_list <- factor(case_when(
  grepl("^HB_", colnames(data_mat)) ~ "HB",
  grepl("^HI_", colnames(data_mat)) ~ "HI",
  grepl("^H_",  colnames(data_mat)) ~ "H",
  grepl("^NB_", colnames(data_mat)) ~ "NB",
  grepl("^NI_", colnames(data_mat)) ~ "NI",
  grepl("^N_",  colnames(data_mat)) ~ "N"
))

design <- model.matrix(~0 + group_list)
colnames(design) <- levels(group_list)

fit <- lmFit(data_mat, design)

# 设置对比：N460 vs HeLa
contrast <- makeContrasts(
  Ctrl    = N - H,
  BrIS    = NB - HB,
  iRT     = NI - HI,
  levels  = design
)

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# 5. 提取差异分析结果
get_volcano_df <- function(fit_obj, coef_idx, p_thresh = 0.05, lfc_thresh = 0.58) {
  tt <- topTable(fit_obj, coef = coef_idx, number = Inf)
  tt$Protein <- rownames(tt)
  tt$diff <- "No.sig"
  tt$diff[tt$P.Value < p_thresh & tt$logFC > lfc_thresh] <- "Up"
  tt$diff[tt$P.Value < p_thresh & tt$logFC < -lfc_thresh] <- "Down"
  return(tt)
}

res_Control <- get_volcano_df(fit2, "Ctrl")
res_BrIS    <- get_volcano_df(fit2, "BrIS")
res_iRT     <- get_volcano_df(fit2, "iRT")

fontsize <- 6.5
# 6. 火山图绘制函数
plot_volcano <- function(df, title_text, p_thresh = 0.05, lfc_thresh = 0.58) {
  ggplot(df, aes(x = logFC, y = -log10(P.Value), color = diff)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_manual(values = c("Up" = "#B2182B", "Down" = "navy", "No.sig" = "grey")) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "black", linewidth = 0.25) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "black", linewidth = 0.25) +
    labs(
      title = title_text, 
      x = expression(log[2]~"(FC)"),
      y = expression(-log[10]~"(Pvalue)")
      ) +
    scale_x_continuous(limits = c(-5, 10), breaks = c(-5, 0, 5, 10), expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5), expand = expansion(mult = c(0.05, 0.05)))+
    theme(
      plot.margin = margin(1,1,1,1, "pt"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.25),
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
      axis.text.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 2, b = 0)),
      axis.text.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 2, l = 0)),
      axis.title.x = element_text(size = fontsize, colour = 'black', margin = margin(t = 1, b = 0)),
      axis.title.y = element_text(size = fontsize, colour = 'black', margin = margin(r = 1, l = 0)),
      plot.title = element_text(size = fontsize, colour = 'black', margin = margin(t = 0, b = 1), hjust = 0.5),
      legend.position = c(0.8, 0.2),
      legend.title = element_blank(),
      legend.text = element_text(size = 6.5, colour = 'black'),
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "pt")
    )
}

p1 <- plot_volcano(res_Control, "NCI-H460 vs HeLa (Ctrl)")
p2 <- plot_volcano(res_BrIS,    "NCI-H460 vs HeLa (BrIS)")
p3 <- plot_volcano(res_iRT,     "NCI-H460 vs HeLa (iRT)")

de_ctrl <- res_Control$Protein[res_Control$diff != "No.sig"]
de_bris <- res_BrIS$Protein[res_BrIS$diff != "No.sig"]
de_irt  <- res_iRT$Protein[res_iRT$diff != "No.sig"]

# save figure
current_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
output_path <- file.path(current_path, "Main_Figures", "Figure6")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

volcano_list <- list(
  "Fig.6e_Volcano_NC_vs_HC.pdf" = p1,
  "Fig.6f_Volcano_NB_vs_HB.pdf" = p2,
  "Fig.6g_Volcano_NI_vs_HI.pdf" = p3
)

for (file_name in names(volcano_list)) {
  print(volcano_list[[file_name]])
  ggsave(
    filename = file.path(output_path, file_name),
    plot = volcano_list[[file_name]],
    device = cairo_pdf,
    width = 1.6, 
    height = 1.7, 
    family = "Arial"
  )
}

venn_list <- list(Ctrl = de_ctrl, BrIS = de_bris, iRT = de_irt)
venn_obj <- venn.diagram(
  x = venn_list,
  filename = NULL,
  col = c("navy", "#D95F02", "#B2182B"),
  alpha = 1.0,
  fill = c("white", "white", "white"),
  lwd = 0.25, 
  cat.cex = fontsize/12,
  cex = fontsize/12,
  cat.dist = c(0.08, 0.08, 0.05),
  main = "Overlap of DEPs (NCI-H460 vs HeLa)",
  main.cex = fontsize/12,
  main.pos = c(0.5, 1.1)
)

grid.newpage()
grid.draw(venn_obj)

cairo_pdf(
  file.path(output_path, "Fig.6h_diff_proteins_venn_HeLa_NCIH460.pdf"), 
  width = 1.68, 
  height = 1.68,
  family = "Arial"
  )
grid.newpage()
grid.draw(venn_obj)
dev.off()