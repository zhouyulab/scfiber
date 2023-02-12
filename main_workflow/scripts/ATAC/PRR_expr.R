library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

cluster_tpm <- read_delim("analysis_v3/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
CPM_TP <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

PRR_cotton <- read_delim("data/genome/HAU/PRR.cotton.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
PRR_info_df <- PRR_cotton[,c("Gid", "PRR")]
PRR_info_df <- PRR_info_df[order(PRR_info_df$PRR, PRR_info_df$Gid),]
PRR_scRNA <- left_join(PRR_info_df, cluster_tpm)
PRR_scRNA_mat <- as.matrix(PRR_scRNA[,3:ncol(PRR_scRNA)])
rownames(PRR_scRNA_mat) <- PRR_scRNA$Gid
PRR_scRNA_mat[rowMaxs(PRR_scRNA_mat)<1,] <- NA
norm_PRR_scRNA_mat <- PRR_scRNA_mat / apply(PRR_scRNA_mat, 1, function(x){return(max(x, na.rm = T))})

PRR_TP_df <- left_join(PRR_info_df, CPM_TP)
PRR_TP_mat <- as.matrix(PRR_TP_df[,3:ncol(PRR_TP_df)])
rownames(PRR_TP_mat) <- PRR_TP_df$Gid
PRR_TP_mat[rowMaxs(PRR_TP_mat)<1,] <- NA
norm_PRR_TP_mat <- PRR_TP_mat / apply(PRR_TP_mat, 1, function(x){return(max(x, na.rm = T))})

scRNA_ht <- Heatmap(
  norm_PRR_scRNA_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  split = PRR_info_df$PRR,
  show_row_names = T, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
TP_ht <- Heatmap(
  norm_PRR_TP_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  cluster_rows = F,
  split = PRR_info_df$PRR,
  show_row_names = T, 
  column_split = c(rep("WT", 9), rep("FL", 9)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_DE_peak2C3/PRR.scRNA.pdf", width=16/inche_cm, height=10/inche_cm)
print(scRNA_ht+TP_ht)
dev.off()
