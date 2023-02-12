library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

cluster_tpm <- read_delim("analysis_v3/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

Gbox_NC_df <- read_delim("data/Gbox_TF/Gbox.NC.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Gbox_PP_df <- read_delim("data/Gbox_TF/Gbox.PP.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Gbox_PP_df <- na.omit(Gbox_PP_df)

Gbox_NC_df$Source <- "NC"
Gbox_PP_df$Source <- "PP"

Gbox_TF_df <- rbind(Gbox_NC_df, Gbox_PP_df)

Gbox_scRNA <- left_join(Gbox_TF_df, cluster_tpm)
Gbox_scRNA_mat <- as.matrix(Gbox_scRNA[,5:ncol(Gbox_scRNA)])
rownames(Gbox_scRNA_mat) <- Gbox_scRNA$Gid
Gbox_scRNA_mat[rowMaxs(Gbox_scRNA_mat)<1,] <- NA
norm_Gbox_scRNA_mat <- Gbox_scRNA_mat / apply(Gbox_scRNA_mat, 1, function(x){return(max(x, na.rm = T))})

Gbox_NC_mat <- norm_Gbox_scRNA_mat[Gbox_TF_df$Source=="NC",]
Gbox_PP_mat <- norm_Gbox_scRNA_mat[Gbox_TF_df$Source=="PP",]

NC_scRNA_ht <- Heatmap(
  Gbox_NC_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  split = Gbox_NC_df$Symbol,
  show_row_names = T, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_DE_peak2C3/Gbox.NC.scRNA.pdf", width=8/inche_cm, height=6/inche_cm)
print(NC_scRNA_ht)
dev.off()

PP_scRNA_ht <- Heatmap(
  Gbox_PP_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  split = Gbox_PP_df$Symbol,
  show_row_names = T, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_DE_peak2C3/Gbox.PP.scRNA.pdf", width=8/inche_cm, height=15/inche_cm)
print(PP_scRNA_ht)
dev.off()
