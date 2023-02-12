library(dplyr)
library(readr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

args = commandArgs(TRUE)
if (length(args) == 3) {
  f_merker_gene <- args[1]
  f_tp_fpkm <- args[2]
  f_plot <- args[3]
} else {
  q()
}

# f_merker_gene <- "/home/sqreb/ias/data/CottonSingleCell/analysis/stat/marker_gene/Marker.gene.tsv"
# f_tp_fpkm <- "/home/sqreb/ias/data/CottonSingleCell/analysis/RNA_seq/merge_expr/TP.FPKM.tsv"

marker_df <- read_delim(f_merker_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
tp_fpkm <- read_delim(f_tp_fpkm, "\t", escape_double = FALSE, trim_ws = TRUE)

tp_fpkm_mat <- as.matrix(tp_fpkm[,2:ncol(tp_fpkm)])
ngs_expr_gid <- tp_fpkm$Gid[apply(tp_fpkm_mat, 1, min)>0]

compute_fc <- function(mat){
  mat[is.na(mat)] <- 0
  x <- t(apply(mat, 1, function(x) sort(x, decreasing = TRUE)))
  return(x[,1]/x[,2])
}

filter_marker_df <- function(block){
  WT_mat <- as.matrix(block[,21:25])
  FL_mat <- as.matrix(block[,26:30])
  block$WT_fc <- compute_fc(WT_mat)
  clu <- block$Cluster[1]
  if(clu=="C3"){
    block$FL_fc <- NA
    block$EnrichScore <- block$WT_fc
  }else{
    block$FL_fc <- compute_fc(FL_mat)
    block$EnrichScore <- apply(block[,c("WT_fc", "FL_fc")], 1, min)
  }
  if(clu!="C4"){
  block <- block[block[,sprintf("WT_%s_ExprRatio", clu)]>0.1,]
  }
  block <- block[order(block$EnrichScore, decreasing = TRUE),]
  return(block[1:10,])
}

#marker_df <- marker_df[marker_df$Gid %in% ngs_expr_gid,]
filtered_marker_df <- plyr::ddply(marker_df, "Cluster", filter_marker_df)
filtered_marker_df[21:30,26:30] <- NA
cutoff <- 2
filtered_marker_df <- filtered_marker_df[filtered_marker_df$EnrichScore > cutoff | filtered_marker_df$Cluster=="C4",]
filtered_marker_df <- filtered_marker_df[!is.na(filtered_marker_df$Cluster),]
WT_mat <- as.matrix(filtered_marker_df[,21:25])
scaled_WT_mat <- WT_mat / apply(WT_mat, 1, function(x) max(x, na.rm = TRUE))

FL_mat <- as.matrix(filtered_marker_df[,26:30])
scaled_FL_mat <- FL_mat / apply(FL_mat, 1, function(x){if(all(is.na(x)))return(NA)else{return( max(x, na.rm = TRUE))}})

scaled_tpm_mat <- cbind(scaled_WT_mat, scaled_FL_mat)
rownames(scaled_tpm_mat) <- filtered_marker_df$Gid

marker_tp_fpkm <- data.frame(Gid=rownames(scaled_tpm_mat))
marker_tp_fpkm <- left_join(marker_tp_fpkm, tp_fpkm, by="Gid")

marker_tp_fpkm_mat <- as.matrix(marker_tp_fpkm[,2:ncol(marker_tp_fpkm)])
rownames(marker_tp_fpkm_mat) <- marker_tp_fpkm$Gid
marker_tp_fpkm_ratio_mat <- marker_tp_fpkm_mat / apply(marker_tp_fpkm_mat, 1, max)

ht <- Heatmap(scaled_tpm_mat, name="Norm. Expr.",
              col =  colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_dend_width = unit(5, "mm"),
              column_dend_height = unit(5, "mm"),
              show_row_dend = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 6),
                title_gp = gpar(fontsize = 6),
                grid_width = unit(2, "mm"),
                grid_height = unit(2, "mm"))
)

ht_tp <- Heatmap(marker_tp_fpkm_ratio_mat, name="Norm. Expr.",
                 col =  colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
                 row_names_gp = gpar(fontsize = 6),
                 column_names_gp = gpar(fontsize = 6),
                 row_title_gp = gpar(fontsize = 6),
                 row_dend_width = unit(5, "mm"),
                 column_dend_height = unit(5, "mm"),
                 show_row_dend = FALSE,
                 show_row_names = TRUE,
                 cluster_columns = FALSE,
                 # cluster_rows = FALSE,
                 heatmap_legend_param = list(
                   labels_gp = gpar(fontsize = 6),
                   title_gp = gpar(fontsize = 6),
                   grid_width = unit(2, "mm"),
                   grid_height = unit(2, "mm"))
)

inche_cm=2.54
pdf(f_plot, width=20/inche_cm, height=10/inche_cm, family="ArialMT")
draw(ht + ht_tp)
dev.off()
