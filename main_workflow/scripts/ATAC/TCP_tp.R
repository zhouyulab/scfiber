library(readr)
library(SingleCellExperiment)
library(ComplexHeatmap)

TCP_gene <- read_delim("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_like_gene <- read_delim("analysis_v3/ATAC/ATAC_kmer_cluster2gene/ACCCT_motif/ACCCT_motif.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

TCP_gene_df <- data.frame(Gid=unique(c(TCP_gene$Gid, TCP_like_gene$Gid)))
TCP_gene_df$IsTCP <- TCP_gene_df$Gid %in% TCP_gene$Gid
TCP_gene_df$IsTCPLike <- TCP_gene_df$Gid %in% TCP_like_gene$Gid

TP_FPKM <- read_delim("analysis_v3/RNA_seq/taco_merge_expr/TP.FPKM.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TP_FPKM$`5DPA` <- NULL

load("analysis_v3/final_cluster/mergeSCE.RData")
cnt_mat <- counts(mergeSCE)
rm(mergeSCE)
gc()
norm_mat <- 1000 * t(t(cnt_mat) / colSums(cnt_mat))

sc_expr_num <- rowSums(cnt_mat>1)
hist(log10(sc_expr_num), breaks = 50)
expr_cnt_mat <- cnt_mat[sc_expr_num>100,]

overlap_sc_mat <- expr_cnt_mat[rownames(expr_cnt_mat) %in% TCP_gene_df$Gid,]
rm(expr_cnt_mat)
overlap_sc_norm_mat <- norm_mat[rownames(norm_mat) %in% rownames(overlap_sc_mat),]

cluster_df <- read_delim("analysis_v3/final_cluster/cluster.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
cluster_df$BarcodeIndx <- 1:nrow(cluster_df)
cluster_df <- cluster_df[cluster_df$Cluster%in%c("C1", "C2", "C3"),]

scRNA_pt <- read_delim("analysis_v3/Trajectory/pt.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
scRNA_pt$PseudoTimeBin <- as.integer(scRNA_pt$PseudoTime*2) / 2
scRNA_pt$Rep <- sapply(strsplit(scRNA_pt$Sample, "[.]"), function(x){return(x[2])})
scRNA_pt$Sample <- sapply(strsplit(scRNA_pt$Sample, "[.]"), function(x){return(x[1])})
scRNA_pt <- inner_join(scRNA_pt, cluster_df)
max(scRNA_pt$PseudoTimeBin)
all_tp_bin <- sort(unique(scRNA_pt$PseudoTimeBin))
all_gene_li <- rownames(overlap_sc_norm_mat)
WT_TCP_tp_mat <- matrix(NA, nrow = length(all_gene_li), ncol = length(all_tp_bin))
FL_TCP_tp_mat <- matrix(NA, nrow = length(all_gene_li), ncol = length(all_tp_bin))
rownames(WT_TCP_tp_mat) <- all_gene_li
rownames(FL_TCP_tp_mat) <- all_gene_li
colnames(WT_TCP_tp_mat) <- sprintf("WT.%s", all_tp_bin)
colnames(FL_TCP_tp_mat) <- sprintf("FL.%s", all_tp_bin)
for(tp_indx in 1:length(all_tp_bin)){
  tp_bin <- all_tp_bin[tp_indx]
  tmp_cell_df <- scRNA_pt[scRNA_pt$PseudoTimeBin==tp_bin,]
  WT_cell_df <- tmp_cell_df[tmp_cell_df$Sample=="WT",]
  FL_cell_df <- tmp_cell_df[tmp_cell_df$Sample=="FL",]
  if(nrow(WT_cell_df)>1){
    tmp_WT_expr_mat <- overlap_sc_norm_mat[,WT_cell_df$BarcodeIndx]
    WT_TCP_tp_mat[, tp_indx] <- rowMeans(tmp_WT_expr_mat)
  }
  if(nrow(FL_cell_df)>1){
    tmp_FL_expr_mat <- overlap_sc_norm_mat[,FL_cell_df$BarcodeIndx]
    FL_TCP_tp_mat[, tp_indx] <- rowMeans(tmp_FL_expr_mat)
  }
}

merge_tcp_mat <- cbind(WT_TCP_tp_mat, FL_TCP_tp_mat)
norm_merge_tcp_mat <- merge_tcp_mat / rowMaxs(merge_tcp_mat, na.rm = T)
TCP_ht <- Heatmap(
  norm_merge_tcp_mat,
  name = "Norm. Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
  cluster_columns=FALSE,
  column_split=c(rep("WT", ncol(merge_tcp_mat)/2), rep("FL", ncol(merge_tcp_mat)/2)),
  show_row_names=FALSE,
  clustering_method_rows = "ward.D2",
  row_title_gp = gpar(foutsize = 5),
  column_title_gp = gpar(foutsize = 5),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5), 
    gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

tcp_info <- data.frame(Gid=rownames(norm_merge_tcp_mat))
tcp_info <- left_join(tcp_info, TCP_gene_df)
col_li <- c("TRUE"="black", "FALSE"="white")
ha <- rowAnnotation(
  TCP=tcp_info$IsTCP,
  TCP_like=tcp_info$IsTCPLike,
  annotation_name_gp=gpar(fontsize = 5),
  width=unit(1, "cm"),
  border=FALSE,
  col=list(
    TCP=col_li,
    TCP_like=col_li
  )
)

inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP.Trajectory.pdf", width=14/inche_cm, height=8/inche_cm)
print(TCP_ht + ha)
dev.off()

