library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

cluster_tpm <- read_delim("analysis_v3/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

CPM_TP <- read_delim("../CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


blast2go <- read_delim("analysis_v3/annotation/blast2GO/blast2go.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_info <- blast2go[grep("TCP[0-9]+", blast2go$Description),]
TCP_info$Gid <- sapply(strsplit(TCP_info$SeqName, "[.]"), function(x){return(x[1])})

TCP_indx <- regexpr("TCP[0-9\\-like]+", TCP_info$Description)
TCP_len <- attr(TCP_indx,"match.length")
TCP_indx <- as.integer(TCP_indx)
TCP_info_df <- data.frame(Gid=TCP_info$Gid, TCP=substr(TCP_info$Description, TCP_indx, TCP_indx+TCP_len))
TCP_info_df <- TCP_info_df[!duplicated(TCP_info_df$Gid),]

TCP_scRNA <- left_join(TCP_info_df, cluster_tpm)
TCP_scRNA_mat <- as.matrix(TCP_scRNA[,3:ncol(TCP_scRNA)])
rownames(TCP_scRNA_mat) <- TCP_scRNA$Gid
expr_TCP_scRNA_mat <- TCP_scRNA_mat[rowMaxs(TCP_scRNA_mat, na.rm = T)>1,]
norm_TCP_scRNA_mat <- expr_TCP_scRNA_mat / apply(expr_TCP_scRNA_mat, 1, function(x){return(max(x, na.rm = T))})
expr_TCP_info_df <- TCP_info_df[TCP_info_df$Gid%in%rownames(norm_TCP_scRNA_mat),]
scRNA_ht <- Heatmap(
  norm_TCP_scRNA_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  split = expr_TCP_info_df$TCP,
  show_row_names = T, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
scRNA_ht

inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP.scRNA.blast2go.pdf", width=8/inche_cm, height=10/inche_cm)
print(scRNA_ht)
dev.off()

TCP_family <- read_delim("data/genome/HAU/TCP_family.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_family$TCP <- sapply(strsplit(TCP_family$Name, "[-]"), function(x){return(x[1])})
names(TCP_family)[3] <- "Gid"

TCP_scRNA <- left_join(TCP_family, cluster_tpm)
TCP_scRNA_mat <- as.matrix(TCP_scRNA[,5:ncol(TCP_scRNA)])
rownames(TCP_scRNA_mat) <- TCP_scRNA$Gid
expr_TCP_scRNA_mat <- TCP_scRNA_mat[rowMaxs(TCP_scRNA_mat, na.rm = T)>1,]
norm_TCP_scRNA_mat <- expr_TCP_scRNA_mat / apply(expr_TCP_scRNA_mat, 1, function(x){return(max(x, na.rm = T))})
expr_TCP_family <- TCP_family[TCP_family$Gid%in%rownames(norm_TCP_scRNA_mat),]

TCP_CPM_TP <- left_join(TCP_family, CPM_TP)
TCP_CPM_TP_mat <- as.matrix(TCP_CPM_TP[,5:ncol(TCP_CPM_TP)])
rownames(TCP_CPM_TP_mat) <- TCP_CPM_TP$Gid
expr_TCP_CPM_TP_mat <- TCP_CPM_TP_mat[rownames(TCP_CPM_TP_mat) %in% rownames(norm_TCP_scRNA_mat),]
norm_TCP_CPM_TP_mat <- expr_TCP_CPM_TP_mat / apply(expr_TCP_CPM_TP_mat, 1, function(x){return(max(x, na.rm = T))})

WT_expr_TCP_CPM_TP_mat <- expr_TCP_CPM_TP_mat[,1:(ncol(expr_TCP_CPM_TP_mat)/2)]
FL_expr_TCP_CPM_TP_mat <- expr_TCP_CPM_TP_mat[,(ncol(expr_TCP_CPM_TP_mat)/2+1):ncol(expr_TCP_CPM_TP_mat)]
WT_log2FC_mat <- log2(WT_expr_TCP_CPM_TP_mat[,2:ncol(WT_expr_TCP_CPM_TP_mat)] / WT_expr_TCP_CPM_TP_mat[,1:(ncol(WT_expr_TCP_CPM_TP_mat)-1)])
WT_log2FC_mat[is.infinite(WT_log2FC_mat)] <- NA
FL_log2FC_mat <- log2(FL_expr_TCP_CPM_TP_mat[,2:ncol(FL_expr_TCP_CPM_TP_mat)] / FL_expr_TCP_CPM_TP_mat[,1:(ncol(FL_expr_TCP_CPM_TP_mat)-1)])
FL_log2FC_mat[is.infinite(FL_log2FC_mat)] <- NA
merge_log2FC_mat <- cbind(WT_log2FC_mat, FL_log2FC_mat)

scRNA_ht <- Heatmap(
  norm_TCP_scRNA_mat, 
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
  split = expr_TCP_family$TCP,
  show_row_names = T, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

RNA_TP_ht <- Heatmap(
  norm_TCP_CPM_TP_mat, 
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
  split = expr_TCP_family$TCP,
  show_row_names = T, 
  column_split = c(rep("WT", 9), rep("FL", 9)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
RNA_TP_ht

inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP.scRNA.pdf", width=16/inche_cm, height=10/inche_cm)
print(scRNA_ht+RNA_TP_ht)
dev.off()

RNA_TP_log2FC_ht <- Heatmap(
  merge_log2FC_mat, 
  name="log2FC",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  split = expr_TCP_family$TCP,
  show_row_names = T, 
  column_split = c(rep("WT", 8), rep("FL", 8)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP.TP.log2FC.pdf", width=10/inche_cm, height=10/inche_cm)
print(RNA_TP_log2FC_ht)
dev.off()



select_TCP <- c("GhTCP14a", "GhTCP14b", "GhTCP14c", "GhTCP7a", "GhTCP21")
select_TCP_cpm_TP <- TCP_CPM_TP[TCP_CPM_TP$TCP%in%select_TCP,]
melt_select_TCP_cpm_TP <- melt(select_TCP_cpm_TP, id.vars = c("Name", "The original ID", "Gid", "TCP"), variable.name = "label", value.name = "CPM")
melt_select_TCP_cpm_TP$Sample <- sapply(strsplit(as.character(melt_select_TCP_cpm_TP$label), "_"), function(x){return(x[1])})
melt_select_TCP_cpm_TP$Time <- sapply(strsplit(as.character(melt_select_TCP_cpm_TP$label), "_"), function(x){return(x[2])})
time_df <- data.frame(Time=c("n48h", "n36h", "n24h", "n12h", "0h", "12h", "24h", "36h", "48h"), TimeInt=c(-48, -36, -24, -12, 0, 12, 24, 36, 48))
melt_select_TCP_cpm_TP <- left_join(melt_select_TCP_cpm_TP, time_df)
p <- ggplot(melt_select_TCP_cpm_TP, aes(x=TimeInt, y=CPM, color=Sample)) +
  geom_line(size=0.3) +
  geom_point(size=0.7) +
  annotate("point", x=0, y=0, color=NA) +
  scale_color_manual(values = c("WT"="red", "FL"="blue")) +
  scale_x_continuous(breaks = time_df$TimeInt, labels = time_df$Time) +
  facet_grid(Gid~., scales = "free_y") +
  labs(x="Time") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size=5),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        legend.title = element_text(family="ArialMT", size=6),
        legend.text = element_text(family="ArialMT", size=5),
        legend.key.size = unit(3, "mm"),
        panel.grid = element_blank()
  )
ggsave("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP7_14_21.ExprTp.pdf", p, width = 7, height = 12, units = "cm")
