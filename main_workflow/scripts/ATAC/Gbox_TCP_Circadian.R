library(readr)
library(ComplexHeatmap)
library(circlize)

CPM_merge <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
CircadianGene_Group <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CircadianGene.Group.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
CircadianGene_Group$CircadianState <- sapply(strsplit(CircadianGene_Group$GroupID, "[_]"), function(x){return(x[1])})
CircadianGene_Group$CircadianType <- sapply(strsplit(CircadianGene_Group$GroupID, "[_]"), function(x){return(x[2])})
CircadianGene_Group$GroupID <- NULL
motif_hit_df <- read_delim("analysis_v3/ATAC/motif/summary/MotifExprTargetGene.HitMat.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

Gbox <- data.frame(Gid=motif_hit_df$Gid[motif_hit_df$Gbox])

cluster_tpm <- read_delim("analysis_v3/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Gbox_topGO <- read_delim("analysis_v3/ATAC/motif/Gbox_motif/topGO/topGO.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

Gbox_CPM <- left_join(Gbox, CPM_merge)
Gbox_CPM_mat <- as.matrix(Gbox_CPM[,2:ncol(Gbox_CPM)])
Gbox_CPM_mat[apply(Gbox_CPM_mat, 1, max)<=1,] <- NA
rownames(Gbox_CPM_mat) <- Gbox_CPM$Gid

cluster_tpm <- left_join(Gbox, cluster_tpm)
cluster_tpm_mat <- as.matrix(cluster_tpm[,2:ncol(cluster_tpm)])
cluster_tpm_mat[apply(cluster_tpm_mat, 1, function(x){return(max(x, na.rm = T))})<=1,] <- NA
rownames(cluster_tpm_mat) <- cluster_tpm$Gid

NA_indx <- is.na(Gbox_CPM_mat[,1]) | is.na(cluster_tpm_mat[,1])
Gbox_CPM_mat <- Gbox_CPM_mat[!NA_indx,]
cluster_tpm_mat <- cluster_tpm_mat[!NA_indx,]
Gbox <- Gbox[!NA_indx,]

norm_cluster_tpm_mat <- cluster_tpm_mat / apply(cluster_tpm_mat, 1, function(x){return(max(x, na.rm = T))})

Gbox_CPM_WT_mat <- Gbox_CPM_mat[,1:9]
Gbox_CPM_FL_mat <- Gbox_CPM_mat[,10:18]
Gbox_CPM_WT_diff_mat <- t(apply(Gbox_CPM_WT_mat, 1, function(x){return(diff(x)>0)}))
Gbox_CPM_FL_diff_mat <- t(apply(Gbox_CPM_FL_mat, 1, function(x){return(diff(x)>0)}))

Gbox_df <- data.frame(Gid=rownames(Gbox_CPM_mat))
Gbox_df <- left_join(Gbox_df, CircadianGene_Group)
Gbox_df$CircadianState[is.na(Gbox_df$CircadianState)] <- "NC"
Gbox_df$CircadianType[is.na(Gbox_df$CircadianType)] <- "NC"

Gbox_CPM_WT_diff_mat_merge <- cbind(t(apply(Gbox_CPM_WT_mat, 1, diff)), t(apply(Gbox_CPM_FL_mat, 1, diff)))
colnames(Gbox_CPM_WT_diff_mat_merge) <- c(
  sprintf("%s.%s", colnames(Gbox_CPM_WT_mat)[1:8], colnames(Gbox_CPM_WT_mat)[2:9]),
  sprintf("%s.%s", colnames(Gbox_CPM_FL_mat)[1:8], colnames(Gbox_CPM_FL_mat)[2:9])
)
scaled_Gbox_CPM_WT_diff_mat_merge <- t(scale(t(Gbox_CPM_WT_diff_mat_merge)))

Gbox_CPM_WT_diff_mat_merge_bool <- Gbox_CPM_WT_diff_mat_merge>0
diff_dist <- dist(Gbox_CPM_WT_diff_mat_merge_bool)
diff_expr_dist <- dist(scaled_Gbox_CPM_WT_diff_mat_merge)
scRNA_dist <- dist(norm_cluster_tpm_mat)
diff_dist[is.na(diff_dist)] <- max(diff_dist, na.rm = T)
diff_expr_dist[is.na(diff_expr_dist)] <- max(diff_expr_dist, na.rm = T)
scRNA_dist[is.na(scRNA_dist)] <- max(scRNA_dist, na.rm = T)

merge_dist <- diff_dist + diff_expr_dist / 2 + scRNA_dist / 10

clu <- hclust(merge_dist, "ward.D2")

PRR_cotton <- read_delim("data/genome/HAU/PRR.cotton.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
circadian <- read_delim("data/genome/HAU/circadian.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
interest_gene <- data.frame(Gid=PRR_cotton$Gid, GeneSymbol=PRR_cotton$PRR)
interest_gene <- rbind(interest_gene, circadian)
interest_gene <- interest_gene[!duplicated(interest_gene$Gid),]
interest_gene_df <- left_join(data.frame(Gid=rownames(norm_cluster_tpm_mat)), interest_gene)
interest_gene_df$Label <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$GeneSymbol)
all(rownames(norm_cluster_tpm_mat) == interest_gene_df$Gid)
all(rownames(scaled_Gbox_CPM_WT_diff_mat_merge) == interest_gene_df$Gid)

subindx <- which(!is.na(interest_gene_df$GeneSymbol))
label_li <- interest_gene_df$Label[subindx]

norm_scRNA_ht <- Heatmap(
  norm_cluster_tpm_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  show_row_names = FALSE, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

scaled_TP_ht <- Heatmap(
  scaled_Gbox_CPM_WT_diff_mat_merge, 
  name="z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  row_order = clu$order,
  cluster_columns = FALSE, 
  show_row_names = F, 
  split = sprintf("%s.%s", Gbox_df$CircadianState, Gbox_df$CircadianType),
  column_split = c(rep("WT", 8), rep("FL", 8)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

ha <- rowAnnotation(
  link = row_anno_link(
    at = subindx,
    labels = label_li,
    labels_gp=gpar(fontsize = 6),
    lines_gp=gpar(lwd = 0.5),
    link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6))
  )
x <- scaled_TP_ht+norm_scRNA_ht+ha
save(x, file = "analysis_v3/ATAC/ATAC_DE_peak2C3/ATAC_DE_peak_target_gene.Gbox.Circadian.RData")
inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_DE_peak2C3/ATAC_DE_peak_target_gene.Gbox.Circadian.pdf", width=15/inche_cm, height=8/inche_cm)
print(scaled_TP_ht+norm_scRNA_ht+ha)
dev.off()
write_tsv(Gbox_df, "analysis_v3/ATAC/ATAC_DE_peak2C3/ATAC_DE_peak_target_gene.Gbox.Circadian.tsv")


MYB_G4 <- data.frame(Gid=motif_hit_df$Gid[motif_hit_df$MYB_G4])
MYB_G4_CPM <- left_join(MYB_G4, CPM_merge)
MYB_G4_CPM_mat <- as.matrix(MYB_G4_CPM[,2:ncol(MYB_G4_CPM)])
MYB_G4_CPM_mat[apply(MYB_G4_CPM_mat, 1, max)<=1,] <- NA
rownames(MYB_G4_CPM_mat) <- MYB_G4_CPM$Gid

cluster_tpm <- left_join(MYB_G4, cluster_tpm)
cluster_tpm_mat <- as.matrix(cluster_tpm[,2:ncol(cluster_tpm)])
cluster_tpm_mat[apply(cluster_tpm_mat, 1, function(x){return(max(x, na.rm = T))})<=1,] <- NA

NA_indx <- is.na(MYB_G4_CPM_mat[,1]) | is.na(cluster_tpm_mat[,1])
MYB_G4_CPM_mat <- MYB_G4_CPM_mat[!NA_indx,]
cluster_tpm_mat <- cluster_tpm_mat[!NA_indx,]
MYB_G4 <- MYB_G4[!NA_indx,]

norm_cluster_tpm_mat <- cluster_tpm_mat / apply(cluster_tpm_mat, 1, function(x){return(max(x, na.rm = T))})

MYB_G4_CPM_WT_mat <- MYB_G4_CPM_mat[,1:9]
MYB_G4_CPM_FL_mat <- MYB_G4_CPM_mat[,10:18]
MYB_G4_CPM_WT_diff_mat <- t(apply(MYB_G4_CPM_WT_mat, 1, function(x){return(diff(x)>0)}))
MYB_G4_CPM_FL_diff_mat <- t(apply(MYB_G4_CPM_FL_mat, 1, function(x){return(diff(x)>0)}))


MYB_G4_df <- data.frame(Gid=rownames(MYB_G4_CPM_mat))
WtWeakCircadianTransNum_li <- do.call(rbind, apply(MYB_G4_CPM_WT_diff_mat, 1, compute_circadian_trans_num))
MYB_G4_df$WtWeakCircadianTransType <- WtWeakCircadianTransNum_li$MaxType
MYB_G4_df$WtWeakCircadianTransNum <- WtWeakCircadianTransNum_li$MaxLen
FlWeakCircadianTransNum_li <-do.call(rbind, apply(MYB_G4_CPM_FL_diff_mat, 1, compute_circadian_trans_num))
MYB_G4_df$FlWeakCircadianTransType <- FlWeakCircadianTransNum_li$MaxType
MYB_G4_df$FlWeakCircadianTransNum <- FlWeakCircadianTransNum_li$MaxLen

MYB_G4_df$CircadianState <- "NC"
MYB_G4_df$CircadianState[(MYB_G4_df$WtWeakCircadianTransNum>=3 | MYB_G4_df$FlWeakCircadianTransNum>=3) & (MYB_G4_df$WtWeakCircadianTransType==MYB_G4_df$FlWeakCircadianTransType)] <- "NotEnriched"
MYB_G4_df$CircadianState[MYB_G4_df$Gid %in% CircadianGene_Group$Gid] <- "Enriched"
MYB_G4_df$CircadianType <- MYB_G4_df$WtWeakCircadianTransType
MYB_G4_df$CircadianType[MYB_G4_df$FlWeakCircadianTransNum>MYB_G4_df$WtWeakCircadianTransNum] <- MYB_G4_df$FlWeakCircadianTransType[MYB_G4_df$FlWeakCircadianTransNum>MYB_G4_df$WtWeakCircadianTransNum]
MYB_G4_df$CircadianType[MYB_G4_df$CircadianState=="NC"] <- 0

MYB_G4_CPM_WT_diff_mat_merge <- cbind(t(apply(MYB_G4_CPM_WT_mat, 1, diff)), t(apply(MYB_G4_CPM_FL_mat, 1, diff)))
colnames(MYB_G4_CPM_WT_diff_mat_merge) <- c(
  sprintf("%s.%s", colnames(MYB_G4_CPM_WT_mat)[1:8], colnames(MYB_G4_CPM_WT_mat)[2:9]),
  sprintf("%s.%s", colnames(MYB_G4_CPM_FL_mat)[1:8], colnames(MYB_G4_CPM_FL_mat)[2:9])
)
scaled_MYB_G4_CPM_WT_diff_mat_merge <- t(scale(t(MYB_G4_CPM_WT_diff_mat_merge)))

MYB_G4_CPM_WT_diff_mat_merge_bool <- MYB_G4_CPM_WT_diff_mat_merge>0
diff_dist <- dist(MYB_G4_CPM_WT_diff_mat_merge_bool)
diff_expr_dist <- dist(scaled_MYB_G4_CPM_WT_diff_mat_merge)
scRNA_dist <- dist(norm_cluster_tpm_mat)
diff_dist[is.na(diff_dist)] <- max(diff_dist, na.rm = T)
diff_expr_dist[is.na(diff_expr_dist)] <- max(diff_expr_dist, na.rm = T)
scRNA_dist[is.na(scRNA_dist)] <- max(scRNA_dist, na.rm = T)

merge_dist <- diff_dist + diff_expr_dist / 2 + scRNA_dist / 10

clu <- hclust(merge_dist, "ward.D2")

norm_scRNA_ht <- Heatmap(
  norm_cluster_tpm_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  show_row_names = FALSE, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

scaled_TP_ht <- Heatmap(
  scaled_MYB_G4_CPM_WT_diff_mat_merge, 
  name="z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  row_order = clu$order,
  cluster_columns = FALSE, 
  show_row_names = F, 
  split = sprintf("%s.%s", MYB_G4_df$CircadianState, MYB_G4_df$CircadianType),
  column_split = c(rep("WT", 8), rep("FL", 8)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
scaled_TP_ht
table(sprintf("%s.%s", MYB_G4_df$CircadianState, MYB_G4_df$CircadianType))

TCP <- read_delim("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_like <- read_delim("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_like_motif/TCP_like_motif.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_df <- data.frame(Gid=unique(c(TCP$Gid, TCP_like$Gid)))
TCP_df$IsTcpMotif <- TCP_df$Gid %in% TCP$Gid
TCP_df$IsTcpLikeMotif <- TCP_df$Gid %in% TCP_like$Gid

TCP_CPM <- left_join(TCP_df, CPM_merge)
expr_TCP_CPM <- TCP_CPM[apply(TCP_CPM[,4:ncol(TCP_CPM)], 1, max)>1,]
expr_TCP_CPM_mat <- as.matrix(expr_TCP_CPM[,4:ncol(expr_TCP_CPM)])
rownames(expr_TCP_CPM_mat) <- expr_TCP_CPM$Gid

expr_TCP_CPM_WT_mat <- expr_TCP_CPM_mat[,1:9]
expr_TCP_CPM_FL_mat <- expr_TCP_CPM_mat[,10:18]
expr_TCP_CPM_WT_diff_mat <- t(apply(expr_TCP_CPM_WT_mat, 1, function(x){return(diff(x)>0)}))
expr_TCP_CPM_FL_diff_mat <- t(apply(expr_TCP_CPM_FL_mat, 1, function(x){return(diff(x)>0)}))


ExprTCP_df <- data.frame(Gid=rownames(expr_TCP_CPM_mat))
WtWeakCircadianTransNum_li <- apply(expr_TCP_CPM_WT_diff_mat, 1, compute_circadian_trans_num)
ExprTCP_df$WtWeakCircadianTransType <- sapply(WtWeakCircadianTransNum_li, function(x){return(x[1])})
ExprTCP_df$WtWeakCircadianTransNum <- sapply(WtWeakCircadianTransNum_li, function(x){return(x[2])})
FlWeakCircadianTransNum_li <- t(apply(expr_TCP_CPM_FL_diff_mat, 1, compute_circadian_trans_num))
ExprTCP_df$FlWeakCircadianTransType <- FlWeakCircadianTransNum_li[,1]
ExprTCP_df$FlWeakCircadianTransNum <- FlWeakCircadianTransNum_li[,2]
ExprTCP_df$WtWeakCircadianTransNum[is.na(ExprTCP_df$WtWeakCircadianTransNum)] <- 0

ExprTCP_df$CircadianState <- "NC"
ExprTCP_df$CircadianState[(ExprTCP_df$WtWeakCircadianTransNum>=4 | ExprTCP_df$FlWeakCircadianTransNum>=4) & (ExprTCP_df$WtWeakCircadianTransType==ExprTCP_df$FlWeakCircadianTransType)] <- "NotEnriched"
ExprTCP_df$CircadianState[ExprTCP_df$Gid %in% CircadianGene_Group$Gid] <- "Enriched"
ExprTCP_df$CircadianType <- ExprTCP_df$WtWeakCircadianTransType
ExprTCP_df$CircadianType[ExprTCP_df$FlWeakCircadianTransNum>ExprTCP_df$WtWeakCircadianTransNum] <- ExprTCP_df$FlWeakCircadianTransType[ExprTCP_df$FlWeakCircadianTransNum>ExprTCP_df$WtWeakCircadianTransNum]
ExprTCP_df$CircadianType[ExprTCP_df$CircadianState=="NC"] <- 0
ExprTCP_df <- left_join(ExprTCP_df, TCP_df)

expr_TCP_CPM_WT_diff_mat_merge <- cbind(t(apply(expr_TCP_CPM_WT_mat, 1, diff)), t(apply(expr_TCP_CPM_FL_mat, 1, diff)))
colnames(expr_TCP_CPM_WT_diff_mat_merge) <- c(
  sprintf("%s.%s", colnames(expr_TCP_CPM_WT_mat)[1:8], colnames(expr_TCP_CPM_WT_mat)[2:9]),
  sprintf("%s.%s", colnames(expr_TCP_CPM_FL_mat)[1:8], colnames(expr_TCP_CPM_FL_mat)[2:9])
)
scaled_expr_TCP_CPM_WT_diff_mat_merge <- t(scale(t(expr_TCP_CPM_WT_diff_mat_merge)))

expr_TCP_CPM_WT_diff_mat_merge_bool <- expr_TCP_CPM_WT_diff_mat_merge>0
diff_dist <- dist(expr_TCP_CPM_WT_diff_mat_merge_bool)
diff_expr_dist <- dist(scaled_expr_TCP_CPM_WT_diff_mat_merge)
merge_dist <- diff_dist + diff_expr_dist / 2
clu <- hclust(merge_dist, "ward.D2")

interest_gene_df <- left_join(data.frame(Gid=ExprTCP_df$Gid), interest_gene)
interest_gene_df$Label <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$GeneSymbol)
interest_gene_df$Label[!is.na(interest_gene_df$Family)] <- sprintf("%s (%s: %s)", interest_gene_df$Gid, interest_gene_df$Family, interest_gene_df$GeneSymbol)[!is.na(interest_gene_df$Family)]
interest_gene_df$Label[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$Family)[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]
interest_gene_df$Label[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- interest_gene_df$Gid[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]

subindx <- which(!is.na(interest_gene_df$GeneSymbol))
label_li <- interest_gene_df$Label[subindx]

scaled_TP_ht <- Heatmap(
  scaled_expr_TCP_CPM_WT_diff_mat_merge, 
  name="z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  row_order = clu$order,
  cluster_columns = FALSE, 
  show_row_names = F, 
  split = sprintf("%s.%s", ExprTCP_df$CircadianState, ExprTCP_df$CircadianType),
  column_split = c(rep("WT", 8), rep("FL", 8)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

col_li <- c("TRUE"="black", "FALSE"="white")
ha <- rowAnnotation(
  IsTCP=ExprTCP_df$IsTcpMotif,
  IsTcpLike=ExprTCP_df$IsTcpLikeMotif,
  # link = row_anno_link(
  #   at = subindx,
  #   labels = label_li,
  #   labels_gp=gpar(fontsize = 6),
  #   lines_gp=gpar(lwd = 0.5),
  #   link_width = unit(4, "mm")),
  # width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)),
  width = unit(4, "mm"),
  col=list(
    IsTCP=col_li,
    IsTcpLike=col_li
  )
  )

inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.Circadian.pdf", width=10/inche_cm, height=8/inche_cm)
print(scaled_TP_ht+ha)
dev.off()
write_tsv(ExprTCP_df, "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.Circadian.tsv")

expr_CPM <- CPM_merge[apply(CPM_merge[,2:ncol(CPM_merge)], 1, max)>1,]
expr_CPM_mat <- as.matrix(expr_CPM[,2:ncol(expr_CPM)])
rownames(expr_CPM_mat) <- expr_CPM$Gid

expr_CPM_WT_mat <- expr_CPM_mat[,1:9]
expr_CPM_FL_mat <- expr_CPM_mat[,10:18]
expr_CPM_WT_diff_mat <- t(apply(expr_CPM_WT_mat, 1, function(x){return(diff(x)>0)}))
expr_CPM_FL_diff_mat <- t(apply(expr_CPM_FL_mat, 1, function(x){return(diff(x)>0)}))



Gbox <- data.frame(Gid=motif_hit_df$Gid[motif_hit_df$Gbox])
MYB_G4 <- data.frame(Gid=motif_hit_df$Gid[motif_hit_df$MYB_G4])
MYB_WT <- data.frame(Gid=motif_hit_df$Gid[motif_hit_df$MYB_WT])
AP2 <- data.frame(Gid=motif_hit_df$Gid[motif_hit_df$AP2])
TCP <- read_delim("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_like <- read_delim("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_like_motif/TCP_like_motif.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
CircadianGene_Group <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CircadianGene.Group.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

do_test <- function(gids1, gids2, total_gene_num, label1, label2){
  gid1_num <- length(unique(gids1))
  gid2_num <- length(unique(gids2))
  overlap_num <- length(unique(intersect(gids1, gids2)))
  exp_overlap_num <- total_gene_num * (gid1_num/total_gene_num) * (gid2_num/total_gene_num)
  gene_num_mat <- matrix(0, nrow = 2, ncol = 2)
  gene_num_mat[1,1] <- overlap_num
  gene_num_mat[2,1] <- gid1_num - overlap_num
  gene_num_mat[1,2] <- gid2_num - overlap_num
  gene_num_mat[2,2] <- total_gene_num - gid1_num - gid2_num + overlap_num
  pval <- fisher.test(gene_num_mat)$p.value
  res <- data.frame(Set1=label1, Set2=label2, Set1Num=gid1_num, Set2Num=gid2_num, OverlapNum=overlap_num, ExpOverlapNum=exp_overlap_num, pval=pval)
  return(res)
}
gene_li <- list(Gbox=Gbox$Gid, MYB_G4=MYB_G4$Gid, MYB_WT=MYB_WT$Gid, AP2=AP2$Gid, TCP=TCP$Gid, TCP_like=TCP_like$Gid, CircadianGeneG13=CircadianGene_Group$Gid[CircadianGene_Group$GroupID!=4], CircadianGeneG4=CircadianGene_Group$Gid[CircadianGene_Group$GroupID==4])
total_gene_num <- nrow(read_delim("analysis_v3/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE))
res <- data.frame()
for(label1 in names(gene_li)){
  for(label2 in names(gene_li)){
    res <- rbind(res, do_test(gene_li[[label1]], gene_li[[label2]], total_gene_num, label1, label2))
  }
}
write_tsv(res, "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_Gbox_Circadian.test.tsv")
res$Label <- sprintf("p=%.1E", res$pval)
TCP_Gbox_test <- res[res$Set1%in%c("TCP", "TCP_like", "Gbox", "MYB_G4", "MYB_WT", "AP2") & res$Set2%in%c("TCP", "TCP_like", "Gbox", "MYB_G4", "MYB_WT", "AP2"),]

TCP_Gbox_melt_df <- melt(TCP_Gbox_test[,c("Set1", "Set2", "OverlapNum", "ExpOverlapNum")], c("Set1", "Set2"), variable.name = "Stat", value.name = "GeneNum")
TCP_Gbox_melt_df$Stat <- factor(TCP_Gbox_melt_df$Stat, levels = c("OverlapNum", "ExpOverlapNum"), labels = c("Overlap gene", "Exp."))
p <- ggplot(TCP_Gbox_melt_df, aes(x=Stat, y=GeneNum, fill=Stat)) +
  geom_bar(stat="identity", position = position_dodge(0.7), width = 0.7) +
  geom_text(mapping = aes(label=GeneNum), vjust=0, size=1.3, color="black", position = position_dodge(0.7)) +
  geom_text(data = TCP_Gbox_test, mapping = aes(x=1.5, y=0, label=Label, fill=NA), vjust=0, size=1.3, color="black") +
  facet_wrap(Set1~Set2, scales = "free_y") +
  scale_fill_manual(values = c("Overlap gene"="black", "Exp."="grey70")) +
  labs(y="#Genes") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
p
ggsave("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_Gbox.test.pdf", p, width = 15, height = 15, limitsize = FALSE, units = "cm")

TCP_Gbox_Circadian_test <- res[res$Set1%in%c("TCP", "TCP_like", "Gbox", "MYB_G4", "MYB_WT", "AP2") & res$Set2%in%c("CircadianGeneG13", "CircadianGeneG4"),]
TCP_Gbox_Circadian_melt_df <- melt(TCP_Gbox_Circadian_test[,c("Set1", "Set2", "OverlapNum", "ExpOverlapNum")], c("Set1", "Set2"), variable.name = "Stat", value.name = "GeneNum")
TCP_Gbox_Circadian_melt_df$Stat <- factor(TCP_Gbox_Circadian_melt_df$Stat, levels = c("OverlapNum", "ExpOverlapNum"), labels = c("Overlap gene", "Exp."))
p <- ggplot(TCP_Gbox_Circadian_melt_df, aes(x=Stat, y=GeneNum, fill=Stat)) +
  geom_bar(stat="identity", position = position_dodge(0.7), width = 0.7) +
  geom_text(mapping = aes(label=GeneNum), vjust=0, size=1.3, color="black", position = position_dodge(0.7)) +
  geom_text(data = TCP_Gbox_Circadian_test, mapping = aes(x=1.5, y=0, label=Label, fill=NA), vjust=0, size=1.3, color="black") +
  facet_wrap(Set1~Set2, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Overlap gene"="black", "Exp."="grey70")) +
  labs(y="#Genes") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
p
ggsave("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_Gbox_Circadian.test.pdf", p, width = 14, height = 6, limitsize = FALSE, units = "cm")
