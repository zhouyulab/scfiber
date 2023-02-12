library(readr)
library(VennDiagram)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(DelayedMatrixStats)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
setwd("~/mu01/project/CottonSingleCell")
f_TCP_gene <- "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.gene.tsv"
f_TCP_like_motif_gene <- "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_like_motif/TCP_like_motif.gene.tsv"
f_Marker_gene <- "analysis_v3/stat/marker_gene/Marker.gene.tsv"
f_scRNA_CPM <- "analysis_v3/stat/cluster_cnt/cluster.tpm.tsv"
f_TP_CPM <- "~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv"
f_interest <- "~/mu01/project/CottonSCE_TP/data/InterestGeneID.tsv"
f_TCP_like_topGO <- "~/mu01/project/CottonSingleCell/analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_like_motif/topGO/topGO.tsv"
f_TCP_topGO <- "~/mu01/project/CottonSingleCell/analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/topGO/topGO.tsv"


TCP_gene <- read_delim(f_TCP_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_like_motif_gene <- read_delim(f_TCP_like_motif_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
Marker_gene <- read_delim(f_Marker_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
cluster_tpm <- read_delim(f_scRNA_CPM, "\t", escape_double = FALSE, trim_ws = TRUE)
expr_tmp <- cluster_tpm[apply(cluster_tpm[,2:ncol(cluster_tpm)], 1, function(x){return(max(x, na.rm = TRUE))})>1,]
expr_tmp$IsC3Max <- apply(expr_tmp[2:ncol(expr_tmp)], 1, which.max)==3
table(expr_tmp$IsC3Max)
TP_CPM_merge <- read_delim(f_TP_CPM, "\t", escape_double = FALSE, trim_ws = TRUE)

C3_marker_gene <- Marker_gene[Marker_gene$Cluster=="C3",]
pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_TCP_like_motif.C3_marker.venn.pdf", width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(
  TCP_like_motif=TCP_like_motif_gene$Gid, 
  TCPmotif=TCP_gene$Gid,
  C3_marker_gene=C3_marker_gene$Gid
  ), filename = NULL,
  fill = brewer.pal(3, "Dark2"),
  borders = FALSE,
  sub.cex = 0.3,
  lwd=0.5)
grid.draw(p)
dev.off()

pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_TCP_like_motif.C3_max.venn.pdf", width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(
  TCP_like_motif=TCP_like_motif_gene$Gid, 
  TCPmotif=TCP_gene$Gid,
  C3_max_gene=expr_tmp$Gid[expr_tmp$IsC3Max]
), filename = NULL,
fill = brewer.pal(3, "Dark2"),
borders = FALSE,
sub.cex = 0.3,
lwd=0.5)
grid.draw(p)
dev.off()

TCP_like_motif_num <- length(unique(TCP_like_motif_gene$Gid))
TCP_motif_num <- length(unique(TCP_gene$Gid))
both_motif_num <- length(unique(intersect(TCP_like_motif_gene$Gid, TCP_gene$Gid)))
C3_max_gene_num <- sum(expr_tmp$IsC3Max)
C3_max_gene_df <- data.frame(
  Source=c("TCP motif", "TCP-like motif", "Both motif", "Other"),
  GeneNum=c(TCP_motif_num-both_motif_num, TCP_like_motif_num-both_motif_num, both_motif_num, C3_max_gene_num-TCP_like_motif_num-TCP_motif_num+both_motif_num)
)
C3_max_gene_df$Ratio <- C3_max_gene_df$GeneNum / C3_max_gene_num
C3_max_gene_df <- C3_max_gene_df[order(C3_max_gene_df$GeneNum, decreasing = T),]
C3_max_gene_df$Source <- factor(C3_max_gene_df$Source, levels = C3_max_gene_df$Source)
C3_max_gene_df$Label <- sprintf("%s (%s, %.1f%%)", C3_max_gene_df$Source, C3_max_gene_df$GeneNum, C3_max_gene_df$Ratio*100)
p <- ggplot(C3_max_gene_df, aes(x="", y=-Ratio, fill=Source)) +
  theme_bw() +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(y = sum(-Ratio)+Ratio/2-c(0, cumsum(-Ratio)[-length(Ratio)]), 
                label = Label), size=1.2, color="white", fontface = "bold") +
  scale_fill_brewer(palette = "Set1") +
  theme(
    text = element_text(size = 6),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(0, units = "mm")
  )
ggsave("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_TCP_like_motif.C3_max.pie.pdf", p, width = 4, height = 4, units = "cm")

# delta_phi_cutoff <- 0.2
# FDR_cutoff <- 0.05
# read_num_cutoff <- 20
# AS_type_li <- c("A3SS", "A5SS", "RI", "SE")
# AS_tp_li <- c("n48h", "n36h", "n24h", "n12h", "0h", "12h", "24h", "36h", "48h")
# AS_sample_li <- c("WT", "FL")
# AS_input <- "~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/rMATS_turbo"
# all_rmats_res_df <- data.frame()
# for(i in 1:(length(AS_tp_li))){
#   for(j in 1:i){
#     if(i==j) next()
#     for(s in AS_sample_li){
#       last_time <- AS_tp_li[j]
#       this_time <- AS_tp_li[i]
#       for(AS_type in AS_type_li){
#         fname <- file.path(AS_input, s, sprintf("%s_vs_%s", last_time, this_time), sprintf("%s.MATS.JC.txt", AS_type))
#         all_rmats_res_df <- rbind(all_rmats_res_df, load_data(fname, s, last_time, this_time, AS_type))
#       }
#     }
#   }
#   
# }
# 
# all_rmats_res_df$AllSample1IJC <- sapply(strsplit(all_rmats_res_df$IJC_SAMPLE_1, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
# all_rmats_res_df$AllSample1SJC <- sapply(strsplit(all_rmats_res_df$SJC_SAMPLE_1, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
# all_rmats_res_df$AllSample2IJC <- sapply(strsplit(all_rmats_res_df$IJC_SAMPLE_2, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
# all_rmats_res_df$AllSample2SJC <- sapply(strsplit(all_rmats_res_df$SJC_SAMPLE_2, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
# all_rmats_res_df$AllSample1 <- all_rmats_res_df$AllSample1IJC + all_rmats_res_df$AllSample1SJC
# all_rmats_res_df$AllSample2 <- all_rmats_res_df$AllSample2IJC + all_rmats_res_df$AllSample2SJC
# 
# all_rmats_res_df$AveIncLevel1 <- sapply(strsplit(all_rmats_res_df$IncLevel1, ","), function(x){return(mean(as.numeric(x), na.rm = T))})
# all_rmats_res_df$AveIncLevel2 <- sapply(strsplit(all_rmats_res_df$IncLevel2, ","), function(x){return(mean(as.numeric(x), na.rm = T))})
# 
# 
# all_rmats_res_df$DE <- "NC"
# all_rmats_res_df$DE[
#   all_rmats_res_df$IncLevelDifference>delta_phi_cutoff & all_rmats_res_df$FDR < FDR_cutoff &
#     all_rmats_res_df$AllSample1>read_num_cutoff & all_rmats_res_df$AllSample2>read_num_cutoff
#   ] <- "Up"
# all_rmats_res_df$DE[
#   (-1*all_rmats_res_df$IncLevelDifference)>delta_phi_cutoff & all_rmats_res_df$FDR < FDR_cutoff &
#     all_rmats_res_df$AllSample1>read_num_cutoff & all_rmats_res_df$AllSample2>read_num_cutoff
#   ] <- "Down"
# 
# write_tsv(all_rmats_res_df, "analysis_v3/ATAC/ATAC_kmer_cluster2gene/All.TP.AS.SourceData.tsv")
# DE_AS_df <- all_rmats_res_df[all_rmats_res_df$DE!="NC",]
# all_as_info <- DE_AS_df %>% group_by(Sample, LastTime, ThisTime) %>% summarise(DEGeneNum=length(unique(GeneID)))
# all_as_info$DEGeneRatio <- all_as_info$DEGeneNum / 70199
# write_tsv(all_as_info, "analysis_v3/ATAC/ATAC_kmer_cluster2gene/All.AS.Ratio.tsv")
# 
# TCP_AS_df <- DE_AS_df[DE_AS_df$GeneID%in%TCP_gene$Gid,]
# TCP_as_info <- TCP_AS_df %>% group_by(Sample, LastTime, ThisTime) %>% summarise(DEGeneNum=length(unique(GeneID)))
# TCP_as_info$DEGeneRatio <- TCP_as_info$DEGeneNum / nrow(TCP_gene)
# write_tsv(TCP_as_info, "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP.AS.Ratio.tsv")
# 
# TCP_like_motif_AS_df <- DE_AS_df[DE_AS_df$GeneID%in%TCP_like_motif_gene$Gid,]
# TCP_like_motif_as_info <- TCP_like_motif_AS_df %>% group_by(Sample, LastTime, ThisTime) %>% summarise(DEGeneNum=length(unique(GeneID)))
# TCP_like_motif_as_info$DEGeneRatio <- TCP_like_motif_as_info$DEGeneNum / nrow(TCP_like_motif_gene)
# write_tsv(TCP_like_motif_as_info, "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_like_motif.AS.Ratio.tsv")

TCP_df <- data.frame(Gid=unique(c(TCP_gene$Gid, TCP_like_motif_gene$Gid)))
TCP_scRNA <- left_join(TCP_df, cluster_tpm)
TCP_TP <- left_join(TCP_df, TP_CPM_merge)
any_expr <- apply(TCP_scRNA[,2:ncol(TCP_scRNA)], 1, function(x){return(max(x, na.rm = T)>1)}) & apply(TCP_TP[,2:ncol(TCP_TP)], 1, function(x){return(max(x, na.rm = T)>1)})
any_expr[is.na(any_expr)] <- FALSE
TCP_df$IsTCP <- TCP_df$Gid %in% TCP_gene$Gid
TCP_df$IsTCPlike <- TCP_df$Gid %in% TCP_like_motif_gene$Gid
TCP_df <- TCP_df[any_expr,]
TCP_scRNA <- TCP_scRNA[any_expr,]
TCP_TP <- TCP_TP[any_expr,]

scRNA_mat <- as.matrix(TCP_scRNA[,2:ncol(TCP_scRNA)])
rownames(scRNA_mat) <- TCP_scRNA$Gid
norm_scRNA_mat <- scRNA_mat / rowMaxs(scRNA_mat, na.rm = T)

TP_mat <- as.matrix(TCP_TP[,2:ncol(TCP_TP)])
rownames(TP_mat) <- TCP_TP$Gid
norm_TP_mat <- TP_mat / rowMaxs(TP_mat, na.rm = T)

scRNA_dist <- dist(norm_scRNA_mat)
TP_dist <- dist(norm_TP_mat)
scRNA_dist[is.na(scRNA_dist)] <- max(scRNA_dist, na.rm = T)
TP_dist[is.na(TP_dist)] <- max(TP_dist, na.rm = T)
merge_dist <- scRNA_dist + TP_dist
clu <- hclust(merge_dist, method = "ward.D2")

TCP_top_GO <- read_delim(f_TCP_topGO, "\t", escape_double = FALSE, trim_ws = TRUE)
TCPlike_top_GO <- read_delim(f_TCP_like_topGO, "\t", escape_double = FALSE, trim_ws = TRUE)
top_GO <- rbind(TCPlike_top_GO, TCP_top_GO)

ribosome_gid_li <- unique(unlist(strsplit(top_GO$Gids[top_GO$GO.ID=="GO:0005840"], " ")))
mitochondrion_gid_li <- unique(unlist(strsplit(top_GO$Gids[top_GO$GO.ID=="GO:0005739"], " ")))
RNA_binding_gid_li <- unique(unlist(strsplit(top_GO$Gids[top_GO$GO.ID=="GO:0003723"], " ")))
purine_ribonucleotide_metabolic_process_gid_li <- unique(unlist(strsplit(top_GO$Gids[top_GO$GO.ID=="GO:0009150"], " ")))

TCP_df$IsRibo <- TCP_df$Gid %in% ribosome_gid_li
TCP_df$IsMitochondrion <- TCP_df$Gid %in% mitochondrion_gid_li
TCP_df$IsPurine <- TCP_df$Gid %in% purine_ribonucleotide_metabolic_process_gid_li
TCP_df$IsRnaBinding <- TCP_df$Gid %in% RNA_binding_gid_li

interest_gene <- read_delim(f_interest, "\t", escape_double = FALSE, trim_ws = TRUE)
interest_gene <- interest_gene[!duplicated(interest_gene[,c("Gid")]),]

interest_df <- left_join(TCP_df, interest_gene)
interest_df$Label <- sprintf("%s (%s)", interest_df$Gid, interest_df$GeneSymbol)
subindx <- which(!is.na(interest_df$GeneSymbol))
label_li <- interest_df$Label[subindx]

ha_col <- c("black", "white")
names(ha_col) <- as.character(c(T, F))
ha <- rowAnnotation(
  IsTCP = TCP_df$IsTCP,
  IsTCPLike = TCP_df$IsTCPlike,
  isRibosome = TCP_df$IsRibo,
  isMitochondrion = TCP_df$IsMitochondrion,
  isPurineMatabolic = TCP_df$IsPurine,
  isRnaBinding = TCP_df$IsRnaBinding,
  link = row_anno_link(
    at = subindx,
    labels = label_li,
    labels_gp=gpar(fontsize = 6),
    lines_gp=gpar(lwd = 0.5),
    link_width = unit(4, "mm")),
  col=list(
    IsTCP=ha_col,
    IsTCPLike=ha_col,
    isRibosome=ha_col,
    isMitochondrion=ha_col,
    isRnaBinding=ha_col,
    isPurineMatabolic=ha_col
  ),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6))
)

scRNA_CPM_ht <- Heatmap(
  norm_scRNA_mat,
  name = "Norm. scRNA CPM",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  show_row_names=FALSE,
  cluster_columns=FALSE,
  row_order=clu$order,
  column_split=c(rep("WT", 5), rep("FL", 5)),
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

TP_CPM_ht <- Heatmap(
  norm_TP_mat,
  name = "Norm. TP. CPM",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  show_row_names=FALSE,
  cluster_columns=FALSE,
  row_order=clu$order,
  column_split=c(rep("WT", 9), rep("FL", 9)),
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
x <- scRNA_CPM_ht + TP_CPM_ht + ha
save(x, file = "TCP_TCPlike.expr.RData")

inche_cm=2.54
pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_TCPlike_motif.expr.ht.pdf", width=25/inche_cm, height=12/inche_cm, family="ArialMT", colormodel = "cmyk")
scRNA_CPM_ht + TP_CPM_ht + ha
dev.off()
