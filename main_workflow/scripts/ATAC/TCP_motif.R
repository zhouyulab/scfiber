library(readr)
library(dplyr)
library(ggplot2)
library(igraph)
library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(VennDiagram)
library(RColorBrewer)

f_ATAC_kmer <- "~/mu01/project/CottonSingleCell/analysis_v3/ATAC/ATAC_kmer_cluster2gene/kmer_clu.filter.M2.C3Max.SourceData.tsv"
f_TP_CPM <- "~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv"
f_interest <- "~/mu01/project/CottonSCE_TP/data/InterestGeneID.tsv"
f_topGO <- "~/mu01/project/CottonSingleCell/analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/topGO/topGO.tsv"
f_marker_gene <- "~/mu01/project/CottonSingleCell/analysis_v3/stat/marker_gene/Marker.gene.tsv"
f_time_course_gene <- "~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/TimeCourse/ImpulseDE2_enrich_gene.CPM.tsv"


kmer_clu_filter_M2_C3Max <- read_delim(f_ATAC_kmer, "\t", escape_double = FALSE, trim_ws = TRUE)

TP_CPM_merge <- read_delim(f_TP_CPM, "\t", escape_double = FALSE, trim_ws = TRUE)

TPC_indx <- unique(c(grep("TGGGC[CT]", kmer_clu_filter_M2_C3Max$Motif), grep("TGGGC[CT]", kmer_clu_filter_M2_C3Max$RevMotif)))
TCP_df <- kmer_clu_filter_M2_C3Max[TPC_indx,]
TCP_df$Is_TCP_Fwd_motif <- FALSE
TCP_df$Is_TCP_Fwd_motif[grep("TGGGC[CT]", TCP_df$Motif)] <- TRUE
TCP_df$IsTCPFwdDir <- xor(TCP_df$CisMotif == "Motif", TCP_df$Is_TCP_Fwd_motif==FALSE)
TCP_df$TCP <- TCP_df$RevMotif
TCP_df$TCP[TCP_df$Is_TCP_Fwd_motif] <- TCP_df$Motif[TCP_df$Is_TCP_Fwd_motif]
table(TCP_df$Is_TCP_Fwd_motif)
TCP_gene_info <- TCP_df %>% group_by(Gid) %>% summarise(TCP=paste(sort(unique(TCP)), collapse = " "), IsTCPFwdDir=paste(sort(unique(IsTCPFwdDir)), collapse = " "))
TCP_gene_info$IsTCPFwdDir <- factor(TCP_gene_info$IsTCPFwdDir, levels = c("TRUE", "FALSE", "FALSE TRUE"), labels = c("Forward", "Reverse", "Both"))
table(TCP_gene_info$IsTCPFwdDir)
TCP_dir_info <- TCP_gene_info %>% group_by(IsTCPFwdDir) %>% summarise(GeneNum=n())
TCP_dir_info$label <- sprintf("%s\nn=%s", TCP_dir_info$IsTCPFwdDir, TCP_dir_info$GeneNum)
TCP_dir_info$NumLabel <- sprintf("n=%s", TCP_dir_info$GeneNum)
p <- ggplot(TCP_dir_info, aes(x="", y=-GeneNum, fill=IsTCPFwdDir)) +
  theme_bw() +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(y = sum(-GeneNum)+GeneNum/2-c(0, cumsum(-GeneNum)[-length(GeneNum)]), 
                label = label), size=1.2, color="white", fontface = "bold") +
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
ggsave("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.motif.strand.pdf", p, width = 4, height = 4, units = "cm")

scRNA_CPM_df <- kmer_clu_filter_M2_C3Max[!duplicated(kmer_clu_filter_M2_C3Max$Gid), c(1:4, 16:25)]
TCP_gene_info <- left_join(TCP_gene_info, scRNA_CPM_df)
TCP_gene_info <- left_join(TCP_gene_info, TP_CPM_merge)
TCP_gene_info$WT_C3_vs_OtherMax <- TCP_gene_info$WT.C3 / apply(TCP_gene_info[,c("WT.C1", "WT.C2", "WT.C4", "WT.C5", "FL.C1", "FL.C2", "FL.C4", "FL.C5")], 1, max)

p <- ggplot(TCP_gene_info, aes(x=IsTCPFwdDir, y=log2(WT_C3_vs_OtherMax), fill=IsTCPFwdDir)) +
  geom_boxplot(size=0.2, outlier.colour = NA) +
  geom_text(data = TCP_dir_info, aes(y=3, label=NumLabel), size=1.2, color="black") +
  scale_fill_brewer(palette = "Set1") +
  labs(x="TCP motif direction", y="log2(WT.C3 / OtherMax)") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.motif.strand.expr.pdf", p, width = 4, height = 4, units = "cm")


TCP_gene_info$TcpMotifCluNum <- sapply(TCP_gene_info$Gid, function(gid){
  tmp_TCP_df <- TCP_df[TCP_df$Gid==gid,]
  start <- apply(tmp_TCP_df[,c("MotifStart", "MotifEnd")], 1, min)
  end <- apply(tmp_TCP_df[,c("MotifStart", "MotifEnd")], 1, max)
  if(nrow(tmp_TCP_df)==1){
    clu_num <- 1
  }else{
    start_li <- c()
    end_li <- c()
    for(i in 1:nrow(tmp_TCP_df)){
      for(j in 1:i){
        if(start[i]>=end[j] | start[j]>=end[i]){
          next()
        }else{
          start_li <- c(start_li, i)
          end_li <- c(end_li, j)
        }
      }
    }
    edge_df <- data.frame(from=start_li, to=end_li)
    graph <- graph_from_data_frame(edge_df, directed=FALSE)
    clu <- igraph::components(graph)
    clu_num <- length(igraph::groups(clu))
  }
  return(clu_num)
})
TCP_gene_info$TcpMotifCluNum <- factor(TCP_gene_info$TcpMotifCluNum, levels = sort(unique(TCP_gene_info$TcpMotifCluNum)), labels = as.character(sort(unique(TCP_gene_info$TcpMotifCluNum))))
TCP_motif_num_info <- TCP_gene_info %>% 
  group_by(TcpMotifCluNum) %>% 
  summarise(
    GeneNum=n(), 
    Pval=wilcox.test(log2(WT_C3_vs_OtherMax), log2(TCP_gene_info$WT_C3_vs_OtherMax[TCP_gene_info$TcpMotifCluNum=="1"]))$p.value
    )
TCP_motif_num_info$NumLabel <- sprintf("n=%s\nP=%.2E", TCP_motif_num_info$GeneNum, TCP_motif_num_info$Pval)
p <- ggplot(TCP_gene_info, aes(x=TcpMotifCluNum, y=log2(WT_C3_vs_OtherMax), fill=TcpMotifCluNum)) +
  geom_boxplot(size=0.2, outlier.colour = NA) +
  geom_text(data = TCP_motif_num_info, aes(y=3, label=NumLabel), size=1.2, color="black") +
  scale_fill_brewer(palette = "Set1") +
  labs(x="#TCP motif cluster", y="log2(WT.C3 / OtherMax)") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.motif.num.expr.pdf", p, width = 4, height = 4, units = "cm")

write_tsv(TCP_gene_info, "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.gene.tsv")

scRNA_CPM_df <- kmer_clu_filter_M2_C3Max[!duplicated(kmer_clu_filter_M2_C3Max$Gid), c(1:4, 16:25)]
TCP_gene_info <- left_join(TCP_gene_info, scRNA_CPM_df)
TCP_gene_info <- left_join(TCP_gene_info, TP_CPM_merge)
TCP_scRNA_CPM_df <- TCP_gene_info[,c("Gid", "WT.C1", "WT.C2", "WT.C3", "WT.C4", "WT.C5", "FL.C1", "FL.C2", "FL.C3", "FL.C4", "FL.C5")]
TCP_scRNA_CPM_mat <- as.matrix(TCP_scRNA_CPM_df[,2:ncol(TCP_scRNA_CPM_df)])
norm_TCP_scRNA_CPM_mat <- TCP_scRNA_CPM_mat / rowMaxs(TCP_scRNA_CPM_mat, na.rm = T)
TCP_TP_CPM_df <- TCP_gene_info[,names(TP_CPM_merge)]
TCP_TP_CPM_mat <- as.matrix(TCP_TP_CPM_df[,2:ncol(TCP_TP_CPM_df)])
norm_TCP_TP_CPM_mat <- TCP_TP_CPM_mat / rowMaxs(TCP_TP_CPM_mat, na.rm = T)

scRNA_dist <- dist(norm_TCP_scRNA_CPM_mat)
TP_dist <- dist(norm_TCP_TP_CPM_mat)
scRNA_dist[is.na(scRNA_dist)] <- max(scRNA_dist, na.rm = T)
TP_dist[is.na(TP_dist)] <- max(TP_dist, na.rm = T)
merge_dist <- scRNA_dist + TP_dist
clu <- hclust(merge_dist, method = "ward.D2")

interest_gene <- read_delim(f_interest, "\t", escape_double = FALSE, trim_ws = TRUE)
interest_gene <- interest_gene[!duplicated(interest_gene[,c("ID", "Gid")]),]

interest_TCP_df <- left_join(TCP_TP_CPM_df, interest_gene)
interest_TCP_df$Label <- sprintf("%s (%s)", interest_TCP_df$Gid, interest_TCP_df$GeneSymbol)
subindx <- which(!is.na(interest_TCP_df$GeneSymbol))
label_li <- interest_TCP_df$Label[subindx]

top_GO <- read_delim(f_topGO, "\t", escape_double = FALSE, trim_ws = TRUE)
ribosome_gid_li <- strsplit(top_GO$Gids[top_GO$GO.ID=="GO:0005840"], " ")[[1]]
mitochondrion_gid_li <- strsplit(top_GO$Gids[top_GO$GO.ID=="GO:0005739"], " ")[[1]]
RNA_binding_gid_li <- strsplit(top_GO$Gids[top_GO$GO.ID=="GO:0003723"], " ")[[1]]
purine_ribonucleotide_metabolic_process_gid_li <- strsplit(top_GO$Gids[top_GO$GO.ID=="GO:0009150"], " ")[[1]]
ha_col <- c("black", "white")
names(ha_col) <- as.character(c(T, F))
ha <- rowAnnotation(
  isRibosome = TCP_TP_CPM_df$Gid %in% ribosome_gid_li,
  isMitochondrion = TCP_TP_CPM_df$Gid %in% mitochondrion_gid_li,
  isPurineMatabolic = TCP_TP_CPM_df$Gid %in% purine_ribonucleotide_metabolic_process_gid_li,
  isRnaBinding = TCP_TP_CPM_df$Gid %in% RNA_binding_gid_li,
  link = row_anno_link(
    at = subindx,
    labels = label_li,
    labels_gp=gpar(fontsize = 6),
    lines_gp=gpar(lwd = 0.5),
    link_width = unit(4, "mm")),
  col=list(
    isRibosome=ha_col,
    isMitochondrion=ha_col,
    isRnaBinding=ha_col,
    isPurineMatabolic=ha_col
  ),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6))
  )

TCP_scRNA_CPM_ht <- Heatmap(
  norm_TCP_scRNA_CPM_mat,
  name = "Norm. scRNA CPM",
  col = colorRamp2(c(0, 1), c("white", "red")),
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

TCP_TP_CPM_ht <- Heatmap(
  norm_TCP_TP_CPM_mat,
  name = "Norm. TP. CPM",
  col = colorRamp2(c(0, 1), c("white", "red")),
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

inche_cm=2.54
pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.expr.ht.pdf", width=25/inche_cm, height=12/inche_cm, family="ArialMT", colormodel = "cmyk")
TCP_scRNA_CPM_ht + TCP_TP_CPM_ht + ha
dev.off()

marker_gene <- read_delim(f_marker_gene, "\t", escape_double = FALSE, trim_ws = TRUE) 
C3_marker_gene <- marker_gene[marker_gene$Cluster=="C3",]
mean(C3_marker_gene$Gid %in% TCP_gene_info$Gid)


time_cource_gene <- read_delim(f_time_course_gene, "\t", escape_double = FALSE, trim_ws = TRUE) 
late_expr_gene <- time_cource_gene[time_cource_gene$Group=="Late",]
mean(late_expr_gene$Gid %in% TCP_gene_info$Gid)
sum(late_expr_gene$Gid %in% TCP_gene_info$Gid)

pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.C3_marker.venn.pdf", width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(TCP_motif_gene=TCP_gene_info$Gid, C3_marker_gene=C3_marker_gene$Gid), filename = NULL,
     fill = brewer.pal(3, "Dark2")[1:2],
     borders = FALSE,
     sub.cex = 0.3,
     lwd=0.5)
grid.draw(p)
dev.off()

pdf("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.late_expr.venn.pdf", width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(TCP_motif_gene=TCP_gene_info$Gid, late_expr_gene=late_expr_gene$Gid), filename = NULL,
                  fill = brewer.pal(3, "Dark2")[c(1, 3)],
                  borders = FALSE,
                  sub.cex = 0.3,
                  lwd=0.5)
grid.draw(p)
dev.off()


interest_GO_df <- top_GO[top_GO$GO.ID %in% c("GO:0005840", "GO:0005739", "GO:0003723", "GO:0009150"),]
interest_GO_df <- interest_GO_df[order(interest_GO_df$Significant, decreasing = T),]
interest_GO_df$x <- sprintf("%s: %s", interest_GO_df$GO.ID, interest_GO_df$Term)
interest_GO_df$x <- factor(interest_GO_df$x, levels = rev(interest_GO_df$x))
interest_GO_df$label <- sprintf("%s / %s\nP=%s", interest_GO_df$Significant, nrow(TCP_gene_info), interest_GO_df$classicFisher)
p <- ggplot(interest_GO_df, aes(x=x, y=Significant, fill=x)) +
  geom_bar(stat = "identity") +
  geom_text(mapping=aes(label=label), size=1.2, color="black", hjust=1) +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  labs(y="#Target gene") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave("analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.term.gene_num.pdf", p, width = 10, height = 3.5, units = "cm")

