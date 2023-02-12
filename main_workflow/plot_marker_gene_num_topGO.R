library(readr)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(topGO)
library(RColorBrewer)

pval_cutoff <- 0.05
TF_gene_list <- read_delim("D:/cottonSingleCell/annotation/TF/TF.gene_list.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

Marker_gene <- read_delim("D:/cottonSingleCell/stat/marker_gene/cluster/Marker.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Marker_gene$IsTF <- Marker_gene$Gid %in% TF_gene_list$Gid
marker_info <- Marker_gene %>% group_by(Cluster, IsTF) %>% summarise(Num=n())
p <- ggplot(marker_info, aes(x=Cluster, y=Num, fill=IsTF)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=Num), size=2, vjust=-0.1) +
  labs(y="#Marker gene") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        panel.grid = element_blank())
p
ggsave("D:/cottonSingleCell/stat/marker_gene/cluster/marker_gene.num.pdf", p, width = 4.5, height = 4.5, units = "cm")

graph_score <- read_delim("D:/cottonSingleCell/stat/Trajectory/Graph_enrich_gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
graph_score <- graph_score[,c("ID", "morans_test_statistic", "morans_I", "p_value", "q_value")]
names(graph_score)[1] <- "Gid"
graph_score <- left_join(graph_score, Marker_gene)
graph_score$IsTF <- graph_score$Gid %in% TF_gene_list$Gid
graph_score$Cluster[is.na(graph_score$Cluster)] <- "Other"
col_li <- c(brewer.pal(8, "Dark2")[1:5], "grey70")
names(col_li) <- c("C1", "C2", "C3", "C4", "C5", "Other")
graph_score$IsTF <- factor(graph_score$IsTF, levels = c(F, T), labels = c("non-TF", "TF"))
p <- ggplot(graph_score, aes(x=Cluster, y=morans_I, fill=Cluster)) +
  geom_boxplot(outlier.colour = NA, size=0.3) +
  facet_grid(~IsTF) +
  lims(y=c(min(graph_score$morans_I), 1)) +
  labs(y="moran's I") +
  theme_bw() +
  scale_fill_manual(values = col_li) +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.grid = element_blank())
p
ggsave("D:/cottonSingleCell/stat/marker_gene/cluster/marker_gene.graph_score.pdf", p, width = 7, height = 4, units = "cm")
graph_score_pval_df <- data.frame(Cluster=c("C1", "C2", "C3", "C4", "C5"), nTF_pval=NA, TF_pval=NA)
nTF_other <- graph_score$morans_I[graph_score$Cluster=="Other" & graph_score$IsTF=="non-TF"]
TF_other <- graph_score$morans_I[graph_score$Cluster=="Other" & graph_score$IsTF=="TF"]
for(clu in c("C1", "C2", "C3", "C4", "C5")){
  nTF_clu <- graph_score$morans_I[graph_score$Cluster==clu & graph_score$IsTF=="non-TF"]
  TF_clu <- graph_score$morans_I[graph_score$Cluster==clu & graph_score$IsTF=="TF"]
  if(length(nTF_clu)>1){
    graph_score_pval_df$nTF_pval[graph_score_pval_df$Cluster==clu] <- wilcox.test(nTF_other, nTF_clu)$p.value
  }
  if(length(TF_clu)>1){
    graph_score_pval_df$TF_pval[graph_score_pval_df$Cluster==clu] <- wilcox.test(TF_other, TF_clu)$p.value
  }
}
write_tsv(graph_score_pval_df, "D:/cottonSingleCell/stat/marker_gene/cluster/marker_gene.graph_score.pval.tsv")

cluster_li <- c("C1", "C2", "C3", "C4", "C5")
source_li <- c("TF", "nTF")

top_df <- data.frame()
enrigh_term_df <- data.frame()
for (tmp_clu in cluster_li) {
  for (tmp_source in source_li) {
    tmp_file <- file.path("D:/cottonSingleCell/stat/marker_gene/cluster", sprintf("topGO_%s", tmp_source), tmp_clu, "topGO.tsv")
    if(!file.exists(tmp_file)) next()
    tmp_topGO <- read_delim(tmp_file, "\t", escape_double = FALSE, col_types = cols(classicFisher = col_double()), trim_ws = TRUE)
    tmp_topGO$classicFisher[is.na(tmp_topGO$classicFisher)] <- 1e-30
    tmp_topGO$Label <- sprintf("%s (%d)", tmp_topGO$Term, tmp_topGO$Significant)
    tmp_topGO$Term <- factor(tmp_topGO$Term, levels = rev(tmp_topGO$Term[!duplicated(tmp_topGO$Term)]))
    tmp_topGO$Cluster <- tmp_clu
    tmp_topGO$TF <- tmp_source
    tmp_topGO$classicKS <- NULL
    tmp_topGO$elimKS <- NULL
    enrigh_term_df <- rbind(enrigh_term_df, tmp_topGO[tmp_topGO$classicFisher<pval_cutoff,])
    tmp_topGO <- plyr::ddply(tmp_topGO, "Ontology", function(block){return(block[1:10,])})
    tmp_topGO$Indx <- 1:nrow(tmp_topGO)
    top_df <- rbind(top_df, tmp_topGO)
  }
}
write_tsv(enrigh_term_df, "D:/cottonSingleCell/stat/marker_gene/cluster/Marker.enrich.term.tsv")

top_df <- na.omit(top_df)
p <- ggplot(top_df, aes(x=-Indx, fill=Ontology, y=-log10(classicFisher))) +
  theme_bw() +
  facet_grid(Cluster~TF) +
  geom_bar(stat="identity", width = 0.6) +
  geom_text(aes(label=Label), size=2, hjust=0) +
  scale_fill_brewer(palette = "Dark2") +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 10, 20, 30), labels = c("1E0", "1E-10", "1E-20", "<1E-30"), limits = c(0, 60)) +
  labs(y="p-value") +
  theme(
    axis.text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
  )
ggsave("D:/cottonSingleCell/stat/marker_gene/cluster/Marker.enrich.term.top.pdf", p, width = 15*2, height = 8*5, units = "cm")

