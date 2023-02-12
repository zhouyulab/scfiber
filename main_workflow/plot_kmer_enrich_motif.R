library(monocle3)
library(SingleCellExperiment)
library(ggplot2)
library(Matrix)
library(readr)
library(dplyr)
library(transport)
library(BiocParallel)
library(igraph)
library(ggtree)

setwd("~/project/CottonSCE")
load("~/project/CottonSCE/data/final_cluster/umapSCE.RData")
load("~/project/CottonSCE/analysis/kmer_graph_enrich/kmer.8.cell_RNA_mat.RData")
kmer_8_graph_enrich_score <- read_delim("analysis/kmer_graph_enrich/kmer.8.graph_enrich_score.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
gene_num_mI_cutoff <- quantile(kmer_8_graph_enrich_score$ExprGeneNum_morans_I, 0.95)
UMI_mI_cutoff <- quantile(kmer_8_graph_enrich_score$UMI_morans_I, 0.95)
kmer_8_graph_enrich_score$ExprGeneNumPass <- kmer_8_graph_enrich_score$ExprGeneNum_morans_I > gene_num_mI_cutoff
kmer_8_graph_enrich_score$UMIPass <- kmer_8_graph_enrich_score$UMI_morans_I > UMI_mI_cutoff
pass_df <- kmer_8_graph_enrich_score %>% group_by(ExprGeneNumPass, UMIPass) %>% summarise(Number=n())
pass_df$x <- 0.1
pass_df$y <- 0.1
pass_df$x[pass_df$ExprGeneNumPass] <- 0.8
pass_df$y[pass_df$UMIPass] <- 0.8
pass_df$Label <- sprintf("n=%d\n(%.1f%%)", pass_df$Number, pass_df$Number / sum(pass_df$Number) * 100)

p <- ggplot(kmer_8_graph_enrich_score, aes(x=ExprGeneNum_morans_I, y=UMI_morans_I)) +
  geom_bin2d(bins=50) +
  geom_hline(yintercept = UMI_mI_cutoff, size=0.3, color="black") +
  geom_vline(xintercept = gene_num_mI_cutoff, size=0.3, color="black") +
  geom_text(pass_df, mapping=aes(x=x, y=y, label=Label), size=1.7) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x="Morans'I for # expressed gene", y="Morans'I for # UMI", fill="# Motifs") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.line = element_line(color="black", size=0.3),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.title = element_text(family="ArialMT", color = "black", size = 6),
    legend.text = element_text(family="ArialMT", color = "black", size = 6)
  )

ggsave("analysis/kmer_graph_enrich/kmer.8.graph_enrich_score.pdf", p, width = 7, height = 4.5, units = "cm")
enriched_motif <- kmer_8_graph_enrich_score[kmer_8_graph_enrich_score$ExprGeneNumPass&kmer_8_graph_enrich_score$UMIPass,]
enriched_motif <- enriched_motif[order(enriched_motif$ExprGeneNum_morans_I, decreasing = T),]
write_tsv(enriched_motif, "analysis/kmer_graph_enrich/kmer.8.graph_enrich_score.enriched.tsv")
df <- data.frame(
  Cluster=corSCE$cluster,
  Sample=corSCE$SampleName,
  Rep=corSCE$Rep,
  UMAP1=corSCE@int_colData@listData$reducedDims$UMAP[,1],
  UMAP2=corSCE@int_colData@listData$reducedDims$UMAP[,2]
)

for (i in 1:nrow(enriched_motif)) {
  motif <- enriched_motif$Motif[i]
  rev_motif <- enriched_motif$RevMotif[i]
  UMI_morans_I <- enriched_motif$UMI_morans_I[i]
  ExprGeneNum_morans_I <- enriched_motif$ExprGeneNum_morans_I[i]
  tmp_df <- df
  tmp_df$UMI <- kmer_cell_li$cnt_mat[motif,]
  tmp_df$ExprGeneNum <- kmer_cell_li$bool_mat[motif,]
  rep2_df <- tmp_df[tmp_df$Rep=="rep2",]
  rep3_df <- tmp_df[tmp_df$Rep=="rep3",]
  p1 <- ggplot() +
    geom_point(data=rep2_df, mapping = aes(x=UMAP1, y=UMAP2, color=log10(ExprGeneNum+1)), alpha=0.2, size=0.1) +
    labs(x="UMAP1", y="UMAP2", color="#ExprGene", title=sprintf("Motif: %s, Reverse motif: %s V2", motif, rev_motif)) +
    scale_color_gradientn(
      breaks=log10(c(0, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)+1),
      labels=c("0", "1E0", "1E1", "1E2", "1E3", "1E4", "1E5", "1E6", "1E7", "1E8"),
      colors=c("blue", "white", "red")) +
    facet_grid(~Sample) +
    theme_bw() +
    scale_alpha_continuous(range = c(0, 0.5)) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  p2 <- ggplot() +
    geom_point(data=rep3_df, mapping = aes(x=UMAP1, y=UMAP2, color=log10(ExprGeneNum+1)), alpha=0.2, size=0.1) +
    labs(x="UMAP1", y="UMAP2", color="#ExprGene", title=sprintf("Motif: %s, Reverse motif: %s V3", motif, rev_motif)) +
    scale_color_gradientn(
      breaks=log10(c(0, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)+1),
      labels=c("0", "1E0", "1E1", "1E2", "1E3", "1E4", "1E5", "1E6", "1E7", "1E8"),
      colors=c("blue", "white", "red")) +
    facet_grid(~Sample) +
    theme_bw() +
    scale_alpha_continuous(range = c(0, 0.5)) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  inche_cm <- 2.54
  pdf(sprintf("analysis/kmer_graph_enrich/UMAP/kmer_8/kmer.8.%s.ExprGene.pdf", motif), width=16/inche_cm, height=6/inche_cm)
  print(p1)
  print(p2)
  dev.off()

  p1 <- ggplot() +
    geom_point(data=rep2_df, mapping = aes(x=UMAP1, y=UMAP2, color=log10(UMI+1)), alpha=0.2, size=0.1) +
    labs(x="UMAP1", y="UMAP2", color="UMI", title=sprintf("Motif: %s, Reverse motif: %s V2", motif, rev_motif)) +
    scale_color_gradientn(
      breaks=log10(c(0, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)+1),
      labels=c("0", "1E0", "1E1", "1E2", "1E3", "1E4", "1E5", "1E6", "1E7", "1E8"),
      colors=c("blue", "white", "red")) +
    facet_grid(~Sample) +
    theme_bw() +
    scale_alpha_continuous(range = c(0, 0.5)) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  p2 <- ggplot() +
    geom_point(data=rep3_df, mapping = aes(x=UMAP1, y=UMAP2, color=log10(UMI+1)), alpha=0.2, size=0.1) +
    labs(x="UMAP1", y="UMAP2", color="UMI", title=sprintf("Motif: %s, Reverse motif: %s V3", motif, rev_motif)) +
    scale_color_gradientn(
      breaks=log10(c(0, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)+1),
      labels=c("0", "1E0", "1E1", "1E2", "1E3", "1E4", "1E5", "1E6", "1E7", "1E8"),
      colors=c("blue", "white", "red")) +
    facet_grid(~Sample) +
    theme_bw() +
    scale_alpha_continuous(range = c(0, 0.5)) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  pdf(sprintf("analysis/kmer_graph_enrich/UMAP/kmer_8/kmer.8.%s.UMI.pdf", motif), width=16/inche_cm, height=6/inche_cm)
  print(p1)
  print(p2)
  dev.off()
}

umi_kmer_mat <- kmer_cell_li$cnt_mat[enriched_motif$Motif,]
bool_kmer_mat <- kmer_cell_li$bool_mat[enriched_motif$Motif,]
rm(kmer_cell_li)
gc()

compute_kmer_distance <- function(mat, umap_df, resolution=0.1){
  # umap_df$UMAP1_bin <- floor(umap_df$UMAP1 / resolution) * resolution
  # umap_df$UMAP2_bin <- floor(umap_df$UMAP2 / resolution) * resolution
  
  kmer_num <- nrow(mat)
  dist_mat <- matrix(0, nrow = kmer_num, ncol = kmer_num)
  rownames(dist_mat) <- rownames(mat)
  colnames(dist_mat) <- rownames(mat)
  indx_li <- list()
  for(i in 2:kmer_num){
    for(j in 1:i){
      indx_li[[sprintf("%d_%d", i, j)]] <- c(i, j)
    }
  }
  dist_li <- bplapply(indx_li, function(x, mat, umap_df){
    i <- x[1]
    j <- x[2]
    cnt_i <- mat[i,]
    cnt_j <- mat[j,]
    tmp_df <- data.frame(rep=umap_df$rep, sample=umap_df$sample, cluster=umap_df$cluster, UMAP1=umap_df$UMAP1, UMAP2=umap_df$UMAP2, sample1=cnt_i, sample2=cnt_j)
    tmp_clu_df <- tmp_df %>% 
      group_by(rep) %>%
      mutate(sample1_ratio=sample1/sum(sample1), sample2_ratio=sample2/sum(sample2)) %>% 
      group_by(sample, cluster) %>% 
      summarise(sample1_ratio=mean(sample1_ratio), sample2_ratio=mean(sample2_ratio))
    tmp_clu_df <- tmp_clu_df[!(tmp_clu_df$sample == "C3" & tmp_clu_df$cluster %in% c("FL.rep2", "FL.rep3")),]
    tmp_dist <- dist(t(as.matrix(tmp_clu_df[,c("sample1_ratio", "sample2_ratio")])))
    return(tmp_dist)
  }, BPPARAM=MulticoreParam(40), mat=mat, umap_df=umap_df)
  
  for(i in 2:kmer_num){
    for(j in 1:i){
      dist_mat[i, j] <- dist_li[[sprintf("%d_%d", i, j)]]
      dist_mat[j, i] <- dist_li[[sprintf("%d_%d", i, j)]]
    }
  }
  dist_mat <- as.dist(dist_mat)
  return(dist_mat)
}

bool_dist <- compute_kmer_distance(bool_kmer_mat, df)
save(bool_dist, file="analysis/kmer_graph_enrich/kmer.8.graph_enrich.bool.dist.RData")
umi_dist <- compute_kmer_distance(umi_kmer_mat, df)
save(umi_dist, file="analysis/kmer_graph_enrich/kmer.8.graph_enrich.umi.dist.RData")


# build_knn <- function(dist_mat, k=5){
#   dist_mat <- as.matrix(dist_mat)
#   n_sample <- nrow(dist_mat)
#   kmer <- rownames(dist_mat)
#   
#   from_indx <- c()
#   to_indx <- c()
#   for (tmp_from_indx in 1:n_sample) {
#     target_dist <- dist_mat[tmp_from_indx,]
#     target_indx <- (1:n_sample)[order(target_dist)][1:(k+1)]
#     target_indx <- target_indx[target_indx!=tmp_from_indx]
#     from_indx <- c(from_indx, rep(kmer[tmp_from_indx], k))
#     to_indx <- c(to_indx, kmer[target_indx])
#   }
#   edge_df <- data.frame(from=from_indx, to=to_indx)
#   g <- graph_from_data_frame(edge_df, directed=FALSE)
#   return(g)
# }
# 
# umi_clu_knn <- build_knn(umi_dist, k=5)
# umi_clu <- cluster_louvain(umi_clu_knn)
# 
# bool_clu_knn <- build_knn(bool_dist, k=5)
# bool_clu <- cluster_louvain(bool_clu_knn)


umi_cluster <- hclust(umi_dist, method = "ward.D2")
plot(umi_cluster)
rect.hclust(umi_cluster,k=2)  

kmer_clu_df <- data.frame(Motif=umi_cluster$labels, Cluster=cutree(umi_cluster, k=2))
kmer_clu_df$Cluster <- factor(kmer_clu_df$Cluster, levels = c(1, 2), labels = c("M1", "M2"))
kmer_clu_df$Motif <- as.character(kmer_clu_df$Motif)

groupInfo <- split(kmer_clu_df$Motif, kmer_clu_df$Cluster)
motif_ggtree <- ggtree(umi_cluster, aes(color=group))
motif_ggtree <- groupOTU(motif_ggtree, groupInfo)
motif_ggtree <- motif_ggtree + 
  scale_color_manual(values = c("M1"="blue", "M2"="red")) +
  theme(
    line = element_line(size=0.3)
  )
ggsave("analysis/kmer_graph_enrich/kmer.8.graph_enrich.umi.cluster.pdf", motif_ggtree, height = 6, width = 10, units = "cm", colormodel = "cmyk")



write_tsv(kmer_clu_df, "analysis/kmer_graph_enrich/kmer.8.graph_enrich.umi.cluster.tsv")

enriched_motif_clu <- left_join(enriched_motif, kmer_clu_df)
enriched_motif_clu$ExprGeneNumPass <- NULL
enriched_motif_clu$UMIPass <- NULL
write_tsv(enriched_motif_clu, "analysis/kmer_graph_enrich/kmer.8.graph_enrich_score.enriched.tsv")

df <- data.frame(
  barcode=corSCE$Barcode,
  sample=corSCE$Sample,
  rep=corSCE$Rep,
  cluster=corSCE$cluster,
  subClu=corSCE$SubClu,
  UMAP1=corSCE@int_colData@listData$reducedDims$UMAP[,1], 
  UMAP2=corSCE@int_colData@listData$reducedDims$UMAP[,2])

all(colnames(umi_kmer_mat) == df$barcode)

for (clu in unique(kmer_clu_df$Cluster)) {
  tmp_df <- df
  tmp_kmer_motif <- kmer_clu_df$Motif[kmer_clu_df$Cluster==clu]
  tmp_expr_mat <- umi_kmer_mat[tmp_kmer_motif,]
  tmp_norm_expr <- colMeans(log10(tmp_expr_mat))
  tmp_df$expr <- tmp_norm_expr
  rep2_df <- tmp_df[tmp_df$rep=="rep2",]
  rep2_df <- rep2_df[sample(1:nrow(rep2_df), nrow(rep2_df), replace = FALSE),]
  p1 <- ggplot() +
    geom_point(data=rep2_df, mapping = aes(x=UMAP1, y=UMAP2, color=expr), alpha=0.2, size=0.1) +
    labs(x="UMAP1", y="UMAP2", color="Average log10 UMI for target gene", title=sprintf("Motif cluster %s", clu)) +
    scale_color_gradientn(
      # breaks=log10(c(0, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)+1),
      # labels=c("0", "1E0", "1E1", "1E2", "1E3", "1E4", "1E5", "1E6", "1E7", "1E8"),
      colors=c("blue", "white", "red")) +
    facet_grid(~sample) +
    theme_bw() +
    scale_alpha_continuous(range = c(0, 0.5)) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  rep3_df <- tmp_df[tmp_df$rep=="rep3",]
  rep3_df <- rep3_df[sample(1:nrow(rep3_df), nrow(rep3_df), replace = FALSE),]
  p2 <- ggplot() +
    geom_point(data=rep3_df, mapping = aes(x=UMAP1, y=UMAP2, color=expr), alpha=0.2, size=0.1) +
    labs(x="UMAP1", y="UMAP2", color="Average log10 UMI for target gene", title=sprintf("Motif cluster %s", clu)) +
    scale_color_gradientn(
      # breaks=log10(c(0, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)+1),
      # labels=c("0", "1E0", "1E1", "1E2", "1E3", "1E4", "1E5", "1E6", "1E7", "1E8"),
      colors=c("blue", "white", "red")) +
    facet_grid(~sample) +
    theme_bw() +
    scale_alpha_continuous(range = c(0, 0.5)) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  inche_cm <- 2.54
  pdf(sprintf("analysis/kmer_graph_enrich/cluster_UMAP/kmer_8/kmer.8.cluster.%s.UMI.pdf", clu), width=16/inche_cm, height=6/inche_cm)
  print(p1)
  print(p2)
  dev.off()
}
