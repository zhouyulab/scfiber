library(readr)
library(ggplot2)
library(RColorBrewer)
library(scater)
library(scran)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(moments)
library(circlize)
library(VennDiagram)

pt <- read_delim("D:/cottonSingleCell/stat/Trajectory/pt.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
load("D:/cottonSingleCell/final_cluster/mergeSCE.RData")
load("D:/cottonSingleCell/final_cluster/umapSCE.RData")
all(corSCE$Barcode == pt$Barcode)
pt$Cluster <- corSCE$cluster
pt$PseudoTime[pt$Cluster %in% c("C4", "C5")] <- NA

df <- data.frame(
  Barcode=corSCE$Barcode,
  Sample=corSCE$SampleName,
  Tech=corSCE$Tech,
  Cluster=corSCE$cluster,
  UMAP1=corSCE@int_colData@listData$reducedDims$UMAP[,1],
  UMAP2=corSCE@int_colData@listData$reducedDims$UMAP[,2],
  PseudoTime=pt$PseudoTime
)

p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=PseudoTime)) +
  geom_point(size=0.1) +
  scale_color_continuous(low="blue", high="orange", na.value="grey70") +
  labs(x="UMAP1", y="UMAP2", fill="Pseudotime") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = c(0.2, 0.8),
    panel.grid = element_blank()
  )
p
ggsave("D:/cottonSingleCell/PseudoTime.pdf", p, height = 10, width = 12, units = "cm", colormodel = "cmyk")

Graph_enrich_gene <- read_delim("D:/cottonSingleCell/stat/Trajectory/Graph_enrich_gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Graph_enrich_gene <- na.omit(Graph_enrich_gene)

TP_FPKM <- read_delim("D:/cottonSingleCell/stat/Trajectory/TP.FPKM.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
cluster_tpm <- read_delim("D:/cottonSingleCell/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

fiber_markers <- findMarkers(mergeSCE, mergeSCE$cluster %in% c("C1", "C2", "C3"), pval.type="all", lfc=0, direction="up")
fiber_df <- fiber_markers@listData$`TRUE`
fiber_df <- as.data.frame(fiber_df)
fiber_df$Gid <- rownames(fiber_df)
fiber_df <- fiber_df[fiber_df$FDR < 1e-3,]
hist(log10(Graph_enrich_gene$morans_test_statistic), breaks=100)
hist(log10(Graph_enrich_gene$morans_I), breaks=100)
mean(Graph_enrich_gene$morans_I<1e-3)
hist(-log10(Graph_enrich_gene$q_value), breaks = 1000)
mean(Graph_enrich_gene$q_value<1e-9)
mean(Graph_enrich_gene$morans_test_statistic>10)



morans_test_statistic_cutoff <- quantile(Graph_enrich_gene$morans_test_statistic, 0.95)
enriched_Graph_enrich_gene <- Graph_enrich_gene[Graph_enrich_gene$morans_test_statistic>morans_test_statistic_cutoff,]
enriched_Graph_enrich_gene <- enriched_Graph_enrich_gene[
  order(enriched_Graph_enrich_gene$morans_test_statistic, decreasing = T), 
  c("ID", "morans_test_statistic", "morans_I", "p_value", "q_value")
  ]
names(enriched_Graph_enrich_gene)[1] <- "Gid"
enriched_Graph_enrich_gene <- enriched_Graph_enrich_gene[enriched_Graph_enrich_gene$Gid %in% rownames(fiber_df),]


compute_bin_expr <- function(start_time, step_time){
  cell_indx <- which((df$PseudoTime >= start_time) & (df$PseudoTime < (start_time + step_time)))
  tmp_norm_sc_mat <- mergeSCE@assays@data@listData$normcounts[rownames(mergeSCE@assays@data@listData$normcounts) %in% enriched_Graph_enrich_gene$Gid,cell_indx]
  tmp_mean_expr <- rowMeans(tmp_norm_sc_mat)
  res <- data.frame(Gid=names(tmp_mean_expr), Start=start_time, scExpr=tmp_mean_expr)
  return(res)
}
step_time=0.5
bin_expr <- lapply(seq(0, max(df$PseudoTime, na.rm = T), step_time), compute_bin_expr, step_time=step_time)
bin_expr <- do.call(rbind, bin_expr)
bin_expr_df <- dcast(bin_expr, Gid~Start)
bin_expr_mat <- as.matrix(bin_expr_df[,2:ncol(bin_expr_df)])
rownames(bin_expr_mat) <- bin_expr_df$Gid
max_value <- rowMaxs(bin_expr_mat)
mid_value <- rowMedians(bin_expr_mat)
hist(log10(max_value/mid_value), breaks=100)
10^0.5
mean((max_value/mid_value)>3)
max_mid_cutoff <- 3
bin_expr_mat <- bin_expr_mat[(max_value/mid_value)>3,]
norm_bin_expr_mat <- (bin_expr_mat) / rowMaxs(bin_expr_mat)


wt_ave <- norm_bin_expr_mat %*% 1:ncol(norm_bin_expr_mat) / ncol(norm_bin_expr_mat)
norm_bin_expr_mat <- norm_bin_expr_mat[order(unlist(wt_ave)),]

bin_2d <- rowSds(norm_bin_expr_mat)
hist(bin_2d, breaks = 50)
bin_3d <- skewness(t(norm_bin_expr_mat))
hist(bin_3d, breaks = 50)
bin_4d <- kurtosis(t(norm_bin_expr_mat))
hist(bin_4d, breaks = 100)

df_stat <- data.frame(D2=bin_2d/max(bin_2d), D3=bin_3d/max(bin_3d), D4=bin_4d/max(bin_4d))
stat_mat <- as.matrix(df_stat)
rownames(stat_mat)

ht <- Heatmap(
  norm_bin_expr_mat,
  # cluster_rows  = F,
  cluster_columns  = F,
  show_row_names  = F
  )
ht

tp_df <- data.frame(Gid=rownames(norm_bin_expr_mat))
tp_df <- left_join(tp_df, TP_FPKM, by="Gid")
tp_mat <- as.matrix(tp_df[,2:ncol(tp_df)])
norm_tp_mat <- (tp_mat) / rowMaxs(tp_mat)

ht2 <- Heatmap(
  norm_tp_mat,
  cluster_columns  = F,
  show_row_names  = F
)

ht3 <- Heatmap(
  stat_mat,
  cluster_columns  = F,
  show_row_names  = F
)


ht + ht2 + ht3

filtered_bin_expr_mat <- bin_expr_mat[(bin_2d>(0.5*max(bin_2d))) & (bin_4d<(0.25*max(bin_4d))),]
filtered_norm_bin_expr_mat <- norm_bin_expr_mat[(bin_2d>(0.5*max(bin_2d))) & (bin_4d<(0.25*max(bin_4d))),]
ht <- Heatmap(
  filtered_norm_bin_expr_mat,
  # col =  colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  cluster_columns  = F,
  show_row_names  = F
)

filtered_tp_df <- data.frame(Gid=rownames(filtered_norm_bin_expr_mat))
filtered_tp_df <- left_join(filtered_tp_df, TP_FPKM, by="Gid")
filtered_tp_mat <- as.matrix(filtered_tp_df[,2:ncol(filtered_tp_df)])
filtered_tp_mat[rowMaxs(filtered_tp_mat)<1,] <- NA
filtered_norm_tp_mat <- (filtered_tp_mat) / rowMaxs(filtered_tp_mat)

ht2 <- Heatmap(
  filtered_norm_tp_mat,
  # col =  colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  cluster_columns  = F,
  show_row_names  = F
)
ht + ht2
inche_cm=2.54
pdf("D:/cottonSingleCell/graph_enrich_gene.pdf", width=10/inche_cm, height=4/inche_cm, family="ArialMT")
draw(ht + ht2)
dev.off()



res_df <- data.frame(Gid=rownames(filtered_bin_expr_mat))
all(rownames(filtered_bin_expr_mat) == rownames(filtered_tp_mat))
res_df <- cbind(res_df, as.data.frame(filtered_bin_expr_mat))
res_df <- cbind(res_df, as.data.frame(filtered_tp_mat))
res_df$MaxTP <- colnames(filtered_bin_expr_mat)[apply(filtered_bin_expr_mat, 1, which.max)]
write_tsv(res_df, "D:/cottonSingleCell/graph_enrich_gene.filter.tsv")

TF_gene_list <- read_delim("D:/cottonSingleCell/TF.gene_list.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
tf_res_df <- inner_join(res_df, TF_gene_list, by="Gid")

df_plot <- df
df_plot$Expr <- mergeSCE@assays@data$normcounts[tf_res_df$Gid[1],]
p <- ggplot(df_plot, aes(x=UMAP1, y=UMAP2, alpha=Expr)) +
  theme_bw() +
  geom_point(size=0.1, alpha=0.01, color="grey70") +
  geom_point(size=0.1, color="red") +
  labs(x="UMAP1", y="UMAP2", alpha="UMI", title=tf_res_df$Gid[1]) +
  scale_alpha_continuous(range = c(0, 0.5)) +
  theme(
    text = element_text(family="ArialMT", size = 6),
    title = element_text(family="ArialMT", size = 6),
    legend.title = element_text(family="ArialMT", size = 6),
    legend.text = element_text(family="ArialMT", size = 6),
    panel.grid = element_blank(),
    legend.background = element_blank()
  )
p

df_plot$Expr2 <- mergeSCE@assays@data$normcounts[tf_res_df$Gid[2],]
p <- ggplot(df_plot, aes(x=UMAP1, y=UMAP2, alpha=Expr2)) +
  theme_bw() +
  geom_point(size=0.1, alpha=0.01, color="grey70") +
  geom_point(size=0.1, color="red") +
  labs(x="UMAP1", y="UMAP2", alpha="UMI", title=tf_res_df$Gid[2]) +
  scale_alpha_continuous(range = c(0, 0.5)) +
  theme(
    text = element_text(family="ArialMT", size = 6),
    title = element_text(family="ArialMT", size = 6),
    legend.title = element_text(family="ArialMT", size = 6),
    legend.text = element_text(family="ArialMT", size = 6),
    panel.grid = element_blank(),
    legend.background = element_blank()
  )
p

ModularityEnrichGene_group <- read_delim("D:/cottonSingleCell/modularity_gene_cluster/ModularityEnrichGene.group.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
df_mod <- left_join(res_df, ModularityEnrichGene_group)

Modularity_filter <- read_delim("D:/cottonSingleCell/modularity_gene_cluster/modularity/Modularity.filter.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
df_mod <- left_join(res_df, Modularity_filter)
df_mod <- na.omit(df_mod)
df_mod <- df_mod[order(as.numeric(df_mod$MaxExprPseudoTime)),]
