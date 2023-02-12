library(RColorBrewer)
library(scran)
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(readr)

args = commandArgs(TRUE)
if (length(args) == 5) {
  f_SNN <- args[1]
  f_mergeSCE <- args[2]
  f_corSCE <- args[3]
  f_cluster <- args[4]
  f_out <- args[5]
} else {
  q()
}

# f_SNN <- "/home/sqreb/ias/data/CottonSingleCell/analysis/find_params/cluster_params/expr/SNN/correction_k_mnn_30/SNN_D_50_SNN_K_50_SNN_TYPE_number.RData"
# f_mergeSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/find_params/correction_params/MNN_norm_HVG_correction/k_mnn_30/sigma_mnn_0.3/mergeSCE.RData"
# f_corSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/find_params/correction_params/MNN_norm_HVG_correction/k_mnn_30/sigma_mnn_0.3/corSCE.RData"
# f_cluster <- "/home/sqreb/ias/data/CottonSingleCell/analysis/find_params/cluster_params/expr/louvain_cluster/correction_k_mnn_30/SNN_D_50_SNN_K_50_SNN_TYPE_number.cluster.tsv"
# f_out <- "~/test"

load(f_corSCE)
load(f_SNN)
cluster <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)

if(!all(cluster$Cell == colnames(corSCE))) stop()

corSCE$cluster <- cluster$Cluster


## Comutre Modularity
cluster_modularity_ratio <- clusterModularity(SNN_graph, cluster$Cluster, as.ratio=TRUE)
MASS::write.matrix(cluster_modularity_ratio, file.path(f_out, "Modularity.ratio.tsv"), sep="\t")
inche_cm=2.54
pdf(file.path(f_out, "Modularity.ratio.pdf"), width=8/inche_cm, height=7/inche_cm, family="ArialMT", colormodel = "cmyk")
p <- pheatmap(log10(cluster_modularity_ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
         col=rev(heat.colors(100)))
print(p)
dev.off()

cluster_modularity <- clusterModularity(SNN_graph, cluster$Cluster)
MASS::write.matrix(cluster_modularity, file.path(f_out, "Modularity.tsv"), sep="\t")
inche_cm=2.54
pdf(file.path(f_out, "Modularity.pdf"), width=8/inche_cm, height=7/inche_cm, family="ArialMT", colormodel = "cmyk")
p <- pheatmap(log10(cluster_modularity+1), cluster_cols=FALSE, cluster_rows=FALSE,
              col=rev(heat.colors(100)))
print(p)
dev.off()


## Cluster UMAP
uniq_clu <- sort(unique(cluster$Cluster))
clu_num <- length(uniq_clu)

if(clu_num <= 12){
  cb_palette <- brewer.pal(12, "Paired")
}else{
  cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                  "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                  "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                  "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                  "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                  "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
}


if(clu_num>length(cb_palette)) stop()

col_li <- cb_palette[1:clu_num]
names(col_li) <- uniq_clu

plot_cluster_umap  <- function(umapSCE, col_li, label){
  df <- data.frame(
    x=umapSCE@int_colData$reducedDims$UMAP[,1], 
    y=umapSCE@int_colData$reducedDims$UMAP[,2], 
    cluster=umapSCE$cluster
    )
  cell_num <- nrow(df)
  rand_indx <- sample(1:cell_num, cell_num)
  df <- df[rand_indx,]
  p <- ggplot(df, aes(x=x, y=y, color=cluster)) +
    theme_bw() +
    geom_point(size=0.1, alpha=0.2) +
    labs(x="UMAP1", y="UMAP2", alpha="UMI", title=label) +
    scale_color_manual(values = col_li) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
}

plot_gene_expr_umap <- function(umapSCE, SCE, gid, label){
  mat <- counts(SCE)
  expr <- mat[which(rownames(SCE)==gid),]
  if(! all(names(expr)==colnames(umapSCE))) stop()
  df <- data.frame(x=umapSCE@int_colData$reducedDims$UMAP[,1], y=umapSCE@int_colData$reducedDims$UMAP[,2], expr=expr)
  p <- ggplot(df, aes(x=x, y=y, alpha=expr)) +
    theme_bw() +
    geom_point(size=0.1, alpha=0.01, color="grey70") +
    geom_point(size=0.1, color="red") +
    labs(x="UMAP1", y="UMAP2", alpha="UMI", title=label) +
    scale_alpha_continuous(range = c(0, 0.5)) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  return(p)
}

ps_clu <- list()
p_clu_all <- plot_cluster_umap(corSCE, col_li, "All")
ps_clu[["All"]] <- p_clu_all

samples <- sort(unique(corSCE$SampleName))
for(s in samples){
  tmp_SCE <- corSCE[,corSCE$SampleName==s]
  ps_clu[[s]] <- plot_cluster_umap(tmp_SCE, col_li, s)
}

inche_cm=2.54
pdf(file.path(f_out, "Cluster.UMAP.pdf"), width=8*length(ps_clu)/inche_cm, height=6/inche_cm, family="ArialMT", colormodel = "cmyk")
grid.arrange(grobs=ps_clu,nrow=1,padding = unit(0, "mm"))
dev.off()

merker_gene <- c(
  "Ghir_A01G016200",
  "Ghir_A01G021600",
  "Ghir_A01G021620",
  "Ghir_A02G007720",
  "Ghir_A03G014630",
  "Ghir_A05G025590",
  "Ghir_A05G040320",
  "Ghir_D13G024880",
  "Ghir_A01G000540",
  "Ghir_D01G015820",
  "Ghir_D05G007410",
  "Ghir_D10G020630",
  "Ghir_A05G040540",
  "Ghir_A09G018540",
  "Ghir_D08G005520"
  )
load(f_mergeSCE)
if(!all(cluster$Cell == colnames(mergeSCE))) stop()
mergeSCE <- mergeSCE[rownames(mergeSCE) %in% merker_gene,]

ps_gene_marker <- ps_clu
for(gid in merker_gene){
  label <- sprintf("%s: %s", "All", gid)
  ps_gene_marker[[label]] <-  plot_gene_expr_umap(corSCE, mergeSCE, gid, label)
  for(s in samples){
    tmpUMAP <- corSCE[,corSCE$SampleName==s]
    tmpSCE <- mergeSCE[,mergeSCE$SampleName==s]
    label <- sprintf("%s: %s", s, gid)
    ps_gene_marker[[label]] <-  plot_gene_expr_umap(tmpUMAP, tmpSCE, gid, label)
    rm(tmpUMAP)
    rm(tmpSCE)
    gc()
  }
  gc()
}

inche_cm=2.54
pdf(file.path(f_out, "GeneMarker.UMAP.pdf"), width=8*3/inche_cm, height=6*ceiling(length(ps_gene_marker)/3)/inche_cm, family="ArialMT", colormodel = "cmyk")
grid.arrange(grobs=ps_gene_marker,ncol=3,padding = unit(0, "mm"))
dev.off()
