library(readr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(RColorBrewer)
library(SingleCellExperiment)


args = commandArgs(TRUE)
if (length(args) == 5) {
  f_cluster <- args[1]
  f_merge <- args[2]
  f_corSCE <- args[3]
  f_mergeSCE <- args[4]
  f_out <- args[5]
} else {
  q()
}

# f_cluster <- "~/data/CottonSingleCell/analysis/base_map/overcluster.tsv"
# f_merge <- "~/data/CottonSingleCell/analysis/base_map/merge.tsv"
# f_out <- "~/data/CottonSingleCell/analysis/final_cluster"

load(f_corSCE)
load(f_mergeSCE)
cluster <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)
merge_df <- read_delim(f_merge, "\t", escape_double = FALSE, trim_ws = TRUE)
cluster <- full_join(cluster, merge_df, by="Cluster")

if(!all(cluster$Barcode == colnames(corSCE))) stop()

corSCE$cluster <- cluster$Merge
corSCE$SubClu <- cluster$SubClu
mergeSCE$cluster <- cluster$Merge
mergeSCE$SubClu <- cluster$SubClu
corSCE <- corSCE[,!is.na(corSCE$cluster)]
mergeSCE <- mergeSCE[,!is.na(mergeSCE$cluster)]
save(corSCE, file=file.path(f_out, "umapSCE.RData"))
save(mergeSCE, file=file.path(f_out, "mergeSCE.RData"))

plot_umap  <- function(umapSCE, col_li, label){
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
    geom_point(size=0.1, alpha=0.6) +
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

plot_cluster_umap  <- function(umapSCE, cluster, color, label){
  df <- data.frame(
    x=umapSCE@int_colData$reducedDims$UMAP[,1], 
    y=umapSCE@int_colData$reducedDims$UMAP[,2], 
    cluster=umapSCE$cluster
  )
  df_cluster <- df[df$cluster==cluster,]

  p <- ggplot(data = NULL, aes(x=x, y=y)) +
    theme_bw() +
    geom_point(data=df, size=0.1, alpha=0.3, color="grey70") +
    geom_point(data=df_cluster, size=0.1, alpha=0.3, color=color) +
    labs(x="UMAP1", y="UMAP2", title=label) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
}


clu_li <- sort(unique(corSCE$cluster))
col_li <- brewer.pal(length(clu_li), "Dark2")
names(col_li) <- clu_li


sample_li <- sort(unique(corSCE$SampleName))
umapSCE_li <- list("All"=corSCE)
ps_li <- list()

label_text <- "All: All cluster"
ps_li[[label_text]] <- plot_umap(corSCE, col_li, label_text)
for(s in sample_li){
  tmp_umap <- corSCE[,corSCE$SampleName==s]
  umapSCE_li[[s]] <- tmp_umap
  label_text <- sprintf("%s: All cluster (#Cell=%d)", s, ncol(tmp_umap))
  ps_li[[label_text]] <- plot_umap(umapSCE_li[[s]], col_li, label_text)
}

for(clu in clu_li){
  color <- col_li[names(col_li)==clu]
  label_text <- sprintf("All: %s (#Cell=%d)", clu, sum(corSCE$cluster==clu))
  ps_li[[label_text]] <- plot_cluster_umap(corSCE, clu, color, label_text)
  
  for(s in sample_li){
    label_text <- sprintf("%s: %s (#Cell=%d)", s, clu, sum(umapSCE_li[[s]]$cluster==clu))
    ps_li[[label_text]] <- plot_cluster_umap(umapSCE_li[[s]], clu, color, label_text)
  }
}

inche_cm=2.54
pdf(file.path(f_out, "cluster.umap.pdf"), width=8*3/inche_cm, height=6+6.5*(ceiling(length(ps_li)/3)-1)/inche_cm, family="ArialMT", colormodel = "cmyk")
grid.arrange(grobs=ps_li,ncol=3,padding = unit(0, "mm"), heights=c(1, rep(1.33, ceiling(length(ps_li)/3)-1)))
dev.off()

df <- data.frame(Sample=corSCE$SampleName, Rep=corSCE$Rep, Cluster=corSCE$cluster, SubClu=corSCE$SubClu, Barcode=colnames(corSCE))
write_tsv(df, file.path(f_out, "cluster.tsv"))
