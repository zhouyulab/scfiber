library(RColorBrewer)
library(scran)
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(readr)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_mergeSCE <- args[1]
  f_corSCE <- args[2]
  f_cluster <- args[3]
  f_out <- args[4]
} else {
  q()
}

# f_cluster = "~/mu01/project/CottonSingleCell/analysis_v3/final_cluster/cluster.tsv"
# f_out = "~/mu01/project/CottonSingleCell/analysis_v3/final_cluster/all_umap"

load(f_corSCE)
cluster <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)

if(!all(cluster$Cell == colnames(corSCE))) stop()

corSCE$cluster <- cluster$Cluster


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

plot_gene_expr_umap <- function(umapSCE, mat, gid, label){
  expr <- mat[which(rownames(mat)==gid),]
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

plot_gene_expr_umap_all_sample <- function(umapSCE, mat, gid, label, sample){
  expr <- mat[which(rownames(mat)==gid),]
  if(! all(names(expr)==colnames(umapSCE))) stop()
  df <- data.frame(x=umapSCE@int_colData$reducedDims$UMAP[,1], y=umapSCE@int_colData$reducedDims$UMAP[,2], expr=expr, sample=sample)
  
  
  p <- ggplot(df, aes(x=x, y=y, alpha=expr, color=sample)) +
    theme_bw() +
    geom_point(size=0.1, alpha=0.01, color="grey70") +
    geom_point(size=0.1) +
    scale_color_manual(values = c("WT"="red", "FL"="blue")) +
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

samples <- sort(unique(corSCE$SampleName))

load(f_mergeSCE)
if(!all(cluster$Cell == colnames(mergeSCE))) stop()

all_gene_li <- rownames(mergeSCE)
cnt_mat <- counts(mergeSCE)
cnt_mat_li <- list()
for(s in samples){
    tmpSCE <- mergeSCE[,mergeSCE$SampleName==s]
    cnt_mat_li[[s]] <- counts(tmpSCE)
}

# all_gene_li <- c("Ghir_A11G028960", "Ghir_D11G029140", "Ghir_A02G004990", "Ghir_A11G020500", "Ghir_D11G020520", "Ghir_D10G010290", "Ghir_A10G008990")

for(gid in all_gene_li){
  print(gid)
  all_cnt <- sum(cnt_mat[which(rownames(mergeSCE)==gid),])
  if(all_cnt<100) next()
  prefix <- substr(gid, 1, 8)
  if(prefix==gid){
    prefix="Other"
  }
  out_dir <- file.path(f_out, prefix)
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }

  f_plot <- file.path(out_dir, sprintf("%s.pdf", gid))
  # if(file.exists(f_plot)) next()
  ps_gene_marker <- list()
  label <- sprintf("%s: %s", "All", gid)
  ps_gene_marker[[label]] <-  plot_gene_expr_umap_all_sample(corSCE, cnt_mat, gid, label, corSCE$SampleName)
  for(s in samples){
    tmpUMAP <- corSCE[,corSCE$SampleName==s]
    tmpSCE <- mergeSCE[,mergeSCE$SampleName==s]
    label <- sprintf("%s: %s", s, gid)
    ps_gene_marker[[label]] <-  plot_gene_expr_umap(tmpUMAP, cnt_mat_li[[s]], gid, label)
    rm(tmpUMAP)
    rm(tmpSCE)
    gc()
  }
  inche_cm=2.54
  pdf(f_plot, width=8*3/inche_cm, height=6*ceiling(length(ps_gene_marker)/3)/inche_cm, family="ArialMT", colormodel = "cmyk")
  grid.arrange(grobs=ps_gene_marker,ncol=3,padding = unit(0, "mm"))
  dev.off()
  gc()
}

