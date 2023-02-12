library(readr)
library(ggplot2)
library(gridExtra)
library(Seurat)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_corSCE <- args[1]
  f_noCorSCE <- args[2]
  f_marker_gene <- args[3]
  f_out <- args[4]
} else {
  q()
}

# f_corSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/v2_analysis/correction_params/Seurat_norm_HVG_correction/n_feature_2000/cca_dim_50/corSCE.RData"
# f_noCorSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/v2_analysis/correction_params/Seurat_norm_HVG_correction/n_feature_2000/cca_dim_50/noCorSCE.RData"
# f_marker_gene <- "/home/sqreb/ias/data/CottonSingleCell/data/protential_marker_gene.seurat.txt"

merker_gene <- as.character(unlist(read_delim(f_marker_gene, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)))

load(f_corSCE)
load(f_noCorSCE)

plot_gene_expr_umap <- function(SCE, gid, label){
  mat <-SCE@assays$RNA@counts
  if(!gid %in% rownames(SCE)){
    return(ggplot())
  }
  expr <- mat[which(rownames(SCE)==gid),]
  if(! all(names(expr)==colnames(SCE))) stop()
  
  
  df <- data.frame(x=SCE@reductions$umap@cell.embeddings[,1], y=SCE@reductions$umap@cell.embeddings[,2], expr=expr)
  
  p <- ggplot(df, aes(x=x, y=y, alpha=expr)) +
    theme_bw() +
    geom_point(size=0.1, alpha=0.01, color="grey70") +
    geom_point(size=0.1, color="red") +
    labs(x="UMAP1", y="UMAP2", alpha="UMI", title=sprintf("%s: %s", label, gid)) +
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

sample_li <- sort(unique(corSCE$Batch))
inche_cm=2.54
pdf(f_out, width=30/inche_cm, height=12/inche_cm, family="ArialMT", colormodel = "cmyk")

for(gid in merker_gene){
  ps <- list()
  
  ## No Correction
  for(s in sample_li){
    tmpSCE <- noCorSCE[,noCorSCE$Batch==s]
    label <- sprintf("No correction:%s", s)
    ps[[label]] <- plot_gene_expr_umap(tmpSCE, gid, label)
  }
  
  ## Correction
  for(s in sample_li){
    tmpSCE <- corSCE[,corSCE$Batch==s]
    label <- sprintf("Correction:%s", s)
    ps[[label]] <- plot_gene_expr_umap(tmpSCE, gid, label)
  }
  grid.arrange(grobs=ps,nrow=2,padding = unit(0, "mm"))
}
dev.off()
