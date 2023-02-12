library(readr)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(reshape2)
library(BiocParallel)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_umapSCE <- args[1]
  f_mergeSCE <- args[2]
  chrom <- args[3]
  f_out <- args[4]
} else {
  q()
}

load(f_umapSCE)
load(f_mergeSCE)

sample_li <- c(rep("WT", 5), rep("FL", 5))
cluster_li <- rep(c("C1", "C2", "C3", "C4", "C5"), 2)

umap_df <- data.frame(
  x=corSCE@int_colData$reducedDims$UMAP[,1], 
  y=corSCE@int_colData$reducedDims$UMAP[,2], 
  cluster=corSCE$cluster,
  Sample=corSCE$SampleName,
  rep=corSCE$Tech
)
cnt_mat <- normcounts(mergeSCE)
rm(mergeSCE)
gc()


plot_umap <- function(gid, umap_df, cnt_mat, f_out){
  all_umap <- umap_df
  all_umap$Expr=cnt_mat[rownames(cnt_mat)==gid,]
  if(sum(all_umap$Expr)<10){
    return(NULL)
  }
  all_umap$Sample <- factor(all_umap$Sample, levels = c("WT", "FL"))
  
  zero_umap <- all_umap[all_umap$Expr==0,]
  expr_umap <- all_umap[all_umap$Expr>0,]
  cell_num <- nrow(expr_umap)
  rand_indx <- sample(1:cell_num, cell_num)
  expr_umap <- expr_umap[rand_indx,]
  p <- ggplot(mapping=aes(x=x, y=y, color=Sample, alpha=Expr)) 
    
  if(any(all_umap$Expr==0)){
    p <- p + geom_point(data=zero_umap, color="grey70", size=0.1, alpha=0.2)
  }

  p <- p +
    geom_point(data=expr_umap, size=0.1) 
    labs(x="UMAP1", y="UMAP2", alpha="norm UMI", title=gid) +
    scale_color_brewer(palette = "Set1") +
    scale_alpha_continuous(trans=scales::log10_trans()) +
      theme_bw() +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  ggsave(file.path(f_out, sprintf("%s.pdf", gid)), width = 11, height = 8, units = "cm")
}

gene_li <- rownames(cnt_mat)
gene_li <- gene_li[grep(chrom, gene_li)]
bplapply(gene_li, plot_umap, umap_df=umap_df, cnt_mat=cnt_mat, f_out=f_out, BPPARAM=MulticoreParam(28))



