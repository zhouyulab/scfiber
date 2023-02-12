library(readr)
library(batchelor)
library(ggplot2)
library(gridExtra)

args = commandArgs(TRUE)
if (length(args) == 5) {
  f_countSCE <- args[1]
  f_corSCE <- args[2]
  f_noCorSCE <- args[3]
  f_marker_gene <- args[4]
  f_out <- args[5]
} else {
  q()
}

# f_countSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/v2_analysis/correction_params/MNN_norm_HVG_correction/k_mnn_20/sigma_mnn_0.02/countSCE.RData"
# f_corSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/v2_analysis/correction_params/MNN_norm_HVG_correction/k_mnn_20/sigma_mnn_0.02/corSCE.RData"
# f_noCorSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/v2_analysis/correction_params/MNN_norm_HVG_correction/k_mnn_20/sigma_mnn_0.02/noCorSCE.RData"
# f_marker_gene <- "/home/sqreb/ias/data/CottonSingleCell/data/protential_marker_gene.txt"

merker_gene <- as.character(unlist(read_delim(f_marker_gene, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)))

load(f_countSCE)
load(f_corSCE)
load(f_noCorSCE)

plot_gene_expr_umap <- function(umapSCE, SCE, gid, label){
  mat <- counts(SCE)
  expr <- mat[which(rownames(SCE)==gid),]
  if(! all(names(expr)==colnames(umapSCE))) stop()
  df <- data.frame(x=umapSCE@int_colData$reducedDims$UMAP[,1], y=umapSCE@int_colData$reducedDims$UMAP[,2], expr=expr)
  cell_num <- nrow(df)
  df <- df[sample(1:cell_num, cell_num, replace = FALSE),]
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

no_correction <- function(...){
  return(correctExperiments(..., PARAM=NoCorrectParam()))
}

mergeSCE <- do.call(no_correction, countSCE)

sample_li <- sort(unique(mergeSCE$Sample))

inche_cm=2.54
pdf(f_out, width=30/inche_cm, height=12/inche_cm, family="ArialMT", colormodel = "cmyk")

for(gid in merker_gene){
  ps <- list()
  
  ## No Correction
  for(s in sample_li){
    tmpSCE <- mergeSCE[,mergeSCE$Sample==s]
    tmpUMAP <- noCorSCE[,noCorSCE$Sample==s]
    label <- sprintf("No correction:%s", s)
    ps[[label]] <- plot_gene_expr_umap(tmpUMAP, tmpSCE, gid, label)
  }
  
  ## Correction
  for(s in sample_li){
    tmpSCE <- mergeSCE[,mergeSCE$Sample==s]
    tmpUMAP <- corSCE[,corSCE$Sample==s]
    label <- sprintf("Correction:%s", s)
    ps[[label]] <- plot_gene_expr_umap(tmpUMAP, tmpSCE, gid, label)
  }
  grid.arrange(grobs=ps,nrow=2,padding = unit(0, "mm"))
}
dev.off()

