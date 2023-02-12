library(readr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(SingleCellExperiment)

args = commandArgs(TRUE)
if (length(args) == 6) {
  f_in <- args[1]
  f_mergeSCE <- args[2]
  f_corSCE <- args[3]
  f_cluster <- args[4]
  ratio_cutoff <- as.numeric(args[5])
  f_out <- args[6]
} else {
  q()
}

# f_in <- "/home/sqreb/ias/data/CottonSingleCell/analysis_v2/recluster/modularity/modularity.tsv"
# ratio_cutoff <- 12

modularity_df <- read_delim(f_in, "\t", escape_double = FALSE, trim_ws = TRUE)
modularity_df <- modularity_df[order(modularity_df$ModularityRatio, decreasing = TRUE),]
modularity_df$Rank <- 1:nrow(modularity_df)
modularity_df$Pass <- modularity_df$ModularityRatio>ratio_cutoff

modularity_df$Pass <- factor(modularity_df$Pass, levels = c(TRUE, FALSE), labels = c("Pass", "Failed"))
label_text <- sprintf("Pass: %d\nFailed: %d", sum(modularity_df$Pass=="Pass"), sum(modularity_df$Pass=="Failed"))

p <- ggplot(modularity_df, aes(x=Rank, y=log10(ModularityRatio), color=Pass)) +
  theme_bw() +
  geom_point(size=0.3) +
  geom_hline(yintercept = log10(ratio_cutoff), color="red") +
  annotate("text", x=2500, y=2.0, label=label_text, size=1.4, hjust=0) +
  annotate("text", x=2500, y=log10(ratio_cutoff), label=sprintf("Cutoff = %s", ratio_cutoff), size=1.4, hjust=0, vjust=-1) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), labels = c("1", "10", "100", "1,000")) +
  scale_color_manual(values = c("Pass"="black", "Failed"="grey70")) +
  labs(x="Rank", y="Modularity ( Obs. / Exp. )") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        panel.grid = element_blank())
ggsave(file.path(f_out, "Modularity.Rank.pdf"), p, height = 4, width = 6, units = "cm", colormodel = "cmyk")

p <- ggplot(modularity_df, aes(x=log10(ModularityRatio), y=ExprCellRatio, color=Pass)) +
  theme_bw() +
  geom_point(size=0.3) +
  annotate("text", x=0.85*max(log10(modularity_df$ModularityRatio)), y=0.95*max(modularity_df$ExprCellRatio), label=label_text, size=1.4, hjust=0, vjust=1) +
  annotate("text", x=log10(ratio_cutoff), y=0.95*max(modularity_df$ExprCellRatio), label=sprintf("Cutoff = %s", ratio_cutoff), size=1.4, hjust=-0.2, vjust=1) +
  geom_vline(xintercept = log10(ratio_cutoff), color="red") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20%", "40%", "60%", "80%", "100%")) +
  scale_x_continuous(breaks = c(0, 1, 2, 3), labels = c("1", "10", "100", "1,000")) +
  scale_color_manual(values = c("Pass"="black", "Failed"="grey70")) +
  labs(x="Modularity ( Obs. / Exp. )", y="Expressed cell") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        panel.grid = element_blank())
ggsave(file.path(f_out, "Modularity.Ratio.pdf"), p, height = 4, width = 6, units = "cm", colormodel = "cmyk")

out_df <- modularity_df[modularity_df$Pass=="Pass",]
out_df$Pass <- NULL
write_tsv(out_df, file.path(f_out, "Modularity.filter.tsv"))


load(f_mergeSCE)
load(f_corSCE)
mergeSCE <- mergeSCE[rownames(mergeSCE) %in% out_df$Gid,]

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

UMAP_folder <- file.path(f_out, "MarkerGene_UMAP")
if(!file.exists(UMAP_folder)){
  dir.create(UMAP_folder)
}

for(indx in 1:nrow(out_df)){
  ps_gene_marker <- list()
  gid <- out_df$Gid[indx]
  ratio <- out_df$ModularityRatio[indx]
  label <- sprintf("%s: %s (Modularity ratio:%.2f)", "All", gid, ratio)
  ps_gene_marker[[label]] <-  plot_gene_expr_umap(corSCE, mergeSCE, gid, label)
  for(s in sort(unique(corSCE$SampleName))){
    tmpUMAP <- corSCE[,corSCE$SampleName==s]
    tmpSCE <- mergeSCE[,mergeSCE$SampleName==s]
    label <- sprintf("%s: %s (Modularity ratio:%.2f)", s, gid, ratio)
    ps_gene_marker[[label]] <-  plot_gene_expr_umap(tmpUMAP, tmpSCE, gid, label)
    rm(tmpUMAP)
    rm(tmpSCE)
    gc()
  }
  gc()
  inche_cm=2.54
  pdf(file.path(UMAP_folder, sprintf("%s.pdf", gid)), width=8*3/inche_cm, height=6*ceiling(length(ps_gene_marker)/3)/inche_cm, family="ArialMT", colormodel = "cmyk")
  grid.arrange(grobs=ps_gene_marker,ncol=3,padding = unit(0, "mm"))
  dev.off()
}

