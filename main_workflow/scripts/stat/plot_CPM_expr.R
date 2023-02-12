library(ggplot2)
library(dplyr)
library(readr)
library(reshape2)
library(BiocParallel)
args = commandArgs(TRUE)
if (length(args) == 2) {
  f_CPM_df <- args[1]
  f_out <- args[2]
} else {
  q()
}

f_CPM_df <- "analysis_v3/stat/cluster_cnt/cluster.tpm.tsv"
f_out <- "analysis_v3/stat/cluster_cnt/cluster_expr_plot"

CPM_df <- read_delim(f_CPM_df, "\t", escape_double = FALSE, trim_ws = TRUE)
CPM_df <- CPM_df[apply(CPM_df[2:ncol(CPM_df)], 1, function(x){return(max(x, na.rm = T))})>0,]
melt_CPM_df <- melt(CPM_df, id.vars = "Gid", variable.name = "Label", value.name = "CPM")
melt_CPM_df$Sample <- sapply(strsplit(as.character(melt_CPM_df$Label), "[.]"), function(x){return(x[1])})
melt_CPM_df$Cluster <- sapply(strsplit(as.character(melt_CPM_df$Label), "[.]"), function(x){return(x[2])})
melt_CPM_df$Sample <- factor(melt_CPM_df$Sample, levels = c("WT", "FL"))
melt_CPM_df$Cluster <- factor(melt_CPM_df$Cluster, levels = c("C1", "C2", "C3", "C4", "C5"))
melt_CPM_df$Chrom <- substr(melt_CPM_df$Gid, 1, 8)
melt_CPM_df$Chrom[!melt_CPM_df$Chrom%in%c(sprintf("Ghir_A%02d", 1:13), sprintf("Ghir_D%02d", 1:13))] <- "Other"

chrom_li <- unique(melt_CPM_df$Chrom)
plot_expr <- function(chrom, f_out, melt_CPM_df){
  tmp_dir <- file.path(f_out, chrom)
  if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir, recursive = T)
  }
  tmp_CPM_df <- melt_CPM_df[melt_CPM_df$Chrom == chrom,]
  gid_li <- unique(tmp_CPM_df$Gid)
  for(gid in gid_li){
    tmp_gene_cpm_df <- tmp_CPM_df[tmp_CPM_df$Gid==gid,]
    p <- ggplot(tmp_gene_cpm_df, aes(x=Cluster, y=Sample, color=Sample, size=CPM)) +
      geom_point() +
      scale_color_brewer(palette = "Set1") +
      scale_size_continuous(range = c(0, 6)) +
      theme_bw() +
      theme(
        text = element_text(family="ArialMT", size=5),
        title = element_text(family="ArialMT", size=5),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(family="ArialMT", size=5),
        legend.background = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(family="ArialMT", size=5),
        panel.grid = element_blank()
      )
    ggsave(file.path(tmp_dir, sprintf("%s.pdf", gid)), p, width = 10, height = 4, units = "cm")
  }
}
bplapply(chrom_li, plot_expr, f_out=f_out, melt_CPM_df=melt_CPM_df, BPPARAM=MulticoreParam(40))
