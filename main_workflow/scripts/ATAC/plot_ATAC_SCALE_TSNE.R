library(ggplot2)
library(readr)
load_tsne <- function(f_tsne, f_clu, f_plot){
  tsne_df <- read_table2(f_tsne, col_names = FALSE)
  names(tsne_df) <- c("TSNE1", "TSNE2")
  cluster_df <- read_delim(f_clu, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  names(cluster_df) <- c("UMI", "cluster")
  cluster_df$Sample <- sapply(strsplit(cluster_df$UMI, "-"), function(x){return(x[2])})
  cluster_df$Sample <- factor(cluster_df$Sample, levels = c("1", "2"), labels = c("WT", "FL"))
  plot_df <- cbind(cluster_df, tsne_df)
  p <- ggplot(plot_df, aes(x=TSNE1, y=TSNE2, color=Sample)) +
    geom_point(size=0.1) +
    scale_color_manual(values = c("WT"="red", "FL"="blue")) +
    theme_bw() +
    theme(text = element_text(family="ArialMT", color = "black", size = 5),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.title = element_text(family="ArialMT", color = "black", size = 5),
          legend.key.size = unit(3, "mm"),
          panel.grid = element_blank())
  ggsave(f_plot, p, width = 6.5, height = 4.5, units = "cm")
}

load_tsne(
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/WT/tsne.txt",
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/WT/cluster_assignments.txt",
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/WT/ATAC.tsne.pdf"
  )

load_tsne(
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/FL/tsne.txt",
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/FL/cluster_assignments.txt",
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/FL/ATAC.tsne.pdf"
)

load_tsne(
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/merge/tsne.txt",
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/merge/cluster_assignments.txt",
  "~/mu01/project/CottonSingleCell/analysis_v3/SCALE_ATAC/merge/ATAC.tsne.pdf"
)
