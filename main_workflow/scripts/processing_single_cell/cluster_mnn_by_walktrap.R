library(readr)
library(scran)
library(scater)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_SNN_graph <- args[1]
  f_corSCE <- args[2]
  walktrap_step <- as.integer(args[3])
  f_out <- args[4]
} else {
  q()
}

# walktrap_step = 4

load(f_SNN_graph)

cluster <- igraph::cluster_walktrap(SNN_graph, steps=walktrap_step)$membership 

cell_num <- table(cluster)
clu_num <- as.data.frame(table(cluster))
clu_num <- clu_num[order(clu_num$Freq, decreasing = TRUE),]
clu_num$label <- sprintf("C%d", seq(1: nrow(clu_num)))
cluster <- factor(cluster, levels = clu_num$cluster, labels = clu_num$label)

load(f_corSCE)
df <- data.frame(Cell=colnames(corSCE), Cluster=cluster)
write_tsv(df, f_out)