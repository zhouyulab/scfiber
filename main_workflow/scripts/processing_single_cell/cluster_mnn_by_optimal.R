library(readr)
library(scran)
library(scater)

args = commandArgs(TRUE)
if (length(args) == 3) {
  f_SNN_graph <- args[1]
  f_corSCE <- args[2]
  f_out <- args[3]
} else {
  q()
}


load(f_SNN_graph)
cluster <- igraph::cluster_optimal(SNN_graph)$membership 

cell_num <- table(cluster)
clu_num <- as.data.frame(table(cluster))
clu_num <- clu_num[order(clu_num$Freq, decreasing = TRUE),]
clu_num$label <- sprintf("C%d", seq(1: nrow(clu_num)))
cluster <- factor(cluster, levels = clu_num$cluster, labels = clu_num$label)

load(f_corSCE)
df <- data.frame(Cell=colnames(corSCE), Cluster=cluster)
write_tsv(df, f_out)