library(scran)
library(scater)

args = commandArgs(TRUE)
if (length(args) == 5) {
  f_corSCE <- args[1]
  npcs <- as.integer(args[2])
  SNN_k <- as.integer(args[3])
  SNN_type <- args[4]
  f_out <- args[5]
} else {
  q()
}

# SNN_k = 5
# SNN_type = "rank"

load(f_corSCE)
corSCE <- runPCA(corSCE, exprs_values="corrected", ncomponents=npcs)
SNN_graph <- buildSNNGraph(corSCE, use.dimred="PCA", k=SNN_k, type=SNN_type)
save(SNN_graph, file=f_out)