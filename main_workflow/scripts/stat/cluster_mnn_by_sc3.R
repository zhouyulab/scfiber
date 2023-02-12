library(readr)
library(SC3)
library(scater)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_corSCE <- args[1]
  k_sc3 <- as.integer(args[2])
  n_cores <- as.integer(args[3])
  f_out <- args[4]
} else {
  q()
}

load(f_corSCE)

build_cluster <- function(sc_expr, k_sc3, n_cores){
  cnt_mat <- normcounts(sc_expr)
  log_cnt_mat <- logcounts(sc_expr)                                                                  
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(cnt_mat),
      logcounts = as.matrix(log_cnt_mat)
    ),
    colData = colData(sc_expr)
  )
  colnames(sce) <- sprintf("CELL%d", 1:ncol(cnt_mat))
  rowData(sce)$feature_symbol <- rownames(cnt_mat)
  sce <- sc3(sce, ks = k_sc3, n_cores=n_cores)
  cluster <- colData(sce)[,sprintf("sc3_%d_clusters", k_sc3)]
  cell_num <- table(cluster)
  clu_num <- as.data.frame(table(cluster))
  clu_num <- clu_num[order(clu_num$Freq, decreasing = TRUE),]
  clu_num$label <- sprintf("C%d", seq(1: nrow(clu_num)))
  cluster <- factor(cluster, levels = clu_num$cluster, labels = clu_num$label)
  return(cluster)
}

cluster <- build_cluster(corSCE, k_sc3, n_cores)
df <- data.frame(Sample=corSCE$SampleName, Rep=corSCE$Rep, Cluster=cluster, Barcode=colnames(corSCE))
write_tsv(df, f_out)