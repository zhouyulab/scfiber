library(readr)
library(scran)

args = commandArgs(TRUE)
if (length(args) == 7) {
  f_SNN <- args[1]
  f_mergeSCE <- args[2]
  expr_cutoff <- as.numeric(args[3])
  min_cell_cutoff <- as.integer(args[4])
  max_cell_ratio <- as.numeric(args[5])
  ncore <- as.integer(args[6])
  f_out <- args[7]
} else {
  q()
}

# f_SNN <- "/home/sqreb/ias/data/CottonSingleCell/analysis/find_params/cluster_params/expr/SNN/correction_k_mnn_30/SNN_D_50_SNN_K_50_SNN_TYPE_number.RData"
# f_mergeSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/find_params/correction_params/MNN_norm_HVG_correction/k_mnn_30/sigma_mnn_0.3/mergeSCE.RData"
# expr_cutoff <- 2
# min_cell_cutoff <- 100
# max_cell_ratio <- 0.5

load(f_mergeSCE)
cnt_mat <- normcounts(mergeSCE)
rm(mergeSCE)
gc()

expr_cell_num <- rowSums(cnt_mat>expr_cutoff)
expr_cell_ratio <- expr_cell_num / ncol(cnt_mat)
filter_cnt_mat <- cnt_mat[(expr_cell_num>min_cell_cutoff) & (expr_cell_ratio<max_cell_ratio),]
filter_cnt_mat <- as.matrix(filter_cnt_mat)
rm(cnt_mat)
gc()

load(f_SNN)

compute_modularity <- function(gid, filter_cnt_mat, expr_cutoff, SNN_graph){
  is_expr <- filter_cnt_mat[which(rownames(filter_cnt_mat)==gid),] > expr_cutoff
  modularity_ratio <- clusterModularity(SNN_graph, is_expr, as.ratio=TRUE)
  df <- data.frame(Gid=gid, ModularityRatio=modularity_ratio[2,2], ExprCellNum=sum(is_expr), ExprCellRatio=mean(is_expr))
  return(df)
}

print("Compute modularity")
res <- lapply(rownames(filter_cnt_mat), compute_modularity, filter_cnt_mat=filter_cnt_mat, expr_cutoff=expr_cutoff, SNN_graph=SNN_graph)
res_df <- plyr::ldply(res, function(x) x)
write_tsv(res_df, f_out)