library(monocle3)
library(SingleCellExperiment)
library(ggplot2)
library(Matrix)
library(readr)
library(dplyr)

setwd("~/project/CottonSCE")

load("~/project/CottonSCE/data/final_cluster/umapSCE.RData")
load("~/project/CottonSCE/data/final_cluster/mergeSCE.RData")
out_dir <- "~/project/CottonSCE/analysis/kmer_graph_enrich"
RNA_cnt_mat <- counts(mergeSCE)
RNA_bool_mat <- RNA_cnt_mat>1
RNA_cnt_mat <- as(RNA_cnt_mat, "dgCMatrix")
RNA_bool_mat <- as(RNA_bool_mat, "dgCMatrix")


kmer_li <- c(8)

kmer2cds <- function(kmer_mat){
  cell_meta <- as.data.frame(colData(corSCE))
  cell_meta$Cell <- sprintf("%s_%s_%s", cell_meta$Sample, cell_meta$Rep, cell_meta$Barcode)
  rownames(cell_meta) <- cell_meta$Cell
  
  kmer_meta <- data.frame(Symbol=rownames(kmer_mat), gene_short_name=rownames(kmer_mat))
  rownames(kmer_meta) <- rownames(kmer_mat)
  
  expr_mat <- as.matrix(kmer_mat)
  colnames(expr_mat) <- cell_meta$Cell
  rm(kmer_mat)
  gc()
  
  print("Create CDS")
  cds <- new_cell_data_set(expr_mat,
                           cell_metadata = cell_meta,
                           gene_metadata = kmer_meta)
  rm(expr_mat)
  
  gc()
  print("Preprocess CDS")
  cds <- preprocess_cds(cds, num_dim = 50)
  print("Reduce dimension")
  cds <- reduce_dimension(cds, reduction_method="UMAP")
  cds@int_colData@listData$reducedDims$UMAP[,1] <- corSCE@int_colData@listData$reducedDims$UMAP[,1]
  cds@int_colData@listData$reducedDims$UMAP[,2] <- corSCE@int_colData@listData$reducedDims$UMAP[,2]
  print("Cluster cells")
  cds <- cluster_cells(cds)
  print("Learn graph")
  cds <- learn_graph(cds)
  return(cds)
}

for (kmer in kmer_li) {
  kmer_cnt_gene <- read_delim(sprintf("data/kmer_motif_cnt_MACS2/kmer.%s.cnt.gene.tsv", kmer), "\t", escape_double = FALSE, trim_ws = TRUE)
  kmer_gene_mat <- as.matrix(kmer_cnt_gene[,3:ncol(kmer_cnt_gene)])
  kmer_df <- kmer_cnt_gene[,c("Motif", "RevMotif")]
  write_tsv(kmer_df, file.path("analysis", "kmer_graph_enrich", sprintf("kmer.%s.tsv", kmer)))
  rownames(kmer_gene_mat) <- kmer_df$Motif
  kmer_gene_bool_mat <- kmer_gene_mat>0
  kmer_gene_bool_mat <- as(kmer_gene_bool_mat, "dgCMatrix")
  tmp_RNA_cnt_mat <- RNA_cnt_mat[colnames(kmer_gene_bool_mat),]
  kmer_cell_RNA_cnt_mat <- kmer_gene_bool_mat %*% tmp_RNA_cnt_mat
  tmp_RNA_bool_mat <- RNA_bool_mat[colnames(kmer_gene_bool_mat),]
  kmer_cell_RNA_bool_mat <- kmer_gene_bool_mat %*% tmp_RNA_bool_mat
  kmer_cell_li <- list("cnt_mat"=kmer_cell_RNA_cnt_mat, "bool_mat"=kmer_cell_RNA_bool_mat)
  save(kmer_cell_li, file=file.path("analysis", "kmer_graph_enrich", sprintf("kmer.%s.cell_RNA_mat.RData", kmer)))

  kmer_cell_cds_li <- lapply(kmer_cell_li, kmer2cds)
  save(kmer_cell_cds_li, file=file.path("analysis", "kmer_graph_enrich", sprintf("kmer.%s.cell_RNA_mat.cds.RData", kmer)))
  load(file.path("analysis", "kmer_graph_enrich", sprintf("kmer.%s.cell_RNA_mat.cds.RData", kmer)))
  enrich_graph_li <- lapply(kmer_cell_cds_li, graph_test)
  cnt_graph_res <- enrich_graph_li$cnt_mat
  cnt_graph_res <- cnt_graph_res[,c("Symbol", "morans_I")]
  names(cnt_graph_res) <- c("Motif", "UMI_morans_I")
  write_tsv(cnt_graph_res, file.path("analysis", "kmer_graph_enrich", sprintf("kmer.%s.expr_graph_res.tsv", kmer)))
  bool_graph_res <- enrich_graph_li$bool_mat
  bool_graph_res <- bool_graph_res[,c("Symbol", "morans_I")]
  names(bool_graph_res) <- c("Motif", "ExprGeneNum_morans_I")
  write_tsv(bool_graph_res, file.path("analysis", "kmer_graph_enrich", sprintf("kmer.%s.bool_graph_res.tsv", kmer)))
  graph_enrich_df <- left_join(kmer_df, cnt_graph_res)
  graph_enrich_df <- left_join(graph_enrich_df, bool_graph_res)
  write_tsv(graph_enrich_df, file.path("analysis", "kmer_graph_enrich", sprintf("kmer.%s.graph_enrich_score.tsv", kmer)))
}


