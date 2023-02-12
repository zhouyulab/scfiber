library(igraph)
library(readr)
library(dplyr)
library(scPred)
library(tidyverse)
library(Matrix)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_corSCE <- args[1]
  f_mergeSCE <- args[2]
  f_gene_group <- args[3]
  min_expr_cutoff <- as.numeric(args[4])
} else {
  q()
}

gene_group <- read_delim(f_gene_group, "\t", escape_double = FALSE, trim_ws = TRUE)
gene_group <- gene_group[gene_group$Group %in% c("G2", "G3", "G5", "G6"),]

load(f_mergeSCE)
rownames(mergeSCE)
print(gene_group$Gid)
print(rownames(mergeSCE) %in% gene_group$Gid)
mergeSCE <- mergeSCE[rownames(mergeSCE) %in% gene_group$Gid,]
gc()

load(f_corSCE)

norm_mat <- normcounts(mergeSCE)
expr_mat <- norm_mat > min_expr_cutoff

corSCE$Label <- "Not Sure"
cell_groups <- sort(unique(gene_group$Group))
for(g in cell_groups){
  gid_li <- gene_group$Gid[gene_group$Group==g]
  corSCE$Label[which(colSums(expr_mat[rownames(expr_mat) %in% gid_li,]) > 0 )] <- g
}

set.seed(0)
SCE_metadata <- as.data.frame(colData(corSCE))

labeled_corSCE <- corSCE[,corSCE$Label!="Not Sure"]
labeled_norm_mat <- normcounts(labeled_corSCE)
labeled_SCE_metadata <- as.data.frame(colData(labeled_corSCE))

i <- createDataPartition(labeled_SCE_metadata$Label, p = 0.70, list = FALSE)
train_data <-  as.matrix(labeled_norm_mat)[, i]
test_data <- as.matrix(labeled_norm_mat)[, -i]
train_info <- labeled_SCE_metadata[i, , drop = FALSE]
test_info <- labeled_SCE_metadata[-i, , drop = FALSE]

scp <- eigenDecompose(train_data, n = 10)
scPred::metadata(scp) <- train_info

scp <- getFeatureSpace(scp, pVar = "Label")
scp@features
plotEigen(scp, group = "Label")
scp <- trainModel(scp)
res <- getTrainResults(scp)

plotTrainProbs(scp)

scp <- scPredict(scp, newData = test_data, threshold = 0.7)
getPredictions(scp)
scp@predMeta <- test_info
crossTab(scp, true = "Label")
