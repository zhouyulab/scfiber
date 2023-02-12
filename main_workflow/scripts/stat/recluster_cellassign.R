library(readr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(cellassign)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_marker <- args[1]
  f_mergeSCE <- args[2]
  f_corSCE <- args[3]
  f_out <- args[4]
} else {
  q()
}
data(example_marker_mat)

# f_marker <- "/home/sqreb/ias/data/CottonSingleCell/data/modularity_marker_gene.tsv"
marker_df <- read_delim(f_marker, "\t", escape_double = FALSE, trim_ws = TRUE)
marker_mat <- as.matrix(marker_df[,2:ncol(marker_df)])
rownames(marker_mat) <- marker_df$Gene

load(f_mergeSCE)

size_factor <- colSums(SummarizedExperiment::assay(mergeSCE, "counts"))

mergeSCE <- mergeSCE[rownames(marker_mat),]

load(f_corSCE)

fit <- cellassign(mergeSCE, marker_gene_info=marker_mat, s = size_factor, return_SCE=TRUE)
save(fit, file = f_out)