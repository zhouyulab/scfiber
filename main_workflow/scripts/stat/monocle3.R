library(monocle3)
library(SingleCellExperiment)
library(readr)
library(dplyr)
library(ggplot2)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_corSCE <- args[1]
  f_cluster <- args[2]
  f_data <- args[3]
  f_plot <- args[4]
} else {
  q()
}

f_corSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/corSCE.RData"
f_cluster <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/cluster.tsv"


load(f_corSCE)
clusterSCE <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)
corSCE$CellCluster <- clusterSCE$Cluster

cell_meta <- as.data.frame(colData(corSCE))
cell_meta$Cell <- sprintf("%s_%s_%s", cell_meta$Sample, cell_meta$Rep, cell_meta$Barcode)
rownames(cell_meta) <- cell_meta$Cell

gene_meta <- as.data.frame(rowData(corSCE))
gene_meta$gene_short_name <- gene_meta$Symbol
expr_mat <- as.matrix(corSCE@assays@data@listData$corrected)
colnames(expr_mat) <- cell_meta$Cell

cds <- new_cell_data_set(expr_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_meta)


cds <- preprocess_cds(cds, num_dim = 50, norm_method="none")
cds <- reduce_dimension(cds, reduction_method="UMAP", umap.n_neighbors = 25, umap.min_dist=0.3)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)


plot_cells(cds,
           reduction_method = "UMAP",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

p <- plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "CellCluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

ggsave(f_plot, p, height = 10, width = 10, units = "cm")
save(cds, file=f_data)
