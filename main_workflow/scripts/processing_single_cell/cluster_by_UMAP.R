library(monocle3)
library(SingleCellExperiment)
library(readr)
library(dplyr)
library(ggplot2)



# f_corSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/corSCE.RData"

load(f_corSCE)
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
cds@int_colData@listData$reducedDims@listData$UMAP <- corSCE@int_colData@listData$reducedDims@listData$UMAP[,1:2]
cds <- cluster_cells(cds)

p <- plot_cells(cds,
           reduction_method = "UMAP",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

ggsave(f_plot, p, height = 10, width = 10, units = "cm")

df <- data.frame(Sample=cds$SampleName, Rep=cds$Rep, Cluster=clusters(cds), Barcode=cds$Barcode)
write_tsv(df, f_out)
