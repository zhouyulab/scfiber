library(monocle3)
library(SingleCellExperiment)
library(readr)


cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
pt <- pseudotime(cds)
pseudotime(cds)[cds$cluster %in% c("C4", "C5")] <- NA
df <- data.frame(Sample=cds$Sample, Barcode=cds$Barcode, PseudoTime=pt)
readr::write_tsv(df, "~/pt.tsv")


load("/home/sqreb/data/CottonSingleCell/analysis/final_cluster/mergeSCE.RData")

cell_meta <- as.data.frame(colData(mergeSCE))
cell_meta$Cell <- sprintf("%s_%s_%s", cell_meta$Sample, cell_meta$Rep, cell_meta$Barcode)
rownames(cell_meta) <- cell_meta$Cell

gene_meta <- as.data.frame(rowData(mergeSCE))
gene_meta$gene_short_name <- gene_meta$Symbol
expr_mat <- counts(mergeSCE)
colnames(expr_mat) <- cell_meta$Cell

cds2 <- new_cell_data_set(expr_mat,
                          cell_metadata = cell_meta,
                          gene_metadata = gene_meta)
rm(expr_mat)
gc()
cds3 <- cds2
cds2 <- preprocess_cds(cds2, num_dim = 50)
cds2@reduce_dim_aux <- cds@reduce_dim_aux
cds2@principal_graph_aux <- cds@principal_graph_aux
cds2@principal_graph <- cds@principal_graph
cds2@clusters <- cds@clusters
cds2@int_colData <- cds@int_colData

x <- graph_test(cds2, neighbor_graph = "principal_graph", verbose = TRUE, cores = 4)


write_tsv(x, "Graph_enrich_gene.tsv")