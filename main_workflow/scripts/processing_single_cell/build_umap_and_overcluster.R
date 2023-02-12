library(monocle3)
library(SingleCellExperiment)
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)

args = commandArgs(TRUE)
if (length(args) == 9) {
  f_mergeSCE <- args[1]
  f_corSCE <- args[2]
  k <- as.integer(args[3])
  sigma <- as.double(args[4])
  f_umapSCE <- args[5]
  f_monoSCE <- args[6]
  f_clu <- args[7]
  f_plot <- args[8]
  f_marker <- args[9]
} else {
  q()
}

  # base_dir = "/home/sqreb/data/CottonSingleCell/analysis/find_params/correction/k_mnn_25_sigma_mnn_0.3"
  # f_mergeSCE = file.path(base_dir, "mergeSCE.RData")
  # f_corSCE = file.path(base_dir, "corSCE.RData")
  # k = 25
  # sigma = 0.3
  # f_umapSCE = file.path(base_dir, "umapSCE.RData")
  # f_monoSCE = file.path(base_dir, "monoSCE.RData")
  # f_clu = file.path(base_dir, "overcluster.tsv")
  # f_plot = file.path(base_dir, "umap.pdf")
  # f_marker = file.path(base_dir, "marker.pdf")

load(f_mergeSCE)
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
cds <- reduce_dimension(cds, reduction_method="UMAP", umap.n_neighbors = k, umap.min_dist=sigma)
cds <- cluster_cells(cds, k=k, num_iter=5)
cds <- learn_graph(cds)


# 
# inche_cm=2.54
# pdf("~/test.pdf", width=12/inche_cm, height=9/inche_cm, family="ArialMT", colormodel = "cmyk")
# plot_cells(cds,
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE)
# dev.off()

inche_cm=2.54
pdf(f_plot, width=12/inche_cm, height=9/inche_cm, family="ArialMT", colormodel = "cmyk")
plot_cells(cds,
           reduction_method = "UMAP",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "SampleName",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by = "Rep",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()

UMAP <- cds@int_colData@listData$reducedDims@listData$UMAP

corSCE@int_colData@listData$reducedDims@listData$UMAP <- UMAP
save(corSCE, file=f_umapSCE)
over_cluster <- data.frame(Sample=cds$SampleName, Rep=cds$Rep, Cluster=clusters(cds), Barcode=cds$Barcode)
write_tsv(over_cluster, f_clu)
dim(UMAP)

save(cds, file=f_monoSCE)

c(
  "Ghir_D13G002250", "Ghir_A02G002840", "Ghir_A01G021230", "Ghir_D05G008850", 
  "Ghir_D08G000460", "Ghir_D08G000460", "Ghir_D08G016080", "Ghir_D05G006460", 
  "Ghir_D05G019650", "Ghir_D05G038790", "Ghir_D10G004420"
  )

df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_A01G021230",]>2)
p1 <- ggplot(df, aes(x=x, y=y)) +
  labs(title="Ghir_A01G021230", x="UMAP1", y="UMAP2") +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3, size=0.2) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3, size=0.2)
p1

df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_D08G000460",]>2)
p2 <- ggplot(df, aes(x=x, y=y)) +
  labs(title="Ghir_D08G000460", x="UMAP1", y="UMAP2") +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3, size=0.2) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3, size=0.2)

df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_D05G019650",]>2)
p3<- ggplot(df, aes(x=x, y=y, color=cnt)) +
  labs(title="Ghir_D05G019650", x="UMAP1", y="UMAP2") +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3, size=0.2) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3, size=0.2)

df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_D10G004420",]>2)
p4<- ggplot(df, aes(x=x, y=y, color=cnt)) +
  labs(title="Ghir_D10G004420", x="UMAP1", y="UMAP2") +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3, size=0.2) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3, size=0.2)

df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_D10G012320",]>2)
p5<- ggplot(df, aes(x=x, y=y, color=cnt)) +
  labs(title="Ghir_D10G012320", x="UMAP1", y="UMAP2") +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3, size=0.2) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3, size=0.2)

df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_A03G015850",]>2)
p6<- ggplot(df, aes(x=x, y=y, color=cnt)) +
  labs(title="Ghir_A03G015850", x="UMAP1", y="UMAP2") +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3, size=0.2) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3, size=0.2)

ps <- list(p2=p2, p3=p3, p4=p4, p5=p5, p6=p6)
inche_cm=2.54
pdf(f_marker, width=12*3/inche_cm, height=9*ceiling(length(ps)/3)/inche_cm, family="ArialMT", colormodel = "cmyk")
grid.arrange(grobs=ps,ncol=3,padding = unit(0, "mm"))
dev.off()


df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_D05G006460",]>2)
p<- ggplot(df, aes(x=x, y=y, color=cnt)) +
  labs(title="Ghir_D05G006460") +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3)
  
p



df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_D05G038790",]>2)
p<- ggplot(df, aes(x=x, y=y, color=cnt)) +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3)
p
df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_D07G007580",]>2)
p<- ggplot(df, aes(x=x, y=y, color=cnt)) +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3)


p


df <- data.frame(x=UMAP[,1], y=UMAP[,2], cnt=mergeSCE@assays@data@listData$normcounts["Ghir_D04G020730",]>1.5)
p<- ggplot(df, aes(x=x, y=y, color=cnt)) +
  labs(title="Ghir_D04G020730", x="UMAP1", y="UMAP2") +
  geom_point(data=df[!df$cnt,], color="grey70", alpha=0.3, size=0.5) +
  geom_point(data=df[df$cnt,], color="red", alpha=0.3, size=2)
p


