library(Seurat)
library(dplyr)
library(argparse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input")
parser$add_argument("-s", "--sample", nargs="*", required=TRUE, type="character", dest = "sample", metavar="sample")
parser$add_argument("-r", "--rep", nargs="*", required=TRUE, type="character", dest = "rep", metavar="rep")
parser$add_argument("-t", "--tech", nargs="*", required=TRUE, type="character", dest = "tech", metavar="tech")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
parser$add_argument("--n-feature", required=FALSE, type="integer", dest = "n_feature", metavar="n_feature", default=2000)
parser$add_argument("--caa-dim", required=FALSE, type="integer", dest = "cca_dim", metavar="cca_dim", default=50)
args <- commandArgs(TRUE)

# args <- c(
#   "-i",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/WT_rep2/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/WT_rep3/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/FL_rep2/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/FL_rep3/filter_feature_bc_matrix",
#   "-s", "WT", "WT", "FL", "FL",
#   "-r", "rep2", "rep3", "rep2", "rep3",
#   "-t", "v2_lib", "v3_lib", "v2_lib", "v3_lib",
#   "-o", "~/"
# )

args <- parser$parse_args(args) 
file_num <- length(args$input)
if((length(args$sample)!=file_num) | (length(args$rep)!=file_num) | (length(args$tech)!=file_num)) stop()

sample_name <- sprintf("%s.%s", args$sample, args$rep)

print("Loading files ...")
sc_count_expr_li <- list()
for(indx in 1:file_num){
  tmp_name <- sample_name[indx]
  sc_cnt <- Read10X(data.dir = args$input[indx])
  sc_count_expr_li[[tmp_name]] <- CreateSeuratObject(counts = sc_cnt, project = tmp_name)
  sc_count_expr_li[[tmp_name]]$Tech <- args$tech[indx]
  sc_count_expr_li[[tmp_name]]$Sample <- args$sample[indx]
  sc_count_expr_li[[tmp_name]]$Rep <- args$rep[indx]
  sc_count_expr_li[[tmp_name]]$Batch <- tmp_name
  # select_indx <- sample(1:ncol(sc_count_expr_li[[tmp_name]]), 2000)
  # sc_count_expr_li[[tmp_name]] <- sc_count_expr_li[[tmp_name]][,select_indx]
  rm(sc_cnt)
  gc()
}

print("Log-normalization and find HVG")
norm_find_hvg <- function(sc){
  sc <- NormalizeData(sc, verbose = FALSE)
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = args$n_feature, verbose = FALSE)
  return(sc)
}
sc_count_expr_li <- lapply(sc_count_expr_li, norm_find_hvg)


print("Find integration anchors")
anchor <- FindIntegrationAnchors(
  object.list = sc_count_expr_li, dims = 1:args$cca_dim
)

print("Integrate")
integrated <- IntegrateData(anchorset = anchor, dims = 1:args$cca_dim)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)

print("runUMAP")
integrated <- RunPCA(integrated, npcs = args$cca_dim, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims=1:args$cca_dim)

corSCE <- integrated
save(corSCE, file=file.path(args$output, "corSCE.RData"))
rm(corSCE)
gc()

col_li <- list()
col_li[["Batch"]] <- brewer.pal(length(sample_name), "Set2")
names(col_li[["Batch"]]) <- sample_name
uniq_tech <- sort(unique(args$tech))
col_li[["Tech"]] <- brewer.pal(max(3, length(uniq_tech)), "Set1")[1:length(uniq_tech)]
names(col_li[["Tech"]]) <- uniq_tech
uniq_sample <- sort(unique(args$sample))
col_li[["Sample"]] <- brewer.pal(max(3, length(uniq_sample)), "Dark2")[1:length(uniq_sample)]
names(col_li[["Sample"]]) <- uniq_sample
uniq_rep <- sort(unique(args$rep))
col_li[["Rep"]] <- brewer.pal(max(3, length(uniq_rep)), "Accent")[1:length(uniq_rep)]
names(col_li[["Rep"]]) <- uniq_rep


p1 <- DimPlot(integrated, reduction = "umap", group.by = "Batch", pt.size=0.1, cols = col_li[["Batch"]])
p2 <- DimPlot(integrated, reduction = "umap", group.by = "Sample", pt.size=0.1, cols = col_li[["Sample"]])
p3 <- DimPlot(integrated, reduction = "umap", group.by = "Tech", pt.size=0.1, cols = col_li[["Tech"]])
p4 <- DimPlot(integrated, reduction = "umap", group.by = "Rep", pt.size=0.1, cols = col_li[["Rep"]])


inche_cm=2.54
pdf(file.path(args$output, "Corrected.pdf"), width=36/inche_cm, height=30/inche_cm)
plot_grid(p1, p2, p3, p4, nrow=2)
dev.off()


print("Merge no correction data")
no_cor_SCE <- sc_count_expr_li[[sample_name[1]]]
for(indx in 2:file_num){
  no_cor_SCE <-  merge(x=no_cor_SCE, y=sc_count_expr_li[[sample_name[indx]]])
}
no_cor_SCE <- FindVariableFeatures(no_cor_SCE, selection.method = "vst", nfeatures = args$n_feature, verbose = FALSE)
no_cor_SCE <- ScaleData(no_cor_SCE, verbose = FALSE)
no_cor_SCE <- RunPCA(no_cor_SCE, npcs = args$cca_dim, verbose = FALSE)
no_cor_SCE <- RunUMAP(no_cor_SCE, reduction = "pca", dims=1:args$cca_dim)

noCorSCE <- no_cor_SCE
save(noCorSCE, file=file.path(args$output, "noCorSCE.RData"))
rm(noCorSCE)
gc()

p1 <- DimPlot(no_cor_SCE, reduction = "umap", group.by = "Batch", pt.size=0.1, cols = col_li[["Batch"]])
p2 <- DimPlot(no_cor_SCE, reduction = "umap", group.by = "Sample", pt.size=0.1, cols = col_li[["Sample"]])
p3 <- DimPlot(no_cor_SCE, reduction = "umap", group.by = "Tech", pt.size=0.1, cols = col_li[["Tech"]])
p4 <- DimPlot(no_cor_SCE, reduction = "umap", group.by = "Rep", pt.size=0.1, cols = col_li[["Rep"]])

inche_cm=2.54
pdf(file.path(args$output, "NoCorrected.pdf"), width=36/inche_cm, height=30/inche_cm)
plot_grid(p1, p2, p3, p4, nrow=2)
dev.off()
