library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
setwd("~/project/CottonSCE")

ATAC_neaberby_promoter <- read_delim("data/scATAC_MACS2_cnt/ATAC.neaberby_promoter.100.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_num <- ATAC_neaberby_promoter %>% group_by(ATAC) %>% summarise(ATAC_num=n())
gene_num <- ATAC_neaberby_promoter %>% group_by(Gid) %>% summarise(Gid_num=n())
ATAC_neaberby_promoter <- left_join(ATAC_neaberby_promoter, ATAC_num)
ATAC_neaberby_promoter <- left_join(ATAC_neaberby_promoter, gene_num)
ATAC_neaberby_promoter <- ATAC_neaberby_promoter[ATAC_neaberby_promoter$ATAC_num==1 & ATAC_neaberby_promoter$Gid_num==1,]
ATAC_neaberby_promoter <- ATAC_neaberby_promoter[,c("ATAC", "Gid")]
names(ATAC_neaberby_promoter)[1] <- "Peak"
load_data <- function(f, label, ATAC_neaberby_promoter){
  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)
  df <- left_join(ATAC_neaberby_promoter, df)
  mat <- as.matrix(df[,3:ncol(df)])
  rownames(mat) <- df$Gid
  obj <- CreateSeuratObject(
    mat,
    project = label,
    assay = "ATAC"
  )
  obj$Sample <- label
  return(obj)
}

sample_li <- c("WT", "FL")

gene_li <- list()
for (s in sample_li) {
  gene_li[[s]] <- load_data(sprintf("data/scATAC_MACS2_cnt/%s.peak.cnt.tsv", s), s, ATAC_neaberby_promoter)
  gc()
  # peak_li[[s]] <- load_data(sprintf("data/scATAC_MACS2_cnt/%s.peak.cnt.tsv", s), s)
}

for(s in sample_li){
  gene_li[[s]] <- NormalizeData(gene_li[[s]])
  gene_li[[s]] <- FindVariableFeatures(
    gene_li[[s]], selection.method = "vst",
    nfeatures = 500, verbose = FALSE)
  VariableFeaturePlot(gene_li[[s]])
  gene_li[[s]] <- ScaleData(gene_li[[s]])
  gene_li[[s]] <- RunPCA(gene_li[[s]], features = VariableFeatures(object = gene_li[[s]]))
}

DimPlot(gene_li[["WT"]], reduction = "pca") + DimPlot(gene_li[["FL"]], reduction = "pca")

for(s in sample_li){
  gene_li[[s]] <- JackStraw(gene_li[[s]], num.replicate = 100)
  gene_li[[s]] <- ScoreJackStraw(gene_li[[s]], dims = 1:20)
  JackStrawPlot(gene_li[[s]], dims = 1:15)
  p <- ElbowPlot(gene_li[[s]])
  ggsave(
    plot = p, 
    filename = file.path("analysis", "ATAC_preprocess", sprintf("ATAC.PC.%s.pdf", s)),
    width = 5,
    height = 4,
    units = "cm"
    )
}

for(s in sample_li){
  gene_li[[s]] <- RunUMAP(gene_li[[s]], reduction = "pca", dims = 1:3, n.neighbors=50, min.dist = 0.5)
  gene_li[[s]] <- FindNeighbors(gene_li[[s]], reduction = "umap", dims = 1:2)
  gene_li[[s]] <- FindClusters(gene_li[[s]], resolution = 0.01)
}

p <- DimPlot(gene_li[["WT"]], reduction = "pca") + 
  DimPlot(gene_li[["FL"]], reduction = "pca") + 
  DimPlot(gene_li[["WT"]], reduction = "umap") + 
  DimPlot(gene_li[["FL"]], reduction = "umap")
ggsave(
  plot = p, 
  filename = file.path("analysis", "ATAC_preprocess", "ATAC.gene.umap.pdf"),
  width = 10,
  height = 8,
  units = "cm"
)


FeaturePlot(gene_li[["WT"]], features = c("Ghir-A01G021620", "Ghir-D01G023150", "Ghir-A01G012170", "Ghir-D01G012940","Ghir-D10G024960", "Ghir-A10G022520"), reduction = "umap", pt.size = 0.1, combine = FALSE)

for(s in sample_li){
  gene_li[[s]]$logCnt <- log10(gene_li[[s]]$nCount_RNA + 1)
}

df <- data.frame(x=gene_li[[s]]@reductions$umap@cell.embeddings[,1], y=gene_li[[s]]@reductions$umap@cell.embeddings[,2], logCnt=log10(scRNA_li[[s]]$nCount_RNA + 1))


save(gene_li, file=file.path("analysis", "ATAC_preprocess", "ATAC.gene_li.RData"))

load("data/final_cluster/mergeSCE.RData")

SCE2seurat <- function(mergeSCE, label){
  SCE <- mergeSCE[,mergeSCE$SampleName == label]
  scRNA_seurat <- CreateSeuratObject(counts(SCE), assay = "RNA", project = label)
  scRNA_seurat$Barcode <- SCE$Barcode
  scRNA_seurat$Sample <- SCE$Sample
  scRNA_seurat$SampleName <- SCE$SampleName
  scRNA_seurat$cluster <- SCE$cluster
  scRNA_seurat$SubClu <- SCE$SubClu
  return(scRNA_seurat)
}

scRNA_li <- list()
for(s in sample_li){
  scRNA_li[[s]] <- SCE2seurat(mergeSCE, s)
}

ATAC_RNA_anchor_li <- list()
for(s in sample_li){
  ATAC_RNA_anchor_li[[s]] <- FindTransferAnchors(
    reference = scRNA_li[[s]] , query = gene_li[[s]], features =  VariableFeatures(object = gene_li[[s]]), 
    reduction = "cca")
}

ATAC_RNA_prediction_li <- list()
for(s in sample_li){
  ATAC_RNA_prediction_li[[s]] <- TransferData(
    anchorset = ATAC_RNA_anchor_li[[s]], 
    refdata = scRNA_li[[s]]$SubClu,
    weight.reduction = gene_li[[s]][["umap"]]
    )
}

table(scRNA_li[["WT"]]$SubClu)
table(ATAC_RNA_prediction_li[["WT"]]$predicted.id)
table(scRNA_li[["FL"]]$SubClu)
table(ATAC_RNA_prediction_li[["FL"]]$predicted.id)

simu_rand_scATAC <- function(SCE){
  cnt_per_gene <- colSums(SCE@assays$RNA@counts)
  cnt_per_cell <- SCE$nCount_RNA
  total_reads <- sum(cnt_per_gene)
  mu_gene <- cnt_per_gene / total_reads
  mat <- matrix(0, nrow = length(cnt_per_gene), ncol = length(cnt_per_cell))
  for(cell_indx in 1:length(cnt_per_cell)){
    cell_cnt <- cnt_per_cell[cell_indx]
    cell_simu <- sapply(mu_gene*cell_cnt, function(x){return(rpois(1, x))})
    mat[,cell_indx] <- cell_simu
  }
  colnames(mat) <- names(cnt_per_cell)
  rownames(mat) <- names(cnt_per_gene)
  obj <- CreateSeuratObject(
    mat,
    project = "Simu",
    assay = "ATAC"
  )
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(
    obj, selection.method = "vst",
    nfeatures = 500, verbose = FALSE)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:3, n.neighbors=50, min.dist = 0.5)
  obj <- FindNeighbors(obj, reduction = "umap", dims = 1:2)
  obj <- FindClusters(obj, resolution = 0.01)
  return(obj)
}

simu_rand_scATAC_li <- lapply(gene_li, simu_rand_scATAC)
save(simu_rand_scATAC_li, file=file.path("analysis", "ATAC_preprocess", "ATAC.gene_li.rand.RData"))

p <- DimPlot(simu_rand_scATAC_li[["WT"]], reduction = "pca") + 
  DimPlot(simu_rand_scATAC_li[["FL"]], reduction = "pca") + 
  DimPlot(simu_rand_scATAC_li[["WT"]], reduction = "umap") + 
  DimPlot(simu_rand_scATAC_li[["FL"]], reduction = "umap")
ggsave(
  plot = p, 
  filename = file.path("analysis", "ATAC_preprocess", "ATAC.rand.umap.pdf"),
  width = 10,
  height = 8,
  units = "cm"
)

load_data <- function(f, label){
  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)
  mat <- as.matrix(df[,2:ncol(df)])
  rownames(mat) <- df$Peak
  obj <- CreateSeuratObject(
    mat,
    project = label,
    assay = "ATAC"
  )
  obj$Sample <- label
  return(obj)
}

peak_li <- list()
for (s in sample_li) {
  peak_li[[s]] <- load_data(sprintf("data/scATAC_MACS2_cnt/%s.peak.cnt.tsv", s), s)
  gc()
}

for(s in sample_li){
  peak_li[[s]] <- NormalizeData(peak_li[[s]])
  peak_li[[s]] <- FindVariableFeatures(
    peak_li[[s]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
  VariableFeaturePlot(peak_li[[s]])
  peak_li[[s]] <- ScaleData(peak_li[[s]])
  peak_li[[s]] <- RunPCA(peak_li[[s]], features = VariableFeatures(object = peak_li[[s]]))
}

for(s in sample_li){
  peak_li[[s]] <- JackStraw(peak_li[[s]], num.replicate = 100)
  peak_li[[s]] <- ScoreJackStraw(peak_li[[s]], dims = 1:20)
  JackStrawPlot(peak_li[[s]], dims = 1:15)
  p <- ElbowPlot(peak_li[[s]])
  ggsave(
    plot = p, 
    filename = file.path("analysis", "ATAC_preprocess", sprintf("ATAC.PC.peak.%s.pdf", s)),
    width = 5,
    height = 4,
    units = "cm"
  )
}

for(s in sample_li){
  peak_li[[s]] <- RunUMAP(peak_li[[s]], reduction = "pca", dims = 1:3, n.neighbors=50, min.dist = 0.5)
  peak_li[[s]] <- FindNeighbors(peak_li[[s]], reduction = "umap", dims = 1:2)
  peak_li[[s]] <- FindClusters(peak_li[[s]], resolution = 0.01)
}

for(s in sample_li){
  peak_li[[s]] <- RunTSNE(peak_li[[s]])
}

ggsave(
  plot = DimPlot(peak_li[["WT"]], reduction = "pca", cols=c("red", "red")), 
  filename = file.path("analysis", "ATAC_preprocess", "ATAC.WT.PCA.pdf"),
  width = 10,
  height = 8,
  units = "cm"
)
ggsave(
  plot = DimPlot(peak_li[["FL"]], reduction = "pca", cols=c("blue", "blue")), 
  filename = file.path("analysis", "ATAC_preprocess", "ATAC.FL.PCA.pdf"),
  width = 10,
  height = 8,
  units = "cm"
)
ggsave(
  plot = DimPlot(peak_li[["WT"]], reduction = "tsne", cols=c("red", "red")), 
  filename = file.path("analysis", "ATAC_preprocess", "ATAC.WT.tSNE.pdf"),
  width = 10,
  height = 8,
  units = "cm"
)
ggsave(
  plot = DimPlot(peak_li[["FL"]], reduction = "tsne", cols=c("blue", "blue")), 
  filename = file.path("analysis", "ATAC_preprocess", "ATAC.FL.tSNE.pdf"),
  width = 10,
  height = 8,
  units = "cm"
)
ggsave(
  plot = DimPlot(peak_li[["WT"]], reduction = "umap", cols=c("red", "red")), 
  filename = file.path("analysis", "ATAC_preprocess", "ATAC.WT.UMAP.pdf"),
  width = 10,
  height = 8,
  units = "cm"
)
ggsave(
  plot = DimPlot(peak_li[["FL"]], reduction = "umap", cols=c("blue", "blue")), 
  filename = file.path("analysis", "ATAC_preprocess", "ATAC.FL.UMAP.pdf"),
  width = 10,
  height = 8,
  units = "cm"
)

save(peak_li, file=file.path("analysis", "ATAC_preprocess", "ATAC.peak_li.RData"))
