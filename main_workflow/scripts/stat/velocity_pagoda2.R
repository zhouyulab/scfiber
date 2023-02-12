library(BiocGenerics)
library(readr)
library(dplyr)
library(RColorBrewer)
library(argparse)
library(SingleCellExperiment)
library(velocyto.R)
library(pagoda2)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input")
parser$add_argument("-s", "--sample", nargs="*", required=TRUE, type="character", dest = "sample", metavar="sample")
parser$add_argument("-r", "--rep", nargs="*", required=TRUE, type="character", dest = "rep", metavar="rep")
parser$add_argument("-c", "--cluster", required=TRUE, type="character", dest = "cluster", metavar="cluster.tsv")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")

parser$add_argument("--nPcs", required=FALSE, type="integer", dest = "nPcs", metavar="nPcs", default=100)
parser$add_argument("--k_KNN", required=FALSE, type="integer", dest = "k_KNN", metavar="k_KNN", default=30)
parser$add_argument("--gam-k", required=FALSE, type="integer", dest = "gam.k", metavar="gam.k", default=10)
parser$add_argument("--n-hvg", required=FALSE, type="integer", dest = "n_hvg", metavar="n_hvg", default=3000)

parser$add_argument("--deltaT", required=FALSE, type="double", dest = "deltaT", metavar="deltaT", default=1)
parser$add_argument("--kCells", required=FALSE, type="integer", dest = "kCells", metavar="kCells", default=25)
parser$add_argument("--fit-quantile", required=FALSE, type="double", dest = "fitQuantile", metavar="fitQuantile", default=0.02)
parser$add_argument("--plotN", required=FALSE, type="integer", dest = "plotN", metavar="plotN", default=100)


args <- commandArgs(TRUE)

setwd("~/ias/data/CottonSingleCell")
args <- strsplit("-i analysis/cellranger/count/WT_rep2/velocyto/WT_rep2.loom analysis/cellranger/count/FL_rep2/velocyto/FL_rep2.loom -s WT FL -r rep2 rep2 --cluster analysis/base_map/cluster.tsv -o /home/sqreb", " ")[[1]]

args <- parser$parse_args(args) 

file_num <- length(args$input)
if((length(args$sample)!=file_num) | (length(args$rep)!=file_num)) stop()

cluster <- read_delim(args$cluster, "\t", escape_double = FALSE, trim_ws = TRUE)
cluster$LoomName <- sprintf("%s_%s:%sx", cluster$Sample, cluster$Rep, sapply(strsplit(cluster$Barcode, "-"), function(x) x[[1]]) )

looms <- list()
for(indx in 1:file_num){
  loom_name <- sprintf("%s_%s", args$sample[indx], args$rep[indx])
  looms[[loom_name]] <- read.loom.matrices(args$input[indx])
}

filter_loom <- function(looms, cluster){
  loom_names <- sort(names(looms))
  new_looms <- list()
  for(indx in 1:length(loom_names)){
    name <- loom_names[indx]
    tmp_spliced_mat <- looms[[name]]$spliced
    tmp_unspliced_mat <- looms[[name]]$unspliced
    tmp_ambiguous_mat <- looms[[name]]$ambiguous
    tmp_barcode_name <- colnames(tmp_spliced_mat)
    selected_indx <- tmp_barcode_name %in% cluster$LoomName
    
    tmp_spliced_mat <- tmp_spliced_mat[, selected_indx]
    tmp_unspliced_mat <- tmp_unspliced_mat[, selected_indx]
    tmp_ambiguous_mat <- tmp_ambiguous_mat[, selected_indx]
    
    tmp_barcode_name <- data.frame(LoomName=colnames(tmp_spliced_mat))
    tmp_info <- left_join(tmp_barcode_name, cluster, by="LoomName")
    tmp_loom_obj <- list(
      "spliced"=tmp_spliced_mat, "unspliced"=tmp_unspliced_mat, "ambiguous"=tmp_ambiguous_mat,
      "Sample"=tmp_info$Sample, "Rep"=tmp_info$Rep, "Cluster"=tmp_info$Cluster, "Barcode"=tmp_info$Barcode
    )
    new_looms[[name]] <- tmp_loom_obj
  }
  return(new_looms)
}

filtered_looms <- filter_loom(looms, cluster)
rm(looms)
gc()

velocity_estimation <- function(
  loom, gam.k=args$gam.k, nPcs=args$nPcs, n_hvg=args$n_hvg, k_KNN=args$k_KNN, 
  deltaT=args$deltaT, kCells=args$kCells, fit.quantile=args$fitQuantile
  ){
  emat <- loom$spliced
  nmat <- loom$unspliced
  cluster <- loom$Cluster
  
  r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
  r$adjustVariance(plot=T,do.par=T,gam.k=gam.k)
  r$calculatePcaReduction(nPcs=nPcs,n.odgenes=n_hvg,maxit=300)
  r$makeKnnGraph(k=k_KNN,type='PCA',center=T,distance='cosine')
  r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
  r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
  
  select_indx <- colnames(emat) %in% rownames(r$counts)
  emat <- emat[,select_indx]
  nmat <- nmat[,select_indx]
  cluster <- cluster[select_indx]
  
  cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))
  names(cluster) <- colnames(emat)
  
  emat <- filter.genes.by.cluster.expression(emat, cluster, min.max.cluster.average = 0.2)
  nmat <- filter.genes.by.cluster.expression(nmat, cluster, min.max.cluster.average = 0.05)
  
  rvel.cd <- gene.relative.velocity.estimates(emat,nmat,n.cores=20,deltaT=deltaT,kCells=kCells,cell.dist=cell.dist,fit.quantile=fit.quantile)
  rvel.cd$tSNE <- r$embeddings$PCA$tSNE

  rvel.cd$Cluster <- cluster
  return(rvel.cd)
}

velocity_li <- lapply(filtered_looms, velocity_estimation)


inche_cm=2.54
pdf(file.path(args$output, "Sample.velocity.tSNE.pdf"), width=20/inche_cm, height=20/inche_cm, family="ArialMT", colormodel = "cmyk")
for(n in names(velocity_li)){
  rvel.cd <- velocity_li[[n]]
  col_li <- factor(rvel.cd$Cluster, levels = c("C1", "C2", "C3", "C4", "C5"), labels = brewer.pal(5, "Dark2"))
  
  show.velocity.on.embedding.cor(
    rvel.cd$tSNE, rvel.cd , 
    cell.colors=ac(col_li,alpha=0.5), 
    n=args$plotN, scale='sqrt', cex=0.8,arrow.scale=2,
    show.grid.flow=TRUE,min.grid.cell.mass=0.5,
    grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1,
    xlab = "tSNE1", ylab = "tSNE2", main=n)
}
dev.off()

merge_loom <- function(looms_li){
  loom_names <- names(looms_li)
  for(indx in 1:length(sample_names)){
    loom_name <- loom_names[indx]
    loom <- looms_li[[loom_name]]
    if(indx == 1){
      spliced <- loom$spliced
      unspliced <- loom$unspliced
      Cluster <- loom$Cluster
    }else{
      spliced <- cbind(spliced, loom$spliced)
      unspliced <- cbind(unspliced, loom$unspliced)
      Cluster <- c(Cluster, loom$Cluster)
    }
  }
  return(list("spliced"=spliced, "unspliced"=unspliced, "Cluster"=Cluster))
}

merged_loom <- merge_loom(filtered_looms)
merged_velocity <- velocity_estimation(merged_loom)

inche_cm=2.54
pdf(file.path(args$output, "Merged.velocity.tSNE.pdf"), width=20/inche_cm, height=20/inche_cm, family="ArialMT", colormodel = "cmyk")
col_li <- factor(merged_velocity$Cluster, levels = c("C1", "C2", "C3", "C4", "C5"), labels = brewer.pal(5, "Dark2"))

show.velocity.on.embedding.cor(
  merged_velocity$UMAP, merged_velocity, 
  cell.colors=ac(col_li,alpha=0.5), n.cores=20,
  n=args$plotN, scale='sqrt', cex=0.8, arrow.scale=2,
  show.grid.flow=TRUE, min.grid.cell.mass=0.5,
  grid.n=40, arrow.lwd=1, do.par=F, cell.border.alpha = 0.1,
  xlab = "tSNE1", ylab = "tSNE2", main="Merged")

dev.off()
