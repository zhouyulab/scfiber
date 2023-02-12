library(BiocGenerics)
library(readr)
library(dplyr)
library(RColorBrewer)
library(argparse)
library(SingleCellExperiment)
library(velocyto.R)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input")
parser$add_argument("-s", "--sample", nargs="*", required=TRUE, type="character", dest = "sample", metavar="sample")
parser$add_argument("-r", "--rep", nargs="*", required=TRUE, type="character", dest = "rep", metavar="rep")
parser$add_argument("-c", "--cluster", required=TRUE, type="character", dest = "cluster", metavar="cluster.tsv")
parser$add_argument("--corSCE", required=TRUE, type="character", dest = "corSCE", metavar="corSCE.RData")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
parser$add_argument("--deltaT", required=FALSE, type="double", dest = "deltaT", metavar="deltaT", default=1)
parser$add_argument("--kCells", required=FALSE, type="integer", dest = "kCells", metavar="kCells", default=25)
parser$add_argument("--fit-quantile", required=FALSE, type="double", dest = "fitQuantile", metavar="fitQuantile", default=0.02)
parser$add_argument("--plotN", required=FALSE, type="integer", dest = "plotN", metavar="plotN", default=100)


args <- commandArgs(TRUE)
args <- parser$parse_args(args) 

file_num <- length(args$input)
if((length(args$sample)!=file_num) | (length(args$rep)!=file_num)) stop()

load(args$corSCE)

cluster <- read_delim(args$cluster, "\t", escape_double = FALSE, trim_ws = TRUE)
cluster$LoomName <- sprintf("%s_%s:%sx", cluster$Sample, cluster$Rep, sapply(strsplit(cluster$Barcode, "-"), function(x) x[[1]]) )

looms <- list()
for(indx in 1:file_num){
  loom_name <- sprintf("%s_%s", args$sample[indx], args$rep[indx])
  looms[[loom_name]] <- read.loom.matrices(args$input[indx])
}

corSCE$LoomName <- sprintf("%s_%s:%sx", corSCE$SampleName, corSCE$Rep, sapply(strsplit(colnames(corSCE), "-"), function(x) x[[1]]))

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

velocity_estimation <- function(loom, corSCE, deltaT=args$deltaT, kCells=args$kCells, fit.quantile=args$fitQuantile){
  emat <- loom$spliced
  nmat <- loom$unspliced
  cluster <- loom$Cluster
  umapSCE <- corSCE[,corSCE$LoomName %in% colnames(emat)]
  loom_order <- order(colnames(emat))
  emat <- emat[,loom_order]
  nmat <- nmat[,loom_order]
  cluster <- cluster[loom_order]
  SCE_order <- order(umapSCE$LoomName)
  umapSCE <- umapSCE[,SCE_order]
  if(! all(umapSCE$LoomName == colnames(emat))) stop()
  
  cell.dist <- as.dist(1-armaCor(t(umapSCE@int_colData@listData$reducedDims$UMAP)))
  cell.dist <- as.matrix(cell.dist)
  colnames(cell.dist) <- colnames(emat)
  rownames(cell.dist) <- colnames(emat)
  cell.dist <- as.dist(cell.dist)

  names(cluster) <- colnames(emat)
  
  emat <- filter.genes.by.cluster.expression(emat, cluster, min.max.cluster.average = 0.2)
  nmat <- filter.genes.by.cluster.expression(nmat, cluster, min.max.cluster.average = 0.05)
  
  rvel.cd <- gene.relative.velocity.estimates(emat,nmat,n.cores=20,deltaT=deltaT,kCells=kCells,cell.dist=cell.dist,fit.quantile=fit.quantile)
  rvel.cd$UMAP <- umapSCE@int_colData@listData$reducedDims$UMAP
  rownames(rvel.cd$UMAP) <- colnames(emat)
  rvel.cd$Cluster <- cluster
  return(rvel.cd)
}

# velocity_li <- lapply(filtered_looms, velocity_estimation, corSCE=corSCE)


# inche_cm=2.54
# pdf(file.path(args$output, "Sample.velocity.UMAP.pdf"), width=20/inche_cm, height=20/inche_cm, family="ArialMT", colormodel = "cmyk")
# for(n in names(velocity_li)){
#   rvel.cd <- velocity_li[[n]]
#   col_li <- factor(rvel.cd$Cluster, levels = c("C1", "C2", "C3", "C4", "C5"), labels = brewer.pal(5, "Dark2"))
# 
#   show.velocity.on.embedding.cor(
#     rvel.cd$UMAP, rvel.cd ,
#     cell.colors=ac(col_li,alpha=0.5),
#     n=args$plotN, scale='sqrt', cex=0.8,arrow.scale=1,
#     show.grid.flow=TRUE,min.grid.cell.mass=0.5,
#     grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1,
#     xlab = "UMAP1", ylab = "UMAP2", main=n)
# }
# dev.off()

merge_loom <- function(looms_li){
  loom_names <- names(looms_li)
  for(indx in 1:length(loom_names)){
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
merged_velocity <- velocity_estimation(merged_loom, corSCE=corSCE)

inche_cm=2.54
pdf(file.path(args$output, "Merged.velocity.UMAP.pdf"), width=20/inche_cm, height=20/inche_cm, family="ArialMT", colormodel = "cmyk")
col_li <- factor(merged_velocity$Cluster, levels = c("C1", "C2", "C3", "C4", "C5"), labels = brewer.pal(5, "Dark2"))
  
show.velocity.on.embedding.cor(
  merged_velocity$UMAP, merged_velocity, 
  cell.colors=ac(col_li,alpha=0.5), n.cores=20,
  n=args$plotN, scale='sqrt', cex=0.8, arrow.scale=1,
  show.grid.flow=TRUE, min.grid.cell.mass=0.5,
  grid.n=40, arrow.lwd=1, do.par=F, cell.border.alpha = 0.1,
  xlab = "UMAP1", ylab = "UMAP2", main="Merged")

dev.off()
