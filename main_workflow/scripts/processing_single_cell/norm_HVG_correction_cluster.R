library(ggplot2)
library(DropletUtils)
library(scran)
library(batchelor)
library(scater)
library(readr)
library(argparse)
library(dplyr)
library(SC3)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input")
parser$add_argument("-s", "--sample", nargs="*", required=TRUE, type="character", dest = "sample", metavar="sample")
parser$add_argument("-r", "--rep", nargs="*", required=TRUE, type="character", dest = "rep", metavar="rep")
parser$add_argument("-b", "--batch", nargs="*", required=TRUE, type="character", dest = "batch", metavar="batch")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
parser$add_argument("--hvg-fdr", required=FALSE, type="double", dest = "hvg_fdr", metavar="hvg_fdr", default=1e-6)
parser$add_argument("--k-mnn", required=FALSE, type="integer", dest = "k_mnn", metavar="k_mnn", default=20)
parser$add_argument("--d-mnn", required=FALSE, type="integer", dest = "d_mnn", metavar="d_mnn", default=50)
parser$add_argument("--ndist-mnn", required=FALSE, type="integer", dest = "ndist_mnn", metavar="ndist_mnn", default=3)
parser$add_argument("--ntop-umap", required=FALSE, type="integer", dest = "ntop_umap", metavar="ntop_umap", default=500)
parser$add_argument("--nneighbors-umap", required=FALSE, type="integer", dest = "nneighbors_umap", metavar="nneighbors_umap", default=15)
parser$add_argument("--k-sc3", required=FALSE, type="integer", dest = "k_sc3", metavar="k_sc3", default=10)
parser$add_argument("--min-cell-per-cluster", required=FALSE, type="integer", dest = "min_cell_per_cluster", metavar="min_cell_per_cluster", default=100)
parser$add_argument("--min-cell-ratio-cluster", required=FALSE, type="double", dest = "min_cell_ratio_cluster", metavar="min_cell_ratio_cluster", default=0.01)
parser$add_argument("--n-cores", required=FALSE, type="integer", dest = "n_cores", metavar="n_cores", default=8)

args <- commandArgs(TRUE)

# args <- c(
#   "-i",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/WT_rep1/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/WT_rep2/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/FL_rep1/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/FL_rep2/filter_feature_bc_matrix",
#   "-s", "WT", "WT", "FL", "FL",
#   "-r", "rep1", "rep2", "rep1", "rep2",
#   "-o", "~/"
# )
# 
# args <- c(
#   "-i",
#   # "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/WT_rep1/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/WT_rep2/filter_feature_bc_matrix",
#   # "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/FL_rep1/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/FL_rep2/filter_feature_bc_matrix",
#   "-s", "WT", "FL",
#   "-r",  "rep2", "rep2",
#   "--max-hvg-num", "500",
#   "-o", "~/"
# )

args <- parser$parse_args(args) 
file_num <- length(args$input)
if((length(args$sample)!=file_num) | (length(args$rep)!=file_num) | (length(args$batch)!=file_num)) stop()

sample_name <- sprintf("%s.%s", args$sample, args$rep)

print("Loading files ...")
sc_count_expr_li <- list()
for(indx in 1:file_num){
  tmp_name <- sample_name[indx]
  sc_count_expr_li[[tmp_name]] <- read10xCounts(args$input[indx], col.names=TRUE, sample.names=tmp_name)
  sc_count_expr_li[[tmp_name]]$TechBatch <- args$batch[indx]
}

my_multiBatchNorm <- function(...){
  return(multiBatchNorm(..., norm.args=list(use_altexps=FALSE)))
}

print("multiBatchNorm ...")
sc_count_expr_li <- do.call(my_multiBatchNorm, sc_count_expr_li)

print("modelGeneVar ...")
dec_li <- lapply(sc_count_expr_li, modelGeneVar)
combined.dec <- do.call(combineVar, dec_li)

df <- combined.dec[,c("mean", "total", "FDR", "bio")]
df$pass <- df$bio>0 & df$FDR<args$hvg_fdr
df$pass[is.na(df$pass)] <- FALSE
df <- as.data.frame(df)
df$pass <- factor(df$pass, levels = c(TRUE, FALSE), labels = c("True", "False"))

p <- ggplot(df, aes(x=mean, y=total, color=pass)) +
  theme_bw() +
  geom_point(size=0.3, alpha=0.3) +
  scale_color_manual(
    values = c("True"="blue", "False"="grey70"), 
    breaks=c("True", "False"), 
    labels=c(
      sprintf("True: %d", sum(df$pass=="True")),
      sprintf("False: %d", sum(df$pass=="False"))
    )
  ) +
  labs(x="Mean log-expression", y="Variance", color="Chosen HVG") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.15, 0.87),
        panel.grid = element_blank())
ggsave(file.path(args$output, "HVG_plot.png"), p, height = 6, width = 8, units = "cm", dpi = 600)

chosen.hvgs <- getTopHVGs(combined.dec, fdr.threshold=args$hvg_fdr)

countSCE <- sc_count_expr_li

no_correction <- function(...){
  return(correctExperiments(..., subset.row=chosen.hvgs, PARAM=NoCorrectParam()))
}

batch_correction <- function(...){
  return(correctExperiments(..., subset.row=chosen.hvgs, PARAM=FastMnnParam(k=args$k_mnn, ndist=args$ndist_mnn, d=args$d_mnn)))
}

build_cluster <- function(sc_expr, corrected){
  cnt_mat <- counts(sc_expr)
  log_cnt_mat <- logcounts(sc_expr)                                                                        
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(cnt_mat),
      logcounts = as.matrix(log_cnt_mat)
    ),
    colData = colData(sc_expr)
  )
  colnames(sce) <- sprintf("CELL%d", 1:ncol(cnt_mat))
  rowData(sce)$feature_symbol <- rownames(cnt_mat)
  sce <- sc3(sce, ks = args$k_sc3, n_cores=args$n_cores)
  cluster <- colData(sce)[,sprintf("sc3_%d_clusters", args$k_sc3)]
  cell_num <- table(cluster)
  fail_clu <- as.integer(names(cell_num)[(cell_num<args$min_cell_per_cluster) | ((cell_num/sum(cell_num))<args$min_cell_ratio_cluster)])
  # dbl <- doubletCluster(sc_expr, cluster)
  # doublet_indx <- isOutlier(dbl$N, log=TRUE, type="lower") | (dbl$lib.size1 < 1 & dbl$lib.size2 < 1 )
  # doublet_clu <- as.integer(rownames(dbl)[doublet_indx])
  cluster[cluster %in% fail_clu] <- NA
  # cluster[cluster %in% doublet_clu] <- NA
  clu_num <- as.data.frame(table(cluster))
  clu_num <- clu_num[order(clu_num$Freq, decreasing = TRUE),]
  clu_num$label <- sprintf("C%d", seq(1: nrow(clu_num)))
  cluster <- factor(cluster, levels = clu_num$cluster, labels = clu_num$label)
  return(cluster)
}

print("No correction ...")
no_correction_sc_expr <- do.call(no_correction, sc_count_expr_li)
no_correction_cluster <- build_cluster(no_correction_sc_expr, FALSE)
no_correction_sc_expr$cluster <- no_correction_cluster
no_correction_sc_expr <- runUMAP(no_correction_sc_expr, subset_row=chosen.hvgs, ntop=args$ntop_umap, n_neighbors=args$nneighbors_umap)

inche_cm=2.54
pdf(file.path(args$output, "NoCorrected.pdf"), width=10/inche_cm, height=8/inche_cm)
plotUMAP(no_correction_sc_expr, colour_by="batch", point_alpha=0.1, point_size=0.1)
plotUMAP(no_correction_sc_expr, colour_by="TechBatch", point_alpha=0.1, point_size=0.1)
plotUMAP(no_correction_sc_expr, colour_by="cluster", point_size=0.1)
dev.off()

print("Correction ...")
correction_sc_expr <- do.call(batch_correction, sc_count_expr_li)
correction_cluster <- build_cluster(correction_sc_expr, TRUE)
correction_sc_expr$cluster <- correction_cluster
correction_sc_expr <- runUMAP(correction_sc_expr, subset_row=chosen.hvgs, ntop=args$ntop_umap, n_neighbors=args$nneighbors_umap, dimred="corrected")

pdf(file.path(args$output, "Corrected.pdf"), width=10/inche_cm, height=8/inche_cm)
plotUMAP(correction_sc_expr, colour_by="batch", point_alpha=0.1, point_size=0.1)
plotUMAP(correction_sc_expr, colour_by="TechBatch", point_alpha=0.1, point_size=0.1)
plotUMAP(correction_sc_expr, colour_by="cluster", point_size=0.1)
dev.off()

umapSCE <- correction_sc_expr[,!is.na(correction_sc_expr$cluster)]
umapSCE$SampleName <- sapply(strsplit(as.character(umapSCE$batch), "[.]"), function(x) x[[1]])
umapSCE$Rep <- sapply(strsplit(as.character(umapSCE$batch), "[.]"), function(x) x[[2]])
save(umapSCE, file=file.path(args$output, "umapSCE.RData"))

df_cluster <- data.frame(Sample=umapSCE$SampleName, Rep=umapSCE$Rep, Cluster=umapSCE$cluster, Barcode=umapSCE$Barcode)
df_cluster <- df_cluster[order(df_cluster$Sample, df_cluster$Rep, df_cluster$Cluster, df_cluster$Barcode),]
write_tsv(df_cluster, file.path(args$output, "BarcodeCluster.tsv"))

pdf(file.path(args$output, "Filtered.pdf"), width=10/inche_cm, height=8/inche_cm)
plotUMAP(umapSCE, colour_by="batch", point_alpha=0.1, point_size=0.1)
plotUMAP(umapSCE, colour_by="TechBatch", point_alpha=0.1, point_size=0.1)
plotUMAP(umapSCE, colour_by="SampleName", point_alpha=0.1, point_size=0.1)
plotUMAP(umapSCE, colour_by="Rep", point_alpha=0.1, point_size=0.1)
plotUMAP(umapSCE, colour_by="cluster", point_size=0.1)
dev.off()

sample_li <- unique(args$sample)
for(tmp_sample in sample_li){
  sample_umapSCE <- umapSCE[,umapSCE$SampleName==tmp_sample]
  save(umapSCE, file=file.path(args$output, sprintf("umapSCE.%s.RData", tmp_sample)))
  pdf(file.path(args$output, sprintf("Filtered.%s.pdf", tmp_sample)), width=10/inche_cm, height=8/inche_cm)
  print(plotTSNE(sample_umapSCE, colour_by="batch", point_alpha=0.1))
  print(plotTSNE(sample_umapSCE, colour_by="cluster"))
  print(plotUMAP(sample_umapSCE, colour_by="batch", point_alpha=0.1))
  print(plotUMAP(sample_umapSCE, colour_by="cluster"))
  dev.off()
}

mergeSCE <- do.call(cbind, countSCE)
mergeSCE$SampleName <- sapply(strsplit(as.character(mergeSCE$Sample), "[.]"), function(x) x[[1]])
mergeSCE$Rep <- sapply(strsplit(as.character(mergeSCE$Sample), "[.]"), function(x) x[[2]])
cluster_uniq <- sprintf("%s.%s.%s", df_cluster$Sample, df_cluster$Rep, df_cluster$Barcode)
sce_uniq <- sprintf("%s.%s.%s", mergeSCE$SampleName, mergeSCE$Rep, mergeSCE$Barcode)
mergeSCE <- mergeSCE[,sce_uniq %in% cluster_uniq]

merge_barcode_df <- data.frame(Barcode=mergeSCE$Barcode, Sample=mergeSCE$SampleName, Rep=mergeSCE$Rep)
merge_barcode_df <- left_join(merge_barcode_df, df_cluster, by=c("Barcode", "Sample", "Rep"))
mergeSCE$cluster <- merge_barcode_df$Cluster
save(mergeSCE, file=file.path(args$output, "mergeSCE.RData"))
