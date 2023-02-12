library(ggplot2)
library(DropletUtils)
library(scran)
library(batchelor)
library(scater)
library(readr)
library(argparse)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(gridExtra)
library(BiocParallel)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input")
parser$add_argument("-s", "--sample", nargs="*", required=TRUE, type="character", dest = "sample", metavar="sample")
parser$add_argument("-r", "--rep", nargs="*", required=TRUE, type="character", dest = "rep", metavar="rep")
parser$add_argument("-t", "--tech", nargs="*", required=TRUE, type="character", dest = "tech", metavar="tech")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
parser$add_argument("--hvg-fdr", required=FALSE, type="double", dest = "hvg_fdr", metavar="hvg_fdr", default=0.05)
parser$add_argument("--k-mnn", required=FALSE, type="integer", dest = "k_mnn", metavar="k_mnn", default=25)
parser$add_argument("--sigma-mnn", required=FALSE, type="double", dest = "sigma_mnn", metavar="sigma_mnn", default=0.3)
args <- commandArgs(TRUE)

# args <- c(
#   "-i",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/WT_rep2/filter_feature_bc_matrix",
#   # "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/WT_rep3/filter_feature_bc_matrix",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/FL_rep2/filter_feature_bc_matrix",
#   # "/home/sqreb/ias/data/CottonSingleCell/analysis/filter/filter_count/FL_rep3/filter_feature_bc_matrix",
#   "-s", "WT", "FL", 
#   "-r", "rep2",  "rep2", 
#   "-t", "v2_lib", "v2_lib", 
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
  sc_count_expr_li[[tmp_name]] <- read10xCounts(args$input[indx], col.names=TRUE, sample.names=tmp_name)
  sc_count_expr_li[[tmp_name]]$Tech <- args$tech[indx]
  sc_count_expr_li[[tmp_name]]$SampleName <- sapply(strsplit(as.character(sc_count_expr_li[[tmp_name]]$Sample), "[.]"), function(x) x[[1]])
  sc_count_expr_li[[tmp_name]]$Rep <- sapply(strsplit(as.character(sc_count_expr_li[[tmp_name]]$Sample), "[.]"), function(x) x[[2]])
  # select_indx <- sample(1:ncol(sc_count_expr_li[[tmp_name]]), 2000)
  # sc_count_expr_li[[tmp_name]] <- sc_count_expr_li[[tmp_name]][,select_indx]
}

my_multiBatchNorm1 <- function(...){
  return(multiBatchNorm(..., normalize.all=TRUE, norm.args = list(log=FALSE)))
}
my_multiBatchNorm2 <- function(...){
  return(multiBatchNorm(..., normalize.all=TRUE))
}

print("multiBatchNorm ...")
sc_count_expr_li <- do.call(my_multiBatchNorm1, sc_count_expr_li)
sc_count_expr_li <- do.call(my_multiBatchNorm2, sc_count_expr_li)

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
save(chosen.hvgs, file=file.path(args$output, "chosen.hvgs.RData"))

countSCE <- sc_count_expr_li
save(countSCE, file=file.path(args$output, "countSCE.RData"))
rm(countSCE)
gc()

no_correction <- function(...){
  return(correctExperiments(..., subset.row = chosen.hvgs, PARAM=NoCorrectParam()))
}
no_correction2 <- function(...){
  return(correctExperiments(..., PARAM=NoCorrectParam()))
}

mergeSCE <- do.call(no_correction2, sc_count_expr_li)
save(mergeSCE, file=file.path(args$output, "mergeSCE.RData"))

cnt_mat <- normcounts(mergeSCE)
expr_cell_num <- rowSums(cnt_mat>2)
expr_cell_ratio <- expr_cell_num / ncol(cnt_mat)
filter_cnt_mat <- cnt_mat[(expr_cell_num>100) & (expr_cell_ratio<0.5) & (expr_cell_ratio>0.01),]
filter_cnt_mat <- as.matrix(filter_cnt_mat)
rm(cnt_mat)
gc()


compute_modularity <- function(gid, filter_cnt_mat, expr_cutoff, SNN_graph){
  library(scran)
  is_expr <- filter_cnt_mat[which(rownames(filter_cnt_mat)==gid),] > expr_cutoff
  modularity_ratio <- clusterModularity(SNN_graph, is_expr, as.ratio=TRUE)
  df <- data.frame(Gid=gid, ModularityRatio=modularity_ratio[2,2], ExprCellNum=sum(is_expr), ExprCellRatio=mean(is_expr))
  return(df)
}

SNN_graph <- buildSNNGraph(mergeSCE[names(mergeSCE) %in% chosen.hvgs,], assay.type="logcounts", k=args$k_mnn, type="rank")

print("Compute modularity")
res <- bplapply(rownames(filter_cnt_mat), compute_modularity, filter_cnt_mat=filter_cnt_mat, expr_cutoff=2, SNN_graph=SNN_graph, BPPARAM=MulticoreParam(8))

res_df <- plyr::ldply(res, function(x) x)
write_tsv(res_df, file.path(args$output, "modularity.tsv"))
modularity_enrich_gene <- as.character(res_df$Gid[res_df$ModularityRatio>10])

chosen.hvgs <- unique(c(chosen.hvgs, modularity_enrich_gene))

rm(mergeSCE)
gc()

batch_correction <- function(...){
  return(correctExperiments(..., subset.row = chosen.hvgs, PARAM=ClassicMnnParam(k=args$k_mnn, sigma=args$sigma_mnn)))
}

plot_umap <- function(umapSCE, color_by, color_label, dim1=1, dim2=2, label=NULL, col_val=NULL, alpha=0.3){
  df <- data.frame(x=umapSCE@int_colData$reducedDims$UMAP[,dim1], y=umapSCE@int_colData$reducedDims$UMAP[,dim2], Type=umapSCE[[color_by]])
  p <- ggplot(df, aes(x=x, y=y, color=Type)) +
    theme_bw() +
    geom_point(size=0.1, alpha=alpha) +
    labs(x=sprintf("UMAP%s", dim1), y=sprintf("UMAP%s", dim2), color=color_label) +
    # scale_color_brewer(palette="Set2") +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  if(!is.null(label)) p <- p + labs(title=label)
  if(!is.null(col_val)) p <- p + scale_color_manual(values = col_val)
  return(p)
}

plot_umap_panal <- function(umapSCE, dim1, dim2, col_li){
  ps <- list()
  ps[["Batch"]] <- plot_umap(umapSCE, "batch", "Batch", col_val = col_li[["Batch"]], dim1=dim1, dim2=dim2)
  ps[["Tech"]] <- plot_umap(umapSCE, "Tech", "Tech", col_val = col_li[["Tech"]], dim1=dim1, dim2=dim2)
  ps[["Sample"]] <- plot_umap(umapSCE, "SampleName", "Sample", col_val = col_li[["Sample"]], dim1=dim1, dim2=dim2)
  ps[["Rep"]] <- plot_umap(umapSCE, "Rep", "Rep", col_val = col_li[["Rep"]], dim1=dim1, dim2=dim2)
  return(ps)
}

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

print("No correction ...")
no_correction_sc_expr <- do.call(no_correction, sc_count_expr_li)
no_correction_sc_expr <- runUMAP(no_correction_sc_expr, subset_row=chosen.hvgs, ncomponents=4)
noCorSCE <- no_correction_sc_expr
save(noCorSCE, file=file.path(args$output, "noCorSCE.RData"))
rm(noCorSCE)
gc()

inche_cm=2.54
pdf(file.path(args$output, "NoCorrected.pdf"), width=20/inche_cm, height=16/inche_cm)
grid.arrange(grobs=plot_umap_panal(no_correction_sc_expr, 1, 2, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(no_correction_sc_expr, 1, 3, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(no_correction_sc_expr, 1, 4, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(no_correction_sc_expr, 2, 3, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(no_correction_sc_expr, 2, 4, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(no_correction_sc_expr, 3, 4, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))

dev.off()

print("Correction ...")
correction_sc_expr <- do.call(batch_correction, sc_count_expr_li)
correction_sc_expr <- runUMAP(correction_sc_expr, subset_row=chosen.hvgs, exprs_values="corrected", ncomponents=4)
corSCE <- correction_sc_expr
save(corSCE, file=file.path(args$output, "corSCE.RData"))
rm(corSCE)
gc()

pdf(file.path(args$output, "Corrected.pdf"), width=20/inche_cm, height=16/inche_cm)
grid.arrange(grobs=plot_umap_panal(correction_sc_expr, 1, 2, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(correction_sc_expr, 1, 3, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(correction_sc_expr, 1, 4, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(correction_sc_expr, 2, 3, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(correction_sc_expr, 2, 4, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
grid.arrange(grobs=plot_umap_panal(correction_sc_expr, 3, 4, col_li), ncol=2, nrow=2, widths=c(1,1), heights=c(1,1), padding = unit(1, "mm"))
dev.off()
