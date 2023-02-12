library(readr)
library(argparse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(BiocParallel)

parser <- ArgumentParser()
parser$add_argument("-m", "--motif", nargs="*", required=TRUE, type="character", dest = "motif", metavar="motif")
parser$add_argument("--motif-base", required=TRUE, type="character", dest = "motif_base", metavar="motif_base")
parser$add_argument("--ATAC-peak-expr", required=TRUE, type="character", dest = "ATAC_peak_expr", metavar="ATAC_peak_expr.tsv")
parser$add_argument("--ATAC-peak-DE", required=TRUE, type="character", dest = "ATAC_peak_DE", metavar="ATAC_peak_DE.tsv")
parser$add_argument("--ATAC-gene-expr", required=TRUE, type="character", dest = "ATAC_gene_expr", metavar="ATAC_gene_expr.tsv")
parser$add_argument("--ATAC-gene-DE", required=TRUE, type="character", dest = "ATAC_gene_DE", metavar="ATAC_gene_DE.tsv")
parser$add_argument("--scRNA-expr", required=TRUE, type="character", dest = "scRNA_expr", metavar="scRNA_expr.tsv")
parser$add_argument("--RNA-TP-expr", required=TRUE, type="character", dest = "RNA_TP_expr", metavar="RNA_TP_expr.tsv")
parser$add_argument("--RNA-TP-DE", required=TRUE, type="character", dest = "RNA_TP_DE", metavar="RNA_TP_DE.tsv")
parser$add_argument("--RNA-TP-circadian", required=TRUE, type="character", dest = "RNA_TP_circadian", metavar="RNA_TP_circadian.tsv")
parser$add_argument("--tf-ann", required=TRUE, type="character", dest = "tf_ann", metavar="tf_ann.tsv")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
args <- commandArgs(TRUE)
args <- parser$parse_args(args) 

args <- c(
  "-m", "TCP", "TCP_like", "Gbox", "MYB_WT", "MYB_G4", "AP2", "WRKY",
  "--motif-base", "analysis_v3/ATAC/motif",
  "--ATAC-peak-expr", "analysis_v3/stat/ATAC_MACS2/ATAC.peak.tsv",
  "--ATAC-peak-DE", "analysis_v3/stat/ATAC_MACS2/ATAC.peak.DE.tsv",
  "--ATAC-gene-expr", "analysis_v3/stat/ATAC_MACS2/ATAC.gene.tsv",
  "--ATAC-gene-DE", "analysis_v3/stat/ATAC_MACS2/ATAC.gene.DE.tsv",
  "--scRNA-expr", "analysis_v3/stat/cluster_cnt/cluster.tpm.tsv",
  "--RNA-TP-expr", "../CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv",
  "--RNA-TP-DE",  "../CottonSCE_TP/analysis/RNA_seq_TP/DE_same_time/DEGs.type.tsv",
  "--RNA-TP-circadian", "~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CircadianGene.Group.tsv",
  "--tf-ann", "analysis_v3/annotation/TF/TF.gene_list.tsv",
  "-o", "analysis_v3/ATAC/motif/summary"
)
args <- parser$parse_args(args)

ATAC_peak_expr <- read_delim(args$ATAC_peak_expr, "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_gene_expr <- read_delim(args$ATAC_gene_expr, "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_peak_DE <- read_delim(args$ATAC_peak_DE, "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_gene_DE <- read_delim(args$ATAC_gene_DE, "\t", escape_double = FALSE, trim_ws = TRUE)

scRNA_expr <- read_delim(args$scRNA_expr, "\t", escape_double = FALSE, trim_ws = TRUE)
RNA_TP_expr <- read_delim(args$RNA_TP_expr, "\t", escape_double = FALSE, trim_ws = TRUE)
RNA_TP_DE <- read_delim(args$RNA_TP_DE, "\t", escape_double = FALSE, trim_ws = TRUE)
RNA_TP_circadian <- read_delim(args$RNA_TP_circadian, "\t", escape_double = FALSE, trim_ws = TRUE)

load_motif_target_peak <- function(motif, motif_base, ATAC_peak_expr){
  f_motif <- file.path(motif_base, sprintf("%s_motif", motif), sprintf("%s_motif.ATAC.bed", motif))
  motif_bed <- read_delim(f_motif, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
  motif_bed <- motif_bed[,1:3]
  names(motif_bed) <- c("Chrom", "Start", "End")
  res <- ATAC_peak_expr[,c("Chrom", "Start", "End")]
  res$HasMotif <- unlist(bplapply(1:nrow(res), function(indx){
    peak_chrom <- res$Chrom[indx]
    peak_start <- res$Start[indx]
    peak_end <- res$End[indx]
    cond1 <- motif_bed$Chrom == peak_chrom
    cond2 <- peak_start <= motif_bed$Start
    cond3 <- peak_end >= motif_bed$End
    return(any(cond1 & cond2 & cond3))
  }, BPPARAM = MulticoreParam(40)))
  names(res)[4] <- motif
  return(res)
}

motif_target_li <- lapply(args$motif, load_motif_target_peak, motif_base=args$motif_base, ATAC_peak_expr=ATAC_peak_expr)
motif_target_peak_df <- ATAC_peak_expr[,c("Chrom", "Start", "End")]
for(i in 1:length(motif_target_li)){
  motif_target_peak_df <- left_join(motif_target_peak_df, motif_target_li[[i]])
}
write_tsv(motif_target_peak_df, file.path(args$output, "motif_target_peak.tsv"))

load_motif_target_gene <- function(motif, motif_base, gene_expr){
  f_motif <- file.path(motif_base, sprintf("%s_motif", motif), sprintf("%s_motif.promoter.ATAC.bed", motif))
  motif_bed <- read_delim(f_motif, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
  motif_gid <- unlist(motif_bed[,4])
  
  res <- data.frame(Gid=gene_expr$Gid, HasMotif=gene_expr$Gid %in% motif_gid)
  names(res)[2] <- motif
  return(res)
}
motif_target_gene_li <- lapply(args$motif, load_motif_target_gene, motif_base=args$motif_base, gene_expr=scRNA_expr)
motif_target_gene_df <- data.frame(Gid=scRNA_expr$Gid)
for(i in 1:length(motif_target_gene_li)){
  motif_target_gene_df <- left_join(motif_target_gene_df, motif_target_gene_li[[i]])
}
write_tsv(motif_target_gene_df, file.path(args$output, "motif_target_gene.tsv"))

do.test <- function(bool_li1, bool_li2, label1, label2){
  assertthat::assert_that(length(bool_li1) == length(bool_li2))
  total_num <- length(bool_li1)
  feature1_num <- sum(bool_li1 == TRUE)
  feature2_num <- sum(bool_li2 == TRUE)
  overlap_num <- sum(bool_li1 & bool_li2)
  exp_overlap_num <- total_num * (feature1_num/total_num) * (feature2_num/total_num)
  num_mat <- matrix(0, nrow = 2, ncol = 2)
  num_mat[1,1] <- overlap_num
  num_mat[2,1] <- feature1_num - overlap_num
  num_mat[1,2] <- feature2_num - overlap_num
  num_mat[2,2] <- total_num - feature1_num - feature2_num + overlap_num
  pval <- fisher.test(num_mat)$p.value
  res <- data.frame(Set1=label1, Set2=label2, Set1Num=feature1_num, Set2Num=feature2_num, Set1Ratio=feature1_num/total_num, Set2Ratio=feature2_num/total_num, OverlapNum=overlap_num, ExpOverlapNum=exp_overlap_num, OverlapRatio=overlap_num/total_num, Log2EnrichRatio=log2(overlap_num/exp_overlap_num), pval=pval)
  return(res)
}

plot_test_res <- function(res_df, prefix, pval_cutoff=1e-2){
  res_df$x_indx <- as.numeric(res_df$Set1)
  res_df$y_indx <- as.numeric(res_df$Set2)
  x_feature_num <- length(unique(res_df$Set1))
  y_feature_num <- length(unique(res_df$Set2))
  p <- ggplot(res_df)+
    geom_rect(mapping = aes(xmin=x_indx-0.5, xmax=x_indx+0.5, ymin=y_indx-0.5, ymax=y_indx+0.5, fill=log10(OverlapNum))) +
    geom_text(mapping = aes(x=x_indx, y=y_indx, label=OverlapNum), color="black", size=1.2) +
    scale_fill_continuous(low="#FFFF11", high="red") +
    scale_x_continuous(breaks = 1:length(levels(res_df$Set1)), labels = levels(res_df$Set1)) +
    scale_y_continuous(breaks = 1:length(levels(res_df$Set2)), labels = levels(res_df$Set2)) +
    theme_bw() +
    theme(
      text = element_text(family="ArialMT", color = "black", size = 5),
      axis.title = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(file.path(prefix, "overlap.num.ht.pdf"), p, width = 1+0.6*x_feature_num, height = 1+0.6*y_feature_num, limitsize = FALSE, units = "cm")
  res_df$Log2EnrichRatio[res_df$pval>pval_cutoff] <- NA
  res_df$Log2EnrichRatio[res_df$Log2EnrichRatio>1] <- 1
  res_df$Log2EnrichRatio[res_df$Log2EnrichRatio< -1] <- -1
  p <- ggplot(res_df)+
    geom_point(mapping = aes(x=x_indx, y=y_indx, color=Log2EnrichRatio, size=-1*log10(pval+1e-20))) +
    geom_point(data=data.frame(x_indx=1, y_indx=1, Log2EnrichRatio=c(1, -1), pval=1), mapping = aes(x=x_indx, y=y_indx, color=Log2EnrichRatio), size=NA) +
    geom_text(mapping = aes(x=x_indx, y=y_indx, label=OverlapNum), color="black", size=1.2) +
    scale_color_gradient2(low="blue", high="red", mid = "white", midpoint = 0, breaks=c(-1, -0.5, 0, 0.5, 1), labels=c("<=-1", "-0.5", "0", "0.5", ">=1")) +
    scale_size_continuous(range = c(0, 6)) +
    scale_x_continuous(breaks = 1:length(levels(res_df$Set1)), labels = levels(res_df$Set1)) +
    scale_y_continuous(breaks = 1:length(levels(res_df$Set2)), labels = levels(res_df$Set2)) +
    labs(color="log2 enrichment", size="ppval") +
    theme_bw() +
    theme(
      text = element_text(family="ArialMT", color = "black", size = 5),
      axis.title = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.key.size = unit(2, "mm"),
      legend.title = element_text(family="ArialMT", color = "black", size = 5),
      legend.text = element_text(family="ArialMT", color = "black", size = 5)
    )
  ggsave(file.path(prefix, "overlap.num.test.pdf"), p, width = 4+0.6*x_feature_num, height = 1+0.6*y_feature_num, limitsize = FALSE, units = "cm")
}

motif_peak_test_li <- list()
for (label1 in args$motif) {
  for (label2 in args$motif) {
    motif_peak_test_li[[sprintf("%s.%s", label1, label2)]] <- do.test(unlist(motif_target_peak_df[,label1]), unlist(motif_target_peak_df[,label2]), label1, label2)
  }
}
motif_peak_test_df <- do.call(rbind, motif_peak_test_li)
motif_peak_test_df$Set1 <- factor(motif_peak_test_df$Set1, levels = args$motif)
motif_peak_test_df$Set2 <- factor(motif_peak_test_df$Set2, levels = args$motif)
f_motif_target_peak_test <- file.path(args$output, "motif_target_peak_test")
if(!dir.exists(f_motif_target_peak_test)){
  dir.create(f_motif_target_peak_test)
}
write_tsv(motif_peak_test_df, file.path(f_motif_target_peak_test, "motif_target_peak.test.tsv"))
plot_test_res(motif_peak_test_df, f_motif_target_peak_test, pval_cutoff=1e-3)

ATAC_expr_peak <- apply(ATAC_peak_expr[,c("WT", "FL")], 1, function(x){return(max(x, na.rm=T))}) > 1
motif_expr_peak_test_li <- list()
for (label1 in args$motif) {
  for (label2 in args$motif) {
    motif_expr_peak_test_li[[sprintf("%s.%s", label1, label2)]] <- do.test(unlist(motif_target_peak_df[,label1])[ATAC_expr_peak], unlist(motif_target_peak_df[,label2])[ATAC_expr_peak], label1, label2)
  }
}
f_motif_expr_target_peak_test <- file.path(args$output, "motif_expr_target_peak")
if(!dir.exists(f_motif_expr_target_peak_test)){
  dir.create(f_motif_expr_target_peak_test)
}
motif_expr_peak_test_df <- do.call(rbind, motif_expr_peak_test_li)
write_tsv(motif_expr_peak_test_df, file.path(f_motif_expr_target_peak_test, "motif_expr_target_peak.test.tsv"))
plot_test_res(motif_expr_peak_test_df, f_motif_expr_target_peak_test, pval_cutoff=1e-3)

motif_gene_test_li <- list()
for (label1 in args$motif) {
  for (label2 in args$motif) {
    motif_gene_test_li[[sprintf("%s.%s", label1, label2)]] <- do.test(unlist(motif_target_gene_df[,label1]), unlist(motif_target_gene_df[,label2]), label1, label2)
  }
}
f_motif_gene_test <- file.path(args$output, "motif_target_gene")
if(!dir.exists(f_motif_gene_test)){
  dir.create(f_motif_gene_test)
}
motif_gene_test_df <- do.call(rbind, motif_gene_test_li)
write_tsv(motif_gene_test_df, file.path(f_motif_gene_test, "motif_target_gene.test.tsv"))
plot_test_res(motif_gene_test_df, f_motif_gene_test, pval_cutoff=1e-3)

scRNA_expr_gene <- apply(scRNA_expr[,2:ncol(scRNA_expr)], 1, function(x){return(max(x, na.rm=T))}) > 1
motif_expr_gene_test_li <- list()
for (label1 in args$motif) {
  for (label2 in args$motif) {
    motif_expr_gene_test_li[[sprintf("%s.%s", label1, label2)]] <- do.test(unlist(motif_target_gene_df[,label1][scRNA_expr_gene]), unlist(motif_target_gene_df[,label2][scRNA_expr_gene]), label1, label2)
  }
}
f_motif_expr_gene_test <- file.path(args$output, "motif_expr_target_gene")
if(!dir.exists(f_motif_expr_gene_test)){
  dir.create(f_motif_expr_gene_test)
}
motif_expr_gene_test_df <- do.call(rbind, motif_expr_gene_test_li)
write_tsv(motif_expr_gene_test_df, file.path(f_motif_expr_gene_test, "motif_expr_target_gene.test.tsv"))
plot_test_res(motif_expr_gene_test_df, f_motif_expr_gene_test, pval_cutoff=1e-3)

plot_expr_distribution <- function(expr_li, bool_df, prefix, ylab, min_expr=-5, max_expr=5, hline=c(-1, 0, 1)){
  expr_df <- cbind(data.frame(Expr=expr_li), bool_df)
  expr_df <- na.omit(expr_df)
  expr_df <- expr_df[order(expr_df$Expr),]
  expr_df$Order <- 1:nrow(expr_df)
  expr_df$Expr[expr_df$Expr<min_expr] <- min_expr
  expr_df$Expr[expr_df$Expr>max_expr] <- max_expr
  min_expr <- min(expr_df$Expr)
  max_expr <- max(expr_df$Expr)
  expr_range <- max_expr - min_expr
  expand_len <- max(10, nrow(expr_df)/100)
  # tmp_li <- list()
  # for(indx in 1:ncol(bool_df)){
  #   feature <- names(bool_df)[indx]
  #   tmp_li[[feature]] <- expr_df[,c(feature, "Order")]
  #   names(tmp_li[[feature]])[1] <- "Feature"
  #   tmp_li[[feature]]$AveRatio <- sapply(tmp_li[[feature]]$Order, function(indx){return(mean(tmp_li[[feature]]$Feature[abs(tmp_li[[feature]]$Order-indx)<=expand_len]))})
  #   tmp_li[[feature]]$FeatureIndex <- indx
  # }
  # tmp_df <- do.call(rbind, tmp_li)
  # tmp_df$ymin <- min_expr - tmp_df$FeatureIndex * ((0.2*expr_range)) + 0.2*0.35*expr_range
  # tmp_df$ymax <- min_expr - tmp_df$FeatureIndex * ((0.2*expr_range)) - 0.2*0.35*expr_range
  # label_df <- data.frame(Index=1:ncol(bool_df), Feature=names(bool_df))
  # label_df$y <-  min_expr - label_df$Index * ((0.2*expr_range))
  # 
  # p <- ggplot() +
  #   geom_line(data=expr_df, mapping = aes(x=Order, y=Expr), size=0.3, color="black") +
  #   geom_rect(data=tmp_df, mapping = aes(xmin=Order-0.5, xmax=Order+0.5, ymin=ymin, ymax=ymax, color=AveRatio)) +
  #   geom_text(data=label_df, mapping = aes(x=-1*expand_len, y=y, label=Feature), hjust=1, color="black", size=1.2) +
  #   geom_hline(yintercept = hline, color="black", size=0.3) +
  #   scale_color_gradient(low="white", high="black") +
  #   labs(y=ylab) +
  #   theme_bw() +
  #   theme(
  #     text = element_text(family="ArialMT", color = "black", size = 5),
  #     axis.title = element_blank(),
  #     panel.background = element_blank(),
  #     panel.grid = element_blank(),
  #     legend.key.size = unit(2, "mm"),
  #     legend.title = element_text(family="ArialMT", color = "black", size = 5),
  #     legend.text = element_text(family="ArialMT", color = "black", size = 5)
  #   )
  #   
  # ggsave(file.path(prefix, "feature.expr.distribution.pdf"), p, width = 10, height = 4+0.3*ncol(bool_df), limitsize = FALSE, units = "cm")
  expr_df <- data.frame()
  for(indx in 1:ncol(bool_df)){
    feature <- names(bool_df)[indx]
    tmp_df <- data.frame(Feature=feature, Expr=expr_li[unlist(bool_df[,indx])])
    expr_df <- rbind(expr_df, tmp_df)
  }
  expr_df$Expr[expr_df$Expr<min_expr] <- min_expr
  expr_df$Expr[expr_df$Expr>max_expr] <- max_expr
  expr_df <- na.omit(expr_df)
  expr_info <- expr_df %>% group_by(Feature) %>% summarise(N=n(), Mu=mean(Expr), Sd=sd(Expr), Pval=wilcox.test(Expr)$p.value)
  expr_info$Label <- sprintf("n=%s\nmu=%.2f\nsd=%.2f\np=%.2E", expr_info$N, expr_info$Mu, expr_info$Sd, expr_info$Pval)
  p <- ggplot(expr_df, aes(x=Feature, y=Expr, fill=Feature)) +
    geom_boxplot(outlier.colour = NA, size=0.3) +
    geom_text(data=expr_info, mapping = aes(x=Feature, y=max_expr, label=Label), size=1.2, color="black") +
    geom_hline(yintercept = hline, color="black", size=0.3) +
    labs(y=ylab) +
    theme_bw() +
    theme(
      text = element_text(family="ArialMT", color = "black", size = 5),
      axis.title.x = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(file.path(prefix, "feature.expr.box_plot.pdf"), p, width = 2+0.6*nrow(expr_info), height = 5, limitsize = FALSE, units = "cm")
  feature1_li <- c()
  feature2_li <- c()
  pval_li <- c()
  for (indx1 in 2:nrow(expr_info)) {
    for (indx2 in 1:indx1) {
      feature1 <- expr_info$Feature[indx1]
      feature2 <- expr_info$Feature[indx2]
      expr1 <- expr_df$Expr[expr_df$Feature==feature1]
      expr2 <- expr_df$Expr[expr_df$Feature==feature2]
      pval <- wilcox.test(expr1, expr2)$p.value
      feature1_li <- c(feature1_li, feature1)
      feature2_li <- c(feature2_li, feature2)
      pval_li <- c(pval_li, pval)
    }
  }
  wilcox_test_df <- data.frame(Feature1=feature1_li, Feature2=feature2_li, Pval=pval_li)
  write_tsv(wilcox_test_df, "feature.expr.wilcox_test.tsv")
}

plot_expr_distribution(log2( (ATAC_peak_expr$WT+1E-3) / (ATAC_peak_expr$FL+1E-3) ), motif_target_peak_df[,4:ncol(motif_target_peak_df)], f_motif_target_peak_test, "log2 ATAC peak signal (WT/FL)")
plot_expr_distribution(log2( (ATAC_peak_expr$WT+1E-3) / (ATAC_peak_expr$FL+1E-3) )[ATAC_expr_peak], motif_target_peak_df[ATAC_expr_peak, 4:ncol(motif_target_peak_df)], f_motif_expr_target_peak_test, "log2 ATAC peak signal (WT/FL)")

plot_expr_distribution(log2( (ATAC_gene_expr$WT+1E-3) / (ATAC_gene_expr$FL+1E-3) ), motif_target_gene_df[, 2:ncol(motif_target_gene_df)], f_motif_gene_test, "log2 ATAC gene signal (WT/FL)")
plot_expr_distribution(log2( (ATAC_gene_expr$WT+1E-3) / (ATAC_gene_expr$FL+1E-3) )[scRNA_expr_gene], motif_target_gene_df[scRNA_expr_gene, 2:ncol(motif_target_gene_df)], f_motif_expr_gene_test, "log2 ATAC gene signal (WT/FL)")

expr_ATAC_peak_info_df <- cbind(data.frame(log2FC=log2( (ATAC_peak_expr$WT+1E-3) / (ATAC_peak_expr$FL+1E-3) )[ATAC_expr_peak]), motif_target_peak_df[ATAC_expr_peak, ])
melt_expr_ATAC_peak_info_df <- melt(expr_ATAC_peak_info_df, id.vars = c("Chrom", "Start", "End", "log2FC"), variable.name = "Feature", value.name = "IsHit")
melt_expr_ATAC_peak_info_df$Tag <- ">1"
melt_expr_ATAC_peak_info_df$Tag[melt_expr_ATAC_peak_info_df$log2FC<1] <- "0~1"
melt_expr_ATAC_peak_info_df$Tag[melt_expr_ATAC_peak_info_df$log2FC<0] <- "-1~0"
melt_expr_ATAC_peak_info_df$Tag[melt_expr_ATAC_peak_info_df$log2FC< -1] <- "<-1"

expr_ATAC_peak_ave_ratio <- melt_expr_ATAC_peak_info_df %>% group_by(Feature) %>% summarise(AveRatio=mean(IsHit), TotalHit=sum(IsHit), TotalNotHit=sum(!IsHit))
melt_expr_ATAC_peak_info_stat <- melt_expr_ATAC_peak_info_df %>% group_by(Feature, Tag) %>% summarise(Num=n(), ObsHit=sum(IsHit), ObsNotHit=sum(!IsHit))
melt_expr_ATAC_peak_info_stat <- left_join(melt_expr_ATAC_peak_info_stat, expr_ATAC_peak_ave_ratio)
melt_expr_ATAC_peak_info_stat$ExpHit <- melt_expr_ATAC_peak_info_stat$Num*melt_expr_ATAC_peak_info_stat$AveRatio
melt_expr_ATAC_peak_info_stat$Pval <- apply(melt_expr_ATAC_peak_info_stat[,c("ObsHit", "ObsNotHit", "TotalHit", "TotalNotHit")], 1, function(x){return(fisher.test(matrix(x, nrow = 2))$p.value)})
write_tsv(melt_expr_ATAC_peak_info_stat, file.path(f_motif_expr_target_peak_test, "Feature.FC_level.fisher_test.tsv"))
melt_expr_ATAC_peak_info_stat$PvalLabel <- sprintf("p=%.1E", melt_expr_ATAC_peak_info_stat$Pval)
plot_df <- data.frame(
  Feature=rep(melt_expr_ATAC_peak_info_stat$Feature, 2), 
  Tag=rep(melt_expr_ATAC_peak_info_stat$Tag, 2),
  PeakNum=c(melt_expr_ATAC_peak_info_stat$ObsHit, melt_expr_ATAC_peak_info_stat$ExpHit),
  Data=c(rep("Obs", nrow(melt_expr_ATAC_peak_info_stat)), rep("Exp", nrow(melt_expr_ATAC_peak_info_stat)))
  )
plot_df$Feature <- factor(plot_df$Feature, levels = args$motif)
plot_df$Tag <- factor(plot_df$Tag, levels = c("<-1", "-1~0", "0~1", ">1"))
plot_df$Data <- factor(plot_df$Data, levels = c("Obs", "Exp"))
p <- ggplot(plot_df, aes(x=Data, y=PeakNum, fill=Data)) +
  geom_bar(stat="identity") +
  geom_text(mapping = aes(label=PeakNum), size=1.2, color="black") +
  geom_text(data=melt_expr_ATAC_peak_info_stat, mapping=aes(label=PvalLabel, x=1.5, y=0, fill=NA), size=1.2, color="black") +
  facet_wrap(Feature~Tag, scales = "free_y", nrow = length(args$motif)) +
  scale_fill_manual(values = c("Obs"="black", "Exp"="grey70")) +
  labs(y="#ATAC peaks") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(f_motif_expr_target_peak_test, "Feature.FC_level.fisher_test.pdf"), p, width = 7, height = 2*length(args$motif), limitsize = FALSE, units = "cm")

f_C3_test <- file.path(args$output, "motif_C3_expr")
if(!dir.exists(f_C3_test)){
  dir.create(f_C3_test)
}
mu_C3_expr <- mean(log10(scRNA_expr$WT.C3+1e-3)[scRNA_expr_gene])
sd_C3_expr <- sd(log10(scRNA_expr$WT.C3+1e-3)[scRNA_expr_gene])
plot_expr_distribution(log10(scRNA_expr$WT.C3+1e-3)[scRNA_expr_gene], motif_target_gene_df[scRNA_expr_gene, 2:ncol(motif_target_gene_df)], f_C3_test, "log10 WT.C3 expr", hline = c(mu_C3_expr, mu_C3_expr-sd_C3_expr, mu_C3_expr+sd_C3_expr))


motif_target_gene_df
scRNA_expr$IsC3Max <- apply(scRNA_expr[,2:11], 1, function(x){return(x[3]==max(x, na.rm = T))})
expr_target_gene <- inner_join(motif_target_gene_df, scRNA_expr[,c("Gid", "IsC3Max")])
write_tsv(expr_target_gene, file.path(args$output, "MotifExprTargetGene.HitMat.tsv"))

TF_gene_list <- read_delim(args$tf_ann, "\t", escape_double = FALSE, trim_ws = TRUE)
expr_target_gene_tf <- left_join(expr_target_gene, TF_gene_list)
all_gene_num <- length(unique(expr_target_gene_tf$Gid))
TF_num <- expr_target_gene_tf %>% group_by(Family) %>% summarise(AllFamilyNum=n(), NotFamilyNum=all_gene_num-AllFamilyNum, AllFamilyRatio=AllFamilyNum/all_gene_num)
MotifHitNum <- melt(expr_target_gene, id.vars = "Gid", variable.name = "Motif", value.name="IsHit") %>% group_by(Motif) %>% summarise(AllMotifNum=n(), AllMotifHitNum=sum(IsHit), AllMotifNotHitNum=sum(!IsHit), MotifHitRatio=mean(IsHit))

TF_target_gene <- expr_target_gene_tf[!is.na(expr_target_gene_tf$Family),]
melt_TF_target_gene <- melt(TF_target_gene[,c("Gid", args$motif, "Family")], id.vars = c("Gid", "Family"), variable.name = "Motif", value.name="IsHit")
TF_target_gene_info <- melt_TF_target_gene %>% group_by(Motif, Family) %>% summarise(AllTfNum=n(), TfTargetNum=sum(IsHit), TfNotTargetNum=sum(!IsHit))
TF_target_gene_info <- left_join(TF_target_gene_info, TF_num)
TF_target_gene_info <- left_join(TF_target_gene_info, MotifHitNum)
TF_target_gene_info$ExpHitNum <- TF_target_gene_info$AllTfNum * TF_target_gene_info$MotifHitRatio
TF_target_gene_info$Pval <- apply(TF_target_gene_info[,c("TfTargetNum", "TfNotTargetNum", "AllMotifHitNum", "AllMotifNotHitNum")], 1, function(x){return(fisher.test(matrix(x, nrow = 2))$p.value)})
write_tsv(TF_target_gene_info, file.path(args$output, "MotifExprTargetGene.TF.test.tsv"))
plot_df <- data.frame(Motif=rep(TF_target_gene_info$Motif, 2), Family=rep(TF_target_gene_info$Family, 2), GeneNum=c(TF_target_gene_info$TfTargetNum, TF_target_gene_info$ExpHitNum), Stat=c(rep("Obs", nrow(TF_target_gene_info)), rep("Exp", nrow(TF_target_gene_info))))
plot_df$Stat <- factor(plot_df$Stat, levels = c("Obs", "Exp"))
plot_df$Motif <- factor(plot_df$Motif, levels = args$motif)
p <- ggplot(plot_df, aes(x=Family, y=GeneNum, fill=Stat)) +
  geom_bar(stat="identity", position = position_dodge(0.7), width = 0.7) +
  facet_grid(Motif~., scales = "free") +
  scale_fill_manual(values = c("Obs"="black", "Exp"="grey70")) +
  labs(y="#TFs") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(family="ArialMT", color = "black", size = 5),
    legend.text = element_text(family="ArialMT", color = "black", size = 5),
    legend.key.size = unit(2, "mm")
  )
ggsave(file.path(args$output, "MotifExprTargetGene.TF.pdf"), p, width = 20, height = 2+2*length(args$motif), limitsize = FALSE, units = "cm")
enriched_test <- TF_target_gene_info[TF_target_gene_info$Pval<1e-3 & TF_target_gene_info$TfTargetNum>TF_target_gene_info$ExpHitNum,]
enriched_test2 <- TF_target_gene_info[TF_target_gene_info$Pval<1e-2 & TF_target_gene_info$Pval>=1e-3 & TF_target_gene_info$TfTargetNum>TF_target_gene_info$ExpHitNum,]
enriched_test3 <- TF_target_gene_info[TF_target_gene_info$Pval<0.05 & TF_target_gene_info$Pval>=1e-2 & TF_target_gene_info$TfTargetNum>TF_target_gene_info$ExpHitNum,]
MYB_G4_target <- melt_TF_target_gene[melt_TF_target_gene$Motif=="MYB_G4" & melt_TF_target_gene$IsHit,]
MYB_G4_target_MYB <- MYB_G4_target[MYB_G4_target$Family%in%c("MYB"),]
MYB_G4_target_MYB <- left_join(MYB_G4_target_MYB, InterestGeneID)

expr_scRNA <- scRNA_expr[scRNA_expr_gene,]
log2scRNA_df <- data.frame(
  Gid=expr_scRNA$Gid, 
  C1=log2((expr_scRNA$WT.C1+1e-3)/(expr_scRNA$FL.C1+1e-3)),
  C2=log2((expr_scRNA$WT.C2+1e-3)/(expr_scRNA$FL.C2+1e-3)),
  C4=log2((expr_scRNA$WT.C4+1e-3)/(expr_scRNA$FL.C4+1e-3)),
  C5=log2((expr_scRNA$WT.C5+1e-3)/(expr_scRNA$FL.C5+1e-3))
  )
melt_log2scRNA_df <- melt(log2scRNA_df, id.vars = "Gid", variable.name = "Cluster", value.name = "log2FC")

MYB_G4_log2scRNA_df <- left_join(melt_log2scRNA_df, data.frame(Motif="MYB_G4", Gid=expr_target_gene$Gid, IsHit=expr_target_gene$MYB_G4))
Gbox_log2scRNA_df <- left_join(melt_log2scRNA_df, data.frame(Motif="Gbox", Gid=expr_target_gene$Gid, IsHit=expr_target_gene$Gbox))

MYB_G4_Gbox_log2scRNA_df <- rbind(MYB_G4_log2scRNA_df, Gbox_log2scRNA_df)
p <- ggplot(MYB_G4_Gbox_log2scRNA_df, aes(x=Cluster, y=log2FC, fill=IsHit)) +
  geom_boxplot(outlier.colour = NA) +
  facet_grid(Motif~.) +
  coord_cartesian(ylim = c(-2, 2))
p
