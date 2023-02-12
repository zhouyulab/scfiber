library(argparse)
library(dplyr)
library(readr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(matrixStats)
library(circlize)
# library(splineTimeR)
library(ImpulseDE2)
# library(lmms)

parser <- ArgumentParser()
parser$add_argument("--tp-fc", required=TRUE, type="character", dest = "tp_fc", metavar="tp_fc.tsv")
parser$add_argument("--tp-cpm", required=TRUE, type="character", dest = "tp_cpm", metavar="tp_cpm.tsv")
parser$add_argument("--tp-cpm-rep", required=TRUE, type="character", dest = "tp_cpm_rep", metavar="tp_cpm_rep.tsv")
parser$add_argument("--ImpulseDE2-PADJ-cutoff", type="double", dest = "ImpulseDE2_padj_cutoff", metavar="ImpulseDE2_padj_cutoff", default=0.01)
parser$add_argument("--ImpulseDE2-mean-cutoff", type="double", dest = "ImpulseDE2_mean_cutoff", metavar="ImpulseDE2_mean_cutoff", default=1)
parser$add_argument("--interest-gene", required=TRUE, type="character", dest = "interest_gene", metavar="interest_gene.txt")
parser$add_argument("--tf-gene", required=TRUE, type="character", dest = "tf_gene", metavar="tf_gene.txt")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")

args <- commandArgs(TRUE)
args <- parser$parse_args(args)

args <- c(
  "--tp-fc",
  "analysis/RNA_seq_TP/merge_expr/feature_counts.merge.tsv",
  "--tp-cpm",
  "analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv",
  "--tp-cpm-rep",
  "analysis/RNA_seq_TP/circadian_gene/CPM.tsv",
  "--interest-gene",
  "data/InterestGeneID.tsv",
  "--tf-gene",
  "data/TF.gene_list.tsv",
  "-o", "analysis/RNA_seq_TP/TimeCourse"
)
args <- parser$parse_args(args)

fc_tp_df <- read_delim(args$tp_fc, "\t", escape_double = FALSE, trim_ws = TRUE)
cpm_tp_df <- read_delim(args$tp_cpm, "\t", escape_double = FALSE, trim_ws = TRUE)
cpm_tp_rep_df <- read_delim(args$tp_cpm_rep, "\t", escape_double = FALSE, trim_ws = TRUE)
tf_gene <- read_delim(args$tf_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
tf_gene <- tf_gene %>% group_by(Gid) %>% summarise(Family=Family[1])
interest_gene <- read_delim(args$interest_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
interest_gene <- interest_gene[!duplicated(interest_gene[,c("ID", "Gid")]),]
tp_rep <- c(-48, -48, -36, -36, -24, -24, -12, -12, 0, 0, 12, 12, 24, 24, 36, 36, 48, 48) + 48

## ImpulseDE2
## WT + FL (Not good)
fc_cnt_mat <- as.matrix(fc_tp_df[,2:ncol(fc_tp_df)])
rownames(fc_cnt_mat) <- fc_tp_df$Gid
# dfAnnotation <- data.frame(
#   Sample=colnames(fc_cnt_mat), 
#   Condition=c(rep("control", ncol(fc_cnt_mat)/2), rep("case", ncol(fc_cnt_mat)/2)),
#   Time=rep(tp_rep, 2),
#   Batch="B1"
# )
# ImpulseDE2_res_all_time <- runImpulseDE2(
#   matCountData    = fc_cnt_mat, 
#   dfAnnotation    = dfAnnotation,
#   boolCaseCtrl    = TRUE,
#   scaNProc        = 40 )
# ImpulseDE2_res_all_time_ht_case <- plotHeatmap(
#   objectImpulseDE2       = ImpulseDE2_res_all_time,
#   strCondition           = "case",
#   boolIdentifyTransients = FALSE,
#   scaQThres              = 0.01)
# draw(ImpulseDE2_res_all_time_ht_case$complexHeatmapRaw)
# ImpulseDE2_res_all_time_ht_control <- plotHeatmap(
#   objectImpulseDE2       = ImpulseDE2_res_all_time,
#   strCondition           = "control",
#   boolIdentifyTransients = FALSE,
#   scaQThres              = 0.01)
# draw(ImpulseDE2_res_all_time_ht_control$complexHeatmapRaw)

## WT only
fc_WT_cnt_mat <- fc_cnt_mat[,1:(ncol(fc_cnt_mat)/2)]
WT_dfAnnotation <- data.frame(
  Sample=colnames(fc_WT_cnt_mat), 
  Condition="case",
  Time=tp_rep,
  Batch="B1"
)
ImpulseDE2_res_WT_all_time <- runImpulseDE2(
  matCountData    = fc_WT_cnt_mat, 
  dfAnnotation    = WT_dfAnnotation,
  scaNProc        = 40 )
ImpulseDE2_res_WT_all_time_ht <- plotHeatmap(
  objectImpulseDE2       = ImpulseDE2_res_WT_all_time,
  strCondition           = "case",
  boolIdentifyTransients = FALSE,
  scaQThres              = 0.01)
draw(ImpulseDE2_res_WT_all_time_ht$complexHeatmapRaw)


fc_night_WT_cnt_mat <- fc_WT_cnt_mat[,tp_rep%%24==12]
night_WT_dfAnnotation <- data.frame(
  Sample=colnames(fc_night_WT_cnt_mat), 
  Condition="case",
  Time=tp_rep[tp_rep%%24==12],
  Batch="B1"
)
ImpulseDE2_res_WT_night_time <- runImpulseDE2(
  matCountData    = fc_night_WT_cnt_mat, 
  dfAnnotation    = night_WT_dfAnnotation,
  scaNProc        = 40 )
fc_light_WT_cnt_mat <- fc_WT_cnt_mat[,tp_rep%%24==0]
light_WT_dfAnnotation <- data.frame(
  Sample=colnames(fc_light_WT_cnt_mat), 
  Condition="case",
  Time=tp_rep[tp_rep%%24==0],
  Batch="B1"
)
ImpulseDE2_res_WT_light_time <- runImpulseDE2(
  matCountData    = fc_light_WT_cnt_mat, 
  dfAnnotation    = light_WT_dfAnnotation,
  scaNProc        = 40 )
ImpulseDE2_res_WT_night_time_ht <- plotHeatmap(
  objectImpulseDE2       = ImpulseDE2_res_WT_night_time,
  strCondition           = "case",
  boolIdentifyTransients = FALSE,
  scaQThres              = 0.01)
ImpulseDE2_res_WT_light_time_ht <- plotHeatmap(
  objectImpulseDE2       = ImpulseDE2_res_WT_light_time,
  strCondition           = "case",
  boolIdentifyTransients = FALSE,
  scaQThres              = 0.01)
draw(ImpulseDE2_res_WT_night_time_ht$complexHeatmapRaw)
draw(ImpulseDE2_res_WT_light_time_ht$complexHeatmapRaw)


filter_ImpulseDE2_res <- function(ImpulseDE2_obj, label, padj_cutoff=args$ImpulseDE2_padj_cutoff, mean_cutoff=args$ImpulseDE2_mean_cutoff){
  DE_res <- ImpulseDE2_obj$dfImpulseDE2Results
  filtered_res <- DE_res[!is.na(DE_res$padj) & DE_res$padj < padj_cutoff & DE_res$mean>mean_cutoff,]
  filtered_res$Label <- label
  return(filtered_res)
}

compute_continue_expr_time_num <- function(zscore_mat, label){
  expr_mat <- zscore_mat > 0
  continue_expr_num <- apply(expr_mat, 1, function(x){
    expr_indx <- sort(which(x))
    if(length(expr_indx)==0){
      return(0)
    }
    max_num <- 0
    tmp_num <- 0
    for(i in 1:length(expr_indx)){
      if(i == 1){
        tmp_num <- 1
      }else{
        last_indx <- expr_indx[i-1]
        this_indx <- expr_indx[i]
        if((this_indx-last_indx)==1){
          tmp_num <- tmp_num + 1
        }else{
          max_num <- max(max_num, tmp_num)
          tmp_num <- 1
        }
      }
    }
    max_num <- max(max_num, tmp_num)
    return(max_num)
  })
  res <- data.frame(Gid=rownames(zscore_mat), ContinueExprNum=continue_expr_num)
  names(res) <- c("Gid", label)
  return(res)
}

all_cpm_mat <- as.matrix(cpm_tp_df[,2:ncol(cpm_tp_df)])
WT_all_cpm_mat <- all_cpm_mat[,1:(ncol(all_cpm_mat)/2)]
rownames(WT_all_cpm_mat) <- cpm_tp_df$Gid
WT_all_zscore_mat <- t(scale(t(WT_all_cpm_mat)))
WT_night_cpm_mat <- WT_all_cpm_mat[,seq(2, ncol(WT_all_zscore_mat), 2)]
WT_night_zscore_mat <- t(scale(t(WT_night_cpm_mat)))
WT_light_cpm_mat <- WT_all_cpm_mat[,seq(1, ncol(WT_all_cpm_mat), 2)]
WT_light_zscore_mat <- t(scale(t(WT_light_cpm_mat)))
WT_all_continue_expr_num <- compute_continue_expr_time_num(WT_all_zscore_mat, "WtAllContinueExprNum")
WT_night_continue_expr_num <- compute_continue_expr_time_num(WT_night_zscore_mat, "WtNightContinueExprNum")
WT_light_continue_expr_num <- compute_continue_expr_time_num(WT_light_zscore_mat, "WtLightContinueExprNum")

ImpulseDE2_res_WT_all_time_fillter_df <- filter_ImpulseDE2_res(ImpulseDE2_res_WT_all_time, "WT All Time", padj_cutoff=args$ImpulseDE2_padj_cutoff, mean_cutoff=args$ImpulseDE2_mean_cutoff)
ImpulseDE2_res_WT_night_fillter_df <- filter_ImpulseDE2_res(ImpulseDE2_res_WT_night_time, "WT Night", padj_cutoff=args$ImpulseDE2_padj_cutoff, mean_cutoff=args$ImpulseDE2_mean_cutoff)
ImpulseDE2_res_WT_light_fillter_df <- filter_ImpulseDE2_res(ImpulseDE2_res_WT_light_time, "WT Light", padj_cutoff=args$ImpulseDE2_padj_cutoff, mean_cutoff=args$ImpulseDE2_mean_cutoff)
all_ImpulseDE2_res <- rbind(ImpulseDE2_res_WT_all_time_fillter_df, ImpulseDE2_res_WT_night_fillter_df)
all_ImpulseDE2_res <- rbind(all_ImpulseDE2_res, ImpulseDE2_res_WT_light_fillter_df)
names(all_ImpulseDE2_res)[1] <- "Gid"
all_ImpulseDE2_res <- left_join(all_ImpulseDE2_res, WT_all_continue_expr_num)
all_ImpulseDE2_res <- left_join(all_ImpulseDE2_res, WT_night_continue_expr_num)
all_ImpulseDE2_res <- left_join(all_ImpulseDE2_res, WT_light_continue_expr_num)
all_ImpulseDE2_res$Pass <- FALSE
all_ImpulseDE2_res$Pass[all_ImpulseDE2_res$Label=="WT All Time" & all_ImpulseDE2_res$WtAllContinueExprNum>=4] <- TRUE
all_ImpulseDE2_res$Pass[all_ImpulseDE2_res$Label=="WT Night" & all_ImpulseDE2_res$WtNightContinueExprNum>=2] <- TRUE
all_ImpulseDE2_res$Pass[all_ImpulseDE2_res$Label=="WT Light" & all_ImpulseDE2_res$WtLightContinueExprNum>=2] <- TRUE
all_ImpulseDE2_res <- all_ImpulseDE2_res[all_ImpulseDE2_res$Pass,]
table(all_ImpulseDE2_res$Label)
write_tsv(all_ImpulseDE2_res, file.path(args$output, "ImpulseDE2.enrich.SourceData.tsv"))

ImpulseDE2_enrich_gene_df <- data.frame(Gid=unique(all_ImpulseDE2_res$Gid))
ImpulseDE2_enrich_gene_df$IsAllTimeChangeGene <- ImpulseDE2_enrich_gene_df$Gid %in% all_ImpulseDE2_res$Gid[all_ImpulseDE2_res$Label=="WT All Time"]
ImpulseDE2_enrich_gene_df$IsLightChangeGene <- ImpulseDE2_enrich_gene_df$Gid %in% all_ImpulseDE2_res$Gid[all_ImpulseDE2_res$Label=="WT Light"]
ImpulseDE2_enrich_gene_df$IsNightChangeGene <- ImpulseDE2_enrich_gene_df$Gid %in% all_ImpulseDE2_res$Gid[all_ImpulseDE2_res$Label=="WT Night"]

ImpulseDE2_enrich_gene_source_info <- ImpulseDE2_enrich_gene_df %>% group_by(IsAllTimeChangeGene, IsLightChangeGene, IsNightChangeGene) %>% summarise(GeneNum=n())

compute_gene_order <- function(cpm_mat, zscore_weight=0.4){
  weight_light_time <- rowSums(t(1:ncol(cpm_mat) * t(cpm_mat))) / rowSums(cpm_mat)
  time_dist <- dist(weight_light_time)
  time_dist[is.na(time_dist)] <- max(time_dist, na.rm = T)
  zscore_mat <- t(scale(t(cpm_mat)))
  zscore_dist <- dist(zscore_mat)
  zscore_dist[is.na(zscore_dist)] <- max(zscore_dist, na.rm = T)
  merge_dist <- zscore_dist / median(zscore_dist) * zscore_weight + time_dist / median(time_dist)
  gene_clu <- hclust(merge_dist, method = "ward.D2")
  gene_order <- gene_clu$order
  return(gene_order)
}

ImpulseDE2_enrich_gene_cpm <- left_join(ImpulseDE2_enrich_gene_df, cpm_tp_df)
ImpulseDE2_enrich_gene_cpm_mat <- as.matrix(ImpulseDE2_enrich_gene_cpm[,5:ncol(ImpulseDE2_enrich_gene_cpm)])
rownames(ImpulseDE2_enrich_gene_cpm_mat) <- ImpulseDE2_enrich_gene_cpm$Gid

ImpulseDE2_enrich_gene_WT_cpm_mat <- ImpulseDE2_enrich_gene_cpm_mat[,1:(ncol(ImpulseDE2_enrich_gene_cpm_mat)/2)]
ImpulseDE2_enrich_gene_FL_cpm_mat <- ImpulseDE2_enrich_gene_cpm_mat[,(ncol(ImpulseDE2_enrich_gene_cpm_mat)/2+1):ncol(ImpulseDE2_enrich_gene_cpm_mat)]
gene_order <- compute_gene_order(ImpulseDE2_enrich_gene_WT_cpm_mat)

ImpulseDE2_enrich_gene_cpm <- ImpulseDE2_enrich_gene_cpm[gene_order,]
ImpulseDE2_enrich_gene_cpm_mat <- ImpulseDE2_enrich_gene_cpm_mat[gene_order,]
ImpulseDE2_enrich_gene_WT_cpm_mat <- ImpulseDE2_enrich_gene_WT_cpm_mat[gene_order,]
ImpulseDE2_enrich_gene_FL_cpm_mat <- ImpulseDE2_enrich_gene_FL_cpm_mat[gene_order,]

ImpulseDE2_enrich_gene_zscore_mat <- t(scale(t(ImpulseDE2_enrich_gene_cpm_mat)))
ImpulseDE2_enrich_gene_zscore_mat[rowMaxs(ImpulseDE2_enrich_gene_cpm_mat)<1,] <- NA
ImpulseDE2_enrich_gene_WT_zscore_mat <- ImpulseDE2_enrich_gene_zscore_mat[,1:(ncol(ImpulseDE2_enrich_gene_zscore_mat)/2)]
ImpulseDE2_enrich_gene_FL_zscore_mat <- ImpulseDE2_enrich_gene_zscore_mat[,(ncol(ImpulseDE2_enrich_gene_zscore_mat)/2+1):ncol(ImpulseDE2_enrich_gene_zscore_mat)]

ImpulseDE2_enrich_gene_group_df <- data.frame(Gid=rownames(ImpulseDE2_enrich_gene_WT_cpm_mat), Group="NotSure")
WT_left_score <- rowMeans(ImpulseDE2_enrich_gene_WT_cpm_mat[,1:4])
WT_left_score[is.na(WT_left_score)] <- 0
WT_right_score <- rowMeans(ImpulseDE2_enrich_gene_WT_cpm_mat[,(ncol(ImpulseDE2_enrich_gene_WT_cpm_mat)-4):ncol(ImpulseDE2_enrich_gene_WT_cpm_mat)])
WT_right_score[is.na(WT_right_score)] <- 0
FL_left_score <- rowMeans(ImpulseDE2_enrich_gene_FL_cpm_mat[,1:4])
FL_left_score[is.na(FL_left_score)] <- 0
FL_right_score <- rowMeans(ImpulseDE2_enrich_gene_FL_cpm_mat[,(ncol(ImpulseDE2_enrich_gene_FL_cpm_mat)-4):ncol(ImpulseDE2_enrich_gene_FL_cpm_mat)])
FL_right_score[is.na(FL_right_score)] <- 0

ImpulseDE2_enrich_gene_group_df$Group <- as.character(ImpulseDE2_enrich_gene_group_df$Group)
ImpulseDE2_enrich_gene_group_df$Group[WT_left_score > (WT_right_score * 1.1) & WT_left_score > 0] <- "Early"
ImpulseDE2_enrich_gene_group_df$Group[(WT_left_score * 1.1) < WT_right_score & WT_right_score > 0] <- "Late"
table(ImpulseDE2_enrich_gene_group_df$Group)
ImpulseDE2_enrich_gene_group_df <- ImpulseDE2_enrich_gene_group_df[ImpulseDE2_enrich_gene_group_df$Group %in% c("Early", "Late"),]
ImpulseDE2_enrich_gene_group_df$Group <- factor(ImpulseDE2_enrich_gene_group_df$Group, levels=c("Early", "Late"))

all_ImpulseDE2_res <- all_ImpulseDE2_res[all_ImpulseDE2_res$Gid %in% ImpulseDE2_enrich_gene_group_df$Gid,]
ImpulseDE2_enrich_gene_df <- data.frame(Gid=unique(all_ImpulseDE2_res$Gid))
ImpulseDE2_enrich_gene_df$IsAllTimeChangeGene <- ImpulseDE2_enrich_gene_df$Gid %in% all_ImpulseDE2_res$Gid[all_ImpulseDE2_res$Label=="WT All Time"]
ImpulseDE2_enrich_gene_df$IsLightChangeGene <- ImpulseDE2_enrich_gene_df$Gid %in% all_ImpulseDE2_res$Gid[all_ImpulseDE2_res$Label=="WT Light"]
ImpulseDE2_enrich_gene_df$IsNightChangeGene <- ImpulseDE2_enrich_gene_df$Gid %in% all_ImpulseDE2_res$Gid[all_ImpulseDE2_res$Label=="WT Night"]

ImpulseDE2_enrich_gene_source_info <- ImpulseDE2_enrich_gene_df %>% group_by(IsAllTimeChangeGene, IsLightChangeGene, IsNightChangeGene) %>% summarise(GeneNum=n())

ImpulseDE2_enrich_gene_cpm <- left_join(ImpulseDE2_enrich_gene_df, cpm_tp_df)
ImpulseDE2_enrich_gene_cpm_mat <- as.matrix(ImpulseDE2_enrich_gene_cpm[,5:ncol(ImpulseDE2_enrich_gene_cpm)])
rownames(ImpulseDE2_enrich_gene_cpm_mat) <- ImpulseDE2_enrich_gene_cpm$Gid

ImpulseDE2_enrich_gene_WT_cpm_mat <- ImpulseDE2_enrich_gene_cpm_mat[,1:(ncol(ImpulseDE2_enrich_gene_cpm_mat)/2)]
ImpulseDE2_enrich_gene_FL_cpm_mat <- ImpulseDE2_enrich_gene_cpm_mat[,(ncol(ImpulseDE2_enrich_gene_cpm_mat)/2+1):ncol(ImpulseDE2_enrich_gene_cpm_mat)]
gene_order <- compute_gene_order(ImpulseDE2_enrich_gene_WT_cpm_mat)

ImpulseDE2_enrich_gene_cpm <- ImpulseDE2_enrich_gene_cpm[gene_order,]
ImpulseDE2_enrich_gene_cpm_mat <- ImpulseDE2_enrich_gene_cpm_mat[gene_order,]
ImpulseDE2_enrich_gene_WT_cpm_mat <- ImpulseDE2_enrich_gene_WT_cpm_mat[gene_order,]
ImpulseDE2_enrich_gene_FL_cpm_mat <- ImpulseDE2_enrich_gene_FL_cpm_mat[gene_order,]

ImpulseDE2_enrich_gene_zscore_mat <- t(scale(t(ImpulseDE2_enrich_gene_cpm_mat)))
ImpulseDE2_enrich_gene_zscore_mat[rowMaxs(ImpulseDE2_enrich_gene_cpm_mat)<1,] <- NA
ImpulseDE2_enrich_gene_WT_zscore_mat <- ImpulseDE2_enrich_gene_zscore_mat[,1:(ncol(ImpulseDE2_enrich_gene_zscore_mat)/2)]
ImpulseDE2_enrich_gene_FL_zscore_mat <- ImpulseDE2_enrich_gene_zscore_mat[,(ncol(ImpulseDE2_enrich_gene_zscore_mat)/2+1):ncol(ImpulseDE2_enrich_gene_zscore_mat)]

ImpulseDE2_enrich_gene_group_df <- data.frame(Gid=rownames(ImpulseDE2_enrich_gene_WT_cpm_mat), Group="NotSure")
WT_left_score <- rowMeans(ImpulseDE2_enrich_gene_WT_cpm_mat[,1:4])
WT_left_score[is.na(WT_left_score)] <- 0
WT_right_score <- rowMeans(ImpulseDE2_enrich_gene_WT_cpm_mat[,(ncol(ImpulseDE2_enrich_gene_WT_cpm_mat)-4):ncol(ImpulseDE2_enrich_gene_WT_cpm_mat)])
WT_right_score[is.na(WT_right_score)] <- 0
FL_left_score <- rowMeans(ImpulseDE2_enrich_gene_FL_cpm_mat[,1:4])
FL_left_score[is.na(FL_left_score)] <- 0
FL_right_score <- rowMeans(ImpulseDE2_enrich_gene_FL_cpm_mat[,(ncol(ImpulseDE2_enrich_gene_FL_cpm_mat)-4):ncol(ImpulseDE2_enrich_gene_FL_cpm_mat)])
FL_right_score[is.na(FL_right_score)] <- 0

ImpulseDE2_enrich_gene_group_df$Group <- as.character(ImpulseDE2_enrich_gene_group_df$Group)
ImpulseDE2_enrich_gene_group_df$Group[WT_left_score > (WT_right_score * 1.1) & WT_left_score > 0] <- "Early"
ImpulseDE2_enrich_gene_group_df$Group[(WT_left_score * 1.1) < WT_right_score & WT_right_score > 0] <- "Late"
table(ImpulseDE2_enrich_gene_group_df$Group)
ImpulseDE2_enrich_gene_group_df <- ImpulseDE2_enrich_gene_group_df[ImpulseDE2_enrich_gene_group_df$Group %in% c("Early", "Late"),]
ImpulseDE2_enrich_gene_group_df$Group <- factor(ImpulseDE2_enrich_gene_group_df$Group, levels=c("Early", "Late"))


wt_zscore_ht <- Heatmap(
  ImpulseDE2_enrich_gene_WT_zscore_mat,
  name="WT expr. z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  split = ImpulseDE2_enrich_gene_group_df$Group,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
    )
)
wt_zscore_ht

fl_zscore_ht <- Heatmap(
  ImpulseDE2_enrich_gene_FL_zscore_mat,
  name="FL expr. z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  split = ImpulseDE2_enrich_gene_group_df$Group,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
fl_zscore_ht

interest_gene_df <- left_join(ImpulseDE2_enrich_gene_group_df, interest_gene)
interest_gene_df <- left_join(interest_gene_df, tf_gene)
interest_gene_df$Label <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$GeneSymbol)
interest_gene_df$Label[!is.na(interest_gene_df$Family)] <- sprintf("%s (%s: %s)", interest_gene_df$Gid, interest_gene_df$Family, interest_gene_df$GeneSymbol)[!is.na(interest_gene_df$Family)]
interest_gene_df$Label[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$Family)[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]
interest_gene_df$Label[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- interest_gene_df$Gid[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]

subindx <- which(!is.na(interest_gene_df$GeneSymbol))
label_li <- interest_gene_df$Label[subindx]
ha <- rowAnnotation(
  link = row_anno_link(
  at = subindx,
  labels = label_li,
  labels_gp=gpar(fontsize = 6),
  lines_gp=gpar(lwd = 0.5),
  link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)))

inche_cm <- 2.54
pdf(file.path(args$output, "ImpulseDE2_enrich_gene.zscore.all.pdf"), width=15/inche_cm, height=20/inche_cm)
print(wt_zscore_ht + fl_zscore_ht + ha)
dev.off()


ImpulseDE2_enrich_gene_WT_zscore_light_mat <- ImpulseDE2_enrich_gene_WT_zscore_mat[,seq(1, ncol(ImpulseDE2_enrich_gene_WT_zscore_mat), 2)]
ImpulseDE2_enrich_gene_FL_zscore_light_mat <- ImpulseDE2_enrich_gene_FL_zscore_mat[,seq(1, ncol(ImpulseDE2_enrich_gene_FL_zscore_mat), 2)]

ImpulseDE2_enrich_gene_WT_cpm_light_mat <- ImpulseDE2_enrich_gene_WT_cpm_mat[,seq(1, ncol(ImpulseDE2_enrich_gene_WT_cpm_mat), 2)]
all(rownames(ImpulseDE2_enrich_gene_WT_cpm_light_mat) == rownames(ImpulseDE2_enrich_gene_WT_zscore_light_mat))
light_gene_order <- compute_gene_order(ImpulseDE2_enrich_gene_WT_cpm_light_mat, zscore_weight=0.4)

ImpulseDE2_enrich_gene_WT_zscore_light_mat <- ImpulseDE2_enrich_gene_WT_zscore_light_mat[light_gene_order,]
ImpulseDE2_enrich_gene_FL_zscore_light_mat <- ImpulseDE2_enrich_gene_FL_zscore_light_mat[light_gene_order,]

wt_zscore_light_ht <- Heatmap(
  ImpulseDE2_enrich_gene_WT_zscore_light_mat,
  name="WT expr. z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  split = ImpulseDE2_enrich_gene_group_df$Group[light_gene_order],
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
wt_zscore_light_ht

fl_zscore_light_ht <- Heatmap(
  ImpulseDE2_enrich_gene_FL_zscore_light_mat,
  name="FL expr. z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  split = ImpulseDE2_enrich_gene_group_df$Group[light_gene_order],
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
fl_zscore_light_ht

light_interest_gene_df <- interest_gene_df[light_gene_order,]
light_subindx <- which(!is.na(light_interest_gene_df$GeneSymbol))
light_label_li <- light_interest_gene_df$Label[light_subindx]
light_ha <- rowAnnotation(link = row_anno_link(
  at = light_subindx,
  labels = light_label_li,
  labels_gp=gpar(fontsize = 6),
  lines_gp=gpar(lwd = 0.5),
  link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(light_label_li, gp = gpar(fontsize = 6)))

pdf(file.path(args$output, "ImpulseDE2_enrich_gene.zscore.light.pdf"), width=12/inche_cm, height=20/inche_cm)
print(wt_zscore_light_ht + fl_zscore_light_ht + light_ha)
dev.off()

ImpulseDE2_enrich_gene_WT_zscore_night_mat <- ImpulseDE2_enrich_gene_WT_zscore_mat[,seq(2, ncol(ImpulseDE2_enrich_gene_WT_zscore_mat), 2)]
ImpulseDE2_enrich_gene_FL_zscore_night_mat <- ImpulseDE2_enrich_gene_FL_zscore_mat[,seq(2, ncol(ImpulseDE2_enrich_gene_FL_zscore_mat), 2)]
ImpulseDE2_enrich_gene_WT_cpm_night_mat <- ImpulseDE2_enrich_gene_WT_cpm_mat[,seq(2, ncol(ImpulseDE2_enrich_gene_WT_cpm_mat), 2)]
all(rownames(ImpulseDE2_enrich_gene_WT_cpm_night_mat) == rownames(ImpulseDE2_enrich_gene_WT_zscore_night_mat))
night_gene_order <- compute_gene_order(ImpulseDE2_enrich_gene_WT_cpm_night_mat, zscore_weight=0.4)

ImpulseDE2_enrich_gene_WT_zscore_night_mat <- ImpulseDE2_enrich_gene_WT_zscore_night_mat[night_gene_order,]
ImpulseDE2_enrich_gene_FL_zscore_night_mat <- ImpulseDE2_enrich_gene_FL_zscore_night_mat[night_gene_order,]

wt_zscore_night_ht <- Heatmap(
  ImpulseDE2_enrich_gene_WT_zscore_night_mat,
  name="WT expr. z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  split = ImpulseDE2_enrich_gene_group_df$Group[night_gene_order],
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
wt_zscore_night_ht

fl_zscore_night_ht <- Heatmap(
  ImpulseDE2_enrich_gene_FL_zscore_night_mat,
  name="FL expr. z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  split = ImpulseDE2_enrich_gene_group_df$Group[night_gene_order],
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
fl_zscore_night_ht

night_interest_gene_df <- interest_gene_df[night_gene_order,]
night_subindx <- which(!is.na(night_interest_gene_df$GeneSymbol))
night_label_li <- night_interest_gene_df$Label[night_subindx]
night_ha <- rowAnnotation(link = row_anno_link(
  at = night_subindx,
  labels = night_label_li,
  labels_gp=gpar(fontsize = 6),
  lines_gp=gpar(lwd = 0.5),
  link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(night_label_li, gp = gpar(fontsize = 6)))

pdf(file.path(args$output, "ImpulseDE2_enrich_gene.zscore.night.pdf"), width=10/inche_cm, height=20/inche_cm)
print(wt_zscore_night_ht + fl_zscore_night_ht + night_ha)
dev.off()

expr_time_region <- apply( t(scale(t(ImpulseDE2_enrich_gene_WT_cpm_light_mat))), 1, function(x){
  x[is.na(x)] <- FALSE
  is_expr <- x>0
  for(i in 1:length(x)){
    if(i==1){
      if(!is_expr[i+1]){
        is_expr[i] <- FALSE
      }
    }else{
      if(i==length(x)){
        if(!is_expr[i-1]){
          is_expr[i] <- FALSE
        }
      }else{
        if((!is_expr[i-1]) & (!is_expr[i+1])){{
          is_expr[i] <- FALSE
        }}
      }
    }
  }
  expr_indx <- which(is_expr)
  start_time <- NA
  end_time <- NA
  if(length(expr_indx)>0 & !(is_expr[1] & is_expr[length(x)])){
    start_time <- names(x)[min(expr_indx)]
    end_time <- names(x)[max(expr_indx)]
  }
  res <- data.frame(LightStartTime=start_time, LightEndTime=end_time)
  return(res)
})
expr_time_region <- do.call(rbind, expr_time_region)
expr_time_region$Gid <- rownames(expr_time_region)

ImpulseDE2_enrich_gene_cpm_group <- left_join(ImpulseDE2_enrich_gene_cpm, ImpulseDE2_enrich_gene_group_df)
ImpulseDE2_enrich_gene_cpm_group <- left_join(ImpulseDE2_enrich_gene_cpm_group, expr_time_region)
write_tsv(ImpulseDE2_enrich_gene_cpm_group, file.path(args$output, "ImpulseDE2_enrich_gene.CPM.tsv"))

expr_time_info <- expr_time_region %>% group_by(LightStartTime, LightEndTime) %>% summarise(GeneNum=n())
expr_time_info$LightStartTime <- as.character(expr_time_info$LightStartTime)
expr_time_info$LightEndTime <- as.character(expr_time_info$LightEndTime)
expr_time_info[is.na(expr_time_info)] <- "NotLightTimeCouseGene"
expr_time_info$LightStartTime <- factor(expr_time_info$LightStartTime, levels = c(colnames(ImpulseDE2_enrich_gene_WT_cpm_light_mat), "NotLightTimeCouseGene"))
expr_time_info$LightEndTime <- factor(expr_time_info$LightEndTime, levels = c(colnames(ImpulseDE2_enrich_gene_WT_cpm_light_mat), "NotLightTimeCouseGene"))
expr_time_info <- expr_time_info[order(expr_time_info$LightStartTime, expr_time_info$LightEndTime),]
write_tsv(expr_time_info, file.path(args$output, "ImpulseDE2_enrich_gene.light_expr_time.num.tsv"))
expr_time_info$Label <- sprintf("%s - %s", expr_time_info$LightStartTime, expr_time_info$LightEndTime)
expr_time_info$Label[expr_time_info$LightStartTime=="NotLightTimeCouseGene"] <- "NotLightTimeCouseGene"
expr_time_info$Label <- factor(expr_time_info$Label, levels = rev(expr_time_info$Label))
p <-  ggplot(expr_time_info, aes(x=Label, y=GeneNum)) +
  geom_bar(stat="identity", fill="black") +
  geom_text(mapping = aes(label=GeneNum), color="white", size=1.4, hjust=1) +
  coord_flip() +
  labs(x="Expr. Time (WT light)", y="#Genes") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 5),
        axis.title = element_text(family="ArialMT", color = "black", size = 5),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.grid = element_blank())
ggsave(file.path(args$output, "ImpulseDE2_enrich_gene.light_expr_time.num.pdf"), p, width = 6, height = 4, units = "cm")


ImpulseDE2_enrich_gene_cpm_group_tf <- inner_join(tf_gene, ImpulseDE2_enrich_gene_cpm_group)
write_tsv(ImpulseDE2_enrich_gene_cpm_group_tf, file.path(args$output, "ImpulseDE2_enrich_gene.TF.CPM.tsv"))
ImpulseDE2_enrich_gene_cpm_group_tf_info <- ImpulseDE2_enrich_gene_cpm_group_tf %>% group_by(Family, Group) %>% summarise(GeneNum=n())
ImpulseDE2_enrich_gene_cpm_group_tf_info_total <- ImpulseDE2_enrich_gene_cpm_group_tf %>% group_by(Family) %>% summarise(GeneNum=n())
ImpulseDE2_enrich_gene_cpm_group_tf_info_total <- ImpulseDE2_enrich_gene_cpm_group_tf_info_total[order(ImpulseDE2_enrich_gene_cpm_group_tf_info_total$GeneNum, decreasing = TRUE),]
ImpulseDE2_enrich_gene_cpm_group_tf_info$Family <- factor(ImpulseDE2_enrich_gene_cpm_group_tf_info$Family, levels = ImpulseDE2_enrich_gene_cpm_group_tf_info_total$Family)
ImpulseDE2_enrich_gene_cpm_group_tf_info$Group <- factor(ImpulseDE2_enrich_gene_cpm_group_tf_info$Group, levels = c("Late", "Early"))
p <- ggplot(ImpulseDE2_enrich_gene_cpm_group_tf_info, aes(x=Family, y=GeneNum, fill=Group)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("Early"="red", "Late"="blue")) +
  theme_bw() +
  labs(y="#TF genes") +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = c(0.8, 0.8),
        legend.key.size = unit(4, "mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(family="ArialMT", size=6),
        panel.grid = element_blank()
  )
ggsave(file.path(args$output, "ImpulseDE2_enrich_gene.TF_gene.num.pdf"), p, width = 10, height = 5, units = "cm")

ImpulseDE2_enrich_gene_group_tf_stat <- ImpulseDE2_enrich_gene_cpm_group_tf_info %>% group_by(Family) %>% summarise(TotalGeneNum=sum(GeneNum), EarlyGeneNum=sum(GeneNum[Group=="Early"]), EarlyGeneRatio= EarlyGeneNum/TotalGeneNum)
ImpulseDE2_enrich_gene_group_tf_stat$TF_group <- "Both expr."
ImpulseDE2_enrich_gene_group_tf_stat$TF_group[ImpulseDE2_enrich_gene_group_tf_stat$EarlyGeneRatio>=0.75] <- "Early expr."
ImpulseDE2_enrich_gene_group_tf_stat$TF_group[ImpulseDE2_enrich_gene_group_tf_stat$EarlyGeneRatio<=0.25] <- "Late expr."
p <- ggplot(ImpulseDE2_enrich_gene_group_tf_stat, aes(x=TotalGeneNum, y=EarlyGeneRatio, color=TF_group)) +
  geom_hline(yintercept = c(0.25, 0.75), color="black", size=0.1, lty=2) +
  geom_point(size=0.2) +
  geom_text(mapping = aes(label=Family), color="black", size=1.2, vjust=2) +
  scale_color_manual(values = c("Early expr."="red", "Late expr."="blue", "Both expr."="grey70")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%"), limits = c(0, 1)) +
  theme_bw() +
  labs(x="#TF genes", y="%Early TF genes") +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(4, "mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(family="ArialMT", size=6),
        panel.grid = element_blank()
  )
ggsave(file.path(args$output, "ImpulseDE2_enrich_gene.TF_gene.expr_state.pdf"), p, width = 8, height = 5, units = "cm")

FL_WT_log2FC_mat <- log2((ImpulseDE2_enrich_gene_FL_cpm_mat+1e-3) / (ImpulseDE2_enrich_gene_WT_cpm_mat+1e-3))
FL_WT_log2FC_df <- as.data.frame(FL_WT_log2FC_mat)
names(FL_WT_log2FC_df) <- sapply(strsplit(names(FL_WT_log2FC_df), "_"), function(x){return(x[2])})
tp_li <- names(FL_WT_log2FC_df)
FL_WT_log2FC_df$Gid <- rownames(FL_WT_log2FC_df)
melt_FL_WT_log2FC_df <- melt(FL_WT_log2FC_df, id.vars = "Gid", variable.name = "TP", value.name = "log2FC")
melt_FL_WT_log2FC_df <- left_join(melt_FL_WT_log2FC_df, ImpulseDE2_enrich_gene_group_df)
melt_FL_WT_log2FC_df$TP <- factor(melt_FL_WT_log2FC_df$TP, levels = tp_li)

FL_WT_log2FC_info <- melt_FL_WT_log2FC_df %>% 
  group_by(TP) %>% 
  summarise(
    MuEarly=mean(log2FC[Group=="Early"]),
    MuLate=mean(log2FC[Group=="Late"]),
    EarlyPval=t.test(log2FC[Group=="Early"])$p.value,
    LatePval=t.test(log2FC[Group=="Late"])$p.value
    )
FL_WT_log2FC_info$Label <- sprintf("MuEarly=%.2f\nMuLate=%.2f\nEarlyPval=%.1E\nLatePval=%.1E", FL_WT_log2FC_info$MuEarly, FL_WT_log2FC_info$MuLate, FL_WT_log2FC_info$EarlyPval, FL_WT_log2FC_info$LatePval)

p <- ggplot(melt_FL_WT_log2FC_df, aes(x=TP, y=log2FC)) +
  geom_hline(yintercept = 0, color="black", size=0.1, lty=2) +
  geom_boxplot(mapping = aes(fill=Group), outlier.colour = NA, size=0.2) +
  geom_text(data=FL_WT_log2FC_info, aes(y=2.3, label=Label), size=1.2, color="black", vjust=0) +
  coord_cartesian(ylim = c(-3, 3)) +
  scale_fill_manual(values = c("Early"="red", "Late"="blue")) +
  theme_bw() +
  labs(y="log2(FL expr. / WT expr.)") +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(family="ArialMT", size=6),
        panel.grid = element_blank()
  )
ggsave(file.path(args$output, "ImpulseDE2_enrich_gene.FL_WT_expr_FC.pdf"), p, width = 12, height = 5, units = "cm")

wilcox_pval <- function(x){
  if(length(x)<=2){return(NA)}
  return(t.test(x)$p.value)
}
tf_FL_WT_log2FC_df <- inner_join(melt_FL_WT_log2FC_df, tf_gene)
tf_FL_WT_log2FC_info <- tf_FL_WT_log2FC_df %>% 
  group_by(Family, TP) %>% 
  summarise(
    nEarly=sum(Group=="Early"),
    nLate=sum(Group=="Late"),
    MuEarly=mean(log2FC[Group=="Early"]),
    MuLate=mean(log2FC[Group=="Late"]),
    EarlyPval=wilcox_pval(log2FC[Group=="Early"]),
    LatePval=wilcox_pval(log2FC[Group=="Late"])
  )
tf_FL_WT_log2FC_info$Label <- sprintf(
  "nEarly=%d\nnLate=%d\nMuEarly=%.2f\nMuLate=%.2f\nEarlyPval=%.1E\nLatePval=%.1E", 
  tf_FL_WT_log2FC_info$nEarly,
  tf_FL_WT_log2FC_info$nLate,
  tf_FL_WT_log2FC_info$MuEarly, 
  tf_FL_WT_log2FC_info$MuLate, 
  tf_FL_WT_log2FC_info$EarlyPval, 
  tf_FL_WT_log2FC_info$LatePval
  )


p <- ggplot(tf_FL_WT_log2FC_df, aes(x=TP, y=log2FC)) +
  geom_hline(yintercept = 0, color="black", size=0.1, lty=2) +
  geom_boxplot(mapping = aes(fill=Group), outlier.colour = NA, size=0.2) +
  geom_text(data=tf_FL_WT_log2FC_info, aes(y=2.5, label=Label), size=1.2, color="black", vjust=0) +
  coord_cartesian(ylim = c(-4, 4)) +
  scale_fill_manual(values = c("Early"="red", "Late"="blue")) +
  facet_wrap(~Family, ncol = 6) +
  theme_bw() +
  labs(y="log2(FL expr. / WT expr.)") +
  theme(text = element_text(family="ArialMT", size=6),
        title = element_text(family="ArialMT", size=6),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(family="ArialMT", size=6),
        panel.grid = element_blank()
  )
ggsave(file.path(args$output, "ImpulseDE2_enrich_gene.FL_WT_expr_FC.each_TF_family.pdf"), p, width = 60, height = 40, units = "cm")

