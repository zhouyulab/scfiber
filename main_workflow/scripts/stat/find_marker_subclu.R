library(scater)
library(scran)
library(dplyr)
library(readr)
library(RColorBrewer)
library(gridExtra)
library(reshape2)

args = commandArgs(TRUE)
if (length(args) == 11) {
  f_mergeSCE <- args[1]
  f_corSCE <- args[2]
  f_cluster <- args[3]
  f_cnt <- args[4]
  f_tf <- args[5]
  fdr_cutoff <- as.numeric(args[6])
  expr_cutoff <- as.integer(args[7])
  bg_ratio <- as.numeric(args[8])
  ratio_fc <- as.numeric(args[9])
  mix_expr_ratio_curoff <- as.numeric(args[10])
  f_out <- args[11]
} else {
  q()
}

# f_mergeSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/mergeSCE.RData"
# f_corSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/corSCE.RData"
# f_cluster <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/cluster.tsv"
# f_cnt <- "/home/sqreb/ias/data/CottonSingleCell/analysis/cluster_cnt/data/cluster.cnt.tsv"
# fdr_cutoff <- 0.05
# expr_cutoff <- 2

load(f_mergeSCE)
load(f_corSCE)

clusterSCE <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)
mergeSCE$cluster <- clusterSCE$SubClu
corSCE$cluster <- clusterSCE$SubClu

cnt <- read_delim(f_cnt, "\t", escape_double = FALSE, trim_ws = TRUE)
cnt$FL.C3.1 <- NA
cnt$FL.C3.2 <- NA

melt_cnt <- melt(cnt, id.vars="Gid", variable.name="Group", value.name="TPM")
melt_cnt$Sample <-sapply(strsplit(as.character(melt_cnt$Group), "[.]"), function(x) x[1])
melt_cnt$Cluster <-sapply(strsplit(as.character(melt_cnt$Group), "[.]"), function(x) x[2])
melt_cnt$Group <- NULL

case_cnt <- dcast(melt_cnt, Gid+Sample~Cluster, value.var="TPM")
cnt_mat <- as.matrix(case_cnt[3:ncol(case_cnt)])
na_indx <- which(apply(cnt_mat, 2, function(x) any(is.na(x))))
case_cnt$MaxClu <- apply(cnt_mat, 1, which.max)

expr_info <- case_cnt %>% group_by(Gid) %>% summarise(Same=length(unique(MaxClu))==1, Miss=any(MaxClu %in% na_indx))
pass_gid <- expr_info$Gid[expr_info$Same | expr_info$Miss]

WT_indx <- mergeSCE$SampleName == "WT"
FL_indx <- (mergeSCE$SampleName == "FL") & (!mergeSCE$cluster %in% c("C3.1", "C3.2"))
markers <- findMarkers(mergeSCE, mergeSCE$cluster, pval.type="all", lfc=0, direction="up")
markers_WT <- findMarkers(mergeSCE[,WT_indx], mergeSCE[,WT_indx]$cluster, pval.type="all", lfc=0, direction="up")
markers_FL <- findMarkers(mergeSCE[,FL_indx], mergeSCE[,FL_indx]$cluster, pval.type="all", lfc=0, direction="up")

filter_markers <- function(block){
  logFC <- block[,4:ncol(block)]
  selected_id <- (block$FDR<fdr_cutoff)
  logFC_mat <- as.matrix(block[,3:ncol(block)])
  df <- data.frame(
    Gid=rownames(block)[selected_id], 
    FDR=block$FDR[selected_id],
    minLogFC=apply(logFC_mat, 1, min)[selected_id]
    )
  return(df)
}

selected_marker <- plyr::ldply(markers, filter_markers, .id="Cluster")
selected_marker_WT <- plyr::ldply(markers_WT, filter_markers, .id="Cluster")
selected_marker_FL <- plyr::ldply(markers_FL, filter_markers, .id="Cluster")
names(selected_marker)[3:4] <- c('FDR', 'MinLogFC')
names(selected_marker_WT)[3:4] <- c('FDR_WT', 'MinLogFC_WT')
names(selected_marker_FL)[3:4] <- c('FDR_FL', 'MinLogFC_FL')

selected_marker_gene_df <- inner_join(selected_marker, selected_marker_WT, by=c("Gid", "Cluster"))
selected_marker_gene_df <- left_join(selected_marker_gene_df, selected_marker_FL, by=c("Gid", "Cluster"))
selected_marker_gene_df <- selected_marker_gene_df[(selected_marker_gene_df$Cluster=="C3")| (!is.na(selected_marker_gene_df$MinLogFC_FL)),]
selected_marker_gene_df$FC_score <- apply(selected_marker_gene_df[,c("MinLogFC_WT", "MinLogFC_FL")], 1, function(x) min(x, na.rm = T))


norm_cnt_mat <- normcounts(mergeSCE)
norm_cnt_mat <- norm_cnt_mat[rownames(norm_cnt_mat)%in%selected_marker_gene_df$Gid,]
expr_mat <- norm_cnt_mat > expr_cutoff
rm(norm_cnt_mat)
gc()

expr_ratio_info <- data.frame(Gid=rownames(expr_mat))
sample_li <- sort(unique(mergeSCE$SampleName))
cluster_li <- sort(unique(mergeSCE$cluster))
for(sample in sample_li){
  for(cluster in cluster_li){
    if(sample=="FL" & cluster %in% c("C3.1", "C3.2")) next()
    part_expr_mat <- expr_mat[,mergeSCE$SampleName==sample & mergeSCE$cluster==cluster]
    expr_ratio_info[[sprintf("%s_%s_ExprRatio", sample, cluster)]] <- rowMeans(part_expr_mat)
  }
}

selected_marker_gene_df <- left_join(selected_marker_gene_df, expr_ratio_info, by="Gid")

stat_clu_ratio <- function(line){
  clu_ratio <- as.numeric(line[10:(length(line))])
  clu <- as.character(sapply(names(line)[10:(length(line))], function(x) strsplit(x, "_")[[1]][2]))
  sample_name <- as.character(sapply(names(line)[10:(length(line))], function(x) strsplit(x, "_")[[1]][1]))
  cluster <- as.character(line["Cluster"])
  cluster_indx <- (clu == cluster)
  return(data.frame(MinClusterRatio=min(clu_ratio[cluster_indx]), MaxNotClusterRatio=max(clu_ratio[!cluster_indx])))
}

summary_expr_ratio <- plyr::ddply(selected_marker_gene_df, "Gid", stat_clu_ratio)
selected_marker_gene_df <- left_join(selected_marker_gene_df, summary_expr_ratio, by="Gid")

TF_gene_list <- read_delim(f_tf, "\t", escape_double = FALSE, trim_ws = TRUE)

selected_marker_gene_df <- selected_marker_gene_df[selected_marker_gene_df$MaxNotClusterRatio<bg_ratio,]

x <- c("Ghir_A06G010490", "Ghir_D13G017770", "Ghir_D08G014350", "Ghir_A01G012170")
print(x)
print(x %in% TF_gene_list$Gid)
print(x %in% selected_marker_gene_df$Gid)
tf_marker_df <- selected_marker_gene_df[selected_marker_gene_df$Gid %in% TF_gene_list$Gid,]
# tf_marker_df <- tf_marker_df[tf_marker_df$MinClusterRatio>tf_marker_df$MaxNotClusterRatio,]

other_marker_df <- selected_marker_gene_df[!selected_marker_gene_df$Gid %in% TF_gene_list$Gid,]
other_marker_df <- other_marker_df[other_marker_df$MinClusterRatio > (ratio_fc * other_marker_df$MaxNotClusterRatio), ]
other_marker_df <- other_marker_df[other_marker_df$MinClusterRatio > mix_expr_ratio_curoff, ]

selected_marker_gene_df <- rbind(tf_marker_df, other_marker_df)
selected_marker_gene_df <- left_join(selected_marker_gene_df, cnt, by="Gid")
selected_marker_gene_df <- selected_marker_gene_df[order(selected_marker_gene_df$Cluster, -selected_marker_gene_df$FC_score),]
selected_marker_gene_df <- selected_marker_gene_df[selected_marker_gene_df$Gid %in% pass_gid,]
write_tsv(selected_marker_gene_df, file.path(f_out, "Marker.gene.tsv"))

## Cluster UMAP
uniq_clu <- sort(unique(clusterSCE$Cluster))
clu_num <- length(uniq_clu)

if(clu_num <= 8){
  cb_palette <- brewer.pal(8, "Dark2")
}else{
  cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                  "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                  "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                  "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                  "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                  "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
}


if(clu_num>length(cb_palette)) stop()

col_li <- cb_palette[1:clu_num]
names(col_li) <- uniq_clu

plot_cluster_umap  <- function(umapSCE, col_li, label){
  df <- data.frame(
    x=umapSCE@int_colData$reducedDims$UMAP[,1], 
    y=umapSCE@int_colData$reducedDims$UMAP[,2], 
    cluster=umapSCE$cluster
  )
  cell_num <- nrow(df)
  rand_indx <- sample(1:cell_num, cell_num)
  df <- df[rand_indx,]
  p <- ggplot(df, aes(x=x, y=y, color=cluster)) +
    theme_bw() +
    geom_point(size=0.1, alpha=0.2) +
    labs(x="UMAP1", y="UMAP2", alpha="UMI", title=label) +
    scale_color_manual(values = col_li) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
}

plot_gene_expr_umap <- function(umapSCE, SCE, gid, label){
  mat <- counts(SCE)
  expr <- mat[which(rownames(SCE)==gid),]
  if(! all(names(expr)==colnames(umapSCE))) stop()
  df <- data.frame(x=umapSCE@int_colData$reducedDims$UMAP[,1], y=umapSCE@int_colData$reducedDims$UMAP[,2], expr=expr)
  p <- ggplot(df, aes(x=x, y=y, alpha=expr)) +
    theme_bw() +
    geom_point(size=0.1, alpha=0.01, color="grey70") +
    geom_point(size=0.1, color="red") +
    labs(x="UMAP1", y="UMAP2", alpha="UMI", title=label) +
    scale_alpha_continuous(range = c(0, 0.5)) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
  return(p)
}

#ps_clu <- list()
#p_clu_all <- plot_cluster_umap(corSCE, col_li, "All")
#ps_clu[["All"]] <- p_clu_all

#samples <- sort(unique(corSCE$SampleName))
#for(s in samples){
#  tmp_SCE <- corSCE[,corSCE$SampleName==s]
#  ps_clu[[s]] <- plot_cluster_umap(tmp_SCE, col_li, s)
#}

#inche_cm=2.54

#for(clu in sort(unique(selected_marker_gene_df$Cluster))){
#  if(clu=="C3") next()
#  tmp_df <- selected_marker_gene_df[selected_marker_gene_df$Cluster==clu,]
#  UMAP_folder <- file.path(f_out, sprintf("UMAP_%s", clu))
#  if(!file.exists(UMAP_folder)){
#    dir.create(UMAP_folder, recursive = TRUE)
#  }
#  for(indx in 1:nrow(tmp_df)){
#    ps_gene_marker <- ps_clu
#    gid <- tmp_df$Gid[indx]
#    clu <- tmp_df$Cluster[indx]
#    label <- sprintf("%s: %s-%s", "All", clu, gid)
#    ps_gene_marker[[label]] <-  plot_gene_expr_umap(corSCE, mergeSCE, gid, label)
#    for(s in samples){
#      tmpUMAP <- corSCE[,corSCE$SampleName==s]
#      tmpSCE <- mergeSCE[,mergeSCE$SampleName==s]
#      label <- sprintf("%s: %s-%s", s, clu, gid)
#      ps_gene_marker[[label]] <-  plot_gene_expr_umap(tmpUMAP, tmpSCE, gid, label)
#      rm(tmpUMAP)
#      rm(tmpSCE)
#      gc()
#    }
#    gc()
#    pdf(file.path(UMAP_folder, sprintf("%s.UMAP.pdf", gid)), width=8*3/inche_cm, height=6*ceiling(length(ps_gene_marker)/3)/inche_cm, family="ArialMT", colormodel = "cmyk")
#    grid.arrange(grobs=ps_gene_marker,ncol=3,padding = unit(0, "mm"))
#    dev.off()
#  }
  
#}

#selected_marker_gene_strict_df <- selected_marker_gene_df[selected_marker_gene_df$FC_score>1.5,]
#for(clu in sort(unique(selected_marker_gene_strict_df$Cluster))){
#  tmp_df <- selected_marker_gene_strict_df[selected_marker_gene_strict_df$Cluster==clu,]
#  UMAP_folder <- file.path(f_out, sprintf("UMAP_strict_%s", clu))
#  if(!file.exists(UMAP_folder)){
#    dir.create(UMAP_folder, recursive = TRUE)
#  }
#  for(indx in 1:nrow(tmp_df)){
#    ps_gene_marker <- ps_clu
#    gid <- tmp_df$Gid[indx]
#    clu <- tmp_df$Cluster[indx]
#    label <- sprintf("%s: %s-%s", "All", clu, gid)
#    ps_gene_marker[[label]] <-  plot_gene_expr_umap(corSCE, mergeSCE, gid, label)
#    for(s in samples){
#      tmpUMAP <- corSCE[,corSCE$SampleName==s]
#      tmpSCE <- mergeSCE[,mergeSCE$SampleName==s]
#      label <- sprintf("%s: %s-%s", s, clu, gid)
#      ps_gene_marker[[label]] <-  plot_gene_expr_umap(tmpUMAP, tmpSCE, gid, label)
#      rm(tmpUMAP)
#      rm(tmpSCE)
#      gc()
#    }
#    gc()
#    pdf(file.path(UMAP_folder, sprintf("%s.UMAP.pdf", gid)), width=8*3/inche_cm, height=6*ceiling(length(ps_gene_marker)/3)/inche_cm, family="ArialMT", colormodel = "cmyk")
#    grid.arrange(grobs=ps_gene_marker,ncol=3,padding = unit(0, "mm"))
#    dev.off()
#  }
#}
