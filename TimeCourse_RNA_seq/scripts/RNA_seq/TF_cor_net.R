library(argparse)
library(dplyr)
library(readr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(igraph)
library(BiocParallel)
library(circlize)
parser <- ArgumentParser()
parser$add_argument("--tp-input", required=TRUE, type="character", dest = "tp_input", metavar="tp_fc.tsv")
parser$add_argument("--other-input", nargs="*", required=TRUE, type="character", dest = "other_input", metavar="other_fc.tsv")
parser$add_argument("--merge-cpm", required=TRUE, type="character", dest = "merge_cpm", metavar="merge_cpm.tsv")
parser$add_argument("--tf", required=TRUE, type="character", dest = "tf", metavar="tf.tsv")
parser$add_argument("--expr-cpm-cutoff", required=FALSE, type="double", dest = "expr_cpm_cutoff", metavar="expr_cpm_cutoff", default=0.2)
parser$add_argument("--interest-gene", required=TRUE, type="character", dest = "interest_gene", metavar="interest_gene.txt")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")


# args <- commandArgs(TRUE)
# args <- parser$parse_args(args)

# setwd("~/mu01/project/CottonSCE_TP")
args <- c(
  "--tp-input",
  "analysis/RNA_seq_TP/merge_expr/feature_counts.merge.tsv",
  "--other-input",
  "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA2.merge.tsv",
  "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA4_part1.merge.tsv",
  "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA5.merge.tsv",
  "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA_web_part1.merge.tsv",
  "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA_web_part2.merge.tsv",
  "--merge-cpm",
  "analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv",
  "--tf", "data/TF.gene_list.tsv",
  "--interest-gene",
  "data/InterestGeneID.tsv",
  "-o",
  "analysis/RNA_seq_TP/TfCorNet"
)
args <- parser$parse_args(args)

fc_tp_df <- read_delim(args$tp_input, "\t", escape_double = FALSE, trim_ws = TRUE)
samples <- unique(sapply(strsplit(names(fc_tp_df)[2:ncol(fc_tp_df)], "[.]") , function(x){return(x[1])}))
tps <- unique(sapply(strsplit(names(fc_tp_df)[2:ncol(fc_tp_df)], "[.]") , function(x){return(x[2])}))
reps <- unique(sapply(strsplit(names(fc_tp_df)[2:ncol(fc_tp_df)], "[.]") , function(x){return(x[3])}))

fc2cpm <- function(x){
  fc <- as.matrix(x[,2:ncol(x)])
  cpm <- 1e6 * t(t(fc) / colSums(fc))
  x[,2:ncol(x)] <- cpm
  return(x)
}
cpm_tp_df <- fc2cpm(fc_tp_df)

tf_gene <- read_delim(args$tf, "\t", escape_double = FALSE, trim_ws = TRUE)

tp_tf_cpm_df <- left_join(data.frame(Gid=unique(tf_gene$Gid)), cpm_tp_df)
all_tf_cpm_df <- tp_tf_cpm_df
for(fname in args$other_input){
  tmp_fc_df <- read_delim(fname, "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp_cpm_df <- fc2cpm(tmp_fc_df)
  all_tf_cpm_df <- left_join(all_tf_cpm_df, tmp_cpm_df)
}
tp_tf_cpm_mat <- as.matrix(tp_tf_cpm_df[,2:ncol(tp_tf_cpm_df)])
rownames(tp_tf_cpm_mat) <- tp_tf_cpm_df$Gid
all_tf_cpm_mat <- as.matrix(all_tf_cpm_df[,2:ncol(all_tf_cpm_df)])
rownames(all_tf_cpm_mat) <- all_tf_cpm_df$Gid

expr_tp_tf_cpm_mat <- tp_tf_cpm_mat[rowMeans(tp_tf_cpm_mat)>args$expr_cpm_cutoff,]
expr_all_tf_cpm_mat <- all_tf_cpm_mat[rowMeans(all_tf_cpm_mat)>args$expr_cpm_cutoff,]

filter_de_tf <- function(cpm_mat, DE_ratio_cutoff=0.5, cv_cutoff=0.5){
  max_expr <- apply(cpm_mat, 1, max)
  min_expr <- apply(cpm_mat, 1, min)
  diff_expr <- max_expr - min_expr
  diff_expr_ratio <- diff_expr / max_expr
  cv <- apply(cpm_mat, 1, sd) / apply(cpm_mat, 1, mean)
  return(cpm_mat[diff_expr_ratio>DE_ratio_cutoff & cv>cv_cutoff,])
}

expr_DE_tp_tf_cpm_mat <- filter_de_tf(expr_tp_tf_cpm_mat)
expr_DE_all_tf_cpm_mat <- filter_de_tf(expr_all_tf_cpm_mat)

expr_tp_tf_pearson_cor_mat <- cor(t(log10(expr_DE_tp_tf_cpm_mat+1e-4)), method = "pearson")
expr_tp_tf_spearman_cor_mat <- cor(t(log10(expr_DE_tp_tf_cpm_mat+1e-4)), method = "spearman")
expr_tp_all_pearson_cor_mat <- cor(t(log10(expr_DE_all_tf_cpm_mat+1e-4)), method = "pearson")
expr_tp_all_spearman_cor_mat <- cor(t(log10(expr_DE_all_tf_cpm_mat+1e-4)), method = "spearman")
# cor_mat_li <- list(
#   "TP_pearson" = expr_tp_tf_pearson_cor_mat,
#   "TP_spearman" = expr_tp_tf_spearman_cor_mat,
#   "all_pearson" = expr_tp_all_pearson_cor_mat,
#   "all_spearman" = expr_tp_all_spearman_cor_mat
# )

cor_mat_li <- list(
  "TP_pearson" = expr_tp_tf_pearson_cor_mat,
  "TP_spearman" = expr_tp_tf_spearman_cor_mat
)

cluster_tf_by_hclust <- function(cor_mat){
  clu <- hclust(as.dist(1-cor_mat), method = "ward.D")
  return(clu)
}
hclust_clu_li <- bplapply(cor_mat_li, cluster_tf_by_hclust, BPPARAM=MulticoreParam(length(cor_mat_li)))


filter_tf_clu_by_num <- function(clu_df, min_num=5){
  clu_info <- clu_df %>% group_by(Cluster) %>% summarise(GeneNum=n())
  enriched_clu_info <- clu_info[clu_info$GeneNum >= min_num,]
  # enriched_clu_info <- enriched_clu_info[-grep("N0.7_", enriched_clu_info$Cluster),]
  enriched_clu <- clu_df[clu_df$Cluster %in% enriched_clu_info$Cluster,]
  return(enriched_clu)
}

compute_ave_cor <- function(cor_mat, clu_df){
  all_clu <- sort(unique(clu_df$Cluster))
  clu_num <- length(all_clu)
  ave_clu_cor_mat <- matrix(0, nrow = clu_num, ncol = clu_num)
  rownames(ave_clu_cor_mat) <- all_clu
  colnames(ave_clu_cor_mat) <- all_clu
  for(indx1 in 1:clu_num){
    for(indx2 in 1:indx1){
      gid1 <- as.character(clu_df$Gid[clu_df$Cluster==all_clu[indx1]])
      gid2 <- as.character(clu_df$Gid[clu_df$Cluster==all_clu[indx2]])
      ave_cor <- mean(cor_mat[gid1, gid2])
      ave_clu_cor_mat[indx1, indx2] <- ave_cor
      ave_clu_cor_mat[indx2, indx1] <- ave_cor
    }
  }
  return(ave_clu_cor_mat)
}


cor_mat2clu <- function(cor_mat, cor_cutoff=0.7){
  high_cor_mat <- cor_mat>0.9
  high_cor_graph <- graph_from_adjacency_matrix(high_cor_mat, mode="upper")
  high_cor_clique <- max_cliques(high_cor_graph, min=5)
  cliqie_graph <- induced_subgraph(high_cor_graph, unique(names(unlist(high_cor_clique))))
  components_li <- groups(components(cliqie_graph))
  group_mat <- matrix(0, nrow = nrow(high_cor_mat), ncol = ncol(high_cor_mat))
  rownames(group_mat) <- rownames(high_cor_mat)
  colnames(group_mat) <- colnames(high_cor_mat)
  for(component_indx in 1:length(components_li)){
    group_mat[rownames(group_mat)%in%components_li[[component_indx]], rownames(group_mat)%in%components_li[[component_indx]]] <- TRUE
  }
  group_mat[cor_mat<cor_cutoff] <- FALSE
  for(indx in 1:nrow(cor_mat)){
    gid <- colnames(cor_mat)[indx]
    for(component_indx in 1:length(components_li)){
      component <- components_li[[component_indx]]
      if(gid %in% component) next()
      tmp_cor <- cor_mat[gid, component]
      high_cor_line <- tmp_cor>cor_cutoff
      
      if(all(high_cor_line)){
        group_mat[rownames(group_mat)==gid, rownames(group_mat)%in%component[high_cor_line]] <- T
        group_mat[rownames(group_mat)%in%component[high_cor_line], rownames(group_mat)==gid] <- T
      }
    }
  }
  group_mat[cor_mat<cor_cutoff] <- FALSE
  group_graph <- graph_from_adjacency_matrix(group_mat, mode="upper", weighted=TRUE)
  group_clu <- groups(igraph::components(group_graph))
  group_clu_df <- data.frame()
  for(indx in 1:length(group_clu)){
    if(length(group_clu[[indx]])>1){
      group_clu_df <- rbind(group_clu_df, data.frame(Gid=group_clu[[indx]], Cluster=sprintf("C%s", indx)))
    }
  }
  
  for(tmp_cutoff in seq(0.9, cor_cutoff, -0.05)){
    print(tmp_cutoff)
    no_group_mat <- matrix(0, nrow = nrow(group_mat), ncol = ncol(group_mat))
    rownames(no_group_mat) <- rownames(group_mat)
    colnames(no_group_mat) <- colnames(group_mat)
    
    no_group_indx <- !rownames(group_mat) %in% group_clu_df$Gid
    no_group_nodes <- rownames(group_mat)[no_group_indx]
    
    rest_cor_mat <- cor_mat[no_group_indx, no_group_indx]
    rest_high_cor_mat <- rest_cor_mat > tmp_cutoff
    
    rest_high_cor_graph <- graph_from_adjacency_matrix(rest_high_cor_mat, mode="upper")
    rest_high_cor_clique <- max_cliques(rest_high_cor_graph, min=5)
    if(length(rest_high_cor_clique)==0) next()
    name_li <- rownames(no_group_mat)
    for(clique_indx in 1:length(rest_high_cor_clique)){
      clique <- names(unlist(rest_high_cor_clique[[clique_indx]]))
      no_group_mat[name_li%in%clique, name_li%in%clique] <- TRUE
      for(gid in no_group_nodes){
        if(all(cor_mat[name_li==gid, name_li%in%clique]>cor_cutoff)){
          no_group_mat[name_li==gid, name_li%in%clique] <- TRUE
          no_group_mat[name_li%in%clique, name_li==gid] <- TRUE
        }
      }
    }
    tmp_group_cor_mat <- cor_mat
    tmp_group_cor_mat[!no_group_mat] <- 0
    tmp_group_cor_mat <- tmp_group_cor_mat[rowSums(tmp_group_cor_mat)>0, rowSums(tmp_group_cor_mat)>0]
    tmp_graph <- graph_from_adjacency_matrix(tmp_group_cor_mat, mode="upper", weighted=TRUE)
    components_li <- groups(components(tmp_graph))
    tmp_clu_df <- data.frame()
    for(component_indx in 1:length(components_li)){
      print(component_indx)
      sub_graph <- induced_subgraph(tmp_graph, components_li[[component_indx]])
      tmp_adj_mat <- as_adjacency_matrix(sub_graph)
      tmp_clu <- MCL::mcl(tmp_adj_mat, addLoops = FALSE)
      tmp_clu_df <- rbind(tmp_clu_df, data.frame(Gid=rownames(tmp_adj_mat), Cluster=sprintf("N%s_%s_%s", tmp_cutoff, component_indx, tmp_clu$Cluster)))
    }
    
    clu_num <- tmp_clu_df %>% group_by(Cluster) %>% summarise(GeneNum=n())
    tmp_clu_df <- tmp_clu_df[tmp_clu_df$Cluster%in%clu_num$Cluster[clu_num$GeneNum>5],]
    group_clu_df <- rbind(group_clu_df, tmp_clu_df)
    no_group_mat[!rownames(no_group_mat) %in% tmp_clu_df$Gid, !rownames(no_group_mat) %in% tmp_clu_df$Gid] <- 0
    group_mat <- group_mat + no_group_mat
  }

  return(group_clu_df)
}

# raw_clu_li <- bplapply(cor_mat_li, cor_mat2clu, BPPARAM=MulticoreParam(length(cor_mat_li)))
raw_clu_li <- lapply(cor_mat_li, cor_mat2clu)
# save(raw_clu_li, file = file.path(args$output, "raw_clu.RData"))

clu_li <- bplapply(raw_clu_li, filter_tf_clu_by_num, min_num=5, BPPARAM=MulticoreParam(length(raw_clu_li)))
clu_cor_li <- list()
for(s in names(clu_li)){
  clu_cor_li[[s]] <- compute_ave_cor(cor_mat_li[[s]], clu_li[[s]])
}

merge_clu <- function(clu_df, clu_cor_mat, max_diff=0.15, ave_diff=0.05){
  from_li <- c()
  to_li <- c()
  for(indx in 1:ncol(clu_cor_mat)){
    cond1 <- apply(clu_cor_mat - clu_cor_mat[indx,], 2, max) < max_diff
    cond2 <- abs(apply(clu_cor_mat - clu_cor_mat[indx,], 2, mean)) < ave_diff
    cond <- cond1 & cond2
    cond[indx] <- FALSE
    if(any(cond)){
      from_name <- rownames(clu_cor_mat)[indx]
      to_name_li <- rownames(clu_cor_mat)[cond]
      from_li <- c(from_li, rep(from_name, length(to_name_li)))
      to_li <- c(to_li, to_name_li)
    }
  }
  
  if(length(from_li)>0){
    clu_df$Cluster <- as.character(clu_df$Cluster)
    edge_df <- data.frame(from=from_li, to=to_li)
    graph <- igraph::graph_from_data_frame(edge_df)
    group_clu <- igraph::components(graph)$membership
    for(indx in unique(group_clu)){
      clu_df$Cluster[clu_df$Cluster %in% names(group_clu)[indx]] <- sprintf("M%s", indx)
    }
  }
  return(clu_df)
}

merged_clu_li <- list()
for(s in names(clu_li)){
  merged_clu_li[[s]] <- merge_clu(clu_li[[s]], clu_cor_li[[s]])
}

ave_clu_cor_li <- list()
for(s in names(merged_clu_li)){
  ave_clu_cor_li[[s]] <- compute_ave_cor(cor_mat_li[[s]], merged_clu_li[[s]])
}
ave_clu_cor_clu_li <- bplapply(ave_clu_cor_li, function(x){
  return(hclust(dist(x), method = "ward.D"))
}, BPPARAM=MulticoreParam(length(ave_clu_cor_li)))


renamed_clu_li <- list()
for(s in names(merged_clu_li)){
  ordered_df <- data.frame(Cluster=ave_clu_cor_clu_li[[s]]$labels[ave_clu_cor_clu_li[[s]]$order])
  ordered_df$NewCluster <- sprintf("G%s", 1:nrow(ordered_df))
  ordered_df$NewCluster <- factor(ordered_df$NewCluster, levels = ordered_df$NewCluster)
  ave_clu_cor_li[[s]] <- ave_clu_cor_li[[s]][ave_clu_cor_clu_li[[s]]$order, ave_clu_cor_clu_li[[s]]$order]
  colnames(ave_clu_cor_li[[s]]) <- ordered_df$NewCluster
  rownames(ave_clu_cor_li[[s]]) <- ordered_df$NewCluster
  renamed_clu_df <- left_join(merged_clu_li[[s]], ordered_df)
  renamed_clu_df <- renamed_clu_df[,c("Gid", "NewCluster")]
  names(renamed_clu_df) <- c("Gid", "Cluster")
  renamed_clu_li[[s]] <- renamed_clu_df
}



compute_connection_ratio <- function(cor_mat, clu_df, cor_cutoff){
  all_clu <- sort(unique(clu_df$Cluster))
  clu_num <- length(all_clu)
  connection_ratio_mat <- matrix(0, nrow = clu_num, ncol = clu_num)
  rownames(connection_ratio_mat) <- all_clu
  colnames(connection_ratio_mat) <- all_clu
  for(indx1 in 1:clu_num){
    for(indx2 in 1:indx1){
      gid1 <- clu_df$Gid[clu_df$Cluster==all_clu[indx1]]
      gid2 <- clu_df$Gid[clu_df$Cluster==all_clu[indx2]]
      connection_ratio <- mean(cor_mat[rownames(cor_mat) %in% gid1, rownames(cor_mat) %in% gid2]>cor_cutoff)
      connection_ratio_mat[indx1, indx2] <- connection_ratio
      connection_ratio_mat[indx2, indx1] <- connection_ratio
    }
  }
  return(connection_ratio_mat)
}

interest_gene <- read_delim(args$interest_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
interest_gene <- interest_gene[!duplicated(interest_gene[,c("ID", "Gid")]),]
interest_gene <- left_join(interest_gene, tf_gene)
merge_CPM <- read_delim(args$merge_cpm, "\t", escape_double = FALSE, trim_ws = TRUE)
tp_num_df <- data.frame(
  TP=c("n48h", "n36h", "n24h", "n12h", "0h", "12h", "24h", "36h", "48h"),
  TP_Num=c(-48, -36, -24, -12, 0, 12, 24, 36, 48)
  )

for(s in names(renamed_clu_li)){
  tmp_dir <- file.path(args$output, s)
  if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir)
  }
  orderd_cor_mat <- cor_mat_li[[s]][hclust_clu_li[[s]]$order, hclust_clu_li[[s]]$order]
  # group_df <- left_join(data.frame(Gid=rownames(orderd_cor_mat)), raw_clu_li[[s]])
  group_df <- left_join(data.frame(Gid=rownames(orderd_cor_mat)), renamed_clu_li[[s]])
  group_df <- group_df[!is.na(group_df$Cluster),]
  orderd_cor_mat <- orderd_cor_mat[rownames(orderd_cor_mat) %in% group_df$Gid, rownames(orderd_cor_mat) %in% group_df$Gid]
  write_tsv(group_df, file.path(tmp_dir, "TF.CorNetGroup.tsv"))
  ave_clu_cor_mat <- ave_clu_cor_li[[s]]
  connection_ratio_mat_0.5 <- compute_connection_ratio(cor_mat_li[[s]], renamed_clu_li[[s]], cor_cutoff=0.5)
  connection_ratio_mat_0.75 <- compute_connection_ratio(cor_mat_li[[s]], renamed_clu_li[[s]], cor_cutoff=0.75)
  connection_ratio_mat_0.9 <- compute_connection_ratio(cor_mat_li[[s]], renamed_clu_li[[s]], cor_cutoff=0.9)
  group1_li <- c()
  group2_li <- c()
  group_cor_li <- c()
  group_connection_ratio_0.5_li <- c()
  group_connection_ratio_0.75_li <- c()
  group_connection_ratio_0.9_li <- c()
  for(indx1 in 1:nrow(ave_clu_cor_mat)){
    for(indx2 in 1:indx1){
      group1_li <- c(group1_li, rownames(ave_clu_cor_mat)[indx1])
      group2_li <- c(group2_li, rownames(ave_clu_cor_mat)[indx2])
      group_cor_li <- c(group_cor_li, ave_clu_cor_mat[indx1, indx2])
      group_connection_ratio_0.5_li <- c(group_connection_ratio_0.5_li, connection_ratio_mat_0.5[indx1, indx2])
      group_connection_ratio_0.75_li <- c(group_connection_ratio_0.75_li, connection_ratio_mat_0.75[indx1, indx2])
      group_connection_ratio_0.9_li <- c(group_connection_ratio_0.9_li, connection_ratio_mat_0.9[indx1, indx2])
    }
  }
  group_connection_df <- data.frame(
    Group1=group1_li, Group2=group2_li, AveCor=group_cor_li, 
    ConnectionRatio0.5=group_connection_ratio_0.5_li, ConnectionRatio0.75=group_connection_ratio_0.75_li, ConnectionRatio0.9=group_connection_ratio_0.9_li
    )
  write_tsv(group_connection_df, file.path(tmp_dir, "TF.Group.Connection.stat.tsv"))
  
  interest_gene_df <- left_join(group_df, interest_gene)
  interest_gene_df$Label <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$GeneSymbol)
  interest_gene_df$Label[!is.na(interest_gene_df$Family)] <- sprintf("%s (%s: %s)", interest_gene_df$Gid, interest_gene_df$Family, interest_gene_df$GeneSymbol)[!is.na(interest_gene_df$Family)]
  interest_gene_df$Label[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$Family)[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]
  interest_gene_df$Label[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- interest_gene_df$Gid[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]
  
  subindx <- which(!is.na(interest_gene_df$GeneSymbol))
  label_li <- interest_gene_df$Label[subindx]
  ht <- Heatmap(
    orderd_cor_mat, 
    name="Cor",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    row_names_gp = gpar(fontsize = 6),
    row_title_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    column_split = group_df$Cluster, 
    row_split = group_df$Cluster, 
    show_row_names = F, 
    show_column_names = F,
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  
  inche_cm <- 2.54
  if(length(label_li)>0){
    ha <- rowAnnotation(
      link = row_anno_link(
        at = subindx,
        labels = label_li,
        labels_gp=gpar(fontsize = 6),
        lines_gp=gpar(lwd = 0.5),
        link_width = unit(4, "mm")),
      width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)))
    
    pdf(file.path(tmp_dir, "TF.gene.cor.pdf"), width=13.5/inche_cm, height=8/inche_cm)
    print(ht+ha)
    dev.off()
  }else{
    pdf(file.path(tmp_dir, "TF.gene.cor.pdf"), width=9.5/inche_cm, height=8/inche_cm)
    print(ht)
    dev.off()
  }
  
  ht <- Heatmap(
    ave_clu_cor_mat, 
    name="Cor",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    row_names_gp = gpar(fontsize = 6),
    row_title_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  pdf(file.path(tmp_dir, "TF.group.cor.pdf"), width=9.5/inche_cm, height=8/inche_cm)
  print(ht)
  dev.off()
  
  ht <- Heatmap(
    connection_ratio_mat_0.5, 
    name="Connection",
    col = colorRamp2(c(0, 1), c("white", "red")),
    row_names_gp = gpar(fontsize = 6),
    row_title_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  pdf(file.path(tmp_dir, "TF.group.connection_0.5.pdf"), width=9.5/inche_cm, height=8/inche_cm)
  print(ht)
  dev.off()
  
  ht <- Heatmap(
    connection_ratio_mat_0.75, 
    name="Connection",
    col = colorRamp2(c(0, 1), c("white", "red")),
    row_names_gp = gpar(fontsize = 6),
    row_title_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  pdf(file.path(tmp_dir, "TF.group.connection_0.75.pdf"), width=9.5/inche_cm, height=8/inche_cm)
  print(ht)
  dev.off()
  
  select_cpm <- left_join(group_df, cpm_tp_df)
  select_cpm_mat <- as.matrix(select_cpm[,3:(ncol(select_cpm)-1)])
  norm_select_cpm_mat <- select_cpm_mat / apply(select_cpm_mat, 1, max)
  norm_select_cpm_mat[is.nan(norm_select_cpm_mat)] <- 0 
  ht <- Heatmap(
    norm_select_cpm_mat, 
    name="NormExpr.",
    col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
    row_names_gp = gpar(fontsize = 6),
    row_title_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    # cluster_rows = FALSE, 
    cluster_columns = FALSE,
    split=select_cpm$Cluster,
    column_split = sapply(strsplit(colnames(norm_select_cpm_mat), "[.]"), function(x){return(x[1])}),
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  if(length(label_li)>0){
    ha <- rowAnnotation(
      link = row_anno_link(
        at = subindx,
        labels = label_li,
        labels_gp=gpar(fontsize = 6),
        lines_gp=gpar(lwd = 0.5),
        link_width = unit(4, "mm")),
      width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)))
    
    pdf(file.path(tmp_dir, "TF.gene.expr.pdf"), width=13.5/inche_cm, height=15/inche_cm)
    print(ht+ha)
    dev.off()
  }else{
    pdf(file.path(tmp_dir, "TF.gene.expr.pdf"), width=9.5/inche_cm, height=15/inche_cm)
    print(ht)
    dev.off()
  }
  
  group_CPM <- left_join(group_df, merge_CPM)
  group_norm_expr <- group_CPM
  group_norm_expr[,3:ncol(group_norm_expr)] <- group_norm_expr[,3:ncol(group_norm_expr)] / apply(group_norm_expr[,3:ncol(group_norm_expr)], 1, max)
  melt_group_norm_expr <- melt(group_norm_expr, id.vars = c("Gid", "Cluster"), variable.name = "Label", value.name = "NormExpr")
  melt_group_norm_expr$Sample <- sapply(strsplit(as.character(melt_group_norm_expr$Label), "_"), function(x){return(x[1])})
  melt_group_norm_expr$TP <- sapply(strsplit(as.character(melt_group_norm_expr$Label), "_"), function(x){return(x[2])})
  melt_group_norm_expr <- left_join(melt_group_norm_expr, tp_num_df)
  melt_group_norm_expr$Sample <- factor(melt_group_norm_expr$Sample, levels = c("WT", "FL"))
  p <- ggplot(melt_group_norm_expr, aes(x=TP_Num, y=NormExpr, group=Gid, color=Sample)) +
    geom_line(size=0.1, alpha=0.3) +
    facet_grid(Cluster~Sample) +
    scale_x_continuous(breaks = tp_num_df$TP_Num, labels = tp_num_df$TP) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    scale_color_brewer(palette = "Set1") +
    labs(y="Norm. Expr.") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", size=5),
          title = element_text(family="ArialMT", size=5),
          axis.title.x = element_blank(),
          axis.text = element_text(color = "black"),
          legend.position = "none",
          panel.grid = element_blank()
    )
  ggsave(file.path(tmp_dir, "TF.gene.expr.group.curve.pdf"), p, width = 6, height = 16, units = "cm")
  
  melt_group_norm_expr_info <- melt_group_norm_expr %>% group_by(Cluster, TP_Num, Sample) %>% summarise(Mu=mean(NormExpr), SD=sd(NormExpr))
  melt_group_norm_expr_info$ymin <- melt_group_norm_expr_info$Mu - melt_group_norm_expr_info$SD
  melt_group_norm_expr_info$ymax <- melt_group_norm_expr_info$Mu + melt_group_norm_expr_info$SD
  p <- ggplot(melt_group_norm_expr_info, aes(x=TP_Num, y=Mu, ymin=ymin, ymax=ymax, color=Sample)) +
    geom_smooth(stat="identity", size=0.2) +
    facet_wrap(~Cluster) +
    coord_cartesian(ylim=c(0, 1)) +
    scale_x_continuous(breaks = tp_num_df$TP_Num, labels = tp_num_df$TP) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_color_brewer(palette = "Set1") +
    labs(y="Norm. Expr.") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", size=5),
          title = element_text(family="ArialMT", size=5),
          axis.title.x = element_blank(),
          axis.text = element_text(color = "black"),
          legend.key.size = unit(4, "mm"),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(family="ArialMT", size=6),
          panel.grid = element_blank()
    )
  ggsave(file.path(tmp_dir, "TF.gene.expr.group.merge.curve.pdf"), p, width = 15, height = 12, units = "cm")
}
