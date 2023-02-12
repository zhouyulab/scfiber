library(readr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(igraph)
library(ggtree)
library(BiocParallel)

args = commandArgs(TRUE)
if (length(args) == 8) {
  f_in <- args[1]
  f_mergeSCE <- args[2]
  f_corSCE <- args[3]
  min_expr_cutoff <- as.numeric(args[4])
  k_nearest <- as.integer(args[5])
  tree_cutoff <- as.numeric(args[6])
  min_group_gene <- as.integer(args[7])
  f_out <- args[8]
} else {
  q()
}

# f_in <- "/home/sqreb/ias/data/CottonSingleCell/analysis_v2/recluster/modularity/Modularity.filter.tsv"
# min_expr_cutoff <- 2
# tree_cutoff <- 0.9
# min_group_gene <- 5
# k_nearest <- 5

load(f_mergeSCE)
marker_df <- read_delim(f_in, "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(mergeSCE)[1:3]
print(marker_df$Gid)[1:3]
print(sum(rownames(mergeSCE) %in% marker_df$Gid))
mergeSCE <- mergeSCE[rownames(mergeSCE) %in% marker_df$Gid,]
gc()

load(f_corSCE)


# find_similar_cell <- function(cell_indx, SNN, k_nearest, inner_ratio_cutoff=0.3){
#   one_order_neighhor_cell <- unique(unlist(ego(SNN, order = k_nearest, cell_indx)))
#   two_order_neighbor_cell <- ego(SNN, order = k_nearest, one_order_neighhor_cell)
#   inner_ratio <- sapply(two_order_neighbor_cell, function(x) mean(x %in% cell_indx))
#   inner_cell <- one_order_neighhor_cell[inner_ratio>inner_ratio_cutoff]
#   return(inner_cell)
# }
# 
# compute_bool_matrix_dist <- function(mat){
#   mat <- mat[,colSums(mat) > 0]
#   gene_num <- nrow(mat)
#   dist_mat <- matrix(0, nrow = gene_num, ncol = gene_num)
#   rownames(dist_mat) <- rownames(mat)
#   colnames(dist_mat) <- rownames(mat)
#   for(i in 2:gene_num){
#     for(j in 1:i){
#       bool1 <- mat[i,]
#       bool2 <- mat[j,]
#       all_bool <- (bool1 + bool2) >0
#       tmp_dist <- sum(bool1 != bool2) / sum(all_bool)
#       dist_mat[i, j] <- tmp_dist
#       dist_mat[j, i] <- tmp_dist
#     }
#   }
#   return(dist_mat)
# }

compute_gene_distance <- function(UMAP_mat, cell_indx_li1, cell_indx_li2, k_nearest=5, drop_ratio=0.05){
  dim1 <- length(cell_indx_li1)
  dim2 <- length(cell_indx_li2)
  
  mat1 <- t(UMAP_mat[cell_indx_li1,])
  mat2 <- t(UMAP_mat[cell_indx_li2,])
  
  mat1_dist <- matrix(diag(as.matrix(crossprod(mat1, mat1))), nrow=dim1, ncol=dim2, byrow=FALSE)
  mat2_dist <- matrix(diag(as.matrix(crossprod(mat2, mat2))), nrow=dim1, ncol=dim2, byrow=TRUE)
  cross_dist <- as.matrix(t(mat1) %*% mat2)
  
  dist_mat <- sqrt(mat1_dist + mat2_dist - 2 * cross_dist)
  
  dim1_dist <- apply(dist_mat, 1, function(x){mean(sort(x)[1:k_nearest])})
  dim2_dist <- apply(dist_mat, 2, function(x){mean(sort(x)[1:k_nearest])})
  
  dim1_retain <- ceiling((1-drop_ratio) * length(dim1_dist))
  dim2_retain <- ceiling((1-drop_ratio) * length(dim2_dist))
  dim1_dist <- sort(dim1_dist)[1:dim1_retain]
  dim2_dist <- sort(dim2_dist)[1:dim2_retain]
  
  return(mean(c(dim1_dist, dim2_dist)))
}

norm_mat <- normcounts(mergeSCE)
expr_mat <- norm_mat > min_expr_cutoff

UMAP_mat <- corSCE@int_colData$reducedDims$UMAP

compute_dist_mat <- function(pos, UMAP_mat, expr_mat, k_nearest=5){
  i <- pos[1]
  j <- pos[2]
  cell_indx_li1 <- which(expr_mat[i,])
  cell_indx_li2 <- which(expr_mat[j,])
  tmp_dist <- compute_gene_distance(UMAP_mat, cell_indx_li1, cell_indx_li2, k_nearest=5)
  return(tmp_dist)
}


gene_num <- nrow(expr_mat)
index_pair <- list()
for(i in 2:gene_num){
  for(j in 1:i){
    index_pair[[paste(i, j)]] <- c(i, j)
  }
}

res <- lapply(index_pair, compute_dist_mat, UMAP_mat=UMAP_mat, expr_mat=expr_mat, k_nearest=k_nearest)

res <- bplapply(index_pair, compute_dist_mat, UMAP_mat=UMAP_mat, expr_mat=expr_mat, k_nearest=k_nearest, BPPARAM=MulticoreParam(workers=30))
gene_dist_mat <- matrix(0, nrow = gene_num, ncol = gene_num)
colnames(gene_dist_mat) <- rownames(expr_mat)
rownames(gene_dist_mat) <- rownames(expr_mat)
for(key in names(res)){
  pos <- as.integer(strsplit(key, " ")[[1]])
  gene_dist_mat[pos[1], pos[2]] <- res[[key]]
  gene_dist_mat[pos[2], pos[1]] <- res[[key]]
}

gene_clu <- hclust(as.dist(gene_dist_mat))
plot(gene_clu)

gene_group <- data.frame(Gid=gene_clu$labels, Group=cutree(gene_clu, h = tree_cutoff))
gene_group_info <- as.data.frame(table(gene_group$Group))
gene_group$Group[gene_group$Group%in%gene_group_info$Var1[gene_group_info$Freq<min_group_gene]] <- NA 
gene_group_info <- as.data.frame(table(gene_group$Group))
gene_group_info <- gene_group_info[order(gene_group_info$Freq, decreasing = TRUE),]
names(gene_group_info)[1] <- "Group"
gene_group_info$Group <- as.integer(as.character(gene_group_info$Group))
gene_group_info$Label <- sprintf("G%d", 1:nrow(gene_group_info))
gene_group <- left_join(gene_group, gene_group_info, by="Group")
gene_group$Label[is.na(gene_group$Label)] <- "Outer"
if("Outer" %in% gene_group$Label){
  gene_group$Label <- factor(gene_group$Label, levels = c(gene_group_info$Label, "Outer"))
}else{
  gene_group$Label <- factor(gene_group$Label, levels = c(gene_group_info$Label))
}

gene_group$Gid <- as.character(gene_group$Gid)

clu_num <- nrow(gene_group_info)
if(clu_num <= 12){
  cb_palette <- brewer.pal(12, "Paired")
}else{
  cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                  "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                  "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                  "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                  "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                  "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
}

if("Outer" %in% gene_group$Label){
  col_li <- c(cb_palette[1:clu_num], "black", "black")
  names(col_li) <- c(gene_group_info$Label, "Outer", "0")
}else{
  col_li <- c(cb_palette[1:clu_num], "black")
  names(col_li) <- c(gene_group_info$Label, "0")
}


groupInfo <- split(gene_group$Gid, gene_group$Label)
gene_group_ggtree <- ggtree(gene_clu, aes(color=group), layout='circular')
gene_group_ggtree <- groupOTU(gene_group_ggtree, groupInfo)
p <- gene_group_ggtree + 
  geom_tiplab(size=1.3) + 
  xlim(0, max(gene_dist_mat)) +
  scale_color_manual(values = col_li) +
  theme(
    text = element_text(family="ArialMT", size = 6)
)
ggsave(file.path(f_out, "ModularityEnrichGene.group.pdf"), p, height = 10, width = 13, units = "cm", colormodel = "cmyk")

gene_group <- gene_group[,c("Gid", "Label")]
names(gene_group)[2] <- "Group"
write_tsv(gene_group, file.path(f_out, "ModularityEnrichGene.group.tsv"))

plot_UMAP <- function(group_expr_gene_num, umapSCE, dim1, dim2, label){
  df <- data.frame(x=umapSCE@int_colData$reducedDims$UMAP[,dim1], y=umapSCE@int_colData$reducedDims$UMAP[,dim2], expr_gene_num=group_expr_gene_num)
  df_highlight <- df[df$expr_gene_num>0,]
  df_bg <- df[df$expr_gene_num==0,]
  
  p <- ggplot(data=NULL, aes(x=x, y=y)) +
    theme_bw() +
    geom_point(data= df_bg, size=0.1, alpha=0.1, color="grey70") +
    geom_point(data= df_highlight, color="red", size=0.1, alpha=0.1) +
    labs(x=sprintf("UMAP%d", dim1), y=sprintf("UMAP%d", dim2), title=label) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.position = "none",
      panel.grid = element_blank()
    )
  return(p)
}

plot_dir <- file.path(f_out, "GeneGroupUMAP")
if(!file.exists(plot_dir)){
  dir.create(plot_dir, recursive = TRUE)
}
inche_cm=2.54
for(g in gene_group_info$Label){
  gid_li <- gene_group$Gid[gene_group$Group==g]
  group_expr_gene_num <- colSums(expr_mat[rownames(expr_mat) %in% gid_li,])
  p <- plot_UMAP(group_expr_gene_num, corSCE, 1, 2, sprintf("%s (#Gene=%d, #Cell=%d)", g, length(gid_li), sum(group_expr_gene_num>0)))
  png(file.path(plot_dir, sprintf("%s.UMAP.png", g)), width=6, height=6, units="cm", family="ArialMT", res=400)
  print(p)
  dev.off()
}


