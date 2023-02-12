library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(readr)
library(BiocParallel)

f_out <- "analysis_v3/Nanopore/stat"
UMAP_df <- as.data.frame(corSCE@int_colData$reducedDims$UMAP)
names(UMAP_df) <- c("UMAP1", "UMAP2")
UMAP_df$Cluster <- corSCE$cluster
UMAP_df$Rep <- corSCE$Rep
UMAP_df$Barcode <- corSCE$Barcode
bin_size <- 0.5
UMAP_df$UMAP1_bin <- as.integer(UMAP_df$UMAP1 / bin_size)
UMAP_df$UMAP2_bin <- as.integer(UMAP_df$UMAP2 / bin_size)
UMAP_df$index <- 1:nrow(UMAP_df)
UMAP_df$Sample <- corSCE$SampleName
UMAP_df$Label <- sprintf("%s_%s_%s", UMAP_df$Sample, UMAP_df$UMAP1_bin, UMAP_df$UMAP2_bin)

scRNA_cnt <- counts(mergeSCE)
fc_merge <- read_delim("data/Nanopore/feature_counts.all.merge.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
fc_merge <- fc_merge[fc_merge$Gid%in%rownames(scRNA_cnt),]
names(fc_merge) <- c("Gid", sprintf("WT.Fiber.%s", 1:24), "ddH2O.1", sprintf("FL.Ovule.%s", 1:24), "ddH2O.2")
fc_merge <- fc_merge[,c("Gid", sprintf("WT.Fiber.%s", 1:24), sprintf("FL.Ovule.%s", 1:24), "ddH2O.1", "ddH2O.2")]

Nanopore_read_per_sample <- apply(fc_merge[,2:ncol(fc_merge)], 2, sum)
read_num_cutoff <- 200000
fc_merge <- fc_merge[,c("Gid", names(Nanopore_read_per_sample)[Nanopore_read_per_sample>read_num_cutoff], "ddH2O.1", "ddH2O.2")]
fc_mat <- as.matrix(fc_merge[, 2:ncol(fc_merge)])
rownames(fc_mat) <- fc_merge$Gid
fc_cpm_mat <- 1e6 * t(t(fc_mat) / colSums(fc_mat))
sample_flag <- sapply(strsplit(colnames(fc_cpm_mat), "[.]"), function(x){return(x[1])}) 
expr_gids <- rownames(fc_cpm_mat)[(apply(fc_cpm_mat[,sample_flag=="WT"], 1, mean) > 1) | (apply(fc_cpm_mat[,sample_flag=="FL"], 1, mean) > 1)]
expr_fc_df <- left_join(data.frame(Gid=expr_gids), fc_merge)
expr_fc_mat <- as.matrix(expr_fc_df[, 2:ncol(expr_fc_df)])
rownames(expr_fc_mat) <- expr_fc_df$Gid
write_tsv(expr_fc_df, file.path(f_out, "Nanopore.Expr.tsv"))


all_cor_df <- data.frame()
for(umap_label in unique(UMAP_df$Label)){
  tmp_umap_df <- UMAP_df[UMAP_df$Label==umap_label,]
  if(nrow(tmp_umap_df)==1){
    next()
  }
  tmp_cnt_mat <- scRNA_cnt[,tmp_umap_df$index]
  tmp_gene_umi <- rowSums(tmp_cnt_mat)
  tmp_gene_umi_df <- data.frame(Gid=names(tmp_gene_umi), scRNA_UMI=tmp_gene_umi)
  expr_tmp_gene_umi_df <- left_join(data.frame(Gid=expr_gids), tmp_gene_umi_df)
  tmp_spearman_cor <- cor(expr_tmp_gene_umi_df$scRNA_UMI, expr_fc_mat, method = "spearman")
  tmp_pearson_cor <- cor(log10(expr_tmp_gene_umi_df$scRNA_UMI+1), log10(expr_fc_mat+1), method = "pearson")
  tmp_cor_df <- data.frame(
    Sample=tmp_umap_df$Sample[1], 
    UMAP1_bin=tmp_umap_df$UMAP1_bin[1],
    UMAP2_bin=tmp_umap_df$UMAP2_bin[1],
    NanoporeSample=colnames(tmp_spearman_cor),
    SpearmanCor=tmp_spearman_cor[1,],
    PearsonCor=tmp_pearson_cor[1,]
    )
  all_cor_df <- rbind(all_cor_df, tmp_cor_df)
}
write_tsv(all_cor_df, file.path(f_out, "Nanopore.scRNA.cor.tsv"))
all_cor_df$UMAP1_bin
p <- ggplot(all_cor_df, aes(xmin=UMAP1_bin*0.5, xmax=UMAP1_bin*0.5+0.5, ymin=UMAP2_bin*0.5, ymax=UMAP2_bin*0.5+0.5, fill=SpearmanCor)) +
  geom_rect() +
  facet_grid(NanoporeSample~Sample) +
  scale_fill_gradient2(low="blue", high="red", midpoint = 0.45, mid="grey90") +
  labs(x="UMAP1", y="UMAP2") +
  theme_bw() +
  theme(
    axis.text = element_text(family="ArialMT", color = "black", size = 6),
    panel.grid = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.title = element_text(family="ArialMT", color = "black", size = 6),
    legend.text = element_text(family="ArialMT", color = "black", size = 6)
  )
ggsave(file.path(f_out, "Nanopore.scRNA.cor.pdf"), p, width = 15, height = 2+ 5*length(unique(all_cor_df$NanoporeSample)), limitsize = FALSE, units = "cm")


Fiber_NGS_df <- read_delim("data/Nanopore/Fiber.DPA5.feature_counts.merge.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Fiber_fc_mat <- as.matrix(Fiber_NGS_df[, 2:ncol(Fiber_NGS_df)])
rownames(Fiber_fc_mat) <- Fiber_NGS_df$Gid
Fiber_fc_cpm_mat <- 1e6 * t(t(Fiber_fc_mat) / colSums(Fiber_fc_mat))
Fiber_expr_gids <- rownames(Fiber_fc_mat)[apply(Fiber_fc_cpm_mat, 1, mean) > 1]
Fiber_expr_fc_df <- left_join(data.frame(Gid=Fiber_expr_gids), Fiber_NGS_df)
Fiber_expr_fc_mat <- as.matrix(Fiber_expr_fc_df[, 2:ncol(Fiber_expr_fc_df)])
rownames(Fiber_expr_fc_mat) <- Fiber_expr_fc_df$Gid
write_tsv(Fiber_expr_fc_df, file.path(f_out, "Fiber.NGS.Expr.tsv"))

Fiber_cor_df <- data.frame()
for(umap_label in unique(UMAP_df$Label)){
  tmp_umap_df <- UMAP_df[UMAP_df$Label==umap_label,]
  if(nrow(tmp_umap_df)==1){
    next()
  }
  tmp_cnt_mat <- scRNA_cnt[,tmp_umap_df$index]
  tmp_gene_umi <- rowSums(tmp_cnt_mat)
  tmp_gene_umi_df <- data.frame(Gid=names(tmp_gene_umi), scRNA_UMI=tmp_gene_umi)
  expr_tmp_gene_umi_df <- left_join(data.frame(Gid=Fiber_expr_gids), tmp_gene_umi_df)
  tmp_spearman_cor <- cor(expr_tmp_gene_umi_df$scRNA_UMI, Fiber_expr_fc_mat, method = "spearman")
  tmp_pearson_cor <- cor(log10(expr_tmp_gene_umi_df$scRNA_UMI+1), log10(Fiber_expr_fc_mat+1), method = "pearson")
  tmp_cor_df <- data.frame(
    Sample=tmp_umap_df$Sample[1], 
    UMAP1_bin=tmp_umap_df$UMAP1_bin[1],
    UMAP2_bin=tmp_umap_df$UMAP2_bin[1],
    FiberNgsRep=colnames(tmp_spearman_cor),
    SpearmanCor=tmp_spearman_cor[1,],
    PearsonCor=tmp_pearson_cor[1,]
  )
  Fiber_cor_df <- rbind(Fiber_cor_df, tmp_cor_df)
}
write_tsv(Fiber_cor_df, file.path(f_out, "FiberNGS.scRNA.cor.tsv"))

p <- ggplot(Fiber_cor_df, aes(xmin=UMAP1_bin*0.5, xmax=UMAP1_bin*0.5+0.5, ymin=UMAP2_bin*0.5, ymax=UMAP2_bin*0.5+0.5, fill=SpearmanCor)) +
  geom_rect() +
  facet_grid(FiberNgsRep~Sample) +
  scale_fill_gradient2(low="blue", high="red", midpoint = 0.45, mid="grey90") +
  labs(x="UMAP1", y="UMAP2") +
  theme_bw() +
  theme(
    axis.text = element_text(family="ArialMT", color = "black", size = 6),
    panel.grid = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.title = element_text(family="ArialMT", color = "black", size = 6),
    legend.text = element_text(family="ArialMT", color = "black", size = 6)
  )
ggsave(file.path(f_out, "FiberNGS.scRNA.cor.pdf"), p, width = 15, height = 2+ 5*length(unique(Fiber_cor_df$FiberNgsRep)), limitsize = FALSE, units = "cm")



