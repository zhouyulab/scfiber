library(readr)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(reshape2)

f_out <- "data/stat"

cluster_rep_cnt <- read_delim("data/cluster.rep.cnt.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Nanopore_feature_cnt <- read_delim("data/feature_counts.all.merge.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

all_expr_df <- inner_join(cluster_rep_cnt, Nanopore_feature_cnt)
cluster_rep_cnt <- all_expr_df[,c("Gid", "WT.rep2.C1", "WT.rep3.C1", "WT.rep2.C2", "WT.rep3.C2", "WT.rep2.C3", "WT.rep3.C3", "WT.rep2.C4", "WT.rep3.C4", "WT.rep2.C5", "WT.rep3.C5", "FL.rep2.C1", "FL.rep3.C1", "FL.rep2.C2", "FL.rep3.C2", "FL.rep2.C4", "FL.rep3.C4", "FL.rep2.C5", "FL.rep3.C5")]
Nanopore_feature_cnt <- all_expr_df[,names(Nanopore_feature_cnt)]

cor_raw_map <- cor(Nanopore_feature_cnt[,2:ncol(Nanopore_feature_cnt)])

cluster_rep_mat <- as.matrix(cluster_rep_cnt[,2:ncol(cluster_rep_cnt)])
rownames(cluster_rep_mat) <- cluster_rep_cnt$Gid
Nanopore_mat <- as.matrix(Nanopore_feature_cnt[,2:ncol(Nanopore_feature_cnt)])
rownames(Nanopore_mat) <- Nanopore_feature_cnt$Gid
Nanopore_read_per_sample <- colSums(Nanopore_mat)
Nanopore_read_per_sample_df <- data.frame(Sample=names(Nanopore_read_per_sample), TotalReadNum=Nanopore_read_per_sample)
Nanopore_read_per_sample_df$Sample <- factor(Nanopore_read_per_sample_df$Sample, levels = Nanopore_read_per_sample_df$Sample)
Nanopore_read_per_sample_df$Type <- sapply(strsplit(as.character(Nanopore_read_per_sample_df$Sample) , "[.]"), function(x){return(x[1])})
H2O_read_num <- Nanopore_read_per_sample_df$TotalReadNum[Nanopore_read_per_sample_df$Type=="ddH2O"]
read_num_cutoff <- 100000
Nanopore_read_per_sample_df$IsSignal <- Nanopore_read_per_sample_df$TotalReadNum>=read_num_cutoff
Nanopore_read_per_sample_df$IsSignal <- factor(Nanopore_read_per_sample_df$IsSignal, levels = c(TRUE, FALSE), labels = c("Signal", "Background"))
write_tsv(Nanopore_read_per_sample_df, file.path(f_out, "Nanopore.ReadPerSample.tsv"))

p <- ggplot(Nanopore_read_per_sample_df, aes(x=Sample, y=TotalReadNum, fill=IsSignal)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = read_num_cutoff, color="black", lty=2, size=0.5) +
  facet_grid(~Type, scales = "free_x") +
  scale_fill_manual(values = c("Signal"="red", "Background"="grey70")) +
  labs(y="#Mapped Nanopore reads") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(2, "mm"),
        panel.grid = element_blank())
ggsave(file.path(f_out, "Nanopore.ReadPerSample.pdf"), p, width = 20, height = 6, units = "cm")

signal_sample_df <- Nanopore_read_per_sample_df[Nanopore_read_per_sample_df$IsSignal=="Signal",]
signal_sample_df$Sample <- as.character(signal_sample_df$Sample)
Nanopore_feature_cnt <- Nanopore_feature_cnt[,c("Gid", signal_sample_df$Sample)]
Nanopore_mat <- as.matrix(Nanopore_feature_cnt[,2:ncol(Nanopore_feature_cnt)])
rownames(Nanopore_mat) <- Nanopore_feature_cnt$Gid


cluster_rep_cpm_mat <- 1e6 * t(t(cluster_rep_mat) / colSums(cluster_rep_mat))
Nanopore_cpm_mat <- 1e6 * t(t(Nanopore_mat) / colSums(Nanopore_mat))

min_scRNA_cpm <- min(cluster_rep_cpm_mat[cluster_rep_cpm_mat!=0])
min_Nanopore_cpm <- min(Nanopore_cpm_mat[Nanopore_cpm_mat!=0])

scRNA_expr_index <- (apply(cluster_rep_cpm_mat[,1:10], 1, mean)>1) | (apply(cluster_rep_cpm_mat[,11:18], 1, mean)>1) | (apply(cluster_rep_cpm_mat[,c("WT.rep2.C3", "WT.rep3.C3")], 1, mean)>1)
Nanopore_expr_index <- (apply(Nanopore_cpm_mat[,signal_sample_df$Sample[signal_sample_df$Type=="WT"]], 1, mean)>1) | (apply(Nanopore_cpm_mat[,signal_sample_df$Sample[signal_sample_df$Type=="FL"]], 1, mean)>1)
expr_index <- scRNA_expr_index | Nanopore_expr_index

expr_cluster_rep_cpm_mat <- cluster_rep_cpm_mat[expr_index,]
expr_Nanopore_cpm_mat <- Nanopore_cpm_mat[expr_index,]

spearman_cor_mat <- cor(expr_cluster_rep_cpm_mat, expr_Nanopore_cpm_mat, method = "spearman")
pearson_cor_mat <- cor(log10(expr_cluster_rep_cpm_mat+min_scRNA_cpm), log10(expr_Nanopore_cpm_mat+min_Nanopore_cpm), method = "pearson")

ht <- Heatmap(
  spearman_cor_mat,
  name = "Spearman Cor.",
  col=colorRamp2(c(0.3, 0.5, 0.7), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = c(rep("WT", 10), rep("FL", 8)),
  column_split = signal_sample_df$Type,
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5)
  ))

inche_cm=2.54
pdf(file.path(f_out, "Nanopore.scRNA.SpearmanCor.ExprGene.pdf"), width=5.5/inche_cm, height=8/inche_cm)
print(ht)
dev.off()


ht <- Heatmap(
  pearson_cor_mat,
  name = "Pearson Cor.",
  col=colorRamp2(c(0.3, 0.5, 0.7), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = c(rep("WT", 10), rep("FL", 8)),
  column_split = signal_sample_df$Type,
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5)
  ))
pdf(file.path(f_out, "Nanopore.scRNA.PearsonCor.ExprGene.pdf"), width=5.5/inche_cm, height=8/inche_cm)
print(ht)
dev.off()

Nanopore_spearman_cor_mat <- cor(expr_Nanopore_cpm_mat, method = "spearman")
Nanopore_pearson_cor_mat <- cor(log10(expr_Nanopore_cpm_mat+min_Nanopore_cpm), method = "pearson")
ht <- Heatmap(
  Nanopore_spearman_cor_mat,
  # col=colorRamp2(c(0.3, 0.5, 0.7), c("blue", "white", "red")),
  name = "Spearman Cor.",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = signal_sample_df$Type,
  column_split = signal_sample_df$Type,
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5)
  ))
inche_cm=2.54
pdf(file.path(f_out, "Nanopore.SpearmanCor.ExprGene.pdf"), width=5.5/inche_cm, height=4/inche_cm)
print(ht)
dev.off()



ht <- Heatmap(
  Nanopore_pearson_cor_mat,
  # col=colorRamp2(c(0.25, 0.5, 0.75), c("blue", "white", "red")),
  name = "Pearson Cor.",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = signal_sample_df$Type,
  column_split = signal_sample_df$Type,
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5)
  ))
inche_cm=2.54
pdf(file.path(f_out, "Nanopore.SpearmanCor.ExprGene.pdf"), width=5.5/inche_cm, height=4/inche_cm)
print(ht)
dev.off()

marker_gene <- read_delim("data/marker_gene.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Nanopore_feature_cnt <- Nanopore_feature_cnt[,c("Gid", "WT.Fiber.rep1", "WT.Fiber.rep3")]
Nanopore_feature_cnt <- left_join(Nanopore_feature_cnt, marker_gene)
all_marker_gene <- na.omit(Nanopore_feature_cnt)

melt_marker_gene <- melt(all_marker_gene, id.vars = c("Gid", "MarkerGene"), variable.name = "Sample", value.name = "UMI")
marker_gene_stat <- melt_marker_gene %>% group_by(MarkerGene, Sample) %>% summarise(AllMarkerGeneNum=n(), ExprMarkerGeneNum=sum(UMI>0), pval=wilcox.test(UMI, melt_marker_gene$UMI[melt_marker_gene$MarkerGene=="C3" & melt_marker_gene$Sample==Sample[1]])$p.value)
marker_gene_stat$Label <- sprintf("%d / %d\np=%.2e", marker_gene_stat$ExprMarkerGeneNum, marker_gene_stat$AllMarkerGeneNum, marker_gene_stat$pval)
p <- ggplot(melt_marker_gene, aes(x=MarkerGene, y=log10(UMI+1), fill=MarkerGene)) +
  geom_boxplot(outlier.colour = NA, size=0.3) +
  geom_text(data=marker_gene_stat, aes(x=MarkerGene, y=log10(max(melt_marker_gene$UMI+1)), label=Label), color="black", size=1.2) +
  facet_grid(~Sample) +
  labs(y="#Nanopore reads", x="Marker gene") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.grid = element_blank())
ggsave(file.path(f_out, "Nanopore.MarkerGene.pdf"), p, width = 10, height = 6, units = "cm")


