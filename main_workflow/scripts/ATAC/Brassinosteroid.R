library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(entropy)

Hormones_transduction <- read_delim("data/Hormones/Transduction/Hormones.transduction.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Brassinosteroid <- Hormones_transduction[Hormones_transduction$Hormones=="Brassinosteroid",]
Gbox_WT_enrich_gene <- read_csv("analysis_v3/ATAC/ATAC_DE_peak2C3/ATAC_DE_peak_target_gene.Gbox.WT_enrich.gid.txt")
Brassinosteroid_scRNA_mat <- as.matrix(Brassinosteroid[,5:ncol(Brassinosteroid)])
rownames(Brassinosteroid_scRNA_mat) <- sprintf("%s (%s)", Brassinosteroid$Gid, Brassinosteroid$Symbol)
Brassinosteroid_info <- Brassinosteroid[,1:4]
Brassinosteroid_info$IsGboxTargetGene <- Brassinosteroid_info$Gid %in% Gbox_WT_enrich_gene
CPM_merge <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Brassinosteroid_TP <- left_join(Brassinosteroid_info, CPM_merge)
Brassinosteroid_TP_mat <- as.matrix(Brassinosteroid_TP[,6:ncol(Brassinosteroid_TP)])
rownames(Brassinosteroid_TP_mat) <- sprintf("%s (%s)", Brassinosteroid_TP$Gid, Brassinosteroid_TP$Symbol)

norm_Brassinosteroid_scRNA_mat <- Brassinosteroid_scRNA_mat / apply(Brassinosteroid_scRNA_mat, 1, function(x){return(max(x, na.rm = T))})
norm_Brassinosteroid_TP_mat <- Brassinosteroid_TP_mat / apply(Brassinosteroid_TP_mat, 1, max)


Brassinosteroid_TP$TP_cv <- apply(Brassinosteroid_TP_mat, 1, function(x){return(sd(x)/mean(x))})
Brassinosteroid_TP$scRNA_cv <- apply(Brassinosteroid_scRNA_mat, 1, function(x){return(sd(x, na.rm = T)/mean(x, na.rm = T))})


scRNA_dist <- dist(norm_Brassinosteroid_scRNA_mat)
TP_dist <- dist(norm_Brassinosteroid_TP_mat)
TP_dist[is.na(TP_dist)] <- max(TP_dist, na.rm = T)
merge_dist <- scRNA_dist + TP_dist
clu <- hclust(merge_dist, method = "ward.D2")

ht_scRNA <- Heatmap(
  norm_Brassinosteroid_scRNA_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  split = Brassinosteroid_info$Symbol,
  column_split = c(rep("WT", 5), rep("FL", 5)),
  cluster_columns = FALSE, 
  row_order = clu$order,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

ht_TP <- Heatmap(
  norm_Brassinosteroid_TP_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  split = Brassinosteroid_info$Symbol,
  column_split = c(rep("WT", 9), rep("FL", 9)),
  cluster_columns = FALSE, 
  row_order = clu$order,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

inche_cm <- 2.54
pdf("analysis_v3/ATAC/ATAC_DE_peak2C3/Brassinosteroid.expr.ht.pdf", width=15.5/inche_cm, height=8/inche_cm)
print(ht_scRNA+ht_TP)
dev.off()

Brassinosteroid_cv_df <- melt(Brassinosteroid_TP[,c("KO", "Symbol", "Gid", "TP_cv",  "scRNA_cv")], c("KO", "Symbol", "Gid"), variable.name = "Stat", value.name = "CV")
Brassinosteroid_cv_df <- na.omit(Brassinosteroid_cv_df)
Brassinosteroid_cv_df$Stat <- factor(Brassinosteroid_cv_df$Stat, levels = c("scRNA_cv", "TP_cv"), labels = c("scRNA-seq", "RNA-seq"))
Brassinosteroid_cv_df$Symbol <- factor(Brassinosteroid_cv_df$Symbol, levels = Brassinosteroid$Symbol[!duplicated(Brassinosteroid$Symbol)])
p <- ggplot(Brassinosteroid_cv_df, mapping = aes(x=Symbol, y=CV, fill=Stat)) +
  geom_boxplot(color="black", size=0.3, outlier.colour = NA) +
  geom_point(color="black", position = position_jitter(width = 0.2), size=0.4) +
  facet_grid(~Stat) +
  labs(y="Coefficient of Variation") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
ggsave("analysis_v3/ATAC/ATAC_DE_peak2C3/Brassinosteroid.cv.pdf", p, width = 7, height = 5, limitsize = FALSE, units = "cm")
