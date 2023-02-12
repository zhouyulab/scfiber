library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(VennDiagram)
library(RColorBrewer)

refGene <- read_delim("data/genome/HAU/Ghirsutumv1.1_gene_model.gene.bed6", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(refGene) <- c("Chrom", "Start", "End", "Gid", "Score", "Strand")
ATAC_peak_DE <- read_delim("analysis_v3/stat/ATAC_MACS2/ATAC.peak.DE.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
cluster_tpm <- read_delim("analysis_v3/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
marker_gene <- read_delim("analysis_v3/stat/marker_gene/Marker.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
CPM_TP <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
expr_time_df <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/TimeCourse/ImpulseDE2_enrich_gene.CPM.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
DiffExpr_C2_Down <- read_delim("~/mu01/project/CottonSingleCell/analysis_v3/stat/diff_expr/DiffExpr.C2.Down.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
public_RNA_seq <- read_delim("analysis_v3/RNA_seq/merge_expr/TP.FPKM.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
public_RNA_seq$`5DPA` <- NULL


motif_hit_df <- read_delim("analysis_v3/ATAC/motif/summary/MotifExprTargetGene.HitMat.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

Gbox <- data.frame(Gid=motif_hit_df$Gid[motif_hit_df$Gbox])

f_out <- "analysis_v3/ATAC/ATAC_DE_peak2C3"



ATAC_peak_DE$Peak <- sprintf("%s_%s_%s", ATAC_peak_DE$Chrom, ATAC_peak_DE$Start, ATAC_peak_DE$End)
WT_peak_DE <- ATAC_peak_DE[ATAC_peak_DE$Tag=="WT",]
refGene$WT_peak_DE <- apply(refGene, 1, function(line, expand=1000){
  gene_chrom <- line[1]
  gene_start <- as.integer(line[2])
  gene_end <- as.integer(line[3])
  gene_strand <- line[6]
  if(gene_strand == "+"){
    start_pos <- gene_start - expand
    end_pos <- gene_end
  }else{
    start_pos <- gene_start
    end_pos <- gene_end + expand
  }
  tmp_WT_peak_DE <- WT_peak_DE[WT_peak_DE$Chrom==gene_chrom,]
  res <- any(!((tmp_WT_peak_DE$Start >= end_pos) | (start_pos >= tmp_WT_peak_DE$End)))
  return(res)
})
refGene$GboxPromoter <- refGene$Gid %in% Gbox$Gid
Gbox_cnt_mat <- table(refGene[,c("WT_peak_DE", "GboxPromoter")])
chisq.test(Gbox_cnt_mat)

WT_DE_ATAC_peak <- refGene[refGene$WT_peak_DE,]

expr_tmp <- cluster_tpm[apply(cluster_tpm[,2:ncol(cluster_tpm)], 1, function(x){return(max(x, na.rm = TRUE))})>1,]
expr_tmp$IsC3Max <- apply(expr_tmp[2:ncol(expr_tmp)], 1, which.max)==3

C3_marker_gene <- marker_gene[marker_gene$Cluster=="C3",]
C1_marker_gene <- marker_gene[marker_gene$Cluster=="C1",]
C2_marker_gene <- marker_gene[marker_gene$Cluster=="C2",]
C2_Down_gene <- DiffExpr_C2_Down$Gid

pdf(file.path(f_out, "ATAC_DE_peak_target_gene.C3_marker_C3_max.venn.pdf"), width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(
  ATAC_DE_WT_target_gene=unique(WT_DE_ATAC_peak$Gid), 
  C3_Max_Expr_Gene=expr_tmp$Gid[expr_tmp$IsC3Max],
  C3_marker_gene=C3_marker_gene$Gid
), filename = NULL,
fill = brewer.pal(3, "Dark2"),
borders = FALSE,
sub.cex = 0.3,
lwd=0.5)
grid.draw(p)
dev.off()

pdf(file.path(f_out, "ATAC_DE_peak_target_gene.C3_marker.venn.pdf"), width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(
  ATAC_DE_WT_target_gene=unique(WT_DE_ATAC_peak$Gid), 
  C3_marker_gene=C3_marker_gene$Gid,
  C2_WT_DEGs=C2_Down_gene
), filename = NULL,
fill = brewer.pal(3, "Dark2"),
borders = FALSE,
sub.cex = 0.3,
lwd=0.5)
grid.draw(p)
dev.off()

pdf(file.path(f_out, "ATAC_DE_peak_target_gene.C3_max.venn.pdf"), width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(
  ATAC_DE_WT_target_gene=unique(WT_DE_ATAC_peak$Gid), 
  C3_Max_Expr_Gene=expr_tmp$Gid[expr_tmp$IsC3Max]
), filename = NULL,
fill = brewer.pal(3, "Dark2")[c(1, 2)],
borders = FALSE,
sub.cex = 0.3,
lwd=0.5)
grid.draw(p)
dev.off()

pdf(file.path(f_out, "ATAC_DE_peak_target_gene.marker.venn.pdf"), width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(
  ATAC_DE_WT_target_gene=unique(WT_DE_ATAC_peak$Gid), 
  C1_marker_gene=C1_marker_gene$Gid,
  C2_marker_gene=C2_marker_gene$Gid,
  C3_marker_gene=C3_marker_gene$Gid,
  C2_WT_DEGs=C2_Down_gene
), filename = NULL,
fill = brewer.pal(5, "Dark2"),
borders = FALSE,
sub.cex = 0.3,
lwd=0.5)
grid.draw(p)
dev.off()

upset_df <- data.frame(Gid=refGene$Gid)
upset_df$ATAC_DE_WT_target_gene <- as.integer(upset_df$Gid %in% WT_DE_ATAC_peak$Gid)
upset_df$C1_marker_gene <- as.integer(upset_df$Gid %in% C1_marker_gene$Gid)
upset_df$C2_marker_gene <- as.integer(upset_df$Gid %in% C2_marker_gene$Gid)
upset_df$C3_marker_gene <- as.integer(upset_df$Gid %in% C3_marker_gene$Gid)
upset_df$C2_WT_DEGs <- as.integer(upset_df$Gid %in% C2_Down_gene)
upset_df$WT_DE_Gbox_promoter_gene <- as.integer(upset_df$Gid %in% Gbox$Gid)
upset_df$WT_DE_Gbox_promoter_gene[upset_df$ATAC_DE_WT_target_gene==0] <- 0
# upset_df$ATAC_DE_WT_target_gene <- NULL
p <- upset(upset_df, nsets =10, order.by = c("freq"), sets.x.label="#Gene")
inche_cm <- 2.54
pdf("ATAC_DE_peak_target_gene.marker.upset.pdf", width=20/inche_cm, height=20/inche_cm)
print(p)
dev.off()

chisq_df <-data.frame(
  DataSet=c(rep("ATAC_DE_WT_target_gene", 4), rep("WT_DE_Gbox_promoter_gene", 4)),
  GeneSet=rep(c("C1_marker_gene", "C2_marker_gene", "C3_marker_gene", "C2_WT_DEGs"), 2))
chisq_df$Pval <- apply(chisq_df, 1, function(x){return(chisq.test(table(upset_df[,x]))$p.value)})

upset_df$Gid[upset_df$C2_marker_gene&upset_df$WT_DE_Gbox_promoter_gene]


WT_C2_Down_ATAC_DE_gene <- DiffExpr_C2_Down[DiffExpr_C2_Down$Gid %in% WT_DE_ATAC_peak$Gid,]
write_tsv(WT_C2_Down_ATAC_DE_gene, file.path(f_out, "WT_C2_Down_ATAC_DE_gene.tsv"))

melt_expr_CPM <- melt(expr_tmp, id.vars = c("Gid", "IsC3Max"), variable.name = "Label", value.name = "CPM")
melt_expr_CPM$Sample <- sapply(strsplit(as.character(melt_expr_CPM$Label), "[.]"), function(x){return(x[1])})
melt_expr_CPM$Cluster <- sapply(strsplit(as.character(melt_expr_CPM$Label), "[.]"), function(x){return(x[2])})
melt_expr_CPM$Sample <- factor(melt_expr_CPM$Sample, levels = c("WT", "FL"))
melt_expr_CPM$Cluster <- factor(melt_expr_CPM$Cluster, levels = c("C1", "C2", "C3", "C4", "C5"))

dcast_expr_CPM <- dcast(melt_expr_CPM[,c("Gid", "Sample", "Cluster", "CPM")], Gid+Cluster~Sample, value.var = "CPM")

ATAC_WT_DE_peak_expr_gene <- melt_expr_CPM[melt_expr_CPM$Gid%in%WT_DE_ATAC_peak$Gid,]
ATAC_WT_DE_peak_expr_gene$Sample <- sapply(strsplit(as.character(ATAC_WT_DE_peak_expr_gene$Label), "[.]"), function(x){return(x[1])})
ATAC_WT_DE_peak_expr_gene$Cluster <- sapply(strsplit(as.character(ATAC_WT_DE_peak_expr_gene$Label), "[.]"), function(x){return(x[2])})
ATAC_WT_DE_peak_expr_gene$Sample <- factor(ATAC_WT_DE_peak_expr_gene$Sample, levels = c("WT", "FL"))
ATAC_WT_DE_peak_expr_gene$Cluster <- factor(ATAC_WT_DE_peak_expr_gene$Cluster, levels = c("C1", "C2", "C3", "C4", "C5"))
write_tsv(ATAC_WT_DE_peak_expr_gene, file.path(f_out, "ATAC_DE_peak_target_gene.expr_gene.CPM.SourceData.tsv"))

ATAC_WT_DE_peak_expr_gene_info <- ATAC_WT_DE_peak_expr_gene %>% 
  group_by(Cluster, Sample) %>% summarise(
    median_CPM=median(CPM)
    )
ATAC_WT_DE_peak_expr_gene_info$Label <- sprintf("Median=%.2f", ATAC_WT_DE_peak_expr_gene_info$median_CPM)
ATAC_WT_DE_peak_expr_gene_pval <- ATAC_WT_DE_peak_expr_gene[ATAC_WT_DE_peak_expr_gene$Cluster!="C3",] %>% 
  group_by(Cluster) %>% summarise(
    Pval=wilcox.test(CPM[Sample=="WT"], CPM[Sample=="FL"])$p.value
  )
ATAC_WT_DE_peak_expr_gene_pval$Label <- sprintf("P-value=%.1E", ATAC_WT_DE_peak_expr_gene_pval$Pval)

p <- ggplot(ATAC_WT_DE_peak_expr_gene, aes(x=Cluster, y=log10(CPM+0.01), fill=Sample)) +
  geom_boxplot(size=0.2, outlier.colour = NA, position = position_dodge(0.7), width=0.5) +
  geom_text(data=ATAC_WT_DE_peak_expr_gene_info, mapping = aes(y=3.5, label=Label), position = position_dodge(0.7), fill=NA, color="black", size=1.2) +
  geom_text(data=ATAC_WT_DE_peak_expr_gene_pval, mapping = aes(y=3.2, label=Label, fill=NA), color="black", size=1.2) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(breaks = log10(c(0, 0.1, 1, 10, 100, 1000)+0.01), labels = c("0", "1E-1", "1E0", "1E1", "1E2", "1E3")) +
  labs(y="CPM") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size=5),
    title = element_text(family="ArialMT", size=5),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(family="ArialMT", size=5),
    legend.background = element_blank(),
    legend.text = element_text(family="ArialMT", size=5),
    panel.grid = element_blank()
  )

ggsave(file.path(f_out, "ATAC_DE_peak_target_gene.CPM.box_plot.pdf"), p, height = 4, width = 8, units = "cm", colormodel = "cmyk")

ATAC_WT_DE_peak_dcast_expr_CPM <- dcast_expr_CPM[dcast_expr_CPM$Gid%in%WT_DE_ATAC_peak$Gid,]
ATAC_WT_DE_peak_dcast_expr_CPM$Log2FC <- log2((ATAC_WT_DE_peak_dcast_expr_CPM$FL+0.01) / (ATAC_WT_DE_peak_dcast_expr_CPM$WT+0.01))
ATAC_WT_DE_peak_dcast_expr_CPM <- na.omit(ATAC_WT_DE_peak_dcast_expr_CPM)
ATAC_WT_DE_peak_dcast_expr_CPM$DE <- "NC"
ATAC_WT_DE_peak_dcast_expr_CPM$DE[ATAC_WT_DE_peak_dcast_expr_CPM$Log2FC>1] <- "FL"
ATAC_WT_DE_peak_dcast_expr_CPM$DE[ATAC_WT_DE_peak_dcast_expr_CPM$Log2FC< -1] <- "WT"
ATAC_WT_DE_peak_dcast_expr_CPM$DE[ATAC_WT_DE_peak_dcast_expr_CPM$WT<1 & ATAC_WT_DE_peak_dcast_expr_CPM$FL<1] <- "NC"
ATAC_WT_DE_peak_dcast_expr_CPM_info <- ATAC_WT_DE_peak_dcast_expr_CPM %>% group_by(Cluster) %>% summarise(Label=sprintf("#WT=%s\n#FL=%s\n#NC=%s", sum(DE=="WT"), sum(DE=="FL"), sum(DE=="NC")))
p <- ggplot(ATAC_WT_DE_peak_dcast_expr_CPM, aes(x=log10(WT+0.01), y=log10(FL+0.01), color=DE)) +
  geom_point(size=0.2, alpha=0.05) +
  geom_text(data=ATAC_WT_DE_peak_dcast_expr_CPM_info, mapping = aes(x=-0.5, y=3, label=Label), color="black", size=1.2) +
  scale_color_manual(values = list("WT"="red", "FL"="blue", "NC"="grey70")) +
  facet_grid(~Cluster) +
  scale_x_continuous(breaks = log10(c(0, 0.1, 1, 10, 100, 1000)+0.01), labels = c("0", "1E-1", "1E0", "1E1", "1E2", "1E3")) +
  scale_y_continuous(breaks = log10(c(0, 0.1, 1, 10, 100, 1000)+0.01), labels = c("0", "1E-1", "1E0", "1E1", "1E2", "1E3")) +
  labs(x="WT CPM", y="FL CPM") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size=5),
    title = element_text(family="ArialMT", size=5),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(family="ArialMT", size=5),
    legend.background = element_blank(),
    legend.text = element_text(family="ArialMT", size=5),
    panel.grid = element_blank()
  )
ggsave(file.path(f_out, "ATAC_DE_peak_target_gene.CPM.WT_FL.pdf"), p, height = 4, width = 14, units = "cm", colormodel = "cmyk")

write_tsv(ATAC_WT_DE_peak_dcast_expr_CPM, file.path(f_out, "ATAC_DE_peak_target_gene.expr_gene.CPM_by_clu.SourceData.tsv"))

CPM_TP_mat <- as.matrix(CPM_TP[,2:ncol(CPM_TP)])
rownames(CPM_TP_mat) <- CPM_TP$Gid
norm_CPM_TP_mat <- CPM_TP_mat / apply(CPM_TP_mat, 1, max)

ATAC_WT_DE_peak_TP_CPM <- CPM_TP[CPM_TP$Gid%in%WT_DE_ATAC_peak$Gid,]
ATAC_WT_DE_peak_TP_CPM_mat <- ATAC_WT_DE_peak_TP_CPM[,2:ncol(ATAC_WT_DE_peak_TP_CPM)]
rownames(ATAC_WT_DE_peak_TP_CPM_mat) <- ATAC_WT_DE_peak_TP_CPM$Gid

expr_TP_stat <- left_join(CPM_TP[,c("Gid")], expr_time_df[,c("Gid", "Group")])
expr_TP_stat$Group[is.na(expr_TP_stat$Group)] <- "Other"
expr_TP_stat$IsWtAtacPeak <- expr_TP_stat$Gid %in% WT_DE_ATAC_peak$Gid

expr_TP_stat_info <- expr_TP_stat %>% group_by(Group, IsWtAtacPeak)

DE_ATAC_peak_expr_time_mat <- as.matrix(table(expr_TP_stat$IsWtAtacPeak, expr_TP_stat$Group))
chi_test <- chisq.test(DE_ATAC_peak_expr_time_mat)

topGO <- read_delim("analysis_v3/ATAC/ATAC_DE_peak2C3/topGO/topGO.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TF_activate_gids <- unique(unlist(strsplit(topGO$Gids[topGO$GO.ID=="GO:0003700"], "[ ]")))

DE_TF_scRNA_df <- left_join(data.frame(Gid=TF_activate_gids), cluster_tpm)
DE_TF_TP_df <- left_join(data.frame(Gid=TF_activate_gids), CPM_TP)
DE_TF_public_TP_df <- left_join(data.frame(Gid=TF_activate_gids), public_RNA_seq)

DE_TF_scRNA_mat <- as.matrix(DE_TF_scRNA_df[,2:ncol(DE_TF_scRNA_df)])
rownames(DE_TF_scRNA_mat) <- DE_TF_scRNA_df$Gid
norm_DE_TF_scRNA_mat <- DE_TF_scRNA_mat / apply(DE_TF_scRNA_mat, 1, function(x){return(max(x, na.rm = T))})

DE_TF_TP_mat <- as.matrix(DE_TF_TP_df[,2:ncol(DE_TF_TP_df)])
rownames(DE_TF_TP_mat) <- DE_TF_TP_df$Gid
norm_DE_TF_TP_mat <- DE_TF_TP_mat / apply(DE_TF_TP_mat, 1, max)

DE_TF_public_TP_mat <- as.matrix(DE_TF_public_TP_df[,2:ncol(DE_TF_public_TP_df)])
rownames(DE_TF_public_TP_mat) <- DE_TF_public_TP_df$Gid
DE_TF_public_TP_mat[rowSums(DE_TF_public_TP_mat)==0,] <- NA
norm_DE_TF_public_TP_mat <- DE_TF_public_TP_mat / apply(DE_TF_public_TP_mat, 1, function(x){return(max(x, na.rm = T))})

tf_gene <- read_delim("analysis_v3/annotation/TF/TF.gene_list.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
interest_gene <- read_delim("data/InterestGeneID.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
interest_gene <- interest_gene[!duplicated(interest_gene[,c("Gid")]),]
interest_gene <- left_join(interest_gene, tf_gene)

interest_gene_df <- left_join(data.frame(Gid=TF_activate_gids), interest_gene)
interest_gene_df$Label <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$GeneSymbol)
interest_gene_df$Label[!is.na(interest_gene_df$Family)] <- sprintf("%s (%s: %s)", interest_gene_df$Gid, interest_gene_df$Family, interest_gene_df$GeneSymbol)[!is.na(interest_gene_df$Family)]
interest_gene_df$Label[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$Family)[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]
interest_gene_df$Label[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- interest_gene_df$Gid[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]

subindx <- which(!is.na(interest_gene_df$GeneSymbol))
label_li <- interest_gene_df$Label[subindx]

scRNA_dist <- dist(norm_DE_TF_scRNA_mat)
TP_dist <- dist(norm_DE_TF_TP_mat)
public_TP_dist <- dist(DE_TF_public_TP_mat)
scRNA_dist[is.na(scRNA_dist)] <- max(scRNA_dist, na.rm = T)
TP_dist[is.na(TP_dist)] <- max(TP_dist, na.rm = T)
public_TP_dist[is.na(public_TP_dist)] <- max(public_TP_dist, na.rm = T)
merge_dist <- 1*scRNA_dist + 1*TP_dist + 1 * public_TP_dist
clu <- hclust(merge_dist, method = "ward.D2")

scRNA_ht <- Heatmap(
  norm_DE_TF_scRNA_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  row_order = clu$order,
  show_row_names = F, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

TP_ht <- Heatmap(
  norm_DE_TF_TP_mat, 
  name="Norm.Expr.",
  # col = colorRamp2(c(0, 1), c("white", "red")),
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  row_order = clu$order,
  cluster_columns = FALSE, 
  show_row_names = F, 
  column_split = c(rep("WT", 9), rep("FL", 9)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

public_TP_ht <- Heatmap(
  norm_DE_TF_public_TP_mat, 
  name="Norm.Expr.",
  # col = colorRamp2(c(0, 1), c("white", "red")),
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  row_order = clu$order,
  cluster_columns = FALSE, 
  show_row_names = F, 
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

ha <- rowAnnotation(
  link = row_anno_link(
    at = subindx,
    labels = label_li,
    labels_gp=gpar(fontsize = 6),
    lines_gp=gpar(lwd = 0.5),
    link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)))

inche_cm <- 2.54
pdf(file.path(f_out, "ATAC_DE_peak_target_gene.TF_GO.ht.pdf"), width=15.5/inche_cm, height=8/inche_cm)
print(scRNA_ht+TP_ht+public_TP_ht+ha)
dev.off()


expr_clu_tpm <- cluster_tpm[apply(cluster_tpm[2:ncol(cluster_tpm)], 1, function(x){return(max(x, na.rm = T))})>1,]
max_expr_cell <- data.frame(Gid=expr_clu_tpm$Gid, MaxExprCell=apply(expr_clu_tpm[,2:ncol(cluster_tpm)], 1, function(x){
  max_expr_cell <- strsplit(names(cluster_tpm)[2:ncol(cluster_tpm)][which.max(x)], "[.]")[[1]][2]
  return(max_expr_cell)
}
  ))

Gbox_max_expr_cell <- inner_join(Gbox, max_expr_cell)
Gbox_max_expr_cell_info <- Gbox_max_expr_cell %>% group_by(MaxExprCell) %>% summarise(Data="Gbox", Indx=1, GeneNum=n())
Gbox_gene_num <- nrow(Gbox_max_expr_cell)
all_gene_num <- nrow(max_expr_cell)
Gbox_rand_max_expr_cell_info <- data.frame()
for (i in 1:10000) {
  tmp_max_expr_cell <- max_expr_cell[sample(1:all_gene_num, Gbox_gene_num, replace = FALSE),]
  tmp_max_expr_cell <- tmp_max_expr_cell %>% group_by(MaxExprCell) %>% summarise(Data="Random", Indx=i, GeneNum=n())
  Gbox_rand_max_expr_cell_info <- rbind(Gbox_rand_max_expr_cell_info, tmp_max_expr_cell)
}
Gbox_rand_max_expr_cell_info <- left_join(Gbox_rand_max_expr_cell_info, data.frame(MaxExprCell=Gbox_max_expr_cell_info$MaxExprCell, GboxGeneNum=Gbox_max_expr_cell_info$GeneNum))
Gbox_rand_max_expr_cell_df <- Gbox_rand_max_expr_cell_info %>% group_by(MaxExprCell, Data) %>% summarise(Mu=mean(GeneNum), Sd=sd(GeneNum), Pval=min(mean(GeneNum<GboxGeneNum), mean(GeneNum>GboxGeneNum))*2)
Gbox_max_expr_cell_df <- data.frame(MaxExprCell=Gbox_max_expr_cell_info$MaxExprCell, Data="Gbox", Mu=Gbox_max_expr_cell_info$GeneNum, Sd=0, Pval=NA)
Gbox_merge_df <- rbind(as.data.frame(Gbox_max_expr_cell_df), as.data.frame(Gbox_rand_max_expr_cell_df))
Gbox_merge_df <- Gbox_merge_df[!is.na(Gbox_merge_df$MaxExprCell),]
write_tsv(Gbox_merge_df, file.path(f_out, "Gbox.MaxExprCell.Expr.test.tsv"))
Gbox_merge_df$Data <- factor(Gbox_merge_df$Data, levels = c("Gbox", "Random"))
Gbox_merge_df$MaxExprCell <- factor(Gbox_merge_df$MaxExprCell, levels = c("C1", "C2", "C3", "C4", "C5"))
Gbox_merge_df$MuRatio <- Gbox_merge_df$Mu / Gbox_gene_num
Gbox_merge_df$SdRatio <- Gbox_merge_df$Sd / Gbox_gene_num
Gbox_merge_df$Label <- sprintf("p=%.1E", Gbox_merge_df$Pval)
Gbox_merge_df$Label[is.na(Gbox_merge_df$Pval)] <- NA
p <- ggplot(Gbox_merge_df) +
  geom_bar(stat = "identity", mapping = aes(x=MaxExprCell, y=MuRatio, fill=Data), position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(mapping = aes(x=MaxExprCell, ymin=MuRatio-SdRatio, ymax=MuRatio+SdRatio, group=Data), color="black", size=0.5, position = position_dodge(width = 0.7), width = 0.4) +
  geom_text(mapping = aes(x=MaxExprCell, y=MuRatio+SdRatio+0.02, label=Label), color="black", size=1.4) +
  scale_fill_manual(values = c("Gbox"="red", "Random"="grey70")) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4), labels = c("0%", "20%", "40%")) +
  labs(x="Max. Expr. Cell", y="Ratio") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size=5),
    title = element_text(family="ArialMT", size=5),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(family="ArialMT", size=5),
    legend.position = c(0.8, 0.8),
    panel.grid = element_blank()
  )
ggsave(file.path(f_out, "Gbox.ExprRatio.pdf"), p, height = 4, width = 5, units = "cm", colormodel = "cmyk")



MYB_G4_max_expr_cell <- inner_join(MYB_G4, max_expr_cell)
MYB_G4_max_expr_cell_info <- MYB_G4_max_expr_cell %>% group_by(MaxExprCell) %>% summarise(Data="MYB_G4", Indx=1, GeneNum=n())
MYB_G4_gene_num <- nrow(MYB_G4_max_expr_cell)
all_gene_num <- nrow(max_expr_cell)
MYB_G4_rand_max_expr_cell_info <- data.frame()
for (i in 1:10000) {
  tmp_max_expr_cell <- max_expr_cell[sample(1:all_gene_num, MYB_G4_gene_num, replace = FALSE),]
  tmp_max_expr_cell <- tmp_max_expr_cell %>% group_by(MaxExprCell) %>% summarise(Data="Random", Indx=i, GeneNum=n())
  MYB_G4_rand_max_expr_cell_info <- rbind(MYB_G4_rand_max_expr_cell_info, tmp_max_expr_cell)
}
MYB_G4_rand_max_expr_cell_info <- left_join(MYB_G4_rand_max_expr_cell_info, data.frame(MaxExprCell=MYB_G4_max_expr_cell_info$MaxExprCell, MYB_G4GeneNum=MYB_G4_max_expr_cell_info$GeneNum))
MYB_G4_rand_max_expr_cell_df <- MYB_G4_rand_max_expr_cell_info %>% group_by(MaxExprCell, Data) %>% summarise(Mu=mean(GeneNum), Sd=sd(GeneNum), Pval=min(mean(GeneNum<MYB_G4GeneNum), mean(GeneNum>MYB_G4GeneNum))*2)
MYB_G4_max_expr_cell_df <- data.frame(MaxExprCell=MYB_G4_max_expr_cell_info$MaxExprCell, Data="MYB_G4", Mu=MYB_G4_max_expr_cell_info$GeneNum, Sd=0, Pval=NA)
MYB_G4_merge_df <- rbind(as.data.frame(MYB_G4_max_expr_cell_df), as.data.frame(MYB_G4_rand_max_expr_cell_df))
MYB_G4_merge_df <- MYB_G4_merge_df[!is.na(MYB_G4_merge_df$MaxExprCell),]
write_tsv(MYB_G4_merge_df, file.path(f_out, "MYB_G4.MaxExprCell.Expr.test.tsv"))
MYB_G4_merge_df$Data <- factor(MYB_G4_merge_df$Data, levels = c("MYB_G4", "Random"))
MYB_G4_merge_df$MaxExprCell <- factor(MYB_G4_merge_df$MaxExprCell, levels = c("C1", "C2", "C3", "C4", "C5"))
MYB_G4_merge_df$MuRatio <- MYB_G4_merge_df$Mu / MYB_G4_gene_num
MYB_G4_merge_df$SdRatio <- MYB_G4_merge_df$Sd / MYB_G4_gene_num
MYB_G4_merge_df$Label <- sprintf("p=%.1E", MYB_G4_merge_df$Pval)
MYB_G4_merge_df$Label[is.na(MYB_G4_merge_df$Pval)] <- NA
p <- ggplot(MYB_G4_merge_df) +
  geom_bar(stat = "identity", mapping = aes(x=MaxExprCell, y=MuRatio, fill=Data), position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(mapping = aes(x=MaxExprCell, ymin=MuRatio-SdRatio, ymax=MuRatio+SdRatio, group=Data), color="black", size=0.5, position = position_dodge(width = 0.7), width = 0.4) +
  geom_text(mapping = aes(x=MaxExprCell, y=MuRatio+SdRatio+0.02, label=Label), color="black", size=1.4) +
  scale_fill_manual(values = c("MYB_G4"="red", "Random"="grey70")) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4), labels = c("0%", "20%", "40%")) +
  labs(x="Max. Expr. Cell", y="Ratio") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size=5),
    title = element_text(family="ArialMT", size=5),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(family="ArialMT", size=5),
    legend.position = c(0.8, 0.8),
    panel.grid = element_blank()
  )
ggsave(file.path(f_out, "MYB_G4.ExprRatio.pdf"), p, height = 4, width = 5, units = "cm", colormodel = "cmyk")


kmer_8_graph_enrich_cluster_motif_pos <- read_delim("analysis_v3/ATAC/ATAC_kmer_cluster2gene/kmer_clu.filter.expr.SourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
kmer_8_graph_enrich_cluster_motif_pos$IsTCP <- FALSE
kmer_8_graph_enrich_cluster_motif_pos$IsTCP[grep("TGGGC[CT]", kmer_8_graph_enrich_cluster_motif_pos$Motif)] <- TRUE
kmer_8_graph_enrich_cluster_motif_pos$IsTCP[grep("TGGGC[CT]", kmer_8_graph_enrich_cluster_motif_pos$RevMotif)] <- TRUE
kmer_8_graph_enrich_cluster_motif_pos$IsTCPLike <- FALSE
kmer_8_graph_enrich_cluster_motif_pos$IsTCPLike[grep("TAGGGT", kmer_8_graph_enrich_cluster_motif_pos$Motif)] <- TRUE
kmer_8_graph_enrich_cluster_motif_pos$IsTCPLike[grep("TAGGGT", kmer_8_graph_enrich_cluster_motif_pos$RevMotif)] <- TRUE

TCP_gids <- unique(kmer_8_graph_enrich_cluster_motif_pos$Gid[kmer_8_graph_enrich_cluster_motif_pos$IsTCP])
TCP_like_gids <- unique(kmer_8_graph_enrich_cluster_motif_pos$Gid[kmer_8_graph_enrich_cluster_motif_pos$IsTCPLike])

TCP_WT_enrich_max_expr_cell <- inner_join(data.frame(Gid=TCP_gids), max_expr_cell)
TCP_WT_enrich_max_expr_cell_info <- TCP_WT_enrich_max_expr_cell %>% group_by(MaxExprCell) %>% summarise(Data="TCP", Indx=1, GeneNum=n())
TCP_gene_num <- nrow(TCP_WT_enrich_max_expr_cell)
all_gene_num <- nrow(max_expr_cell)
TCP_rand_max_expr_cell_info <- data.frame()
for (i in 1:10000) {
  tmp_max_expr_cell <- max_expr_cell[sample(1:all_gene_num, TCP_gene_num, replace = FALSE),]
  tmp_max_expr_cell <- tmp_max_expr_cell %>% group_by(MaxExprCell) %>% summarise(Data="Random", Indx=i, GeneNum=n())
  TCP_rand_max_expr_cell_info <- rbind(TCP_rand_max_expr_cell_info, tmp_max_expr_cell)
}
TCP_rand_max_expr_cell_info <- left_join(TCP_rand_max_expr_cell_info, data.frame(MaxExprCell=TCP_WT_enrich_max_expr_cell_info$MaxExprCell, TCPGeneNum=TCP_WT_enrich_max_expr_cell_info$GeneNum))
TCP_rand_max_expr_cell_df <- TCP_rand_max_expr_cell_info %>% group_by(MaxExprCell, Data) %>% summarise(Mu=mean(GeneNum), Sd=sd(GeneNum), Pval=min(mean(GeneNum<TCPGeneNum), mean(GeneNum>TCPGeneNum))*2)
TCP_max_expr_cell_df <- data.frame(MaxExprCell=TCP_WT_enrich_max_expr_cell_info$MaxExprCell, Data="TCP", Mu=TCP_WT_enrich_max_expr_cell_info$GeneNum, Sd=0, Pval=NA)
TCP_merge_df <- rbind(as.data.frame(TCP_max_expr_cell_df), as.data.frame(TCP_rand_max_expr_cell_df))
TCP_merge_df <- TCP_merge_df[!is.na(TCP_merge_df$MaxExprCell),]
write_tsv(TCP_merge_df, file.path(f_out, "TCP.WT_enrich.MaxExprCell.Expr.test.tsv"))
TCP_merge_df$Data <- factor(TCP_merge_df$Data, levels = c("TCP", "Random"))
TCP_merge_df$MaxExprCell <- factor(TCP_merge_df$MaxExprCell, levels = c("C1", "C2", "C3", "C4", "C5"))
TCP_merge_df$MuRatio <- TCP_merge_df$Mu / TCP_gene_num
TCP_merge_df$SdRatio <- TCP_merge_df$Sd / TCP_gene_num
TCP_merge_df$Label <- sprintf("p=%.1E", TCP_merge_df$Pval)
TCP_merge_df$Label[is.na(TCP_merge_df$Pval)] <- NA
p <- ggplot(TCP_merge_df) +
  geom_bar(stat = "identity", mapping = aes(x=MaxExprCell, y=MuRatio, fill=Data), position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(mapping = aes(x=MaxExprCell, ymin=MuRatio-SdRatio, ymax=MuRatio+SdRatio, group=Data), color="black", size=0.5, position = position_dodge(width = 0.7), width = 0.4) +
  geom_text(mapping = aes(x=MaxExprCell, y=MuRatio+SdRatio+0.02, label=Label), color="black", size=1.4) +
  scale_fill_manual(values = c("TCP"="red", "Random"="grey70")) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4), labels = c("0%", "20%", "40%")) +
  labs(x="Max. Expr. Cell", y="Ratio") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size=5),
    title = element_text(family="ArialMT", size=5),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(family="ArialMT", size=5),
    legend.position = c(0.8, 0.8),
    panel.grid = element_blank()
  )
ggsave(file.path(f_out, "ATAC_DE_peak_target_gene.TCP.ExprRatio.pdf"), p, height = 4, width = 5, units = "cm", colormodel = "cmyk")



TCP_like_WT_enrich_max_expr_cell <- inner_join(data.frame(Gid=TCP_like_gids), max_expr_cell)
TCP_like_WT_enrich_max_expr_cell_info <- TCP_like_WT_enrich_max_expr_cell %>% group_by(MaxExprCell) %>% summarise(Data="TCP_like", Indx=1, GeneNum=n())
TCP_like_gene_num <- nrow(TCP_like_WT_enrich_max_expr_cell)
all_gene_num <- nrow(max_expr_cell)
TCP_like_rand_max_expr_cell_info <- data.frame()
for (i in 1:10000) {
  tmp_max_expr_cell <- max_expr_cell[sample(1:all_gene_num, TCP_like_gene_num, replace = FALSE),]
  tmp_max_expr_cell <- tmp_max_expr_cell %>% group_by(MaxExprCell) %>% summarise(Data="Random", Indx=i, GeneNum=n())
  TCP_like_rand_max_expr_cell_info <- rbind(TCP_like_rand_max_expr_cell_info, tmp_max_expr_cell)
}
TCP_like_rand_max_expr_cell_info <- left_join(TCP_like_rand_max_expr_cell_info, data.frame(MaxExprCell=TCP_like_WT_enrich_max_expr_cell_info$MaxExprCell, TCP_likeGeneNum=TCP_like_WT_enrich_max_expr_cell_info$GeneNum))
TCP_like_rand_max_expr_cell_df <- TCP_like_rand_max_expr_cell_info %>% group_by(MaxExprCell, Data) %>% summarise(Mu=mean(GeneNum), Sd=sd(GeneNum), Pval=min(mean(GeneNum<TCP_likeGeneNum), mean(GeneNum>TCP_likeGeneNum))*2)
TCP_like_max_expr_cell_df <- data.frame(MaxExprCell=TCP_like_WT_enrich_max_expr_cell_info$MaxExprCell, Data="TCP_like", Mu=TCP_like_WT_enrich_max_expr_cell_info$GeneNum, Sd=0, Pval=NA)
TCP_like_merge_df <- rbind(as.data.frame(TCP_like_max_expr_cell_df), as.data.frame(TCP_like_rand_max_expr_cell_df))
TCP_like_merge_df <- TCP_like_merge_df[!is.na(TCP_like_merge_df$MaxExprCell),]
write_tsv(TCP_like_merge_df, file.path(f_out, "TCP_like.WT_enrich.MaxExprCell.Expr.test.tsv"))
TCP_like_merge_df$Data <- factor(TCP_like_merge_df$Data, levels = c("TCP_like", "Random"))
TCP_like_merge_df$MaxExprCell <- factor(TCP_like_merge_df$MaxExprCell, levels = c("C1", "C2", "C3", "C4", "C5"))
TCP_like_merge_df$MuRatio <- TCP_like_merge_df$Mu / TCP_like_gene_num
TCP_like_merge_df$SdRatio <- TCP_like_merge_df$Sd / TCP_like_gene_num
TCP_like_merge_df$Label <- sprintf("p=%.1E", TCP_like_merge_df$Pval)
TCP_like_merge_df$Label[is.na(TCP_like_merge_df$Pval)] <- NA
p <- ggplot(TCP_like_merge_df) +
  geom_bar(stat = "identity", mapping = aes(x=MaxExprCell, y=MuRatio, fill=Data), position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(mapping = aes(x=MaxExprCell, ymin=MuRatio-SdRatio, ymax=MuRatio+SdRatio, group=Data), color="black", size=0.5, position = position_dodge(width = 0.7), width = 0.4) +
  geom_text(mapping = aes(x=MaxExprCell, y=MuRatio+SdRatio+0.02, label=Label), color="black", size=1.4) +
  scale_fill_manual(values = c("TCP_like"="red", "Random"="grey70")) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4), labels = c("0%", "20%", "40%")) +
  labs(x="Max. Expr. Cell", y="Ratio") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size=5),
    title = element_text(family="ArialMT", size=5),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(family="ArialMT", size=5),
    legend.position = c(0.8, 0.8),
    panel.grid = element_blank()
  )
ggsave(file.path(f_out, "ATAC_DE_peak_target_gene.TCP_like.ExprRatio.pdf"), p, height = 4, width = 5, units = "cm", colormodel = "cmyk")





Gbox_topGO <- read_delim("analysis_v3/ATAC/ATAC_DE_peak2C3/Gbox_WT_enrich_GO/topGO.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

TF_activate_gids <- unique(unlist(strsplit(Gbox_topGO$Gids[Gbox_topGO$GO.ID=="GO:0006355"], "[ ]")))

DE_TF_scRNA_df <- left_join(data.frame(Gid=TF_activate_gids), cluster_tpm)
DE_TF_TP_df <- left_join(data.frame(Gid=TF_activate_gids), CPM_TP)
DE_TF_public_TP_df <- left_join(data.frame(Gid=TF_activate_gids), public_RNA_seq)

DE_TF_scRNA_mat <- as.matrix(DE_TF_scRNA_df[,2:ncol(DE_TF_scRNA_df)])
rownames(DE_TF_scRNA_mat) <- DE_TF_scRNA_df$Gid
norm_DE_TF_scRNA_mat <- DE_TF_scRNA_mat / apply(DE_TF_scRNA_mat, 1, function(x){return(max(x, na.rm = T))})

DE_TF_TP_mat <- as.matrix(DE_TF_TP_df[,2:ncol(DE_TF_TP_df)])
rownames(DE_TF_TP_mat) <- DE_TF_TP_df$Gid
norm_DE_TF_TP_mat <- DE_TF_TP_mat / apply(DE_TF_TP_mat, 1, max)

DE_TF_public_TP_mat <- as.matrix(DE_TF_public_TP_df[,2:ncol(DE_TF_public_TP_df)])
rownames(DE_TF_public_TP_mat) <- DE_TF_public_TP_df$Gid
DE_TF_public_TP_mat[rowSums(DE_TF_public_TP_mat)==0,] <- NA
norm_DE_TF_public_TP_mat <- DE_TF_public_TP_mat / apply(DE_TF_public_TP_mat, 1, function(x){return(max(x, na.rm = T))})

interest_gene_df <- left_join(data.frame(Gid=TF_activate_gids), interest_gene)
interest_gene_df$Label <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$GeneSymbol)
interest_gene_df$Label[!is.na(interest_gene_df$Family)] <- sprintf("%s (%s: %s)", interest_gene_df$Gid, interest_gene_df$Family, interest_gene_df$GeneSymbol)[!is.na(interest_gene_df$Family)]
interest_gene_df$Label[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- sprintf("%s (%s)", interest_gene_df$Gid, interest_gene_df$Family)[!is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]
interest_gene_df$Label[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)] <- interest_gene_df$Gid[is.na(interest_gene_df$Family) & is.na(interest_gene_df$GeneSymbol)]

subindx <- which(!is.na(interest_gene_df$GeneSymbol))
label_li <- interest_gene_df$Label[subindx]

scRNA_dist <- dist(norm_DE_TF_scRNA_mat)
TP_dist <- dist(norm_DE_TF_TP_mat)
public_TP_dist <- dist(DE_TF_public_TP_mat)
scRNA_dist[is.na(scRNA_dist)] <- max(scRNA_dist, na.rm = T)
TP_dist[is.na(TP_dist)] <- max(TP_dist, na.rm = T)
public_TP_dist[is.na(public_TP_dist)] <- max(public_TP_dist, na.rm = T)
merge_dist <- 1*scRNA_dist + 0.1*TP_dist
clu <- hclust(merge_dist, method = "ward.D2")

scRNA_ht <- Heatmap(
  norm_DE_TF_scRNA_mat, 
  name="Norm.Expr.",
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  cluster_columns = FALSE, 
  row_order = clu$order,
  show_row_names = F, 
  column_split = c(rep("WT", 5), rep("FL", 5)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

TP_ht <- Heatmap(
  norm_DE_TF_TP_mat, 
  name="Norm.Expr.",
  # col = colorRamp2(c(0, 1), c("white", "red")),
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  row_order = clu$order,
  cluster_columns = FALSE, 
  show_row_names = F, 
  column_split = c(rep("WT", 9), rep("FL", 9)),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

public_TP_ht <- Heatmap(
  norm_DE_TF_public_TP_mat, 
  name="Norm.Expr.",
  # col = colorRamp2(c(0, 1), c("white", "red")),
  col = colorRamp2(c(0, 0.5, 1), c("magenta", "black", "yellow")),
  row_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  row_order = clu$order,
  cluster_columns = FALSE, 
  show_row_names = F, 
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

ha <- rowAnnotation(
  link = row_anno_link(
    at = subindx,
    labels = label_li,
    labels_gp=gpar(fontsize = 6),
    lines_gp=gpar(lwd = 0.5),
    link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)))

inche_cm <- 2.54
pdf(file.path(f_out, "ATAC_DE_peak_target_gene.Gbox.TF_GO.ht.pdf"), width=15.5/inche_cm, height=8/inche_cm)
print(scRNA_ht+TP_ht+public_TP_ht+ha)
dev.off()
