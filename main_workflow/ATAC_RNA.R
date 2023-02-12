library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(Gviz)
library(GenomicFeatures)
library(ComplexHeatmap)
library(circlize)

setwd("D:/cottonSingleCell")
ATAC_gene_expr <- read_delim("stat/ATAC_MACS2/ATAC.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_gene_expr <- ATAC_gene_expr[,c("Gid", "WT", "FL")]
ATAC_gene_DE <- read_delim("stat/ATAC_MACS2/ATAC.gene.DE.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RNA_tpm <-  read_delim("stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_peak_DE <- read_delim("stat/ATAC_MACS2/ATAC.peak.DE.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
refGene <- read_delim("Ghirsutumv1.1_gene_model.gene.bed6", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(refGene) <- c("Chrom", "Start", "End", "Name", "Score", "Strand")
WT_peak_DE <- ATAC_peak_DE[ATAC_peak_DE$Tag=="WT",]
FL_peak_DE <- ATAC_peak_DE[ATAC_peak_DE$Tag=="FL",]
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
refGene$FL_peak_DE <- apply(refGene, 1, function(line, expand=1000){
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
  tmp_FL_peak_DE <- FL_peak_DE[FL_peak_DE$Chrom==gene_chrom,]
  res <- any(!((tmp_FL_peak_DE$Start >= end_pos) | (start_pos >= tmp_FL_peak_DE$End)))
  return(res)
})
RNA_tpm$IsExpr <- apply(RNA_tpm[,2:ncol(RNA_tpm)], 1, function(x){return(max(x, na.rm = T)>1)})
refGene$IsExpr <- refGene$Name %in% RNA_tpm$Gid[RNA_tpm$IsExpr]
refGene$expr_WT_peak_DE <- refGene$WT_peak_DE & refGene$IsExpr
refGene$expr_FL_peak_DE <- refGene$FL_peak_DE & refGene$IsExpr

wt_target_gene <- refGene$Name[refGene$expr_WT_peak_DE]


melt_RNA_tpm <- melt(RNA_tpm, id.vars="Gid", variable.name="RNA_Cluster", value.name = "RNA_tpm")
melt_ATAC_expr <- melt(ATAC_gene_expr, id.vars="Gid", variable.name="ATAC_Sample", value.name = "ATAC_Expr")

ATAC_RNA_expr <- inner_join(melt_RNA_tpm, melt_ATAC_expr)
ATAC_RNA_expr <- ATAC_RNA_expr[(ATAC_RNA_expr$RNA_tpm+ATAC_RNA_expr$ATAC_Expr)>0,]
ATAC_RNA_expr <- na.omit(ATAC_RNA_expr)

p <- ggplot(ATAC_RNA_expr, aes(y=ATAC_Expr, x=RNA_tpm)) +
  geom_point(size=0.6, alpha=0.1) +
  facet_grid(ATAC_Sample~RNA_Cluster) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        # axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        # legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        panel.grid = element_blank())
ggsave("stat/ATAC_MACS2/atac_RNA.png", p, width = 50, height = 10, units = "cm", dpi=600)

hist(log10(rowMeans(ATAC_gene_expr[,-1])), breaks = 50)

ATAC_log2FC <- data.frame(Gid=ATAC_gene_expr$Gid, ATAC_log2FC=log2((ATAC_gene_expr$FL+1e-3)/(ATAC_gene_expr$WT+1e-3)))
ATAC_log2FC <- ATAC_log2FC[(ATAC_gene_expr$FL>0.1)|(ATAC_gene_expr$WT>0.1),]
RNA_WT_mat <- as.matrix(RNA_tpm[,c("WT.C1", "WT.C2", "WT.C4", "WT.C5")])
RNA_FL_mat <- as.matrix(RNA_tpm[,c("FL.C1", "FL.C2", "FL.C4", "FL.C5")])
colnames(RNA_WT_mat) <- c("C1", "C2", "C4", "C5")
colnames(RNA_FL_mat) <- c("C1", "C2", "C4", "C5")
hist(log10(rowMeans(RNA_WT_mat)), breaks = 50)

RNA_max_mat <- pmax(RNA_FL_mat, RNA_WT_mat)

RNA_log2FC <- log2((RNA_FL_mat+1e-3)/(RNA_WT_mat+1e-3))
RNA_log2FC[RNA_max_mat<1] <- NA
RNA_log2FC <- cbind(data.frame(Gid=RNA_tpm$Gid), RNA_log2FC)
melt_RNA_log2FC <- melt(RNA_log2FC, id.vars="Gid", variable.name="RNA_Cluster", value.name = "RNA_log2FC")
melt_RNA_log2FC <- na.omit(melt_RNA_log2FC)
ATAC_RNA_log2FC <- inner_join(ATAC_log2FC, melt_RNA_log2FC)
names(ATAC_RNA_log2FC)
ATAC_RNA_log2FC$ATAC_DE <- ATAC_RNA_log2FC$Gid %in% ATAC_gene_DE$Gid
x <- ATAC_RNA_log2FC[ATAC_RNA_log2FC$ATAC_DE,]
ATAC_RNA_log2FC %>% group_by(RNA_Cluster, ATAC_DE) %>% summarise(R=cor(ATAC_log2FC, RNA_log2FC, method = "spearman"))

ATAC_RNA_log2FC$ATAC_log2FC[ATAC_RNA_log2FC$ATAC_log2FC < -5] <- -5
ATAC_RNA_log2FC$ATAC_log2FC[ATAC_RNA_log2FC$ATAC_log2FC > 5] <- 5
ATAC_RNA_log2FC$RNA_log2FC[ATAC_RNA_log2FC$RNA_log2FC < -10] <- -10
ATAC_RNA_log2FC$RNA_log2FC[ATAC_RNA_log2FC$RNA_log2FC > 10] <- 10
p <- ggplot(mapping=aes(y=ATAC_log2FC, x=RNA_log2FC)) +
  geom_point(data=ATAC_RNA_log2FC[!ATAC_RNA_log2FC$ATAC_DE,], size=0.6, alpha=0.1, color="grey") +
  geom_point(data=ATAC_RNA_log2FC[ATAC_RNA_log2FC$ATAC_DE & ATAC_RNA_log2FC$ATAC_log2FC<0,], size=0.5, color="red", alpha=0.5) +
  geom_point(data=ATAC_RNA_log2FC[ATAC_RNA_log2FC$ATAC_DE & ATAC_RNA_log2FC$ATAC_log2FC>0,], size=0.5, color="blue", alpha=0.5) +
  geom_hline(yintercept = 0, size=0.3, lty=2) +
  geom_vline(xintercept = 0, size=0.3, lty=2) +
  # scale_fill_continuous(trans = scales::log10_trans()) +
  facet_grid(~RNA_Cluster) +
  # lims(x=c(-15, 15), y=c(-10, 10)) +
  labs(x="log2(FL/WT) for scRNA-seq", y="log2(FL/WT) for scATAC-seq", fill="Number") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        # legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        panel.grid = element_blank())
ggsave("stat/ATAC_MACS2/atac_RNA.log2FC.pdf", p, width = 15, height = 4.8, units = "cm")

gene_DE_num <- ATAC_gene_DE %>% group_by(Tag) %>% summarise(Num=n())
gene_DE_num$Stat <- "Gene"
peak_DE_num <- ATAC_peak_DE %>% group_by(Tag) %>% summarise(Num=n())
peak_DE_num$Stat <- "Peak"
peak_target_gene_DE_num <- data.frame(
  Tag=c("WT", "FL"),
  Num=c(sum(refGene$WT_peak_DE),  sum(refGene$FL_peak_DE)),
  Stat=c("Target gene", "Target gene")
)

DE_num <- rbind(gene_DE_num, peak_DE_num, peak_target_gene_DE_num)
DE_num$Tag <- factor(DE_num$Tag, levels = c("WT", "FL"), labels = c("Down", "Up"))
DE_num$Stat <- factor(DE_num$Stat, levels = c("Peak", "Target gene", "Gene"))
p <- ggplot(DE_num, aes(x=Tag, y=Num, fill=Tag)) +
  geom_bar(stat="identity") +
  geom_point(aes(y=Num*1.05), color=NA) +
  geom_text(aes(label=Num), size=1.7, angle=-90, vjust=-0.2, color="black") +
  coord_flip() +
  theme_bw() +
  labs(y="Number") +
  facet_grid(~Stat, scales = "free_x") +
  scale_fill_manual(values = c("Down"="blue", "Up"="red")) +
  theme(text = element_text(family="ArialMT", size=7),
        title = element_text(family="ArialMT", size=8),
        axis.text = element_text(color = "black"),
        # axis.text.y = element_text(angle = 45, hjust = 1, size=7),
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()
  )
p
ggsave("stat/ATAC_MACS2/atac_DE_num.pdf", p, width = 6, height = 2.5, units = "cm")

gene_DE <- ATAC_gene_DE
gene_DE <- left_join(gene_DE, RNA_tpm)
gene_DE <- gene_DE[order(gene_DE$Log2FC, decreasing = TRUE),]
gene_DE_WT <- as.matrix(gene_DE[,c("WT.C1", "WT.C2", "WT.C4", "WT.C5")])
gene_DE_FL <- as.matrix(gene_DE[,c("FL.C1", "FL.C2", "FL.C4", "FL.C5")])
gene_DE_log2FC <- log2((gene_DE_FL+1e-3) / (gene_DE_WT+1e-3))

gene_DE_AVE_tpm <- (gene_DE_WT+gene_DE_FL) / 2

# gene_DE_log2FC[(gene_DE_AVE_tpm) < 0.1] <- NA
colnames(gene_DE_log2FC) <- c("C1", "C2", "C4", "C5")
colnames(gene_DE_AVE_tpm) <- c("C1", "C2", "C4", "C5")

gene_DE$Log2FC[is.infinite(gene_DE$Log2FC) & gene_DE$Log2FC<0] <- -15
rownames(gene_DE_log2FC) <- gene_DE$Gid
rownames(gene_DE_AVE_tpm) <- gene_DE$Gid

DiffExpr <- read_delim("stat/diff_expr/DiffExpr.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

log_gene_DE_AVE_tpm <- log10(gene_DE_AVE_tpm)
log_gene_DE_AVE_tpm[is.infinite(log_gene_DE_AVE_tpm)] <- NA

ha <- rowAnnotation(
  ATAC_log2FC=anno_barplot(gene_DE$Log2FC, width=unit(5, "mm"), gp=gpar(fontsize = 5)),
  annotation_name_gp=gpar(fontsize = 5),
  annotation_legend_param=list(
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm"),
    title_gp = gpar(fontsize = 6),
    labels_gp=gpar(fontsize = 5)
  )
)

gene_DE_log2FC[is.na(log_gene_DE_AVE_tpm)] <- NA
ht <- Heatmap(gene_DE_log2FC, name="scRNA-seq log2FC",
              column_title = "log2FC",
              col =  colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_dend_width = unit(5, "mm"),
              column_dend_height = unit(5, "mm"),
              show_row_dend = TRUE,
              show_row_names = TRUE,
              row_title = NULL,
              row_title_gp = gpar(fontsize = 6),
              column_title_gp = gpar(fontsize = 6),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              split = gene_DE$Log2FC<0,
              width=unit(1, "cm"),
              left_annotation = ha,
              cell_fun = function(j, i, x, y, w, h, col) {
                gid <- rownames(gene_DE_log2FC)[i]
                clu <- colnames(gene_DE_log2FC)[j]
                if(sum(DiffExpr$Gid == gid & DiffExpr$Cluster==clu)){
                  label = "*"
                }else{
                  label = ""
                }
                grid.text(label, x, y)
              },
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 6),
                title_gp = gpar(fontsize = 6),
                grid_width = unit(2, "mm"),
                grid_height = unit(2, "mm"))
)

ht2 <- Heatmap(log_gene_DE_AVE_tpm, name="scRNA-seq log10 average TPM",
               column_title = "log10 ave. TPM",
               col =  colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_dend_width = unit(5, "mm"),
              column_dend_height = unit(5, "mm"),
              show_row_dend = TRUE,
              show_row_names = TRUE,
              row_title = NULL,
              row_title_gp = gpar(fontsize = 6),
              column_title_gp = gpar(fontsize = 6),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              split = gene_DE$Log2FC<0,
              width=unit(1, "cm"),
              cell_fun = function(j, i, x, y, w, h, col) {
                gid <- rownames(gene_DE_log2FC)[i]
                clu <- colnames(gene_DE_log2FC)[j]
                if(sum(DiffExpr$Gid == gid & DiffExpr$Cluster==clu)){
                  label = "*"
                }else{
                  label = ""
                }
                grid.text(label, x, y)
              },
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 6),
                title_gp = gpar(fontsize = 6),
                grid_width = unit(2, "mm"),
                grid_height = unit(2, "mm"))
)

inche_cm <- 2.54
pdf("stat/ATAC_MACS2/ATAC.RNA.DE.ht.pdf", width=10/inche_cm, height=40/inche_cm)
print(ht + ht2)
dev.off()

mask_log_gene_DE_AVE_tpm <- log_gene_DE_AVE_tpm
mask_log_gene_DE_AVE_tpm[mask_log_gene_DE_AVE_tpm<0] <- NA
mask_gene_DE_log2FC <- gene_DE_log2FC
mask_gene_DE_log2FC[is.na(mask_log_gene_DE_AVE_tpm)] <- NA

low_expr_indx <- rowSums(!is.na(mask_log_gene_DE_AVE_tpm)) == 0
log_gene_DE_AVE_tpm <- log_gene_DE_AVE_tpm[!low_expr_indx,]
gene_DE_log2FC <- gene_DE_log2FC[!low_expr_indx>0,]
mask_log_gene_DE_AVE_tpm <- mask_log_gene_DE_AVE_tpm[!low_expr_indx,]
mask_gene_DE_log2FC <- mask_gene_DE_log2FC[!low_expr_indx>0,]


ha <- rowAnnotation(
  ATAC_log2FC=anno_barplot(gene_DE$Log2FC[!low_expr_indx], width=unit(5, "mm"), gp=gpar(fontsize = 5)),
  annotation_name_gp=gpar(fontsize = 5),
  annotation_legend_param=list(
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm"),
    title_gp = gpar(fontsize = 6),
    labels_gp=gpar(fontsize = 5)
  )
)
ht <- Heatmap(gene_DE_log2FC, name="scRNA-seq log2FC",
              column_title = "log2FC",
              col =  colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_dend_width = unit(5, "mm"),
              column_dend_height = unit(5, "mm"),
              show_row_dend = TRUE,
              show_row_names = TRUE,
              row_title = NULL,
              row_title_gp = gpar(fontsize = 6),
              column_title_gp = gpar(fontsize = 6),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              split = (gene_DE$Log2FC<0)[!low_expr_indx],
              width=unit(1, "cm"),
              left_annotation = ha,
              cell_fun = function(j, i, x, y, w, h, col) {
                gid <- rownames(gene_DE_log2FC)[i]
                clu <- colnames(gene_DE_log2FC)[j]
                if(sum(DiffExpr$Gid == gid & DiffExpr$Cluster==clu)){
                  label = "*"
                }else{
                  label = ""
                }
                grid.text(label, x, y)
              },
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 6),
                title_gp = gpar(fontsize = 6),
                grid_width = unit(2, "mm"),
                grid_height = unit(2, "mm"))
)

ht2 <- Heatmap(log_gene_DE_AVE_tpm, name="scRNA-seq log10 average TPM",
               column_title = "log10 ave. TPM",
               col =  colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
               row_names_gp = gpar(fontsize = 6),
               column_names_gp = gpar(fontsize = 6),
               row_dend_width = unit(5, "mm"),
               column_dend_height = unit(5, "mm"),
               show_row_dend = TRUE,
               show_row_names = TRUE,
               row_title = NULL,
               row_title_gp = gpar(fontsize = 6),
               column_title_gp = gpar(fontsize = 6),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               split = (gene_DE$Log2FC<0)[!low_expr_indx],
               width=unit(1, "cm"),
               cell_fun = function(j, i, x, y, w, h, col) {
                 gid <- rownames(gene_DE_log2FC)[i]
                 clu <- colnames(gene_DE_log2FC)[j]
                 if(sum(DiffExpr$Gid == gid & DiffExpr$Cluster==clu)){
                   label = "*"
                 }else{
                   label = ""
                 }
                 grid.text(label, x, y)
               },
               heatmap_legend_param = list(
                 labels_gp = gpar(fontsize = 6),
                 title_gp = gpar(fontsize = 6),
                 grid_width = unit(2, "mm"),
                 grid_height = unit(2, "mm"))
)

ht3 <- Heatmap(mask_gene_DE_log2FC, name="scRNA-seq log2FC",
              column_title = "log2FC",
              col =  colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_dend_width = unit(5, "mm"),
              column_dend_height = unit(5, "mm"),
              show_row_dend = TRUE,
              show_row_names = TRUE,
              row_title = NULL,
              row_title_gp = gpar(fontsize = 6),
              column_title_gp = gpar(fontsize = 6),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              split = (gene_DE$Log2FC<0)[!low_expr_indx],
              width=unit(1, "cm"),
              left_annotation = ha,
              cell_fun = function(j, i, x, y, w, h, col) {
                gid <- rownames(gene_DE_log2FC)[i]
                clu <- colnames(gene_DE_log2FC)[j]
                if(sum(DiffExpr$Gid == gid & DiffExpr$Cluster==clu)){
                  label = "*"
                }else{
                  label = ""
                }
                grid.text(label, x, y)
              },
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 6),
                title_gp = gpar(fontsize = 6),
                grid_width = unit(2, "mm"),
                grid_height = unit(2, "mm"))
)

ht4 <- Heatmap(mask_log_gene_DE_AVE_tpm, name="scRNA-seq log10 average TPM",
               column_title = "log10 ave. TPM",
               col =  colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
               row_names_gp = gpar(fontsize = 6),
               column_names_gp = gpar(fontsize = 6),
               row_dend_width = unit(5, "mm"),
               column_dend_height = unit(5, "mm"),
               show_row_dend = TRUE,
               show_row_names = TRUE,
               row_title = NULL,
               row_title_gp = gpar(fontsize = 6),
               column_title_gp = gpar(fontsize = 6),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               split = (gene_DE$Log2FC<0)[!low_expr_indx],
               width=unit(1, "cm"),
               cell_fun = function(j, i, x, y, w, h, col) {
                 gid <- rownames(gene_DE_log2FC)[i]
                 clu <- colnames(gene_DE_log2FC)[j]
                 if(sum(DiffExpr$Gid == gid & DiffExpr$Cluster==clu)){
                   label = "*"
                 }else{
                   label = ""
                 }
                 grid.text(label, x, y)
               },
               heatmap_legend_param = list(
                 labels_gp = gpar(fontsize = 6),
                 title_gp = gpar(fontsize = 6),
                 grid_width = unit(2, "mm"),
                 grid_height = unit(2, "mm"))
)

inche_cm <- 2.54
pdf("stat/ATAC_MACS2/ATAC.RNA.DE.expr.ht.pdf", width=20/inche_cm, height=30/inche_cm)
print(ht + ht2 + ht3 + ht4)
dev.off()

ref_gene <- read_delim("Ghirsutumv1.1_gene_model.gene.bed6", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(ref_gene) <- c("Chrom", "Start", "End", "Gid", "Score", "Strand")
ref_gene <- ref_gene[,c("Chrom", "Start", "End", "Gid", "Strand")]

find_cloest_gene <- function(line){
  chrom <- line["Chrom"]
  start <- as.integer(line["Start"])
  end <- as.integer(line["End"])
  gene_with_same_chrom <- ref_gene[ref_gene$Chrom==chrom,]
  gene_with_same_chrom$Pos <- "Hit"
  gene_with_same_chrom$Pos[gene_with_same_chrom$End<start] <- "Left"
  gene_with_same_chrom$Pos[gene_with_same_chrom$Start>end] <- "Right"
  gene_with_same_chrom$Dist <- 0
  gene_with_same_chrom$Dist[gene_with_same_chrom$Pos=="Left"] <- start - gene_with_same_chrom$End[gene_with_same_chrom$Pos=="Left"]
  gene_with_same_chrom$Dist[gene_with_same_chrom$Pos=="Right"] <- gene_with_same_chrom$Start[gene_with_same_chrom$Pos=="Right"] - end
  
  min_dist <- min(gene_with_same_chrom$Dist)
  clost_gene_id <- gene_with_same_chrom$Gid[gene_with_same_chrom$Dist==min_dist]
  gene_1k <- gene_with_same_chrom[gene_with_same_chrom$Dist <= 1000,]
  upstream_gene <- gene_1k$Gid[(gene_1k$Pos=="Left" & gene_1k$Strand=="+") | (gene_1k$Pos=="Right" & gene_1k$Strand=="-") | gene_1k$Pos=="Hit"]
  downstream_gene <- gene_1k$Gid[(gene_1k$Pos=="Left" & gene_1k$Strand=="-") | (gene_1k$Pos=="Right" & gene_1k$Strand=="+") | gene_1k$Pos=="Hit"]
  
  return(data.frame(
    Chrom=chrom,
    Start=start,
    End=end,
    CloestGid=paste(clost_gene_id, collapse = ","),
    ClosetGeneDist=min_dist,
    GeneUpStream1K=paste(upstream_gene, collapse = ","),
    GeneDownStream1K=paste(downstream_gene, collapse = ",")
  ))
}
ATAC_peak_DE <- ATAC_peak_DE[ATAC_peak_DE$Chrom %in% unique(ref_gene$Chrom),]
clost_gene_df <- apply(ATAC_peak_DE, 1, find_cloest_gene)
clost_gene_df <- do.call(rbind, clost_gene_df)
clost_gene_df <- right_join(ATAC_peak_DE, clost_gene_df)
write_tsv(clost_gene_df, "stat/ATAC_MACS2/ATAC.DE_peak.cloest_gene.tsv")

mean(clost_gene_df$ClosetGeneDist>1000)

p <- ggplot(clost_gene_df, aes(x=log10(ClosetGeneDist+1))) +
  geom_histogram(bins = 50, fill="black") +
  facet_grid(~Tag) +
  theme_bw() +
  scale_x_continuous(breaks = c(0, log10(1e2+1), log10(1e4+1), log10(1e6+1)), labels = c("0", "1E2", "1E4", "1E6")) +
  labs(y="#DE ATAC peaks", x="Distense to cloest gene") +
  theme(text = element_text(family="ArialMT", size=7),
        title = element_text(family="ArialMT", size=8),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid = element_blank()
  )
ggsave("stat/ATAC_MACS2/atac_DE_peak_gene_dist.pdf", p, width = 7, height = 4, units = "cm")



Marker_gene <- read_delim("stat/marker_gene/cluster/Marker.gene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_gene_expr_C3_marker <- ATAC_gene_expr[ATAC_gene_expr$Gid %in% Marker_gene$Gid[Marker_gene$Cluster=="C3"],]
ATAC_gene_expr_C3_marker <- left_join(ATAC_gene_expr_C3_marker, ATAC_gene_DE)
ATAC_gene_expr_C3_marker$Tag
p <- ggplot(ATAC_gene_expr_C3_marker, aes(x=log2(WT), y=log2(FL))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)
p

ATAC_gene_expr_C2_marker <- ATAC_gene_expr[ATAC_gene_expr$Gid %in% Marker_gene$Gid[Marker_gene$Cluster=="C2"],]
ATAC_gene_expr_C2_marker <- left_join(ATAC_gene_expr_C2_marker, ATAC_gene_DE)
ATAC_gene_expr_C2_marker$Tag
p <- ggplot(ATAC_gene_expr_C2_marker, aes(x=log2(WT), y=log2(FL))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)
p

