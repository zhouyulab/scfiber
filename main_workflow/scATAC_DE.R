library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
setwd("~/project/CottonSCE")
load("analysis/ATAC_preprocess/ATAC.gene_li.RData")

merge_gene_ATAC <- merge(gene_li[["WT"]], gene_li[["FL"]])
merge_gene_ATAC <- NormalizeData(merge_gene_ATAC, normalization.method="LogNormalize")
save(merge_gene_ATAC, file=file.path("analysis", "ATAC_preprocess", "ATAC.gene.merge.RData"))
Idents(merge_gene_ATAC) <- merge_gene_ATAC$Sample
ave_expr <- AverageExpression(merge_gene_ATAC, verbose = FALSE)$RNA
ave_expr$Gid <- gsub("-", "_", rownames(ave_expr))
ave_expr <- ave_expr[,c("Gid", "WT", "FL")]
write_tsv(ave_expr, "analysis/ATAC_DE/ATAC.gene.tsv")

ATAC_DE <- FindAllMarkers(merge_gene_ATAC, logfc.threshold=0.1, min.pct=0.01, return.thresh=0.2, only.pos = TRUE)
ave_expr$gene <- rownames(ave_expr)
ATAC_DE_df <- left_join(ATAC_DE, ave_expr)
ATAC_DE_df <- ATAC_DE_df[,c("Gid", "WT", "FL", "pct.1", "pct.2", "cluster")]
names(ATAC_DE_df)[6] <- "Tag"
ATAC_DE_df$Log2FC <- log2(ATAC_DE_df$FL / ATAC_DE_df$WT)
ATAC_DE_df$Log2pctFC <- log2(ATAC_DE_df$pct.1 / ATAC_DE_df$pct.2)
ATAC_DE_df <- ATAC_DE_df[abs(ATAC_DE_df$Log2FC)>log2(2),]
write_tsv(ATAC_DE_df, "analysis/ATAC_DE/ATAC.gene.DE.tsv")

load("analysis/ATAC_preprocess/ATAC.peak_li.RData")
merge_peak_ATAC <- merge(peak_li[["WT"]], peak_li[["FL"]])

merge_peak_ATAC <- NormalizeData(merge_peak_ATAC, scale.factor = 200000, normalization.method="LogNormalize")
rm(peak_li)
gc()

save(merge_peak_ATAC, file=file.path("analysis", "ATAC_preprocess", "ATAC.peak.merge.RData"))
Idents(merge_peak_ATAC) <- merge_peak_ATAC$Sample
ave_expr <- AverageExpression(merge_peak_ATAC, verbose = FALSE)$ATAC
ave_expr$Chrom <- sapply(strsplit(rownames(ave_expr), "-"), function(x){return(
  if(length(x)==4){
    sprintf("%s_%s", x[1], x[2])
  }else{
    x[1]
  })})
ave_expr$Start <- sapply(strsplit(rownames(ave_expr), "-"), function(x){return(as.integer(x[length(x)-1]))})
ave_expr$End <- sapply(strsplit(rownames(ave_expr), "-"), function(x){return(as.integer(x[length(x)]))})
ave_expr <- ave_expr[,c("Chrom", "Start", "End", "WT", "FL")]
write_tsv(ave_expr, "analysis/ATAC_DE/ATAC.peak.tsv")

ATAC_peak_DE <- FindAllMarkers(merge_peak_ATAC, logfc.threshold=0.1, min.pct=0.01, return.thresh=0.2, only.pos = TRUE)
ave_expr$gene <- rownames(ave_expr)
ATAC_peak_DE_df <- left_join(ATAC_peak_DE, ave_expr)
ATAC_peak_DE_df <- ATAC_peak_DE_df[,c("Chrom", "Start", "End", "WT", "FL", "pct.1", "pct.2", "cluster")]
names(ATAC_peak_DE_df)[8] <- "Tag"
ATAC_peak_DE_df$Log2FC <- log2(ATAC_peak_DE_df$FL / ATAC_peak_DE_df$WT)
ATAC_peak_DE_df$Log2pctFC <- log2(ATAC_peak_DE_df$pct.1 / ATAC_peak_DE_df$pct.2)
write_tsv(ATAC_peak_DE_df, "analysis/ATAC_DE/ATAC.peak.DE.SourceData.tsv")
ATAC_peak_weak_DF_df <- ATAC_peak_DE_df[abs(ATAC_peak_DE_df$Log2FC)>log2(1.5) & abs(ATAC_peak_DE_df$Log2FC)<=log2(2),]
write_tsv(ATAC_peak_weak_DF_df, "analysis/ATAC_DE/ATAC.peak.weakDE.tsv")
ATAC_peak_DE_df <- ATAC_peak_DE_df[abs(ATAC_peak_DE_df$Log2FC)>log2(2),]
write_tsv(ATAC_peak_DE_df, "analysis/ATAC_DE/ATAC.peak.DE.tsv")
