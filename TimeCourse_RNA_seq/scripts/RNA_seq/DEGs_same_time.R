library(dplyr)
library(DESeq2)
library(readr)
library(reshape2)
library(ggplot2)

args = commandArgs(TRUE)
if (length(args) == 5) {
  f_featureCount <- args[1]
  f_interest_gene <- args[2]
  padj_cutoff <- as.numeric(args[3])
  lg2fc_cutoff <- as.numeric(args[4])
  f_out <- args[5]
} else {
  q()
}

f_featureCount <- "analysis/RNA_seq_TP/merge_expr/feature_counts.merge.tsv"
padj_cutoff <- 1e-3
lg2fc_cutoff <- 1
f_out <- "analysis/RNA_seq_TP/DE_same_time"
f_interest_gene <- "data/InterestGeneID.tsv"

featureCount_df <- read_delim(f_featureCount, "\t", escape_double = FALSE, trim_ws = TRUE)
sample_df <- data.frame(Label=names(featureCount_df)[2:ncol(featureCount_df)])
sample_df$Label <- as.character(sample_df$Label)
sample_df$Sample <- sapply(strsplit(sample_df$Label, "[.]"), function(x){return(x[1])})
sample_df$TP <- sapply(strsplit(sample_df$Label, "[.]"), function(x){return(x[2])})
sample_df$Rep <- sapply(strsplit(sample_df$Label, "[.]"), function(x){return(x[3])})
tp_df <- data.frame(
  TP_num=c(-48, -36, -24, -12, 0, 12, 24, 36, 48),
  TP=c("n48h", "n36h", "n24h", "n12h", "0h", "12h", "24h", "36h", "48h")
)
sample_df <- left_join(sample_df, tp_df)
sample_df$TP <- factor(sample_df$TP, levels = tp_df$TP)
tp_df$TP <- factor(tp_df$TP, levels = tp_df$TP)


compute_neighbor_DEGs <- function(tp, padj_cutoff=1e-3, lg2fc_cutoff=1){
  tmp_sample_df <- sample_df[sample_df$TP==tp,]
  DEGs_li <- list()
  DEGs_state <- featureCount_df[,"Gid"]
  wt_sample <- tmp_sample_df[tmp_sample_df$Sample == "WT",]
  fl_sample <- tmp_sample_df[tmp_sample_df$Sample == "FL",]
  cnt_mat <- as.matrix(featureCount_df[,c(wt_sample$Label, fl_sample$Label)])
  rownames(cnt_mat) <- featureCount_df$Gid
  cond_df <- data.frame(Sample=c(rep("WT", nrow(wt_sample)), rep("FL", nrow(fl_sample))))
  cond_df$Sample <- factor(cond_df$Sample, levels = c("WT", "FL"))
  dds <- DESeqDataSetFromMatrix(countData = cnt_mat, colData=cond_df, design=~Sample)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  res$Gid <- rownames(res)
  res$TP <- tp
  res$DE <- "NC"
  res$DE[(res$padj<padj_cutoff) & (res$log2FoldChange>lg2fc_cutoff)] <- "FL"
  res$DE[(res$padj<padj_cutoff) & (res$log2FoldChange< (-1 * lg2fc_cutoff))] <- "WT"
  return(res)
}

res_li <- list()
for(tp in unique(sample_df$TP)){
  res_li[[tp]] <- compute_neighbor_DEGs(tp, padj_cutoff, lg2fc_cutoff)
}
res_df <- do.call(rbind, res_li)
write_tsv(res_df, file.path(f_out, "RNA_seq.SameTime.DE.SourceData.tsv"))
DEGs <- res_df[res_df$DE!="NC",]
write_tsv(res_df, file.path(f_out, "RNA_seq.SameTime.DEGs.tsv"))

tps <- sample_df$TP[!duplicated(sample_df$TP)]
DE_info <- DEGs %>% group_by(TP, DE) %>% summarise(GeneNum=n())
DE_info$DE <- factor(DE_info$DE, levels = c("WT", "FL"), labels = c("WT enriched", "FL enriched"))
DE_info$TP <- factor(DE_info$TP, levels=tps)
p <- ggplot(DE_info, aes(x=TP, y=GeneNum, fill=DE)) +
  geom_bar(stat="identity", position = position_dodge(0.7), width = 0.7) +
  geom_text(mapping = aes(label=GeneNum), size=1.3, color="black", position = position_dodge(0.7), vjust=0) +
  scale_fill_manual(values = c("WT enriched"="red", "FL enriched"="blue")) +
  labs(x="Time", y="#Genes") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.title = element_text(family="ArialMT", color = "black", size = 5)
  )
p
ggsave(file.path(f_out, "DEGs.num.pdf"), p, width = 7, height = 4, limitsize = FALSE, units = "cm")

all_DEGs_gids <- unique(DEGs$Gid)
DEGs_mat <- matrix("NC", nrow = length(all_DEGs_gids), ncol = length(tps))
for(indx in 1:nrow(DEGs)){
  DEGs_mat[which(all_DEGs_gids==DEGs$Gid[indx]), which(tps==DEGs$TP[indx])] <- DEGs$DE[indx]
}
rownames(DEGs_mat) <- all_DEGs_gids
colnames(DEGs_mat) <- tps
DEGs_gene_info <- data.frame(
  Gid=all_DEGs_gids, 
  WT_enrich_num=apply(DEGs_mat, 1, function(x){return(sum(x=="WT"))}),
  FL_enrich_num=apply(DEGs_mat, 1, function(x){return(sum(x=="FL"))})
  )
DEGs_gene_info$Tag <- "BothEnrich"
DEGs_gene_info$Tag[DEGs_gene_info$FL_enrich_num==0] <- "WtEnrich"
DEGs_gene_info$Tag[DEGs_gene_info$WT_enrich_num==0] <- "FlEnrich"
DEGs_gene_info_num <- DEGs_gene_info %>% group_by(Tag) %>% summarise(GeneNum=n())
DEGs_gene_info_num$Tag <- factor(DEGs_gene_info_num$Tag, levels = c("WtEnrich", "FlEnrich", "BothEnrich"))
p <- ggplot(DEGs_gene_info_num, aes(x=Tag, y=GeneNum, fill=Tag)) +
  geom_bar(stat="identity") +
  geom_text(mapping = aes(label=GeneNum), size=1.3, color="black", vjust=0) +
  scale_fill_manual(values = c("WtEnrich"="red", "FlEnrich"="blue", "BothEnrich"="grey70")) +
  labs(y="#Genes") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  )
ggsave(file.path(f_out, "DEGs.type.num.pdf"), p, width = 4, height = 4, limitsize = FALSE, units = "cm")
write_tsv(DEGs_gene_info, file.path(f_out, "DEGs.type.tsv"))
