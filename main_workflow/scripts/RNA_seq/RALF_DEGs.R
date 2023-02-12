library(dplyr)
library(ggplot2)
library(readr)
library(DESeq2)

args = commandArgs(TRUE)
if (length(args) == 4) {
  f_featureCount <- args[1]
  padj_cutoff <- as.numeric(args[2])
  lg2fc_cutoff <- as.numeric(args[3])
  f_out <- args[4]
} else {
  q()
}


# padj_cutoff <- 0.05
# lg2fc_cutoff <- 1
# f_out <- "analysis_v3/RALF_RNA_seq/DEGs"
# f_featureCount <- "analysis_v3/RALF_RNA_seq/merge_expr/feature_counts.merge.tsv"

RALF_fc_df <- read_delim(f_featureCount, "\t", escape_double = FALSE, trim_ws = TRUE)
cnt_mat <- as.matrix(RALF_fc_df[,c("WT.rep1", "WT.rep2", "RALF.rep1", "RALF.rep2")])
rownames(cnt_mat) <- RALF_fc_df$Gid
cond_df <- data.frame(Sample=c("WT", "WT", "RALF", "RALF"))
cond_df$Sample <- factor(cond_df$Sample, levels = c("WT", "RALF"))
dds <- DESeqDataSetFromMatrix(countData = cnt_mat, colData=cond_df, design=~Sample)
dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res$Gid <- rownames(res)
res$DE <- "NC"
res$DE[(res$padj<padj_cutoff) & (res$log2FoldChange>lg2fc_cutoff)] <- "Up"
res$DE[(res$padj<padj_cutoff) & (res$log2FoldChange< (-1 * lg2fc_cutoff))] <- "Down"
DE_info <- res %>% group_by(DE) %>% summarise(Num=n())
DE_info <- DE_info[DE_info$DE!="NC",]
p <- ggplot(res, aes(x=log2FoldChange, y=-1*log10(padj+1e-10), color=DE)) +
  geom_point(size=0.3, alpha=0.3) +
  annotate("text", x=-5, y=8, label=sprintf("Up=%s\nDown=%s", DE_info$Num[DE_info$DE=="Up"], DE_info$Num[DE_info$DE=="Down"]), color="black", size=2) + 
  labs(x="log2 fold change", y="PADJ") +
  scale_y_continuous(breaks = c(0, 5, 10), labels = c("1E0", "1E-5", "1E-10")) +
  scale_color_manual(values = c("Up"="red", "Down"="blue", "NC"="grey70")) +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(f_out, "RALF.DEGs.vocano.pdf"), p, width = 6, height = 5, limitsize = FALSE, units = "cm")

p <- ggplot(DE_info, aes(x=DE, y=Num, fill=DE)) +
  geom_bar(stat="identity") +
  geom_text(mapping = aes(label=Num), vjust=0, size=1.3, color="black") +
  labs(y="#Genes") +
  scale_fill_manual(values = c("Up"="red", "Down"="blue")) +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(f_out, "RALF.DEGs.gene_num.pdf"), p, width = 4, height = 5, limitsize = FALSE, units = "cm")

write_tsv(res, file.path(f_out, "RALF.DEGs.tsv"))
