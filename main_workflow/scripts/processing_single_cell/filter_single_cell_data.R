# https://doi.org/10.1016/j.cell.2019.10.003

library(ggplot2)
library(DropletUtils)

args = commandArgs(TRUE)
if (length(args) == 6) {
  f_in_folder <- args[1]
  FDR_cutoff <- as.numeric(args[2])
  UMI_cutoff <- as.numeric(args[3])
  f_filter_cnt <- args[4]
  f_barcode_rank_plot <- args[5]
  f_filter_plot <- args[6]
} else {
  q()
}

# f_in_folder <- "/home/sqreb/ias/data/CottonSingleCell/analysis/cellranger/count/FL_rep2/outs/raw_feature_bc_matrix"
# FDR_cutoff <- 1e-3
# UMI_cutoff <- 5000

sc_count_expr <- read10xCounts(f_in_folder)
sc_count <- counts(sc_count_expr)

barcode_rand_res <- barcodeRanks(sc_count)
o <- order(barcode_rand_res$rank)
df <- data.frame(Rank=barcode_rand_res$rank, Total=barcode_rand_res$total)
df_filter <- data.frame(Rank=barcode_rand_res$rank[o], Total=barcode_rand_res$fitted[o])
p <- ggplot(mapping = aes(x=Rank, y=Total)) +
  theme_bw() +
  geom_point(data = df, color="black", size=0.3) +
  geom_line(data=df_filter, color="red", size=0.1) +
  geom_hline(yintercept = metadata(barcode_rand_res)$knee, color="dodgerblue", lty=2) +
  geom_hline(yintercept = metadata(barcode_rand_res)$inflection, color="forestgreen", lty=2) +
  geom_hline(yintercept = UMI_cutoff, color="black", lty=2) +
  annotate("text", x=1, y=metadata(barcode_rand_res)$knee, label="knee", vjust=1.2, hjust = 0, size=1.3) +
  annotate("text", x=1, y=metadata(barcode_rand_res)$inflection, label="inflection", vjust=1.2, hjust = 0, size=1.3) +
  annotate("text", x=1, y=UMI_cutoff, label="UMI cutoff", vjust=1.2, hjust = 0, size=1.3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x="Rank", y="Total") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank())
ggsave(f_barcode_rank_plot, p, height = 4, width = 5, units = "cm", dpi = 600)

barcodes <- sc_count_expr@colData$Barcode
empty_drops_res <- emptyDrops(sc_count)
is.cell <- (empty_drops_res$FDR <= FDR_cutoff) & (empty_drops_res$Total >= UMI_cutoff)
is.cell[is.na(is.cell)] <- FALSE
cell_sc_count <- sc_count[,is.cell]
cell_barcodes <- barcodes[is.cell]

log10_total_umi <- log10(colSums(cell_sc_count))
median_log10_total_umi <- median(log10_total_umi)
mad_log10_total_umi <- mad(log10_total_umi)
umi_cutoff <- median_log10_total_umi - 2 * mad_log10_total_umi
umi_pass <- log10_total_umi > umi_cutoff

log10_gene_num <- log10(colSums(cell_sc_count>0))
median_log10_gene_num <- median(log10_gene_num)
mad_log10_gene_num <- mad(log10_gene_num)
gene_num_cutoff <- median_log10_gene_num - 2 * mad_log10_gene_num
gene_num_pass <- log10_gene_num > gene_num_cutoff

filtered_count <- cell_sc_count[,umi_pass&gene_num_pass]
filtered_barcodes <- cell_barcodes[umi_pass&gene_num_pass]

all_log10_total_umi <- log10(colSums(sc_count))
all_log10_gene_num <- log10(colSums(sc_count>0))
df <- data.frame(logTotal=log10(empty_drops_res$Total+1), pProb=-empty_drops_res$LogProb, isCell=is.cell, clean=is.cell&(all_log10_total_umi>umi_cutoff)&(all_log10_gene_num>gene_num_cutoff))
df <- na.omit(df)
df$Type <- "Empty"
df$Type[df$isCell] <- "LQ"
df$Type[df$isCell&df$clean] <- "HQ"
df$Type <- factor(df$Type, levels = c("HQ", "LQ", "Empty"), labels = c("HQ", "LQ", "Empty"))
HQ_num <- sum(df$Type=="HQ")
LQ_num <- sum(df$Type=="LQ")
all_num <- length(is.cell)
empty_num <- all_num - sum(is.cell, na.rm = TRUE)

p <- ggplot(df, aes(x=logTotal, y=pProb, color=Type)) +
  theme_bw() +
  geom_point(size=0.3, alpha=0.4) +
  scale_color_manual(
    values = c("Empty"="black", "HQ"="red", "LQ"="blue"), 
    breaks=c("HQ", "LQ", "Empty"), 
    labels=c(
      sprintf("HQ, %d, %.2f%%", HQ_num, HQ_num/all_num*100),
      sprintf("LQ, %d, %.2f%%", LQ_num, LQ_num/all_num*100),
      sprintf("Empty, %d, %.2f%%", empty_num, empty_num/all_num*100)
      )
    ) +
  labs(x="Log total UMI count", y="-Log Probability") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.25, 0.9),
        panel.grid = element_blank())

ggsave(f_filter_plot, p, height = 6, width = 8, units = "cm", dpi = 600)


GeneNum <- colSums(filtered_count>0) 
ExprGeneNum <- colSums(filtered_count>1)
gene_num_cutoff <- quantile(GeneNum, c(0.025, 0.975))
expr_gene_num_cutoff <- quantile(ExprGeneNum, c(0.025, 0.975))

pass_id <- (GeneNum>=gene_num_cutoff[1]) & (GeneNum<=gene_num_cutoff[2]) & (ExprGeneNum>=expr_gene_num_cutoff[1]) & (ExprGeneNum<=expr_gene_num_cutoff[2])
filtered_count <- filtered_count[,pass_id]
filtered_barcodes <- filtered_barcodes[pass_id]

write10xCounts(f_filter_cnt, filtered_count, barcodes = filtered_barcodes)
