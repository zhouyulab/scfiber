library(RColorBrewer)
library(SingleCellExperiment)
library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

load("D:/cottonSingleCell/final_cluster/mergeSCE.RData")

TP_FPKM <- read_delim("D:/cottonSingleCell/TP.FPKM.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TP_FPKM$`5DPA` <- NULL
tp_mat <- as.matrix(TP_FPKM[,2:ncol(TP_FPKM)])
rownames(tp_mat) <- TP_FPKM$Gid
expr_mat <- tp_mat[rowMax(tp_mat) > 1,]
dim(expr_mat)
diff_ratio <- apply(
  tp_mat, 1, function(x){
    s_x <- sort(x, decreasing = TRUE)
    return((s_x[1]+1e-3) / (s_x[2]+1e-3))
  }
  )
tp_spe_gene <- rownames(tp_mat)[diff_ratio>2]

TP_FPKM <- TP_FPKM[TP_FPKM$Gid %in% tp_spe_gene,]

clu_UMAP <- read_delim("D:/cottonSingleCell/clu_UMAP.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
clu_UMAP$x <- as.integer(clu_UMAP$UMAP1 * 5)/5
clu_UMAP$y <- as.integer(clu_UMAP$UMAP2 * 5)/5
clu_UMAP$cell_indx <- 1:nrow(clu_UMAP)

cnt_mat <- counts(mergeSCE)
rm(mergeSCE)
gc()

norm_mat <- 1000 * t(t(cnt_mat) / colSums(cnt_mat))

sc_expr_num <- rowSums(cnt_mat>1)
hist(log10(sc_expr_num), breaks = 50)
expr_cnt_mat <- cnt_mat[sc_expr_num>100,]

overlap_sc_mat <- expr_cnt_mat[rownames(expr_cnt_mat) %in% tp_spe_gene,]
rm(expr_cnt_mat)
dim(overlap_sc_mat)
overlap_sc_norm_mat <- norm_mat[rownames(norm_mat) %in% rownames(overlap_sc_mat),]
overlap_tp_df <- data.frame(Gid=rownames(overlap_sc_mat))
overlap_tp_df <- left_join(overlap_tp_df, TP_FPKM, by="Gid")
overlap_tp_mat <- as.matrix(overlap_tp_df[,2:ncol(overlap_tp_df)])
rownames(overlap_tp_mat) <- overlap_tp_df$Gid

overlap_tp_max_indx <- apply(overlap_tp_mat, 1, which.max)
overlap_tp_num <- table(overlap_tp_max_indx)
overlap_tp_wt <- 1/overlap_tp_num[overlap_tp_max_indx]
overlap_tp_wt <- overlap_tp_wt / sum(overlap_tp_wt)

umap_bin_info <- clu_UMAP %>% group_by(x, y) %>% summarise(Num=n())
compute_cor <- function(line){
  cell_indx <- clu_UMAP$cell_indx[clu_UMAP$x==as.numeric(line["x"]) & clu_UMAP$y==as.numeric(line["y"])]
  tmp_norm_sc_mat <- as.matrix(overlap_sc_norm_mat[,cell_indx])
  tmp_mean_expr <- rowMeans(tmp_norm_sc_mat)
  tmp_combine_mat <- cbind(tmp_mean_expr, overlap_tp_mat)
  tmp_cov <- cov.wt(tmp_combine_mat, wt=as.numeric(overlap_tp_wt))
  tmp_cor_mat <- cov2cor(tmp_cov$cov)
  tmp_cor <- tmp_cor_mat[1,2:ncol(tmp_cor_mat)]
  res <- c(as.vector(line), tmp_cor)
  names(res) <- c(names(line), names(tmp_cor))
  return(res)
}
res <- apply(umap_bin_info, 1, compute_cor)
res_df <- as.data.frame(t(res))
res_df$MaxCorTag <- apply(res_df[,colnames(overlap_tp_mat)], 1, function(line){return(names(line)[which.max(line)])})
res_df$MaxCorTag <- factor(res_df$MaxCorTag, levels = colnames(overlap_tp_mat))
max_cor <- rowMaxs(as.matrix(res_df[,colnames(overlap_tp_mat)]))
hist(max_cor, breaks=50)
res_df$MaxCorTag[max_cor<0.25] <- NA
write_tsv(res_df, "D:/cottonSingleCell/sc_tp.umap.bin.tsv")

ps <- list()
p <- ggplot(res_df) +
  geom_rect(aes(xmin=x-0.1, xmax=x+0.1, ymin=y-0.1, ymax=y+0.1, fill=MaxCorTag)) +
  scale_fill_brewer(palette="RdYlBu", direction=-1, na.value="grey70") +
  labs(x="UMAP1", y="UMAP2", fill="Time") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = c(0.2, 0.8),
    panel.grid = element_blank()
  )
ps[["All"]] <- p

max_cor <- max(res_df[,colnames(overlap_tp_mat)])

p <- ggplot(res_df) +
  geom_rect(aes(xmin=x-0.1, xmax=x+0.1, ymin=y-0.1, ymax=y+0.1, fill=m3DPA)) +
  scale_fill_continuous(low="white", high="red", limits=c(0, max_cor)) +
  labs(x="UMAP1", y="UMAP2", fill="Cor. with -3DPA") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = c(0.2, 0.8),
    panel.grid = element_blank()
  )
ps[["-3"]] <- p

p <- ggplot(res_df) +
  geom_rect(aes(xmin=x-0.1, xmax=x+0.1, ymin=y-0.1, ymax=y+0.1, fill=m1DPA)) +
  scale_fill_continuous(low="white", high="red", limits=c(0, max_cor)) +
  labs(x="UMAP1", y="UMAP2", fill="Cor. with -1DPA") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = c(0.2, 0.8),
    panel.grid = element_blank()
  )
ps[["-1"]] <- p

p <- ggplot(res_df) +
  geom_rect(aes(xmin=x-0.1, xmax=x+0.1, ymin=y-0.1, ymax=y+0.1, fill=`0DPA`)) +
  scale_fill_continuous(low="white", high="red", limits=c(0, max_cor)) +
  labs(x="UMAP1", y="UMAP2", fill="Cor. with 0DPA") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = c(0.2, 0.8),
    panel.grid = element_blank()
  )
ps[["0"]] <- p
p <- ggplot(res_df) +
  geom_rect(aes(xmin=x-0.1, xmax=x+0.1, ymin=y-0.1, ymax=y+0.1, fill=`1DPA`)) +
  scale_fill_continuous(low="white", high="red", limits=c(0, max_cor)) +
  labs(x="UMAP1", y="UMAP2", fill="Cor. with 1DPA") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = c(0.2, 0.8),
    panel.grid = element_blank()
  )
ps[["1"]] <- p

p <- ggplot(res_df) +
  geom_rect(aes(xmin=x-0.1, xmax=x+0.1, ymin=y-0.1, ymax=y+0.1, fill=`3DPA`)) +
  scale_fill_continuous(low="white", high="red", limits=c(0, max_cor)) +
  labs(x="UMAP1", y="UMAP2", fill="Cor. with 3DPA") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.background = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.position = c(0.2, 0.8),
    panel.grid = element_blank()
  )
ps[["3"]] <- p

inche_cm=2.54
pdf("D:/cottonSingleCell/sc_tp.umap.bin.pdf", width=5*3/inche_cm, height=4*2/inche_cm, family="ArialMT", colormodel = "cmyk")
grid.arrange(grobs=ps,ncol=3,nrow=2,widths=c(1,1,1), heights=c(1,1),padding = unit(0, "mm"))
dev.off()
