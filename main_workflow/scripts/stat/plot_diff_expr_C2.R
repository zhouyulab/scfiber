library(ggplot2)
library(dplyr)
library(readr)
library(ComplexHeatmap)
library(circlize)
args = commandArgs(TRUE)
if (length(args) == 5) {
  f_diff_expr <- args[1]
  f_tf <- args[2]
  f_tpm <- args[3]
  f_tp_fpkm <- args[4]
  f_out <- args[5]
} else {
  q()
}

# f_tpm <- "/home/sqreb/ias/data/CottonSingleCell/analysis/cluster_cnt/data/cluster.tpm.tsv"
# f_tf <- "/home/sqreb/ias/data/CottonSingleCell/analysis/annotation/TF/TF.gene_list.tsv"
# f_tp_fpkm <- "/home/sqreb/ias/data/CottonSingleCell/analysis/RNA_seq/merge_expr/TP.FPKM.tsv"
# f_diff_expr <- "/home/sqreb/ias/data/CottonSingleCell/analysis/stat/diff_expr/DiffExpr.tsv"

diff_expr_df <- read_delim(f_diff_expr, "\t", escape_double = FALSE, trim_ws = TRUE)
tf_df <- read_delim(f_tf, "\t", escape_double = FALSE, trim_ws = TRUE)
tpm <- read_delim(f_tpm, "\t", escape_double = FALSE, trim_ws = TRUE)
tp_fpkm <- read_delim(f_tp_fpkm, "\t", escape_double = FALSE, trim_ws = TRUE)

diff_expr_num <- diff_expr_df %>% group_by(Cluster, Reg) %>% summarise(GeneNum=n())
df <- data.frame(Cluster=c("C1", "C1", "C2", "C2", "C4", "C4", "C5", "C5"), Reg=rep(c("Up", "Down"), 4))
df <- left_join(df, diff_expr_num)
df$GeneNum[is.na(df$GeneNum)] <- 0
df$Reg <- factor(df$Reg, levels = c("Up", "Down"))
p <- ggplot(df, aes(x=Cluster, y=GeneNum, fill=Reg)) +
  theme_bw() +
  geom_bar(stat="identity", position = position_dodge(0.7), width = 0.7) +
  geom_text(aes(label=GeneNum), position = position_dodge(0.7), vjust=-0.2, size=1.2) +
  scale_fill_manual(values = c("Up"="#E41A1C", "Down"="#377EB8")) +
  labs(y="#Gene", fill="FL/WT") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())

ggsave(file.path(f_out, "DiffExpr.stat.pdf"), p, height = 4, width = 4.8, units = "cm", colormodel = "cmyk")

tf_summary <- tf_df %>% group_by(Gid) %>% summarise(TF=paste(sort(unique(Family)), collapse = ", "))

c2_down_df <- diff_expr_df[diff_expr_df$Cluster=="C2" & diff_expr_df$Reg=="Down",]
c2_down_df <- c2_down_df[order(c2_down_df$logFC),]
c2_down_df$FC_rank <- 1:nrow(c2_down_df)
c2_down_df <- left_join(c2_down_df, tf_summary, by="Gid")
c2_down_df$Label <- sprintf("%s (%s)", c2_down_df$Gid, c2_down_df$TF)
c2_down_df$IsTF <- ! is.na(c2_down_df$TF)
c2_down_df$Label[!c2_down_df$IsTF] <- ""
c2_down_df$IsTF <- factor(c2_down_df$IsTF, levels=c(T, F), labels=c("T", "F"))

# c2_down_df <- left_join(c2_down_df, tpm, by="Gid")

c2_down_tp_fpkm <- data.frame(Gid=c2_down_df$Gid)
c2_down_tp_fpkm <- left_join(c2_down_tp_fpkm, tp_fpkm, by="Gid")

names(c2_down_df)
tpm_mat <- as.matrix(c2_down_df[,(ncol(c2_down_df)-13):(ncol(c2_down_df)-4)])
rownames(tpm_mat) <- c2_down_df$Gid
colnames(tpm_mat)

tp_fpkm_mat <- as.matrix(c2_down_tp_fpkm[,2:ncol(c2_down_tp_fpkm)])
rownames(tp_fpkm_mat) <- c2_down_tp_fpkm$Gid
tp_fpkm_ratio_mat <- tp_fpkm_mat / apply(tp_fpkm_mat, 1, max)

ratio_mat <- tpm_mat / tpm_mat[,"WT.C2"]
is_wt_c3_max <- apply(ratio_mat, 1, which.max) == 3
is_wt_c2_max <- apply(ratio_mat, 1, which.max) == 2
df_ha <- data.frame(Type=rep("Other", nrow(ratio_mat)))
df_ha$Type <- as.character(df_ha$Type)
df_ha$Type[is_wt_c3_max] <- "WT.C3 enrich"
df_ha$Type[is_wt_c2_max] <- "WT.C2 enrich"
rownames(df_ha) <- rownames(ratio_mat)
ha <- rowAnnotation(
  df = df_ha, 
  col = list(Type = c("WT.C2 enrich" =  "#408e35", "WT.C3 enrich" = "#ca7e1a", "Other"="white")),
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm"))
  )

ht <- Heatmap(log2(ratio_mat+1e-3), name="log2 (Expr / WT.C2)",
              col =  colorRamp2(c(-3, 0, 3), c("blue", "grey95", "red")),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_dend_width = unit(5, "mm"),
              column_dend_height = unit(5, "mm"),
              show_row_dend = FALSE,
              show_row_names = FALSE,
              cluster_columns = FALSE,
              row_title_gp = gpar(fontsize = 6),
              split = df_ha$Type,
              # cluster_rows = FALSE,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 6),
                title_gp = gpar(fontsize = 6),
                grid_width = unit(2, "mm"),
                grid_height = unit(2, "mm"))
)

ht_tp <- Heatmap(tp_fpkm_ratio_mat, name="Norm, Expr.",
                 col =  colorRamp2(c(0, 0.25, 1), c("magenta", "black", "yellow")),
                 row_names_gp = gpar(fontsize = 6),
                 column_names_gp = gpar(fontsize = 6),
                 row_title_gp = gpar(fontsize = 6),
                 row_dend_width = unit(5, "mm"),
                 column_dend_height = unit(5, "mm"),
                 show_row_dend = FALSE,
                 show_row_names = FALSE,
                 cluster_columns = FALSE,
                 split = df_ha$Type,
                 # cluster_rows = FALSE,
                 heatmap_legend_param = list(
                   labels_gp = gpar(fontsize = 6),
                   title_gp = gpar(fontsize = 6),
                   grid_width = unit(2, "mm"),
                   grid_height = unit(2, "mm"))
)


ha_tf <- rowAnnotation(link = row_anno_link(
  at = which(c2_down_df$IsTF=="T"), 
  labels = c2_down_df$Label[c2_down_df$IsTF=="T"], 
  labels_gp=gpar(fontsize = 6),
  lines_gp=gpar(lwd = 0.5),
  link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(c2_down_df$Label[c2_down_df$IsTF=="T"], gp = gpar(fontsize = 6)))

inche_cm=2.54
pdf(file.path(f_out, "DiffExpr.C2.Down.TPM2C2.ht.pdf"), width=18/inche_cm, height=4/inche_cm, family="ArialMT")
draw(ht + ht_tp + ha + ha_tf)
dev.off()

c2_down_df <- left_join(c2_down_df, tp_fpkm, by="Gid")

c2_down_df$Type <- df_ha$Type
p <- ggplot(c2_down_df, aes(x=FC_rank, y=-logFC)) +
  theme_bw() +
  geom_line(size=0.1) +
  geom_point(data=c2_down_df[c2_down_df$IsTF=="F",], color="grey30", size=0.1, alpha=0.6) +
  geom_point(data=c2_down_df[c2_down_df$IsTF=="T",], aes(color=Type), size=0.5) +
  geom_text(aes(label=Label), hjust=-0.05, vjust=-0.2, size=1.2) +
  scale_y_continuous(breaks = c(2, 4, 6, 8), labels = c(4, 16, 64, 256)) +
  scale_color_manual(values = c("WT.C2 enrich" =  "#408e35", "WT.C3 enrich" = "#ca7e1a", "Other"="black")) +
  labs(x="Fold change rank", y="Fold change") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        panel.grid = element_blank())
ggsave(file.path(f_out, "DiffExpr.C2.Down.TF.pdf"), p, height = 4.4, width = 7, units = "cm", colormodel = "cmyk")

res <- c2_down_df
res$IsTF <- NULL
res$Label <- NULL
write_tsv(res, file.path(f_out, "DiffExpr.C2.Down.tsv"))
