library(dplyr)
library(ggplot2)
library(readr)
library(DESeq2)
library(reshape2)

f_out <- "analysis_v3/RALF_RNA_seq/DEGs"
f_RALF_DEGs <- "analysis_v3/RALF_RNA_seq/DEGs/RALF.DEGs.tsv"
f_RNA_TP <- "~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv"
f_cluster_tpm <- "analysis_v3/stat/cluster_cnt/cluster.tpm.tsv"
f_Gbox <- "analysis_v3/ATAC/ATAC_DE_peak2C3/ATAC_DE_peak_target_gene.Gbox.WT_enrich.gid.txt"
f_CircadianGene_Group <- "~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CircadianGene.Group.tsv"
f_TCP <- "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.gene.tsv"
f_TCP_like <- "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_like_motif/TCP_like_motif.gene.tsv"
f_AP2 <- "analysis_v3/ATAC/ATAC_kmer_cluster2gene/AP2/topGO/AP2.gene.txt"
f_marker_gene <- "analysis_v3/stat/marker_gene/Marker.gene.tsv"
f_DiffExpr_C2_Down <- "analysis_v3/stat/diff_expr/DiffExpr.C2.Down.tsv"
f_RNA_TP_same_time_DEGs <- "~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/DE_same_time/DEGs.type.tsv"
f_PRR <- "data/genome/HAU/PRR.cotton.tsv"
f_AUXIN <- "data/Hormones/AUXIN.tsv"

RALF_DEGs_df <- read_delim(f_RALF_DEGs, "\t", escape_double = FALSE, trim_ws = TRUE)
RALF_DEGs_df <- RALF_DEGs_df[RALF_DEGs_df$DE!="NC", c("Gid", "DE")]
RALF_up_gids <- RALF_DEGs_df$Gid[RALF_DEGs_df$DE=="Up"]
RALF_down_gids <- RALF_DEGs_df$Gid[RALF_DEGs_df$DE=="Down"]

RNA_TP_CPM_df <- read_delim(f_RNA_TP, "\t", escape_double = FALSE, trim_ws = TRUE)
cluster_tpm <- read_delim(f_cluster_tpm, "\t", escape_double = FALSE, trim_ws = TRUE)
C3_max_gene <- cluster_tpm$Gid[apply(cluster_tpm[,2:ncol(cluster_tpm)], 1, function(x){ names(x)[which.max(x)]=="WT.C3" & max(x, na.rm = T)>1 })]
CircadianGene_Group <- read_delim(f_CircadianGene_Group, "\t", escape_double = FALSE, trim_ws = TRUE)

RNA_TP_same_time_DEGs <- read_delim(f_RNA_TP_same_time_DEGs, "\t", escape_double = FALSE, trim_ws = TRUE)

Marker_gene <- read_delim(f_marker_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
DiffExpr_C2_Down <- read_delim(f_DiffExpr_C2_Down, "\t", escape_double = FALSE, trim_ws = TRUE)

Gbox <- read_csv(f_Gbox)
names(Gbox) <- "Gid"

TCP_gene <- read_delim(f_TCP, "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_like_gene <- read_delim(f_TCP_like, "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_gene_df <- data.frame(Gid=unique(c(TCP_gene$Gid, TCP_like_gene$Gid)))
TCP_gene_df$IsTCP <- TCP_gene_df$Gid %in% TCP_gene$Gid
TCP_gene_df$IsTCPLike <- TCP_gene_df$Gid %in% TCP_like_gene$Gid
AP2_gene <- read_delim(f_AP2, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
names(AP2_gene) <- "Gid"

PRR <- read_delim(f_PRR, "\t", escape_double = FALSE, trim_ws = TRUE)
PRR <- PRR[,c("PRR", "Gid")]
AUXIN <- read_delim(f_AUXIN, "\t", escape_double = FALSE, trim_ws = TRUE)
names(AUXIN) <- c("AUXIN", "Gid")

do_test <- function(gids1, gids2, total_gene_num, label1, label2){
  gid1_num <- length(unique(gids1))
  gid2_num <- length(unique(gids2))
  overlap_num <- length(unique(intersect(gids1, gids2)))
  exp_overlap_num <- total_gene_num * (gid1_num/total_gene_num) * (gid2_num/total_gene_num)
  gene_num_mat <- matrix(0, nrow = 2, ncol = 2)
  gene_num_mat[1,1] <- overlap_num
  gene_num_mat[2,1] <- gid1_num - overlap_num
  gene_num_mat[1,2] <- gid2_num - overlap_num
  gene_num_mat[2,2] <- total_gene_num - gid1_num - gid1_num + overlap_num
  pval <- fisher.test(gene_num_mat)$p.value
  res <- data.frame(Set1=label1, Set2=label2, Set1Num=gid1_num, Set2Num=gid2_num, OverlapNum=overlap_num, ExpOverlapNum=exp_overlap_num, pval=pval)
  return(res)
}
RALF_gene_li <- list("RALF_Up"=RALF_up_gids, "RALF_Down"=RALF_down_gids)
test_gene_li <- list(
  "C3_Max"=C3_max_gene,
  "TCP"=TCP_gene$Gid,
  "AP2"=AP2_gene$Gid,
  "Gbox"=Gbox$Gid,
  "Circadian"=CircadianGene_Group$Gid,
  "S_Light"=CircadianGene_Group$Gid[CircadianGene_Group$GroupID=="SignificantCircadian_Light"],
  "S_Night"=CircadianGene_Group$Gid[CircadianGene_Group$GroupID=="SignificantCircadian_Night"],
  "IS_Light"=CircadianGene_Group$Gid[CircadianGene_Group$GroupID=="WeakCircadian_Light"],
  "IS_Night"=CircadianGene_Group$Gid[CircadianGene_Group$GroupID=="WeakCircadian_Night"],
  "C1_Marker"=Marker_gene$Gid[Marker_gene$Cluster=="C1"],
  "C2_Marker"=Marker_gene$Gid[Marker_gene$Cluster=="C2"],
  "C3_Marker"=Marker_gene$Gid[Marker_gene$Cluster=="C3"],
  "C4_Marker"=Marker_gene$Gid[Marker_gene$Cluster=="C4"],
  "C5_Marker"=Marker_gene$Gid[Marker_gene$Cluster=="C5"],
  "C2Down"=DiffExpr_C2_Down$Gid,
  "PRR"=PRR$Gid,
  "AUXIN"=AUXIN$Gid,
  "RNA_TP_WT_enrich"=RNA_TP_same_time_DEGs$Gid[RNA_TP_same_time_DEGs$Tag=="WtEnrich"],
  "RNA_TP_FL_enrich"=RNA_TP_same_time_DEGs$Gid[RNA_TP_same_time_DEGs$Tag=="FlEnrich"],
  "RNA_TP_both_enrich"=RNA_TP_same_time_DEGs$Gid[RNA_TP_same_time_DEGs$Tag=="BothEnrich"]
)
total_gene_num <- nrow(RNA_TP_CPM_df)
test_li <- list()
for(RALF_label in names(RALF_gene_li)){
  for(test_label in names(test_gene_li)){
    test_li[[sprintf("%s.%s", RALF_label, test_label)]] <- do_test(RALF_gene_li[[RALF_label]], test_gene_li[[test_label]], total_gene_num, RALF_label, test_label)
  }
}
test_df <- do.call(rbind, test_li)
test_df$Label <- sprintf("p=%.1E", test_df$pval)
write_tsv(test_df, file.path(f_out, "feature_test.SourceData.tsv"))
melt_df <- melt(test_df[,c("Set1", "Set2", "OverlapNum", "ExpOverlapNum")], c("Set1", "Set2"), variable.name = "Stat", value.name = "GeneNum")
melt_df$Stat <- factor(melt_df$Stat, levels = c("OverlapNum", "ExpOverlapNum"), labels = c("Overlap gene", "Exp."))
p <- ggplot(melt_df, aes(x=Stat, y=GeneNum, fill=Stat)) +
  geom_bar(stat="identity", position = position_dodge(0.7), width = 0.7) +
  geom_text(mapping = aes(label=GeneNum), vjust=0, size=1.3, color="black", position = position_dodge(0.7)) +
  geom_text(data = test_df, mapping = aes(x=1.5, y=0, label=Label, fill=NA), vjust=0, size=1.3, color="black") +
  facet_grid(Set1~Set2, scales = "free_y") +
  scale_fill_manual(values = c("Overlap gene"="black", "Exp."="grey70")) +
  labs(y="#Genes") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
p
ggsave(file.path(f_out, "RALF.DEGs.feature_test.pdf"), p, width = 20, height = 6, limitsize = FALSE, units = "cm")

RALF_up_G3_FL_enrich_df <- data.frame(
  Gid=unique(c(RALF_gene_li[["RALF_Up"]], test_gene_li[["S_Light"]], test_gene_li[["S_Night"]], test_gene_li[["IS_Light"]], test_gene_li[["IS_Night"]], test_gene_li[["RNA_TP_FL_enrich"]]))
  )
RALF_up_G3_FL_enrich_df$IsRalfUp <- RALF_up_G3_FL_enrich_df$Gid %in% RALF_gene_li[["RALF_Up"]]
RALF_up_G3_FL_enrich_df$IsLight <- RALF_up_G3_FL_enrich_df$Gid %in% c(test_gene_li[["S_Light"]], test_gene_li[["IS_Light"]])
RALF_up_G3_FL_enrich_df$IsNight <- RALF_up_G3_FL_enrich_df$Gid %in% c(test_gene_li[["S_Night"]], test_gene_li[["IS_Night"]])
RALF_up_G3_FL_enrich_df$IsFlEnrich <- RALF_up_G3_FL_enrich_df$Gid %in% test_gene_li[["RNA_TP_FL_enrich"]]

pdf(file.path(f_out, "RALF_up_G3_FL_enrich_df.pdf"), width = 48/25.4, height = 48/25.4, family="ArialMT", colormodel = "cmyk")
p <- venn.diagram(list(
  RalfUp=RALF_gene_li[["RALF_Up"]], 
  IsLight=unique(c(test_gene_li[["S_Light"]], test_gene_li[["IS_Light"]])),
  IsNight=unique(c(test_gene_li[["S_Night"]], test_gene_li[["IS_Night"]])),
  FlEnrich=test_gene_li[["RNA_TP_FL_enrich"]]
), filename = NULL,
fill = brewer.pal(4, "Dark2"),
borders = FALSE,
sub.cex = 0.3,
lwd=0.5)
grid.draw(p)
dev.off()

RALF_up_G3_FL_enrich_df <- RALF_up_G3_FL_enrich_df[RALF_up_G3_FL_enrich_df$IsRalfUp,]

blast2go <- read_delim("analysis_v3/annotation/blast2GO/blast2go.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
blast2go$X1 <-NULL
names(blast2go)[2] <- "Tid"
blast2go$Gid <- sapply(strsplit(blast2go$Tid, "[.]"), function(x){return(x[1])})
blast2go <- blast2go[!duplicated(blast2go$Gid),]

At <- read_delim("data/genome/HAU/Gh_PrimaryPeptide_vs_Arabidopsis_Annotation.txt", "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE, skip = 1)
At  <- na.omit(At)
At <- At[,c(1, 2)]
names(At) <- c("Tid", "TairID")
At$Gid <- sapply(strsplit(At$Tid, "[.]"), function(x){return(x[1])})
At <- At[!duplicated(At$Gid),]
At$Tid <- NULL

RALF_CircadianGene_FL_enrich <- left_join(RALF_up_G3_FL_enrich_df, blast2go)
RALF_CircadianGene_FL_enrich <- left_join(RALF_CircadianGene_FL_enrich, At)
write_tsv(RALF_CircadianGene_FL_enrich, file.path(f_out, "RALF_CircadianGene_FL_enrich.tsv"))

RALF_CircadianGene <- inner_join(RALF_DEGs_df, CircadianGene_Group)
RALF_CircadianGene <- left_join(RALF_CircadianGene, blast2go)
RALF_CircadianGene <- left_join(RALF_CircadianGene, At)
write_tsv(RALF_CircadianGene, file.path(f_out, "RALF.DEGs.CircadianGene.tsv"))


RALF_Gbox <- inner_join(RALF_DEGs_df, Gbox)
RALF_Gbox <- left_join(RALF_Gbox, blast2go)
RALF_Gbox <- left_join(RALF_Gbox, At)
write_tsv(RALF_Gbox, file.path(f_out, "RALF.DEGs.Gbox.tsv"))

RALF_Marker <- inner_join(RALF_DEGs_df, Marker_gene)
RALF_Marker <- left_join(RALF_Marker, blast2go)
RALF_Marker <- left_join(RALF_Marker, At)
write_tsv(RALF_Marker, file.path(f_out, "RALF.DEGs.Marker.tsv"))

circadian_gene_li <- list(
  "S_Light"=CircadianGene_Group$Gid[CircadianGene_Group$GroupID=="SignificantCircadian_Light"],
  "S_Night"=CircadianGene_Group$Gid[CircadianGene_Group$GroupID=="SignificantCircadian_Night"],
  "IS_Light"=CircadianGene_Group$Gid[CircadianGene_Group$GroupID=="WeakCircadian_Light"],
  "IS_Night"=CircadianGene_Group$Gid[CircadianGene_Group$GroupID=="WeakCircadian_Night"]
)
RNA_TP_DEGs_gene_li <- list(
  "RNA_TP_WT_enrich"=RNA_TP_same_time_DEGs$Gid[RNA_TP_same_time_DEGs$Tag=="WtEnrich"],
  "RNA_TP_FL_enrich"=RNA_TP_same_time_DEGs$Gid[RNA_TP_same_time_DEGs$Tag=="FlEnrich"],
  "RNA_TP_both_enrich"=RNA_TP_same_time_DEGs$Gid[RNA_TP_same_time_DEGs$Tag=="BothEnrich"]
)
test_li <- list()
for(circadian_label in names(circadian_gene_li)){
  for(test_label in names(RNA_TP_DEGs_gene_li)){
    test_li[[sprintf("%s.%s", circadian_label, test_label)]] <- do_test(circadian_gene_li[[circadian_label]], RNA_TP_DEGs_gene_li[[test_label]], total_gene_num, circadian_label, test_label)
  }
}
test_df <- do.call(rbind, test_li)
test_df$Label <- sprintf("p=%.1E", test_df$pval)
write_tsv(test_df, file.path(f_out, "RNA_TP_DEGs.Circadian.feature_test.SourceData.tsv"))
melt_df <- melt(test_df[,c("Set1", "Set2", "OverlapNum", "ExpOverlapNum")], c("Set1", "Set2"), variable.name = "Stat", value.name = "GeneNum")
melt_df$Stat <- factor(melt_df$Stat, levels = c("OverlapNum", "ExpOverlapNum"), labels = c("Overlap gene", "Exp."))
p <- ggplot(melt_df, aes(x=Stat, y=GeneNum, fill=Stat)) +
  geom_bar(stat="identity", position = position_dodge(0.7), width = 0.7) +
  geom_text(mapping = aes(label=GeneNum), vjust=0, size=1.3, color="black", position = position_dodge(0.7)) +
  geom_text(data = test_df, mapping = aes(x=1.5, y=0, label=Label, fill=NA), vjust=0, size=1.3, color="black") +
  facet_grid(Set1~Set2, scales = "free_y") +
  scale_fill_manual(values = c("Overlap gene"="black", "Exp."="grey70")) +
  labs(y="#Genes") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
p
ggsave(file.path(f_out, "RNA_TP_DEGs.Circadian.feature_test.pdf"), p, width = 9, height = 9, limitsize = FALSE, units = "cm")

