library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(e1071)
library(ComplexHeatmap)
library(circlize)

setwd("~/project/CottonSCE")

kmer_cluster_pos_df <- read_delim("data/kmer_cluster/kmer.8.graph_enrich.cluster.motif_pos.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
peak_pos_df <- read_delim("data/kmer_cluster/peak_pos.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
peak_pos_df$PeakMid <- (peak_pos_df$PeakStart + peak_pos_df$PeakEnd) / 2

bg_df <- data.frame()
for (i in 1:100) {
  peak_pos_df$RandMotifMid <- apply(peak_pos_df[,c("PeakStart", "PeakEnd")], 1, function(x){return(sample(x[1]:x[2], 1))})
  peak_pos_df$RandMotif2TSS <- peak_pos_df$RandMotifMid - peak_pos_df$TSS
  peak_pos_df$RandMotif2TSS[peak_pos_df$Strand=="-"] <- -1 * peak_pos_df$RandMotif2TSS[peak_pos_df$Strand=="-"]
  peak_pos_df$RandMotif2Peak <- peak_pos_df$RandMotifMid - peak_pos_df$PeakMid
  peak_pos_df$RandMotif2Peak[peak_pos_df$Strand=="-"] <- -1 * peak_pos_df$RandMotif2Peak[peak_pos_df$Strand=="-"]
  rand_df <- data.frame(
    RandIndx=i,
    Gid=rep(peak_pos_df$Gid, 2),
    Cis="Rand",
    Stat=c(rep("Motif2TSS", nrow(peak_pos_df)), rep("Motif2Peak", nrow(peak_pos_df))),
    Value=c(peak_pos_df$RandMotif2TSS, peak_pos_df$RandMotif2Peak)
  )
  bg_df <- rbind(bg_df, rand_df)
}

tmp_motif <- kmer_cluster_pos_df$Motif[2]
tmp_motif_df <- kmer_cluster_pos_df[kmer_cluster_pos_df$Motif==tmp_motif,]

stat_motif_pos <- function(block, bg_df, f_plot){
  block$MotifMid <- (block$MotifStart + block$MotifEnd) / 2
  block$PeakMid <- (block$PeakStart + block$PeakEnd) / 2
  block$Motif2TSS <- block$MotifMid - block$TSS
  block$Motif2TSS[block$Strand=="-"] <- -1 * block$Motif2TSS[block$Strand=="-"]
  block$Motif2Peak <- block$MotifMid - block$PeakMid
  block$Motif2Peak[block$Strand=="-"] <- -1 * block$Motif2Peak[block$Strand=="-"]
  
  signal_df <- data.frame(
    Gid=rep(block$Gid, 2),
    Cis=rep(block$CisMotif, 2),
    Stat=c(rep("Motif2TSS", nrow(block)), rep("Motif2Peak", nrow(block))),
    Value=c(block$Motif2TSS, block$Motif2Peak)
  )
  
  merge_df <- rbind(signal_df, bg_df)
  
  motif <- block$Motif[1]
  rev_motif <- block$RevMotif[1]
  motif_num <- sum(block$CisMotif=="Motif")
  rev_motif_num <- sum(block$CisMotif=="RevMotif")
  merge_df$Cis <- factor(merge_df$Cis, levels = c("Motif", "RevMotif", "Rand"), labels = c("Motif", "Reverse motif", "Random"))
  
  p <- ggplot(merge_df, aes(x=Value, color=Cis)) +
    geom_vline(xintercept = 0, size=0.3, lty=2) +
    geom_density(size=0.3) +
    facet_grid(~Stat, scales = "free") +
    coord_cartesian(xlim = c(-2000, 2000)) +
    labs(x="Position", y="Density", title = sprintf("Motif: %s (%d) Reverse motif: %s (%d)", motif, motif_num, rev_motif, rev_motif_num)) +
    scale_color_manual(values = c("Motif"="red", "Reverse motif"="blue", "Random"="grey70")) +
    theme_bw() +
    theme(
      text = element_text(family="ArialMT", color = "black", size = 6),
      axis.line = element_line(color="black", size=0.3),
      panel.grid = element_blank(),
      legend.key.size = unit(3, "mm"),
      legend.title = element_blank(),
      legend.text = element_text(family="ArialMT", color = "black", size = 6)
    )
  ggsave(file.path(f_plot, sprintf("%s.pos.pdf", motif)), p, height = 4.5, width = 12, units = "cm", colormodel = "cmyk")
  
  
  res <- signal_df %>% 
    group_by(Cis, Stat) %>% 
    summarise(Mu=mean(Value), Std=sd(Value), Medain=median(Value), Mad=mad(Value), Skew=skewness(Value), Kurt=kurtosis(Value))
  return(res)
}

pos_stat <- plyr::ddply(kmer_cluster_pos_df, "Motif", stat_motif_pos, bg_df=bg_df, f_plot=file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_pos_distribution"))
write_tsv(pos_stat, file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_pos.shape.tsv"))

RandStat <- bg_df %>% 
  group_by(RandIndx, Cis, Stat) %>% 
  summarise(Mu=mean(Value), Std=sd(Value), Medain=median(Value), Mad=mad(Value), Skew=skewness(Value), Kurt=kurtosis(Value))


motif2TSS_signal_mad <- dcast(pos_stat[pos_stat$Stat=="Motif2TSS", c("Motif", "Cis", "Mad")], Motif~Cis, value.var="Mad")
names(motif2TSS_signal_mad) <- c("Motif", "SameStrand", "DiffStrand")
motif2TSS_signal_mad$StrongMad <- apply(motif2TSS_signal_mad[,c("SameStrand", "DiffStrand")], 1, min)
motif2TSS_signal_mad$WeakMad <- apply(motif2TSS_signal_mad[,c("SameStrand", "DiffStrand")], 1, max)

motif_hit_num <- kmer_cluster_pos_df %>% group_by(Motif, CisMotif) %>% summarise(HitNum=n())
all_hit_num <- sort(unique(motif_hit_num$HitNum))

ave_rand_mad <- c()
sd_rand_mad <- c()
for(hit in all_hit_num){
  tmp_mad <- c()
  for(i in 1:100){
    tmp_pos_df <- peak_pos_df[sample(1:nrow(peak_pos_df), hit),]
    tmp_pos_df$RandMotifMid <- apply(tmp_pos_df[,c("PeakStart", "PeakEnd")], 1, function(x){return(sample(x[1]:x[2], 1))})
    tmp_pos_df$RandMotif2TSS <- tmp_pos_df$RandMotifMid - tmp_pos_df$TSS
    tmp_pos_df$RandMotif2TSS[tmp_pos_df$Strand=="-"] <- -1 * tmp_pos_df$RandMotif2TSS[tmp_pos_df$Strand=="-"]
    tmp_mad <- c(tmp_mad, mad(tmp_pos_df$RandMotif2TSS))
  }
  ave_rand_mad <- c(ave_rand_mad, mean(tmp_mad))
  sd_rand_mad <- c(sd_rand_mad, sd(tmp_mad))
}
RandMad_df <- data.frame(HitNum=all_hit_num, RandAveMad=ave_rand_mad, RandSdMad=sd_rand_mad)
motif_hit_num <- left_join(motif_hit_num, RandMad_df)
write_tsv(motif_hit_num, file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_pos.RandMad.tsv"))
motif_hit_num$MadCutoff <- motif_hit_num$RandAveMad - 2*motif_hit_num$RandSdMad
case_motif_hit_num <- dcast(motif_hit_num[, c("Motif", "CisMotif", "MadCutoff")], Motif~CisMotif, value.var="MadCutoff")
names(case_motif_hit_num) <- c("Motif", "SameStrandMadCutoff", "DiffStrandMadCutoff")

motif2TSS_signal_mad <- left_join(motif2TSS_signal_mad, case_motif_hit_num)
motif2TSS_signal_mad$SamePass <- motif2TSS_signal_mad$SameStrand<motif2TSS_signal_mad$SameStrandMadCutoff
motif2TSS_signal_mad$DiffPass <- motif2TSS_signal_mad$DiffStrand<motif2TSS_signal_mad$DiffStrandMadCutoff
motif2TSS_signal_mad$PassNum <- rowSums(motif2TSS_signal_mad[,c("SamePass", "DiffPass")]) 
motif2TSS_signal_mad$Tag <- factor(motif2TSS_signal_mad$PassNum, levels = c(0, 1, 2), labels = c("Decentralized", "1-centralized", "2-centralized"))

write_tsv(motif2TSS_signal_mad, file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_pos.Shape.SourceData.tsv"))

motif2TSS_signal_mad_info <- motif2TSS_signal_mad %>% group_by(Tag) %>% summarise(Num=n())
motif2TSS_signal_mad_info$Label <- sprintf("%s (n=%d)", motif2TSS_signal_mad_info$Tag, motif2TSS_signal_mad_info$Num)

p <- ggplot(motif2TSS_signal_mad, aes(x=StrongMad, y=WeakMad-StrongMad, color=Tag)) +
  geom_point(size=0.3) +
  scale_color_manual(values = c("Decentralized"="grey70", "1-centralized"="red", "2-centralized"="blue"), breaks = motif2TSS_signal_mad_info$Tag, labels=motif2TSS_signal_mad_info$Label) +
  labs(x="Shape MAD of sharp one bewteen motif and reverse motif (bp)", y="Difference of shape MAD between\nmotif and reverse motif (bp)") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.line = element_line(color="black", size=0.3),
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.position = c(0.8, 0.85),
    legend.key.size = unit(3, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(family="ArialMT", color = "black", size = 6)
  )
ggsave(file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_pos.Shape.1.pdf"), p, height = 6, width = 8, units = "cm", colormodel = "cmyk")

p <- ggplot(motif2TSS_signal_mad, aes(x=StrongMad, y=WeakMad, color=Tag)) +
  geom_point(size=0.3) +
  scale_color_manual(values = c("Decentralized"="grey70", "1-centralized"="red", "2-centralized"="blue"), breaks = motif2TSS_signal_mad_info$Tag, labels=motif2TSS_signal_mad_info$Label) +
  labs(x="Shape MAD of sharp one bewteen motif and reverse motif (bp)", y="Shape MAD of broad one bewteen\nmotif and reverse motif (bp)") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.line = element_line(color="black", size=0.3),
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.position = c(0.8, 0.15),
    legend.key.size = unit(3, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(family="ArialMT", color = "black", size = 6)
  )
ggsave(file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_pos.Shape.2.pdf"), p, height = 6, width = 8, units = "cm", colormodel = "cmyk")

TSS_stat <- pos_stat[pos_stat$Stat=="Motif2TSS",]
cond1 <- (TSS_stat$Motif %in% motif2TSS_signal_mad$Motif[motif2TSS_signal_mad$SamePass]) & (TSS_stat$Cis=="Motif")
cond2 <- (TSS_stat$Motif %in% motif2TSS_signal_mad$Motif[motif2TSS_signal_mad$DiffPass]) & (TSS_stat$Cis=="RevMotif")
target_motif <- TSS_stat[cond1|cond2,]
names(target_motif)[2] <- "CisMotif"
target_motif <- target_motif[,c("Motif", "CisMotif", "Medain", "Mad")]

kmer_cluster_pos_df$MotifMid <- (kmer_cluster_pos_df$MotifStart + kmer_cluster_pos_df$MotifEnd) / 2
kmer_cluster_pos_df$Motif2TSS <- kmer_cluster_pos_df$MotifMid - kmer_cluster_pos_df$TSS
kmer_cluster_pos_df$Motif2TSS[kmer_cluster_pos_df$Strand=="-"] <- -1 * kmer_cluster_pos_df$Motif2TSS[kmer_cluster_pos_df$Strand=="-"]
kmer_cluster_pos_df <- inner_join(kmer_cluster_pos_df, target_motif)
pass_kmer_cluster <- na.omit(kmer_cluster_pos_df)
pass_kmer_cluster <- pass_kmer_cluster[abs(pass_kmer_cluster$Motif2TSS-pass_kmer_cluster$Medain) < pass_kmer_cluster$Mad,]
write_tsv(pass_kmer_cluster, file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_clu.filter.SourceData.tsv"))


motif_cluster <- read_delim("analysis/kmer_graph_enrich/kmer.8.graph_enrich.umi.cluster.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
cluster_tpm <- read_delim("data/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

expr_tmp <- cluster_tpm[apply(cluster_tpm[,2:ncol(cluster_tpm)], 1, function(x){return(max(x, na.rm = TRUE))})>1,]
write_tsv(pass_kmer_cluster, file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_clu.filter.expr.SourceData.tsv"))


C3_max_gene <- expr_tmp[apply(expr_tmp[2:ncol(expr_tmp)], 1, which.max)==3,]
all_gene_num <- nrow(expr_tmp)
all_C3_max_gene_num <- nrow(C3_max_gene)

expr_pass_motif_gene <- pass_kmer_cluster[pass_kmer_cluster$Gid %in% expr_tmp$Gid,]
expr_pass_motif_gene <- left_join(expr_pass_motif_gene, motif_cluster)
expr_pass_motif_gene_C3_max <- expr_pass_motif_gene[!duplicated(expr_pass_motif_gene[,c("Gid", "Motif", "Cluster")]),] %>% 
  group_by(Motif, RevMotif, Cluster) %>% 
  summarise(
    GeneNum=n(), 
    C3MaxNum=sum(Gid %in% C3_max_gene$Gid), 
    C3MaxRatio=mean(Gid %in% C3_max_gene$Gid), 
    pval=fisher.test(matrix(c(all_gene_num, all_C3_max_gene_num, GeneNum, C3MaxNum), nrow=2))$p.value
    )

write_tsv(expr_pass_motif_gene_C3_max, file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_clu.filter.C3Max.SourceData.tsv"))

expr_pass_motif_gene_C3_max$log2FC <- log2(expr_pass_motif_gene_C3_max$C3MaxRatio/ mean(expr_tmp$Gid %in% C3_max_gene$Gid))

C3_ratio_cutoff <- c(1 / 1.5 * mean(expr_tmp$Gid %in% C3_max_gene$Gid), 1.5*mean(expr_tmp$Gid %in% C3_max_gene$Gid))
expr_pass_motif_gene_C3_max$LabelName <- sprintf("%s-%s", expr_pass_motif_gene_C3_max$Motif, expr_pass_motif_gene_C3_max$RevMotif)
expr_pass_motif_gene_C3_max$Label <- ""
expr_pass_motif_gene_C3_max$Label[grep("CACGTG",expr_pass_motif_gene_C3_max$LabelName)] <- expr_pass_motif_gene_C3_max$LabelName[grep("CACGTG",expr_pass_motif_gene_C3_max$LabelName)]
# expr_pass_motif_gene_C3_max$Label[grep("TGGGCT",expr_pass_motif_gene_C3_max$LabelName)] <- expr_pass_motif_gene_C3_max$LabelName[grep("TGGGCT",expr_pass_motif_gene_C3_max$LabelName)]
# expr_pass_motif_gene_C3_max$Label[grep("TGGGCC",expr_pass_motif_gene_C3_max$LabelName)] <- expr_pass_motif_gene_C3_max$LabelName[grep("TGGGCC",expr_pass_motif_gene_C3_max$LabelName)]
# expr_pass_motif_gene_C3_max$Label[grep("TAGGGT",expr_pass_motif_gene_C3_max$LabelName)] <- expr_pass_motif_gene_C3_max$LabelName[grep("TAGGGT",expr_pass_motif_gene_C3_max$LabelName)]

p <- ggplot(expr_pass_motif_gene_C3_max, aes(x=log2FC, y=-1*log10(pval+1e-10), color=Cluster)) +
  geom_point(data=expr_pass_motif_gene_C3_max[expr_pass_motif_gene_C3_max$Label=="" & expr_pass_motif_gene_C3_max$Cluster=="M1",], size=0.4, color="grey70") +
  geom_point(data=expr_pass_motif_gene_C3_max[expr_pass_motif_gene_C3_max$Cluster=="M2",], size=0.4) +
  geom_point(data=expr_pass_motif_gene_C3_max[expr_pass_motif_gene_C3_max$Label!="",], size=0.6) +
  geom_text(aes(label=Label), size=1.2, color="black") +
  geom_hline(yintercept = -1*log10(0.05), size=0.3, lty=2) +
  geom_vline(xintercept = log2(c(1.5, 1/1.5)), size=0.3, lty=2) +
  geom_vline(xintercept = 0, size=0.3) +
  labs(x="log2(Obs./Exp.)", y="p-value") +
  facet_grid(~Cluster) +
  scale_color_manual(values = c("M1"="blue", "M2"="red")) +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size = 6),
    title = element_text(family="ArialMT", size = 6),
    legend.position = "none",
    panel.grid = element_blank()
  )
ggsave(file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_clu.filter.C3Max.volcano.pdf"), p, height = 4.8, width = 10, units = "cm", colormodel = "cmyk")


p <- ggplot(expr_pass_motif_gene_C3_max, aes(x=GeneNum, y=C3MaxRatio*100, color=Cluster)) +
  geom_point(size=0.4) +
  geom_text(aes(label=Label), size=1.2, color="black") +
  geom_hline(yintercept = mean(expr_tmp$Gid %in% C3_max_gene$Gid)*100, size=0.3) +
  geom_hline(yintercept = C3_ratio_cutoff*100, size=0.3, lty=2) +
  annotate("text", x=10000, y=mean(expr_tmp$Gid %in% C3_max_gene$Gid)*100, label="Background", size=1.4, vjust=1.1) +
  annotate("text", x=10000, y=C3_ratio_cutoff*100, label="Cutoff (1.5 x background)", size=1.4, vjust=1.1) +
  scale_x_log10() +
  labs(x="# target gene", y="% C3 maximum expression genes") +
  facet_grid(~Cluster) +
  scale_color_manual(values = c("M1"="blue", "M2"="red")) +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", size = 6),
    title = element_text(family="ArialMT", size = 6),
    legend.position = "none",
    panel.grid = element_blank()
  )
p
ggsave(file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_clu.filter.C3Max.pdf"), p, height = 4.8, width = 10, units = "cm", colormodel = "cmyk")

M2_select_motif <- expr_pass_motif_gene_C3_max[expr_pass_motif_gene_C3_max$Cluster=="M2" & expr_pass_motif_gene_C3_max$C3MaxRatio>C3_ratio_cutoff[2],]

M2_C3_target_genes <- pass_kmer_cluster[pass_kmer_cluster$Motif %in% M2_select_motif$Motif & pass_kmer_cluster$Gid %in% C3_max_gene$Gid,]

uniq_M2_C3_target_genes <- M2_C3_target_genes[!duplicated(M2_C3_target_genes[,c("Gid", "Motif")]),]
uniq_M2_C3_target_genes <- left_join(uniq_M2_C3_target_genes, cluster_tpm)

motif_overlap_gene_matrix <- matrix(0, nrow = nrow(M2_select_motif), ncol = nrow(M2_select_motif))
colnames(motif_overlap_gene_matrix) <- M2_select_motif$Label
rownames(motif_overlap_gene_matrix) <- M2_select_motif$Label
for(i in 1:nrow(M2_select_motif)){
  motif1 <- M2_select_motif$Motif[i]
  gene_li1 <- M2_C3_target_genes$Gid[M2_C3_target_genes$Motif==motif1]
  for(j in 1:i){
    motif2 <- M2_select_motif$Motif[j]
    gene_li2 <- M2_C3_target_genes$Gid[M2_C3_target_genes$Motif==motif2]
    overlap_num <- sum(unique(gene_li1) %in% unique(gene_li2))
    motif_overlap_gene_matrix[i, j] <- overlap_num
    motif_overlap_gene_matrix[j, i] <- overlap_num
  }
}
motif_dist <- dist( t(motif_overlap_gene_matrix / sqrt(diag(motif_overlap_gene_matrix)))/ sqrt(diag(motif_overlap_gene_matrix)), method = "manhattan")
motif_clu <- hclust(motif_dist, method = "ward.D2")

plot(motif_clu)
motif_overlap_gene_matrix <- motif_overlap_gene_matrix[motif_clu$order, motif_clu$order]

ht <- Heatmap(
  log10(motif_overlap_gene_matrix+1),
  name = "#target gene (log10+1)",
  col = colorRamp2(log10(c(0, 1, 10, 100, 1000) + 1), colorRamp2(c(1, 5), c("white", "red"))(c(1, 1.1, 1.5, 3, 5))),
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.text(motif_overlap_gene_matrix[i, j], x, y, gp=gpar(fontsize = 4))
  },
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 5)
  )
  )
inche_cm=2.54
pdf(file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_clu.filter.M2.C3Max.motif.ht.pdf"), width=12/inche_cm, height=10/inche_cm)
print(ht)
dev.off()
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
# "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"

M2_C3_target_genes$MotifCluster <- M2_C3_target_genes$Motif
M2_C3_target_genes$MotifCluster[M2_C3_target_genes$Motif %in% c("AAACCCTA", "AACCCTAG", "AAAACCCT")] <- "AAAACCCTAG"
M2_C3_target_genes$MotifCluster[M2_C3_target_genes$Motif %in% c("AAGCCCAA", "AAAGCCCA", "AAAAGCCC")] <- "AAAAGCCCAA"
M2_C3_target_genes$MotifCluster[M2_C3_target_genes$Motif %in% c("AGCCCAAT", "GCCCAATA")] <- "AGCCCAATA"
M2_C3_target_genes$MotifCluster[M2_C3_target_genes$Motif %in% c("AAGGCCCA", "AGGCCCAA", "GGCCCAAA")] <- "AAGGCCCAAA"
write_tsv(M2_C3_target_genes, file.path("analysis", "ATAC_kmer_cluster2gene", "kmer_clu.filter.M2.C3Max.SourceData.tsv"))

for(clu in unique(M2_C3_target_genes$MotifCluster)){
  target_gids <- unique(M2_C3_target_genes$Gid[M2_C3_target_genes$MotifCluster==clu])
  write_csv(data.frame(target_gids), file.path("analysis", "ATAC_kmer_cluster2gene", "cluster_gene_li", sprintf("%s.txt", clu)), col_names = FALSE)
}

bg_gene_li <- C3_max_gene$Gid
write_csv(data.frame(bg_gene_li), file.path("analysis", "ATAC_kmer_cluster2gene", "cluster_gene_li", "C3.background.txt"), col_names = FALSE)
