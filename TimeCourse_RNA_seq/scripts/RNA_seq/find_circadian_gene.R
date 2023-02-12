library(dplyr)
library(DESeq2)
library(readr)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(ggplot2)
library(topGO)

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
padj_cutoff <- 0.05
lg2fc_cutoff <- log2(1.5)
f_out <- "analysis/RNA_seq_TP/circadian_gene"
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


cpm_df <- featureCount_df
for(indx in 2:ncol(featureCount_df)){
  cpm_df[,indx] <- 1e6 * featureCount_df[,indx] / sum(featureCount_df[,indx])
}
write_tsv(cpm_df, file.path(f_out, "CPM.tsv"))

cpm_mat <- as.matrix(cpm_df[,2:ncol(cpm_df)])
expr_cpm <- cpm_mat[rowMeans(cpm_mat) > 1,]
cpm_cor <- cor(log10(expr_cpm+1e-5))

ht <- Heatmap(cpm_cor, name="Correlation",
              col = colorRamp2(c(0.5, 0.75, 1), c("#2c7bb6", "#ffffbf", "#d7191c")),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_dend_width = unit(5, "mm"),
              column_dend_height = unit(5, "mm"),
              show_row_dend = FALSE,
              # show_column_names = FALSE,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 6),
                title_gp = gpar(fontsize = 6),
                grid_width = unit(2, "mm"),
                grid_height = unit(2, "mm"))
)
inche_cm=2.54
pdf(file.path(f_out, "RNA_seq.Cor.pdf"), width=10/inche_cm, height=9.5/inche_cm, family="ArialMT", colormodel = "cmyk")
draw(ht)
dev.off()

cpm_merge <- cpm_df[,"Gid"]
for(s in unique(sample_df$Sample)){
  for(tp in unique(sample_df$TP)){
    labels <- sample_df$Label[sample_df$Sample==s & sample_df$TP==tp]
    tmp_df <- data.frame(Gid=cpm_df$Gid ,Expr=rowMeans(cpm_df[,labels]))
    names(tmp_df) <- c("Gid", sprintf("%s_%s", s, tp))
    cpm_merge <- left_join(cpm_merge, tmp_df)
  }
}
write_tsv(cpm_merge, file.path(f_out, "CPM.merge.tsv"))


CPM_mat <- as.matrix(cpm_merge[,2:ncol(cpm_merge)])
rownames(CPM_mat) <- cpm_merge$Gid
CPM_mat <- CPM_mat[apply(CPM_mat, 1, max)>1,]
CPM_WT_mat <- CPM_mat[,1:9]
CPM_FL_mat <- CPM_mat[,10:18]

cv_df <- data.frame(
  Gid=rownames(CPM_WT_mat), 
  WT_CV=apply(CPM_WT_mat, 1, function(x){if(max(x)==0){return(0)}else{return( sd(x)/mean(x))}}),
  FL_CV=apply(CPM_FL_mat, 1, function(x){if(max(x)==0){return(0)}else{return( sd(x)/mean(x))}})
  )
cv_df$MaxCV <- apply(cv_df[,2:3], 1, max)

CPM_WT_diff_mat <- t(apply(CPM_WT_mat, 1, function(x){return(diff(x)>0)}))
CPM_FL_diff_mat <- t(apply(CPM_FL_mat, 1, function(x){return(diff(x)>0)}))

compute_circadian_trans_num <- function(line){
  if(all(is.na(line))){
    return(c(MaxType=NA, MaxLen=0))
  }
  len <- length(line)
  continue_li <- rep(FALSE, len-1)
  for (indx in 1:(len-1)) {
    if(line[indx]!=line[indx+1]){
      continue_li[indx] <- TRUE
    }
  }
  if(!any(continue_li)){
    return(0)
  }
  continue_indx <- as.integer(which(continue_li == TRUE))
  max_len <- 0
  max_type <- NULL
  this_len <- 0
  last_indx <- NULL
  for(i in continue_indx){
    if(is.null(last_indx)){
      last_indx <- i
      this_len <- 1
    }else{
      if(i-last_indx == 1){
        this_len <- this_len + 1
        last_indx <- i
      }else{
        if(this_len > max_len){
          max_type <- xor(last_indx%%2==1, line[last_indx])
          max_len <- this_len
        }
        last_indx <- i
        this_len <- 1
      }
    }
  }
  if(this_len > max_len){
    max_type <- xor(last_indx%%2==1, line[last_indx])
    max_len <- this_len
  }
  return(data.frame(MaxType=max_type, MaxLen=max_len))
}
weak_circadian_df <- data.frame(Gid=rownames(CPM_mat))
WtWeakCircadianTransNum_li <- do.call(rbind, apply(CPM_WT_diff_mat, 1, compute_circadian_trans_num))
weak_circadian_df$WtWeakCircadianTransType <- WtWeakCircadianTransNum_li$MaxType
weak_circadian_df$WtWeakCircadianTransNum <- WtWeakCircadianTransNum_li$MaxLen
FlWeakCircadianTransNum_li <-do.call(rbind, apply(CPM_FL_diff_mat, 1, compute_circadian_trans_num))
weak_circadian_df$FlWeakCircadianTransType <- FlWeakCircadianTransNum_li$MaxType
weak_circadian_df$FlWeakCircadianTransNum <- FlWeakCircadianTransNum_li$MaxLen

weak_circadian_df <- left_join(weak_circadian_df, cv_df)
x <- c(rep(1, 4), rep(1.5, 4))
cv_cutoff <- sd(x) / mean(x)

time_num_cutoff <- 5
weak_circadian_df$CircadianState <- "NC"
weak_circadian_df$WtWeakCircadianTransType[weak_circadian_df$WtWeakCircadianTransNum<time_num_cutoff] <- -1
weak_circadian_df$FlWeakCircadianTransType[weak_circadian_df$FlWeakCircadianTransNum<time_num_cutoff] <- -1
weak_circadian_df$CircadianState[(weak_circadian_df$WtWeakCircadianTransNum>=time_num_cutoff | weak_circadian_df$FlWeakCircadianTransNum>=time_num_cutoff)] <- "WeakCircadian"
weak_circadian_df$CircadianState[weak_circadian_df$WtWeakCircadianTransType==1 & weak_circadian_df$WtWeakCircadianTransType==0] <- "NC"
weak_circadian_df$CircadianState[weak_circadian_df$WtWeakCircadianTransType==0 & weak_circadian_df$WtWeakCircadianTransType==1] <- "NC"

weak_circadian_df$CircadianType <- weak_circadian_df$WtWeakCircadianTransType
weak_circadian_df$CircadianType[weak_circadian_df$FlWeakCircadianTransNum>weak_circadian_df$WtWeakCircadianTransNum] <- weak_circadian_df$FlWeakCircadianTransType[weak_circadian_df$FlWeakCircadianTransNum>weak_circadian_df$WtWeakCircadianTransNum]
weak_circadian_df$CircadianType <- factor(weak_circadian_df$CircadianType, levels = c(0, 1, -1), labels = c("Night", "Light", "NC"))
weak_circadian_df$CircadianType[weak_circadian_df$CircadianState=="NC"] <- "NC"
table(weak_circadian_df$CircadianState)
table(weak_circadian_df$CircadianType)

all_weak_circadian_gid <- weak_circadian_df$Gid[weak_circadian_df$CircadianState!="NC"]

compute_neighbor_DEGs <- function(s, padj_cutoff=0.05, lg2fc_cutoff=1){
  tmp_sample_df <- sample_df[sample_df$Sample==s,]
  DEGs_li <- list()
  DEGs_state <- featureCount_df[,"Gid"]
  for(indx in 2:nrow(tp_df)){
    tp0 <- as.character(tp_df$TP[indx-1])
    tp1 <- as.character(tp_df$TP[indx])
    tp0_sample <- tmp_sample_df[tmp_sample_df$TP == tp0,]
    tp1_sample <- tmp_sample_df[tmp_sample_df$TP == tp1,]
    cnt_mat <- as.matrix(featureCount_df[,c(tp0_sample$Label, tp1_sample$Label)])
    rownames(cnt_mat) <- featureCount_df$Gid
    cond_df <- data.frame(TP=c(rep(as.character(tp0), nrow(tp0_sample)), rep(as.character(tp1), nrow(tp1_sample))))
    cond_df$TP <- factor(cond_df$TP, levels = c(tp0, tp1))
    compair_label <- sprintf("%s_%s", tp0, tp1)
    dds <- DESeqDataSetFromMatrix(countData = cnt_mat, colData=cond_df, design=~TP)
    dds <- DESeq(dds)
    res <- results(dds)
    res <- as.data.frame(res)
    res$Gid <- rownames(res)
    res$Compare <- compair_label
    res$DE <- "NC"
    res$DE[(res$padj<padj_cutoff) & (res$log2FoldChange>lg2fc_cutoff)] <- "Up"
    res$DE[(res$padj<padj_cutoff) & (res$log2FoldChange< (-1 * lg2fc_cutoff))] <- "Down"
    DEGs_li[[compair_label]] <- res
    res_DE_df <- res[,c("Gid", "DE")]
    names(res_DE_df) <- c("Gid", compair_label)
    DEGs_state <- left_join(DEGs_state, res_DE_df)
  }
  DEGs_state$CircadianTransNum <- apply(DEGs_state, 1, function(line){
    DE_state <- line[2:length(line)]
    up_indx <- as.integer(which(DE_state == "Up"))
    down_indx <- as.integer(which(DE_state == "Down"))
    if( (length(up_indx)==0) | (length(down_indx)==0) ){
      return(0)
    }
    if( length(unique(up_indx %% 2)) != 1 ){
      return(0)
    }
    if( length(unique(down_indx %% 2)) != 1 ){
      return(0)
    }
    if( unique(down_indx %% 2) == unique(up_indx %% 2) ){
      return(0)
    }
    all_indx <- sort(c(up_indx, down_indx))
    diff_indx <- diff(all_indx)
    return(sum(diff_indx == 1))
  })
  res <- list(
    "DEGs_li" = DEGs_li,
    "DEGs_state" = DEGs_state
  )
  return(res)
}

WT_DEGs_li <- compute_neighbor_DEGs("WT", padj_cutoff, lg2fc_cutoff)
FL_DEGs_li <- compute_neighbor_DEGs("FL", padj_cutoff, lg2fc_cutoff)

write_tsv(WT_DEGs_li$DEGs_state, file.path(f_out, "WT.DEGs.DE_state.SourceData.tsv"))
write_tsv(FL_DEGs_li$DEGs_state, file.path(f_out, "FL.DEGs.DE_state.SourceData.tsv"))

WT_DEGs_info <- do.call(rbind, WT_DEGs_li$DEGs_li)
FL_DEGs_info <- do.call(rbind, FL_DEGs_li$DEGs_li)

write_tsv(WT_DEGs_info, file.path(f_out, "WT.DEGs.SourceData.tsv"))
write_tsv(FL_DEGs_info, file.path(f_out, "FL.DEGs.SourceData.tsv"))

WT_log2fc_df <- dcast(WT_DEGs_info[,c("Gid", "Compare", "log2FoldChange")], Gid~Compare, value.var = "log2FoldChange")
WT_log2fc_df <- WT_log2fc_df[,names(WT_DEGs_li$DEGs_state)[1:(ncol(WT_DEGs_li$DEGs_state)-1)]]
WT_log2fc_df <- left_join(data.frame(Gid=WT_DEGs_li$DEGs_state$Gid), WT_log2fc_df)
FL_log2fc_df <- dcast(FL_DEGs_info[,c("Gid", "Compare", "log2FoldChange")], Gid~Compare, value.var = "log2FoldChange")
FL_log2fc_df <- FL_log2fc_df[,names(FL_DEGs_li$DEGs_state)[1:(ncol(FL_DEGs_li$DEGs_state)-1)]]
FL_log2fc_df <- left_join(data.frame(Gid=FL_DEGs_li$DEGs_state$Gid), FL_log2fc_df)

write_tsv(WT_log2fc_df, file.path(f_out, "WT.log2FC.SourceData.tsv"))
write_tsv(FL_log2fc_df, file.path(f_out, "FL.log2FC.SourceData.tsv"))

WT_DEGs_state <- WT_DEGs_li$DEGs_state[,c("Gid", "CircadianTransNum")]
names(WT_DEGs_state) <- c("Gid", "WtCircadianTransNum")
FL_DEGs_state <- FL_DEGs_li$DEGs_state[,c("Gid", "CircadianTransNum")]
names(FL_DEGs_state) <- c("Gid", "FlCircadianTransNum")
DEGs_state <- inner_join(WT_DEGs_state, FL_DEGs_state)
DEGs_state$WT_Circadian <- "NoCircadian"
DEGs_state$WT_Circadian[DEGs_state$WtCircadianTransNum==1] <- "WeakCircadian"
DEGs_state$WT_Circadian[DEGs_state$WtCircadianTransNum>1] <- "StrongCircadian"
DEGs_state$FL_Circadian <- "NoCircadian"
DEGs_state$FL_Circadian[DEGs_state$FlCircadianTransNum==1] <- "WeakCircadian"
DEGs_state$FL_Circadian[DEGs_state$FlCircadianTransNum>1] <- "StrongCircadian"
DEGs_state$WT_Circadian[!DEGs_state$Gid%in%all_weak_circadian_gid] <- "NoCircadian"
DEGs_state$FL_Circadian[!DEGs_state$Gid%in%all_weak_circadian_gid] <- "NoCircadian"

DEGs_state_info <- DEGs_state %>% group_by(WT_Circadian, FL_Circadian) %>% summarise(GeneNum=n())
DEGs_state_info <- DEGs_state_info[order(DEGs_state_info$GeneNum, decreasing = T),]
write_tsv(DEGs_state_info, file.path(f_out, "CircadianGeneNum.tsv"))

circadian_gene <- DEGs_state[DEGs_state$WT_Circadian=="StrongCircadian" | DEGs_state$FL_Circadian=="StrongCircadian", "Gid"]
write_tsv(circadian_gene, file.path(f_out, "AllCircadianGene.txt"), col_names=FALSE)

weak_circadian_df$CircadianState[weak_circadian_df$Gid%in%circadian_gene$Gid] <- "SignificantCircadian"
table(weak_circadian_df$CircadianState)
write_tsv(weak_circadian_df, file.path(f_out, "WeakCircadian.SourceData.tsv"))

weak_circadian_CPM_df <- weak_circadian_df[weak_circadian_df$CircadianState!="NC", c("Gid", "CircadianState", "CircadianType")]
weak_circadian_CPM_df <- left_join(weak_circadian_CPM_df, cpm_merge)
write_tsv(weak_circadian_CPM_df, file.path(f_out, "WeakCircadian.CPM.tsv"))
x <- data.frame(
Gid=c("Ghir_A01G003520", "Ghir_A01G008090", "Ghir_A01G017870", "Ghir_A03G005650", "Ghir_A05G009320",
"Ghir_A05G013360", "Ghir_A05G016440", "Ghir_A07G018790", "Ghir_A10G015260", "Ghir_D01G019360",
"Ghir_D02G021280", "Ghir_D03G004930", "Ghir_D03G013210", "Ghir_D05G013110", "Ghir_D07G009540",
"Ghir_D10G020290", "Ghir_D11G007710", "Ghir_D11G008250", "Ghir_D12G012000"))
x <- left_join(x, weak_circadian_df)

df2DE_mat <- function(df){
  DE_mat <- as.matrix(df[,2:ncol(df)])
  DE_num_mat <- matrix(0, nrow = nrow(DE_mat), ncol = ncol(DE_mat))
  DE_num_mat[DE_mat=="Up"] <- 1
  DE_num_mat[DE_mat=="Down"] <- -1
  rownames(DE_num_mat) <- df$Gid
  colnames(DE_num_mat) <- colnames(DE_mat)
  return(DE_num_mat)
}

cluster_tpm <- read_delim("~/mu01/project/CottonSingleCell/analysis_v3/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

circadian_gene <- data.frame(Gid=weak_circadian_CPM_df$Gid)
circadian_gene_WT_DE <- left_join(circadian_gene, WT_DEGs_li$DEGs_state)
circadian_gene_WT_DE$CircadianTransNum <- NULL
circadian_gene_WT_DE_mat <- df2DE_mat(circadian_gene_WT_DE)
colnames(circadian_gene_WT_DE_mat) <- sprintf("WT_%s", colnames(circadian_gene_WT_DE_mat))
circadian_gene_FL_DE <- left_join(circadian_gene, FL_DEGs_li$DEGs_state)
circadian_gene_FL_DE$CircadianTransNum <- NULL
circadian_gene_FL_DE_mat <- df2DE_mat(circadian_gene_FL_DE)
colnames(circadian_gene_FL_DE_mat) <- sprintf("FL_%s", colnames(circadian_gene_FL_DE_mat))
merge_DE_mat <- cbind(circadian_gene_WT_DE_mat, circadian_gene_FL_DE_mat)

circadian_cluster_tpm <- left_join(circadian_gene, cluster_tpm)
circadian_cluster_tpm_mat <- as.matrix(circadian_cluster_tpm[,2:ncol(circadian_cluster_tpm)])
rownames(circadian_cluster_tpm_mat) <- circadian_cluster_tpm$Gid
cluster_max_cpm <- apply(circadian_cluster_tpm_mat, 1, function(x){return(max(x, na.rm = T))})
cluster_norm_cpm_mat <- circadian_cluster_tpm_mat / cluster_max_cpm

circadian_cpm <- left_join(circadian_gene, cpm_merge)
circadian_cpm_mat <- as.matrix(circadian_cpm[,2:ncol(circadian_cpm)])
rownames(circadian_cpm_mat) <- circadian_cpm$Gid
max_cpm <- rowMax(circadian_cpm_mat)
norm_cpm_mat <- circadian_cpm_mat / max_cpm

circadian_WT_log2fc <- left_join(circadian_gene, WT_log2fc_df)
circadian_WT_log2fc_mat <- as.matrix(circadian_WT_log2fc[,2:ncol(circadian_WT_log2fc)])
rownames(circadian_WT_log2fc_mat) <- circadian_WT_log2fc$Gid
circadian_FL_log2fc <- left_join(circadian_gene, FL_log2fc_df)
circadian_FL_log2fc_mat <- as.matrix(circadian_FL_log2fc[,2:ncol(circadian_FL_log2fc)])
rownames(circadian_FL_log2fc_mat) <- circadian_FL_log2fc$Gid
circadian_log2fc_mat <- cbind(circadian_WT_log2fc_mat, circadian_FL_log2fc_mat)

DEGs_type <- read_delim("analysis/RNA_seq_TP/DE_same_time/DEGs.type.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


circadian_group_df <- left_join(circadian_gene, circadian_group_df)

circadian_group_df <- data.frame(Gid=weak_circadian_df$Gid, GroupID=sprintf("%s_%s", weak_circadian_df$CircadianState, weak_circadian_df$CircadianType))
circadian_group_df <- left_join(circadian_group_df, DEGs_type)
circadian_group_df <- circadian_group_df[,c("Gid", "GroupID", "Tag")]
circadian_group_df$Tag[is.na(circadian_group_df$Tag)] <- "NC"

write_tsv(circadian_group_df, file.path(f_out, "CircadianGene.Group.tsv"))

circadian_group_de_info <- circadian_group_df %>% group_by(GroupID, Tag) %>% summarise(GeneNum=n())
circadian_group_de_info$Time <- sapply(strsplit(circadian_group_de_info$GroupID, "_"), function(x){return(x[2])})
circadian_group_de_info$Significant <- sapply(strsplit(circadian_group_de_info$GroupID, "_"), function(x){return(x[1])})
circadian_group_de_info$Tag <- factor(circadian_group_de_info$Tag, levels = c("WtEnrich", "FlEnrich", "BothEnrich", "NC"), labels = c("WT enriched", "fl enriched", "Both enriched", "NC"))
p <- ggplot(circadian_group_de_info, aes(x=Time, y=GeneNum, fill=Tag)) +
  geom_bar(stat="identity", position = position_dodge(0.7), width = 0.7) +
  geom_text(mapping=aes(label=GeneNum), size=1.2, color="black", position = position_dodge(0.7), vjust=0) +
  facet_grid(Significant~., scales = "free_y") +
  scale_fill_manual(values = c("WT enriched"="#E41A1C", "fl enriched"="#377EB8", "Both enriched"="#4DAF4A", "NC"="grey70")) +
  labs(y="#Circadian genes") +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    title = element_text(family="ArialMT", color = "black", size = 6),
    axis.text = element_text(family="ArialMT", color = "black", size = 6),
    legend.title = element_blank(),
    legend.text = element_text(family="ArialMT", color = "black", size = 6),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color="black", size=0.3),
    axis.ticks = element_line(color="black", size=0.3),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill=NA, color="black", size=0.3)
  )
ggsave(file.path(f_out, "CircadianGene.DEGs.pdf"), p, height = 6, width = 8, units = "cm", colormodel = "cmyk")


DE_dist <- dist(merge_DE_mat)
norm_cpm_dist <- dist(norm_cpm_mat / 0.5)
log2fc_dist <- dist(circadian_log2fc_mat / quantile(abs(circadian_log2fc_mat), 0.90, na.rm=T), method="manhattan")
log2fc_dist[is.na(log2fc_dist)] <- max(log2fc_dist, na.rm = T)
merge_dist <- DE_dist + 0.1 * log2fc_dist + norm_cpm_dist
cl <- hclust(merge_dist, method = "ward.D2")

WT_norm_mat <- norm_cpm_mat[,1:(ncol(norm_cpm_mat)/2)]
FL_norm_mat <- norm_cpm_mat[,(1+ncol(norm_cpm_mat)/2):ncol(norm_cpm_mat)]
WT_large_mat <- WT_norm_mat > FL_norm_mat
light_WT_large_mat <- WT_large_mat[,seq(1, ncol(WT_large_mat), 2)]
night_WT_large_mat <- WT_large_mat[,seq(2, ncol(WT_large_mat), 2)]
circadian_expr_large_df <- data.frame(
  Gid=circadian_gene$Gid,
  lightWtLargeNum=rowSums(light_WT_large_mat), 
  lightFlLargeNum=rowSums(!light_WT_large_mat),
  nightWtLargeNum=rowSums(night_WT_large_mat), 
  nightFlLargeNum=rowSums(!night_WT_large_mat)
  )
circadian_expr_large_df <- left_join(circadian_expr_large_df, circadian_group_df)
circadian_expr_large_df$circadian_time <- sapply(strsplit(as.character(circadian_expr_large_df$GroupID), "_"), function(x){return(x[2])})
circadian_expr_large_df$WtLargeNum <- circadian_expr_large_df$lightWtLargeNum
circadian_expr_large_df$FlLargeNum <- circadian_expr_large_df$lightFlLargeNum
circadian_expr_large_df$WtLargeNum[circadian_expr_large_df$circadian_time=="Night"] <- circadian_expr_large_df$nightWtLargeNum[circadian_expr_large_df$circadian_time=="Night"]
circadian_expr_large_df$FlLargeNum[circadian_expr_large_df$circadian_time=="Night"] <- circadian_expr_large_df$nightFlLargeNum[circadian_expr_large_df$circadian_time=="Night"]
table(circadian_expr_large_df$WtLargeNum, circadian_expr_large_df$FlLargeNum)
circadian_expr_large_df$DE <- "NC"
circadian_expr_large_df$DE[circadian_expr_large_df$WtLargeNum<=1] <- "FL enrich"
circadian_expr_large_df$DE[circadian_expr_large_df$FlLargeNum<=1] <- "WT enrich"

PRR_cotton <- read_delim("~/mu01/project/CottonSingleCell/data/genome/HAU/PRR.cotton.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
PRR_cotton <- PRR_cotton[,c("Gid", "PRR")]
PRR_cotton <- PRR_cotton[!duplicated(PRR_cotton$Gid),]
PRR_circadian <- inner_join(PRR_cotton, circadian_expr_large_df[,c("Gid", "GroupID", "Tag", "DE")])

circadian_group_gene_num <- circadian_group_df %>% group_by(GroupID) %>% summarise(GeneNum=n())
write_tsv(circadian_group_gene_num, file.path(f_out, "CircadianGene.GroupNum.tsv"))

de_ht <- Heatmap(
  merge_DE_mat,
  name = "DE",
  cluster_columns=FALSE,
  column_split=c(rep("WT DE", ncol(circadian_gene_WT_DE_mat)), rep("FL DE", ncol(circadian_gene_FL_DE_mat))),
  show_row_names=FALSE,
  row_order = cl$order,
  row_split = circadian_group_df$GroupID,
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
cpm_ht <- Heatmap(
  norm_cpm_mat,
  name = "Norm. Expr.",
  col = colorRamp2(c(0, 0.6, 1), c("blue", "white", "red")),
  cluster_columns=FALSE,
  column_split=c(rep("WT Norm. Expr.", ncol(norm_cpm_mat)/2), rep("FL Norm. Expr.", ncol(norm_cpm_mat)/2)),
  show_row_names=FALSE,
  row_order = cl$order,
  row_split = circadian_group_df$GroupID,
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
log2fc_ht <- Heatmap(
  circadian_log2fc_mat,
  name = "log2FC",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_columns=FALSE,
  column_split=c(rep("WT log2FC", ncol(circadian_log2fc_mat)/2), rep("FL log2FC", ncol(circadian_log2fc_mat)/2)),
  show_row_names=FALSE,
  row_order = cl$order,
  row_split = circadian_group_df$GroupID,
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
cluster_ht <- Heatmap(
  cluster_norm_cpm_mat,
  name = "Norm. TPM",
  col = colorRamp2(c(0, 0.25, 1), c("magenta", "black", "yellow")),
  cluster_columns=FALSE,
  column_split=c(rep("WT scRNA", ncol(cluster_norm_cpm_mat)/2), rep("FL scRNA", ncol(cluster_norm_cpm_mat)/2)),
  show_row_names=FALSE,
  row_order = cl$order,
  row_split = circadian_group_df$GroupID,
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

select_gid <- c("Ghir_D04G020730", "Ghir_A12G028530", "Ghir_D12G028680", "Ghir_A05G021080", "Ghir_A04G012270", "Ghir_D05G028710", "Ghir_D12G012900", "Ghir_A05G028660", "Ghir_D04G016690", "Ghir_A12G012650", "Ghir_D04G017160", "Ghir_D05G013110", "Ghir_A13G007780", "Ghir_D04G004030", "Ghir_D12G000700", "Ghir_D07G002440", "Ghir_A04G012800", "Ghir_A04G012800", "Ghir_A05G005120", "Ghir_A07G008670", "Ghir_A05G038050", "Ghir_D06G016590", "Ghir_D08G018460", "Ghir_D12G017660", "Ghir_D12G017670", "Ghir_A08G017600", "Ghir_A12G017450", "Ghir_D10G011790", "Ghir_A06G015820", "Ghir_D10G025850", "Ghir_D04G005020", "Ghir_A11G003570", "Ghir_A12G028530", "Ghir_A05G039070", "Ghir_D10G012330", "Ghir_D03G016220", "Ghir_D12G028680", "Ghir_A11G028960", "Ghir_D11G003520", "Ghir_D11G029140", "Ghir_D07G008730", "Ghir_A10G023300", "Ghir_D07G015440", "Ghir_A10G015770")
interest_gene <- read_delim(f_interest_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
all_interest_gene <- interest_gene
all_interest_gene <- all_interest_gene[!duplicated(all_interest_gene$Gid),]
interest_gene <- interest_gene[interest_gene$Gid%in%select_gid,]
interest_gene <- rbind(interest_gene, data.frame(
  GeneSymbol=c("qFWPB_D02", "FE1", "FL2"),
  ID=c("_", "_", "_"),
  Gid=c("Ghir_D02G000220", "Ghir_D04G017160", "Ghir_D11G020470")
))
interest_gene <- interest_gene[!duplicated(interest_gene$Gid),]
interest_circadian_df <- left_join(circadian_gene, interest_gene)
interest_circadian_df$Label <- sprintf("%s (%s)", interest_circadian_df$Gid, interest_circadian_df$GeneSymbol)
subindx <- which(!is.na(interest_circadian_df$GeneSymbol))
label_li <- interest_circadian_df$Label[subindx]
ha <- rowAnnotation(link = row_anno_link(
  at = subindx,
  labels = label_li,
  labels_gp=gpar(fontsize = 6),
  lines_gp=gpar(lwd = 0.5),
  link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)))

x <- de_ht + log2fc_ht + cluster_ht + ha
save(x, file = file.path(f_out, "CircadianGene.HeatMap.RData"))

inche_cm=2.54
pdf(file.path(f_out, "CircadianGene.HeatMap.pdf"), width=20/inche_cm, height=40/inche_cm, family="ArialMT", colormodel = "cmyk")
de_ht + log2fc_ht + cluster_ht + ha
dev.off()


significant_circ <- circadian_expr_large_df[circadian_expr_large_df$GroupID %in% c("SignificantCircadian_Night", "SignificantCircadian_Light"),]
sig_merge_DE_mat <- merge_DE_mat[rownames(merge_DE_mat) %in% significant_circ$Gid,]
all(rownames(sig_merge_DE_mat) == significant_circ$Gid)
sig_circadian_log2fc_mat <- circadian_log2fc_mat[rownames(circadian_log2fc_mat) %in% significant_circ$Gid,]
sig_cluster_norm_cpm_mat <- cluster_norm_cpm_mat[rownames(cluster_norm_cpm_mat) %in% significant_circ$Gid,]
sig_norm_cpm_mat <- norm_cpm_mat[rownames(norm_cpm_mat) %in% significant_circ$Gid,]

sig_DE_dist <- dist(sig_merge_DE_mat)
sig_norm_cpm_dist <- dist(sig_norm_cpm_mat / 0.5)
sig_log2fc_dist <- dist(sig_circadian_log2fc_mat / quantile(abs(sig_circadian_log2fc_mat), 0.90, na.rm=T), method="manhattan")
sig_log2fc_dist[is.na(sig_log2fc_dist)] <- max(sig_log2fc_dist, na.rm = T)
sig_merge_dist <- sig_DE_dist + 0.1 * sig_log2fc_dist + sig_norm_cpm_dist
sig_merge_dist <- as.matrix(sig_merge_dist)
sig_merge_dist[significant_circ$GroupID=="SignificantCircadian_Night", significant_circ$GroupID=="SignificantCircadian_Light"] <- max(sig_merge_dist)
sig_merge_dist[significant_circ$GroupID=="SignificantCircadian_Light", significant_circ$GroupID=="SignificantCircadian_Night"] <- max(sig_merge_dist)
sig_merge_dist <- as.dist(sig_merge_dist)
sig_cl <- hclust(sig_merge_dist, method = "ward.D2")

sig_circadian_group <- cutree(sig_cl, k = 4)
significant_circ$Group <- sig_circadian_group
table(significant_circ$Group, significant_circ$GroupID)

significant_circ$Group <- factor(significant_circ$Group, levels = 1:4, labels = c("L1", "L2", "N1", "L3"))
write_tsv(significant_circ[,c("Gid", "GroupID", "Group", "DE")], file.path(f_out, "SigCircadianGene.Group.tsv"))

sig_de_ht <- Heatmap(
  sig_merge_DE_mat,
  name = "DE",
  cluster_columns=FALSE,
  column_split=c(rep("WT DE", ncol(circadian_gene_WT_DE_mat)), rep("FL DE", ncol(circadian_gene_FL_DE_mat))),
  show_row_names=FALSE,
  row_order = sig_cl$order,
  row_split = significant_circ$Group,
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
sig_log2fc_ht <- Heatmap(
  sig_circadian_log2fc_mat,
  name = "log2FC",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_columns=FALSE,
  column_split=c(rep("WT log2FC", ncol(circadian_log2fc_mat)/2), rep("FL log2FC", ncol(circadian_log2fc_mat)/2)),
  show_row_names=FALSE,
  row_order = sig_cl$order,
  row_split = significant_circ$Group,
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
sig_cluster_ht <- Heatmap(
  sig_cluster_norm_cpm_mat,
  name = "Norm. TPM",
  col = colorRamp2(c(0, 0.25, 1), c("magenta", "black", "yellow")),
  cluster_columns=FALSE,
  column_split=c(rep("WT scRNA", ncol(cluster_norm_cpm_mat)/2), rep("FL scRNA", ncol(cluster_norm_cpm_mat)/2)),
  show_row_names=FALSE,
  row_order = sig_cl$order,
  row_split = significant_circ$Group,
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)

select_gid <- c("Ghir_D04G020730", "Ghir_A12G028530", "Ghir_D12G028680", "Ghir_A05G021080", "Ghir_A04G012270", "Ghir_D05G028710", "Ghir_D12G012900", "Ghir_A05G028660", "Ghir_D04G016690", "Ghir_A12G012650", "Ghir_D04G017160", "Ghir_D05G013110", "Ghir_A13G007780", "Ghir_D04G004030", "Ghir_D12G000700", "Ghir_D07G002440", "Ghir_A04G012800", "Ghir_A04G012800", "Ghir_A05G005120", "Ghir_A07G008670", "Ghir_A05G038050", "Ghir_D06G016590", "Ghir_D08G018460", "Ghir_D12G017660", "Ghir_D12G017670", "Ghir_A08G017600", "Ghir_A12G017450", "Ghir_D10G011790", "Ghir_A06G015820", "Ghir_D10G025850", "Ghir_D04G005020", "Ghir_A11G003570", "Ghir_A12G028530", "Ghir_A05G039070", "Ghir_D10G012330", "Ghir_D03G016220", "Ghir_D12G028680", "Ghir_A11G028960", "Ghir_D11G003520", "Ghir_D11G029140", "Ghir_D07G008730", "Ghir_A10G023300", "Ghir_D07G015440", "Ghir_A10G015770")
interest_gene <- read_delim(f_interest_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
all_interest_gene <- interest_gene
all_interest_gene <- all_interest_gene[!duplicated(all_interest_gene$Gid),]
interest_gene <- interest_gene[interest_gene$Gid%in%select_gid,]
interest_gene <- rbind(interest_gene, data.frame(
  GeneSymbol=c("qFWPB_D02", "FE1", "FL2"),
  ID=c("_", "_", "_"),
  Gid=c("Ghir_D02G000220", "Ghir_D04G017160", "Ghir_D11G020470")
))
interest_gene <- rbind(interest_gene, data.frame(GeneSymbol=PRR_cotton$PRR, ID="_", Gid=PRR_cotton$Gid))
interest_gene <- interest_gene[!duplicated(interest_gene$Gid),]
interest_circadian_df <- left_join(significant_circ, interest_gene)
interest_circadian_df$Label <- sprintf("%s (%s)", interest_circadian_df$Gid, interest_circadian_df$GeneSymbol)
subindx <- which(!is.na(interest_circadian_df$GeneSymbol))
label_li <- interest_circadian_df$Label[subindx]
ha <- rowAnnotation(
  DE = significant_circ$DE,
  col=list(
    DE=c("WT enrich"="red", "FL enrich"="blue", "NC"="grey70")
  ),
  link = row_anno_link(
  at = subindx,
  labels = label_li,
  labels_gp=gpar(fontsize = 6),
  lines_gp=gpar(lwd = 0.5),
  link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)))

x <- sig_de_ht + sig_log2fc_ht + sig_cluster_ht + ha
save(x, file = file.path(f_out, "SigCircadianGene.HeatMap.RData"))

inche_cm=2.54
pdf(file.path(f_out, "CircadianGene.HeatMap.pdf"), width=20/inche_cm, height=40/inche_cm, family="ArialMT", colormodel = "cmyk")
de_ht + log2fc_ht + cluster_ht + ha
dev.off()

norm_cpm_df <- as.data.frame(norm_cpm_mat)
norm_cpm_df$Gid <- rownames(norm_cpm_df)
norm_cpm_df$Group <- group
melt_norm_cpm_df <- melt(norm_cpm_df, c("Gid", "Group"), variable.name = "SampleLabel", value.name = "NormExpr")
melt_norm_cpm_df$SampleLabel <- as.character(melt_norm_cpm_df$SampleLabel)
melt_norm_cpm_df$Sample <- sapply(strsplit(melt_norm_cpm_df$SampleLabel, "_"), function(x){return(x[1])})
melt_norm_cpm_df$TP <- sapply(strsplit(melt_norm_cpm_df$SampleLabel, "_"), function(x){return(x[2])})
melt_norm_cpm_df <- left_join(melt_norm_cpm_df, tp_df)
melt_norm_cpm_df$Sample <- factor(melt_norm_cpm_df$Sample, levels = c("WT", "FL"))
p <- ggplot(melt_norm_cpm_df) +
  geom_line(mapping = aes(x=TP_num, y=NormExpr, group=Gid, color=Sample), size=0.1, alpha=0.05) +
  facet_grid(Group~Sample) + 
  scale_x_continuous(breaks = tp_df$TP_num, labels = tp_df$TP) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%")) +
  labs(x="Time points", y="Norm. Expr.") +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.title = element_text(family="ArialMT", color = "black", size = 5)
  )
ggsave(file.path(f_out, "CircadianGene.ExprLine.pdf"), p, width = 8, height = 6, limitsize = FALSE, units = "cm")

f_term <- "data/genome/HAU/topGo.gene.map"
geneID2GO <- readMappings(file = f_term)
geneList <- factor(rep(1, length(geneID2GO)), levels=c(0, 1))
names(geneList) <- names(geneID2GO)
BP_GOdata <- new("topGOdata", ontology="BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

circadian_gene_df <- data.frame(Gid=genesInTerm(BP_GOdata, "GO:0007623")[[1]])
circadian_gene_df <- left_join(circadian_gene_df, weak_circadian_df)
circadian_gene_df$CircadianState[is.na(circadian_gene_df$CircadianState)] <- "NoExpr"
circadian_gene_df$CircadianType <- as.character(circadian_gene_df$CircadianType)
circadian_gene_df$CircadianType[is.na(circadian_gene_df$CircadianType)] <- "NoExpr"
circadian_gene_df$Label <- sprintf("%s_%s", circadian_gene_df$CircadianState, circadian_gene_df$CircadianType)
circadian_gene_df <- left_join(circadian_gene_df, all_interest_gene)
circadian_gene_df$GeneName <- circadian_gene_df$Gid
circadian_gene_df$GeneName[!is.na(circadian_gene_df$GeneSymbol)] <- sprintf("%s (%s)", circadian_gene_df$Gid[!is.na(circadian_gene_df$GeneSymbol)], circadian_gene_df$GeneSymbol[!is.na(circadian_gene_df$GeneSymbol)])

table(circadian_gene_df$CircadianState, circadian_gene_df$CircadianType)
circadian_gene_CPM_df <- left_join(data.frame(Gid=circadian_gene_df$Gid), cpm_merge)
circadian_gene_CPM_mat <- as.matrix(circadian_gene_CPM_df[,2:ncol(circadian_gene_CPM_df)])
rownames(circadian_gene_CPM_mat) <- circadian_gene_df$GeneName
circadian_gene_CPM_norm_mat <- t(scale(t(circadian_gene_CPM_mat)))

clu_dist <- dist(circadian_gene_CPM_norm_mat)
clu_dist[is.na(clu_dist)] <- max(clu_dist, na.rm = T)
cl <- hclust(clu_dist, method = "ward.D2")

ht <- Heatmap(
  circadian_gene_CPM_norm_mat,
  name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_columns=FALSE,
  cluster_rows = F,
  column_split=c(rep("WT Norm. Expr.", ncol(norm_cpm_mat)/2), rep("FL Norm. Expr.", ncol(norm_cpm_mat)/2)),
  show_row_names=T,
  row_order = cl$order,
  row_split = circadian_gene_df$Label,
  row_title_gp = gpar(foutsize = 6),
  column_title_gp = gpar(foutsize = 6),
  row_names_gp = gpar(foutsize = 5),
  column_names_gp = gpar(foutsize = 5),
  heatmap_legend_param = list(
    labels_gp = gpar(foutsize = 5),
    title_gp = gpar(foutsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
inche_cm=2.54
pdf(file.path(f_out, "Circadian.Term.HeatMap.pdf"), width=14/inche_cm, height=30/inche_cm, family="ArialMT", colormodel = "cmyk")
ht
dev.off()


f_plot_dir <- file.path(f_out, "CircadianGeneCPM_plot")
if(dir.exists(f_plot_dir)){
  unlink(f_plot_dir, recursive = TRUE)
}
dir.create(f_plot_dir, recursive = TRUE)
for(g in unique(circadian_group_df$GroupID)){
  f_group_dir <- file.path(f_plot_dir, sprintf("Group_%s", g))
  dir.create(f_group_dir, recursive = TRUE)
  tmp_circadian_group_df <- circadian_group_df[circadian_group_df$GroupID==g,]
  melt_tmp_circadian_group_df <- melt(tmp_circadian_group_df, c("Gid", "GroupID"), variable.name = "Label", value.name = "CPM")
  melt_tmp_circadian_group_df$Label <- as.character(melt_tmp_circadian_group_df$Label)
  melt_tmp_circadian_group_df$Sample <- sapply(strsplit(melt_tmp_circadian_group_df$Label, "[_]"), function(x){return(x[1])})
  melt_tmp_circadian_group_df$TP <- sapply(strsplit(melt_tmp_circadian_group_df$Label, "[_]"), function(x){return(x[2])})
  melt_tmp_circadian_group_df <- left_join(melt_tmp_circadian_group_df, tp_df)
  
  tmp_rep_cmp_group_df <- cpm_df[cpm_df$Gid %in% tmp_circadian_group_df$Gid,]
  melt_tmp_rep_cmp_group_df <- melt(tmp_rep_cmp_group_df, c("Gid"), variable.name = "Label", value.name = "CPM")
  melt_tmp_rep_cmp_group_df$Label <- as.character(melt_tmp_rep_cmp_group_df$Label)
  melt_tmp_rep_cmp_group_df$Sample <- sapply(strsplit(melt_tmp_rep_cmp_group_df$Label, "[.]"), function(x){return(x[1])})
  melt_tmp_rep_cmp_group_df$TP <- sapply(strsplit(melt_tmp_rep_cmp_group_df$Label, "[.]"), function(x){return(x[2])})
  melt_tmp_rep_cmp_group_df <- left_join(melt_tmp_rep_cmp_group_df, tp_df)
  
  tp1_df <- tp_df
  names(tp1_df) <- c("TP1_num", "TP1")
  tp2_df <- tp_df
  names(tp2_df) <- c("TP2_num", "TP2")
  
  tmp_WT_DE_group_df <- circadian_gene_WT_DE[circadian_gene_WT_DE$Gid %in% tmp_circadian_group_df$Gid,]
  melt_tmp_WT_DE_group_df <- melt(tmp_WT_DE_group_df, c("Gid"), variable.name = "DE_TPs", value.name = "DE")
  melt_tmp_WT_DE_group_df$TP1 <- sapply(strsplit(as.character(melt_tmp_WT_DE_group_df$DE_TPs), "[_]"), function(x){return(x[1])})
  melt_tmp_WT_DE_group_df$TP2 <- sapply(strsplit(as.character(melt_tmp_WT_DE_group_df$DE_TPs), "[_]"), function(x){return(x[2])})
  melt_tmp_WT_DE_group_df <- melt_tmp_WT_DE_group_df[melt_tmp_WT_DE_group_df$DE!="NC",]
  melt_tmp_WT_DE_group_df <- left_join(melt_tmp_WT_DE_group_df, tp1_df)
  melt_tmp_WT_DE_group_df <- left_join(melt_tmp_WT_DE_group_df, tp2_df)
  
  tmp_FL_DE_group_df <- circadian_gene_FL_DE[circadian_gene_FL_DE$Gid %in% tmp_circadian_group_df$Gid,]
  melt_tmp_FL_DE_group_df <- melt(tmp_FL_DE_group_df, c("Gid"), variable.name = "DE_TPs", value.name = "DE")
  melt_tmp_FL_DE_group_df$TP1 <- sapply(strsplit(as.character(melt_tmp_FL_DE_group_df$DE_TPs), "[_]"), function(x){return(x[1])})
  melt_tmp_FL_DE_group_df$TP2 <- sapply(strsplit(as.character(melt_tmp_FL_DE_group_df$DE_TPs), "[_]"), function(x){return(x[2])})
  melt_tmp_FL_DE_group_df <- melt_tmp_FL_DE_group_df[melt_tmp_FL_DE_group_df$DE!="NC",]
  melt_tmp_FL_DE_group_df <- left_join(melt_tmp_FL_DE_group_df, tp1_df)
  melt_tmp_FL_DE_group_df <- left_join(melt_tmp_FL_DE_group_df, tp2_df)
  
  for(gid in unique(tmp_circadian_group_df$Gid)){
    tmp_mu_cpm_df <- melt_tmp_circadian_group_df[melt_tmp_circadian_group_df$Gid==gid,]
    tmp_rep_cpm_df <- melt_tmp_rep_cmp_group_df[melt_tmp_rep_cmp_group_df$Gid==gid,]
    tmp_rep_cpm_info <- tmp_rep_cpm_df %>% group_by(Sample, TP_num) %>% summarise(Mu=mean(CPM), SD=sd(CPM))
    tmp_rep_cpm_info$Setoff <- 0.3
    tmp_rep_cpm_info$Setoff[tmp_rep_cpm_info$Sample=="WT"] <- -0.3
    tmp_rep_cpm_info$ErrorBarMax <- tmp_rep_cpm_info$Mu + tmp_rep_cpm_info$SD
    tmp_rep_cpm_info$ErrorBarMin <- tmp_rep_cpm_info$Mu - tmp_rep_cpm_info$SD
    tmp_rep_cpm_info$ErrorBarMin[tmp_rep_cpm_info$ErrorBarMin<0] <- 0
    
    max_scale <- max(tmp_rep_cpm_info$ErrorBarMax)
    
    tmp_WT_DE <- melt_tmp_WT_DE_group_df[melt_tmp_WT_DE_group_df$Gid==gid,]
    if(nrow(tmp_WT_DE)>0){
      tmp_WT_DE$y_setoff <- -1*0.025*max_scale
      tmp_WT_DE$y_setoff[tmp_WT_DE$DE=="Up"] <- 0.025*max_scale
    }
    tmp_FL_DE <- melt_tmp_FL_DE_group_df[melt_tmp_FL_DE_group_df$Gid==gid,]
    if(nrow(tmp_FL_DE)>0){
      tmp_FL_DE$y_setoff <- -1*0.025*max_scale
      tmp_FL_DE$y_setoff[tmp_FL_DE$DE=="Up"] <- 0.025*max_scale
    }
    
    circadian_tp_df <- tp_df
    circadian_tp_df$start_TP <- circadian_tp_df$TP_num - 2
    circadian_tp_df$end_TP <- circadian_tp_df$TP_num + 10
    circadian_tp_df$Label <- "Light"
    circadian_tp_df$Label[(circadian_tp_df$TP_num %% 24)!=0] <- "Dark"

    p <- ggplot() +
      geom_line(tmp_mu_cpm_df, mapping = aes(x=TP_num, y=CPM, color=Sample), size=0.5) +
      geom_point(tmp_rep_cpm_df, mapping = aes(x=jitter(TP_num, factor = 0.2), y=CPM, color=Sample), size=0.5) +
      geom_errorbar(tmp_rep_cpm_info, mapping = aes(x=TP_num+Setoff, ymin=ErrorBarMin, ymax=ErrorBarMax, color=Sample), size=0.2, width=1) +
      annotate("text", x=0, y=-0.03*max_scale, color="black", label="Light", size=2) +
      geom_rect(data = circadian_tp_df, mapping=aes(xmin=start_TP, xmax=end_TP, ymax=-1*0.07*max_scale, ymin=-1*0.1*max_scale, fill=Label)) +
      annotate("text", x=0, y=-0.16*max_scale, color="black", label="WT DE", size=2) +
      annotate("text", x=0, y=-0.32*max_scale, color="black", label="FL DE", size=2) 
    if(nrow(tmp_WT_DE)>0){
      p <- p + geom_segment(
        tmp_WT_DE, 
        mapping = aes(x=TP1_num+1, xend=TP2_num-1, y=-1*0.25*max_scale - y_setoff, yend=-1*0.25*max_scale + y_setoff), 
        arrow = arrow(length = unit(0.1, "cm"), type="closed"),
        size=0.2, color="black") 
    }
    if(nrow(tmp_FL_DE)>0){
      p <- p + geom_segment(
        tmp_FL_DE, 
        mapping = aes(x=TP1_num+1, xend=TP2_num-1, y=-1*0.42*max_scale - y_setoff, yend=-1*0.42*max_scale + y_setoff), 
        arrow = arrow(length = unit(0.1, "cm"), type="closed"),
        size=0.2, color="black")
    }
    p <- p + scale_x_continuous(breaks = tp_df$TP_num, labels = tp_df$TP) +
      scale_color_manual(values = c(WT="red", FL="blue")) +
      scale_fill_manual(values = c(Light="yellow", Dark="grey40")) +
      labs(x="Time points", y="CPM", title=gid) +
      theme_bw() +
      theme(
        text = element_text(family="ArialMT", color = "black", size = 5),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.1, 0.8),
        legend.key.size = unit(2, "mm"),
        legend.spacing = unit(1, "mm"),
        legend.title = element_blank()
      )
    ggsave(file.path(f_group_dir, sprintf("%s.pdf", gid)), p, width = 8, height = 6, limitsize = FALSE, units = "cm")
  }
}


