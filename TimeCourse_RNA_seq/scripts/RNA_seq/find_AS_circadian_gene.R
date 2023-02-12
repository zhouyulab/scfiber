library(argparse)
library(dplyr)
library(readr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", required=TRUE, type="character", dest = "input", metavar="input_dir")
parser$add_argument("-t", "--tp", nargs="*", required=TRUE, type="character", dest = "tp", metavar="tp")
parser$add_argument("-s", "--sample", nargs="*", required=TRUE, type="character", dest = "sample", metavar="sample")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
parser$add_argument("--FDR-cutoff", type="double", dest = "FDR_cutoff", metavar="FDR_cutoff", default=0.05)
parser$add_argument("--delta-phi-cutoff", type="double", dest = "delta_phi_cutoff", metavar="delta_phi_cutoff", default=0.2)
parser$add_argument("--read-num-cutoff", type="integer", dest = "read_num_cutoff", metavar="read_num_cutoff", default=20)
parser$add_argument("--interest-gene", required=TRUE, type="character", dest = "interest_gene", metavar="InterestGeneID.tsv")
parser$add_argument("--cpm", required=TRUE, type="character", dest = "cpm", metavar="cpm.tsv")

args <- commandArgs(TRUE)
args <- parser$parse_args(args)

# args <- c(
#   "-i", "analysis/RNA_seq_TP/rMATS_turbo",
#   "-t", "n48h", "n36h", "n24h", "n12h", "0h", "12h", "24h", "36h", "48h",
#   "-s", "WT", "FL",
#   "-o", "analysis/RNA_seq_TP/AS_circadian_gene",
#   "--interest-gene", "data/InterestGeneID.tsv",
#   "--cpm", "analysis/RNA_seq_TP/circadian_gene/CPM.tsv"
# )
# args <- parser$parse_args(args)

AS_type_li <- c("A3SS", "A5SS", "RI", "SE")

load_data <- function(fname, s, last_time, this_time, AS_type){
  df <- read_delim(fname, "\t", escape_double = FALSE, trim_ws = TRUE)
  df$chr <- sapply(strsplit(df$chr, "chr"), function(x){return(x[2])})
  if(AS_type %in% c("SE")){
    df$Key <- sprintf("%s_%s_%s_%s_%s_%s_%s_%s", df$chr, df$exonStart_0base, df$exonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, df$strand)
  }else{
    if(AS_type %in% c("A3SS", "A5SS")){
      df$Key <- sprintf("%s_%s_%s_%s_%s_%s_%s_%s", df$chr, df$longExonStart_0base, df$longExonEnd, df$shortES, df$shortEE, df$flankingES, df$flankingEE, df$strand)
    }else{
      if(AS_type %in% c("RI")){
        df$Key <- sprintf("%s_%s_%s_%s_%s_%s_%s_%s", df$chr, df$riExonStart_0base, df$riExonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, df$strand)
      }
    }
  }
  df$Sample <- s
  df$LastTime <- last_time
  df$ThisTime <- this_time
  df$AsType <- AS_type
  df <- df[,c(
    "Sample", "LastTime", "ThisTime", "AsType", "GeneID", "Key", 
    "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", 
    "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference"
    )]
  return(df)
}

all_rmats_res_df <- data.frame()
for(s in args$sample){
  for(i in 1:(length(args$tp)-1)){
    last_time <- args$tp[i]
    this_time <- args$tp[i+1]
    for(AS_type in AS_type_li){
      fname <- file.path(args$input, s, sprintf("%s_vs_%s", last_time, this_time), sprintf("%s.MATS.JC.txt", AS_type))
      all_rmats_res_df <- rbind(all_rmats_res_df, load_data(fname, s, last_time, this_time, AS_type))
    }
  }
}

all_rmats_res_df$AllSample1IJC <- sapply(strsplit(all_rmats_res_df$IJC_SAMPLE_1, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
all_rmats_res_df$AllSample1SJC <- sapply(strsplit(all_rmats_res_df$SJC_SAMPLE_1, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
all_rmats_res_df$AllSample2IJC <- sapply(strsplit(all_rmats_res_df$IJC_SAMPLE_2, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
all_rmats_res_df$AllSample2SJC <- sapply(strsplit(all_rmats_res_df$SJC_SAMPLE_2, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
all_rmats_res_df$AllSample1 <- all_rmats_res_df$AllSample1IJC + all_rmats_res_df$AllSample1SJC
all_rmats_res_df$AllSample2 <- all_rmats_res_df$AllSample2IJC + all_rmats_res_df$AllSample2SJC

all_rmats_res_df$AveIncLevel1 <- sapply(strsplit(all_rmats_res_df$IncLevel1, ","), function(x){return(mean(as.numeric(x), na.rm = T))})
all_rmats_res_df$AveIncLevel2 <- sapply(strsplit(all_rmats_res_df$IncLevel2, ","), function(x){return(mean(as.numeric(x), na.rm = T))})

all_rmats_res_df$DE <- "NC"
all_rmats_res_df$DE[
  all_rmats_res_df$IncLevelDifference>args$delta_phi_cutoff & all_rmats_res_df$FDR < args$FDR_cutoff &
    all_rmats_res_df$AllSample1>args$read_num_cutoff & all_rmats_res_df$AllSample2>args$read_num_cutoff
  ] <- "Up"
all_rmats_res_df$DE[
  (-1*all_rmats_res_df$IncLevelDifference)>args$delta_phi_cutoff & all_rmats_res_df$FDR < args$FDR_cutoff &
    all_rmats_res_df$AllSample1>args$read_num_cutoff & all_rmats_res_df$AllSample2>args$read_num_cutoff
  ] <- "Down"
write_tsv(all_rmats_res_df, file.path(args$output, "AS.rMATS.SourceData.tsv"))
DE_rmats_res_df <- all_rmats_res_df[all_rmats_res_df$DE != "NC",]
DE_AS_num_stat <- DE_rmats_res_df %>% group_by(Sample, LastTime, ThisTime, AsType, DE) %>% summarise(DE_AS_num=n())
write_tsv(DE_AS_num_stat, file.path(args$output, "AS.DE.AsNum.tsv"))
DE_AS_num_stat$Times <- sprintf("%s-%s", DE_AS_num_stat$LastTime, DE_AS_num_stat$ThisTime)
time_pair <- c()
for(i in 1:(length(args$tp)-1)){
  last_time <- args$tp[i]
  this_time <- args$tp[i+1]
  time_pair <- c(time_pair, sprintf("%s-%s", last_time, this_time))
}
DE_AS_num_stat$Times <- factor(DE_AS_num_stat$Times, levels = time_pair)

p <- ggplot(DE_AS_num_stat, aes(x=Times, y=DE_AS_num, fill=DE)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_grid(Sample~AsType) +
  scale_fill_manual(values = c("Down"="blue", "Up"="red")) +
  labs(y="#DE AS") +
  theme_bw() +
  theme(
    axis.text = element_text(family="ArialMT", color = "black", size = 5),
    axis.title = element_text(family="ArialMT", color = "black", size = 5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family="ArialMT", color = "black", size = 5),
    panel.grid = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.text = element_text(family="ArialMT", color = "black", size = 5),
    legend.title = element_blank()
  )
ggsave(file.path(args$output, "DE_AD.num.pdf"), p, width = 14, height = 6, limitsize = FALSE, units = "cm")

uniq_DE_rmats_res_df <- DE_rmats_res_df[,c("AsType", "GeneID", "Key")]
uniq_DE_rmats_res_df <- uniq_DE_rmats_res_df[!duplicated(uniq_DE_rmats_res_df),]

all_rmats_res_df$TimePair <- sprintf("%s-%s", all_rmats_res_df$LastTime, all_rmats_res_df$ThisTime)
WT_data <- all_rmats_res_df[all_rmats_res_df$Sample=="WT",]
FL_data <- all_rmats_res_df[all_rmats_res_df$Sample=="FL",]

find_DE_mat <- function(rmats_res_df, rmats_data, time_pair){
  DE_mat <- matrix(0, nrow = nrow(rmats_data), ncol = length(time_pair))
  rownames(DE_mat) <- rmats_data$Key
  colnames(DE_mat) <- time_pair
  for(key_indx in 1:length(rmats_data$Key)){
    key = rmats_data$Key[key_indx]
    for(tp_indx in 1:length(time_pair)){
      tp <- time_pair[tp_indx]
      select_indx <- rmats_res_df$Key==key & rmats_res_df$TimePair==tp
      if(any(select_indx)){
        DE_mat[key_indx, tp_indx] <- rmats_res_df$DE[select_indx]
        # if(rmats_res_df$DE[select_indx] == "Up"){
        #   DE_mat[key_indx, tp_indx] <- 1
        # }else{
        #   if(rmats_res_df$DE[select_indx] == "Down"){
        #     DE_mat[key_indx, tp_indx] <- -1
        #   }
        # }
        
      }
    }
  }
  return(DE_mat)
}

find_delta_phi_mat <- function(rmats_res_df, rmats_data, time_pair){
  delta_phi_mat <- matrix(NA, nrow = nrow(rmats_data), ncol = length(time_pair))
  rownames(delta_phi_mat) <- rmats_data$Key
  colnames(delta_phi_mat) <- time_pair
  for(key_indx in 1:length(rmats_data$Key)){
    key = rmats_data$Key[key_indx]
    for(tp_indx in 1:length(time_pair)){
      tp <- time_pair[tp_indx]
      select_indx <- rmats_res_df$Key==key & rmats_res_df$TimePair==tp
      if(any(select_indx)){
        delta_phi_mat[key_indx, tp_indx] <- rmats_res_df$IncLevelDifference[select_indx]
      }
    }
  }
  return(delta_phi_mat)
}

find_phi_mat <- function(rmats_res_df, rmats_data, time_li){
  phi_mat <- matrix(NA, nrow = nrow(rmats_data), ncol = length(time_li))
  rownames(phi_mat) <- rmats_data$Key
  colnames(phi_mat) <- time_li
  for(key_indx in 1:length(rmats_data$Key)){
    key = rmats_data$Key[key_indx]
    for(tp_indx2 in 2:length(time_li)){
      tp_indx1 <- tp_indx2 - 1
      tp1 <- time_li[tp_indx1]
      tp2 <- time_li[tp_indx2]
      tp <- sprintf("%s-%s", tp1, tp2)
      select_indx <- rmats_res_df$Key==key & rmats_res_df$TimePair==tp
      if(any(select_indx)){
        phi1 <- sapply(strsplit(rmats_res_df$IncLevel1[select_indx], ","), function(x){return(mean(as.numeric(x), na.rm = T))})
        phi2 <- sapply(strsplit(rmats_res_df$IncLevel2[select_indx], ","), function(x){return(mean(as.numeric(x), na.rm = T))})
        phi_mat[key_indx, tp_indx1] <- phi1
        phi_mat[key_indx, tp_indx2] <- phi2
      }
    }
  }
  return(phi_mat)
}

WT_DE_mat <- find_DE_mat(WT_data, uniq_DE_rmats_res_df, time_pair)
FL_DE_mat <- find_DE_mat(FL_data, uniq_DE_rmats_res_df, time_pair)

WT_delta_phi_mat <- find_delta_phi_mat(WT_data, uniq_DE_rmats_res_df, time_pair)
FL_delta_phi_mat <- find_delta_phi_mat(FL_data, uniq_DE_rmats_res_df, time_pair)

WT_phi_mat <- find_phi_mat(WT_data, uniq_DE_rmats_res_df, args$tp)
FL_phi_mat <- find_phi_mat(FL_data, uniq_DE_rmats_res_df, args$tp)

compute_circadian_trans_num <- function(DE_li){
  up_indx <- as.integer(which(DE_li == "Up"))
  down_indx <- as.integer(which(DE_li == "Down"))
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
}

compute_odd_tag <- function(WT_DE_li, FL_DE_li){
  WT_up_indx <- as.integer(which(WT_DE_li == "Up"))
  WT_down_indx <- as.integer(which(WT_DE_li == "Down"))
  
  FL_up_indx <- as.integer(which(FL_DE_li == "Up"))
  FL_down_indx <- as.integer(which(FL_DE_li == "Down"))
  
  up_indx <- sort(unique(c(WT_up_indx, FL_up_indx)))
  down_indx <- sort(unique(c(WT_down_indx, FL_down_indx)))
  
  if( length(unique(up_indx %% 2)) > 1 ){
    return("NoTag")
  }
  if( length(unique(down_indx %% 2)) > 1 ){
    return("NoTag")
  }
  if( (length(up_indx)==0) | (length(down_indx)==0) ){
    return("NC")
  }
  if( unique(down_indx %% 2) == unique(up_indx %% 2) ){
    return("NoTag")
  }
  if(any(unique(up_indx %% 2) == 1)){
    return("Up")
  }else{
    return("Down")
  }
}

uniq_DE_rmats_res_df$WtCircadianTransNum <- apply(WT_DE_mat, 1, compute_circadian_trans_num)
uniq_DE_rmats_res_df$FlCircadianTransNum <- apply(FL_DE_mat, 1, compute_circadian_trans_num)
uniq_DE_rmats_res_df$OddTag <- sapply(1:nrow(uniq_DE_rmats_res_df), function(indx){
  return(compute_odd_tag(WT_DE_mat[indx,], FL_DE_mat[indx,]))
})

uniq_DE_rmats_res_df$WT_Circadian <- "NoCircadian"
uniq_DE_rmats_res_df$WT_Circadian[uniq_DE_rmats_res_df$WtCircadianTransNum==1] <- "WeakCircadian"
uniq_DE_rmats_res_df$WT_Circadian[uniq_DE_rmats_res_df$WtCircadianTransNum>2] <- "StrongCircadian"
uniq_DE_rmats_res_df$FL_Circadian <- "NoCircadian"
uniq_DE_rmats_res_df$FL_Circadian[uniq_DE_rmats_res_df$FlCircadianTransNum==1] <- "WeakCircadian"
uniq_DE_rmats_res_df$FL_Circadian[uniq_DE_rmats_res_df$FlCircadianTransNum>2] <- "StrongCircadian"

DEGs_state_info <- uniq_DE_rmats_res_df %>% group_by(AsType, WT_Circadian, FL_Circadian) %>% summarise(AsNum=n(), GeneNum=length(unique(GeneID)))
DEGs_state_info <- DEGs_state_info[order(DEGs_state_info$GeneNum, decreasing = T),]
write_tsv(DEGs_state_info, file.path(args$output, "AsCircadianDEGsStat.tsv"))

circadian_AS_gene_cond1 <- (uniq_DE_rmats_res_df$WT_Circadian=="StrongCircadian" | uniq_DE_rmats_res_df$FL_Circadian=="StrongCircadian")
circadian_AS_gene_cond2 <- uniq_DE_rmats_res_df$OddTag %in% c("Up", "Down")
circadian_AS_gene_indx <- circadian_AS_gene_cond1 & circadian_AS_gene_cond2

circadian_gene <- uniq_DE_rmats_res_df[circadian_AS_gene_indx, c("AsType", "GeneID", "Key")]
write_tsv(circadian_gene, file.path(args$output, "AllAsCircadianGene.txt"), col_names=FALSE)

circadian_gene_info <- circadian_gene %>% group_by(AsType) %>% summarise(AsNum=n(), GeneNum=length(unique(GeneID)))
write_tsv(circadian_gene_info, file.path(args$output, "AsCircadianGeneNum.tsv"))

circadian_gene_df <- uniq_DE_rmats_res_df[circadian_AS_gene_indx,]
circadian_AS_gene_WT_DE_mat <- WT_DE_mat[circadian_AS_gene_indx,]
circadian_AS_gene_FL_DE_mat <- FL_DE_mat[circadian_AS_gene_indx,]

circadian_AS_gene_WT_delta_phi_mat <- WT_delta_phi_mat[circadian_AS_gene_indx,]
circadian_AS_gene_FL_delta_phi_mat <- FL_delta_phi_mat[circadian_AS_gene_indx,]

circadian_AS_gene_WT_phi_mat <- WT_phi_mat[circadian_AS_gene_indx,]
circadian_AS_gene_FL_phi_mat <- FL_phi_mat[circadian_AS_gene_indx,]


merge_DE_mat <- cbind(circadian_AS_gene_WT_DE_mat, circadian_AS_gene_FL_DE_mat)
merge_DE_mat_num <- matrix(0, nrow = nrow(merge_DE_mat), ncol = ncol(merge_DE_mat))
merge_DE_mat_num[merge_DE_mat=="Up"] <- 1
merge_DE_mat_num[merge_DE_mat=="Down"] <- -1
colnames(merge_DE_mat_num) <- colnames(merge_DE_mat)
rownames(merge_DE_mat_num) <- rownames(merge_DE_mat)
merge_delta_phi_mat <- cbind(circadian_AS_gene_WT_delta_phi_mat, circadian_AS_gene_FL_delta_phi_mat)
merge_phi_mat <- cbind(circadian_AS_gene_WT_phi_mat, circadian_AS_gene_FL_phi_mat)

DE_dist <- dist(merge_DE_mat_num)
delta_phi_dist <- dist(merge_delta_phi_mat)
delta_phi_dist[is.na(delta_phi_dist)] <- median(delta_phi_dist, na.rm = T)
merge_dist <- DE_dist + delta_phi_dist
cl <- hclust(merge_dist, method = "ward.D2")

as_order <- order(circadian_gene_df$OddTag, apply(circadian_gene_df[,c("WtCircadianTransNum", "FlCircadianTransNum")], 1, max), decreasing = T)
de_ht <- Heatmap(
  merge_DE_mat_num,
  name = "DE",
  cluster_columns=FALSE,
  column_split=c(rep("WT DE", ncol(circadian_AS_gene_WT_DE_mat)), rep("FL DE", ncol(circadian_AS_gene_FL_DE_mat))),
  show_row_names=FALSE,
  row_order = cl$order,
  row_split = circadian_gene_df$AsType,
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
de_ht

psi_ht <- Heatmap(
  merge_phi_mat,
  name = "PSI",
  col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
  cluster_columns=FALSE,
  column_split=c(rep("WT PSI", ncol(merge_phi_mat)/2), rep("FL PSI", ncol(merge_phi_mat)/2)),
  show_row_names=FALSE,
  row_order = cl$order,
  row_split = circadian_gene_df$AsType,
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
psi_ht

delta_phi_ht <- Heatmap(
  merge_delta_phi_mat,
  name = "delta PSI",
  cluster_columns=FALSE,
  column_split=c(rep("WT delta PSI", ncol(merge_delta_phi_mat)/2), rep("FL delta PSI", ncol(merge_delta_phi_mat)/2)),
  show_row_names=FALSE,
  row_order = cl$order,
  row_split = circadian_gene_df$AsType,
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
delta_phi_ht

interest_gene <- read_delim(args$interest_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
names(interest_gene)[3] <- "GeneID"
interest_gene <- interest_gene[!duplicated(interest_gene[,c("ID", "GeneID")]),]
interest_circadian_df <- left_join(circadian_gene_df, interest_gene)
interest_circadian_df$Label <- sprintf("%s (%s)", interest_circadian_df$GeneID, interest_circadian_df$GeneSymbol)
subindx <- which(!is.na(interest_circadian_df$GeneSymbol))
label_li <- interest_circadian_df$Label[subindx]
ha <- rowAnnotation(link = row_anno_link(
  at = subindx,
  labels = label_li,
  labels_gp=gpar(fontsize = 6),
  lines_gp=gpar(lwd = 0.5),
  link_width = unit(4, "mm")),
  width = unit(4, "mm") + max_text_width(label_li, gp = gpar(fontsize = 6)))


inche_cm=2.54
pdf(file.path(args$output, "CircadianAsGene.HeatMap.pdf"), width=25/inche_cm, height=12/inche_cm, family="ArialMT", colormodel = "cmyk")
de_ht + delta_phi_ht + psi_ht + ha
dev.off()

psi_time1_df <- data.frame(
  Sample=all_rmats_res_df$Sample, TP=all_rmats_res_df$LastTime, AsType=all_rmats_res_df$AsType, 
  GeneID=all_rmats_res_df$GeneID, Key=all_rmats_res_df$Key,
  rep1=sapply(strsplit(all_rmats_res_df$IncLevel1, ","), function(x){return(as.numeric(x[1]))}),
  rep2=sapply(strsplit(all_rmats_res_df$IncLevel1, ","), function(x){return(as.numeric(x[2]))})
  )
psi_time2_df <- data.frame(
  Sample=all_rmats_res_df$Sample, TP=all_rmats_res_df$ThisTime, AsType=all_rmats_res_df$AsType, 
  GeneID=all_rmats_res_df$GeneID, Key=all_rmats_res_df$Key,
  rep1=sapply(strsplit(all_rmats_res_df$IncLevel2, ","), function(x){return(as.numeric(x[1]))}),
  rep2=sapply(strsplit(all_rmats_res_df$IncLevel2, ","), function(x){return(as.numeric(x[2]))})
)
psi_df <- rbind(psi_time1_df, psi_time2_df)
psi_df <- psi_df[!duplicated(psi_df[,c("Sample", "TP", "AsType", "GeneID", "Key")]),]
melt_psi_df <- melt(psi_df, id.vars = c("Sample", "TP", "AsType", "GeneID", "Key"), variable.name = "Rep", value.name = "PSI")
DE_AS_psi_df <- inner_join(circadian_gene_df, melt_psi_df)
DE_AS_psi_df <- na.omit(DE_AS_psi_df)

cpm_df <- read_delim(args$cpm, "\t", escape_double = FALSE, trim_ws = TRUE)
names(cpm_df)[1] <- "GeneID"
melt_cpm_df <- melt(cpm_df, id.vars = "GeneID", variable.name = "CpmLabel", value.name = "CPM")
melt_cpm_df$CpmLabel <- as.character(melt_cpm_df$CpmLabel)
melt_cpm_df$Sample <- sapply(strsplit(melt_cpm_df$CpmLabel, "[.]"), function(x){return(x[1])})
melt_cpm_df$TP <- sapply(strsplit(melt_cpm_df$CpmLabel, "[.]"), function(x){return(x[2])})
melt_cpm_df$Rep <- sapply(strsplit(melt_cpm_df$CpmLabel, "[.]"), function(x){return(x[3])})
melt_cpm_df$CpmLabel <- NULL
DE_AS_psi_df <- left_join(DE_AS_psi_df, melt_cpm_df)
DE_AS_psi_df$IncludeCPM <- DE_AS_psi_df$CPM * DE_AS_psi_df$PSI
DE_AS_psi_df$ExcludeCPM <- DE_AS_psi_df$CPM * (1 - DE_AS_psi_df$PSI)
write_tsv(DE_AS_psi_df, file.path(args$output, "AsCircadianPsiCpmSourceData.tsv"))

tp_df <- data.frame(
  TP_num=c(-48, -36, -24, -12, 0, 12, 24, 36, 48),
  TP=c("n48h", "n36h", "n24h", "n12h", "0h", "12h", "24h", "36h", "48h")
)

melt_DE_AS_psi_df <- melt(DE_AS_psi_df[,c("AsType", "GeneID", "Key", "Sample", "TP", "Rep", "PSI", "CPM", "IncludeCPM", "ExcludeCPM")], id.vars = c("AsType", "GeneID", "Key", "Sample", "TP", "Rep"), value.name = "Value", variable.name = "Stat")
melt_DE_AS_psi_df <- left_join(melt_DE_AS_psi_df, tp_df)
melt_DE_AS_psi_df$TP <- factor(melt_DE_AS_psi_df$TP, levels = args$tp)
melt_DE_AS_psi_df$Sample <- factor(melt_DE_AS_psi_df$Sample, levels = c("WT", "FL"))
melt_DE_AS_psi_df$Stat <- factor(
  melt_DE_AS_psi_df$Stat, 
  levels = c("PSI", "CPM", "IncludeCPM", "ExcludeCPM"), 
  labels = c("PSI", "Gene CPM", "Include AS CPM", "Exclude AS CPM"))

circadian_tp_df <- tp_df
circadian_tp_df$start_TP <- circadian_tp_df$TP_num - 2
circadian_tp_df$end_TP <- circadian_tp_df$TP_num + 10
circadian_tp_df$Label <- "Light"
circadian_tp_df$Label[(circadian_tp_df$TP_num %% 24)!=0] <- "Dark"

for(key in circadian_gene_df$Key){
  tmp_melt_DE_AS_psi_df <- melt_DE_AS_psi_df[melt_DE_AS_psi_df$Key==key,]
  melt_tmp_DE_AS_psi_info_df <- tmp_melt_DE_AS_psi_df %>% group_by(AsType, GeneID, Key, Sample, TP, TP_num, Stat) %>% summarise(n=n(), Mu=mean(Value), SD=sd(Value))
  melt_tmp_DE_AS_psi_info_df$ErrorBarMin <- melt_tmp_DE_AS_psi_info_df$Mu - melt_tmp_DE_AS_psi_info_df$SD
  melt_tmp_DE_AS_psi_info_df$ErrorBarMax <- melt_tmp_DE_AS_psi_info_df$Mu + melt_tmp_DE_AS_psi_info_df$SD
  melt_tmp_DE_AS_psi_info_df$Setoff <- 0.3
  melt_tmp_DE_AS_psi_info_df$Setoff[melt_tmp_DE_AS_psi_info_df$Sample=="WT"] <- -0.3
  
  p <- ggplot() +
    geom_line(melt_tmp_DE_AS_psi_info_df, mapping = aes(x=TP_num, y=Mu, color=Sample), size=0.5) +
    geom_point(tmp_melt_DE_AS_psi_df, mapping = aes(x=jitter(TP_num, factor = 0.2), y=Value, color=Sample), size=0.5) +
    geom_errorbar(melt_tmp_DE_AS_psi_info_df, mapping = aes(x=TP_num+Setoff, ymin=ErrorBarMin, ymax=ErrorBarMax, color=Sample), size=0.2, width=1) +
    facet_grid(Stat~., scales = "free_y") +
    scale_x_continuous(breaks = tp_df$TP_num, labels = tp_df$TP) +
    scale_color_manual(values = c(WT="red", FL="blue")) +
    scale_fill_manual(values = c(Light="yellow", Dark="grey40")) +
    labs(x="Time points", y="Value", title=sprintf("%s: %s", tmp_melt_DE_AS_psi_df$AsType[1], key)) +
    theme_bw() +
    theme(
      text = element_text(family="ArialMT", color = "black", size = 5),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.background = element_blank(),
      legend.key.size = unit(2, "mm"),
      legend.spacing = unit(1, "mm"),
      legend.title = element_blank()
    )
  f_out_dir <- file.path(args$output, "CircadianAsGene_plot", tmp_melt_DE_AS_psi_df$AsType[1])
  if(!dir.exists(f_out_dir)){
    dir.create(f_out_dir, recursive = T)
  }
  ggsave(file.path(f_out_dir, sprintf("%s.%s.pdf", tmp_melt_DE_AS_psi_df$GeneID[1], key)), p, width = 8, height = 7, limitsize = FALSE, units = "cm")
}



