library(readr)
library(dplyr)
library(DropletUtils)
library(reshape2)
library(zeallot)
library(ggplot2)

args = commandArgs(TRUE)
if (length(args) == 3) {
  f_mergeSCE <- args[1]
  f_cluster <- args[2]
  f_out <- args[3]
} else {
  q()
}
# setwd("~/ias/data/CottonSingleCell")
# f_mergeSCE <- "analysis/base_map/mergeSCE.RData"
# f_cluster <- "analysis/base_map/cluster.tsv"

load(f_mergeSCE)
clusterSCE <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)
mergeSCE$cluster <- clusterSCE$Cluster
cluster_level = c("C1", "C2", "C3", "C4", "C5")

## Cell number
cell_clu_num <- clusterSCE %>% group_by(Cluster) %>% summarise(CellNumber=n())
clusterSCE$Cluster <- factor(clusterSCE$Cluster, levels = cluster_level)

cell_num <- clusterSCE %>% group_by(Sample, Cluster) %>% summarise(CellNumber=n())
cell_num$Cluster <- factor(cell_num$Cluster, levels = cluster_level)
sample_num <- clusterSCE %>% group_by(Sample) %>% summarise(SampleCellNumber=n())
sc_stat <- left_join(cell_num, sample_num, by="Sample")
sc_stat$CellRatio <- sc_stat$CellNumber / sc_stat$SampleCellNumber
sc_stat$CellRatioText <- sprintf("%.2f%%", sc_stat$CellRatio * 100)

## UMI
cell_umi <- colSums(counts(mergeSCE))
cell_umi <- data.frame(Sample=mergeSCE$SampleName, Rep=mergeSCE$Rep, Barcode=names(cell_umi), UMI=cell_umi)
cell_umi <- left_join(clusterSCE, cell_umi, by=c("Sample", "Rep", "Barcode"))
cell_clu_umi <- cell_umi %>% group_by(Sample, Cluster) %>% summarise(TotalUMI=sum(UMI), MedianUMIPerCell=median(UMI), MadUMIPerCell=mad(UMI))
sample_umi <- cell_umi %>% group_by(Sample) %>% summarise(SampleTotalUMI=sum(UMI))
sc_stat <- left_join(sc_stat, cell_clu_umi, by=c("Sample", "Cluster"))
sc_stat <- left_join(sc_stat, sample_umi, by="Sample")
sc_stat$TotalUMIRatio <- sc_stat$TotalUMI / sc_stat$SampleTotalUMI
sc_stat$TotalUMIRatioText <- sprintf("%.2f%%", sc_stat$TotalUMIRatio * 100)

## Gene number
cell_gene_num <- c()
detected_gene_li <- list()
for(tmp_batch in unique(mergeSCE$Sample)){
  tmmSCE <- mergeSCE[, mergeSCE$Sample==tmp_batch]
  tmp_cnt <- normcounts(tmmSCE)
  cell_gene_num <- c(cell_gene_num, colSums(tmp_cnt>0))
  c(tmp_sample, tmp_rep) %<-% strsplit(tmp_batch, "[.]")[[1]]
  tmp_barcode_df <- clusterSCE[which((clusterSCE$Sample==tmp_sample) & (clusterSCE$Rep==tmp_rep)),]
  for(tmp_clu in unique(tmp_barcode_df$Cluster)){
    tmp_barcode_li <- tmp_barcode_df$Barcode[tmp_barcode_df$Cluster==tmp_clu]
    if(length(tmp_barcode_li)>=2){
      tmp_mat <- tmp_cnt[,colnames(tmp_cnt) %in% tmp_barcode_li]
      tmp_detected_gene <- rownames(tmp_mat)[rowSums(tmp_mat)>0]
    }else{
      if(length(tmp_barcode_li)==0){
        tmp_detected_gene <- c()
      }else{
        tmp_detected_gene <- names(tmp_mat)[tmp_mat>0]
      }
    }
    tmp_key <- sprintf("%s.%s", tmp_sample, tmp_clu)
    if(! tmp_key %in% names(detected_gene_li)){
      detected_gene_li[[tmp_key]] <- c()
    }
    detected_gene_li[[tmp_key]] <- unique(c(detected_gene_li[[tmp_key]], tmp_detected_gene))
  }
}
clu_gene_num <- unlist(lapply(detected_gene_li, length))
clu_gene_num <- data.frame(Class=names(clu_gene_num), TotalGeneNum=clu_gene_num)
clu_gene_num$Sample <- sapply(strsplit(as.character(clu_gene_num$Class), "[.]"), function(x) x[[1]])
clu_gene_num$Cluster <- sapply(strsplit(as.character(clu_gene_num$Class), "[.]"), function(x) x[[2]])
clu_gene_num[, "Class"] <- NULL
sc_stat <- left_join(sc_stat, clu_gene_num, by=c("Sample", "Cluster"))

cell_gene_num <- data.frame(Barcode=names(cell_gene_num), GeneNum=cell_gene_num)
cell_gene_num <- left_join(clusterSCE, cell_gene_num, by="Barcode")
cell_clu_gene_num <- cell_gene_num %>% group_by(Sample, Cluster) %>% summarise(MedianGeneNumPerCell=median(GeneNum), MadGeneNumPerCell=mad(GeneNum))
sc_stat <- left_join(sc_stat, cell_clu_gene_num, by=c("Sample", "Cluster"))
sc_stat$Cluster <- factor(sc_stat$Cluster, levels = cluster_level)

umi_per_gene <- table(mergeSCE@assays@data@listData$counts@x)
umi_per_gene_df <- as.data.frame(umi_per_gene)
names(umi_per_gene_df) <- c("UMI", "GeneNum")
umi_per_gene_df$UMI <- as.integer(as.character(umi_per_gene_df$UMI))
umi_per_gene_df <- left_join(data.frame(UMI=1:max(umi_per_gene_df$UMI)), umi_per_gene_df, by="UMI")
umi_per_gene_df$GeneNum[is.na(umi_per_gene_df$GeneNum)] <- 0
p <- ggplot(umi_per_gene_df, aes(x=UMI, y=log10(GeneNum))) +
  geom_bar(color="black", stat = "identity") +
  theme_bw() +
  labs(x="UMI / (cell x gene)", y="Number") +
  xlim(0, 600) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c("1e0", "1e2", "1e4", "1e6", "1e8", "1e10")) +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "UMIPerGeneCell.BarPlot.pdf"), p, height = 4, width = 6, units = "cm", colormodel = "cmyk")

fn <- ecdf(mergeSCE@assays@data@listData$normcounts@x)
cum_umi_per_gene <- data.frame(UMI=1:max(umi_per_gene_df$UMI), cumSum=fn(1:max(umi_per_gene_df$UMI)))
p <- ggplot(cum_umi_per_gene, aes(x=UMI, y=cumSum)) +
  theme_bw() +
  geom_line() +
  labs(x="UMI / (cell x gene)", y="Cumulative proportion") +
  xlim(0, 600) +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "UMIPerGeneCell.CumSum.pdf"), p, height = 4, width = 6, units = "cm", colormodel = "cmyk")

p <- ggplot(sc_stat, aes(x=Cluster, y=CellNumber, fill=Sample)) +
  geom_bar(stat="identity", width=0.8, position=position_dodge(0.8)) +
  geom_text(aes(label=CellNumber), position=position_dodge(0.8), size=1.2, vjust=-0.2) +
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000), labels = c("0", "5,000", "10,000", "15,000", "20,000")) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y="#Cell") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "Cluster.CellNumber.pdf"), p, height = 4, width = 10, units = "cm", colormodel = "cmyk")

p <- ggplot(sc_stat, aes(x=Cluster, y=CellRatio, fill=Sample)) +
  geom_bar(stat="identity", width=0.8, position=position_dodge(0.8)) +
  geom_text(aes(label=CellRatioText), position=position_dodge(0.8), size=1.2, vjust=-0.2) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0%", "20%", "40%", "60%", "80%", "100%")) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y="Cell ratio") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "Cluster.CellRatio.pdf"), p, height = 4, width = 10, units = "cm", colormodel = "cmyk")

p <- ggplot(sc_stat, aes(x=Cluster, y=TotalUMI, fill=Sample)) +
  geom_bar(stat="identity", width=0.8, position=position_dodge(0.8)) +
  geom_text(aes(label=sprintf("%.2f", TotalUMI/1e7)), position=position_dodge(0.8), size=1.2, vjust=-0.2) +
  scale_y_continuous(breaks = seq(0, 20, 2) * 1e7, labels = as.character(seq(0, 20, 2))) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y="Total UMI x1e6") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "Cluster.TotalUMI.pdf"), p, height = 4, width = 10, units = "cm", colormodel = "cmyk")

p <- ggplot(sc_stat, aes(x=Cluster, y=TotalUMIRatio, fill=Sample)) +
  geom_bar(stat="identity", width=0.8, position=position_dodge(0.8)) +
  geom_text(aes(label=TotalUMIRatioText), position=position_dodge(0.8), size=1.2, vjust=-0.2) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0%", "20%", "40%", "60%", "80%", "100%")) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y="UMI ratio") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "Cluster.TotalUMIRatio.pdf"), p, height = 4, width = 10, units = "cm", colormodel = "cmyk")

p <- ggplot(cell_umi, aes(x=Cluster, y=log10(UMI), fill=Sample)) +
  geom_violin(draw_quantiles = 0.5) +
  # geom_text(data = sc_stat, aes(y=max(log10(cell_umi$UMI)) ,label=CellNumber), position=position_dodge(0.8), size=1.2, vjust=-0.2) +
  scale_y_continuous(breaks = c(3, log10(c(2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000)), 4), labels = c("1e3", rep("", 8), "1e4")) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y="UMI per cell") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "Cluster.UMIPerCell.pdf"), p, height = 4, width = 10, units = "cm", colormodel = "cmyk")

p <- ggplot(sc_stat, aes(x=Cluster, y=TotalGeneNum, fill=Sample)) +
  geom_bar(stat="identity", width=0.8, position=position_dodge(0.8)) +
  geom_text(aes(label=TotalGeneNum), position=position_dodge(0.8), size=1.2, vjust=-0.2) +
  scale_y_continuous(breaks = seq(0, 20, 2) * 1e4, labels = as.character(seq(0, 20, 2))) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y="#Detected gene x10,000") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "Cluster.DetectedGeneNumber.pdf"), p, height = 4, width = 10, units = "cm", colormodel = "cmyk")

p <- ggplot(cell_gene_num, aes(x=Cluster, y=log10(GeneNum), fill=Sample)) +
  geom_violin() +
  # geom_text(data = sc_stat, aes(y=max(log10(cell_umi$UMI)) ,label=CellNumber), position=position_dodge(0.8), size=1.2, vjust=-0.2) +
  scale_y_continuous(breaks = c(2, log10(c(200, 500)), 3, log10(c(2000, 5000))), labels = c("100", "200", "500", "1,000", "2,000", "5,000")) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y="#Gene per cell") +
  theme(text = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(f_out, "Cluster.GeneNumPerCell.pdf"), p, height = 4, width = 10, units = "cm", colormodel = "cmyk")

col_names <- c(
  "Sample", "Cluster", 
  "CellNumber", "CellRatio", 
  "TotalUMI", "TotalUMIRatio", 
  "MedianUMIPerCell", 
  "TotalGeneNum",
  "MedianGeneNumPerCell"
)
sc_stat_res <- sc_stat[,col_names]

sample_stat <- sc_stat_res %>% group_by(Sample) %>% summarise(
  CellNumber=sum(CellNumber), CellRatio=1, TotalUMI=sum(TotalUMI), TotalUMIRatio=1
  )

sample_umi <- cell_umi %>% group_by(Sample) %>% summarise(MedianUMIPerCell=median(UMI), MadUMIPerCell=mad(UMI))
sample_stat <- left_join(sample_stat, sample_umi, by="Sample")

sample_detected_gene_li <- list()
all_detected_gene <- c()
for(key in names(detected_gene_li)){
  c(tmp_sample, tmp_cluster) %<-% strsplit(key, "[.]")[[1]]
  if(! tmp_sample %in% names(sample_detected_gene_li)){
    sample_detected_gene_li[[tmp_sample]] <- c()
  }
  tmp_gene <- detected_gene_li[[key]]
  sample_detected_gene_li[[tmp_sample]] <- unique(c(sample_detected_gene_li[[tmp_sample]], tmp_gene))
  all_detected_gene <- unique(c(all_detected_gene, tmp_gene))
}
sample_gene_num <- unlist(lapply(sample_detected_gene_li, length))
sample_gene_num <- data.frame(Sample=names(sample_gene_num), TotalGeneNum=sample_gene_num)
sample_stat <- left_join(sample_stat, sample_gene_num, by="Sample")

sample_gene_num <- cell_gene_num %>% group_by(Sample) %>% summarise(MedianGeneNumPerCell=median(GeneNum), MadGeneNumPerCell=mad(GeneNum))
sample_stat <- left_join(sample_stat, sample_gene_num, by="Sample")
sample_stat$Cluster <- "ALL"
sample_stat <- sample_stat[,col_names]

all_stat <- data.frame(
  Sample="ALL",
  Cluster="ALL",
  CellNumber=sum(sample_stat$CellNumber),
  CellRatio=NA,
  TotalUMI=sum(sample_stat$TotalUMI),
  TotalUMIRatio=NA,
  MedianUMIPerCell=median(cell_umi$UMI),
  MadUMIPerCell=mad(cell_umi$UMI),
  TotalGeneNum=length(all_detected_gene),
  MedianGeneNumPerCell=median(cell_gene_num$GeneNum),
  MadGeneNumPerCell=mad(cell_gene_num$GeneNum)
)
all_stat <- all_stat[,col_names]

sc_stat_res <- as.data.frame(sc_stat_res)
sample_stat <- as.data.frame(sample_stat)
all_stat_res <- rbind(sc_stat_res, sample_stat, all_stat)
write_tsv(all_stat_res, file.path(f_out, "Cluster.stat.tsv"))
