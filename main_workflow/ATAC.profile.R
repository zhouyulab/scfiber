library(readr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
setwd("D:/cottonSingleCell")

ATAC_peak <- read_delim("stat/ATAC_MACS2/merge.peak.bed3", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(ATAC_peak) <- c("Chrom", "Start", "End")
ATAC_peak$Len <- ATAC_peak$End - ATAC_peak$Start

p <- ggplot(ATAC_peak, aes(x=Len)) +
  geom_histogram(fill="black", binwidth = 10) +
  labs(x="scATAC-seq peak length", y="#Peak") +
  theme_bw() +
  lims(x=c(0, 2500)) +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        panel.grid = element_blank())
ggsave("stat/ATAC_MACS2/ATAC.peak.len.pdf", p, width = 5, height = 4, units = "cm")

WT_promoter_profile <- read_delim("stat/ATAC_MACS2/profile/WT.promoter.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
FL_promoter_profile <- read_delim("stat/ATAC_MACS2/profile/FL.promoter.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
all(WT_promoter_profile$Id == FL_promoter_profile$Id)
WT_sig <- rowSums(WT_promoter_profile[,2:ncol(WT_promoter_profile)])
FL_sig <- rowSums(WT_promoter_profile[,2:ncol(FL_promoter_profile)])
profile_order <- order(WT_sig+FL_sig, decreasing = TRUE)
WT_promoter_profile <- WT_promoter_profile[profile_order,]
FL_promoter_profile <- FL_promoter_profile[profile_order,]
WT_sig <- rowSums(WT_promoter_profile[,2:ncol(WT_promoter_profile)])
FL_sig <- rowSums(WT_promoter_profile[,2:ncol(FL_promoter_profile)])
all_zero_indx <- (WT_sig+FL_sig) == 0
WT_promoter_profile <- WT_promoter_profile[!all_zero_indx,]
FL_promoter_profile <- FL_promoter_profile[!all_zero_indx,]

WT_mat <- as.matrix(WT_promoter_profile[,2:ncol(WT_promoter_profile)])
rownames(WT_mat) <- WT_promoter_profile$Id
FL_mat <- as.matrix(FL_promoter_profile[,2:ncol(FL_promoter_profile)])
rownames(FL_mat) <- FL_promoter_profile$Id

rm(WT_promoter_profile)
rm(FL_promoter_profile)
gc()

WT_pos_sum <- colSums(WT_mat)
FL_pos_sum <- colSums(FL_mat)

norm_WT_pos_sum <- length(WT_pos_sum) * WT_pos_sum / sum(WT_pos_sum)
norm_FL_pos_sum <- length(FL_pos_sum) * FL_pos_sum / sum(FL_pos_sum)

df <- data.frame(
  x=rep(as.integer(names(norm_WT_pos_sum)), 2),
  y=c(norm_WT_pos_sum, norm_FL_pos_sum), 
  label=c(rep("WT", length(norm_WT_pos_sum)), rep("FL", length(norm_FL_pos_sum)))
  )
df <- df[!df$x %in% c(-1000, 1000),]

df$label <- factor(df$label, levels = c("WT", "FL"))

p <- ggplot(df, aes(x=x, y=y, color=label)) +
  geom_line(size=0.3) +
  scale_x_continuous(breaks = c(-1000, 0, 1000), labels = c("-1000", "TSS", "1000"))+
  labs(y="ATAC signal") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", color = "black", size = 6),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.8, 0.8),
        panel.grid = element_blank())
ggsave("stat/ATAC_MACS2/ATAC.signal.profile.pdf", p, width = 5, height = 4, units = "cm")

