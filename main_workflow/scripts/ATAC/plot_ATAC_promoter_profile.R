library(ggplot2)
library(readr)
library(dplyr)
library(reshape2)

args = commandArgs(TRUE)
if (length(args) == 3) {
  f_WT <- args[1]
  f_FL <- args[2]
  f_out <- args[3]
} else {
  q()
}

# setwd("~/mu01/project/CottonSingleCell")

profile2df <- function(profile, label){
  profile_mat <- as.matrix(profile[, 3:(ncol(WT_promoter_profile)-1)])
  x <- as.integer(colnames(profile_mat))
  y <- colMeans(profile_mat)
  df <- data.frame(x=x, y=y, label=label)
}

WT_promoter_profile <- read_delim(f_WT, "\t", escape_double = FALSE, trim_ws = TRUE)
FL_promoter_profile <- read_delim(f_FL, "\t", escape_double = FALSE, trim_ws = TRUE)

WT_promoter_profile_df <- profile2df(WT_promoter_profile, "WT")
FL_promoter_profile_df <- profile2df(FL_promoter_profile, "FL")

promoter_df <- rbind(WT_promoter_profile_df, FL_promoter_profile_df)

p <- ggplot(promoter_df, aes(x=x, y=y, color=label)) +
  geom_line(size=0.3) +
  scale_x_continuous(breaks = c(-1000, 0, 1000), labels = c("-1000", "TSS", "1000")) +
  scale_color_manual(values = c("WT"="blue", "FL"="red")) +
  theme_bw() +
  theme_bw() +
  labs(y="ATAC-seq signal") +
  theme(text = element_text(family="ArialMT", color = "black", size = 5),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        panel.grid = element_blank())
ggsave(f_out, p, width = 5, height = 3.8, units = "cm")
