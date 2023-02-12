library(readr)
library(reshape2)
library(ggplot2)
CPM <- read_delim("analysis/RNA_seq_TP/circadian_gene/CPM.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
f_out <- "analysis/RNA_seq_TP/all_expr_plot"

CPM <- CPM[apply(CPM[,2:ncol(CPM)], 1, max) > 0,]

melt_cpm <- melt(CPM, "Gid", variable.name = "Label", value.name = "CPM")
melt_cpm$Label <- as.character(melt_cpm$Label)
melt_cpm$Sample <- sapply(strsplit(melt_cpm$Label, "[.]"), function(x){return(x[1])})
melt_cpm$TP <- sapply(strsplit(melt_cpm$Label, "[.]"), function(x){return(x[2])})
melt_cpm$Rep <- sapply(strsplit(melt_cpm$Label, "[.]"), function(x){return(x[3])})
melt_cpm$Chrom <- substr(melt_cpm$Gid, 1, 8)
melt_cpm$Chrom[!melt_cpm$Chrom%in%c(sprintf("Ghir_A%02d", 1:13), sprintf("Ghir_D%02d", 1:13))] <- "Other"

tp_df <- data.frame(
  TP_num=c(-48, -36, -24, -12, 0, 12, 24, 36, 48),
  TP=c("n48h", "n36h", "n24h", "n12h", "0h", "12h", "24h", "36h", "48h")
)
melt_cpm <- left_join(melt_cpm, tp_df)

chrom_li <- unique(melt_cpm$Chrom)


plot_expr <- function(chrom, f_out, tp_df, melt_cpm){
  tmp_dir <- file.path(f_out, chrom)
  if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir, recursive = T)
  }
  tmp_CPM_df <- melt_cpm[melt_cpm$Chrom == chrom,]
  gid_li <- unique(tmp_CPM_df$Gid)
  for(gid in gid_li){
    tmp_rep_cpm_df <- tmp_CPM_df[tmp_CPM_df$Gid==gid,]
    tmp_rep_cpm_info <- tmp_rep_cpm_df %>% group_by(Sample, TP_num) %>% summarise(Mu=mean(CPM), SD=sd(CPM))
    tmp_rep_cpm_info$Setoff <- 0.3
    tmp_rep_cpm_info$Setoff[tmp_rep_cpm_info$Sample=="WT"] <- -0.3
    tmp_rep_cpm_info$ErrorBarMax <- tmp_rep_cpm_info$Mu + tmp_rep_cpm_info$SD
    tmp_rep_cpm_info$ErrorBarMin <- tmp_rep_cpm_info$Mu - tmp_rep_cpm_info$SD
    tmp_rep_cpm_info$ErrorBarMin[tmp_rep_cpm_info$ErrorBarMin<0] <- 0
    max_scale <- max(tmp_rep_cpm_info$ErrorBarMax)
    
    circadian_tp_df <- tp_df
    circadian_tp_df$start_TP <- circadian_tp_df$TP_num - 2
    circadian_tp_df$end_TP <- circadian_tp_df$TP_num + 10
    circadian_tp_df$Label <- "Light"
    circadian_tp_df$Label[(circadian_tp_df$TP_num %% 24)!=0] <- "Dark"
    p <- ggplot() +
      geom_line(tmp_rep_cpm_info, mapping = aes(x=TP_num, y=Mu, color=Sample), size=0.5) +
      geom_point(tmp_rep_cpm_df, mapping = aes(x=jitter(TP_num, factor = 0.2), y=CPM, color=Sample), size=0.5) +
      geom_errorbar(tmp_rep_cpm_info, mapping = aes(x=TP_num+Setoff, ymin=ErrorBarMin, ymax=ErrorBarMax, color=Sample), size=0.2, width=1) +
      annotate("text", x=0, y=-0.03*max_scale, color="black", label="Light", size=2) +
      geom_rect(data = circadian_tp_df, mapping=aes(xmin=start_TP, xmax=end_TP, ymax=-1*0.07*max_scale, ymin=-1*0.1*max_scale, fill=Label)) +
      annotate("text", x=0, y=-0.16*max_scale, color="black", label="WT DE", size=2) +
      annotate("text", x=0, y=-0.32*max_scale, color="black", label="FL DE", size=2) + 
      scale_x_continuous(breaks = tp_df$TP_num, labels = tp_df$TP) +
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
    ggsave(file.path(tmp_dir, sprintf("%s.pdf", gid)), p, width = 10, height = 4, units = "cm") 
  }
}
bplapply(chrom_li, plot_expr, f_out=f_out, tp_df=tp_df, melt_cpm=melt_cpm, BPPARAM=MulticoreParam(40))
