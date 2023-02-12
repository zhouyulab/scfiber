library(RColorBrewer)
library(scran)
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(readr)

args = commandArgs(TRUE)
if (length(args) == 3) {
  f_corSCE <- args[1]
  f_cluster <- args[2]
  f_out <- args[3]
} else {
  q()
}

load(f_corSCE)
cluster <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)

if(!all(cluster$Cell == colnames(corSCE))) stop()

corSCE$cluster <- cluster$Cluster


## Cluster UMAP
uniq_clu <- sort(unique(cluster$Cluster))
clu_num <- length(uniq_clu)

if(clu_num <= 8){
  cb_palette <- brewer.pal(8, "Dark2")
}else{
  cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                  "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                  "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                  "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                  "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                  "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
}


if(clu_num>length(cb_palette)) stop()

col_li <- cb_palette[1:clu_num]
names(col_li) <- uniq_clu

plot_cluster_umap  <- function(umapSCE, col_li, label){
  df <- data.frame(
    x=umapSCE@int_colData$reducedDims$UMAP[,1], 
    y=umapSCE@int_colData$reducedDims$UMAP[,2], 
    cluster=umapSCE$cluster
  )
  cell_num <- nrow(df)
  rand_indx <- sample(1:cell_num, cell_num)
  df <- df[rand_indx,]
  p <- ggplot(df, aes(x=x, y=y, color=cluster)) +
    theme_bw() +
    geom_point(size=0.1, alpha=0.6) +
    labs(x="UMAP1", y="UMAP2", alpha="UMI", title=label) +
    scale_color_manual(values = col_li) +
    theme(
      text = element_text(family="ArialMT", size = 6),
      title = element_text(family="ArialMT", size = 6),
      legend.title = element_text(family="ArialMT", size = 6),
      legend.text = element_text(family="ArialMT", size = 6),
      panel.grid = element_blank(),
      legend.background = element_blank()
    )
}


ps_clu <- list()
p_clu_all <- plot_cluster_umap(corSCE, col_li, "All")
ps_clu[["All"]] <- p_clu_all

samples <- sort(unique(corSCE$SampleName))
for(s in samples){
  tmp_SCE <- corSCE[,corSCE$SampleName==s]
  ps_clu[[s]] <- plot_cluster_umap(tmp_SCE, col_li, s)
}

inche_cm=2.54
pdf(f_out, width=8*3/inche_cm, height=6*ceiling(length(ps_clu)/3)/inche_cm, family="ArialMT", colormodel = "cmyk")
grid.arrange(grobs=ps_clu,ncol=3,padding = unit(0, "mm"))
dev.off()
