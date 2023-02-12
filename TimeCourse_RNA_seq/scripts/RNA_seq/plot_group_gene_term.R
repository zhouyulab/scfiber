library(readr)
library(argparse)
library(dplyr)
library(ggplot2)
parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="term.tsv")
parser$add_argument("--label", nargs="*", required=TRUE, type="character", dest = "label", metavar="label")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "out_plot", metavar="term.pdf")
args <- commandArgs(TRUE)
args <- parser$parse_args(args) 

# args <- c(
#   "-i",
#   "analysis/RNA_seq_TP/TimeCourse/group_topGO/Early/topGO.tsv",
#   "analysis/RNA_seq_TP/TimeCourse/group_topGO/Late/topGO.tsv",
#   "--label",
#   "Early", "Late", 
#   "-o", "analysis/RNA_seq_TP/TimeCourse/group_topGO/topGO.pdf"
# )
# args <- parser$parse_args(args)

read_term <- function(f, label, pval_cutoff=0.05, top_n=10){
  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)
  df$Label <- label
  df <- df %>% group_by(Ontology) %>% slice_head(n=top_n) %>% mutate(Indx=n():1)
  df$classicFisher[df$classicFisher=="< 1e-30"] <- "1e-30"
  df$classicFisher <- as.numeric(df$classicFisher)
  df <- df[df$classicFisher < pval_cutoff,]
  df <- as.data.frame(df)
  return(df)
}

term_li <- list()
for(i in 1:length(args$input)){
  term_li[[args$label[i]]] <- read_term(args$input[i], args$label[i])
}
all_term_df <- do.call(rbind, term_li)
all_term_df$Label <- factor(all_term_df$Label, levels = args$label)
ymax <- max(-log10(all_term_df$classicFisher))
ybreaks <- seq(0, ymax, floor(ymax/4.))
p <- ggplot(all_term_df, aes(x=Indx, fill=Ontology, y=-log10(classicFisher))) +
  theme_bw() +
  geom_bar(stat="identity", width = 0.6) +
  geom_text(aes(label=Term), size=2, hjust=0) +
  scale_fill_brewer(palette = "Dark2") +
  coord_flip() +
  scale_y_continuous(breaks = ybreaks, limits = c(0, ymax*2), labels = sprintf("%.0E", 10 ** ybreaks)) +
  facet_grid(Ontology~Label) +
  labs(y="p-value") +
  theme(
    axis.text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line()
  )
ggsave(args$out_plot, p, width = 6*length(args$input), height = 12, limitsize = FALSE, units = "cm")

