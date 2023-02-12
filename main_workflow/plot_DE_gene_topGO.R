library(readr)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(topGO)
library(RColorBrewer)

pval_cutoff <- 0.05

clu <- c("C2", "C3")
top_df <- data.frame()
enrigh_term_df <- data.frame()
for(tmp_clu in clu){
  tmp_file <- file.path("D:/cottonSingleCell/stat/diff_expr/DiffExpr.C2.topGO", sprintf("%s.enrich", tmp_clu), "topGO.tsv")
  if(!file.exists(tmp_file)) next()
  tmp_topGO <- read_delim(tmp_file, "\t", escape_double = FALSE, col_types = cols(classicFisher = col_double()), trim_ws = TRUE)
  tmp_topGO$classicFisher[is.na(tmp_topGO$classicFisher)] <- 1e-30
  tmp_topGO$Label <- sprintf("%s (%d)", tmp_topGO$Term, tmp_topGO$Significant)
  tmp_topGO$Term <- factor(tmp_topGO$Term, levels = rev(tmp_topGO$Term[!duplicated(tmp_topGO$Term)]))
  tmp_topGO$Cluster <- tmp_clu
  tmp_topGO$classicKS <- NULL
  tmp_topGO$elimKS <- NULL
  enrigh_term_df <- rbind(enrigh_term_df, tmp_topGO[tmp_topGO$classicFisher<pval_cutoff,])
  tmp_topGO <- plyr::ddply(tmp_topGO, "Ontology", function(block){return(block[1:10,])})
  tmp_topGO$Indx <- 1:nrow(tmp_topGO)
  top_df <- rbind(top_df, tmp_topGO)
}


top_df <- na.omit(top_df)
p <- ggplot(top_df, aes(x=-Indx, fill=Ontology, y=-log10(classicFisher))) +
  theme_bw() +
  facet_grid(~Cluster) +
  geom_bar(stat="identity", width = 0.6) +
  geom_text(aes(label=Label), size=2, hjust=0) +
  scale_fill_brewer(palette = "Dark2") +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 10, 20, 30), labels = c("1E0", "1E-10", "1E-20", "<1E-30"), limits = c(0, 60)) +
  labs(y="p-value") +
  theme(
    axis.text = element_text(family="ArialMT", color = "black", size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
  )
p
ggsave("D:/cottonSingleCell/stat/diff_expr/DiffExpr.C2.topGO.pdf", p, width = 15*2, height = 8, units = "cm")

