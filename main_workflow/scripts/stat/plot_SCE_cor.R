library(readr)
library(reshape2)
library(zeallot)
library(ComplexHeatmap)
library(circlize)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input")
parser$add_argument("--rep1", nargs="*", required=TRUE, type="character", dest = "rep1", metavar="rep1")
parser$add_argument("--rep2", nargs="*", required=TRUE, type="character", dest = "rep2", metavar="rep2")
parser$add_argument("--rep1-name", nargs="*", required=TRUE, type="character", dest = "rep1_name", metavar="rep1_name")
parser$add_argument("--rep2-name", nargs="*", required=TRUE, type="character", dest = "rep2_name", metavar="rep2_name")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
args <- commandArgs(TRUE)

# args <- c(
#   "-i", "~/sample_rep_cnt.tsv", 
#   "--rep1", "rep1", "--rep2", "rep2", 
#   "--rep1-name", "rep1", "--rep2-name", "rep2",
#   "-o", "~/"
#   )

args <- parser$parse_args(args) 

cnt <- read_delim(args$input, "\t", escape_double = FALSE, trim_ws = TRUE)
cnt_mat <- as.matrix(cnt[2:ncol(cnt)])

batchs <- names(cnt)[2:ncol(cnt)]
sample_li <- sapply(strsplit(as.character(batchs), "[.]"), function(x) x[[1]])
uniq_sample <- unique(sample_li)

build_cnt_mat <- function(cnt_mat, sample, select_rep, rep_name){
  batchs <- colnames(cnt_mat)
  sample_li <- sapply(strsplit(as.character(batchs), "[.]"), function(x) x[[1]])
  rep_li <- sapply(strsplit(as.character(batchs), "[.]"), function(x) x[[2]])
  clu_li <- sapply(strsplit(as.character(batchs), "[.]"), function(x) x[[3]])
  uniq_clu <- sprintf("C%d", 1:length(unique(clu_li)))
  if(!all(uniq_clu %in% clu_li)) stop()
  tmp_cnt_mat <- matrix(0, nrow = nrow(cnt_mat), ncol = length(uniq_clu))
  for(indx in 1:length(uniq_clu)){
    select_indx <- (sample_li == sample) & (rep_li %in% select_rep) & (clu_li == uniq_clu[indx])
    tmp_mat <- cnt_mat[,select_indx]
    if(length(select_rep)==1){
      tmp_cnt_mat[,indx] <- tmp_mat
    }else{
      tmp_cnt_mat[,indx] <- rowSums(tmp_mat)
    }
  }
  # colnames(tmp_cnt_mat) <- sprintf("%s.%s", rep_name, uniq_clu)
  colnames(tmp_cnt_mat) <- uniq_clu
  rownames(tmp_cnt_mat) <- rownames(cnt_mat)
  return(tmp_cnt_mat)
}


for(sample in uniq_sample){
  tmp_rep1_mat <- build_cnt_mat(cnt_mat, sample, args$rep1, args$rep1_name)
  tmp_rep2_mat <- build_cnt_mat(cnt_mat, sample, args$rep2, args$rep2_name)
  cor_mat <- cor(log10(tmp_rep2_mat+1), log10(tmp_rep1_mat+1), method="pearson")
  cor_mat_text <- matrix(sprintf("%.2f" ,cor_mat), nrow = nrow(cor_mat))
  ht <- Heatmap(cor_mat, name="Correlation",
                col = colorRamp2(c(0.5, 0.75, 1), c("#2c7bb6", "#ffffbf", "#d7191c")),
                row_title_gp = gpar(fontsize = 6),
                row_names_gp = gpar(fontsize = 6),
                column_title_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_dend_width = unit(5, "mm"),
                column_dend_height = unit(5, "mm"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                column_title = args$rep1_name,
                row_title = args$rep2_name,
                cell_fun = function(j, i, x, y, w, h, col) {
                  grid.text(cor_mat_text[i, j], x, y, gp = gpar(fontsize = 5))
                },
                heatmap_legend_param = list(
                  labels_gp = gpar(fontsize = 6),
                  title_gp = gpar(fontsize = 6),
                  grid_width = unit(2, "mm"),
                  grid_height = unit(2, "mm"))
  )
  inche_cm=2.54
  pdf(file.path(args$output, sprintf("%s.cor.pdf", sample)), width=10/inche_cm, height=8/inche_cm, family="ArialMT", colormodel = "cmyk")
  draw(ht)
  dev.off()
  write.csv(cor_mat,file.path(args$output, sprintf("%s.cor.csv", sample)))
}
