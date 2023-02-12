library(readr)
library(argparse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input.tbl")
parser$add_argument("-t", "--timepoint", nargs="*", required=TRUE, type="character", dest = "timepoint", metavar="timepoint")
parser$add_argument("-s", "--srr", nargs="*", required=TRUE, type="character", dest = "srr", metavar="SRR")
parser$add_argument("--out-data", required=TRUE, type="character", dest = "out_data", metavar="FPKM.sample.tsv")
parser$add_argument("--out-plot", required=TRUE, type="character", dest = "out_plot", metavar="cor.pdf")
parser$add_argument("--out-merge", required=TRUE, type="character", dest = "out_merge", metavar="FPKM.timepoint.tsv")
args <- commandArgs(TRUE)

# args <- c(
#   "-i", 
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/RNA_seq/stringtie/0DPA.SRR1695183.tab",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/RNA_seq/stringtie/0DPA.SRR2081051.tab",
#   "/home/sqreb/ias/data/CottonSingleCell/analysis/RNA_seq/stringtie/1DPA.SRR1695184.tab",
#   "-t", "0DPA", "0DPA", "1DPA",
#   "-s", "SRR1695183", "SRR2081051", "SRR1695184",
#   "--out-plot", "~/test.cor.pdf",
#   "--out-data", "~/test.tsv",
#   "--out-merge", "~/test.merge.tsv"
#   )

args <- parser$parse_args(args) 

file_num <- length(args$input)
if((length(args$timepoint)!=file_num) | (length(args$srr)!=file_num)) stop()


load_data <- function(f, tp, srr){
  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)
  df <- df[,c("Gene ID", "FPKM")]
  names(df) <- c("Gid", sprintf("%s.%s", tp, srr))
  return(df)
}

for(indx in 1:file_num){
  tmp_df <- load_data(args$input[indx], args$timepoint[indx], args$srr[indx])
  if(indx == 1){
    df <- tmp_df
  }else{
    df <- left_join(df, tmp_df, by="Gid")
  }
}

write_tsv(df, args$out_data)

expr_mat <- as.matrix(df[,2:ncol(df)])
cor_mat <- cor(expr_mat, method = "spearman")
cor_mat_text <- matrix(sprintf("%.2f" ,cor_mat), nrow = nrow(cor_mat))

inche_cm=2.54
pdf(args$out_plot, width=20/inche_cm, height=18/inche_cm, family="ArialMT", colormodel = "cmyk")
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
draw(ht)
dev.off()

merge_df <- data.frame(Gid=df$Gid)
for(tp in args$timepoint[!duplicated(args$timepoint)]){
  hit_indx <- args$timepoint==tp
  if(sum(hit_indx)>1){
    merge_df[[tp]] <- rowMeans(expr_mat[,hit_indx])
  }else{
    merge_df[[tp]] <- expr_mat[,hit_indx]
  }
  
}
write_tsv(merge_df, args$out_merge)