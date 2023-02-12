library(readr)
library(argparse)
library(dplyr)
parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input.tbl")
parser$add_argument("--label", nargs="*", required=TRUE, type="character", dest = "label", metavar="label")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "out_data", metavar="FPKM.sample.tsv")
args <- commandArgs(TRUE)
args <- parser$parse_args(args) 

file_num <- length(args$input)
if( (length(args$label)!=file_num) ) stop()

load_data <- function(f, label){
  print(label)
  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
  names(df) <- c("Gid", label)
  return(df)
}

for(indx in 1:file_num){
  tmp_df <- load_data(args$input[indx], args$label[indx])
  if(indx == 1){
    df <- tmp_df
  }else{
    df <- inner_join(df, tmp_df, by=c("Gid"))
  }
}

write_tsv(df, args$out_data)

