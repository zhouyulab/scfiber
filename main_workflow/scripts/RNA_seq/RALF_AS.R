library(argparse)
library(dplyr)
library(readr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", required=TRUE, type="character", dest = "input", metavar="input_dir")
parser$add_argument("--FDR-cutoff", type="double", dest = "FDR_cutoff", metavar="FDR_cutoff", default=0.05)
parser$add_argument("--delta-phi-cutoff", type="double", dest = "delta_phi_cutoff", metavar="delta_phi_cutoff", default=0.2)
parser$add_argument("--read-num-cutoff", type="integer", dest = "read_num_cutoff", metavar="read_num_cutoff", default=20)
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")

args <- commandArgs(TRUE)
args <- parser$parse_args(args)

args <- c(
  "-i", "analysis_v3/RALF_RNA_seq/rMATS_turbo",
  "-o", "analysis_v3/RALF_RNA_seq/AS"
)
args <- parser$parse_args(args)

AS_type_li <- c("A3SS", "A5SS", "RI", "SE")

load_data <- function(fname, AS_type){
  df <- read_delim(fname, "\t", escape_double = FALSE, trim_ws = TRUE)
  df$chr <- sapply(strsplit(df$chr, "chr"), function(x){return(x[2])})
  if(AS_type %in% c("SE")){
    df$Key <- sprintf("%s_%s_%s_%s_%s_%s_%s_%s", df$chr, df$exonStart_0base, df$exonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, df$strand)
  }else{
    if(AS_type %in% c("A3SS", "A5SS")){
      df$Key <- sprintf("%s_%s_%s_%s_%s_%s_%s_%s", df$chr, df$longExonStart_0base, df$longExonEnd, df$shortES, df$shortEE, df$flankingES, df$flankingEE, df$strand)
    }else{
      if(AS_type %in% c("RI")){
        df$Key <- sprintf("%s_%s_%s_%s_%s_%s_%s_%s", df$chr, df$riExonStart_0base, df$riExonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, df$strand)
      }
    }
  }
  df$AsType <- AS_type
  df <- df[,c(
    "AsType", "GeneID", "Key", 
    "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", 
    "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference"
  )]
  return(df)
}

all_rmats_res_df <- data.frame()
for(AS_type in AS_type_li){
  fname <- file.path(args$input, sprintf("%s.MATS.JC.txt", AS_type))
  all_rmats_res_df <- rbind(all_rmats_res_df, load_data(fname, AS_type))
}


all_rmats_res_df$AllSample1IJC <- sapply(strsplit(all_rmats_res_df$IJC_SAMPLE_1, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
all_rmats_res_df$AllSample1SJC <- sapply(strsplit(all_rmats_res_df$SJC_SAMPLE_1, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
all_rmats_res_df$AllSample2IJC <- sapply(strsplit(all_rmats_res_df$IJC_SAMPLE_2, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
all_rmats_res_df$AllSample2SJC <- sapply(strsplit(all_rmats_res_df$SJC_SAMPLE_2, ","), function(x){return(sum(as.numeric(x), na.rm = T))})
all_rmats_res_df$AllSample1 <- all_rmats_res_df$AllSample1IJC + all_rmats_res_df$AllSample1SJC
all_rmats_res_df$AllSample2 <- all_rmats_res_df$AllSample2IJC + all_rmats_res_df$AllSample2SJC

all_rmats_res_df$AveIncLevel1 <- sapply(strsplit(all_rmats_res_df$IncLevel1, ","), function(x){return(mean(as.numeric(x), na.rm = T))})
all_rmats_res_df$AveIncLevel2 <- sapply(strsplit(all_rmats_res_df$IncLevel2, ","), function(x){return(mean(as.numeric(x), na.rm = T))})

all_rmats_res_df$DE <- "NC"
all_rmats_res_df$DE[
  all_rmats_res_df$IncLevelDifference>args$delta_phi_cutoff & all_rmats_res_df$FDR < args$FDR_cutoff &
    all_rmats_res_df$AllSample1>args$read_num_cutoff & all_rmats_res_df$AllSample2>args$read_num_cutoff
  ] <- "Up"
all_rmats_res_df$DE[
  (-1*all_rmats_res_df$IncLevelDifference)>args$delta_phi_cutoff & all_rmats_res_df$FDR < args$FDR_cutoff &
    all_rmats_res_df$AllSample1>args$read_num_cutoff & all_rmats_res_df$AllSample2>args$read_num_cutoff
  ] <- "Down"
write_tsv(all_rmats_res_df, file.path(args$output, "AS.rMATS.SourceData.tsv"))
DE_rmats_res_df <- all_rmats_res_df[all_rmats_res_df$DE != "NC",]
DE_AS_num_stat <- DE_rmats_res_df %>% group_by(AsType, DE) %>% summarise(DE_AS_num=n())
write_tsv(DE_AS_num_stat, file.path(args$output, "AS.DE.AsNum.tsv"))

AsCircadianPsiCpmSourceData <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/AS_circadian_gene/AsCircadianPsiCpmSourceData.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
all_rmats_res_df$IsAsCircadian <- all_rmats_res_df$Key %in% AsCircadianPsiCpmSourceData$Key
table(all_rmats_res_df[,c("DE", "IsAsCircadian")])
