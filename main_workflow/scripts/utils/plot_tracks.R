library(argparse)
library(readr)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(Sushi)
library(BiocParallel)
parser_opts <- function(){
  parser <- ArgumentParser()
  parser$add_argument("--bw", nargs="*", required=FALSE, type="character", dest = "bw", metavar="data.bw")
  parser$add_argument("--bw-label", nargs="*", required=FALSE, type="character", dest = "bw_label", metavar="bw_label")
  parser$add_argument("--bw-color", nargs="*", required=FALSE, type="character", dest = "bw_color", metavar="bw_color")
  parser$add_argument("--bg", nargs="*", required=FALSE, type="character", dest = "bg", metavar="data.bg")
  parser$add_argument("--bg-label", nargs="*", required=FALSE, type="character", dest = "bg_label", metavar="bg_label")
  parser$add_argument("--bg-color", nargs="*", required=FALSE, type="character", dest = "bg_color", metavar="bg_color")
  parser$add_argument("--ref-bed12", nargs="*", required=FALSE, type="character", dest = "ref_bed12", metavar="ref.bed12")
  parser$add_argument("--ref-bed12-label", nargs="*", required=FALSE, type="character", dest = "ref_bed12_label", metavar="ref_bed12_label")
  parser$add_argument("--ref-bed12-color", nargs="*", required=FALSE, type="character", dest = "ref_bed12_color", metavar="ref_bed12_color")
  parser$add_argument("--ref-bed6", nargs="*", required=FALSE, type="character", dest = "ref_bed6", metavar="ref.bed6")
  parser$add_argument("--ref-bed6-label", nargs="*", required=FALSE, type="character", dest = "ref_bed6_label", metavar="ref_bed6_label")
  parser$add_argument("--ref-bed6-color", nargs="*", required=FALSE, type="character", dest = "ref_bed6_color", metavar="ref_bed6_color")
  parser$add_argument("--loc-bed", required=FALSE, type="character", dest = "loc_bed", metavar="loc_bed")
  parser$add_argument("--expand-len", required=TRUE, type="integer", dest = "expand_len", metavar="expand_len")
  parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")
  return(parser)
}

args <- commandArgs(TRUE)
parser <- parser_opts()
# args <- c(
#   "--bw", "analysis_v3/ATAC/bw/WT.bw", "analysis_v3/ATAC/bw/FL.bw",
#   "--bw-label", "WT ATAC", "FL ATAC", "--bw-color", "red", "blue",
#   "--ref-bed12", "data/genome/HAU/Ghirsutumv1.1_gene_model.bed12", 
#   "--ref-bed12-label", "Gene", "--ref-bed12-color", "black",
#   "--ref-bed6", "analysis_v3/ATAC/data_MACS2/ATAC.All.bed6",
#   "--ref-bed6-label", "ATAC peak", "--ref-bed6-color", "black",
#   "--loc-bed", "data/genome/HAU/Ghirsutumv1.1_gene_model.gene.bed6",
#   "--expand-len", "1000", "-o", "analysis_v3/ATAC/track"
# )
args <- parser$parse_args(args)

check_args <- function(f, label, col){
  if (is.null(f)) return(NULL)
  if (length(f) != length(label)) stop()
  if (is.null(col)){
    col <- rep("black", length(f))
  }
  if (length(col)==1){
    col <- rep(col, length(f))
  }
  if (length(f) != length(col)) stop()
  cfg <- list(file=f, label=label, color=col)
  return(cfg)
}

plot_bg <- function(cfg, gr){
  sample_len <- length(cfg$file)
  for(indx in 1:sample_len){
    try({
      f <- cfg$file[indx]
      label <- cfg$label[indx]
      col <- cfg$color[indx]
      bg <- import.bedGraph(f)
      bg <- subsetByOverlaps(bg, gr)
      bg_df <- as.data.frame(bg)
      bg <- data.frame(chrom=bg_df$seqnames, start=bg_df$start, end=bg_df$end, value=abs(bg_df$score))
      chrom <- as.character(gr@seqnames@values)
      start <- gr@ranges@start
      end <- start + gr@ranges@width
      plotBedgraph(bg, chrom, start, end, color = col)
      labelgenome(chrom,start,end,n=3,scale="Kb")
      axis(side=2,las=2,tcl=.2)
      mtext(label,side=2,line=1.75,cex=1,font=2)
    })
  }
}

plot_bw <- function(cfg, gr){
  sample_len <- length(cfg$file)
  for(indx in 1:sample_len){
    try({
      f <- cfg$file[indx]
      label <- cfg$label[indx]
      col <- cfg$color[indx]
      bw <- import.bw(f, selection=BigWigSelection(gr))
      bw_df <- as.data.frame(bw)
      bw <- data.frame(chrom=bw_df$seqnames, start=bw_df$start, end=bw_df$end, value=abs(bw_df$score))
      chrom <- as.character(gr@seqnames@values)
      start <- gr@ranges@start
      end <- start + gr@ranges@width
      plotBedgraph(bw, chrom, start, end, color = col)
      labelgenome(chrom,start,end,n=3,scale="Kb")
      axis(side=2,las=2,tcl=.2)
      mtext(label,side=2,line=1.75,cex=1,font=2)
    })
  }
}

plot_bed6 <- function(cfg, gr){
  sample_len <- length(cfg$file)
  for(indx in 1:sample_len){
    try({
      f <- cfg$file[indx]
      label <- cfg$label[indx]
      col <- cfg$color[indx]
      bed <- import.bed(f)
      bed <- subsetByOverlaps(bed, gr)
      bed_df <- as.data.frame(bed)
      bed <- data.frame(chrom=bed_df$seqnames, start=bed_df$start, end=bed_df$end, name=bed_df$name, score=bed_df$score, strand=bed_df$strand)
      chrom <- as.character(gr@seqnames@values)
      start <- gr@ranges@start
      end <- start + gr@ranges@width
      plotBed(bed, chrom, start, end, color = col, numbins=500)
      labelgenome(chrom,start,end,n=3,scale="Kb")
      mtext(label,side=2,line=1.75,cex=1,font=2)
    })
  }
}

plot_bed12 <- function(cfg, gr){
  sample_len <- length(cfg$file)
  for(indx in 1:sample_len){
    f <- cfg$file[indx]
    label <- cfg$label[indx]
    col <- cfg$color[indx]
    bed <- import.bed(f)
    bed <- subsetByOverlaps(bed, gr)
    bed_df <- as.data.frame(bed)
    
    data_df <- data.frame()
    for(g_idx in 1:length(bed)){
      g_info <- bed[g_idx]
      exons <- blocks(g_info)[[1]]
      exon_df <- as.data.frame(exons)
      exon_df$strand_int <- 1
      exon_df$strand_int[exon_df$strand=="-"] <- -1
      tmp_df <- data.frame(chrom=exon_df$seqnames, start=exon_df$start, stop=exon_df$end, gene=g_info$name, score=g_info$score, strand=exon_df$strand_int, type="exon")
      data_df <- rbind(data_df, tmp_df)
    }
    chrom <- as.character(gr@seqnames@values)
    start <- gr@ranges@start
    end <- start + gr@ranges@width
    plotGenes(data_df, chrom, start, end, types = data_df$type, col = col)
    labelgenome(chrom,start,end,n=3,scale="Kb")
    mtext(label,side=2,line=1.75,cex=1,font=2)
  }
}

plot_track <- function(cfg, bg_cfg, bw_cfg, bed12_cfg, bed6_cfg){
  fname <- cfg[["fname"]]
  gr <- cfg[["gr"]]
  fig_num <- 0
  if(! is.null(bg_cfg)){
    fig_num <- fig_num + length(bg_cfg$file)
  }
  if(! is.null(bw_cfg)){
    fig_num <- fig_num + length(bw_cfg$file)
  }
  if(! is.null(bed12_cfg)){
    fig_num <- fig_num + length(bed12_cfg$file)
  }
  if(! is.null(bed6_cfg)){
    fig_num <- fig_num + length(bed6_cfg$file)
  }
  
  inche_cm=2.54
  pdf(fname, width=20/inche_cm, height=6*fig_num/inche_cm, family="ArialMT", colormodel = "cmyk")
  par(mfrow = c(fig_num, 1))
  if(! is.null(bg_cfg)){
    plot_bg(bg_cfg, gr)
  }
  if(! is.null(bw_cfg)){
    plot_bw(bw_cfg, gr)
  }
  if(! is.null(bed6_cfg)){
    plot_bed6(bed6_cfg, gr)
  }
  if(! is.null(bed12_cfg)){
    plot_bed12(bed12_cfg, gr)
  }
  dev.off()
}

bg_cfg <- check_args(args$bg, args$bg_label, args$bg_color)
bw_cfg <- check_args(args$bw, args$bw_label, args$bw_color)
bed12_cfg <- check_args(args$ref_bed12, args$ref_bed12_label, args$ref_bed12_color)
bed6_cfg <- check_args(args$ref_bed6, args$ref_bed6_label, args$ref_bed6_color)

loc_bed <- read_delim(args$loc_bed, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
loc_bed <- loc_bed[,c(1, 2, 3, 4)]
names(loc_bed) <- c("Chrom", "Start", "End", "Name")
loc_bed$ChromAbbr <- loc_bed$Chrom
loc_bed$ChromAbbr[grep("Scaffold", loc_bed$ChromAbbr)] <- "Scaffold"
for(chrom in unique(loc_bed$ChromAbbr)){
  print(chrom)
  tmp_dir <- file.path(args$output, chrom)
  if(!dir.exists(tmp_dir)){
    dir.create(tmp_dir)
  }
  tmp_loc <- loc_bed[loc_bed$ChromAbbr==chrom,]
  cfg_li <- list()
  for(indx in 1:nrow(tmp_loc)){
    fname <- file.path(tmp_dir, sprintf("%s.pdf", tmp_loc$Name[indx]))
    gr <- GRanges(seqnames=tmp_loc$Chrom[indx], ranges=IRanges(max(0, tmp_loc$Start[indx]-args$expand_len), max(0, tmp_loc$End[indx]+args$expand_len)))
    cfg <- list(fname=fname, gr=gr)
    cfg_li[[tmp_loc$Name[indx]]] <- cfg
    # plot_track(cfg, bg_cfg, bw_cfg, bed12_cfg, bed6_cfg)
  }
  bplapply(cfg_li, plot_track, bg_cfg=bg_cfg, bw_cfg=bw_cfg, bed12_cfg=bed12_cfg, bed6_cfg=bed6_cfg, BPPARAM=MulticoreParam(28))
}


