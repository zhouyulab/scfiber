library(readr)
library(dplyr)
library(argparse)
library(ggplot2)
library(gridExtra)
library(extrafont)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input.tsv")
parser$add_argument("-t", "--tag", nargs="*", required=TRUE, type="character", dest = "tag", metavar="tag")
parser$add_argument("--group-label", nargs="*", required=TRUE, type="character", dest = "group_label", metavar="group_label")
parser$add_argument("--color-label", nargs="*", required=TRUE, type="character", dest = "color_label", metavar="color_label")
parser$add_argument("--out-plot", required=TRUE, type="character", dest = "out_plot", metavar="out_plot.pdf")
parser$add_argument("--out-data", required=TRUE, type="character", dest = "out_data", metavar="out_data.tsv")
parser$add_argument("--title", required=TRUE, type="character", dest = "title", metavar="title")
parser$add_argument("--type", required=TRUE, type="character", dest = "type", metavar="type", choices=c("point", "region"))
parser$add_argument("--filter-cnt", action="store_true", dest = "filter_cnt")
args <- commandArgs(TRUE)
args <- parser$parse_args(args)

if(length(args$input) != length(args$tag)) stop()
if(length(args$input) != length(args$group_label)) stop()
if(length(args$input) != length(args$color_label)) stop()

plot_point <- function(f_li, tag_li, group_li, color_li, title, plot_xlab=TRUE, do.filter=FALSE){
  expr_df <- data.frame()
  for(file_indx in 1:length(f_li)){
    df <- read_delim(f_li[file_indx], "\t", escape_double = FALSE, trim_ws = TRUE)
    expr_mat <- as.matrix(df[,2:ncol(df)])
    if(do.filter){
      mid_val <- apply(expr_mat, 2, median)
      mad_val <- apply(expr_mat, 2, mad)
      cutoff_val <- mid_val + 3 * mad_val
      for(indx in 1:ncol(expr_mat)){
        expr_mat[,indx][expr_mat[,indx]>cutoff_val[indx]] <- cutoff_val[indx]
      }
    }
    expr_li <- colSums(expr_mat)
    expr_li <- 100 * expr_li / sum(expr_li)
    tmp_expr_df <- data.frame(
      Position=as.numeric(names(df)[2:ncol(df)]),
      Expr=expr_li,
      Label=tag_li[file_indx],
      Group=group_li[file_indx],
      Color=color_li[file_indx]
      )
    tmp_expr_df <- tmp_expr_df[2:(nrow(tmp_expr_df)-1),]
    expr_df <- rbind(expr_df, tmp_expr_df)
  }
  
  p <- ggplot(expr_df, aes(x=Position, y=Expr, color=Color)) +
    geom_line(size=0.1) +
    facet_grid(~Group) +
    scale_x_continuous(breaks = c(min(expr_df$Position-1), 0, max(expr_df$Position+1)), labels = c(as.character(min(expr_df$Position)-1), title, as.character(max(expr_df$Position)+1))) +
    theme_bw() +
    theme(
      text = element_text(size=6, family = "ArialMT"),
      axis.text = element_text(color="black"),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size=6),
      legend.key.size = unit(3, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(1,1,1,1,"mm"),
      panel.spacing = unit(0, "mm"),
      panel.background = element_blank()
    )
  if(! plot_xlab){
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  res <- list(data=expr_df, plot=p)
  return(res)
}

plot_scale <- function(f_li, tag_li, group_li, color_li, title, plot_xlab=TRUE, do.filter=FALSE){
  expr_df <- data.frame()
  for(file_indx in 1:length(f_li)){
    df <- read_delim(f_li[file_indx], "\t", escape_double = FALSE, trim_ws = TRUE)
    expr_mat <- as.matrix(df[,2:ncol(df)])
    if(do.filter){
      mid_val <- apply(expr_mat, 2, median)
      mad_val <- apply(expr_mat, 2, mad)
      cutoff_val <- mid_val + 3 * mad_val
      for(indx in 1:ncol(expr_mat)){
        expr_mat[,indx][expr_mat[,indx]>cutoff_val[indx]] <- cutoff_val[indx]
      }
    }
    expr_li <- colSums(expr_mat)
    expr_li <- 100 * expr_li / sum(expr_li)
    tmp_expr_df <- data.frame(
      Location=names(df)[2:ncol(df)],
      Expr=expr_li,
      Label=tag_li[file_indx],
      Group=group_li[file_indx],
      Color=color_li[file_indx]
    )
    tmp_expr_df$AxisLabel <- ""
    tmp_expr_df$AxisLabel[tmp_expr_df$Location=="0.00"] <- "Start site"
    tmp_expr_df$AxisLabel[tmp_expr_df$Location=="0"] <- "End site"
    
    start_site_indx <- which(tmp_expr_df$AxisLabel == "Start site")
    end_site_indx <- which(tmp_expr_df$AxisLabel == "End site")
    
    tmp_expr_df <- tmp_expr_df[(start_site_indx - 100): (end_site_indx + 100),]
    tmp_expr_df$AxisLabel[1] <- as.character(tmp_expr_df$Location[1])
    tmp_expr_df$AxisLabel[nrow(tmp_expr_df)] <- as.character(tmp_expr_df$Location[nrow(tmp_expr_df)])
    tmp_expr_df$x <- 1:nrow(tmp_expr_df)
    label_indx <- which(tmp_expr_df$AxisLabel != "")
    
    expr_df <- rbind(expr_df, tmp_expr_df)
  }
  
  
  p <- ggplot(expr_df, aes(x=x, y=Expr, color=Color)) +
    geom_line(size=0.1) +
    scale_x_continuous(breaks = label_indx, labels = expr_df$AxisLabel[label_indx]) +
    facet_grid(~Group) +
    theme_bw() +
    theme(
      text = element_text(size=6, family = "ArialMT"),
      axis.title.x = element_blank(),
      axis.text = element_text(color="black"),
      legend.title = element_blank(),
      legend.text = element_text(size=6),
      legend.key.size = unit(3, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(1,1,1,1,"mm"),
      panel.spacing = unit(0, "mm"),
      panel.background = element_blank()
    )
  if(! plot_xlab){
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  res <- list(data=expr_df, plot=p)
  return(res)
}

if(args$type == "point"){
  res <- plot_point(args$input, args$tag, args$group_label, args$color_label, args$title, plot_xlab=TRUE, do.filter=args$filter_cnt)
}else{
  res <- plot_scale(args$input, args$tag, args$group_label, args$color_label, args$title, plot_xlab=TRUE, do.filter=args$filter_cnt)
}
write_tsv(res$data, args$out_data)
ggsave(args$out_plot, res$plot, width = 3+5*length(unique(args$group_label)), height = 1+4, units = "cm")