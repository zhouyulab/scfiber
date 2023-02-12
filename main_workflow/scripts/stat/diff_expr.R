library(readr)
library(edgeR)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(SingleCellExperiment)

args = commandArgs(TRUE)
if (length(args) == 8) {
  f_mergeSCE <- args[1]
  f_corSCE <- args[2]
  f_cluster <- args[3]
  f_tpm <- args[4]
  FC_cutoff <- as.numeric(args[5])
  fdr_cutoff <- as.numeric(args[6])
  cell_num_cutoff <- as.integer(args[7])
  f_out <- args[8]
} else {
  q()
}

# f_mergeSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/mergeSCE.RData"
# f_corSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/corSCE.RData"
# f_cluster <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/cluster.tsv"
# f_tpm <- "/home/sqreb/ias/data/CottonSingleCell/analysis/cluster_cnt/data/cluster.rep.tpm.tsv"
# FC_cutoff <- 2
# fdr_cutoff <- 0.05
# cell_num_cutoff <- 100
# f_out <- "~/test"


load(f_mergeSCE)
load(f_corSCE)

clusterSCE <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)
mergeSCE$cluster <- clusterSCE$Cluster
corSCE$cluster <- clusterSCE$Cluster

tpm <- read_delim(f_tpm, "\t", escape_double = FALSE, trim_ws = TRUE)

uniq_sample <- sort(unique(mergeSCE$SampleName), decreasing=T)
uniq_cluster <- c("C1", "C2", "C4", "C5")
if(length(uniq_sample)<2) stop()

downsample_cell <- function(cnt, select_indx){
  dsm <- cnt[,select_indx]
  g_cnt <- rowSums(dsm)
  return(g_cnt)
}

build_downsample_mat <- function(SCE, sample_name, cluster){
  sample_data <- SCE[, (mergeSCE$SampleName==sample_name) & (mergeSCE$cluster==cluster)]
  cnt <- counts(sample_data)
  reps <- sort(unique(SCE$Rep))
  select_index_group <- lapply(reps, function(x){return(which(sample_data$Rep==x))})
  dsm <- sapply(select_index_group, function(x){return(downsample_cell(cnt, x))})
  colnames(dsm) <- sprintf("%s.%s.%s", sample_name, cluster, reps)
  return(dsm)
}

df_res <- data.frame()
all_expr_res <- data.frame()
reps <- sort(unique(mergeSCE$Rep))

for(sample_i in 1:(length(uniq_sample)-1)){
  for(sample_j in (sample_i+1):length(uniq_sample)){
    for(clu in uniq_cluster){
      sample1 <- uniq_sample[sample_i]
      sample2 <- uniq_sample[sample_j]
      sample1_num <- sum((mergeSCE$SampleName==sample1)&(mergeSCE$cluster==clu))
      sample2_num <- sum((mergeSCE$SampleName==sample2)&(mergeSCE$cluster==clu))
      if(min(sample1_num, sample2_num) < cell_num_cutoff) next()
      sample1_mat <- build_downsample_mat(mergeSCE, sample1, clu)
      sample2_mat <- build_downsample_mat(mergeSCE, sample2, clu)
      combine_mat <- cbind(sample1_mat, sample2_mat)
      
      group <- factor(c(rep(1, length(reps)),rep(2, length(reps))))
      obj <- DGEList(counts=combine_mat,group=group)
      keep <- filterByExpr(obj)
      obj <- obj[keep, , keep.lib.sizes=FALSE]
      obj <- calcNormFactors(obj)
      design <- model.matrix(~group)
      obj <- estimateDisp(obj,design)
      fit <- glmQLFit(obj,design)
      qlf <- glmQLFTest(fit,coef=2)
      res <- as.data.frame(topTags(qlf, n=nrow(obj$counts)))
      
      res$Sample1 <- sample1
      res$Sample1Num <- sample1_num
      res$Sample2 <- sample2
      res$Sample2Num <- sample2_num
      res$Cluster <- clu
      res$Reg <- "NC"
      res$Reg[(res$logFC>log2(FC_cutoff)) & (res$FDR<fdr_cutoff)] <- "Up"
      res$Reg[(-res$logFC>log2(FC_cutoff)) & (res$FDR<fdr_cutoff)] <- "Down"
      res$Gid <- rownames(res)
      res <- res[,c("Gid", "Sample1", "Sample2", "Cluster", "Sample1Num", "Sample2Num", "Reg", "logFC", "logCPM", "F", "PValue", "FDR")]
      res$Reg <- factor(res$Reg, levels = c("Up", "NC", "Down"), labels =  c("Up", "NC", "Down"))
      res <- as.data.frame(res)
      all_expr_res <- rbind(all_expr_res, res)
      
      res$ppval <- -log10(res$PValue+1e-3)
      res$pfdr <- -log10(res$FDR+1e-3)
      col_li <- c("#cd2631", "#e9e9e9", "#4b79ae")
      names(col_li) <- c("Up", "NC", "Down")
      
      max_lfc <- max(abs(res$logFC), na.rm = TRUE)
      
      p <- ggplot(res, aes(x=logFC, y=pfdr, color=Reg)) +
        geom_point(size=0.3) +
        labs(x="Log2 fold change", y="FDR", title = sprintf("%s (%d) vs. %s (%d)", sample2, sample2_num, sample1, sample1_num)) +
        scale_y_continuous(
          breaks = c(0, 1, 2, 3), 
          labels = c(
            expression("10"^"-0"), expression("10"^"-1"), expression("10"^"-2"), expression("10"^"-3")
          )) +
        scale_color_manual(values = col_li) +
        coord_cartesian(xlim = c(-max_lfc, max_lfc), ylim = c(0, 3)) +
        theme_bw() +
        theme(text = element_text(family="ArialMT", color = "black", size = 6),
              axis.line = element_line(color = "black"),
              axis.text = element_text(color = "black"),
              axis.ticks = element_line(color = "black"),
              legend.background = element_blank(),
              legend.title = element_blank(),
              legend.key.size = unit(4, "mm"),
              panel.grid = element_blank())
      
      diff_expr_info <- as.data.frame(res %>% group_by(Reg) %>% summarise(Num=n()))
      if(any(res$Reg!="NC")){
        diff_expr_info$Text <- sprintf("N=%d", diff_expr_info$Num)
        diff_expr_info <- diff_expr_info[which(diff_expr_info$Reg!="NC"),]
        diff_expr_info$x <- -max_lfc*0.8
        diff_expr_info$x[diff_expr_info$Reg=="Up"] <- -diff_expr_info$x[diff_expr_info$Reg=="Up"]
        p <- p + geom_text(aes(x=x, y=2.8, label=Text), data=diff_expr_info, size=1.5, color="black")
      }
      
      tmp_dir <- file.path(f_out, "edgeR", sprintf("%s_vs_%s", sample2, sample1))
      if(! dir.exists(tmp_dir)) dir.create(tmp_dir, recursive=TRUE)
      ggsave(file.path(tmp_dir, sprintf("DiffExpr.%s.pdf", clu)), p, height = 9, width = 12, units = "cm", dpi = 600)
    }
  }
}
diff_expr_df <- all_expr_res[which(all_expr_res$Reg!="NC"),]
write_tsv(left_join(diff_expr_df, tpm, by="Gid"), file.path(f_out, "DiffExpr.tsv"))
write_tsv(left_join(all_expr_res, tpm, by="Gid"), file.path(f_out, "DiffExpr.data.tsv"))
