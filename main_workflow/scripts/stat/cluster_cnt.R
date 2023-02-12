library(readr)
library(dplyr)
library(SingleCellExperiment)
args = commandArgs(TRUE)
if (length(args) == 5) {
  f_mergeSCE <- args[1]
  f_cluster <- args[2]
  f_sample_rep_cnt <- args[3]
  f_sample_cnt <- args[4]
  clu_col <- args[5]
} else {
  q()
}
# f_out <- "~/cluster.cnt.tsv"
# f_mergeSCE <- "~/mergeSCE.RData"

load(f_mergeSCE)
clusterSCE <- read_delim(f_cluster, "\t", escape_double = FALSE, trim_ws = TRUE)
mergeSCE$cluster <- unlist(clusterSCE[, clu_col])

batch_naems <- sort(unique(mergeSCE$Sample))
cluster_names <- sort(unique(mergeSCE$cluster))
df <- data.frame()
for(batch_name in batch_naems){
  for(cluster in cluster_names){
    tmpSCE <- mergeSCE[,(mergeSCE$Sample==batch_name) & (mergeSCE$cluster==cluster)]
    tmp_cnt <- normcounts(tmpSCE)
    gene_umi <- rowSums(tmp_cnt)
    tmp_df <- data.frame(Gid=names(gene_umi), UMI=gene_umi)
    names(tmp_df)[2] <- sprintf("%s.%s", batch_name, cluster)
    if(nrow(df)==0){
      df <- tmp_df
    }else{
      if(!all(tmp_df$Gid==df$Gid)) stop()
      df <- left_join(df, tmp_df, by="Gid")
    }
  }
}
write_tsv(df, f_sample_rep_cnt)

df <- data.frame()
sample_names <- unique(mergeSCE$SampleName)
for(sample_name in sample_names){
  for(cluster in cluster_names){
    tmpSCE <- mergeSCE[,(mergeSCE$SampleName==sample_name) & (mergeSCE$cluster==cluster)]
    tmp_cnt <- normcounts(tmpSCE)
    gene_umi <- rowSums(tmp_cnt)
    tmp_df <- data.frame(Gid=names(gene_umi), UMI=gene_umi)
    colName <- sprintf("%s.%s", sample_name, cluster)
    names(tmp_df)[2] <- colName
    if(nrow(df)==0){
      df <- tmp_df
    }else{
      if(!all(tmp_df$Gid==df$Gid)) stop()
      if(colName %in% names(df)){
        df[, colName] <- df[, colName] + tmp_df[, colName]
      }else{
        df <- left_join(df, tmp_df, by="Gid")
      }
    }
  }
}
write_tsv(df, f_sample_cnt)
