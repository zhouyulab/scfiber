motif_hit_df <- read_delim("analysis_v3/ATAC/motif/summary/MotifExprTargetGene.HitMat.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
cluster_tpm <- read_delim("analysis_v3/stat/cluster_cnt/cluster.tpm.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RNA_seq_TP <- read_delim("~/mu01/project/CottonSCE_TP/analysis/RNA_seq_TP/circadian_gene/CPM.merge.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
blast2go <- read_delim("analysis_v3/annotation/blast2GO/blast2go.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
blast2go <- blast2go[,c("SeqName", "Description")]
names(blast2go) <- c("Tid", "BLAST2GO")
blast2go$Gid <- sapply(strsplit(blast2go$Tid, "[.]"), function(x){return(x[1])})
blast2go <- blast2go[!duplicated(blast2go$Gid),]
blast2go <- blast2go[,c("Gid", "BLAST2GO")]
cotton2tair <- read_delim("data/genome/HAU/Gh_PrimaryPeptide_vs_Arabidopsis_Annotation.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 1)
cotton2tair <- cotton2tair[,c(1,4)]
names(cotton2tair) <- c("Tid", "TairID")
cotton2tair$Gid <- sapply(strsplit(cotton2tair$Tid, "[.]"), function(x){return(x[1])})
cotton2tair <- cotton2tair[!duplicated(cotton2tair$Gid),]
cotton2tair <- cotton2tair[,c("Gid", "TairID")]

Gbox <- data.frame(Gid=motif_hit_df$Gid[motif_hit_df$Gbox])

Gbox_CPM <- left_join(Gbox, CPM_merge)
Gbox_CPM_mat <- as.matrix(Gbox_CPM[,2:ncol(Gbox_CPM)])
Gbox_CPM_mat[apply(Gbox_CPM_mat, 1, max)<=1,] <- NA
rownames(Gbox_CPM_mat) <- Gbox_CPM$Gid
Gbox_cluster_tpm <- left_join(Gbox, cluster_tpm)
Gbox_cluster_tpm_mat <- as.matrix(Gbox_cluster_tpm[,2:ncol(Gbox_cluster_tpm)])
Gbox_cluster_tpm_mat[apply(Gbox_cluster_tpm_mat, 1, function(x){return(max(x, na.rm = T))})<=1,] <- NA
rownames(Gbox_cluster_tpm_mat) <- Gbox_cluster_tpm$Gid
NA_indx <- is.na(Gbox_CPM_mat[,1]) | is.na(Gbox_cluster_tpm_mat[,1])
Gbox <- data.frame(Gid=Gbox$Gid[!NA_indx])

Gbox <- left_join(Gbox, cluster_tpm)
Gbox <- left_join(Gbox, RNA_seq_TP)
Gbox <- left_join(Gbox, blast2go)
Gbox <- left_join(Gbox, cotton2tair)

write_tsv(Gbox, "analysis_v3/ATAC/motif/Gbox.info.tsv")

f_TCP_gene <- "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP/TCP.gene.tsv"
f_TCP_like_motif_gene <- "analysis_v3/ATAC/ATAC_kmer_cluster2gene/TCP_like_motif/TCP_like_motif.gene.tsv"
TCP_gene <- read_delim(f_TCP_gene, "\t", escape_double = FALSE, trim_ws = TRUE)
TCP_like_motif_gene <- read_delim(f_TCP_like_motif_gene, "\t", escape_double = FALSE, trim_ws = TRUE)

expr_tmp <- cluster_tpm[apply(cluster_tpm[,2:ncol(cluster_tpm)], 1, function(x){return(max(x, na.rm = TRUE))})>1,]
expr_tmp$IsC3Max <- apply(expr_tmp[2:ncol(expr_tmp)], 1, which.max)==3

TCP_df <- data.frame(Gid=unique(c(TCP_gene$Gid, TCP_like_motif_gene$Gid)))
TCP_df$IsTCP <- TCP_df$Gid %in% TCP_gene$Gid
TCP_df$IsTCP_like <- TCP_df$Gid %in% TCP_like_motif_gene$Gid
TCP_df <- TCP_df[TCP_df$Gid%in%expr_tmp$Gid[expr_tmp$IsC3Max],]

TCP_df <- left_join(TCP_df, cluster_tpm)
TCP_df <- left_join(TCP_df, RNA_seq_TP)
TCP_df <- left_join(TCP_df, blast2go)
TCP_df <- left_join(TCP_df, cotton2tair)
write_tsv(TCP_df, "analysis_v3/ATAC/motif/TCP.info.tsv")
