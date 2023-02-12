library(dplyr)
library(readr)

args = commandArgs(TRUE)
if (length(args) == 3) {
  f_map <- args[1]
  f_tfdb <- args[2]
  f_out <- args[3]
} else {
  q()
}
# f_map <- "~/ias/data/CottonSingleCell/analysis/colinear/prot_pair_blast/HUA2NUA_NBI.tsv"
# f_tfdb <- "~/ias/common/PlantTFDB/tf_list/gosHir/Ghi_TF_list"
gene_map <- read_delim(f_map, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(gene_map) <- c("Gid", "NUA_NBI")
tfdb <- read_delim(f_tfdb, "\t", escape_double = FALSE, trim_ws = TRUE)
tfdb$TF_ID <- NULL
names(tfdb)[1] <- "NUA_NBI"
res <- inner_join(gene_map, tfdb)
write_tsv(res, f_out)
