library(SingleCellExperiment)
library(readr)

args = commandArgs(TRUE)
if (length(args) == 2) {
  f_mergeSCE <- args[1]
  f_out <- args[2]
} else {
  q()
}
# f_mergeSCE <- "/home/sqreb/ias/data/CottonSingleCell/analysis/base_map/mergeSCE.RData"

load(f_mergeSCE)
cnt <- counts(mergeSCE)
norm_cnt <- normcounts(mergeSCE)
size_factor <- colSums(norm_cnt) / colSums(cnt)
df <- data.frame(Sample=mergeSCE$SampleName, Rep=mergeSCE$Rep, Barcode=mergeSCE$Barcode, SizeFactor=size_factor)
write_tsv(df, f_out)