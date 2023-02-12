library(readr)

args = commandArgs(TRUE)
if (length(args) == 2) {
  f_in <- args[1]
  f_out <- args[2]
} else {
  q()
}

cnt <- read_delim(f_in, "\t", escape_double = FALSE, trim_ws = TRUE)
val <- cnt[,2:ncol(cnt)]
cnt_colSum <- colSums(val)
val[,cnt_colSum < 1e5] <- NA
cnt[,2:ncol(cnt)] <- as.data.frame(t(1e6 * t(val) / cnt_colSum))
write_tsv(cnt, f_out)
