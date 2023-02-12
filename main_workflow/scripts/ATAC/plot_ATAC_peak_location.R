library(readr)
library(reshape2)
library(ggplot2)

peak_location <- read_delim("analysis_v3/ATAC/peak_profile/peak.location.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
peak_location$AtacPeakRatio <- peak_location$PeakBaseNum / peak_location$AllBaseNum
peak_location <- peak_location[order(peak_location$AtacPeakRatio, decreasing = T),]
df <- data.frame(
  Stat=c(rep("#Peak base", nrow(peak_location)), rep("Peak base / Genome base", nrow(peak_location))),
  Location=rep(peak_location$Location, 2),
  Value=c(peak_location$PeakBaseNum, peak_location$AtacPeakRatio)
  )

df$Location <- factor(df$Location, levels = rev(peak_location$Location))

p <- ggplot(df, aes(x=Location, y=Value)) +
  geom_bar(fill="black", stat="identity") +
  facet_grid(~Stat, scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(
    text = element_text(family="ArialMT", color = "black", size = 6),
    axis.line = element_line(color="black", size=0.3),
    axis.title = element_blank(),
    panel.grid = element_blank(),
  )
ggsave("analysis_v3/ATAC/peak_profile/peak.location.pdf", p, height = 4.5, width = 9, units = "cm", colormodel = "cmyk")

