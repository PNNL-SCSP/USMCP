library(dplyr)
d1 <- read.csv("MOp_Raw_intensities.csv", header = T, stringsAsFactors = F)

d1 <- d1[!duplicated(d1$PG.Genes), ]

write.csv(d1, "MOp_Raw_Quan_values.csv", row.names = F)

# this code block is discarded.
# reason is need to keep all 3211 protein to reproduce the clustering plot.