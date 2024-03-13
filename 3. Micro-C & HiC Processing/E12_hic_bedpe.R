###
### Conversion of raw HiCPro matrix output to BEDPE format
###

library(HiCcompare)
setwd("./")
options(scipen = 999)

reslist <- c(
  400,
  600,
  800,
  1000,
  2500,
  5000,
  10000,
  20000,
  25000,
  40000,
  50000,
  100000,
  250000,
  500000,
  1000000
)

for (i in reslist) {
  cat("Processing res ", i, "\n", sep = "")
  mat <- read.table(paste("distal.", format(i, scientific = F), ".matrix", sep = ""))
  bed <- read.table(paste("distal.", format(i, scientific = F), "_abs.bed", sep = ""))
  distal <- hicpro2bedpe(mat, bed)
  distal <- distal$cis

  for (chrom in c(1:19, "X")) {
    cat("Processing chrom ", chrom, "\n", sep = "")
    eval(parse(text = (
      paste("distal_chrom = distal$`", chrom, "`", sep = "")
    )))
    write.table(
      distal_chrom[, c(2, 5, 7)],
      paste("bedpe/distal_", i, "_", chrom, "_raw.txt", sep = ""),
      col.names = F,
      row.names = F,
    )
  }
}

for (i in reslist) {
  cat("Processing res ", i, "\n", sep = "")
  mat <- read.table(paste("proximal.", format(i, scientific = F), ".matrix", sep = ""))
  bed <- read.table(paste("proximal.", format(i, scientific = F), "_abs.bed", sep = ""))
  proximal <- hicpro2bedpe(mat, bed)
  proximal <- proximal$cis

  for (chrom in c(1:19, "X")) {
    cat("Processing chrom ", chrom, "\n", sep = "")
    eval(parse(text = (
      paste("proximal_chrom = proximal$`", chrom, "`", sep = "")
    )))
    write.table(
      proximal_chrom[, c(2, 5, 7)],
      paste("bedpe/proximal_", i, "_", chrom, "raw.txt", sep = ""),
      col.names = F,
      row.names = F,
    )
  }
}
