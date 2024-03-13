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
  mat <- read.table(paste("clamped.", format(i, scientific = F), ".matrix", sep = ""))
  bed <- read.table(paste("clamped.", format(i, scientific = F), "_abs.bed", sep = ""))
  clamped <- hicpro2bedpe(mat, bed)
  clamped <- clamped$cis

  for (chrom in c(1:19, "X")) {
    cat("Processing chrom ", chrom, "\n", sep = "")
    eval(parse(text = (
      paste("clamped_chrom = clamped$`", chrom, "`", sep = "")
    )))
    write.table(
      clamped_chrom[, c(2, 5, 7)],
      paste("bedpe/clamped_", i, "_", chrom, "_raw.txt", sep = ""),
      col.names = F,
      row.names = F,
    )
  }
}

for (i in reslist) {
  cat("Processing res ", i, "\n", sep = "")
  mat <- read.table(paste("unclamped.", format(i, scientific = F), ".matrix", sep = ""))
  bed <- read.table(paste("unclamped.", format(i, scientific = F), "_abs.bed", sep = ""))
  unclamped <- hicpro2bedpe(mat, bed)
  unclamped <- unclamped$cis

  for (chrom in c(1:19, "X")) {
    cat("Processing chrom ", chrom, "\n", sep = "")
    eval(parse(text = (
      paste("unclamped_chrom = unclamped$`", chrom, "`", sep = "")
    )))
    write.table(
      unclamped_chrom[, c(2, 5, 7)],
      paste("bedpe/unclamped_", i, "_", chrom, "raw.txt", sep = ""),
      col.names = F,
      row.names = F,
    )
  }
}
