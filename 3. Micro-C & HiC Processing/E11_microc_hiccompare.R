###
### HiCcompare loess normalization on BEDPE matrices
###

library(multiHiCcompare)
library(HiCcompare)
library(BiocParallel)

register(MulticoreParam(workers = 16), default = TRUE)
options(scipen = 999)

setwd("./")

reslist <- c(
  1000000,
  500000,
  250000,
  100000,
  50000,
  25000,
  20000,
  10000,
  5000,
  2500,
  1000,
  800,
  600,
  400
)
chromlist <- c(1:19, "X")

mm39_blacklist <- readRDS("mm39_blacklist.RDS")
for (i in reslist) {
  for (chrom in chromlist) {
    clamped <- read.table(paste("bedpe/clamped_", i, "_", chrom, "_raw.txt", sep = ""))
    unclamped <- read.table(paste("bedpe/unclamped_", i, "_", chrom, "raw.txt", sep = ""))

    joint_table <- create.hic.table(
      clamped,
      unclamped,
      chr = paste("chr", chrom, sep = ""),
      exclude.regions = mm39_blacklist,
      exclude.overlap = 0.2,
      include.zeros = T
    )
    hic.table <-
      hic_loess(
        joint_table,
        Plot = F,
        Plot.smooth = FALSE,
        parallel = T
      )
    hic.table <-
      hic_compare(
        hic.table,
        A.min = NA,
        adjust.dist = TRUE,
        p.method = "fdr",
        Plot = F,
        parallel = T
      )

    write.table(
      hic.table,
      paste("hiccompare/hictable_", i, "_", chrom, ".txt", sep = ""),
      col.names = T,
      row.names = F,
      quote = F
    )
    write.table(
      hic.table[, c("start1", "start2", "adj.IF1")],
      paste(
        "hiccompare/clamped_hiccompare_",
        i,
        "_",
        chrom,
        ".txt",
        sep = ""
      ),
      col.names = F,
      row.names = F,
      quote = F
    )
    write.table(
      hic.table[, c("start1", "start2", "adj.IF2")],
      paste(
        "hiccompare/unclamped_hiccompare_",
        i,
        "_",
        chrom,
        ".txt",
        sep = ""
      ),
      col.names = F,
      row.names = F,
      quote = F
    )
  }
}
