###
### TADcompare boundary calling
###

library(TADCompare)
library(HiCcompare)
library(BiocParallel)

options(scipen = 999)
setwd("./")
dir <- getwd()

numCores <- 14
register(MulticoreParam(workers = numCores), default = TRUE)

reslist <- c(500000, 250000, 100000, 50000, 25000, 10000)
for (i in reslist) {
  for (chrom in c(1:19, "X")) {
    distal <- read.table(paste("distal_hiccompare_", i, "_", chrom, ".txt", sep = ""))
    proximal <- read.table(paste("proximal_hiccompare_", i, "_", chrom, ".txt", sep = ""))

    results <- TADCompare(proximal, distal, resolution = i)
    write.table(
      results$TAD_Frame,
      paste(
        "tadcompare/tadcompare_e12_",
        format(i, scientific = F),
        "_",
        chrom,
        "tadframe",
        sep = ""
      )
    )
    write.table(
      results$Boundary_Scores,
      paste(
        "tadcompare/tadcompare_e12_",
        format(i, scientific = F),
        "_",
        chrom,
        "boundaryscores",
        sep = ""
      )
    )
    saveRDS(
      results,
      paste(
        "tadcompare/tadcompare_e12_",
        format(i, scientific = F),
        "_",
        chrom,
        ".rds",
        sep = ""
      )
    )
  }
}
