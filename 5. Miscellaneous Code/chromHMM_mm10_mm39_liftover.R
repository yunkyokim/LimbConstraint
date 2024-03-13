###
### Liftover code converting mm10 ChromHMM annotations to mm39
### https://doi.org/10.1038/s42003-021-01756-4


library(rtracklayer)
library(rGREAT)

options(scipen = 999)
setwd("./")
dir <- getwd()

chain_file <- "chromhmm/mm10ToMm39.over.chain"
chain_file <- import.chain(chain_file)

bed_file <- "chromhmm/E12_18state_chromhmm.bed"
granges <- import.bed(bed_file)
liftover <- liftOver(granges, chain_file)
saveRDS(liftover, file = "chromhmm/E12_18state_chromhmm.rds")

bed_file <- "chromhmm/E11_18state_chromhmm.bed"
granges <- import.bed(bed_file)
liftover <- liftOver(granges, chain_file)
saveRDS(liftover, file = "chromhmm/E11_18state_chromhmm.rds")

bed_file <- "chromhmm/E12_18state_chromhmm_200.bed"
granges <- import.bed(bed_file)
liftover <- liftOver(granges, chain_file)
saveRDS(liftover, file = "chromhmm/E12_18state_chromhmm_200.rds")

bed_file <- "chromhmm/E11_18state_chromhmm_200.bed"
granges <- import.bed(bed_file)
liftover <- liftOver(granges, chain_file)
saveRDS(liftover, file = "chromhmm/E11_18state_chromhmm_200.rds")
