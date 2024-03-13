###
### P(s) curve calculation from raw HiCPro matrices
###

library(dplyr)
library(HiCExperiment)
library(HiContacts)

options(scipen = 999)
setwd("./")
dir <- getwd()

hf <- HicproFile("distal.1000.matrix", "distal.1000_abs.bed")
distal <- import(hf)
hf <- HicproFile("proximal.1000.matrix", "proximal.1000_abs.bed")
proximal <- import(hf)
rm(hf)

ps_1 <- distanceLaw(distal) |> mutate(sample = "distal")
ps_2 <- distanceLaw(proximal) |> mutate(sample = "proximal")
ps_rep1 <- rbind(ps_1, ps_2)
saveRDS(ps_rep1, "ps_e12_hic.rds")
