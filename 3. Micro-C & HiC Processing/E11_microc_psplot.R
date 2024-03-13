###
### P(s) curve calculation from raw HiCPro matrices
###

library(dplyr)
library(HiCExperiment)
library(HiContacts)

options(scipen = 999)
setwd("./")
dir <- getwd()

hf <- HicproFile("clamped.1000.matrix", "clamped.1000_abs.bed")
clamped <- import(hf)
hf <- HicproFile("unclamped.1000.matrix", "unclamped.1000_abs.bed")
unclamped <- import(hf)
rm(hf)

ps_1 <- distanceLaw(clamped) |> mutate(sample = "clamped")
ps_2 <- distanceLaw(unclamped) |> mutate(sample = "unclamped")
ps_rep1 <- rbind(ps_1, ps_2)
saveRDS(ps_rep1, "ps_e11_microc.rds")
