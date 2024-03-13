###
### Comparative analysis between loop-gene correlations in Micro-C vs Hi-C datasets
###

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyverse)
library(reshape)
library(Seurat)
library(dplyr)
library(MASS)
library(ggplot2)
library(viridis)
library(multiHiCcompare)
library(HiCcompare)
library(ggpubr)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

options(scipen = 999)
setwd("./")
dir <- getwd()

##### Distal vs Proximal Hi-C - basic preprocessing
mat <- read.csv("E12_distprox_counts.csv")
diffmat <- read.csv("E12_distprox_DEG.csv")
mat$comp_diff <- rowMeans(mat[, 3:5]) / rowMeans(mat[, 6:8])
mat$DFL <- rowMeans(mat[, 3:5])
mat$PFL <- rowMeans(mat[, 6:8])
mat <- mat[!is.na(mat$comp_diff), ]
mat <- mat[!is.nan(mat$comp_diff), ]
mat <- mat[is.finite(mat$comp_diff), ]
mat <- mat[abs(mat$comp_diff) > 0, ]
orgmat <- mat

mat <- mat[mat$DFL1 >= 10 &
  mat$DFL2 >= 10 &
  mat$DFL3 >= 10 |
  mat$PFL1 >= 10 &
    mat$PFL2 >= 10 &
    mat$PFL3 >= 10, ]

loop_common <- readRDS("loops_dd_merged_hictable_ann_chip")
loop_common <- loop_common[loop_common$A1_gene != "" | loop_common$A2_gene != "", ]
loop_common <- loop_common[str_detect(loop_common[, "A1_ann"], "Tss") > 0 &
  str_detect(loop_common[, "A2_ann"], "Enh") > 0 &
  loop_common$A1_gene != "" |
  str_detect(loop_common[, "A2_ann"], "Tss") > 0 &
    str_detect(loop_common[, "A1_ann"], "Enh") > 0 &
    loop_common$A2_gene != "", ]

for (i in 1:length(rownames(loop_common))) {
  if (str_detect(loop_common[i, "A1_ann"], "Tss") > 0 & str_detect(loop_common[i, "A2_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(loop_common[i, "A1_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    loop_common[i, "genelist"] <- loop_common[i, "A1_gene"]
  } else if (str_detect(loop_common[i, "A2_ann"], "Tss") > 0 & str_detect(loop_common[i, "A1_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(loop_common[i, "A2_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    loop_common[i, "genelist"] <- loop_common[i, "A2_gene"]
  } else {
    anchor_genes <- c(strsplit(loop_common[i, "A1_gene"], split = "/"), strsplit(loop_common[i, "A2_gene"], split = "/"))
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    loop_common[i, "genelist"] <- paste(loop_common[i, "A1_gene"], "/", loop_common[i, "A2_gene"], sep = "")
  }
  if (length(anchor_genes) > 0) {
    loop_common[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, "comp_diff"])
    loop_common[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, "comp_diff"], length(anchor_genes)))
    loop_common[i, "mean_exp_anchor_CC"] <- mean(mat[mat$symbol %in% anchor_genes, "DFL"])
    loop_common[i, "mean_exp_anchor_UC"] <- mean(mat[mat$symbol %in% anchor_genes, "PFL"])
  }
}

loop_common <- loop_common[!is.na(loop_common$mean_exp_anchor), ]
loop_common <- loop_common[!is.na(loop_common$hic_diff), ]
loop_common$mean_exp_anchor <- log2(loop_common$mean_exp_anchor)
loop_common <- loop_common[is.finite(loop_common$mean_exp_anchor), ]
loop_common$mean_exp_rand <- log2(loop_common$mean_exp_rand)
loop_common <- loop_common[is.finite(loop_common$mean_exp_anchor), ]
loop_common$adj.IF1 <- as.numeric(loop_common$adj.IF1)
loop_common <- loop_common[!is.na(loop_common$adj.IF1), ]
loop_common$adj.IF2 <- as.numeric(loop_common$adj.IF2)
loop_common <- loop_common[!is.na(loop_common$adj.IF2), ]
loop_common$loop_length <- log2(as.numeric(loop_common$end2) - as.numeric(loop_common$start1))
loop_common_dd <- loop_common


##### Clamped vs Unclamped Micro-C - basic preprocessing
mat <- read.csv("E11_clamp_counts.csv")
diffmat <- read.csv("E11_clamp_DEG.csv")
mat$comp_diff <- rowMeans(mat[, 3:4]) / rowMeans(mat[, 5:6])
mat$CC <- rowMeans(mat[, 3:4])
mat$UC <- rowMeans(mat[, 5:6])
mat <- mat[!is.na(mat$comp_diff), ]
mat <- mat[!is.nan(mat$comp_diff), ]
mat <- mat[is.finite(mat$comp_diff), ]
mat <- mat[abs(mat$comp_diff) > 0, ]
orgmat <- mat

mat <- mat[mat$DISTCOM2 >= 10 &
  mat$DISTCOM4 >= 10 |
  mat$DISTUNCOM2 >= 10 &
    mat$DISTUNCOM4 >= 10, ]

loop_common <- readRDS("loops_merged_hictable_ann_chip")
loop_common <- loop_common[loop_common$A1_gene != "" | loop_common$A2_gene != "", ]
loop_common <- loop_common[str_detect(loop_common[, "A1_ann"], "Tss") > 0 &
  str_detect(loop_common[, "A2_ann"], "Enh") > 0 &
  loop_common$A1_gene != "" |
  str_detect(loop_common[, "A2_ann"], "Tss") > 0 &
    str_detect(loop_common[, "A1_ann"], "Enh") > 0 &
    loop_common$A2_gene != "", ]

for (i in 1:length(rownames(loop_common))) {
  if (str_detect(loop_common[i, "A1_ann"], "Tss") > 0 & str_detect(loop_common[i, "A2_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(loop_common[i, "A1_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    loop_common[i, "genelist"] <- loop_common[i, "A1_gene"]
  } else if (str_detect(loop_common[i, "A2_ann"], "Tss") > 0 & str_detect(loop_common[i, "A1_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(loop_common[i, "A2_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    loop_common[i, "genelist"] <- loop_common[i, "A2_gene"]
  } else {
    anchor_genes <- c(strsplit(loop_common[i, "A1_gene"], split = "/"), strsplit(loop_common[i, "A2_gene"], split = "/"))
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    loop_common[i, "genelist"] <- paste(loop_common[i, "A1_gene"], "/", loop_common[i, "A2_gene"], sep = "")
  }
  if (length(anchor_genes) > 0) {
    loop_common[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, 7])
    loop_common[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, 7], length(anchor_genes)))
    loop_common[i, "mean_exp_anchor_CC"] <- mean(mat[mat$symbol %in% anchor_genes, "CC"])
    loop_common[i, "mean_exp_anchor_UC"] <- mean(mat[mat$symbol %in% anchor_genes, "UC"])
  }
}

loop_common <- loop_common[!is.na(loop_common$mean_exp_anchor), ]
loop_common <- loop_common[!is.na(loop_common$hic_diff), ]
loop_common$mean_exp_anchor <- log2(loop_common$mean_exp_anchor)
loop_common <- loop_common[is.finite(loop_common$mean_exp_anchor), ]
loop_common$mean_exp_rand <- log2(loop_common$mean_exp_rand)
loop_common <- loop_common[is.finite(loop_common$mean_exp_anchor), ]
loop_common$adj.IF1 <- as.numeric(loop_common$adj.IF1)
loop_common <- loop_common[!is.na(loop_common$adj.IF1), ]
loop_common$adj.IF2 <- as.numeric(loop_common$adj.IF2)
loop_common <- loop_common[!is.na(loop_common$adj.IF2), ]
loop_common <- loop_common[!is.na(loop_common$mean_exp_anchor_CC), ]
loop_common <- loop_common[!is.na(loop_common$mean_exp_anchor_UC), ]
loop_common$loop_length <- log2(as.numeric(loop_common$end2) - as.numeric(loop_common$start1))
loop_common_comp <- loop_common

saveRDS(loop_common_comp, "loop_common_comp_switch")
saveRDS(loop_common_dd, "loop_common_dd_switch")

##### Shared Loop Annotation
loop_common_comp <- readRDS("loop_common_comp_switch")
loop_common_dd <- readRDS("loop_common_dd_switch")
loop_common_comp$switchID <- 1:nrow(loop_common_comp)
loop_common_dd$switchID <- 1:nrow(loop_common_dd)

##### Annotation loops based on pattern of shared anchors
for (i in 1:nrow(loop_common_dd)) {
  loop_subset <- loop_common_comp[loop_common_comp$chrom1 == loop_common_dd[i, "chrom1"] &
    loop_common_comp$start1 >= loop_common_dd[i, "start1"] &
    loop_common_comp$end1 <= loop_common_dd[i, "end1"] &
    loop_common_comp$start2 >= loop_common_dd[i, "start2"] &
    loop_common_comp$end2 <= loop_common_dd[i, "end2"], ]
  if (nrow(loop_subset) > 0) {
    loop_common_dd[i, "A1A2_shared"] <- "SHARED"
    loop_common_dd[i, "A1A2_shared_strength"] <- mean(loop_subset$hic_diff)
    loop_common_comp[loop_common_comp$switchID %in% loop_subset$switchID, "A1A2_shared"] <- loop_common_dd[i, "switchID"]
  } else {
    loop_common_dd[i, "A1A2_shared"] <- "UNIQUE"
  }
  print(paste(i, "in", nrow(loop_common_dd)))
}

for (i in 1:nrow(loop_common_dd)) {
  loop_subset <- loop_common_comp[loop_common_comp$chrom1 == loop_common_dd[i, "chrom1"] &
    loop_common_comp$start2 >= loop_common_dd[i, "start1"] &
    loop_common_comp$end2 <= loop_common_dd[i, "end1"] &
    loop_common_comp$start1 >= loop_common_dd[i, "start2"] &
    loop_common_comp$end1 <= loop_common_dd[i, "end2"], ]
  if (nrow(loop_subset) > 0) {
    loop_common_dd[i, "A2A1_shared"] <- "SHARED"
    loop_common_dd[i, "A2A1_shared_strength"] <- mean(loop_subset$hic_diff)
    loop_common_comp[loop_common_comp$switchID %in% loop_subset$switchID, "A2A1_shared"] <- loop_common_dd[i, "switchID"]
  } else {
    loop_common_dd[i, "A2A1_shared"] <- "UNIQUE"
  }
  print(paste(i, "in", nrow(loop_common_dd)))
}

for (i in 1:nrow(loop_common_dd)) {
  loop_subset <- loop_common_comp[loop_common_comp$chrom1 == loop_common_dd[i, "chrom1"] &
    loop_common_comp$start1 >= loop_common_dd[i, "start1"] &
    loop_common_comp$end1 <= loop_common_dd[i, "end1"] &
    loop_common_comp$start2 > loop_common_dd[i, "end2"], ]
  if (nrow(loop_subset) > 0) {
    loop_common_dd[i, "A1_shared_greater"] <- "SHARED"
    loop_common_dd[i, "A1_shared_greater_strength"] <- mean(loop_subset$hic_diff)
    loop_common_comp[loop_common_comp$switchID %in% loop_subset$switchID, "A1_shared_greater"] <- loop_common_dd[i, "switchID"]
  } else {
    loop_common_dd[i, "A1_shared_greater"] <- "UNIQUE"
  }
  print(paste(i, "in", nrow(loop_common_dd), ", nrow = ", nrow(loop_subset)))
}

for (i in 1:nrow(loop_common_dd)) {
  loop_subset <- loop_common_comp[loop_common_comp$chrom1 == loop_common_dd[i, "chrom1"] &
    loop_common_comp$start1 >= loop_common_dd[i, "start1"] &
    loop_common_comp$end1 <= loop_common_dd[i, "end1"] &
    loop_common_comp$end2 < loop_common_dd[i, "start2"], ]
  if (nrow(loop_subset) > 0) {
    loop_common_dd[i, "A1_shared_lesser"] <- "SHARED"
    loop_common_dd[i, "A1_shared_lesser_strength"] <- mean(loop_subset$hic_diff)
    loop_common_comp[loop_common_comp$switchID %in% loop_subset$switchID, "A1_shared_lesser"] <- loop_common_dd[i, "switchID"]
  } else {
    loop_common_dd[i, "A1_shared_lesser"] <- "UNIQUE"
  }
  print(paste(i, "in", nrow(loop_common_dd), ", nrow = ", nrow(loop_subset)))
}

for (i in 1:nrow(loop_common_dd)) {
  loop_subset <- loop_common_comp[loop_common_comp$chrom1 == loop_common_dd[i, "chrom1"] &
    loop_common_comp$start2 >= loop_common_dd[i, "start2"] &
    loop_common_comp$end2 <= loop_common_dd[i, "end2"] &
    loop_common_comp$end1 < loop_common_dd[i, "start1"], ]
  if (nrow(loop_subset) > 0) {
    loop_common_dd[i, "A2_shared_greater"] <- "SHARED"
    loop_common_dd[i, "A2_shared_greater_strength"] <- mean(loop_subset$hic_diff)
    loop_common_comp[loop_common_comp$switchID %in% loop_subset$switchID, "A2_shared_greater"] <- loop_common_dd[i, "switchID"]
  } else {
    loop_common_dd[i, "A2_shared_greater"] <- "UNIQUE"
  }
  print(paste(i, "in", nrow(loop_common_dd)))
}

for (i in 1:nrow(loop_common_dd)) {
  loop_subset <- loop_common_comp[loop_common_comp$chrom1 == loop_common_dd[i, "chrom1"] &
    loop_common_comp$start2 >= loop_common_dd[i, "start2"] &
    loop_common_comp$end2 <= loop_common_dd[i, "end2"] &
    loop_common_comp$start1 > loop_common_dd[i, "end1"], ]
  if (nrow(loop_subset) > 0) {
    loop_common_dd[i, "A2_shared_lesser"] <- "SHARED"
    loop_common_dd[i, "A2_shared_lesser_strength"] <- mean(loop_subset$hic_diff)
    loop_common_comp[loop_common_comp$switchID %in% loop_subset$switchID, "A2_shared_lesser"] <- loop_common_dd[i, "switchID"]
  } else {
    loop_common_dd[i, "A2_shared_lesser"] <- "UNIQUE"
  }
  print(paste(i, "in", nrow(loop_common_dd)))
}

##### Plotting Shared Loops - Permissive (Anchors must be shared, but can have other looping interactions)
dd_AA_shared <- loop_common_dd[loop_common_dd$A1A2_shared == "SHARED", ]
threshold <- 20
loop_threshold <- log2(20000)
ggplot(dd_AA_shared[as.numeric(dd_AA_shared$adj.IF1) >= threshold &
  dd_AA_shared$loop_length >= loop_threshold |
  as.numeric(dd_AA_shared$adj.IF2) >= threshold &
    dd_AA_shared$loop_length >= loop_threshold, ], aes(x = hic_diff, y = mean_exp_anchor)) +
  scale_color_viridis(option = "F") + # , limits = c(0,2.5)
  geom_jitter(aes(color = get_density(hic_diff, mean_exp_anchor, n = 500)), size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 1) +
  geom_smooth(method = lm, color = "seagreen2") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(color = "Density", title = "") +
  xlab(bquote(atop(NA, atop(
    textstyle("Loop Intensity Change"),
    textstyle("(Distal/Proximal)")
  )))) +
  ylab(bquote(atop(NA, atop(
    textstyle("log"["2"] ~ "Expresion Fold Change"),
    textstyle("(Distal/Proximal)")
  )))) +
  ylim(-1.0, 1.0) +
  xlim(-1.0, 1.0)


comp_AA_shared <- loop_common[is.finite(loop_common_comp$A1A2_shared), ]
threshold <- 20
loop_threshold <- log2(20000)
ggplot(comp_AA_shared[as.numeric(comp_AA_shared$adj.IF1) >= threshold &
  comp_AA_shared$loop_length >= loop_threshold |
  as.numeric(comp_AA_shared$adj.IF2) >= threshold &
    comp_AA_shared$loop_length >= loop_threshold, ], aes(x = hic_diff, y = mean_exp_anchor)) +
  scale_color_viridis(option = "F") + # , limits = c(0,2.5)
  geom_jitter(aes(color = get_density(hic_diff, mean_exp_anchor, n = 500)), size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 1) +
  geom_smooth(method = lm, color = "seagreen2") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(color = "Density", title = "") +
  xlab(bquote(atop(NA, atop(
    textstyle("Loop Intensity Change"),
    textstyle("(Compressed/Control)")
  )))) +
  ylab(bquote(atop(NA, atop(
    textstyle("log"["2"] ~ "Expresion Fold Change"),
    textstyle("(Compressed/Control)")
  )))) +
  ylim(-1.0, 1.0) +
  xlim(-1.0, 1.0)


##### Plotting Shared Loops - Strict (Anchors must be shared and have no other interactions)
dd_AA_shared_strict <- loop_common_dd[loop_common_dd$A1_shared_greater == "UNIQUE" &
  loop_common_dd$A1_shared_lesser == "UNIQUE" &
  loop_common_dd$A2_shared_greater == "UNIQUE" &
  loop_common_dd$A2_shared_lesser == "UNIQUE" &
  loop_common_dd$A1A2_shared == "SHARED", ]

threshold <- 20
loop_threshold <- log2(20000)
ggplot(dd_AA_shared_strict[as.numeric(dd_AA_shared_strict$adj.IF1) >= threshold &
  dd_AA_shared_strict$loop_length >= loop_threshold |
  as.numeric(dd_AA_shared_strict$adj.IF2) >= threshold &
    dd_AA_shared_strict$loop_length >= loop_threshold, ], aes(x = hic_diff, y = mean_exp_anchor)) +
  scale_color_viridis(option = "F") + # , limits = c(0,2.5)
  geom_jitter(aes(color = get_density(hic_diff, mean_exp_anchor, n = 500)), size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 1) +
  geom_smooth(method = lm, color = "seagreen2") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(color = "Density", title = "") +
  xlab(bquote(atop(NA, atop(
    textstyle("Loop Intensity Change"),
    textstyle("(Distal/Proximal)")
  )))) +
  ylab(bquote(atop(NA, atop(
    textstyle("log"["2"] ~ "Expresion Fold Change"),
    textstyle("(Distal/Proximal)")
  )))) +
  ylim(-1.0, 1.0) +
  xlim(-1.0, 1.0)



comp_AA_shared_strict <- loop_common_comp[is.na(loop_common_comp$A1_shared_greater) &
  is.na(loop_common_comp$A1_shared_lesser) &
  is.na(loop_common_comp$A2_shared_greater) &
  is.na(loop_common_comp$A2_shared_lesser) &
  is.finite(loop_common_comp$A1A2_shared), ]
threshold <- 20
loop_threshold <- log2(20000)
ggplot(comp_AA_shared_strict[as.numeric(comp_AA_shared_strict$adj.IF1) >= threshold &
  comp_AA_shared_strict$loop_length >= loop_threshold |
  as.numeric(comp_AA_shared_strict$adj.IF2) >= threshold &
    comp_AA_shared_strict$loop_length >= loop_threshold, ], aes(x = hic_diff, y = mean_exp_anchor)) +
  scale_color_viridis(option = "F") + # , limits = c(0,2.5)
  geom_jitter(aes(color = get_density(hic_diff, mean_exp_anchor, n = 500)), size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 1) +
  geom_smooth(method = lm, color = "seagreen2") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(color = "Density", title = "") +
  xlab(bquote(atop(NA, atop(
    textstyle("Loop Intensity Change"),
    textstyle("(Compressed/Control)")
  )))) +
  ylab(bquote(atop(NA, atop(
    textstyle("log"["2"] ~ "Expresion Fold Change"),
    textstyle("(Compressed/Control)")
  )))) +
  ylim(-1.0, 1.0) +
  xlim(-1.0, 1.0)
