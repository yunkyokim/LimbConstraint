###
### Correlative loop analysis between Micro-C and RNA-seq for clamped vs unclamped samples
### Analogous methods used for distal vs proximal HiC

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(ggmosaic)
library(ggpol)
library(ggridges)

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

##### Expression Annotation
set.seed(123456789)
mat <- read.csv("deseq_counts_clamp_only.csv")
diffmat <- read.csv("res_clamp_only.csv")
mat$comp_diff <- rowMeans(mat[, 3:4]) / rowMeans(mat[, 5:6])
mat$clamped <- rowMeans(mat[, 3:4])
mat$unclamped <- rowMeans(mat[, 5:6])
mat <- mat[!is.na(mat$comp_diff), ]
mat <- mat[!is.nan(mat$comp_diff), ]
mat <- mat[is.finite(mat$comp_diff), ]
mat <- mat[abs(mat$comp_diff) > 0, ]
orgmat <- mat

mat <- mat[mat$DISTCOM2 >= 10 &
  mat$DISTCOM4 >= 10 |
  mat$DISTUNCOM2 >= 10 &
    mat$DISTUNCOM4 >= 10, ]
mat <- mat[mat$symbol %in% diffmat$symbol, ]

feature_common <- readRDS("loops_merged_hictable_ann_chip") # load in annotated data
feature_common <- feature_common[feature_common$A1_gene != "" | feature_common$A2_gene != "", ]
feature_common <- feature_common[str_detect(feature_common[, "A1_ann"], "Tss") > 0 & # filter for enhancer-promoter interactions
  str_detect(feature_common[, "A2_ann"], "Enh") > 0 &
  feature_common$A1_gene != "" |
  str_detect(feature_common[, "A2_ann"], "Tss") > 0 &
    str_detect(feature_common[, "A1_ann"], "Enh") > 0 &
    feature_common$A2_gene != "", ]

for (i in 1:length(rownames(feature_common))) {
  if (str_detect(feature_common[i, "A1_ann"], "Tss") > 0 & str_detect(feature_common[i, "A2_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common[i, "A1_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common[i, "genelist"] <- feature_common[i, "A1_gene"]
  } else if (str_detect(feature_common[i, "A2_ann"], "Tss") > 0 & str_detect(feature_common[i, "A1_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common[i, "A2_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common[i, "genelist"] <- feature_common[i, "A2_gene"]
  } else {
    anchor_genes <- c(strsplit(feature_common[i, "A1_gene"], split = "/"), strsplit(feature_common[i, "A2_gene"], split = "/"))
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common[i, "genelist"] <- paste(feature_common[i, "A1_gene"], "/", feature_common[i, "A2_gene"], sep = "")
  }
  if (length(anchor_genes) > 0) {
    feature_common[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, 7])
    feature_common[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, "comp_diff"], length(intersect(mat$symbol, anchor_genes))))
    feature_common[i, "mean_exp_anchor_clamped"] <- mean(mat[mat$symbol %in% anchor_genes, "clamped"])
    feature_common[i, "mean_exp_anchor_unclamped"] <- mean(mat[mat$symbol %in% anchor_genes, "unclamped"])
  }
}

feature_common <- feature_common[!is.na(feature_common$mean_exp_anchor), ]
feature_common <- feature_common[!is.na(feature_common$hic_diff), ]
feature_common$mean_exp_anchor <- log2(feature_common$mean_exp_anchor)
feature_common <- feature_common[is.finite(feature_common$mean_exp_anchor), ]
feature_common$mean_exp_rand <- log2(feature_common$mean_exp_rand)
feature_common <- feature_common[is.finite(feature_common$mean_exp_anchor), ]
feature_common$adj.IF1 <- as.numeric(feature_common$adj.IF1)
feature_common <- feature_common[!is.na(feature_common$adj.IF1), ]
feature_common$adj.IF2 <- as.numeric(feature_common$adj.IF2)
feature_common <- feature_common[!is.na(feature_common$adj.IF2), ]
feature_common <- feature_common[!is.na(feature_common$mean_exp_anchor_clamped), ]
feature_common <- feature_common[!is.na(feature_common$mean_exp_anchor_unclamped), ]
feature_common$loop_length <- log2(as.numeric(feature_common$end2) - as.numeric(feature_common$start1))

nrow(feature_common) - nrow(feature_common[feature_common$hic_diff <= -log2(1.25), ]) - nrow(feature_common[feature_common$hic_diff >= log2(1.25), ]) # 12749
nrow(feature_common[feature_common$hic_diff >= log2(1.25), ]) # 443
nrow(feature_common[feature_common$hic_diff <= -log2(1.25), ]) # 904

##### Venn Diagram Plotting
data <- data.frame(
  group = c("DEDOWN", "NONDE", "DEUP"),
  value = c(
    nrow(feature_common[feature_common$hic_diff <= -log2(1.25), ]),
    nrow(feature_common) - nrow(feature_common[feature_common$hic_diff <= -log2(1.25), ]) - nrow(feature_common[feature_common$hic_diff >= log2(1.25), ]),
    nrow(feature_common[feature_common$hic_diff >= log2(1.25), ])
  )
)

sampleorder <- c("DEUP", "NONDE", "DEDOWN")
ggplot(data, aes(x = "", y = value, fill = factor(group, level = sampleorder))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 1.57) +
  theme_void() +
  NoLegend() +
  scale_fill_manual(values = c("NONDE" = "grey90", "DEUP" = "#4D9221FF", "DEDOWN" = "#C51B7DFF"))

##### Correlative analysis between Micro-C and RNA-seq for enhancer-promoter loops
bin_corr <- data.frame()
width <- 20
for (i in 0:14) {
  loop_bin <- feature_common[
    feature_common$adj.IF1 >= 0 + i * width | feature_common$adj.IF2 >= 0 + i * width,
    c("hic_diff", "mean_exp_anchor")
  ]
  loop_cor <- cor.test(loop_bin$hic_diff, loop_bin$mean_exp_anchor, method = c("spearman"))
  bin_corr[i + 1, "bin"] <- paste("> ", 0 + i * width, sep = "")
  bin_corr[i + 1, "value"] <- loop_cor$estimate[[1]]
  bin_corr[i + 1, "pvalue"] <- loop_cor$p.value
  bin_corr[i + 1, "npvalue"] <- -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.position = c(.5, .8),
    legend.title = element_text(size = 10),
    legend.key.width = unit(0.65, "cm"),
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  labs(color = "Density", title = "", fill = expression("-log"["10"] ~ "\n" * italic(P) * "-value")) +
  xlab("Loop Intensity Threshold") +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-Loop Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.15, 0.15) +
  scale_x_discrete(labels = c("> 0", "", "", "> 60", "", "", "> 120", "", "", "> 180", "", "", "> 240", "", ""))


##### Correlative analysis within conditions
##### Clamped condition
bin_corr <- data.frame()
width <- 20
for (i in 0:14) {
  loop_bin <- feature_common[
    feature_common$adj.IF1 >= 0 + i * width | feature_common$adj.IF2 >= 0 + i * width,
    c("adj.IF1", "mean_exp_anchor_clamped")
  ]
  loop_cor <- cor.test(loop_bin$adj.IF1, loop_bin$mean_exp_anchor_clamped, method = c("spearman"))
  bin_corr[i + 1, "bin"] <- paste("> ", 0 + i * width, sep = "")
  bin_corr[i + 1, "value"] <- loop_cor$estimate[[1]]
  bin_corr[i + 1, "pvalue"] <- loop_cor$p.value
  bin_corr[i + 1, "npvalue"] <- -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.position = c(.5, .25),
    legend.title = element_text(size = 10),
    legend.key.width = unit(0.65, "cm"),
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  labs(color = "Density", title = "", fill = expression("-log"["10"] ~ "\n" * italic(P) * "-value")) +
  xlab("Loop Intensity Threshold") +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-Loop Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.15, 0.15) +
  scale_x_discrete(labels = c("> 0", "", "", "> 60", "", "", "> 120", "", "", "> 180", "", "", "> 240", "", "")) # comp_clamped_loops 450x260

##### Unclamped condition
bin_corr <- data.frame()
width <- 20
for (i in 0:14) {
  loop_bin <- feature_common[
    feature_common$adj.IF1 >= 0 + i * width | feature_common$adj.IF2 >= 0 + i * width,
    c("adj.IF2", "mean_exp_anchor_unclamped")
  ]
  loop_cor <- cor.test(loop_bin$adj.IF2, loop_bin$mean_exp_anchor_unclamped, method = c("spearman"))
  bin_corr[i + 1, "bin"] <- paste("> ", 0 + i * width, sep = "")
  bin_corr[i + 1, "value"] <- loop_cor$estimate[[1]]
  bin_corr[i + 1, "pvalue"] <- loop_cor$p.value
  bin_corr[i + 1, "npvalue"] <- -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.position = c(.5, .25),
    legend.title = element_text(size = 10),
    legend.key.width = unit(0.65, "cm"),
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  labs(color = "Density", title = "", fill = expression("-log"["10"] ~ "\n" * italic(P) * "-value")) +
  xlab("Loop Intensity Threshold") +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-Loop Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.15, 0.15) +
  scale_x_discrete(labels = c("> 0", "", "", "> 60", "", "", "> 120", "", "", "> 180", "", "", "> 240", "", ""))


##### Random promoter control correlations
bin_corr <- data.frame()
width <- 20
for (i in 0:14) {
  loop_bin <- feature_common[
    feature_common$adj.IF1 >= 0 + i * width | feature_common$adj.IF2 >= 0 + i * width,
    c("hic_diff", "mean_exp_rand")
  ]
  loop_cor <- cor.test(loop_bin$hic_diff, loop_bin$mean_exp_rand, method = c("spearman"))
  bin_corr[i + 1, "bin"] <- paste("> ", 0 + i * width, sep = "")
  bin_corr[i + 1, "value"] <- loop_cor$estimate[[1]]
  bin_corr[i + 1, "pvalue"] <- loop_cor$p.value
  bin_corr[i + 1, "npvalue"] <- -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.position = c(.5, .8),
    legend.title = element_text(size = 10),
    legend.key.width = unit(0.65, "cm"),
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  labs(color = "Density", title = "", fill = expression("-log"["10"] ~ "\n" * italic(P) * "-value")) +
  xlab("Loop Intensity Threshold") +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-Loop Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.15, 0.15) +
  scale_x_discrete(labels = c("> 0", "", "", "> 60", "", "", "> 120", "", "", "> 180", "", "", "> 240", "", ""))

##### Correlative analysis within conditions, random promoters
##### Clamped condition, random promoters
bin_corr <- data.frame()
width <- 20
for (i in 0:14) {
  loop_bin <- feature_common[
    feature_common$adj.IF1 >= 0 + i * width | feature_common$adj.IF2 >= 0 + i * width,
    c("adj.IF1", "mean_exp_rand")
  ]
  loop_cor <- cor.test(loop_bin$adj.IF1, loop_bin$mean_exp_rand, method = c("spearman"))
  bin_corr[i + 1, "bin"] <- paste("> ", 0 + i * width, sep = "")
  bin_corr[i + 1, "value"] <- loop_cor$estimate[[1]]
  bin_corr[i + 1, "pvalue"] <- loop_cor$p.value
  bin_corr[i + 1, "npvalue"] <- -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.position = c(.5, .8),
    legend.title = element_text(size = 10),
    legend.key.width = unit(0.65, "cm"),
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  labs(color = "Density", title = "", fill = expression("-log"["10"] ~ "\n" * italic(P) * "-value")) +
  xlab("Loop Intensity Threshold") +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-Loop Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.15, 0.15) +
  scale_x_discrete(labels = c("> 0", "", "", "> 60", "", "", "> 120", "", "", "> 180", "", "", "> 240", "", ""))

##### Unclamped condition, random promoters
bin_corr <- data.frame()
width <- 20
for (i in 0:14) {
  loop_bin <- feature_common[
    feature_common$adj.IF1 >= 0 + i * width | feature_common$adj.IF2 >= 0 + i * width,
    c("adj.IF2", "mean_exp_rand")
  ]
  loop_cor <- cor.test(loop_bin$adj.IF2, loop_bin$mean_exp_rand, method = c("spearman"))
  bin_corr[i + 1, "bin"] <- paste("> ", 0 + i * width, sep = "")
  bin_corr[i + 1, "value"] <- loop_cor$estimate[[1]]
  bin_corr[i + 1, "pvalue"] <- loop_cor$p.value
  bin_corr[i + 1, "npvalue"] <- -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.position = c(.5, .8),
    legend.title = element_text(size = 10),
    legend.key.width = unit(0.65, "cm"),
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  labs(color = "Density", title = "", fill = expression("-log"["10"] ~ "\n" * italic(P) * "-value")) +
  xlab("Loop Intensity Threshold") +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-Loop Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.15, 0.15) +
  scale_x_discrete(labels = c("> 0", "", "", "> 60", "", "", "> 120", "", "", "> 180", "", "", "> 240", "", ""))


##### Gene enrichment testing - variable or non-variable loops in DEG set
difflist <- diffmat$symbol
nondifflist <- setdiff(mat$symbol, diffmat$symbol)
nondifflist <- nondifflist[!nondifflist %in% c("a", "T")]

threshold <- 20
feature_common2 <- feature_common[feature_common$adj.IF1 > threshold | feature_common$adj.IF2 > threshold, ]

anchorgenes <- c(strsplit(feature_common2[abs(feature_common2$hic_diff) >= log2(1.25), "genelist"], "/")) # Differential genes
anchorgenes <- unlist(anchorgenes)
anchorgenes <- anchorgenes[nchar(anchorgenes) > 0]
anchorgenes <- unique(anchorgenes)

anchorgenes2 <- c(strsplit(feature_common2[abs(feature_common2$hic_diff) <= log2(1.25), "genelist"], "/")) # Non-Differential genes
anchorgenes2 <- unlist(anchorgenes2)
anchorgenes2 <- anchorgenes2[nchar(anchorgenes2) > 0]
anchorgenes2 <- unique(anchorgenes2)

dat <- data.frame(
  "Variable Loop" = c(length(intersect(difflist, unique(anchorgenes))), length(intersect(setdiff(nondifflist, difflist), unique(anchorgenes)))),
  "Non-Variable Loop" = c(length(intersect(difflist, unique(anchorgenes2))), length(intersect(setdiff(nondifflist, difflist), unique(anchorgenes2)))),
  row.names = c("DEG", "Non-DEG"), stringsAsFactors = FALSE
)

dat <- as.data.frame(t(dat))
mosaicplot(dat, main = "Mosaic plot", color = TRUE)
fisher.test(dat)

dat$variable <- c("Variable Loop", "Non-Variable Loop")
dat_long <- gather(dat, key = "deg", value = "count", -c("variable"))

ggplot(data = dat_long) +
  geom_mosaic(aes(x = variable, fill = deg, weight = count), offset = 0.025, alpha = 1) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent"),
    aspect.ratio = 1.4,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5), line = element_blank()
  ) +
  scale_fill_manual(values = c("deepskyblue3", "grey90")) +
  xlab("") +
  ylab("") +
  NoLegend()


##### Gene enrichment testing - variable or non-variable loops per gene
threshold <- 20
feature_common2 <- feature_common[feature_common$adj.IF1 > threshold | feature_common$adj.IF2 > threshold, ]

loopspergene <- data.frame(symbol = diffmat$symbol)
for (i in 1:nrow(loopspergene)) {
  loop_subset <- feature_common2[str_detect(feature_common2$genelist, paste("\\b", loopspergene[i, "symbol"], "\\b", sep = "")), c("genelist", "hic_diff")]
  if (length(rownames(loop_subset)) > 0) {
    loopspergene[i, "number_diff"] <- nrow(loop_subset[abs(loop_subset$hic_diff) > log2(1.25), ])
    loopspergene[i, "number_nondiff"] <- nrow(loop_subset[abs(loop_subset$hic_diff) < log2(1.25), ])
  } else {
    loopspergene[i, "number_diff"] <- 0
    loopspergene[i, "number_nondiff"] <- 0
  }
}

loopspergene2 <- data.frame(symbol = setdiff(mat$symbol, diffmat$symbol))
for (i in 1:nrow(loopspergene2)) {
  loop_subset <- feature_common2[str_detect(feature_common2$genelist, paste("\\b", loopspergene2[i, "symbol"], "\\b", sep = "")), c("genelist", "hic_diff")]
  if (length(rownames(loop_subset)) > 0) {
    loopspergene2[i, "number_diff"] <- nrow(loop_subset[abs(loop_subset$hic_diff) > log2(1.25), ])
    loopspergene2[i, "number_nondiff"] <- nrow(loop_subset[abs(loop_subset$hic_diff) < log2(1.25), ])
  } else {
    loopspergene2[i, "number_diff"] <- 0
    loopspergene2[i, "number_nondiff"] <- 0
  }
}

loopspergene <- loopspergene[loopspergene$number_nondiff > 0, ]
loopspergene$ind <- loopspergene$number_diff / loopspergene$number_nondiff
loopspergene2 <- loopspergene2[loopspergene2$number_nondiff > 0, ]
loopspergene2$ind <- loopspergene2$number_diff / loopspergene2$number_nondiff

wilcox.test(loopspergene$ind, loopspergene2$ind)

loopspergene$condition <- "DEG"
loopspergene2$condition <- "Non-DEG"

loopspergene <- rbind(loopspergene, loopspergene2)
loopspergene$condition <- factor(loopspergene$condition)
loopspergene$ind2 <- log2(loopspergene$ind + 1)
my_comparisons <- list(c("DEG", "Non-DEG"))

ggplot(loopspergene[loopspergene$ind <= 3, ], aes(x = ind, color = factor(condition, levels = c("Non-DEG", "DEG")))) +
  geom_density(linewidth = 1) +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(color = "Density", title = "") +
  xlim(0, 3) +
  xlab("Variable/Non-Variable Loops") +
  ylab("Density") +
  NoLegend() +
  scale_y_continuous(trans = "sqrt") +
  scale_x_continuous(trans = "sqrt") +
  scale_color_manual(values = c("gray70", "deepskyblue3"))

##### Enhancer-Promoter ChIP-Seq Correlations
threshold <- 40
feature_common2 <- feature_common[feature_common$adj.IF1 >= threshold | feature_common$adj.IF2 >= threshold, ]
chip_list <- c(
  "E11_ATAC_WL_NP", "E12_RNAPII_BP", "E11_CTCF_NP", "E11_SMC1A_NP",
  "E11_H3K27ac_NP", "E11_H3K27me3_BP", "E11_H3K36me3_BP", "E11_H3K4me1_BP",
  "E11_H3K4me2_NP", "E11_H3K4me3_NP"
)
loop_cor_test <- cor.test(feature_common2$hic_diff, feature_common2$mean_exp_anchor, method = ("spearman"))
loop_chip <- data.frame(pvalue = loop_cor_test$p.value, rho = loop_cor_test$estimate)

for (i in 1:length(chip_list)) {
  for (j in 1:nrow(feature_common2)) {
    print(c(i, j))
    if (str_detect(feature_common2[j, "A1_ann"], "Tss") > 0 & str_detect(feature_common2[j, "A2_ann"], "Tss") == 0) {
      feature_common2[j, "chip_temp"] <- feature_common2[j, paste(chip_list[i], "_A1", sep = "")]
    } else if (str_detect(feature_common2[j, "A2_ann"], "Tss") > 0 & str_detect(feature_common2[j, "A1_ann"], "Tss") == 0) {
      feature_common2[j, "chip_temp"] <- feature_common2[j, paste(chip_list[i], "_A2", sep = "")]
    } else {
      feature_common2[j, "chip_temp"] <- rowMeans(feature_common2[j, c(paste(chip_list[i], "_A1", sep = ""), paste(chip_list[i], "_A2", sep = ""))])
    }
  }
  loop_cor_test <- cor.test(feature_common2$hic_diff, feature_common2$chip_temp, method = ("spearman"))
  loop_chip[i + 1, ] <- data.frame(pvalue = loop_cor_test$p.value, rho = loop_cor_test$estimate)
}

loop_chip$label <- c("Diff RNA", "ATAC", "RNAPII", "CTCF", "Smc1a", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3")
loop_chip$npvalue <- -log(loop_chip$pvalue)
label_order <- loop_chip[order(loop_chip$rho), "label"]
# label_order = c("H3K27ac","Diff RNA","CTCF","H3K4me1","ATAC","H3K36me3","RNAPII","H3K4me3","Smc1a","H3K4me2","H3K27me3")

heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight = 0.2)
ggplot(loop_chip, aes(x = reorder(label, rho), y = rho, fill = npvalue)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
  theme(
    aspect.ratio = 0.6,
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.position = c(.28, .8),
    legend.title = element_text(size = 10),
    legend.key.width = unit(0.65, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key = element_rect(fill = "transparent", colour = "transparent")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  labs(color = "Density", title = "", fill = expression("-log"["10"] ~ "\n" * italic(P) * "-value")) +
  xlab("") +
  ylab(bquote(atop(NA, atop(
    textstyle("Loop-ChIP Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.1, 0.1)


##### Enhancer-Promoter Breakdown by Enhancer Subset
##### Active Enhancers
feature_common_active <- readRDS("loops_merged_hictable_ann_chip_merge")
feature_common_active <- feature_common_active[feature_common_active$A1_gene != "" | feature_common_active$A2_gene != "", ]
feature_common_active$A1_ann <- paste0(feature_common_active$A1_ann, "/")
feature_common_active$A2_ann <- paste0(feature_common_active$A2_ann, "/")

feature_common_active <- feature_common_active[str_detect(feature_common_active[, "A1_ann"], "Tss") > 0 &
  str_detect(feature_common_active[, "A2_ann"], "Enh") > 0 &
  str_detect(feature_common_active[, "A2_ann"], "EnhPr") == 0 &
  str_detect(feature_common_active[, "A2_ann"], "EnhLo") == 0 &
  str_detect(feature_common_active[, "A2_ann"], "EnhPois") == 0 &
  feature_common_active$A1_gene != "" |
  str_detect(feature_common_active[, "A2_ann"], "Tss") > 0 &
    str_detect(feature_common_active[, "A1_ann"], "Enh") > 0 &
    str_detect(feature_common_active[, "A1_ann"], "EnhPr") == 0 &
    str_detect(feature_common_active[, "A1_ann"], "EnhLo") == 0 &
    str_detect(feature_common_active[, "A1_ann"], "EnhPois") == 0 &
    feature_common_active$A2_gene != "", ]

for (i in 1:length(rownames(feature_common_active))) {
  if (str_detect(feature_common_active[i, "A1_ann"], "Tss") > 0 & str_detect(feature_common_active[i, "A2_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common_active[i, "A1_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_active[i, "genelist"] <- feature_common_active[i, "A1_gene"]
  } else if (str_detect(feature_common_active[i, "A2_ann"], "Tss") > 0 & str_detect(feature_common_active[i, "A1_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common_active[i, "A2_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_active[i, "genelist"] <- feature_common_active[i, "A2_gene"]
  } else {
    anchor_genes <- c(strsplit(feature_common_active[i, "A1_gene"], split = "/"), strsplit(feature_common_active[i, "A2_gene"], split = "/"))
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_active[i, "genelist"] <- paste(feature_common_active[i, "A1_gene"], "/", feature_common_active[i, "A2_gene"], sep = "")
  }
  if (length(anchor_genes) > 0) {
    feature_common_active[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, 7])
    feature_common_active[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, 7], length(anchor_genes)))
    feature_common_active[i, "mean_exp_anchor_clamped"] <- mean(mat[mat$symbol %in% anchor_genes, "clamped"])
    feature_common_active[i, "mean_exp_anchor_unclamped"] <- mean(mat[mat$symbol %in% anchor_genes, "unclamped"])
  }
}

feature_common_active <- feature_common_active[!is.na(feature_common_active$mean_exp_anchor), ]
feature_common_active <- feature_common_active[!is.na(feature_common_active$hic_diff), ]
feature_common_active$mean_exp_anchor <- log2(feature_common_active$mean_exp_anchor)
feature_common_active <- feature_common_active[is.finite(feature_common_active$mean_exp_anchor), ]
feature_common_active$mean_exp_rand <- log2(feature_common_active$mean_exp_rand)
feature_common_active <- feature_common_active[is.finite(feature_common_active$mean_exp_anchor), ]
feature_common_active$adj.IF1 <- as.numeric(feature_common_active$adj.IF1)
feature_common_active <- feature_common_active[!is.na(feature_common_active$adj.IF1), ]
feature_common_active$adj.IF2 <- as.numeric(feature_common_active$adj.IF2)
feature_common_active <- feature_common_active[!is.na(feature_common_active$adj.IF2), ]
feature_common_active <- feature_common_active[!is.na(feature_common_active$mean_exp_anchor_clamped), ]
feature_common_active <- feature_common_active[!is.na(feature_common_active$mean_exp_anchor_unclamped), ]
feature_common_active$loop_length <- log2(as.numeric(feature_common_active$end2) - as.numeric(feature_common_active$start1))


##### Poised Enhancers
feature_common_poised <- readRDS("loops_merged_hictable_ann_chip_merge")
feature_common_poised <- feature_common_poised[feature_common_poised$A1_gene != "" | feature_common_poised$A2_gene != "", ]
feature_common_poised$A1_ann <- paste0(feature_common_poised$A1_ann, "/")
feature_common_poised$A2_ann <- paste0(feature_common_poised$A2_ann, "/")
feature_common_poised <- feature_common_poised[str_detect(feature_common_poised[, "A1_ann"], "Tss") > 0 &
  str_detect(feature_common_poised[, "A2_ann"], "Enh/") == 0 &
  str_detect(feature_common_poised[, "A2_ann"], "EnhPr") == 0 &
  str_detect(feature_common_poised[, "A2_ann"], "EnhLo") == 0 &
  str_detect(feature_common_poised[, "A2_ann"], "EnhPois") > 0 &
  feature_common_poised$A1_gene != "" |
  str_detect(feature_common_poised[, "A2_ann"], "Tss") > 0 &
    str_detect(feature_common_poised[, "A1_ann"], "Enh/") == 0 &
    str_detect(feature_common_poised[, "A1_ann"], "EnhPr") == 0 &
    str_detect(feature_common_poised[, "A1_ann"], "EnhLo") == 0 &
    str_detect(feature_common_poised[, "A1_ann"], "EnhPois") > 0 &
    feature_common_poised$A2_gene != "", ]

for (i in 1:length(rownames(feature_common_poised))) {
  if (str_detect(feature_common_poised[i, "A1_ann"], "Tss") > 0 & str_detect(feature_common_poised[i, "A2_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common_poised[i, "A1_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_poised[i, "genelist"] <- feature_common_poised[i, "A1_gene"]
  } else if (str_detect(feature_common_poised[i, "A2_ann"], "Tss") > 0 & str_detect(feature_common_poised[i, "A1_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common_poised[i, "A2_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_poised[i, "genelist"] <- feature_common_poised[i, "A2_gene"]
  } else {
    anchor_genes <- c(strsplit(feature_common_poised[i, "A1_gene"], split = "/"), strsplit(feature_common_poised[i, "A2_gene"], split = "/"))
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_poised[i, "genelist"] <- paste(feature_common_poised[i, "A1_gene"], "/", feature_common_poised[i, "A2_gene"], sep = "")
  }
  if (length(anchor_genes) > 0) {
    feature_common_poised[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, 7])
    feature_common_poised[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, 7], length(anchor_genes)))
    feature_common_poised[i, "mean_exp_anchor_clamped"] <- mean(mat[mat$symbol %in% anchor_genes, "clamped"])
    feature_common_poised[i, "mean_exp_anchor_unclamped"] <- mean(mat[mat$symbol %in% anchor_genes, "unclamped"])
  }
}

feature_common_poised <- feature_common_poised[!is.na(feature_common_poised$mean_exp_anchor), ]
feature_common_poised <- feature_common_poised[!is.na(feature_common_poised$hic_diff), ]
feature_common_poised$mean_exp_anchor <- log2(feature_common_poised$mean_exp_anchor)
feature_common_poised <- feature_common_poised[is.finite(feature_common_poised$mean_exp_anchor), ]
feature_common_poised$mean_exp_rand <- log2(feature_common_poised$mean_exp_rand)
feature_common_poised <- feature_common_poised[is.finite(feature_common_poised$mean_exp_anchor), ]
feature_common_poised$adj.IF1 <- as.numeric(feature_common_poised$adj.IF1)
feature_common_poised <- feature_common_poised[!is.na(feature_common_poised$adj.IF1), ]
feature_common_poised$adj.IF2 <- as.numeric(feature_common_poised$adj.IF2)
feature_common_poised <- feature_common_poised[!is.na(feature_common_poised$adj.IF2), ]
feature_common_poised <- feature_common_poised[!is.na(feature_common_poised$mean_exp_anchor_clamped), ]
feature_common_poised <- feature_common_poised[!is.na(feature_common_poised$mean_exp_anchor_unclamped), ]
feature_common_poised$loop_length <- log2(as.numeric(feature_common_poised$end2) - as.numeric(feature_common_poised$start1))

##### Low Enhancers
feature_common_low <- readRDS("loops_merged_hictable_ann_chip_merge")
feature_common_low <- feature_common_low[feature_common_low$A1_gene != "" | feature_common_low$A2_gene != "", ]
feature_common_low$A1_ann <- paste0(feature_common_low$A1_ann, "/")
feature_common_low$A2_ann <- paste0(feature_common_low$A2_ann, "/")
feature_common_low <- feature_common_low[str_detect(feature_common_low[, "A1_ann"], "Tss") > 0 &
  str_detect(feature_common_low[, "A2_ann"], "Enh/") == 0 &
  str_detect(feature_common_low[, "A2_ann"], "EnhPr") == 0 &
  str_detect(feature_common_low[, "A2_ann"], "EnhLo") > 0 &
  str_detect(feature_common_low[, "A2_ann"], "EnhPois") == 0 &
  feature_common_low$A1_gene != "" |
  str_detect(feature_common_low[, "A2_ann"], "Tss") > 0 &
    str_detect(feature_common_low[, "A1_ann"], "Enh/") == 0 &
    str_detect(feature_common_low[, "A1_ann"], "EnhPr") == 0 &
    str_detect(feature_common_low[, "A1_ann"], "EnhLo") > 0 &
    str_detect(feature_common_low[, "A1_ann"], "EnhPois") == 0 &
    feature_common_low$A2_gene != "", ]

for (i in 1:length(rownames(feature_common_low))) {
  if (str_detect(feature_common_low[i, "A1_ann"], "Tss") > 0 & str_detect(feature_common_low[i, "A2_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common_low[i, "A1_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_low[i, "genelist"] <- feature_common_low[i, "A1_gene"]
  } else if (str_detect(feature_common_low[i, "A2_ann"], "Tss") > 0 & str_detect(feature_common_low[i, "A1_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common_low[i, "A2_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_low[i, "genelist"] <- feature_common_low[i, "A2_gene"]
  } else {
    anchor_genes <- c(strsplit(feature_common_low[i, "A1_gene"], split = "/"), strsplit(feature_common_low[i, "A2_gene"], split = "/"))
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_low[i, "genelist"] <- paste(feature_common_low[i, "A1_gene"], "/", feature_common_low[i, "A2_gene"], sep = "")
  }
  if (length(anchor_genes) > 0) {
    feature_common_low[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, 7])
    feature_common_low[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, 7], length(anchor_genes)))
    feature_common_low[i, "mean_exp_anchor_clamped"] <- mean(mat[mat$symbol %in% anchor_genes, "clamped"])
    feature_common_low[i, "mean_exp_anchor_unclamped"] <- mean(mat[mat$symbol %in% anchor_genes, "unclamped"])
  }
}

feature_common_low <- feature_common_low[!is.na(feature_common_low$mean_exp_anchor), ]
feature_common_low <- feature_common_low[!is.na(feature_common_low$hic_diff), ]
feature_common_low$mean_exp_anchor <- log2(feature_common_low$mean_exp_anchor)
feature_common_low <- feature_common_low[is.finite(feature_common_low$mean_exp_anchor), ]
feature_common_low$mean_exp_rand <- log2(feature_common_low$mean_exp_rand)
feature_common_low <- feature_common_low[is.finite(feature_common_low$mean_exp_anchor), ]
feature_common_low$adj.IF1 <- as.numeric(feature_common_low$adj.IF1)
feature_common_low <- feature_common_low[!is.na(feature_common_low$adj.IF1), ]
feature_common_low$adj.IF2 <- as.numeric(feature_common_low$adj.IF2)
feature_common_low <- feature_common_low[!is.na(feature_common_low$adj.IF2), ]
feature_common_low <- feature_common_low[!is.na(feature_common_low$mean_exp_anchor_clamped), ]
feature_common_low <- feature_common_low[!is.na(feature_common_low$mean_exp_anchor_unclamped), ]
feature_common_low$loop_length <- log2(as.numeric(feature_common_low$end2) - as.numeric(feature_common_low$start1))


##### Primed Enhancers
feature_common_primed <- readRDS("loops_merged_hictable_ann_chip_merge")
feature_common_primed <- feature_common_primed[feature_common_primed$A1_gene != "" | feature_common_primed$A2_gene != "", ]
feature_common_primed$A1_ann <- paste0(feature_common_primed$A1_ann, "/")
feature_common_primed$A2_ann <- paste0(feature_common_primed$A2_ann, "/")
feature_common_primed <- feature_common_primed[str_detect(feature_common_primed[, "A1_ann"], "Tss") > 0 &
  str_detect(feature_common_primed[, "A2_ann"], "Enh/") == 0 &
  str_detect(feature_common_primed[, "A2_ann"], "EnhPr") > 0 &
  str_detect(feature_common_primed[, "A2_ann"], "EnhLo") == 0 &
  str_detect(feature_common_primed[, "A2_ann"], "EnhPois") == 0 &
  feature_common_primed$A1_gene != "" |
  str_detect(feature_common_primed[, "A2_ann"], "Tss") > 0 &
    str_detect(feature_common_primed[, "A1_ann"], "Enh/") == 0 &
    str_detect(feature_common_primed[, "A1_ann"], "EnhPr") > 0 &
    str_detect(feature_common_primed[, "A1_ann"], "EnhLo") == 0 &
    str_detect(feature_common_primed[, "A1_ann"], "EnhPois") == 0 &
    feature_common_primed$A2_gene != "", ]

for (i in 1:length(rownames(feature_common_primed))) {
  if (str_detect(feature_common_primed[i, "A1_ann"], "Tss") > 0 & str_detect(feature_common_primed[i, "A2_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common_primed[i, "A1_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_primed[i, "genelist"] <- feature_common_primed[i, "A1_gene"]
  } else if (str_detect(feature_common_primed[i, "A2_ann"], "Tss") > 0 & str_detect(feature_common_primed[i, "A1_ann"], "Tss") == 0) {
    anchor_genes <- strsplit(feature_common_primed[i, "A2_gene"], split = "/")
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_primed[i, "genelist"] <- feature_common_primed[i, "A2_gene"]
  } else {
    anchor_genes <- c(strsplit(feature_common_primed[i, "A1_gene"], split = "/"), strsplit(feature_common_primed[i, "A2_gene"], split = "/"))
    anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
    feature_common_primed[i, "genelist"] <- paste(feature_common_primed[i, "A1_gene"], "/", feature_common_primed[i, "A2_gene"], sep = "")
  }
  if (length(anchor_genes) > 0) {
    feature_common_primed[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, 7])
    feature_common_primed[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, 7], length(anchor_genes)))
    feature_common_primed[i, "mean_exp_anchor_clamped"] <- mean(mat[mat$symbol %in% anchor_genes, "clamped"])
    feature_common_primed[i, "mean_exp_anchor_unclamped"] <- mean(mat[mat$symbol %in% anchor_genes, "unclamped"])
  }
}

feature_common_primed <- feature_common_primed[!is.na(feature_common_primed$mean_exp_anchor), ]
feature_common_primed <- feature_common_primed[!is.na(feature_common_primed$hic_diff), ]
feature_common_primed$mean_exp_anchor <- log2(feature_common_primed$mean_exp_anchor)
feature_common_primed <- feature_common_primed[is.finite(feature_common_primed$mean_exp_anchor), ]
feature_common_primed$mean_exp_rand <- log2(feature_common_primed$mean_exp_rand)
feature_common_primed <- feature_common_primed[is.finite(feature_common_primed$mean_exp_anchor), ]
feature_common_primed$adj.IF1 <- as.numeric(feature_common_primed$adj.IF1)
feature_common_primed <- feature_common_primed[!is.na(feature_common_primed$adj.IF1), ]
feature_common_primed$adj.IF2 <- as.numeric(feature_common_primed$adj.IF2)
feature_common_primed <- feature_common_primed[!is.na(feature_common_primed$adj.IF2), ]
feature_common_primed <- feature_common_primed[!is.na(feature_common_primed$mean_exp_anchor_clamped), ]
feature_common_primed <- feature_common_primed[!is.na(feature_common_primed$mean_exp_anchor_unclamped), ]
feature_common_primed$loop_length <- log2(as.numeric(feature_common_primed$end2) - as.numeric(feature_common_primed$start1))


##### Enhancer Breakdown Plotting
##### Active Enhancers
threshold <- 40
loop_threshold <- log2(20000)
ggplot(feature_common_active[as.numeric(feature_common_active$adj.IF1) >= threshold &
  feature_common_active$loop_length >= loop_threshold |
  as.numeric(feature_common_active$adj.IF2) >= threshold &
    feature_common_active$loop_length >= loop_threshold, ], aes(x = hic_diff, y = mean_exp_anchor)) +
  scale_color_viridis(option = "F") +
  geom_jitter(aes(color = get_density(hic_diff, mean_exp_anchor, n = 500)), size = 0.1) +
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
  ylim(-1.5, 1.5) +
  xlim(-1.5, 1.5)

##### Poised Enhancers
ggplot(feature_common_poised[as.numeric(feature_common_poised$adj.IF1) >= threshold &
  feature_common_poised$loop_length >= loop_threshold |
  as.numeric(feature_common_poised$adj.IF2) >= threshold &
    feature_common_poised$loop_length >= loop_threshold, ], aes(x = hic_diff, y = mean_exp_anchor)) +
  scale_color_viridis(option = "F") +
  geom_jitter(aes(color = get_density(hic_diff, mean_exp_anchor, n = 500)), size = 0.1) +
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
  ylim(-1.5, 1.5) +
  xlim(-1.5, 1.5)

##### Primed Enhancers
ggplot(feature_common_primed[as.numeric(feature_common_primed$adj.IF1) >= threshold &
  feature_common_primed$loop_length >= loop_threshold |
  as.numeric(feature_common_primed$adj.IF2) >= threshold &
    feature_common_primed$loop_length >= loop_threshold, ], aes(x = hic_diff, y = mean_exp_anchor)) +
  scale_color_viridis(option = "F") +
  geom_jitter(aes(color = get_density(hic_diff, mean_exp_anchor, n = 500)), size = 0.1) +
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
  ylim(-1.5, 1.5) +
  xlim(-1.5, 1.5)

##### Low Enhancers
ggplot(feature_common_low[as.numeric(feature_common_low$adj.IF1) >= threshold &
  feature_common_low$loop_length >= loop_threshold |
  as.numeric(feature_common_low$adj.IF2) >= threshold &
    feature_common_low$loop_length >= loop_threshold, ], aes(x = hic_diff, y = mean_exp_anchor)) +
  scale_color_viridis(option = "F") +
  geom_jitter(aes(color = get_density(hic_diff, mean_exp_anchor, n = 500)), size = 0.1) +
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
  ylim(-1.5, 1.5) +
  xlim(-1.5, 1.5)
