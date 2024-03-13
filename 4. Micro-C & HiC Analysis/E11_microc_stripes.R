###
### Correlative stripe analysis between Micro-C and RNA-seq for clamped vs unclamped samples
### Analogous methods used for distal vs proximal HiC

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(Seurat)
library(dplyr)
library(ggplot2)

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

# left stripes
feature_common <- readRDS("stripes_left_merged_hictable_ann_chip")
feature_common <- feature_common[feature_common$A1_gene != "", ]
feature_common <- feature_common[str_detect(feature_common[, "A1_ann"], "Tss") > 0 & feature_common$A1_gene != "", ] # filter anchors for presence of promoter

for (i in 1:length(rownames(feature_common))) {
  anchor_genes <- c(strsplit(feature_common[i, "A1_gene"], split = "/"))
  anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
  feature_common[i, "genelist"] <- paste(feature_common[i, "A1_gene"], sep = "")

  if (length(anchor_genes) > 0) {
    feature_common[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, "comp_diff"])
    feature_common[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, "comp_diff"], length(intersect(mat$symbol, anchor_genes))))
    feature_common[i, "mean_exp_anchor_clamped"] <- mean(mat[mat$symbol %in% anchor_genes, "clamped"])
    feature_common[i, "mean_exp_anchor_unclamped"] <- mean(mat[mat$symbol %in% anchor_genes, "unclamped"])
  }
}
loop_left <- feature_common

# right stripes
feature_common <- readRDS("stripes_right_merged_hictable_ann_chip")
feature_common <- feature_common[feature_common$A2_gene != "", ]
feature_common <- feature_common[str_detect(feature_common[, "A2_ann"], "Tss") > 0 & feature_common$A2_gene != "", ]

for (i in 1:length(rownames(feature_common))) {
  anchor_genes <- c(strsplit(feature_common[i, "A2_gene"], split = "/"))
  anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
  feature_common[i, "genelist"] <- paste(feature_common[i, "A2_gene"], sep = "")

  if (length(anchor_genes) > 0) {
    feature_common[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, "comp_diff"])
    feature_common[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, "comp_diff"], length(intersect(mat$symbol, anchor_genes))))
    feature_common[i, "mean_exp_anchor_clamped"] <- mean(mat[mat$symbol %in% anchor_genes, "clamped"])
    feature_common[i, "mean_exp_anchor_unclamped"] <- mean(mat[mat$symbol %in% anchor_genes, "unclamped"])
  }
}
loop_right <- feature_common

names(loop_right)[names(loop_right) == "A2_ann"] <- "A1_ann" # swap orientation to account for different stripe directions
names(loop_right)[names(loop_right) == "A2_gene"] <- "A1_gene"

feature_common <- rbind(loop_right, loop_left)
feature_common <- feature_common[!is.na(feature_common$mean_exp_anchor), ]
feature_common <- feature_common[!is.na(feature_common$hic_diff), ]
feature_common <- feature_common[is.finite(feature_common$hic_diff), ]
feature_common$mean_exp_anchor <- log2(feature_common$mean_exp_anchor)
feature_common <- feature_common[is.finite(feature_common$mean_exp_anchor), ]
feature_common$mean_exp_rand <- log2(feature_common$mean_exp_rand)
feature_common$adj.IF1 <- as.numeric(feature_common$adj.IF1)
feature_common <- feature_common[!is.na(feature_common$adj.IF1), ]
feature_common$adj.IF2 <- as.numeric(feature_common$adj.IF2)
feature_common <- feature_common[!is.na(feature_common$adj.IF2), ]
feature_common$loop_length <- log2(as.numeric(feature_common$end2) - as.numeric(feature_common$start1))

nrow(feature_common) - nrow(feature_common[feature_common$hic_diff <= -log2(1.25), ]) - nrow(feature_common[feature_common$hic_diff >= log2(1.25), ]) # 8044
nrow(feature_common[feature_common$hic_diff >= log2(1.25), ]) # 1
nrow(feature_common[feature_common$hic_diff <= -log2(1.25), ]) # 107

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

##### Correlative Analysis
bin_corr <- data.frame()
width <- 150
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

##### Anchor analysis
heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 20), label = function(x) sprintf("%.0f", x)) +
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
  xlab("Stripe Intensity Threshold") +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-Stripe Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.25, 0.25) +
  scale_x_discrete(labels = c("> 0", "", "", "> 450", "", "", "> 900", "", "", "> 1350", "", "", "> 1800", "", ""))

##### Random Promoter Controls
bin_corr <- data.frame()
width <- 10
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
  scale_fill_gradientn(colors = heat_colors, na.value = "red3", limits = c(0, 20)) +
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
  xlab("Stripe Intensity Threshold") +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-Loop Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.15, 0.15) +
  scale_x_discrete(labels = c("> 0", "", "", "> 30", "", "", "> 60", "", "", "> 90", "", "", "> 120", "", ""))
