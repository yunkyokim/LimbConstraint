###
### Volcano analysis of all HicCompare Interactions
###


library(rtracklayer)
library(Seurat)
library(dplyr)
library(ggplot2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(plotgardener)
library(EnhancedVolcano)
library(respR)

options(scipen = 999)
setwd("./")
dir <- getwd()

##### Clamped vs Unclamped
comp_hictable <- read.table("E11_microc_hictable_merge_500000.txt", header = T)
hictable_volcano <- data.frame(
  "labels" = rownames(comp_hictable),
  "log2FC" = log2(as.numeric(comp_hictable$adj.IF1) / as.numeric(comp_hictable$adj.IF2)), # clamped over unclamped
  "padj" = as.numeric(comp_hictable$p.adj)
)

hictable_volcano <- na.omit(hictable_volcano)
hictable_volcano <- hictable_volcano[hictable_volcano$padj < 0.99, ]
nrow(hictable_volcano[abs(hictable_volcano$log2FC) > log2(1.25) & hictable_volcano$padj < 0.05, ])

hictable_volcano <- subsample(hictable_volcano, n = 2, plot = F)

keyvals <- ifelse(
  hictable_volcano$log2FC < -log2(1.25) & hictable_volcano$padj < 0.05, "#C51B7DFF",
  ifelse(hictable_volcano$log2FC > log2(1.25) & hictable_volcano$padj < 0.05, "#4D9221FF",
    "grey"
  )
)
keyvals[is.na(keyvals)] <- "grey"
names(keyvals)[keyvals == "#4D9221FF"] <- "high"
names(keyvals)[keyvals == "grey"] <- "mid"
names(keyvals)[keyvals == "#C51B7DFF"] <- "low"

EnhancedVolcano(hictable_volcano,
  lab = NA,
  x = "log2FC",
  y = "padj",
  title = "",
  subtitle = "",
  caption = "",
  pCutoff = 0.05,
  FCcutoff = log2(1.25),
  pointSize = 0.1,
  labSize = 4.0,
  labCol = "black",
  axisLabSize = 11,
  colAlpha = 1,
  colCustom = keyvals,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  xlim = c(-4, 4),
  ylim = c(0, 10),
  boxedLabels = T,
  drawConnectors = TRUE,
  arrowheads = F
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11)
  ) +
  xlab(expression("-log"["2"] ~ "Fold Change")) +
  ylab(expression("log"["10"] ~ italic(P))) +
  NoLegend()

##### Distal vs Proximal
comp_hictable <- read.table("E12_hic_hictable_merge_500000.txt", header = T)
hictable_volcano <- data.frame(
  "labels" = rownames(comp_hictable),
  "log2FC" = log2(as.numeric(comp_hictable$adj.IF1) / as.numeric(comp_hictable$adj.IF2)), # distal over proximal
  "padj" = as.numeric(comp_hictable$p.adj)
)

hictable_volcano <- na.omit(hictable_volcano)
hictable_volcano <- hictable_volcano[hictable_volcano$padj < 0.99, ]
nrow(hictable_volcano[abs(hictable_volcano$log2FC) > log2(1.25) & hictable_volcano$padj < 0.05, ])

keyvals <- ifelse(
  hictable_volcano$log2FC < -log2(1.25) & hictable_volcano$padj < 0.05, "#C51B7DFF",
  ifelse(hictable_volcano$log2FC > log2(1.25) & hictable_volcano$padj < 0.05, "#4D9221FF",
    "grey"
  )
)
keyvals[is.na(keyvals)] <- "grey"
names(keyvals)[keyvals == "#4D9221FF"] <- "high"
names(keyvals)[keyvals == "grey"] <- "mid"
names(keyvals)[keyvals == "#C51B7DFF"] <- "low"

EnhancedVolcano(hictable_volcano,
  lab = NA,
  x = "log2FC",
  y = "padj",
  title = "",
  subtitle = "",
  caption = "",
  pCutoff = 0.05,
  FCcutoff = log2(1.25),
  pointSize = 0.1,
  labSize = 4.0,
  labCol = "black",
  axisLabSize = 11,
  colAlpha = 1,
  colCustom = keyvals,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  xlim = c(-4, 4),
  ylim = c(0, 10),
  boxedLabels = T,
  drawConnectors = TRUE,
  arrowheads = F
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11)
  ) +
  xlab(expression("-log"["2"] ~ "Fold Change")) +
  ylab(expression("log"["10"] ~ italic(P))) +
  NoLegend()
