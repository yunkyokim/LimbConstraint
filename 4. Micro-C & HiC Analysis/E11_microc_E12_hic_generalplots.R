###
### General plotting of contact matrices
###

library(rtracklayer)
library(GenomicRanges)
library(HiCcompare)
library(AnnotationDbi)
library(plotgardener)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(multiHiCcompare)

options(scipen = 999)
setwd("./")

mm39 <- assembly(
  Genome = "mm39",
  TxDb = TxDb.Mmusculus.UCSC.mm39.refGene,
  OrgDb = "org.Mm.eg.db",
  BSgenome = "BSgenome.Mmusculus.UCSC.mm39"
)

##### HoxD Plotting
i <- 5000
clamped_plot <- read.table(paste("hiccompare/clamped_hiccompare_", format(i, scientific = F), "_2.txt", sep = ""))
unclamped_plot <- read.table(paste("hiccompare/unclamped_hiccompare_", format(i, scientific = F), "_2.txt", sep = ""))

start <- 73646853 # Atf2
end <- 75499751 # Hnrnpa3
trange <- 10
pageCreate(width = 16, height = 10, default.units = "cm")
hicPlot <- plotHicSquare(
  data = unclamped_plot,
  zrange = c(0, trange),
  bg = "white",
  chrom = "2", chromstart = start, chromend = end,
  assembly = mm39,
  palette = colorRampPalette(colors = c("white", "#FF2345", "violetred4")),
  colorTrans = "linear",
  x = 0, y = 0, width = 6, height = 6, default.units = "cm"
)
annoGenomeLabel(
  plot = hicPlot,
  fontsize = 8,
  x = 0, y = 6.2, width = 6, height = 2, default.units = "cm"
)
plotGenes(
  assembly = mm39,
  chrom = "chr2", chromstart = start, chromend = end,
  x = 0, y = 6.8, width = 6, height = 1.5, default.units = "cm",
  fill = c("black", "black"),
  fontsize = 8,
  fontcolor = c("darkgreen", "dodgerblue4"),
  stroke = 0.1
)

## 2
hicPlot <- plotHicSquare(
  data = clamped_plot,
  bg = "white",
  zrange = c(0, trange),
  chrom = "2", chromstart = start, chromend = end,
  assembly = mm39,
  palette = colorRampPalette(colors = c("white", "#FF2345", "violetred4")),
  colorTrans = "linear",
  x = 7, y = 0, width = 6, height = 6, default.units = "cm"
)
annoGenomeLabel(
  plot = hicPlot,
  fontsize = 8,
  x = 7, y = 6.2, width = 6, height = 2, default.units = "cm"
)
plotGenes(
  assembly = mm39,
  chrom = "chr2", chromstart = start, chromend = end,
  x = 7, y = 6.8, width = 6, height = 1.5, default.units = "cm",
  fill = c("black", "black"),
  fontsize = 8,
  fontcolor = c("darkgreen", "dodgerblue4"),
  stroke = 0.1
)
annoHeatmapLegend(
  plot = hicPlot,
  x = 13.5, y = 0.6, width = 0.5, height = 2, default.units = "cm",
  fontcolor = "black"
)
pageGuideHide() # comp_hoxd5k, 700x400

i <- 5000
distal_plot <- read.table(paste("hiccompare/distal_hiccompare_", format(i, scientific = F), "_2.txt", sep = ""))
proximal_plot <- read.table(paste("hiccompare/proximal_hiccompare_", format(i, scientific = F), "_2.txt", sep = ""))

trange <- 10
pageCreate(width = 16, height = 10, default.units = "cm")
hicPlot <- plotHicSquare(
  data = proximal_plot,
  zrange = c(0, trange),
  bg = "white",
  chrom = "2", chromstart = start, chromend = end,
  assembly = mm39,
  palette = colorRampPalette(colors = c("white", "#FF2345", "violetred4")),
  colorTrans = "linear",
  x = 0, y = 0, width = 6, height = 6, default.units = "cm"
)
annoGenomeLabel(
  plot = hicPlot,
  fontsize = 8,
  x = 0, y = 6.2, width = 6, height = 2, default.units = "cm"
)
plotGenes(
  assembly = mm39,
  chrom = "chr2", chromstart = start, chromend = end,
  x = 0, y = 6.8, width = 6, height = 1.5, default.units = "cm",
  fill = c("black", "black"),
  fontsize = 8,
  fontcolor = c("darkgreen", "dodgerblue4"),
  stroke = 0.1
)

## 2
hicPlot <- plotHicSquare(
  data = distal_plot,
  bg = "white",
  zrange = c(0, trange),
  chrom = "2", chromstart = start, chromend = end,
  assembly = mm39,
  palette = colorRampPalette(colors = c("white", "#FF2345", "violetred4")),
  colorTrans = "linear",
  x = 7, y = 0, width = 6, height = 6, default.units = "cm"
)
annoGenomeLabel(
  plot = hicPlot,
  fontsize = 8,
  x = 7, y = 6.2, width = 6, height = 2, default.units = "cm"
)
plotGenes(
  assembly = mm39,
  chrom = "chr2", chromstart = start, chromend = end,
  x = 7, y = 6.8, width = 6, height = 1.5, default.units = "cm",
  fill = c("black", "black"),
  fontsize = 8,
  fontcolor = c("darkgreen", "dodgerblue4"),
  stroke = 0.1
)
annoHeatmapLegend(
  plot = hicPlot,
  x = 13.5, y = 0.6, width = 0.5, height = 2, default.units = "cm",
  fontcolor = "black"
)
pageGuideHide() # 12_hoxd5k, 700x400


##### HoxA Plotting
i <- 5000
clamped_plot <- read.table(paste("hiccompare/clamped_hiccompare_", format(i, scientific = F), "_6.txt", sep = ""))
unclamped_plot <- read.table(paste("hiccompare/unclamped_hiccompare_", format(i, scientific = F), "_6.txt", sep = ""))

start <- 50539543 # Atf2
end <- 53055616 # Cycs
trange <- 10
pageCreate(width = 16, height = 10, default.units = "cm")
hicPlot <- plotHicSquare(
  data = unclamped_plot,
  zrange = c(0, trange),
  bg = "white",
  chrom = "6", chromstart = start, chromend = end,
  assembly = mm39,
  palette = colorRampPalette(colors = c("white", "#FF2345", "violetred4")),
  colorTrans = "linear",
  x = 0, y = 0, width = 6, height = 6, default.units = "cm"
)
annoGenomeLabel(
  plot = hicPlot,
  fontsize = 8,
  x = 0, y = 6.2, width = 6, height = 2, default.units = "cm"
)
plotGenes(
  assembly = mm39,
  chrom = "chr6", chromstart = start, chromend = end,
  x = 0, y = 6.8, width = 6, height = 1.5, default.units = "cm",
  fill = c("black", "black"),
  fontsize = 8,
  fontcolor = c("darkgreen", "dodgerblue4"),
  stroke = 0.1
)

## 2
hicPlot <- plotHicSquare(
  data = clamped_plot,
  bg = "white",
  zrange = c(0, trange),
  chrom = "6", chromstart = start, chromend = end,
  assembly = mm39,
  palette = colorRampPalette(colors = c("white", "#FF2345", "violetred4")),
  colorTrans = "linear",
  x = 7, y = 0, width = 6, height = 6, default.units = "cm"
)
annoGenomeLabel(
  plot = hicPlot,
  fontsize = 8,
  x = 7, y = 6.2, width = 6, height = 2, default.units = "cm"
)
plotGenes(
  assembly = mm39,
  chrom = "chr6", chromstart = start, chromend = end,
  x = 7, y = 6.8, width = 6, height = 1.5, default.units = "cm",
  fill = c("black", "black"),
  fontsize = 8,
  fontcolor = c("darkgreen", "dodgerblue4"),
  stroke = 0.1
)
annoHeatmapLegend(
  plot = hicPlot,
  x = 13.5, y = 0.6, width = 0.5, height = 2, default.units = "cm",
  fontcolor = "black"
)
pageGuideHide() # comp_hoxa5k, 700x400

i <- 5000
distal_plot <- read.table(paste("hiccompare/distal_hiccompare_", format(i, scientific = F), "_6.txt", sep = ""))
proximal_plot <- read.table(paste("hiccompare/proximal_hiccompare_", format(i, scientific = F), "_6.txt", sep = ""))

trange <- 10
pageCreate(width = 16, height = 10, default.units = "cm")
hicPlot <- plotHicSquare(
  data = proximal_plot,
  zrange = c(0, trange),
  bg = "white",
  chrom = "6", chromstart = start, chromend = end,
  assembly = mm39,
  palette = colorRampPalette(colors = c("white", "#FF2345", "violetred4")),
  colorTrans = "linear",
  x = 0, y = 0, width = 6, height = 6, default.units = "cm"
)
annoGenomeLabel(
  plot = hicPlot,
  fontsize = 8,
  x = 0, y = 6.2, width = 6, height = 2, default.units = "cm"
)
plotGenes(
  assembly = mm39,
  chrom = "chr6", chromstart = start, chromend = end,
  x = 0, y = 6.8, width = 6, height = 1.5, default.units = "cm",
  fill = c("black", "black"),
  fontsize = 8,
  fontcolor = c("darkgreen", "dodgerblue4"),
  stroke = 0.1
)

## 2
hicPlot <- plotHicSquare(
  data = distal_plot,
  bg = "white",
  zrange = c(0, trange),
  chrom = "6", chromstart = start, chromend = end,
  assembly = mm39,
  palette = colorRampPalette(colors = c("white", "#FF2345", "violetred4")),
  colorTrans = "linear",
  x = 7, y = 0, width = 6, height = 6, default.units = "cm"
)
annoGenomeLabel(
  plot = hicPlot,
  fontsize = 8,
  x = 7, y = 6.2, width = 6, height = 2, default.units = "cm"
)
plotGenes(
  assembly = mm39,
  chrom = "chr6", chromstart = start, chromend = end,
  x = 7, y = 6.8, width = 6, height = 1.5, default.units = "cm",
  fill = c("black", "black"),
  fontsize = 8,
  fontcolor = c("darkgreen", "dodgerblue4"),
  stroke = 0.1
)
annoHeatmapLegend(
  plot = hicPlot,
  x = 13.5, y = 0.6, width = 0.5, height = 2, default.units = "cm",
  fontcolor = "black"
)
pageGuideHide()
