###
### RNA-seq analysis of E12.5 Distal vs Proximal Forelimb
### https://doi.org/10.1101/gad.348934.121

library(DESeq2)
library(sva)
library(tximport)
library(readr)
library(IHW)
library(AnnotationDbi)
library(EnhancedVolcano)
library(Seurat)

setwd("./")

files <- file.path(list.dirs("salmon_distal_proximal", recursive = F), "quant.sf")
names(files) <- list.dirs("salmon_distal_proximal", recursive = F, full.names = F)
samples <- read_csv("sampledata_distal_proximal.csv")
tx2gene <- read_tsv("salmon_tx2gene.tsv", col_names = F)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi_combat <- ComBat_seq(txi$counts, batch = samples$batch)

##### ComBat Batch Correction
dds <- DESeqDataSetFromMatrix(round(txi_combat), colData = samples, design = ~condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
dds$condition <- factor(dds$condition, levels = c("DFL", "PFL"))

dds_norm <- rlog(dds, blind = T)
mat <- assay(dds_norm)
mm <- model.matrix(~condition, colData(dds_norm))
mat <- limma::removeBatchEffect(mat, batch = dds_norm$batch, design = mm)
assay(dds_norm) <- mat

##### DESeq2 Processing
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "DFL", "PFL"))

ens.str <- rownames(res)
res$symbol <- mapIds(org.Mm.eg.db,
  keys = ens.str,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res$entrez <- mapIds(org.Mm.eg.db,
  keys = ens.str,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res$refseq <- mapIds(org.Mm.eg.db,
  keys = ens.str,
  column = "REFSEQ",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res$foldchange <- (2^res$log2FoldChange)

res_clamp <- as.data.frame(res)
res_clamp <- res_clamp[res_clamp$padj < 0.05, ]
res_clamp <- na.omit(res_clamp)
res_clamp <- res_clamp[order(-res_clamp$foldchange), ]
write.csv(res_clamp, "E12_distal_proximal_DEG.csv")


deseq_counts <- as.data.frame(counts(dds, normalized = TRUE))
deseq_counts$ensembl <- rownames(deseq_counts)
deseq_counts$symbol <- unlist(mapIds(org.Mm.eg.db,
  keys = deseq_counts$ensembl,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
))

deseq_counts <- deseq_counts[complete.cases(deseq_counts), ]
deseq_counts$ensembl <- NULL
deseq_counts <- aggregate(. ~ symbol, deseq_counts, sum)
rownames(deseq_counts) <- deseq_counts$symbol
write.csv(deseq_counts, "E12_distal_proximal_counts.csv")

##### Plotting
res <- lfcShrink(dds, contrast = c("condition", "DFL", "PFL"), res = res, type = "normal")
res$symbol <- mapIds(org.Mm.eg.db,
  keys = ens.str,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

keyvals <- ifelse(
  res$log2FoldChange < -log2(1.25) & res$padj < 0.05, "#C51B7DFF",
  ifelse(res$log2FoldChange > log2(1.25) & res$padj < 0.05, "#4D9221FF",
    "grey"
  )
)
keyvals[is.na(keyvals)] <- "grey"
names(keyvals)[keyvals == "#4D9221FF"] <- "high"
names(keyvals)[keyvals == "grey"] <- "mid"
names(keyvals)[keyvals == "#C51B7DFF"] <- "low"

EnhancedVolcano(res,
  lab = NA,
  x = "log2FoldChange",
  y = "padj",
  title = "",
  subtitle = "",
  caption = "",
  pCutoff = 0.05,
  FCcutoff = log2(1.25),
  pointSize = 0.5,
  labSize = 4.0,
  labCol = "black",
  axisLabSize = 11,
  colAlpha = 1,
  colCustom = keyvals,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  xlim = c(-2.5, 2.5),
  ylim = c(0, 290),
  boxedLabels = T,
  drawConnectors = TRUE,
  arrowheads = F
) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11)
  ) +
  xlab(expression("log"["2"] ~ "Fold Change")) +
  ylab(expression("log"["10"] ~ italic(P))) +
  NoLegend()
