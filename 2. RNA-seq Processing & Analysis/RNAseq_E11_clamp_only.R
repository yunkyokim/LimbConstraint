###
### RNA-seq analysis for clamp vs unclamped only, in addition to supplementary plotting
###

library(DESeq2)
library(sva)
library(tximport)
library(readr)
library(IHW)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggplot2)
library(Seurat)
library(tidyr)
library(dplyr)
library(clusterProfiler)
library(ggpubr)

setwd("./")
dir <- getwd()

files <- file.path(list.dirs("salmon", recursive = F), "quant.sf")
names(files) <- list.dirs("salmon", recursive = F, full.names = F)
samples <- read_csv("sampledata.csv")
tx2gene <- read_tsv("salmon_tx2gene.tsv", col_names = F)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi_combat <- ComBat_seq(txi$counts, batch = samples$batch)

dds <- DESeqDataSetFromMatrix(round(txi_combat), colData = samples, design = ~condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
dds$condition <- factor(dds$condition, levels = c("DISTCOM", "DISTUNCOM", "AMNION"))

dds_norm <- rlog(dds, blind = F)
mat <- assay(dds_norm)
mm <- model.matrix(~condition, colData(dds_norm))
mat <- limma::removeBatchEffect(mat, batch = dds_norm$batch, design = mm)
assay(dds_norm) <- mat
dds <- DESeq(dds)

##### Clamp vs Control
res <- results(dds, contrast = c("condition", "DISTCOM", "DISTUNCOM"))

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
res_clamp <- res_clamp[abs(res_clamp$log2FoldChange) > log2(1.25), ]
res_clamp <- na.omit(res_clamp)
res_clamp <- res_clamp[order(-res_clamp$foldchange), ]
write.csv(res_clamp, "res_clamp.csv")

##### Volcano Plotting
res <- lfcShrink(dds, contrast = c("condition", "DISTCOM", "DISTUNCOM"), res = res, type = "normal")
ens.str <- rownames(res)
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


##### ChIP Peak Overlap Analysis
res_clamp <- read.csv("res_clamp.csv")
res_const <- read.csv("res_constraint.csv")

diff_clamp_up <- res_clamp[res_clamp$log2FoldChange > log2(1.25), "symbol"]
diff_clamp_down <- res_clamp[res_clamp$log2FoldChange < -log2(1.25), "symbol"]
diff_const_up <- res_const[res_const$log2FoldChange > log2(1.25), "symbol"]
diff_const_down <- res_const[res_const$log2FoldChange < -log2(1.25), "symbol"]

files <- file.path(list.files("chip_peaks", recursive = F))
names <- sapply(str_split(files, "_.P"), `[`, 1)
files <- paste("chip_peaks/", files, sep = "")

chip_markers <- data.frame()
for (i in 1:length(files)) {
  chip_peaks <- read_tsv(files[i])
  chip_peaks <- chip_peaks[grepl("TSS", chip_peaks[["Annotation"]]), ]

  chip_peaks_clamp_UP <- chip_peaks[chip_peaks$`Gene Name` %in% diff_clamp_up, ]
  chip_peaks_clamp_DOWN <- chip_peaks[chip_peaks$`Gene Name` %in% diff_clamp_down, ]
  chip_peaks_const_UP <- chip_peaks[chip_peaks$`Gene Name` %in% diff_const_up, ]
  chip_peaks_const_DOWN <- chip_peaks[chip_peaks$`Gene Name` %in% diff_const_down, ]

  chip_markers[i, "name"] <- names[i]
  chip_markers[i, "clamp_up_score"] <- mean(chip_peaks_clamp_UP$`Peak Score`)
  chip_markers[i, "clamp_down_score"] <- mean(chip_peaks_clamp_DOWN$`Peak Score`)
  chip_markers[i, "const_up_score"] <- mean(chip_peaks_const_UP$`Peak Score`)
  chip_markers[i, "const_down_score"] <- mean(chip_peaks_const_DOWN$`Peak Score`)

  chip_markers[i, "clamp_up_prop"] <- length(unique(chip_peaks_clamp_UP$`Gene Name`)) / length(diff_clamp_up)
  chip_markers[i, "clamp_down_prop"] <- length(unique(chip_peaks_clamp_DOWN$`Gene Name`)) / length(diff_clamp_down)
  chip_markers[i, "const_up_prop"] <- length(unique(chip_peaks_const_UP$`Gene Name`)) / length(diff_const_up)
  chip_markers[i, "const_down_prop"] <- length(unique(chip_peaks_const_DOWN$`Gene Name`)) / length(diff_const_down)
}

chip_markers <- chip_markers[(chip_markers$name != "E11 H3K9me"), ]
chip_markers_long <- data.frame(name = append(chip_markers$name, c(chip_markers$name, chip_markers$name, chip_markers$name)))
chip_markers_long$score <- append(chip_markers$clamp_up_score, c(chip_markers$clamp_down_score, chip_markers$const_up_score, chip_markers$const_down_score))
chip_markers_long$prop <- append(chip_markers$clamp_up_prop, c(chip_markers$clamp_down_prop, chip_markers$const_up_prop, chip_markers$const_down_prop))
chip_markers_long$sample <- c(rep("Clamp DE UP", 13), rep("Clamp DE DOWN", 13), rep("Constraint DE UP", 13), rep("Constraint DE DOWN", 13))


heat_colors <- (colorRampPalette(c("grey", "red"))(50))

##### Plot ChIP-seq Dotplot
sampleorder <- rev(c("Clamp DE UP", "Clamp DE DOWN", "Constraint DE UP", "Constraint DE DOWN"))
ggplot(chip_markers_long, aes(x = name, y = factor(sample, level = sampleorder), color = score, size = prop)) +
  geom_point() +
  scale_color_gradientn(colors = heat_colors) +
  scale_radius(range = 15 * range(chip_markers_long$prop)) +
  theme_classic() +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    text = element_text(size = 11), axis.text = element_text(size = 12),
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 11, colour = "black", angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 11, colour = "black"), legend.box = "horizontal"
  ) +
  labs(title = "", y = "", x = "", color = "Median Peak \nScore", size = "Proportion")


##### E11 Spatial Signatures
long <- readRDS("long_spatial.rds")
long$value <- as.character(long$value)

res_clamp <- read.csv("res_clamp.csv")
genelist <- res_clamp[, "log2FoldChange"]
names(genelist) <- as.character(res_clamp[, "symbol"])
genelist <- sort(genelist, decreasing = TRUE)
x <- GSEA(genelist, TERM2GENE = long, pvalueCutoff = 1, minGSSize = 0)
write.csv(x, "GSEA_clamp.csv")

res_const <- read.csv("res_constraint.csv")
res_const <- res_const[abs(res_const$log2FoldChange) > log2(1.25), ]
write.csv(res_const, "res_const.csv")
genelist <- res_const[, "log2FoldChange"]
names(genelist) <- as.character(res_const[, "symbol"])
genelist <- sort(genelist, decreasing = TRUE)
x <- GSEA(genelist, TERM2GENE = long, pvalueCutoff = 1, minGSSize = 0)
write.csv(x, "GSEA_const.csv")

GSEA_clamp <- as.data.frame(read.csv("GSEA_clamp.csv"))
GSEA_const <- as.data.frame(read.csv("GSEA_const.csv"))

##### Clamp Plotting
for (i in 1:nrow(GSEA_clamp)) {
  if (GSEA_clamp[i, 6] > 0) {
    GSEA_clamp[i, "fill"] <- "#4D9221FF"
  } else {
    GSEA_clamp[i, "fill"] <- "#C51B7DFF"
  }
}
GSEA_clamp$Description <- gsub("_markers", "", GSEA_clamp$Description)
GSEA_clamp$Description <- gsub("_", " ", GSEA_clamp$Description)
ggbarplot(GSEA_clamp,
  x = "Description", y = "NES",
  fill = "fill",
  color = "fill",
  palette = (unique(GSEA_clamp$fill)),
  sort.val = "asc",
  sort.by.groups = FALSE,
  ylab = "NES",
  xlab = "",
  ggtitle = "Distal Clamped vs Unclamped",
  rotate = TRUE,
  ggtheme = theme_classic() +
    NoLegend() +
    theme(
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black")
    )
)

##### Constraint Plotting
for (i in 1:nrow(GSEA_const)) {
  if (GSEA_const[i, "NES"] > 0) {
    GSEA_const[i, "fill"] <- "#4D9221FF"
  } else {
    GSEA_const[i, "fill"] <- "#C51B7DFF"
  }
}
GSEA_const$Description <- gsub("_markers", "", GSEA_const$Description)
GSEA_const$Description <- gsub("_", " ", GSEA_const$Description)
ggbarplot(GSEA_const,
  x = "Description", y = "NES",
  fill = "fill",
  color = "white",
  palette = c("#4D9221FF", "#C51B7DFF"),
  sort.val = "asc",
  sort.by.groups = FALSE,
  x.text.angle = 90,
  ylab = "NES",
  xlab = "",
  legend.title = "Amnion ON vs Unclamped",
  rotate = TRUE,
  ggtheme = theme_classic() +
    NoLegend() +
    theme(
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black")
    )
)


##### Mechanotransduction Signatures
sig_files <- list.files("mech_sigs") # mSigDB Signatures
sig_data <- data.frame()
for (i in 1:length(sig_files)) {
  gmt <- as.data.frame(read_tsv(paste("mech_sigs/", sig_files[i], sep = "")))
  gmt <- colnames(gmt[, 3:ncol(gmt)])
  gmt <- data.frame(sig = rep(sig_files[i], length(gmt)), genes = gmt)
  sig_data <- rbind(sig_data, gmt)
}

##### Clamp Mech
res_clamp <- read.csv("res_clamp.csv")
genelist <- res_clamp[, "log2FoldChange"]
names(genelist) <- as.character(res_clamp[, "symbol"])
genelist <- sort(genelist, decreasing = TRUE)

GSEA_clamp <- GSEA(genelist, TERM2GENE = sig_data, pvalueCutoff = 1, minGSSize = 0)
GSEA_clamp <- as.data.frame(GSEA_clamp)
for (i in 1:nrow(GSEA_clamp)) {
  if (GSEA_clamp[i, 5] > 0) {
    GSEA_clamp[i, "fill"] <- "#4D9221FF"
  } else {
    GSEA_clamp[i, "fill"] <- "#C51B7DFF"
  }
}

GSEA_clamp$Description <- sapply(strsplit(GSEA_clamp$Description, ".v"), `[`, 1)
GSEA_clamp$Description <- gsub("_", " ", GSEA_clamp$Description)
ggbarplot(GSEA_clamp,
  x = "Description", y = "NES",
  fill = "fill",
  color = "fill",
  palette = (unique(GSEA_clamp$fill)),
  sort.val = "asc",
  sort.by.groups = FALSE,
  ylab = "NES",
  xlab = "",
  ggtitle = "Distal Clamped vs Unclamped",
  rotate = TRUE,
  ggtheme = theme_classic() +
    NoLegend() +
    theme(
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black")
    )
) + ylim(-2.5, 2.5)
GSEA_clamp[GSEA_clamp$p.adjust >= 0.05, "Description"]


##### Constraint Mech
res_clamp <- read.csv("res_const.csv")
genelist <- res_clamp[, "log2FoldChange"]
names(genelist) <- as.character(res_clamp[, "symbol"])
genelist <- sort(genelist, decreasing = TRUE)

GSEA_clamp <- GSEA(genelist, TERM2GENE = sig_data, pvalueCutoff = 1, minGSSize = 0)
GSEA_clamp <- as.data.frame(GSEA_clamp)
for (i in 1:nrow(GSEA_clamp)) {
  if (GSEA_clamp[i, 5] > 0) {
    GSEA_clamp[i, "fill"] <- "#4D9221FF"
  } else {
    GSEA_clamp[i, "fill"] <- "#C51B7DFF"
  }
}

GSEA_clamp$Description <- sapply(strsplit(GSEA_clamp$Description, ".v"), `[`, 1)
GSEA_clamp$Description <- gsub("_", " ", GSEA_clamp$Description)
ggbarplot(GSEA_clamp,
  x = "Description", y = "NES",
  fill = "fill",
  color = "fill",
  palette = (unique(GSEA_clamp$fill)),
  sort.val = "asc",
  sort.by.groups = FALSE,
  ylab = "NES",
  xlab = "",
  ggtitle = "Distal Clamped vs Unclamped",
  rotate = TRUE,
  ggtheme = theme_classic() +
    NoLegend() +
    theme(
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black")
    )
) + ylim(-2.5, 2.5)
GSEA_clamp[GSEA_clamp$p.adjust >= 0.05, "Description"]


##### E9.5-E13.5 Organogenesis Atlas Signatures
traj <- read.csv("DE_gene_sub_trajectory.csv")
traj <- traj[, c("class", "gene_short_name")]

E11trajgene <- traj %>%
  group_by(class) %>%
  mutate(row = row_number()) %>%
  spread(class, gene_short_name) %>%
  select(-row)
E11trajgene <- as.list(E11trajgene)
E11trajgene <- na.omit(E11trajgene)

# Clamp
res_clamp <- read.csv("res_clamp.csv")
clamplist <- res_clamp$log2FoldChange
names(clamplist) <- res_clamp$symbol

clamp_fsea <- GSEA(clamplist, TERM2GENE = traj, pvalueCutoff = 1, minGSSize = 0)
write.csv(clamp_fsea, "GSEA_clamp_atlas.csv")
clamp_fsea <- as.data.frame(read.csv("GSEA_clamp_atlas.csv"))

for (i in 1:nrow(clamp_fsea)) {
  if (clamp_fsea[i, "NES"] > 0) {
    clamp_fsea[i, "fill"] <- "#4D9221FF"
  } else {
    clamp_fsea[i, "fill"] <- "#C51B7DFF"
  }
}
clamp_fsea2 <- clamp_fsea[clamp_fsea$p.adjust < 0.05, ]
ggbarplot(clamp_fsea2,
  x = "Description", y = "NES",
  fill = "fill",
  color = "white",
  palette = c("#4D9221FF", "#C51B7DFF"),
  sort.val = "asc",
  sort.by.groups = FALSE,
  x.text.angle = 90,
  ylab = "NES",
  xlab = "",
  legend.title = "Clamped vs Unclamped",
  rotate = TRUE,
  ggtheme = theme_classic() +
    NoLegend() +
    theme(
      aspect.ratio = 1.2,
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black")
    )
)

# Constraint
res_const <- read.csv("res_const.csv")
res_const <- res_const[-c(1, 1365, 1364, 1363), ]
clamplist <- res_const$log2FoldChange
names(clamplist) <- res_const$symbol

const_fsea <- GSEA(clamplist, TERM2GENE = traj, pvalueCutoff = 1, minGSSize = 0)
write.csv(const_fsea, "GSEA_const_atlas.csv")
const_fsea <- as.data.frame(read.csv("GSEA_const_atlas.csv"))

for (i in 1:nrow(const_fsea)) {
  if (const_fsea[i, "NES"] > 0) {
    const_fsea[i, "fill"] <- "#4D9221FF"
  } else {
    const_fsea[i, "fill"] <- "#C51B7DFF"
  }
}
const_fsea2 <- const_fsea[const_fsea$p.adjust < 0.05, ]
ggbarplot(const_fsea2,
  x = "Description", y = "NES",
  fill = "fill",
  color = "white",
  palette = c("#4D9221FF", "#C51B7DFF"),
  sort.val = "asc",
  sort.by.groups = FALSE,
  x.text.angle = 90,
  ylab = "NES",
  xlab = "",
  legend.title = "Constrained vs Unclamped",
  rotate = TRUE,
  ggtheme = theme_classic() +
    NoLegend() +
    theme(
      aspect.ratio = 1.2,
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black")
    )
)


##### ECDF Analysis
##### CLamp ECDF
deseq_counts <- read.csv("deseq_counts.csv")
res_clamp <- read.csv("res_clamp.csv")

deseq_counts$amnion <- log2(rowMeans(deseq_counts[, 3:4]))
deseq_counts$clamp <- log2(rowMeans(deseq_counts[, 5:6]))
deseq_counts$control <- log2(rowMeans(deseq_counts[, 7:8]))
deseq_counts <- na.omit(deseq_counts)
deseq_counts <- deseq_counts[is.finite(deseq_counts$amnion), ]
deseq_counts <- deseq_counts[is.finite(deseq_counts$control), ]
deseq_counts <- deseq_counts[is.finite(deseq_counts$clamp), ]

ggplot() +
  stat_ecdf(data = deseq_counts[!(deseq_counts$symbol %in% res_clamp[, "symbol"]), ], aes(control), geom = "step", size = 1, color = "grey20") +
  stat_ecdf(data = deseq_counts[deseq_counts$symbol %in% res_clamp[res_clamp$log2FoldChange > log2(1.25), "symbol"], ], aes(control), geom = "step", size = 1, color = "#4D9221FF") +
  stat_ecdf(data = deseq_counts[deseq_counts$symbol %in% res_clamp[res_clamp$log2FoldChange < -(log2(1.25)), "symbol"], ], aes(control), geom = "step", size = 1, color = "#C51B7DFF") +
  theme_classic() +
  theme(
    aspect.ratio = 1, plot.title = element_blank(),
    text = element_text(size = 11), axis.text = element_text(size = 11),
    plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Cumulative Frequency", x = expression("log"["2"] ~ "TMM")) +
  NoLegend() # clamp_control_ecdf
ggplot() +
  stat_ecdf(data = deseq_counts[!(deseq_counts$symbol %in% res_clamp[, "symbol"]), ], aes(amnion), geom = "step", size = 1, color = "grey20", pad = TRUE) +
  stat_ecdf(data = deseq_counts[deseq_counts$symbol %in% res_clamp[res_clamp$log2FoldChange > log2(1.25), "symbol"], ], aes(amnion), geom = "step", size = 1, color = "#4D9221FF", pad = TRUE) +
  stat_ecdf(data = deseq_counts[deseq_counts$symbol %in% res_clamp[res_clamp$log2FoldChange < -(log2(1.25)), "symbol"], ], aes(amnion), geom = "step", size = 1, color = "#C51B7DFF", pad = TRUE) +
  theme_classic() +
  theme(
    aspect.ratio = 1, plot.title = element_blank(),
    text = element_text(size = 11), axis.text = element_text(size = 11),
    plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Cumulative Frequency", x = expression("log"["2"] ~ "TMM")) +
  NoLegend() # clamp_amnion_ecdf


mean(deseq_counts[deseq_counts$symbol %in% res_clamp[res_clamp$log2FoldChange > log2(1.25), "symbol"], "amnion"])
mean(deseq_counts[deseq_counts$symbol %in% res_clamp[res_clamp$log2FoldChange < -(log2(1.25)), "symbol"], "amnion"])
mean(deseq_counts[!(deseq_counts$symbol %in% res_clamp[, "symbol"]), "amnion"])

mean(deseq_counts[deseq_counts$symbol %in% res_clamp[res_clamp$log2FoldChange > log2(1.25), "symbol"], "control"])
mean(deseq_counts[deseq_counts$symbol %in% res_clamp[res_clamp$log2FoldChange < -(log2(1.25)), "symbol"], "control"])
mean(deseq_counts[!(deseq_counts$symbol %in% res_clamp[, "symbol"]), "control"])

##### Constraint ECDF
deseq_counts <- read.csv("deseq_counts.csv")
res_const <- read.csv("res_constraint.csv")
res_const <- res_const[abs(res_const$log2FoldChange) > log2(1.25), ]

deseq_counts$amnion <- log2(rowMeans(deseq_counts[, 3:4]))
deseq_counts$clamp <- log2(rowMeans(deseq_counts[, 5:6]))
deseq_counts$control <- log2(rowMeans(deseq_counts[, 7:8]))
deseq_counts <- na.omit(deseq_counts)
deseq_counts <- deseq_counts[is.finite(deseq_counts$amnion), ]
deseq_counts <- deseq_counts[is.finite(deseq_counts$control), ]
deseq_counts <- deseq_counts[is.finite(deseq_counts$clamp), ]

ggplot() +
  stat_ecdf(data = deseq_counts[!(deseq_counts$symbol %in% res_const[, "symbol"]), ], aes(control), geom = "step", size = 1, color = "grey20") +
  stat_ecdf(data = deseq_counts[deseq_counts$symbol %in% res_const[res_const$log2FoldChange > log2(1.25), "symbol"], ], aes(control), geom = "step", size = 1, color = "#4D9221FF") +
  stat_ecdf(data = deseq_counts[deseq_counts$symbol %in% res_const[res_const$log2FoldChange < -(log2(1.25)), "symbol"], ], aes(control), geom = "step", size = 1, color = "#C51B7DFF") +
  theme_classic() +
  theme(
    aspect.ratio = 1, plot.title = element_blank(),
    text = element_text(size = 11), axis.text = element_text(size = 11),
    plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Cumulative Frequency", x = expression("log"["2"] ~ "TMM")) +
  NoLegend() # const_control_ecdf
ggplot() +
  stat_ecdf(data = deseq_counts[!(deseq_counts$symbol %in% res_const[, "symbol"]), ], aes(amnion), geom = "step", size = 1, color = "grey20", pad = TRUE) +
  stat_ecdf(data = deseq_counts[deseq_counts$symbol %in% res_const[res_const$log2FoldChange > log2(1.25), "symbol"], ], aes(amnion), geom = "step", size = 1, color = "#4D9221FF", pad = TRUE) +
  stat_ecdf(data = deseq_counts[deseq_counts$symbol %in% res_const[res_const$log2FoldChange < -(log2(1.25)), "symbol"], ], aes(amnion), geom = "step", size = 1, color = "#C51B7DFF", pad = TRUE) +
  theme_classic() +
  theme(
    aspect.ratio = 1, plot.title = element_blank(),
    text = element_text(size = 11), axis.text = element_text(size = 11),
    plot.margin = margin(t = 7, r = 0, b = 7, l = 0, unit = "pt"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(title = "", y = "Cumulative Frequency", x = expression("log"["2"] ~ "TMM")) +
  NoLegend() # const_amnion_ecdf
