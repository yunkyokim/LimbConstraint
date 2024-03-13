###
### TAD calling and analysis using GENOVA
###

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(GENOVA)

options(scipen = 999)
setwd("./")
dir <- getwd()


##### Matrix formatting for GENOVA
i <- 50000
bedpe <- read.table("hictable_merge_50000.txt", header = T)
genova_bed <- read.table(paste("matrix/raw/clamped.", i, "_abs.bed", sep = ""))
bedpe$chr1 <- str_replace(bedpe$chr1, "chr", "")
bedpe$chr2 <- str_replace(bedpe$chr2, "chr", "")

bedpe_sub <- data.frame(
  "V1" = genova_bed[match(paste(bedpe$chr1, bedpe$start1), paste(genova_bed$V1, genova_bed$V2)), "V4"],
  "V2" = genova_bed[match(paste(bedpe$chr2, bedpe$start2), paste(genova_bed$V1, genova_bed$V2)), "V4"],
  "V3" = bedpe$IF1
)
bedpe_sub <- bedpe_sub[order(bedpe_sub$V1), ]
bedpe_sub <- na.omit(bedpe_sub)
bedpe_sub$V3 <- formatC(as.numeric(bedpe_sub$V3), format = "f", flag = "0", digits = 10)
write.table(bedpe_sub, paste("genova_input/clamped_", i, "_hictable_genova.matrix", sep = ""), quote = F, col.names = F, row.names = F)

bedpe_sub <- data.frame(
  "V1" = genova_bed[match(paste(bedpe$chr1, bedpe$start1), paste(genova_bed$V1, genova_bed$V2)), "V4"],
  "V2" = genova_bed[match(paste(bedpe$chr2, bedpe$start2), paste(genova_bed$V1, genova_bed$V2)), "V4"],
  "V3" = bedpe$IF2
)
bedpe_sub <- bedpe_sub[order(bedpe_sub$V1), ]
bedpe_sub <- na.omit(bedpe_sub)
bedpe_sub$V3 <- formatC(as.numeric(bedpe_sub$V3), format = "f", flag = "0", digits = 10)
write.table(bedpe_sub, paste("genova_input/unclamped_", i, "_hictable_genova.matrix", sep = ""), quote = F, col.names = F, row.names = F)


##### Load data into GENOVA
clamped <- load_contacts(
  signal_path = "genova_input/clamped_50000_hictable_genova.matrix",
  indices_path = "matrix/raw/clamped.50000_abs.bed",
  sample_name = "Clamped",
  colour = "blue",
  scale_bp = NULL,
  balancing = F,
  verbose = T
)
unclamped <- load_contacts(
  signal_path = "genova_input/unclamped_50000_hictable_genova.matrix",
  indices_path = "matrix/raw/unclamped.50000_abs.bed",
  sample_name = "Control",
  colour = "red",
  scale_bp = NULL,
  balancing = F
)


comp_insulation <- insulation_score(list(unclamped, clamped), window = 25)
TADcalls <- call_TAD_insulation(comp_insulation)

hic_matrixplot(
  exp1 = unclamped,
  exp2 = clamped,
  chrom = "6",
  start = 40000000,
  end = 70000000,
  tads = list(TADcalls$Control, TADcalls$Clamped), # see ATA
  tads.type = list("lower", "upper"), # only plot in lower triangle
  tads.colour = c("green", "grey"), # green TAD-borders
  cut.off = 25
) # upper limit of contacts

ATA_comp <- ATA(
  list(
    "Control" = unclamped,
    "Clamped" = clamped
  ),
  bed = TADcalls$Control
)

visualise(ATA_comp,
  colour_lim = c(0, 50),
  dist_thres = c(200000, Inf),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  focus = 1
) # which entry to use as comparison

clamped <- ATA_comp$signal_raw$Clamped
unclamped <- ATA_comp$signal_raw$Control
loop_common <- as.data.frame(row.names(clamped))

##### Extract data from ATA output
loop_common$chrom <- sapply(strsplit(loop_common$`row.names(clamped)`, ":"), `[`, 1)
loop_common$start <- sapply(strsplit(sapply(strsplit(loop_common$`row.names(clamped)`, ":"), `[`, 2), "-"), `[`, 1)
loop_common$end <- sapply(strsplit(loop_common$`row.names(clamped)`, "-"), `[`, 2)
loop_common$`row.names(clamped)` <- NULL

for (i in 1:nrow(loop_common)) {
  clampedsub <- clamped[i, 1:100, 1:100]
  clampedsub <- clampedsub[order(as.numeric(rownames(clampedsub))), ]
  diag(clampedsub) <- 0

  unclampedsub <- unclamped[i, 1:100, 1:100]
  unclampedsub <- unclampedsub[order(as.numeric(rownames(unclampedsub))), ]
  diag(unclampedsub) <- 0

  loop_common[i, "clampedtad"] <- sum(colSums(clampedsub))
  loop_common[i, "unclampedtad"] <- sum(colSums(unclampedsub))
  print(paste(i, "in", nrow(loop_common)))
}

loop_common$hic_diff <- log2(loop_common$clampedtad / loop_common$unclampedtad)
loop_common$start <- as.numeric(loop_common$start)
loop_common$end <- as.numeric(loop_common$end)
saveRDS(loop_common, "collection_comp_tads")

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
# mat = mat[mat$symbol %in% diffmat$symbol,]


refgene <- readRDS("E:/Compression_MicroC_Full/mustache_analysis/mm39refgene.rds")
# for (i in 1:length(rownames(loop_common))){
#  row_grange <- GRanges(paste0(loop_common[i,"chrom"]), IRanges((loop_common[i,"start"]), (loop_common[i,"end"])))
#  distances = distanceToNearest(row_grange, refgene)
#  refgene.subset = subsetByOverlaps(refgene, row_grange, maxgap = 0)
#  refgene.subset <- refgene.subset[order(distances[subjectHits(findOverlaps(refgene.subset, row_grange))])]
#  loop_common[i,"intragenes"] = paste(unique(as.character(refgene.subset$gene_name)), collapse = "/")
#
#  print(paste("intra ann", i))
# }
# saveRDS(loop_common, "tadcompare_zero/genova_intratad_50000_ann")
loop_common <- readRDS("tadcompare_zero/genova_intratad_50000_ann")


for (i in 1:length(rownames(loop_common))) {
  anchor_genes <- c(strsplit(loop_common[i, "intragenes"], split = "/"))
  anchor_genes <- unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))

  if (length(anchor_genes) > 0) {
    loop_common[i, "mean_exp_anchor"] <- mean(mat[mat$symbol %in% anchor_genes, 7])
    loop_common[i, "mean_exp_rand"] <- mean(sample(mat[is.finite(mat$comp_diff) == T, "comp_diff"], length(intersect(mat$symbol, anchor_genes))))
    loop_common[i, "mean_exp_anchor_clamped"] <- mean(mat[mat$symbol %in% anchor_genes, "clamped"])
    loop_common[i, "mean_exp_anchor_unclamped"] <- mean(mat[mat$symbol %in% anchor_genes, "unclamped"])
  }
}

loop_common <- loop_common[!is.na(loop_common$mean_exp_anchor), ]
loop_common <- loop_common[!is.na(loop_common$hic_diff), ]

loop_common$mean_exp_anchor <- log2(loop_common$mean_exp_anchor)
loop_common <- loop_common[is.finite(loop_common$mean_exp_anchor), ]

loop_common$mean_exp_rand <- log2(loop_common$mean_exp_rand)
loop_common <- loop_common[is.finite(loop_common$mean_exp_anchor), ]

loop_common$size <- log2(as.numeric(loop_common$end) - as.numeric(loop_common$start))

nrow(loop_common) - nrow(loop_common[loop_common$hic_diff <= -log2(1.25), ]) - nrow(loop_common[loop_common$hic_diff >= log2(1.25), ]) # 775
nrow(loop_common[loop_common$hic_diff >= log2(1.25), ]) # 0
nrow(loop_common[loop_common$hic_diff <= -log2(1.25), ]) # 1

data <- data.frame(
  group = c("DEDOWN", "NONDE", "DEUP"),
  value = c(
    nrow(loop_common[loop_common$hic_diff <= -log2(1.25), ]),
    nrow(loop_common) - nrow(loop_common[loop_common$hic_diff <= -log2(1.25), ]) - nrow(loop_common[loop_common$hic_diff >= log2(1.25), ]),
    nrow(loop_common[loop_common$hic_diff >= log2(1.25), ])
  )
)

sampleorder <- c("DEUP", "NONDE", "DEDOWN")
ggplot(data, aes(x = "", y = value, fill = factor(group, level = sampleorder))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 1.57) +
  theme_void() +
  NoLegend() +
  scale_fill_manual(values = c("NONDE" = "grey90", "DEUP" = "#4D9221FF", "DEDOWN" = "#C51B7DFF"))

##### Correlation and Plotting
bin_corr <- data.frame()
width <- 2.5
for (i in 0:14) {
  loop_bin <- loop_common[
    loop_common$clampedtad / 10000 >= 0 + i * width | loop_common$unclampedtad / 10000 >= 0 + i * width,
    c("hic_diff", "mean_exp_anchor")
  ]
  loop_cor <- cor.test(loop_bin$hic_diff, loop_bin$mean_exp_anchor, method = c("spearman"))
  bin_corr[i + 1, "bin"] <- paste("> ", 0 + i * width, sep = "")
  bin_corr[i + 1, "value"] <- loop_cor$estimate[[1]]
  bin_corr[i + 1, "pvalue"] <- loop_cor$p.value
  bin_corr[i + 1, "npvalue"] <- -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("grey95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) + # 3.7 cutoff
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "#red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
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
  xlab(expression(paste("TAD Intensity Threshold (x ", 10^4, ")", sep = ""))) +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-TAD Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.3, 0.3) +
  scale_x_discrete(labels = c("> 0", "", "", "> 9", "", "", "> 18", "", "", "> 27", "", "", "> 36", "", ""))

##### Random Control Correlation
bin_corr <- data.frame()
width <- 2.5
for (i in 0:14) {
  loop_bin <- loop_common[
    loop_common$clampedtad / 10000 >= 0 + i * width | loop_common$unclampedtad / 10000 >= 0 + i * width,
    c("hic_diff", "mean_exp_rand")
  ]
  loop_cor <- cor.test(loop_bin$hic_diff, loop_bin$mean_exp_rand, method = c("spearman"))
  bin_corr[i + 1, "bin"] <- paste("> ", 0 + i * width, sep = "")
  bin_corr[i + 1, "value"] <- loop_cor$estimate[[1]]
  bin_corr[i + 1, "pvalue"] <- loop_cor$p.value
  bin_corr[i + 1, "npvalue"] <- -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("grey95", "red3"))(50))
par(lheight = 0.2)
ggplot(bin_corr, aes(x = factor(bin, level = bin), y = value, fill = npvalue)) + # 3.7 cutoff
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradientn(colors = heat_colors, na.value = "#red3", limits = c(0, 10), label = function(x) sprintf("%.1f", x)) +
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
  xlab(expression(paste("TAD Intensity Threshold (x ", 10^4, ")", sep = ""))) +
  ylab(bquote(atop(NA, atop(
    textstyle("Gene-TAD Correlation"),
    textstyle("(Spearman's " * rho * ")")
  )))) +
  ylim(-0.3, 0.3) +
  scale_x_discrete(labels = c("> 0", "", "", "> 9", "", "", "> 18", "", "", "> 27", "", "", "> 36", "", ""))

saveRDS(loop_common, "comp_tad")
