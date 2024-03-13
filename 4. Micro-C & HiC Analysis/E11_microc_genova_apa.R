###
### APA analysis for feature calls
###

library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(GENOVA)

options(scipen = 999)
setwd("./")
dir <- getwd()


# Matrix formatting for GENOVA
i <- 5000
bedpe <- read.table("hictable_merge_5000.txt", header = T)
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
  signal_path = "genova_input/clamped_5000_hictable_genova.matrix",
  indices_path = "matrix/raw/clamped.5000_abs.bed",
  sample_name = "Clamped",
  colour = "blue",
  scale_bp = NULL,
  balancing = F,
  verbose = T
)
unclamped <- load_contacts(
  signal_path = "genova_input/unclamped_5000_hictable_genova.matrix",
  indices_path = "matrix/raw/unclamped.5000_abs.bed",
  sample_name = "Control",
  colour = "red",
  scale_bp = NULL,
  balancing = F
)

##### APA all loops
feature_common <- readRDS("collection_comp_loops")
feature_apa_bedpe <- feature_common[, 1:6]

feature_apa <- APA(list(" " = clamped),
  size_bin = 25,
  outlier_filter = c(0, 1),
  dist_thres = c(sqrt(2 * (5000 * 25)^2), Inf),
  bedpe = feature_apa_bedpe[, ]
)

visualise(feature_apa,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 1.5))

feature_apa <- APA(list(" " = unclamped),
  size_bin = 25,
  outlier_filter = c(0, 1),
  dist_thres = c(sqrt(2 * (5000 * 25)^2), Inf),
  bedpe = feature_apa_bedpe[, ]
)

visualise(feature_apa,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 1.5))


##### APA upregulated loops
feature_common <- readRDS("collection_comp_loops")
feature_common <- feature_common[is.finite(feature_common$start1), ]
feature_apa_bedpe <- feature_common[feature_common$hic_diff >= log2(1.25), 1:6]
feature_apa_bedpe <- feature_apa_bedpe[is.finite(feature_apa_bedpe$start1), ]

feature_apa <- APA(list(" " = clamped),
  size_bin = 25,
  outlier_filter = c(0, 1),
  dist_thres = c(sqrt(2 * (5000 * 25)^2), Inf),
  bedpe = feature_apa_bedpe[, ]
)

visualise(feature_apa,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 1.5))

feature_apa <- APA(list(" " = unclamped),
  size_bin = 25,
  outlier_filter = c(0, 1),
  dist_thres = c(sqrt(2 * (5000 * 25)^2), Inf),
  bedpe = feature_apa_bedpe[, ]
)

visualise(feature_apa,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 1.5))


##### APA downregulated loops
feature_common <- readRDS("collection_comp_loops")
feature_common <- feature_common[is.finite(feature_common$start1), ]
feature_apa_bedpe <- feature_common[feature_common$hic_diff <= -log2(1.25), 1:6]
feature_apa_bedpe <- feature_apa_bedpe[is.finite(feature_apa_bedpe$start1), ]

feature_apa <- APA(list(" " = clamped),
  size_bin = 25,
  outlier_filter = c(0, 1),
  dist_thres = c(sqrt(2 * (5000 * 25)^2), Inf),
  bedpe = feature_apa_bedpe[, ]
)

visualise(feature_apa,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 1.5))

feature_apa <- APA(list(" " = unclamped),
  size_bin = 25,
  outlier_filter = c(0, 1),
  dist_thres = c(sqrt(2 * (5000 * 25)^2), Inf),
  bedpe = feature_apa_bedpe[, ]
)

visualise(feature_apa,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 1.5))

##### APA all stripes
feature_common <- readRDS("collection_comp_stripes")
for (i in 1:nrow(feature_common)) {
  if (feature_common[i, "dir"] == "left") {
    feature_common[i, "start2"] <- feature_common[i, "start1"]
    feature_common[i, "end2"] <- feature_common[i, "end1"]
  }
}

feature_apa_bedpe <- feature_common[, 1:6]
feature_ara <- ARA(list(" " = clamped),
  size_bin = 100,
  outlier_filter = c(0, 1),
  bed = feature_apa_bedpe
)

visualise(feature_ara,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 2.5))

feature_apa_bedpe <- feature_common[, 1:6]
feature_ara <- ARA(list(" " = unclamped),
  size_bin = 100,
  outlier_filter = c(0, 1),
  bed = feature_apa_bedpe
)

visualise(feature_ara,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 2.5))


##### APA downregulated stripes
feature_apa_bedpe <- feature_common[feature_common$hic_diff <= -log2(1.25), 1:6]
feature_apa_bedpe <- feature_apa_bedpe[is.finite(feature_apa_bedpe$start1), ]

feature_ara <- ARA(list(" " = clamped),
  size_bin = 100,
  outlier_filter = c(0, 1),
  bed = feature_apa_bedpe
)

visualise(feature_ara,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 2.5)) #

feature_ara <- ARA(list(" " = unclamped),
  size_bin = 100,
  outlier_filter = c(0, 1),
  bed = feature_apa_bedpe
)

visualise(feature_ara,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 2.5))


##### APA all boundaries
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

feature_common <- readRDS("collection_comp_boundaries")
feature_common <- feature_common[feature_common$boundarygenes != "", ]

feature_apa_bedpe <- feature_common[, c(8, 1, 1)]
feature_ara <- ARA(list(" " = clamped),
  size_bin = 100,
  outlier_filter = c(0, 1),
  bed = feature_apa_bedpe
)

visualise(feature_ara,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 2.5))


feature_apa_bedpe <- feature_common[, 1:6]
feature_ara <- ARA(list(" " = unclamped),
  size_bin = 100,
  outlier_filter = c(0, 1),
  bed = feature_apa_bedpe
)

visualise(feature_ara,
  title = "",
  colour_lim = c(0, 1.5),
  colour_lim_contrast = c(-1, 1),
  metric = "diff",
  contrast = 1
) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(fill = "Contacts") +
  ggplot2::scale_fill_gradientn(colours = c("white", "#FF2345", "black"), limits = c(0, 2.5))
