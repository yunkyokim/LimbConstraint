###
### Feature Annotation for Stripes Left Calls
### Analogous approach used for E12.5 distal/proximal Hi-C datasets


library(rtracklayer)
library(GenomicRanges)

options(scipen = 999)
setwd("./")

##### HiCcompare Annotation
##### Here we annotate the stripes by the HiCcompare (loess) normalized interaction using 1kb-resolution data
##### stripes left - gate by start anchor
#####
loop_common <- readRDS("stripes_left_merged")
loop_common$loop_num <- seq.int(nrow(loop_common))

hicmerge <- read.table("hictable_merge_1000.txt", header = T) # concatenated file of loess-normalized hictables from all chromosomes
hicmerge[, "chr1"] <- sub("chr", "", hicmerge[, "chr1"])
hicmerge[, "chr2"] <- sub("chr", "", hicmerge[, "chr2"])
hicmerge$start1 <- as.numeric(hicmerge$start1)
hicmerge <- hicmerge[!is.na(hicmerge$start1), ]

for (i in 1:length(rownames(loop_common))) {
  start_subset <- seq(
    from = loop_common[i, "start1"],
    to = loop_common[i, "end1"] - 1000,
    by = 1000
  )
  end_subset <- seq(
    from = loop_common[i, "start1"],
    to = loop_common[i, "end2"] - 1000,
    by = 1000
  ) # take all contacts up to the 2nd "anchor" call by chromosight

  hicmerge_subset <- hicmerge[hicmerge$chr1 == loop_common[i, "chrom1"] &
    hicmerge$start1 %in% start_subset &
    hicmerge$start2 %in% end_subset, c("adj.IF1", "adj.IF2")]
  loop_common[i, "adj.IF1"] <- sum(as.numeric(hicmerge_subset$adj.IF1))
  loop_common[i, "adj.IF2"] <- sum(as.numeric(hicmerge_subset$adj.IF2))
  loop_common[i, "SD.IF1"] <- sd(as.numeric(hicmerge_subset$adj.IF1))
  loop_common[i, "SD.IF2"] <- sd(as.numeric(hicmerge_subset$adj.IF2))
  loop_common[i, "sparsity"] <- length(hicmerge_subset$adj.IF1)

  print(paste("loop", i, "of", length(rownames(loop_common))))
}

saveRDS(loop_common, "stripes_left_merged_hictable")


##### Anchor Genomic Feature Annotation
##### Here we annotate the stripes by the genes located at the anchor, as well as ChromHMM state
refgene <- readRDS("mm39refgene.rds")
E11_18_chromhmm <- readRDS("E11_18state_chromhmm.rds")

loop_common$hic_diff <- log2((as.numeric(loop_common$adj.IF1) + 0) / (as.numeric(loop_common$adj.IF2) +
  0))
loop_common <- loop_common[, !(names(loop_common) %in% c("bin1", "bin2"))]

for (i in 1:length(rownames(loop_common))) {
  row_grange <-
    GRanges(
      paste0("chr", loop_common[i, 1]),
      IRanges(loop_common[i, 2], loop_common[i, 3])
    )
  overlap <- findOverlaps(row_grange, E11_18_chromhmm)
  chromhmm.subset <- unlist(E11_18_chromhmm[subjectHits(overlap)])
  chromhmm.subset <- chromhmm.subset[order(width(chromhmm.subset), decreasing = T)]
  annotation <- chromhmm.subset
  loop_common[i, "A1_ann"] <- paste(unique(annotation$name), collapse = "/")

  row_grange <-
    GRanges(paste0(loop_common[i, 1]), IRanges(loop_common[i, 2], loop_common[i, 3]))
  distances <- distanceToNearest(row_grange, refgene)
  refgene.subset <- subsetByOverlaps(refgene, row_grange, maxgap = (width(row_grange) *
    0))
  refgene.subset <-
    refgene.subset[order(distances[subjectHits(findOverlaps(refgene.subset, row_grange))])]
  loop_common[i, "A1_gene"] <- paste(unique(as.character(refgene.subset$gene_name)), collapse = "/")

  print(paste("stripe", i, "of", length(rownames(loop_common))))
}

##### Blacklisted Region Removal
blacklist <- readRDS("mm39_blacklist.rds")
for (i in 1:length(rownames(loop_common))) {
  row_grange <-
    GRanges(
      paste0("chr", loop_common[i, 1]),
      IRanges(loop_common[i, 2], loop_common[i, 3])
    )
  overlap <- findOverlaps(row_grange, blacklist)
  chromhmm.subset <- blacklist[subjectHits(overlap)]
  chromhmm.subset <- chromhmm.subset[order(width(chromhmm.subset), decreasing = T)]
  annotation <- chromhmm.subset
  loop_common[i, "A1_blacklist"] <- paste(unique(annotation$name), collapse = "/")

  row_grange <-
    GRanges(
      paste0("chr", loop_common[i, 4]),
      IRanges(loop_common[i, 5], loop_common[i, 6])
    )
  overlap <- findOverlaps(row_grange, blacklist)
  chromhmm.subset <- blacklist[subjectHits(overlap)]
  chromhmm.subset <- chromhmm.subset[order(width(chromhmm.subset), decreasing = T)]
  annotation <- chromhmm.subset
  loop_common[i, "A2_blacklist"] <- paste(unique(annotation$name), collapse = "/")

  print(paste("BLACKLIST stripe", i, "of", length(rownames(loop_common))))
}
loop_common <- loop_common[loop_common$A1_blacklist == "" &
  loop_common$A2_blacklist == "", ]

saveRDS(loop_common, "stripes_left_merged_hictable_ann")

##### Anchor Genomic Feature Annotation
##### Here we annotate the anchors by MACS2 significant ChIP-seq peak calls
chip_files <- list.files("chip_files") # contains all ChIP-seq peak files for relevant histone modifications

for (j in 1:length(chip_files)) {
  peaks <- read.table(paste("chip_files/", chip_files[j], sep = ""))

  for (i in 1:length(rownames(loop_common))) {
    peaks_subset <- peaks[peaks$V1 == loop_common[i, "chrom1"] &
      peaks$V2 >= loop_common[i, "start1"] &
      peaks$V3 <= loop_common[i, "end1"], ]
    loop_common[i, chip_files[j]] <- mean(peaks_subset$V7)
    print(paste(chip_files[j], i, "out of", length(rownames(loop_common))))
  }
}

saveRDS(loop_common, "stripes_left_merged_hictable_ann_chip")
