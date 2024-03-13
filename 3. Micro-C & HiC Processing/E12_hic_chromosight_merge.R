###
### Hierarchical merging for loop and stripe features called from chromosight
###


library(rtracklayer)
library(GenomicRanges)

options(scipen = 999)
setwd("./")
dir <- getwd()

res <- c(
  1000000,
  500000,
  250000,
  100000,
  50000,
  40000,
  25000,
  20000,
  10000,
  5000,
  2500
)

##### Loop Merging
i <- 1
highres <- read.table(
  paste("distal_", res[i], "_loops_distalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_proximal <- read.table(
  paste("distal_", res[i], "_loops_proximalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_proximal <- highres_proximal$score
highres$pvalue_proximal <- highres_proximal$pvalue
highres$qvalue_proximal <- highres_proximal$qvalue
highres <- highres[!(is.na(highres$score) &
  is.na(highres$score_proximal)), ]

for (i in 1:(length(res) - 1)) {
  lowres <- read.table(
    paste("distal_", res[i + 1], "_loops_distalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_proximal <- read.table(
    paste("distal_", res[i + 1], "_loops_proximalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres$score_proximal <- lowres_proximal$score
  lowres$pvalue_proximal <- lowres_proximal$pvalue
  lowres$qvalue_proximal <- lowres_proximal$qvalue
  lowres <- lowres[!(is.na(lowres$score) &
    is.na(lowres$score_proximal)), ]

  highres[, 15] <- "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (lowres[b, 1] == highres[a, 1] &
          lowres[b, 2] >= highres[a, 2] &
          lowres[b, 3] <= highres[a, 3] &
          lowres[b, 4] == highres[a, 4] &
          lowres[b, 5] >= highres[a, 5] &
          lowres[b, 6] <= highres[a, 6]) {
          highres[a, 15] <- "Duplicate"
          break
        }
      }
    }
    lowres[, 15] <- "Unique"
    highres <- rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres <- subset(highres, V15 != "Duplicate")
}

distal <- highres
saveRDS(distal, "distal_loops_merged")

i <- 1

highres <- read.table(
  paste("proximal_", res[i], "_loops_distalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_proximal <- read.table(
  paste("proximal_", res[i], "_loops_proximalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_proximal <- highres_proximal$score
highres$pvalue_proximal <- highres_proximal$pvalue
highres$qvalue_proximal <- highres_proximal$qvalue
highres <- highres[!(is.na(highres$score) &
  is.na(highres$score_proximal)), ]

for (i in 1:(length(res) - 1)) {
  lowres <- read.table(
    paste("proximal_", res[i + 1], "_loops_distalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_proximal <- read.table(
    paste("proximal_", res[i + 1], "_loops_proximalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres$score_proximal <- lowres_proximal$score
  lowres$pvalue_proximal <- lowres_proximal$pvalue
  lowres$qvalue_proximal <- lowres_proximal$qvalue
  lowres <- lowres[!(is.na(lowres$score) &
    is.na(lowres$score_proximal)), ]

  highres[, 15] <- "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (lowres[b, 1] == highres[a, 1] &
          lowres[b, 2] >= highres[a, 2] &
          lowres[b, 3] <= highres[a, 3] &
          lowres[b, 4] == highres[a, 4] &
          lowres[b, 5] >= highres[a, 5] &
          lowres[b, 6] <= highres[a, 6]) {
          highres[a, 15] <- "Duplicate"
          break
        }
      }
    }
    lowres[, 15] <- "Unique"
    highres <- rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres <- subset(highres, V15 != "Duplicate")
}
proximal <- highres
saveRDS(proximal, "proximal_loops_merged")

loop_mustache <- rbind(distal, proximal)
loop_mustache <- loop_mustache[!duplicated(loop_mustache[, 1:6]), ]
write.table(
  loop_mustache[, 1:6],
  paste("loops_E11_merged.txt", sep = ""),
  col.names = F,
  row.names = F,
  quote = F,
  sep = "\t"
)

##### Stripes Left Merging
i <- 1
highres <- read.table(
  paste("distal_", res[i], "_stripes_left_distalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_proximal <- read.table(
  paste("distal_", res[i], "_stripes_left_proximalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_proximal <- highres_proximal$score
highres$pvalue_proximal <- highres_proximal$pvalue
highres$qvalue_proximal <- highres_proximal$qvalue
highres <- highres[!(is.na(highres$score) &
  is.na(highres$score_proximal)), ]

for (i in 1:(length(res) - 1)) {
  lowres <- read.table(
    paste("distal_", res[i + 1], "_stripes_left_distalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_proximal <- read.table(
    paste("distal_", res[i + 1], "_stripes_left_proximalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres$score_proximal <- lowres_proximal$score
  lowres$pvalue_proximal <- lowres_proximal$pvalue
  lowres$qvalue_proximal <- lowres_proximal$qvalue
  lowres <- lowres[!(is.na(lowres$score) &
    is.na(lowres$score_proximal)), ]

  highres[, 15] <- "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (lowres[b, 1] == highres[a, 1] &
          lowres[b, 2] >= highres[a, 2] &
          lowres[b, 3] <= highres[a, 3]) {
          # lowres[b,4] == highres[a,4] &
          # lowres[b,5] >= highres[a,5] &
          # lowres[b,6] <= highres[a,6]) {
          highres[a, 15] <- "Duplicate"
          break
        }
      }
    }
    lowres[, 15] <- "Unique"
    highres <- rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres <- subset(highres, V15 != "Duplicate")
}

distal <- highres
saveRDS(distal, "distal_stripes_left_merged")

i <- 1

highres <- read.table(
  paste("proximal_", res[i], "_stripes_left_distalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_proximal <- read.table(
  paste("proximal_", res[i], "_stripes_left_proximalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_proximal <- highres_proximal$score
highres$pvalue_proximal <- highres_proximal$pvalue
highres$qvalue_proximal <- highres_proximal$qvalue
highres <- highres[!(is.na(highres$score) &
  is.na(highres$score_proximal)), ]

for (i in 1:(length(res) - 1)) {
  lowres <- read.table(
    paste("proximal_", res[i + 1], "_stripes_left_distalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_proximal <- read.table(
    paste(
      "proximal_",
      res[i + 1],
      "_stripes_left_proximalquant.tsv",
      sep = ""
    ),
    header = T,
    fill = TRUE
  )
  lowres$score_proximal <- lowres_proximal$score
  lowres$pvalue_proximal <- lowres_proximal$pvalue
  lowres$qvalue_proximal <- lowres_proximal$qvalue
  lowres <- lowres[!(is.na(lowres$score) &
    is.na(lowres$score_proximal)), ]

  highres[, 15] <- "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (lowres[b, 1] == highres[a, 1] &
          lowres[b, 2] >= highres[a, 2] &
          lowres[b, 3] <= highres[a, 3]) {
          # lowres[b,4] == highres[a,4] &
          # lowres[b,5] >= highres[a,5] &
          # lowres[b,6] <= highres[a,6]) {
          highres[a, 15] <- "Duplicate"
          break
        }
      }
    }
    lowres[, 15] <- "Unique"
    highres <- rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres <- subset(highres, V15 != "Duplicate")
}
proximal <- highres
saveRDS(proximal, "proximal_stripes_left_merged")

loop_mustache <- rbind(distal, proximal)
loop_mustache <- loop_mustache[!duplicated(loop_mustache[, 1:6]), ]
write.table(
  loop_mustache[, 1:6],
  paste("stripes_left_E11_merged.txt", sep = ""),
  col.names = F,
  row.names = F,
  quote = F,
  sep = "\t"
)


##### Stripes Right Merging
i <- 1
highres <- read.table(
  paste("distal_", res[i], "_stripes_right_distalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_proximal <- read.table(
  paste("distal_", res[i], "_stripes_right_proximalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_proximal <- highres_proximal$score
highres$pvalue_proximal <- highres_proximal$pvalue
highres$qvalue_proximal <- highres_proximal$qvalue
highres <- highres[!(is.na(highres$score) &
  is.na(highres$score_proximal)), ]

for (i in 1:(length(res) - 1)) {
  lowres <- read.table(
    paste("distal_", res[i + 1], "_stripes_right_distalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_proximal <- read.table(
    paste("distal_", res[i + 1], "_stripes_right_proximalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres$score_proximal <- lowres_proximal$score
  lowres$pvalue_proximal <- lowres_proximal$pvalue
  lowres$qvalue_proximal <- lowres_proximal$qvalue
  lowres <- lowres[!(is.na(lowres$score) &
    is.na(lowres$score_proximal)), ]

  highres[, 15] <- "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if ( # lowres[b,1] == highres[a,1] &
          # lowres[b,2] >= highres[a,2] &
          # lowres[b,3] <= highres[a,3] &
          lowres[b, 4] == highres[a, 4] &
            lowres[b, 5] >= highres[a, 5] &
            lowres[b, 6] <= highres[a, 6]) {
          highres[a, 15] <- "Duplicate"
          break
        }
      }
    }
    lowres[, 15] <- "Unique"
    highres <- rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres <- subset(highres, V15 != "Duplicate")
}

distal <- highres
saveRDS(distal, "distal_stripes_right_merged")

i <- 1

highres <- read.table(
  paste("proximal_", res[i], "_stripes_right_distalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_proximal <- read.table(
  paste("proximal_", res[i], "_stripes_right_proximalquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_proximal <- highres_proximal$score
highres$pvalue_proximal <- highres_proximal$pvalue
highres$qvalue_proximal <- highres_proximal$qvalue
highres <- highres[!(is.na(highres$score) &
  is.na(highres$score_proximal)), ]

for (i in 1:(length(res) - 1)) {
  lowres <- read.table(
    paste("proximal_", res[i + 1], "_stripes_right_distalquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_proximal <- read.table(
    paste(
      "proximal_",
      res[i + 1],
      "_stripes_right_proximalquant.tsv",
      sep = ""
    ),
    header = T,
    fill = TRUE
  )
  lowres$score_proximal <- lowres_proximal$score
  lowres$pvalue_proximal <- lowres_proximal$pvalue
  lowres$qvalue_proximal <- lowres_proximal$qvalue
  lowres <- lowres[!(is.na(lowres$score) &
    is.na(lowres$score_proximal)), ]

  highres[, 15] <- "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if ( # lowres[b,1] == highres[a,1] &
          # lowres[b,2] >= highres[a,2] &
          # lowres[b,3] <= highres[a,3] &
          lowres[b, 4] == highres[a, 4] &
            lowres[b, 5] >= highres[a, 5] &
            lowres[b, 6] <= highres[a, 6]) {
          highres[a, 15] <- "Duplicate"
          break
        }
      }
    }
    lowres[, 15] <- "Unique"
    highres <- rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres <- subset(highres, V15 != "Duplicate")
}
proximal <- highres
saveRDS(proximal, "proximal_stripes_right_merged")

loop_mustache <- rbind(distal, proximal)
loop_mustache <- loop_mustache[!duplicated(loop_mustache[, 1:6]), ]
write.table(
  loop_mustache[, 1:6],
  paste("stripes_right_E11_merged.txt", sep = ""),
  col.names = F,
  row.names = F,
  quote = F,
  sep = "\t"
)
