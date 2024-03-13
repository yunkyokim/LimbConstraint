###
### Hierarchical merging for loop and stripe features called from chromosight
###


library(rtracklayer)
library(GenomicRanges)

options(scipen = 999)
setwd("./")
dir = getwd()

res = c(1000000,
        500000,
        250000,
        100000,
        50000,
        40000,
        25000,
        20000,
        10000,
        5000,
        2500)

##### Loop Merging
i = 1
highres = read.table(
  paste("clamped_", res[i], "_loops_clampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_unclamped = read.table(
  paste("clamped_", res[i], "_loops_unclampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_unclamped = highres_unclamped$score
highres$pvalue_unclamped = highres_unclamped$pvalue
highres$qvalue_unclamped = highres_unclamped$qvalue
highres = highres[!(is.na(highres$score) &
                      is.na(highres$score_unclamped)), ]

for (i in 1:(length(res) - 1)) {
  lowres = read.table(
    paste("clamped_", res[i + 1], "_loops_clampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_unclamped = read.table(
    paste("clamped_", res[i + 1], "_loops_unclampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres$score_unclamped = lowres_unclamped$score
  lowres$pvalue_unclamped = lowres_unclamped$pvalue
  lowres$qvalue_unclamped = lowres_unclamped$qvalue
  lowres = lowres[!(is.na(lowres$score) &
                      is.na(lowres$score_unclamped)), ]
  
  highres[, 15] = "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (lowres[b, 1] == highres[a, 1] &
            lowres[b, 2] >= highres[a, 2] &
            lowres[b, 3] <= highres[a, 3] &
            lowres[b, 4] == highres[a, 4] &
            lowres[b, 5] >= highres[a, 5] &
            lowres[b, 6] <= highres[a, 6]) {
          highres[a, 15] = "Duplicate"
          break
        }
      }
    }
    lowres[, 15] = "Unique"
    highres = rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres = subset(highres, V15 != "Duplicate")
}

clamped = highres
saveRDS(clamped, "clamped_loops_merged")

i = 1

highres = read.table(
  paste("unclamped_", res[i], "_loops_clampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_unclamped = read.table(
  paste("unclamped_", res[i], "_loops_unclampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_unclamped = highres_unclamped$score
highres$pvalue_unclamped = highres_unclamped$pvalue
highres$qvalue_unclamped = highres_unclamped$qvalue
highres = highres[!(is.na(highres$score) &
                      is.na(highres$score_unclamped)), ]

for (i in 1:(length(res) - 1)) {
  lowres = read.table(
    paste("unclamped_", res[i + 1], "_loops_clampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_unclamped = read.table(
    paste("unclamped_", res[i + 1], "_loops_unclampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres$score_unclamped = lowres_unclamped$score
  lowres$pvalue_unclamped = lowres_unclamped$pvalue
  lowres$qvalue_unclamped = lowres_unclamped$qvalue
  lowres = lowres[!(is.na(lowres$score) &
                      is.na(lowres$score_unclamped)), ]
  
  highres[, 15] = "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (lowres[b, 1] == highres[a, 1] &
            lowres[b, 2] >= highres[a, 2] &
            lowres[b, 3] <= highres[a, 3] &
            lowres[b, 4] == highres[a, 4] &
            lowres[b, 5] >= highres[a, 5] &
            lowres[b, 6] <= highres[a, 6]) {
          highres[a, 15] = "Duplicate"
          break
        }
      }
    }
    lowres[, 15] = "Unique"
    highres = rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres = subset(highres, V15 != "Duplicate")
}
unclamped = highres
saveRDS(unclamped, "unclamped_loops_merged")

loop_mustache = rbind(clamped, unclamped)
loop_mustache = loop_mustache[!duplicated(loop_mustache[, 1:6]), ]
write.table(
  loop_mustache[, 1:6],
  paste("loops_E11_merged.txt", sep = ""),
  col.names = F,
  row.names = F,
  quote = F ,
  sep = "\t"
)

##### Stripes Left Merging
i = 1
highres = read.table(
  paste("clamped_", res[i], "_stripes_left_clampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_unclamped = read.table(
  paste("clamped_", res[i], "_stripes_left_unclampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_unclamped = highres_unclamped$score
highres$pvalue_unclamped = highres_unclamped$pvalue
highres$qvalue_unclamped = highres_unclamped$qvalue
highres = highres[!(is.na(highres$score) &
                      is.na(highres$score_unclamped)), ]

for (i in 1:(length(res) - 1)) {
  lowres = read.table(
    paste("clamped_", res[i + 1], "_stripes_left_clampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_unclamped = read.table(
    paste("clamped_", res[i + 1], "_stripes_left_unclampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres$score_unclamped = lowres_unclamped$score
  lowres$pvalue_unclamped = lowres_unclamped$pvalue
  lowres$qvalue_unclamped = lowres_unclamped$qvalue
  lowres = lowres[!(is.na(lowres$score) &
                      is.na(lowres$score_unclamped)), ]
  
  highres[, 15] = "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (lowres[b, 1] == highres[a, 1] &
            lowres[b, 2] >= highres[a, 2] &
            lowres[b, 3] <= highres[a, 3]) {
            #lowres[b,4] == highres[a,4] &
            #lowres[b,5] >= highres[a,5] &
            #lowres[b,6] <= highres[a,6]) {
            highres[a, 15] = "Duplicate"
            break
      }
    }
  }
  lowres[, 15] = "Unique"
  highres = rbind(highres, lowres)
  print(table(highres$V15))
}
highres = subset(highres, V15 != "Duplicate")
}

clamped = highres
saveRDS(clamped, "clamped_stripes_left_merged")

i = 1

highres = read.table(
  paste("unclamped_", res[i], "_stripes_left_clampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_unclamped = read.table(
  paste("unclamped_", res[i], "_stripes_left_unclampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_unclamped = highres_unclamped$score
highres$pvalue_unclamped = highres_unclamped$pvalue
highres$qvalue_unclamped = highres_unclamped$qvalue
highres = highres[!(is.na(highres$score) &
                      is.na(highres$score_unclamped)), ]

for (i in 1:(length(res) - 1)) {
  lowres = read.table(
    paste("unclamped_", res[i + 1], "_stripes_left_clampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_unclamped = read.table(
    paste(
      "unclamped_",
      res[i + 1],
      "_stripes_left_unclampedquant.tsv",
      sep = ""
    ),
    header = T,
    fill = TRUE
  )
  lowres$score_unclamped = lowres_unclamped$score
  lowres$pvalue_unclamped = lowres_unclamped$pvalue
  lowres$qvalue_unclamped = lowres_unclamped$qvalue
  lowres = lowres[!(is.na(lowres$score) &
                      is.na(lowres$score_unclamped)), ]
  
  highres[, 15] = "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (lowres[b, 1] == highres[a, 1] &
            lowres[b, 2] >= highres[a, 2] &
            lowres[b, 3] <= highres[a, 3]) {
            #lowres[b,4] == highres[a,4] &
            #lowres[b,5] >= highres[a,5] &
            #lowres[b,6] <= highres[a,6]) {
            highres[a, 15] = "Duplicate"
            break
      }
    }
  }
  lowres[, 15] = "Unique"
  highres = rbind(highres, lowres)
  print(table(highres$V15))
}
highres = subset(highres, V15 != "Duplicate")
}
unclamped = highres
saveRDS(unclamped, "unclamped_stripes_left_merged")

loop_mustache = rbind(clamped, unclamped)
loop_mustache = loop_mustache[!duplicated(loop_mustache[, 1:6]), ]
write.table(
  loop_mustache[, 1:6],
  paste("stripes_left_E11_merged.txt", sep = ""),
  col.names = F,
  row.names = F,
  quote = F ,
  sep = "\t"
)


##### Stripes Right Merging
i = 1
highres = read.table(
  paste("clamped_", res[i], "_stripes_right_clampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_unclamped = read.table(
  paste("clamped_", res[i], "_stripes_right_unclampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_unclamped = highres_unclamped$score
highres$pvalue_unclamped = highres_unclamped$pvalue
highres$qvalue_unclamped = highres_unclamped$qvalue
highres = highres[!(is.na(highres$score) &
                      is.na(highres$score_unclamped)), ]

for (i in 1:(length(res) - 1)) {
  lowres = read.table(
    paste("clamped_", res[i + 1], "_stripes_right_clampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_unclamped = read.table(
    paste("clamped_", res[i + 1], "_stripes_right_unclampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres$score_unclamped = lowres_unclamped$score
  lowres$pvalue_unclamped = lowres_unclamped$pvalue
  lowres$qvalue_unclamped = lowres_unclamped$qvalue
  lowres = lowres[!(is.na(lowres$score) &
                      is.na(lowres$score_unclamped)), ]
  
  highres[, 15] = "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (#lowres[b,1] == highres[a,1] &
          #lowres[b,2] >= highres[a,2] &
          #lowres[b,3] <= highres[a,3] &
          lowres[b, 4] == highres[a, 4] &
          lowres[b, 5] >= highres[a, 5] &
          lowres[b, 6] <= highres[a, 6]) {
          highres[a, 15] = "Duplicate"
          break
        }
      }
    }
    lowres[, 15] = "Unique"
    highres = rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres = subset(highres, V15 != "Duplicate")
}

clamped = highres
saveRDS(clamped, "clamped_stripes_right_merged")

i = 1

highres = read.table(
  paste("unclamped_", res[i], "_stripes_right_clampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres_unclamped = read.table(
  paste("unclamped_", res[i], "_stripes_right_unclampedquant.tsv", sep = ""),
  header = T,
  fill = TRUE
)
highres$score_unclamped = highres_unclamped$score
highres$pvalue_unclamped = highres_unclamped$pvalue
highres$qvalue_unclamped = highres_unclamped$qvalue
highres = highres[!(is.na(highres$score) &
                      is.na(highres$score_unclamped)), ]

for (i in 1:(length(res) - 1)) {
  lowres = read.table(
    paste("unclamped_", res[i + 1], "_stripes_right_clampedquant.tsv", sep = ""),
    header = T,
    fill = TRUE
  )
  lowres_unclamped = read.table(
    paste(
      "unclamped_",
      res[i + 1],
      "_stripes_right_unclampedquant.tsv",
      sep = ""
    ),
    header = T,
    fill = TRUE
  )
  lowres$score_unclamped = lowres_unclamped$score
  lowres$pvalue_unclamped = lowres_unclamped$pvalue
  lowres$qvalue_unclamped = lowres_unclamped$qvalue
  lowres = lowres[!(is.na(lowres$score) &
                      is.na(lowres$score_unclamped)), ]
  
  highres[, 15] = "Unique"
  if (length(rownames(lowres)) > 0) {
    for (a in 1:length(rownames(highres))) {
      for (b in 1:length(rownames(lowres))) {
        if (#lowres[b,1] == highres[a,1] &
          #lowres[b,2] >= highres[a,2] &
          #lowres[b,3] <= highres[a,3] &
          lowres[b, 4] == highres[a, 4] &
          lowres[b, 5] >= highres[a, 5] &
          lowres[b, 6] <= highres[a, 6]) {
          highres[a, 15] = "Duplicate"
          break
        }
      }
    }
    lowres[, 15] = "Unique"
    highres = rbind(highres, lowres)
    print(table(highres$V15))
  }
  highres = subset(highres, V15 != "Duplicate")
}
unclamped = highres
saveRDS(unclamped, "unclamped_stripes_right_merged")

loop_mustache = rbind(clamped, unclamped)
loop_mustache = loop_mustache[!duplicated(loop_mustache[, 1:6]), ]
write.table(
  loop_mustache[, 1:6],
  paste("stripes_right_E11_merged.txt", sep = ""),
  col.names = F,
  row.names = F,
  quote = F ,
  sep = "\t"
)
