###
### Analysis of borders by identified by TADcompare
### 

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(Seurat)
library(dplyr)
library(ggplot2)


setwd("./")
dir = getwd()

##### Merge boundary calls for all chromosomes
jointtad = read.table(paste("tadcompare_50000_", 1, "boundaryscores", sep = ""))
jointtad$Chromosome = 1
for (i in c(2:19, "X")){
  tad_mat = read.table(paste("tadcompare_50000_", i, "boundaryscores", sep = ""))
  tad_mat$Chromosome = i
  jointtad = rbind(jointtad, tad_mat)                 
}
saveRDS(jointtad, "joint_50000")

test = jointtad = read.table(paste("tadcompare_50000_", 1, "tadframe", sep = ""))

##### Annotation of boundary genes
refgene = readRDS("mm39refgene.rds")

jointtad = readRDS("joint_50000")
for (i in 1:length(rownames(jointtad))){
  row_grange <- GRanges(paste0(jointtad[i,"Chromosome"]), IRanges((jointtad[i,"Boundary"]-50000/2), (jointtad[i,"Boundary"]+50000/2)))
  distances = distanceToNearest(row_grange, refgene)
  refgene.subset = subsetByOverlaps(refgene, row_grange, maxgap = 0)
  refgene.subset <- refgene.subset[order(distances[subjectHits(findOverlaps(refgene.subset, row_grange))])]
  jointtad[i,"boundarygenes"] = paste(unique(as.character(refgene.subset$gene_name)), collapse = "/")
  
  print(paste("joint", i))
}
saveRDS(jointtad, "joint_50000_ann")

##### Integration with clamp/unclamped RNA-seq
set.seed(123456789)
mat = read.csv("deseq_counts_clamp_only.csv")
mat$comp_diff = rowMeans(mat[,3:4])/rowMeans(mat[,5:6])
mat$clamped = rowMeans(mat[,3:4])
mat$unclamped = rowMeans(mat[,5:6])
mat = mat[!is.na(mat$comp_diff),]
mat = mat[!is.nan(mat$comp_diff),]
mat = mat[is.finite(mat$comp_diff),]
mat = mat[abs(mat$comp_diff) > 0,]
orgmat = mat

mat = mat[mat$CLAMPED1>=10 & 
            mat$CLAMPED2>=10 | 
            mat$UNCLAMPED1>=10 & 
            mat$UNCLAMPED2>=10,]
#mat = mat[mat$symbol %in% diffmat$symbol,]

loop_common = readRDS("joint_50000_ann") #load in data
loop_common = loop_common[loop_common$boundarygenes != "",]

for (i in 1:length(rownames(loop_common))){
  anchor_genes = c(strsplit(loop_common[i,"boundarygenes"], split = "/"))
  anchor_genes = unlist(lapply(anchor_genes, function(x) x[nchar(x) >= 1]))
  
  if (length(anchor_genes)>0){
    loop_common[i,"mean_exp_anchor"] = mean(mat[mat$symbol %in% anchor_genes,7])
    loop_common[i,"mean_exp_rand"] = mean(sample(mat[is.finite(mat$comp_diff) == T, "comp_diff"],length(intersect(mat$symbol,anchor_genes)))) #equal number of random genes assigned
    loop_common[i,"mean_exp_anchor_clamped"] = mean(mat[mat$symbol %in% anchor_genes,"clamped"])
    loop_common[i,"mean_exp_anchor_unclamped"] = mean(mat[mat$symbol %in% anchor_genes,"unclamped"])
  }
}

loop_common = loop_common[!is.na(loop_common$mean_exp_anchor),]
loop_common$mean_exp_anchor = log2(loop_common$mean_exp_anchor)
loop_common = loop_common[is.finite(loop_common$mean_exp_anchor),]
loop_common$mean_exp_rand = log2(loop_common$mean_exp_rand)
loop_common = loop_common[is.finite(loop_common$mean_exp_anchor),]
loop_common = loop_common[!is.na(loop_common$mean_exp_anchor_clamped),]
loop_common = loop_common[!is.na(loop_common$mean_exp_anchor_unclamped),]
loop_common$Gap_Score = loop_common$Gap_Score*-1 #want positve values to reflect matrix 2 enrichment, which is Clamped
loop_common = loop_common[loop_common$TAD_Score1 >= 0.5 | loop_common$TAD_Score2 >= 0.5,]

nrow(loop_common) - nrow(loop_common[loop_common$Gap_Score<= -2,]) - nrow(loop_common[loop_common$Gap_Score>= 2,]) #5000
nrow(loop_common[loop_common$Gap_Score>= 2,]) #312
nrow(loop_common[loop_common$Gap_Score<= -2,]) #308

data <- data.frame(
  group=c("DEDOWN","NONDE", "DEUP"),
  value=c(nrow(loop_common[loop_common$Gap_Score<= -2 ,]),
          nrow(loop_common) - nrow(loop_common[loop_common$Gap_Score<= -2,]) - nrow(loop_common[loop_common$Gap_Score>= 2,]),
          nrow(loop_common[loop_common$Gap_Score>= 2,])))

sampleorder = c("DEUP","NONDE", "DEDOWN")
ggplot(data, aes(x="", y=value, fill=factor(group,level = sampleorder))) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=1.57) +
  theme_void() +
  NoLegend()+
  scale_fill_manual(values = c("NONDE" = "grey90","DEUP" =  "#4D9221FF","DEDOWN" =  "#C51B7DFF"))#comp_boundaries_pie, 400x400

##### Correlation Analysis
bin_corr = data.frame()
width = 0.3
for (i in 0:14){
  loop_bin = loop_common[loop_common$TAD_Score1 >= 0+i*width | loop_common$TAD_Score2 >= 0+i*width,
                         c("Gap_Score", "mean_exp_anchor")]
  loop_cor = cor.test(loop_bin$Gap_Score, loop_bin$mean_exp_anchor, method = c("spearman"))
  bin_corr[i+1,"bin"] = paste("> ", 0+i*width, sep = "")
  bin_corr[i+1,"value"] = loop_cor$estimate[[1]]
  bin_corr[i+1,"pvalue"] = loop_cor$p.value
  bin_corr[i+1,"npvalue"] = -log(loop_cor$p.value)
}

##### Plotting
heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight=0.2) 
ggplot(bin_corr, aes(x=factor(bin, level = bin), y=value, fill = npvalue)) +
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_gradientn(colors = heat_colors, na.value = "red3",limits = c(0, 10),label = function(x) sprintf("%.1f", x))+
  theme(aspect.ratio = 1, 
        axis.text=element_text(size=11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust=1, vjust = 1),
        axis.title.y = element_text(hjust=0.5),
        axis.title=element_text(size=11),
        plot.title = element_text(hjust = 0.5),
        legend.direction="horizontal",
        legend.position=c(.5,.8),
        legend.title = element_text(size = 10),
        legend.key.width = unit(0.65, "cm"),
        legend.key.height = unit(0.5, "cm"))+
  guides(fill = guide_colourbar(title.position = "top", title.hjust=0.5))+
  labs(color = "Density", title = "", fill = expression("-log"["10"]~"\n"*italic(P)*"-value")) +
  xlab("Boundary Score Threshold") +
  ylab(bquote(atop(NA, atop(textstyle("Gene-Boundary Correlation"),
                            textstyle("(Spearman's "*rho*")"))))) +
  ylim(-0.15,0.15) +
  scale_x_discrete(labels= c("> 0", "", "","> 0.9","","","> 1.8","","","> 2.7","","","> 3.6","",""))

##### Random Promoters
bin_corr = data.frame()
width = 0.3
for (i in 0:14){
  loop_bin = loop_common[loop_common$TAD_Score1 >= 0+i*width | loop_common$TAD_Score2 >= 0+i*width,
                         c("Gap_Score", "mean_exp_rand")]
  loop_cor = cor.test(loop_bin$Gap_Score, loop_bin$mean_exp_rand, method = c("spearman"))
  bin_corr[i+1,"bin"] = paste("> ", 0+i*width, sep = "")
  bin_corr[i+1,"value"] = loop_cor$estimate[[1]]
  bin_corr[i+1,"pvalue"] = loop_cor$p.value
  bin_corr[i+1,"npvalue"] = -log(loop_cor$p.value)
}

heat_colors <- (colorRampPalette(c("gray95", "red3"))(50))
par(lheight=0.2) 
ggplot(bin_corr, aes(x=factor(bin, level = bin), y=value, fill = npvalue)) + #3.7 cutoff
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_gradientn(colors = heat_colors, na.value = "red3",limits = c(0, 10),label = function(x) sprintf("%.1f", x))+
  theme(aspect.ratio = 1, 
        axis.text=element_text(size=11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust=1, vjust = 1),
        axis.title.y = element_text(hjust=0.5),
        axis.title=element_text(size=11),
        plot.title = element_text(hjust = 0.5),
        legend.direction="horizontal",
        legend.position=c(.5,.25),
        legend.title = element_text(size = 10),
        legend.key.width = unit(0.65, "cm"),
        legend.key.height = unit(0.5, "cm"))+
  guides(fill = guide_colourbar(title.position = "top", title.hjust=0.5))+
  labs(color = "Density", title = "", fill = expression("-log"["10"]~"\n"*italic(P)*"-value")) +
  xlab("Boundary Score Threshold") +
  ylab(bquote(atop(NA, atop(textstyle("Gene-Boundary Correlation"),
                            textstyle("(Spearman's "*rho*")"))))) +
  ylim(-0.15,0.15) +
  scale_x_discrete(labels= c("> 0", "", "","> 0.9","","","> 1.8","","","> 2.7","","","> 3.6","","")) 
