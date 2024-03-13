###
### RNA-seq analysis of E11.5 clamped, unclamped, constrained, unconstrained, and amnion-retained samples
### 

library(DESeq2)
library(sva)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(IHW)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(reshape2) 
library(reshape) 
library(plyr)
library(stringr)
library(ggthemes) 

setwd("")
dir = getwd()

files = file.path(list.dirs("star_E11", recursive = F), "quant.sf")
names(files) = list.dirs("star_E11", recursive = F, full.names = F)
samples =read_csv("sampledata_E11.csv")
tx2gene = read_tsv("salmon_tx2gene.tsv", col_names = F)

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
txi_combat = ComBat_seq(txi$counts, batch = samples$batch)

dds <- DESeqDataSetFromMatrix(round(txi_combat), colData = samples, design = ~ condition)
dds = dds[rowSums(counts(dds) >= 10) >= 10,] 
dds$condition <- factor(dds$condition, levels = c("AMNION", "UNCOM", "CLAMP", "CONT"))

dds_norm <- rlog(dds, blind =  F)
mat <- assay(dds_norm)
mm <- model.matrix(~condition, colData(dds_norm))
mat <- limma::removeBatchEffect(mat, batch=dds_norm$batch, design=mm)
assay(dds_norm) <- mat

ggplot(data = dds_pca, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5)+
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text.x =element_text(colour="black"),
        axis.text.y =element_text(colour="black"),
        text=element_text(size=11), 
        axis.text=element_text(size=11)) + 
  ggforce::geom_mark_ellipse(aes(fill = group, color = group))+
  NoLegend() +
  ylab("PC2 (22% Variance)") +
  xlab("PC1 (58% Variance)") +
  scale_color_manual(values=c("#FF64B2", "grey", "#29B473", "#F6921E"))+ 
  scale_fill_manual(values=c("#FF64B2", "grey", "#29B473", "#F6921E")) #350x350

dds <- DESeq(dds)

##### Sample Correlation Heatmap
vsd <- vst(dds, blind =  F, nsub = 1000)
#vsd <- rlog(dds, blind =  F)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
levels(vsd$condition) = c("Amnion", "Control", "Clamped", "Constrained")
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$batch, sep="-")
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "YlGn")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, 
         cellheight=15, cellwidth = 15, border_color = "white") #rnaseq_cor, 550x450

##### Constraint Comparison
res = results(dds, contrast = c("condition", "CONT","UNCOM"))

ens.str <- rownames(res)
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$refseq <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="REFSEQ",
                     keytype="ENSEMBL",
                     multiVals="first")
res$foldchange = (2^res$log2FoldChange)

res_clamp = as.data.frame(res)
res_clamp = res_clamp[res_clamp$padj<0.05,]
res_clamp = na.omit(res_clamp)
res_clamp = res_clamp[order(-res_clamp$foldchange),]
write.csv(res_clamp, "E11_constraint.csv")
res_clamp = res_clamp[abs(res_clamp$log2FoldChange)>log2(1.25),]
write.csv(res_clamp, "res_const.csv")

deseq_counts = as.data.frame(counts(dds, normalized=TRUE))
deseq_counts$ensembl <- rownames(deseq_counts)
deseq_counts$symbol <- unlist(mapIds(org.Mm.eg.db,
                                     keys=deseq_counts$ensembl,
                                     column="SYMBOL",
                                     keytype="ENSEMBL",
                                     multiVals="first"))


deseq_counts = deseq_counts[complete.cases(deseq_counts),]
deseq_counts$ensembl = NULL
deseq_counts = aggregate(.~symbol, deseq_counts, sum)
rownames(deseq_counts) = deseq_counts$symbol
write.csv(deseq_counts, "deseq_counts_E11_constraint.csv")

##### Constraint vs Control Volcano Plot
res = results(dds, contrast = c("condition", "CONT","UNCOM"))
res <- lfcShrink(dds,contrast = c("condition", "CONT", "UNCOM"), res=res, type = 'normal')
ens.str <- rownames(res)
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

keyvals <- ifelse(
  res$log2FoldChange < -log2(1.25) & res$padj < 0.05, '#C51B7DFF',
  ifelse(res$log2FoldChange > log2(1.25) & res$padj < 0.05, '#4D9221FF',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == '#4D9221FF'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == '#C51B7DFF'] <- 'low'

EnhancedVolcano(res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
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
                xlim = c(-2.5,2.5),
                ylim = c(0, 290),
                boxedLabels = T,
                drawConnectors = TRUE,
                arrowheads = F) +
  theme(aspect.ratio = 1, 
        axis.text=element_text(size=11, color = "black"),
        axis.title=element_text(size=11))+
  xlab(expression("log"["2"]~"Fold Change")) +
  ylab(expression("log"["10"]~italic(P))) +
  NoLegend() 

##### Amnion-Retained vs Control Volcano Plot
res = results(dds, contrast = c("condition", "CONT","UNCOM"))
res <- lfcShrink(dds,contrast = c("condition", "CONT", "UNCOM"), res=res, type = 'normal')
ens.str <- rownames(res)
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

keyvals <- ifelse(
  res$log2FoldChange < -log2(1.25) & res$padj < 0.05, '#C51B7DFF',
  ifelse(res$log2FoldChange > log2(1.25) & res$padj < 0.05, '#4D9221FF',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == '#4D9221FF'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == '#C51B7DFF'] <- 'low'

EnhancedVolcano(res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
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
                xlim = c(-2.5,2.5),
                ylim = c(0, 290),
                boxedLabels = T,
                drawConnectors = TRUE,
                arrowheads = F) +
  theme(aspect.ratio = 1, 
        axis.text=element_text(size=11, color = "black"),
        axis.title=element_text(size=11))+
  xlab(expression("log"["2"]~"Fold Change")) +
  ylab(expression("log"["10"]~italic(P))) +
  NoLegend() 


##### Individual Gene Plotting
##### Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summarized
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- reshape::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

##### HoxA Gene Panel
rownames(deseq_counts) = deseq_counts$symbol
hoxa_diff= deseq_counts[c("Hoxa1", "Hoxa2", "Hoxaas2","Hoxa3", "Hoxaas3","Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7", "Hoxa9", "Hoxa10", "Hoxa11", "Hoxa11os", "Hoxa13", "Hoxd3", "Hoxd4", "Hoxd8","Hoxd9","Hoxd10", "Hoxd11", "Hoxd12","Hoxd13"),]
hoxa_diff <- melt(hoxa_diff, id = c("symbol"), variable.name = "variable")
hoxa_diff = data_summary(hoxa_diff, varname = "value", groupnames = c("variable", "symbol"))
hoxa_diff$zscore <- ave(hoxa_diff$value, hoxa_diff$symbol, FUN=scale) #zscore calculation

geneorder = c("Hoxa1", "Hoxa2", "Hoxaas2","Hoxa3", "Hoxaas3","Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7", "Hoxa9", "Hoxa10", "Hoxa11", "Hoxa11os", "Hoxa13", "Hoxd3", "Hoxd4", "Hoxd8","Hoxd9","Hoxd10", "Hoxd11", "Hoxd12","Hoxd13")
variableorder = c("Amnion","Constrained", "Clamped", "Control" )
palette <- brewer.pal(11, "Purples") 
heat_colors <- (colorRampPalette(palette)(50))

ggplot(hoxa_diff, aes(x=factor(symbol,level = geneorder), y=factor(variable,level = variableorder), fill=zscore))+
  theme_tufte(base_family="Helvetica")+
  geom_tile(color="white", size=0.2)+
  scale_fill_gradientn(colors = heat_colors)+
  coord_equal()+
  theme(plot.title=element_text(size = 11),
        axis.text.x = element_text(face = "italic",angle = 45, hjust=1,vjust = 1, colour = "black"),
        axis.text.y = element_text( color = "black"),
        axis.ticks=element_blank(),
        axis.text=element_text(size=11),
        strip.text.x = element_text(size = 6))+
  labs(fill = expression("Z-Score"))+
  xlab("")+
  ylab("")

##### HoxD Gene Panel
hoxd_diff= deseq_counts[grepl("Hoxd", rownames(deseq_counts)),]
hoxd_diff <- melt(hoxd_diff, id = c("symbol"), variable.name = "variable")
hoxd_diff = data_summary(hoxd_diff, varname = "value", groupnames = c("variable", "symbol"))
hoxd_diff$zscore <- ave(hoxd_diff$value, hoxd_diff$symbol, FUN=scale) #zscore calculation

geneorder = c("Hoxd3", "Hoxd4", "Hoxd8","Hoxd9","Hoxd10", "Hoxd11", "Hoxd12","Hoxd13" )
variableorder = c("Amnion","Constrained", "Clamped", "Control" )
palette <- brewer.pal(11, "Purples") 
heat_colors <- (colorRampPalette(palette)(50))

ggplot(hoxd_diff, aes(x=factor(symbol,level = geneorder), y=factor(variable,level = variableorder), fill=zscore))+
  theme_tufte(base_family="Helvetica")+
  geom_tile(color="white", size=0.25)+
  scale_fill_gradientn(colors = heat_colors)+
  coord_equal()+
  theme(plot.title=element_text(size = 11),
        axis.text.x = element_text(face = "italic",angle = 45, hjust=1,vjust = 1, colour = "black"),
        axis.text.y = element_text( color = "black"),
        axis.ticks=element_blank(),
        axis.text=element_text(size=11),
        strip.text.x = element_text(size = 6))+
  labs(fill = expression("Z-Score"))+
  xlab("")+
  ylab("")


##### Mechanotransduction Genes
rownames(deseq_counts) = deseq_counts$symbol
hoxa_diff= deseq_counts[c("Sfrp1","Gpc3", "Celsr2", "Wnt3","Wnt4", "Wnt6","Wnt7a","Wnt7b","Wnt10a","Wnt10b", "Tgfb1","Ppp2r2b", "Wwc1", "Cdh1", "Ajuba", "Llgl2"),]
hoxa_diff <- melt(hoxa_diff, id = c("symbol"), variable.name = "variable")
hoxa_diff = data_summary(hoxa_diff, varname = "value", groupnames = c("variable", "symbol"))
hoxa_diff$zscore <- ave(hoxa_diff$value, hoxa_diff$symbol, FUN=scale) #zscore calculation

geneorder = c("Sfrp1","Gpc3", "Celsr2", "Wnt3","Wnt4", "Wnt6","Wnt7a","Wnt7b","Wnt10a","Wnt10b", "Tgfb1","Ppp2r2b", "Wwc1", "Cdh1", "Ajuba", "Llgl2")
variableorder = c("Amnion","Constrained", "Clamped", "Control" )
palette <- brewer.pal(11, "Purples") 
heat_colors <- (colorRampPalette(palette)(50))

ggplot(hoxa_diff, aes(x=factor(symbol,level = geneorder), y=factor(variable,level = variableorder), fill=zscore))+
  theme_tufte(base_family="Helvetica")+
  geom_tile(color="white", size=0.2)+
  scale_fill_gradientn(colors = heat_colors)+
  coord_equal()+
  theme(plot.title=element_text(size = 11),
        axis.text.x = element_text(face = "italic",angle = 45, hjust=1,vjust = 1, colour = "black"),
        axis.text.y = element_text( color = "black"),
        axis.ticks=element_blank(),
        axis.text=element_text(size=11),
        strip.text.x = element_text(size = 6))+
  labs(fill = expression("Z-Score"))+
  xlab("")+
  ylab("")
