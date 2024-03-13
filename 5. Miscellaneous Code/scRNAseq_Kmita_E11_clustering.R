###
### scRNA-seq clustering code of E11.5 whole limb
### https://doi.org/10.3390/jdb8040031


library(Seurat)
library(ggplot2)
library(data.table)

setwd("./")
dir <- getwd()

load("GSM4811355_E115.data.R")

limb_E11 <- CreateSeuratObject(counts = E115.data)
VlnPlot(
  limb_E11,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 3,
  pt.size = 0.1
)
limb_E11 <- subset(limb_E11, subset = nFeature_RNA > 1000 &
  nCount_RNA > 500 & nCount_RNA < 60000)

limb_E11 <- NormalizeData(limb_E11)
limb_E11 <-
  FindVariableFeatures(limb_E11, selection.method = "vst", nfeatures = 10000)
all.genes <- rownames(limb_E11)
limb_E11 <- ScaleData(limb_E11, features = all.genes)
limb_E11 <-
  RunPCA(limb_E11, features = VariableFeatures(object = limb_E11))

limb_E11 <- FindNeighbors(limb_E11, dims = 1:20)
limb_E11 <- FindClusters(limb_E11, resolution = 3)
limb_E11 <- RunUMAP(limb_E11, dims = 1:20, min.dist = 1)

DimPlot(limb_E11,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  )

clusterlength <- 40
clusters <- c(1, 2, 6, 7, 12, 15, 17, 20, 22, 23, 26, 27)
label1 <- "Distal"
label2 <- "Other"

markerlabel <-
  function(seuratdata,
           clusters,
           label1,
           label2,
           clusterlength) {
    clusteridents <- data.frame()
    for (x in 1:clusterlength) {
      if ((x - 1) %in% clusters == TRUE) {
        clusteridents[x, 1] <- label1
      } else {
        clusteridents[x, 1] <- label2
      }
    }
    clusteridents <- as.list(clusteridents$V1)
    names(clusteridents) <- levels(seuratdata)
    seuratdata <- RenameIdents(seuratdata, clusteridents)
    print(levels(seuratdata))
    return(seuratdata)
  }

# Distal_Mesenchyme (Hoxa13) 1,2,6,7,12,15,17,20,22,23,26,27
limb_Distal_Mesenchyme <- limb_E11
limb_Distal_Mesenchyme <- markerlabel(
  limb_Distal_Mesenchyme,
  c(1, 2, 6, 7, 12, 15, 17, 20, 22, 23, 26, 27),
  "Distal_Mesenchyme",
  "Other",
  40
)
Distal_Mesenchyme_markers <- FindMarkers(limb_Distal_Mesenchyme, ident.1 = "Distal_Mesenchyme")
Distal_Mesenchyme_markers <- Distal_Mesenchyme_markers[Distal_Mesenchyme_markers$p_val_adj <
  0.05, ]
Distal_Mesenchyme_markers <- Distal_Mesenchyme_markers[Distal_Mesenchyme_markers$avg_log2FC >
  0, ]
write.csv(
  Distal_Mesenchyme_markers,
  "scRNA_E115_Distal_Mesenchyme_markers.csv"
)

DimPlot(
  limb_Distal_Mesenchyme,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Proximal_Mesenchyme (Hoxa11) 0,3,5,8,9,10,11,13,14,16,18,19,23,31,32
limb_Proximal_Mesenchyme <- limb_E11
limb_Proximal_Mesenchyme <- markerlabel(
  limb_Proximal_Mesenchyme,
  c(0, 3, 5, 8, 9, 10, 11, 13, 14, 16, 18, 19, 23, 31, 32),
  "Proximal_Mesenchyme",
  "Other",
  40
)
Proximal_Mesenchyme_markers <- FindMarkers(limb_Proximal_Mesenchyme, ident.1 = "Proximal_Mesenchyme")
Proximal_Mesenchyme_markers <- Proximal_Mesenchyme_markers[Proximal_Mesenchyme_markers$p_val_adj <
  0.05, ]
Proximal_Mesenchyme_markers <- Proximal_Mesenchyme_markers[Proximal_Mesenchyme_markers$avg_log2FC >
  0, ]
write.csv(
  Proximal_Mesenchyme_markers,
  "scRNA_E115_Proximal_Mesenchyme_markers.csv"
)

DimPlot(
  limb_Proximal_Mesenchyme,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("springgreen3", "grey87")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Anterior Mesenchyme (Pax9) 5,9,11,23
limb_Anterior_Mesenchyme <- limb_E11
limb_Anterior_Mesenchyme <- markerlabel(
  limb_Anterior_Mesenchyme,
  c(5, 9, 11, 23),
  "Anterior_Mesenchyme",
  "Other",
  40
)
Anterior_Mesenchyme_markers <- FindMarkers(limb_Anterior_Mesenchyme, ident.1 = "Anterior_Mesenchyme")
Anterior_Mesenchyme_markers <- Anterior_Mesenchyme_markers[Anterior_Mesenchyme_markers$p_val_adj <
  0.05, ]
Anterior_Mesenchyme_markers <- Anterior_Mesenchyme_markers[Anterior_Mesenchyme_markers$avg_log2FC >
  0, ]
write.csv(
  Anterior_Mesenchyme_markers,
  "scRNA_E115_Anterior_Mesenchyme_markers.csv"
)

DimPlot(
  limb_Anterior_Mesenchyme,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Posterior Mesenchyme (Shh) 1,2,6,12
limb_Posterior_Mesenchyme <- limb_E11
limb_Posterior_Mesenchyme <- markerlabel(
  limb_Posterior_Mesenchyme,
  c(1, 2, 6, 12),
  "Posterior_Mesenchyme",
  "Other",
  40
)
Posterior_Mesenchyme_markers <- FindMarkers(limb_Posterior_Mesenchyme, ident.1 = "Posterior_Mesenchyme")
Posterior_Mesenchyme_markers <- Posterior_Mesenchyme_markers[Posterior_Mesenchyme_markers$p_val_adj <
  0.05, ]
Posterior_Mesenchyme_markers <- Posterior_Mesenchyme_markers[Posterior_Mesenchyme_markers$avg_log2FC >
  0, ]
write.csv(
  Posterior_Mesenchyme_markers,
  "scRNA_E115_Posterior_Mesenchyme_markers.csv"
)

DimPlot(
  limb_Posterior_Mesenchyme,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Central Mesenchyme (Grem1) 1,14,15,20,22
limb_Central_Mesenchyme <- limb_E11
limb_Central_Mesenchyme <- markerlabel(
  limb_Central_Mesenchyme,
  c(1, 14, 15, 20, 22),
  "Central_Mesenchyme",
  "Other",
  40
)
Central_Mesenchyme_markers <- FindMarkers(limb_Central_Mesenchyme, ident.1 = "Central_Mesenchyme")
Central_Mesenchyme_markers <- Central_Mesenchyme_markers[Central_Mesenchyme_markers$p_val_adj <
  0.05, ]
Central_Mesenchyme_markers <- Central_Mesenchyme_markers[Central_Mesenchyme_markers$avg_log2FC >
  0, ]
write.csv(
  Central_Mesenchyme_markers,
  "scRNA_E115_Central_Mesenchyme_markers.csv"
)

DimPlot(
  limb_Central_Mesenchyme,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Peripheral Mesenchyme (Msx2) 1,2,5,6,9,11,12,15,23
limb_Peripheral_Mesenchyme <- limb_E11
limb_Peripheral_Mesenchyme <- markerlabel(
  limb_Peripheral_Mesenchyme,
  c(1, 2, 5, 6, 9, 11, 12, 15, 23),
  "Peripheral_Mesenchyme",
  "Other",
  40
)
Peripheral_Mesenchyme_markers <- FindMarkers(limb_Peripheral_Mesenchyme, ident.1 = "Peripheral_Mesenchyme")
Peripheral_Mesenchyme_markers <- Peripheral_Mesenchyme_markers[Peripheral_Mesenchyme_markers$p_val_adj <
  0.05, ]
Peripheral_Mesenchyme_markers <- Peripheral_Mesenchyme_markers[Peripheral_Mesenchyme_markers$avg_log2FC >
  0, ]
write.csv(
  Peripheral_Mesenchyme_markers,
  "scRNA_E115_Peripheral_Mesenchyme_markers.csv"
)

DimPlot(
  limb_Peripheral_Mesenchyme,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Condensates (Sox9) 3,5,7,8,11,16,17,19,25,26,27,29
limb_Condensates <- limb_E11
limb_Condensates <- markerlabel(
  limb_Condensates,
  c(3, 5, 7, 8, 11, 16, 17, 19, 25, 26, 27, 29),
  "Condensates",
  "Other",
  40
)
Condensates_markers <- FindMarkers(limb_Condensates, ident.1 = "Condensates")
Condensates_markers <- Condensates_markers[Condensates_markers$p_val_adj <
  0.05, ]
Condensates_markers <- Condensates_markers[Condensates_markers$avg_log2FC >
  0, ]
write.csv(Condensates_markers, "scRNA_E115_Condensates_markers.csv")

DimPlot(
  limb_Condensates,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Apical Ectodermal Ridge (Fgf8) 35
limb_Apical_Ectodermal_Ridge <- limb_E11
limb_Apical_Ectodermal_Ridge <- markerlabel(
  limb_Apical_Ectodermal_Ridge,
  c(35),
  "Apical_Ectodermal_Ridge",
  "Other",
  40
)
Apical_Ectodermal_Ridge_markers <- FindMarkers(limb_Apical_Ectodermal_Ridge, ident.1 = "Apical_Ectodermal_Ridge")
Apical_Ectodermal_Ridge_markers <- Apical_Ectodermal_Ridge_markers[Apical_Ectodermal_Ridge_markers$p_val_adj <
  0.05, ]
Apical_Ectodermal_Ridge_markers <- Apical_Ectodermal_Ridge_markers[Apical_Ectodermal_Ridge_markers$avg_log2FC >
  0, ]
write.csv(
  Apical_Ectodermal_Ridge_markers,
  "scRNA_E115_Apical_Ectodermal_Ridge_markers.csv"
)

DimPlot(
  limb_Apical_Ectodermal_Ridge,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Dorsal/Ventral Ectoderm (En1/Wnt7a) 21,34,35
limb_Dorsal_Ventral_Ectoderm <- limb_E11
limb_Dorsal_Ventral_Ectoderm <- markerlabel(
  limb_Dorsal_Ventral_Ectoderm,
  c(21, 34, 35),
  "Dorsal_Ventral_Ectoderm",
  "Other",
  40
)
Dorsal_Ventral_Ectoderm_markers <- FindMarkers(limb_Dorsal_Ventral_Ectoderm, ident.1 = "Dorsal_Ventral_Ectoderm")
Dorsal_Ventral_Ectoderm_markers <- Dorsal_Ventral_Ectoderm_markers[Dorsal_Ventral_Ectoderm_markers$p_val_adj <
  0.05, ]
Dorsal_Ventral_Ectoderm_markers <- Dorsal_Ventral_Ectoderm_markers[Dorsal_Ventral_Ectoderm_markers$avg_log2FC >
  0, ]
write.csv(
  Dorsal_Ventral_Ectoderm_markers,
  "scRNA_E115_Dorsal_Ventral_Ectoderm_markers.csv"
)

DimPlot(
  limb_Dorsal_Ventral_Ectoderm,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Early Myogenic Cells (Myf5/Pax3) 28
limb_Early_Myogenic_Cells <- limb_E11
limb_Early_Myogenic_Cells <- markerlabel(
  limb_Early_Myogenic_Cells,
  c(28),
  "Early_Myogenic_Cells",
  "Other",
  40
)
Early_Myogenic_Cells_markers <- FindMarkers(limb_Early_Myogenic_Cells, ident.1 = "Early_Myogenic_Cells")
Early_Myogenic_Cells_markers <- Early_Myogenic_Cells_markers[Early_Myogenic_Cells_markers$p_val_adj <
  0.05, ]
Early_Myogenic_Cells_markers <- Early_Myogenic_Cells_markers[Early_Myogenic_Cells_markers$avg_log2FC >
  0, ]
write.csv(
  Early_Myogenic_Cells_markers,
  "scRNA_E115_Early_Myogenic_Cells_markers.csv"
)

DimPlot(
  limb_Early_Myogenic_Cells,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Late Myogenic Cells (Myod1) 37
limb_Late_Myogenic_Cells <- limb_E11
limb_Late_Myogenic_Cells <- markerlabel(
  limb_Late_Myogenic_Cells,
  c(37),
  "Late_Myogenic_Cells",
  "Other",
  40
)
Late_Myogenic_Cells_markers <- FindMarkers(limb_Late_Myogenic_Cells, ident.1 = "Late_Myogenic_Cells")
Late_Myogenic_Cells_markers <- Late_Myogenic_Cells_markers[Late_Myogenic_Cells_markers$p_val_adj <
  0.05, ]
Late_Myogenic_Cells_markers <- Late_Myogenic_Cells_markers[Late_Myogenic_Cells_markers$avg_log2FC >
  0, ]
write.csv(
  Late_Myogenic_Cells_markers,
  "scRNA_E115_Late_Myogenic_Cells_markers.csv"
)

DimPlot(
  limb_Late_Myogenic_Cells,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

# Endothelium (Pecam1) 36,39
limb_Endothelium <- limb_E11
limb_Endothelium <- markerlabel(limb_Endothelium, c(36, 39), "Endothelium", "Other", 40)
Endothelium_markers <- FindMarkers(limb_Endothelium, ident.1 = "Endothelium")
Endothelium_markers <- Endothelium_markers[Endothelium_markers$p_val_adj <
  0.05, ]
Endothelium_markers <- Endothelium_markers[Endothelium_markers$avg_log2FC >
  0, ]
write.csv(Endothelium_markers, "scRNA_E115_Endothelium_markers.csv")

DimPlot(
  limb_Endothelium,
  reduction = "umap",
  pt.size = 0.5,
  label = TRUE,
  cols = c("grey87", "springgreen3")
) +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = "pt"
    )
  ) +
  NoLegend()

spatial_markers <-
  qpcR:::cbind.na(
    rownames(Distal_Mesenchyme_markers),
    rownames(Proximal_Mesenchyme_markers),
    rownames(Anterior_Mesenchyme_markers),
    rownames(Posterior_Mesenchyme_markers),
    rownames(Central_Mesenchyme_markers),
    rownames(Peripheral_Mesenchyme_markers),
    rownames(Condensates_markers),
    rownames(Apical_Ectodermal_Ridge_markers),
    rownames(Dorsal_Ventral_Ectoderm_markers),
    rownames(Early_Myogenic_Cells_markers),
    rownames(Late_Myogenic_Cells_markers),
    rownames(Endothelium_markers)
  )

colnames(spatial_markers) <- c(
  "Distal_Mesenchyme_markers",
  "Proximal_Mesenchyme_markers",
  "Anterior_Mesenchyme_markers",
  "Posterior_Mesenchyme_markers",
  "Central_Mesenchyme_markers",
  "Peripheral_Mesenchyme_markers",
  "Condensates_markers",
  "Apical_Ectodermal_Ridge_markers",
  "Dorsal_Ventral_Ectoderm_markers",
  "Early_Myogenic_Cells_markers",
  "Late_Myogenic_Cells_markers",
  "Endothelium_markers"
)

library(data.table)
long <- melt(spatial_markers, variable.name = "area")
long[, 1] <- NULL
long <- na.omit(long)

saveRDS(long, "long_E11_spatial.rds")
