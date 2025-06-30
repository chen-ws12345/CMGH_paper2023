setwd("/Volumes/My_Passport/SMO_scData_analysis/SMO_48hr_PHX/")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Hmisc)
library(EnhancedVolcano)
library(magrittr)
library(tidyverse)
library(cowplot)

# Load the dataset
cre_48hr <- readRDS("./savedRDS/Cre_48hr_PHX_4357_PC10.rds")
cre_48hr
head(cre_48hr@meta.data)
cre_48hr@meta.data$orig_seurat_clusters <- paste0("C",cre_48hr@meta.data$seurat_clusters)

luc_48hr <- readRDS("./savedRDS/Luc_48hr_PHX_2955_PC10.rds")
luc_48hr
head(luc_48hr@meta.data)
luc_48hr@meta.data$orig_seurat_clusters <- paste0("L",luc_48hr@meta.data$seurat_clusters)

# Combine datasets
ifnb.list <- lapply(X = list(luc_48hr,cre_48hr), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# Perform integration
# Vision 3.2
ifnb.list <- list(luc_48hr,cre_48hr)
ifnb.list
features <- SelectIntegrationFeatures(object.list = ifnb.list) # select features that are repeatedly variable across datasets for integration
length(features) # 2000
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features) #dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors) # dims = 1:20

# Perform an integrated analysis
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# UMAP and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# Visulization
head(immune.combined@meta.data)
tail(immune.combined@meta.data)
DimPlot(immune.combined, reduction = "umap", group.by = "condition")
DimPlot(immune.combined, reduction = "umap", label = T)

DimPlot(immune.combined, reduction = "umap", group.by = "condition", 
        split.by = "condition")
DimPlot(immune.combined, reduction = "umap", split.by = "condition", 
        group.by = "orig_seurat_clusters",label = TRUE)

saveRDS(immune.combined, file = "./savedRDS/Cre&Luc_48hr_PHX_7312_seurat_integrated.rds")

theme_set(theme_cowplot())
# Find Markers
DefaultAssay(immune.combined) <- "RNA"
genes = c('Mki67','Top2a','Ccnb1','Foxm1','Ccna2',
          'Cdk1')
VlnPlot(immune.combined, features = genes)
VlnPlot(immune.combined, features = 'Axin2', split.by = "condition")
DotPlot(immune.combined, features = genes) + RotatedAxis()
DotPlot(immune.combined, features = genes, 
        cols = c("blue", "red"), dot.scale = 7, 
        split.by = "condition") + 
  RotatedAxis()

nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

# Identify differential expressed genes across conditions
immune.combined$celltype.condition <- paste(Idents(immune.combined), immune.combined$condition, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.condition"
head(immune.combined@meta.data)
markers <- FindMarkers(immune.combined, ident.1 = "1_Cre_0hr", ident.2 = "1_Luc_0hr", verbose = FALSE)
head(markers, n = 10)
dim(markers)
genes = "Cyp2f2"
FeaturePlot(immune.combined, features = genes, split.by = "condition",
            cols = c("grey", "red"))
VlnPlot(immune.combined, features = genes, split.by = "condition", 
        group.by = "celltype", pt.size = 0, combine = FALSE)




