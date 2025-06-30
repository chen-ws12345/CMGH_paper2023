#BiocManager::install("scmap")
#devtools::install_github("pcahan1/singleCellNet")

library(SingleCellExperiment)
library(Seurat)
library(scmap)

setwd("/Volumes/My_Passport/SMO_scData_analysis/SMO_48hr_PHX/savedRDS")
"%!in%" <- Negate("%in%")
luc <- readRDS("./Luc_48hr_PHX_2955_PC10.rds")
luc <- subset(luc, seurat_clusters %in% c("3","4","7"))
luc # 383 samples
logcounts <- luc[['RNA']]@data
logcounts[1:5,1:5]
meta <- luc@meta.data
head(meta)
luc_ann <- data.frame(row.names = rownames(meta),
                      cell_type1 = paste0("L",meta$seurat_clusters))
head(luc_ann)

luc_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)),
                                colData = luc_ann)
rowData(luc_sce)$feature_symbol <- rownames(luc_sce)
luc_sce <- luc_sce[!duplicated(rownames(luc_sce)),]
luc_sce

cre <- readRDS("./Cre_48hr_PHX_4357_PC10.rds")
cre <- subset(cre, seurat_clusters %in% c("3","5"))
cre # 722 samples
logcounts <- cre[['RNA']]@data
logcounts[1:5,1:5]
meta <- cre@meta.data
head(meta)
cre_ann <- data.frame(row.names = rownames(meta),
                      cell_type1 = paste0("C",meta$seurat_clusters))
head(cre_ann)

cre_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)),
                                colData = cre_ann)
rowData(cre_sce)$feature_symbol <- rownames(cre_sce)
cre_sce <- cre_sce[!duplicated(rownames(cre_sce)),]
cre_sce

# Feature Selection
luc_sce <- selectFeatures(luc_sce, suppress_plot = FALSE)
table(rowData(luc_sce)$scmap_features)
cre_sce <- selectFeatures(cre_sce, suppress_plot = FALSE)
table(rowData(cre_sce)$scmap_features)

# scmap - cluster
cre_sce <- indexCluster(cre_sce)
head(metadata(cre_sce)$scmap_cluster_index)

# Projection
scmapCluster_results <- scmapCluster(
  projection = luc_sce,
  index_list = list(
    res = metadata(cre_sce)$scmap_cluster_index
  )
)
head(scmapCluster_results$scmap_cluster_labs)
dim(scmapCluster_results$scmap_cluster_labs) # [1] 4357    1

# Visulization
dev.off()
plot(
  getSankey(
    colData(luc_sce)$cell_type1,
    scmapCluster_results$scmap_cluster_labs[,'res'],
    plot_height = 400
  )
)

# scmap - cell
set.seed(1)
cre_sce <- indexCell(cre_sce)
scmapCell_results <- scmapCell(
  luc_sce,
  list(
    res = metadata(cre_sce)$scmap_cell_index
  )
)
scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results,
  list(
    as.character(colData(cre_sce)$cell_type1)
  )
)
head(scmapCell_clusters$scmap_cluster_labs)
dev.off()
plot(
  getSankey(
    colData(luc_sce)$cell_type1,
    scmapCell_clusters$scmap_cluster_labs[,'res'],
    plot_height = 400
  )
)


