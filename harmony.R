#BiocManager::install("harmony")
library(harmony)

immune.combined[['RNA']]@counts[1:5,1:5]
immune.combined[['RNA']]@data[1:5,1:5]
counts <- immune.combined[['RNA']]@counts
meta <- immune.combined@meta.data

pbmc <- CreateSeuratObject(counts = counts, project = "Harmony_integration")
pbmc
pbmcsca <- NormalizeData(pbmc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
head(pbmcsca@meta.data)
all(rownames(pbmcsca@meta.data) == rownames(meta)) #[1] TRUE
pbmcsca@meta.data$orig_seurat_clusters <- meta$orig_seurat_clusters
pbmcsca@meta.data$condition <- meta$condition
pbmcsca <- RunHarmony(pbmcsca, group.by.vars = "condition")
# Harmony converged after 5 iterations
pbmcsca <- RunUMAP(pbmcsca, reduction = "harmony", dims = 1:20)
pbmcsca <- FindNeighbors(pbmcsca, reduction = "harmony", dims = 1:20) %>% FindClusters()
head(pbmcsca@meta.data)
DimPlot(pbmcsca, group.by = "condition")
DimPlot(pbmcsca, group.by = "orig_seurat_clusters", label = T,
        split.by = "condition")
DimPlot(pbmcsca, group.by = "condition", label = F,
        split.by = "condition")

