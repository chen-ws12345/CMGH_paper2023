setwd("/Volumes/My_Passport/SMO_scData_analysis/")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Hmisc)
library(EnhancedVolcano)
library(magrittr)
library(tidyverse)
library(cowplot)

Cre_48hr.data <- Read10X(data.dir = "./Cre_48hr_PHX/filtered_feature_bc_matrix")
Cre_48hr <- CreateSeuratObject(counts = Cre_48hr.data, min.cells = 3, min.features = 200)
Cre_48hr #16425 features across 85348 samples
Cre_48hr[["percent.mt"]] <- PercentageFeatureSet(Cre_48hr, pattern = "^mt-")
head(Cre_48hr@meta.data, 5)
VlnPlot(Cre_48hr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Cre_48hr
Cre_48hr <- subset(Cre_48hr, subset = nFeature_RNA > 1750 & percent.mt <= 60)
Cre_48hr # 16425 features across 4357 samples
summary(Cre_48hr@meta.data)
head(Cre_48hr@meta.data)
Cre_48hr@meta.data$condition <- "Cre_48hr_PHX"

pbmc <- Cre_48hr
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

saveRDS(pbmc, file = "./savedRDS/Cre_48hr_PHX_4357_PC10.rds")

pbmc <- readRDS("/Users/tianyichen/Desktop/Regeneration_manuscript/Cre_48hr_PHX_4357_PC10.rds")
pbmc
head(pbmc@meta.data)
table(pbmc@meta.data$seurat_clusters)


## output for cellphoneDB ###
meta <- data.frame(Cell = rownames(pbmc@meta.data),
                   cell_type = pbmc@meta.data$seurat_clusters)
head(meta)
write.table(meta, file = "./SMO_48hr_PHX/cellphoneDB/Cre_48_meta.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
data <- pbmc[['RNA']]@data
data <- as.matrix(data)
data[1:5,1:5]
rownames(data) <- toupper(rownames(data))
dim(data) # [1] 16425  4357
write.table(data, file = "./SMO_48hr_PHX/cellphoneDB/Cre_48_counts.txt", sep = "\t",
            row.names = T, col.names = T, quote = F)
#### OVER #######
# read in cells - info from scanpy output
tmp <- read.table("./SMO_48hr_PHX/scvelo/Cre_48hr_PHX_files/smo_cre_phx48_cell_4357.txt", sep = "\t",
                  row.names = 1, header = T)
head(tmp)
dim(tmp)
#`%!in%` = Negate(`%in%`)
#tmp[rownames(tmp)[which(rownames(tmp) %!in% rownames(pbmc@meta.data))],]
#which(rownames(tmp) == "CAGGCCAAGGGCCAAT-1")
all(rownames(tmp) == rownames(pbmc@meta.data)) # [1] TRUE
#pbmc@meta.data$seurat_clusters <- paste0("C",pbmc@meta.data$seurat_clusters)
Idents(pbmc) <- pbmc@meta.data$seurat_clusters
head(pbmc@meta.data)
coldata <- data.frame(seurat_clusters = pbmc@meta.data$seurat_clusters,
                      row.names = rownames(pbmc@meta.data))
head(coldata)
write.table(coldata, file = "./SMO_48hr_PHX/scvelo/Cre_48hr_PHX_files/smo_cre_phx48_cell_4357_seurat_clusters.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)

df <- pbmc[["umap"]]@cell.embeddings
head(df)
write.table(df, file = "./SMO_48hr_PHX/scvelo/Cre_48hr_PHX_files/smo_cre_phx48_cell_4357_umap.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)
df <- pbmc[["pca"]]@cell.embeddings
head(df)
write.table(df, file = "./SMO_48hr_PHX/scvelo/Cre_48hr_PHX_files/smo_cre_phx48_cell_4357_pca50.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)
### finished output for scanpy and scVelo



pbmc[["RNA"]]@counts[1:5,1:5]
counts <- pbmc[["RNA"]]@counts
counts[1:5,1:5]
dim(counts) # [1] 14088  2955
write.table(counts, file = "~/Desktop/RISC/Cre_PHX_48hr_4357_scRNA_counts.txt",
            sep = "\t", row.names = T, col.names = T, quote = F)
meta <- pbmc@meta.data
write.table(meta, file = "~/Desktop/RISC/Cre_PHX_48hr_4357_metadata.txt",
            sep = "\t", row.names = T, col.names = T, quote = F)

all(colnames(pbmc[["RNA"]]@counts) == rownames(pbmc@meta.data)) # [1] TRUE

genes = c('Smo','Shh','Ihh')
genes = c('Tgfbr1','Tgfbr2','Tgfbr3')
genes = c('Mki67','Ccna2','Top2a','Ccnb1')
genes = c('Ihh')
VlnPlot(pbmc, features = genes, ncol = 2)
genes = c('Yap1','Wwtr1')
genes = c('Acsl3','Acsl4','Slc7a11')
genes = c('Fbp1')
genes = c('F2r','Thrb','Bcl2','Mcl1')
DotPlot(pbmc, features = genes) + RotatedAxis()

# "Des","Lrat","Cd52","Ptprc",
# "Apoa1","Hnf4a"
# 'Mki67','Ccna2','Top2a','Ccnb1'
# 'Axin2','Yap1','Sox9','Lgr5','Afp'
# 'Epcam','Krt7','Krt19','Cldn7','Hnf1b'

# Find markers
cluster.markers <- FindMarkers(pbmc, ident.1 = 'C5', min.pct = 0.25,
                               only.pos = TRUE)
head(cluster.markers, n = 5)
tail(cluster.markers)
dim(cluster.markers)
sigDEG <- cluster.markers %>% filter(p_val_adj <= 0.05) %>% filter(avg_log2FC >= 0.35)
dim(sigDEG)
head(sigDEG)
tail(sigDEG)
df <- data.frame(genes = rownames(sigDEG))
head(df)
write.table(df, file = "./sigDEs/C5_vs_others_padj005_035_marker_genes.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)

# return full list of genes
cluster.markers <- FindMarkers(pbmc, ident.1 = 'C5', #ident.2 = c(0,1,2,6), 
                               min.cells.group = 1, 
                               min.cells.feature = 1,
                               min.pct = 0,
                               logfc.threshold = 0,
                               only.pos = FALSE)
head(cluster.markers)
tail(cluster.markers)
dim(cluster.markers) # [1] 13752  5
## Create ranked list for GSEA analysis
cluster.markers$fcSign=sign(cluster.markers$avg_log2FC)
cluster.markers$logP=-log10(cluster.markers$p_val)
cluster.markers$metric=cluster.markers$logP/cluster.markers$fcSign
head(cluster.markers)
cluster.markers$Name <- rownames(cluster.markers)
y<-cluster.markers[,c("Name", "metric")]
head(y)
dim(y) 
ranked_list <- y[order(y$metric, decreasing = T),]
ranked_list$Gene.name <- str_to_upper(rownames(ranked_list))
head(ranked_list)
ranked_list <- ranked_list[,c("Gene.name", "metric")]
rownames(ranked_list) <- NULL
head(ranked_list)
ranked_list$metric[1:2] <- 282
PATH = "/Volumes/My_Passport/SMO_scData_analysis/SMO_48hr_PHX/gsea/"
write.table(ranked_list, file = paste0(PATH, "C5_vs_C0+1+2+6+3_full.rnk"), 
            sep = "\t", row.names = F, quote = F, col.names = F)

write.table(cluster.markers, file = "/Users/tianyichen/Desktop/SMO_scRNA/seurat_output/smo_luc_48hr_PHX/PC10/3_vs_0126.txt",
            sep = "\t", row.names = T, col.names = T, quote = F)

#### slingshot tradeSeq for Heps trajectory
pbmc <- readRDS("./Cre_48hr_PHX_4357_PC10.rds")
pbmc
head(pbmc@meta.data)
"%!in%" <- Negate("%in%")
pbmc@meta.data$seurat_clusters <- paste0("C",pbmc@meta.data$seurat_clusters)
hep <- subset(pbmc, seurat_clusters %!in% c("C4","C7","C8","C9"))
hep # 16425 features across 3745 samples
DimPlot(hep, reduction = "umap", label = T, 
        group.by = "seurat_clusters")
counts <- hep[["RNA"]]@counts
logcounts <- hep[['RNA']]@data
df <- hep[["umap"]]@cell.embeddings %>% as.data.frame()
head(df)
dim(df) # [1] 3745    2

cl <- hep@meta.data$seurat_clusters
head(cl)

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
sce <- SingleCellExperiment(assays = List(counts = counts, logcounts = logcounts))
pca_data <- prcomp(t(logcounts(sce)), rank=50)

PCA=pca_data$x
dim(PCA) # [1] 3745   50
head(PCA)
class(PCA)
df <- as.matrix(df)
head(df)
class(df)
all(rownames(PCA) == rownames(df)) # TRUE
reducedDims(sce) <- list(PCA=pca_data$x, UMAP=df)
colData(sce)$seurat_clusters <- cl

sce
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = 'seurat_clusters',
                 start.clus = 'C1')
sce

saveRDS(sce, file = "~/Desktop/smo_Cre_48hr_PHX.rds")
saveRDS(sce, file = "~/Desktop/updated_smo_Cre_48hr_PHX.rds")

sce <- readRDS("~/Desktop/smo_Cre_48hr_PHX.rds")
sce <- readRDS("~/Desktop/updated_smo_Cre_48hr_PHX.rds")

counts <- as.matrix(counts(sce))
counts[1:5,1:5]
crv <- SlingshotDataSet(sce)
crv

set.seed(5)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, 
                   nGenes = 200, verbose = T)

# Run on HARDAC - since the process takes really long!!!!
set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
head(pseudotime)
cellWeights <- slingCurveWeights(crv)
head(cellWeights)

sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = TRUE)

### Read in HARDAC output of updated sce

sce
table(rowData(sce)$tradeSeq$converged)
# FALSE  TRUE 
# 1163 12925
assoRes <- associationTest(sce, l2fc = log2(2))
head(assoRes)
dim(assoRes) # [1] 14088     4

assoRes_lineage <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
head(assoRes_lineage)
dim(assoRes_lineage) # [1] 14088    13

startRes_lineage <- startVsEndTest(sce, lineages = TRUE)
head(startRes_lineage)
oStart <- order(startRes_lineage$waldStat, decreasing = TRUE)
head(oStart)
sigGeneStart <- names(sce)[oStart[1]]
sigGeneStart
counts <- counts(sce)
plotSmoothers(sce, counts, gene = sigGeneStart)
tmp <- readRDS("~/Desktop/smo_Luc_48hr_PHX.rds")
crv <- SlingshotDataSet(tmp)
crv
plotGeneCount(crv, counts, gene = sigGeneStart)

endRes <- diffEndTest(sce, pairwise=TRUE)
head(endRes)
dim(endRes) # [1] 14088    15
endRes_sub1vs2 <- endRes[,c(4,5,6,13)]
head(endRes_sub1vs2)
o <- order(endRes_sub1vs2$waldStat_1vs2, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
sigGene
plotSmoothers(sce, counts, sigGene)
plotGeneCount(crv, counts, gene = sigGene)

endRes_sub1vs2$padj <- p.adjust(endRes_sub1vs2$pvalue_1vs2, "fdr")
head(endRes_sub1vs2)
sigGenes <- rownames(endRes_sub1vs2)[which(endRes_sub1vs2$padj <= 0.05)]
head(sigGenes)  
length(sigGenes)# [1] 1839
df <- data.frame(sigGenes1vs2 = sigGenes)
head(df)  
write.table(df, file = "~/Desktop/smo_luc_48PHX_tradeSeq_lin1vs2_diffEnd_sigGenes.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)  

patternRes <- patternTest(sce)
head(patternRes)
patternRes$padj <- p.adjust(patternRes$pvalue, "fdr")
genes <- rownames(patternRes)[which(patternRes$padj <= 0.05)]
head(genes)
length(genes) # 1204
plotSmoothers(sce, counts, gene = genes[1])
plotGeneCount(crv, counts, gene = genes[1])
df <- data.frame(sigGenes1vs2 = genes)
head(df)  
write.table(df, file = "~/Desktop/smo_luc_48PHX_tradeSeq_lin1vs2_diffPattern_sigGenes.txt",
            sep = "\t", row.names = F, col.names = T, quote = F) 

#UpSetR::upset(fromList(list(diffEnd = sigGenes, diffPattern = genes)))
idx <- intersect(sigGenes, genes)
length(idx) # [1] 1006

