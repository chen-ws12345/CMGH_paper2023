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
library(homologene)

luc_48hr.data <- Read10X(data.dir = "/Volumes/My_Passport/SMO_scData_analysis/Luc_48hr_PHX/filtered_feature_bc_matrix")
luc_48hr <- CreateSeuratObject(counts = luc_48hr.data, min.cells = 3, min.features = 200)
luc_48hr
luc_48hr[["percent.mt"]] <- PercentageFeatureSet(luc_48hr, pattern = "^mt-")
head(luc_48hr@meta.data, 5)
VlnPlot(luc_48hr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
luc_48hr
luc_48hr <- subset(luc_48hr, subset = percent.mt <= 60)
luc_48hr # 2955 cells
head(luc_48hr@meta.data)
luc_48hr@meta.data$condition <- "Luc_48hr_PHX"

pbmc <- luc_48hr
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
DimPlot(pbmc, reduction = "umap", label = T, 
        group.by = "seurat_clusters")

saveRDS(pbmc, file = "./Luc_48hr_PHX_2955_PC10.rds")

pbmc <- readRDS("/Users/tianyichen/Desktop/Regeneration_manuscript/Luc_48hr_PHX_2955_PC10.rds")
pbmc
head(pbmc@meta.data)
table(pbmc@meta.data$seurat_clusters)

pbmc@meta.data$orig_seurat_clusters <- paste0("L",pbmc@meta.data$seurat_clusters)
genes = c('Fbp1')
DotPlot(pbmc, features = genes, group.by = "orig_seurat_clusters") + RotatedAxis()

## output for cellphoneDB ###
meta <- data.frame(Cell = rownames(pbmc@meta.data),
                   cell_type = pbmc@meta.data$orig_seurat_clusters)
head(meta)
write.table(meta, file = "./SMO_48hr_PHX/cellphoneDB/Luc_48_meta.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
data <- pbmc[['RNA']]@data
data <- as.matrix(data)
data[1:5,1:5]
rownames(data) <- toupper(rownames(data))
dim(data) # [1] 14088  2955
write.table(data, file = "./SMO_48hr_PHX/cellphoneDB/Luc_48_counts.txt", sep = "\t",
            row.names = T, col.names = T, quote = F)
#### OVER #######

pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:30)
DimPlot(pbmc, reduction = "umap", label = T, 
        group.by = "orig_seurat_clusters")

FeatureScatter(pbmc, feature1 = "Lrat", feature2 = "Alb")


# read in cells - info from scanpy output
tmp <- read.table("./SMO_48hr_PHX/scvelo/Luc_48hr_PHX_files/smo_luc_phx48_cell_2955.txt", sep = "\t",
                  row.names = 1, header = T)
head(tmp)
all(rownames(tmp) == rownames(pbmc@meta.data)) # [1] TRUE
head(pbmc@meta.data)
coldata <- data.frame(seurat_clusters = pbmc@meta.data$orig_seurat_clusters,
                      row.names = rownames(pbmc@meta.data))
head(coldata)
write.table(coldata, file = "./SMO_48hr_PHX/scvelo/Luc_48hr_PHX_files/smo_luc_phx48_cell_2955_seurat_clusters.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)

df <- pbmc[["umap"]]@cell.embeddings
head(df)
write.table(df, file = "./SMO_48hr_PHX/scvelo/Luc_48hr_PHX_files/smo_luc_phx48_cell_2955_umap.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)
df <- pbmc[["pca"]]@cell.embeddings
head(df)
write.table(df, file = "./SMO_48hr_PHX/scvelo/Luc_48hr_PHX_files/smo_luc_phx48_cell_2955_pca50.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)
### finished output for scanpy and scVelo


all(colnames(pbmc[["RNA"]]@counts) == rownames(pbmc@meta.data)) # [1] TRUE
pbmc[["RNA"]]@counts[1:5,1:5]
counts <- pbmc[["RNA"]]@counts
counts[1:5,1:5]
dim(counts) # [1] 14088  2955
write.table(counts, file = "~/Desktop/RISC/Luc_PHX_48hr_2955_scRNA_counts.txt",
            sep = "\t", row.names = T, col.names = T, quote = F)
meta <- pbmc@meta.data
write.table(meta, file = "~/Desktop/RISC/Luc_PHX_48hr_2955_metadata.txt",
            sep = "\t", row.names = T, col.names = T, quote = F)

genes = c("Smo")
VlnPlot(pbmc, features = genes)

# 'Mki67','Ccna2','Top2a','Ccnb1'
# 'Axin2','Yap1','Sox9','Lgr5','Afp'
# 'Epcam','Krt7','Krt19','Cldn7','Hnf1b'

# Find markers
Idents(pbmc) <- pbmc@meta.data$orig_seurat_clusters
cluster.markers <- FindMarkers(pbmc, ident.1 = "L4", ident.2 = "L3", min.pct = 0.25,
                               logfc.threshold = 0.35, min.diff.pct = 0.2,
                               only.pos = TRUE)
head(cluster.markers, n = 5)
tail(cluster.markers)
dim(cluster.markers)
df <- data.frame(genes = rownames(cluster.markers))
head(df)
PATH = "/Volumes/My_Passport/SMO_scData_analysis/SMO_48hr_PHX/sigDEs/"
write.table(df, file = paste0(PATH, "L4_vs_L3_markers.txt"),
            sep = "\t", row.names = F, col.names = F, quote = F)

# return full list of genes
cluster.markers <- FindMarkers(pbmc, ident.1 = c("L2","L6"), ident.2 = c("L0","L1"), 
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
PATH = "/Volumes/My_Passport/SMO_scData_analysis/SMO_48hr_PHX/gsea/"
write.table(ranked_list, file = paste0(PATH, "L2+6_vs_L0+1_full.rnk"), 
            sep = "\t", row.names = F, quote = F, col.names = F)

###############################################################################################
############### Calculate Gene Set Signature Score from KEGG or HALLMARK databases ###########
###############################################################################################
source("/Volumes/My_Passport/TC_analysis/Locasale_sc_metabolic_landscape/Single-Cell-Metabolic-Landscape/utils.R")
source("/Volumes/My_Passport/SMO_scData_analysis/GeneSetSignatureScore.R")
kegg_pathways <- gmtPathways("/Volumes/My_Passport/TC_analysis/Locasale_sc_metabolic_landscape/Single-Cell-Metabolic-Landscape/Data/KEGG_metabolism.gmt")
kegg_pathways
View(names(kegg_pathways))
hallmark_pathways <- gmtPathways("/Volumes/My_Passport/TC_analysis/Locasale_sc_metabolic_landscape/Single-Cell-Metabolic-Landscape/Data/h.all.v6.1.symbols.gmt")
hallmark_pathways
View(names(hallmark_pathways))

mat <- GetAssayData(object = pbmc[["RNA"]], slot = "data")
mat[1:5,1:5]
class(mat)
dim(mat)

genes <- hallmark_pathways[['HALLMARK_HEDGEHOG_SIGNALING']]
genes
length(genes)
tmp <- human2mouse(genes)
tmp
input_genes <- tmp$mouseGene
length(which(input_genes %in% rownames(mat)))
genes <- input_genes[which(input_genes %in% rownames(mat))]

df <- read.table("/Volumes/My_Passport/SMO_scData_analysis/SMO_0hr/sigDEs/c3_vs_c0_markerGenes_192.txt",
                 sep = "\t", header = F)
head(df)
input_genes <- df$V1

genes <- input_genes[which(input_genes %in% rownames(mat))]
length(genes)
exprs <- colMeans(mat[genes,],na.rm = T)
#exprs <- colMeans(sub_mat[genes,],na.rm = T)
head(exprs)
all(names(exprs) == colnames(mat))
meta <- pbmc@meta.data
meta$avg_score <- unname(exprs)
head(meta)
scaled_score <- (meta$avg_score - min(meta$avg_score))/(max(meta$avg_score) - min(meta$avg_score))
meta$scaled_avg_score <- scaled_score
head(meta)

range(meta$scaled_avg_score)
range(meta$avg_score)

umap <- as.data.frame(Embeddings(object = pbmc, reduction = "umap"))
head(umap)
all(rownames(umap) == rownames(meta))
umap$avg_score <- meta$avg_score
umap$scaled_avg_score <- meta$scaled_avg_score

ggplot(umap, aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color = scaled_avg_score), size=1.2, alpha=0.5)+
  scale_color_steps(low = "yellow", high = "red") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right") +
  scale_alpha_continuous(range = c(0.5,1.5))



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### slingshot tradeSeq for Heps trajectory #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
"%!in%" <- Negate("%in%")
hep <- subset(pbmc, orig_seurat_clusters %!in% c("L5","L4","L7"))
hep # 14088 features across 2861 samples
DimPlot(hep, reduction = "umap", label = T, 
        group.by = "orig_seurat_clusters")
counts <- hep[["RNA"]]@counts
logcounts <- hep[['RNA']]@data
df <- hep[["umap"]]@cell.embeddings %>% as.data.frame()
head(df)
dim(df) # [1] 2861    2

cl <- data.frame(seurat_clusters = hep@meta.data$orig_seurat_clusters,
                 row.names = rownames(hep@meta.data))
head(cl)

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
sce <- SingleCellExperiment(assays = List(counts = counts, logcounts = logcounts))
df <- as.matrix(df)
head(df)
class(df)
reducedDims(sce) <- list(UMAP=df)
colData(sce)$seurat_clusters <- cl$seurat_clusters

sce
saveRDS(sce, file = "~/Desktop/smo_Luc_48hr_PHX.rds")
saveRDS(sce, file = "~/Desktop/updated_smo_Luc_48hr_PHX.rds")

sce <- readRDS("~/Desktop/smo_Luc_48hr_PHX.rds")
sce <- readRDS("~/Desktop/updated_smo_Luc_48hr_PHX.rds")
sce
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = 'seurat_clusters',
                 start.clus = 'L0')

#counts <- as.matrix(counts(sce))
#counts[1:5,1:5]
crv <- SlingshotDataSet(sce)
crv

#set.seed(5)
#icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, 
#                   nGenes = 200, verbose = T)
#set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
head(pseudotime)

cellWeights <- slingCurveWeights(crv)
head(cellWeights)

pbmc <- readRDS("./Luc_48hr_PHX_2955_PC10.rds")
pbmc
genes <- VariableFeatures(pbmc)
head(genes)
length(genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst")
var_genes <- VariableFeatures(pbmc)
length(var_genes) # 2000

sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
             nknots = 6, verbose = TRUE)
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



