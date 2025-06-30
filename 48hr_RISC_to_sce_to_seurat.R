setwd("/Users/tianyichen/Desktop/Regeneration_manuscript")

library(RColorBrewer)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(tidyverse)


#counts <- readMM("RISC_0+1+2+3+5_smo_0hr_mat.mtx") # raw counts
#counts[1:5,1:5]
logcounts <- readMM("logmat.mtx") # logcounts/RISC corrected logmat
logcounts[1:5,1:5]
genes <- read.table("gene.tsv", sep = "\t",
                    header = T, row.names = 1)
head(genes)
rownames(counts) <- genes$Symbol
rownames(logcounts) <- genes$Symbol
cells <- read.table("cell.tsv", sep = "\t",
                    header = T, row.names = 1, check.names = F)
head(cells)
tail(cells)
colnames(counts) <- cells$scBarcode
colnames(logcounts) <- cells$scBarcode

logcounts[1:5,1:5]
df <- as.data.frame(logcounts)
df[1:5,1:5]
write.table(df, file = "./expression.tsv", sep = "\t", row.names = T, 
            col.names = T, quote = F)
df <- read.table(file = "./expression.tsv", sep = "\t", row.names = 1,
                 header = T)
df[1:5,1:5]


umap <- read.table("umap.tsv", sep = "\t",
                   header = F, check.names = F)
head(umap)
rownames(umap) <- colnames(counts)
umap <- as.matrix(umap)
counts <- as.matrix(counts)
counts[1:5,1:5]
logcounts <- as.matrix(logcounts)
logcounts[1:5,1:5]
dim(logcounts) # [1] 13212  7258

#cl <- paste0("c",cells$seurat_clusters)
cl <- cells$seurat_clusters
length(cl) # 7258
head(cl)
class(cl)

sce <- SingleCellExperiment(assays = List(counts = logcounts, logcounts = logcounts))
reducedDims(sce) <- list(UMAP=umap)
colData(sce)$seurat_clusters <- cl
sce # dim: 13566 5293

seuratfromsce <- as.Seurat(sce, counts = "counts", data = "logcounts")
seuratfromsce
head(seuratfromsce@meta.data)
seuratfromsce@assays$RNA[1:5,1:5]
#seuratfromsce@assays$data <- logcounts
#seuratfromsce@assays$data[1:5,1:5]
Idents(seuratfromsce) <- "seurat_clusters"
DimPlot(seuratfromsce, reduction = "UMAP", group.by = "seurat_clusters", label = T) 
seuratfromsce$condition <- cells$condition
DimPlot(seuratfromsce, reduction = "UMAP", group.by = "condition") 

gene = c('Afp','Sox9')
VlnPlot(seuratfromsce, features = gene)
VlnPlot(seuratfromsce, features = gene, group.by = "condition")
DotPlot(seuratfromsce, features = gene, group.by = "condition")
FeaturePlot(seuratfromsce, features = gene)
FeatureScatter(seuratfromsce, feature1 = "Smo",
                              feature2 = "Epcam")

saveRDS(seuratfromsce, file = "../Luc+Cre_seuratfromsce.rds")

seuratfromsce <- readRDS("./Luc+Cre_seuratfromsce.rds")
seuratfromsce <- readRDS("./Luc+Cre_hepOnly_seuratfromsce.rds")
seuratfromsce
head(seuratfromsce@meta.data)
table(seuratfromsce@meta.data$seurat_clusters)
DimPlot(seuratfromsce, reduction = "UMAP", group.by = "condition")
DimPlot(seuratfromsce, reduction = "UMAP", group.by = "seurat_clusters", label = T)
table(seuratfromsce@meta.data$seurat_clusters)

### Find Markers ####################
hep <- SetIdent(seuratfromsce, value = "condition")
hep <- SetIdent(seuratfromsce, value = "seurat_clusters")

hep # 13212 features across 7258
hep[['RNA']][1:5,1:5]
head(hep@meta.data)

cluster.markers <- FindMarkers(hep, ident.1 = c("C5"), ident.2 = c('L0','L1','L2','L3','L4','L6','L7'), 
                               only.pos = TRUE, logfc.threshold = 0.25) # log(2)
cluster.markers <- FindMarkers(hep, ident.1 = c('Cre_48hr_PHX'), ident.2 = c('Luc_48hr_PHX'), 
                               only.pos = TRUE, logfc.threshold = 0.25)
head(cluster.markers)
cluster.markers <- cluster.markers %>% filter(p_val_adj <= 0.05)
dim(cluster.markers) #[1] 1459    5
tail(cluster.markers)
cluster.markers <- cluster.markers[order(cluster.markers$avg_log2FC, decreasing = T),]
markers <- rownames(cluster.markers)
gs <- rownames(hep)
length(gs) # [1] 13212
library(GeneOverlap)
tmp <- read_csv("./aging_geneSig_l2fc.csv")
head(tmp)
tmp <- as.data.frame(tmp)
tmp_markers <- tmp$Gene_Name
length(tmp_markers) # [1] 497

go.obj <- newGeneOverlap(markers, tmp_markers, genome.size = length(gs))
go.obj
go.obj <- testGeneOverlap(go.obj)
go.obj
print(go.obj)
getIntersection(go.obj)
save_genes <- as.data.frame(getIntersection(go.obj))
head(save_genes)
write.table(save_genes, file = "./geneOverlap_aging30m_smoNeg.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)


table(hep@meta.data$condition)
cluster.markers <- FindMarkers(hep, ident.1 = c('Cre_48hr_PHX'), ident.2 = c('Luc_48hr_PHX'), 
                               min.cells.group = 1, 
                               min.cells.feature = 1,
                               min.pct = 0,
                               logfc.threshold = 0, #log(2),
                               only.pos = FALSE)
cluster.markers <- FindMarkers(hep, ident.1 = c('C5'), ident.2 = c('L0','L1','L2','L3','L4','L6','L7'),
                               min.cells.group = 1, 
                               min.cells.feature = 1,
                               min.pct = 0,
                               logfc.threshold = 0, #log(2),
                               only.pos = FALSE)
head(cluster.markers, n = 5)
tail(cluster.markers)
dim(cluster.markers) # [1] 12913    5
# order list, pull out gene name and log2fc, and convert genes to uppercase
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
ranked_list <- na.omit(ranked_list)  
head(ranked_list)
tail(ranked_list)
ranked_list$metric[1] <- 290
ranked_list$metric[2] <- 290

write.table(ranked_list, file = "./gsea/hep_48hr_PHX_Cre_vs_Luc_AllHeps_gsea.rnk", sep = "\t", row.names = F, quote = F, col.names = F)

## Extract sig DEs for HOMER analysis
DE <- cluster.markers %>% filter(p_val_adj <= 0.05) # avg_log2FC > 0 & 
head(DE)
tail(DE)
dim(DE) # [1] 298   5
DE <- DE[order(DE$avg_log2FC, decreasing = T),]
dim(DE)
posDE <- DE %>% filter(avg_log2FC < 0)
dim(posDE)
df <- data.frame(genes = rownames(posDE))
head(df)
dim(df)
write.table(df, file = "/Users/tianyichen/Desktop/SMO_48hr_PHX/L3_vs_C3_negDE_padj005.txt",
            sep = "\t", row.names = F, quote = F, col.names = F)

### Using EXTEND to predict telomerase activity #####
seuratfromsce # 13212 features across 6553 samples
mat <- seuratfromsce@assays$RNA[1:13212,1:6553]
mat[1:5,1:5]
class(mat)
#write.table(mat, file = "./Luc+Cre_hepOnly_seuratfromsce_data.txt", sep = "\t", row.names = T, 
#            col.names = T, quote = F)
mat[1:5,1:5]
mat <- as.matrix(mat)
class(mat)
source("/Volumes/My_Passport/SMO_scData_analysis/EXTEND/R/ComponentAndMarkerFunction.r")
source("/Volumes/My_Passport/SMO_scData_analysis/EXTEND/R/ComponentOneAndMarkerFunction.r")
source("/Volumes/My_Passport/SMO_scData_analysis/EXTEND/R/ComponentTwoAndMarkerFunction.r")
source("/Volumes/My_Passport/SMO_scData_analysis/EXTEND/R/InputData.r")
source("/Volumes/My_Passport/SMO_scData_analysis/EXTEND/R/IterativeRS.r")
source("/Volumes/My_Passport/SMO_scData_analysis/EXTEND/R/MarkerFunction.r")
source("/Volumes/My_Passport/SMO_scData_analysis/EXTEND/R/RunEXTEND.r")

genes = c("4930503L19Rik","Cep72","Pole2","Mcm4","Hells","Rexo5",
          "Ehmt2","Lin9","Tert","Terc","Paxip1")
genes[which(genes %in% rownames(mat))]
"Terc" %in% rownames(mat)

setwd("/Volumes/My_Passport/SMO_scData_analysis/SMO_48hr_PHX/")
RunEXTEND(mat)
tel_scores <- read.table("TelomeraseScores.txt", header = T, row.names = 1, check.names = F)
head(tel_scores)
rownames(tel_scores) <- sapply(strsplit(as.character(rownames(tel_scores)), "\\."), "[",1)
head(tel_scores)
rownames(tel_scores) <- paste0(rownames(tel_scores),"-1")
head(tel_scores)
meta <- seuratfromsce@meta.data
head(meta)
idx <- intersect(rownames(meta),rownames(tel_scores))
length(idx) # 6553
meta <- meta[idx,]
tel_scores <- tel_scores[idx,]
all(rownames(meta) == rownames(tel_scores)) # [1] TRUE
meta$scores <- tel_scores$NormEXTENDScores
#meta$seurat_clusters <- paste0("c",meta$seurat_clusters)
head(meta)

ggboxplot(meta,x="seurat_clusters",y="scores",size =  0.8, width = 0.6,outlier.shape=NA,color = "seurat_clusters", add = "jitter",add.params= list(size=0.3),legend="none")+
  theme(axis.text.x=element_text(angle=30,size=8,vjust=1,hjust=1))+
  labs(y = "EXTEND Scores") 

my_comparisons1 = list(c("C2","C0"),
                       c("C2","C1"),c("C2","C6"),c("C1","C6"),
                       c("C1","C3"),c("C1","C5"),c("C6","C3"),
                       c("C6","C5"),c("C3","C5"))
# c("L0","L1"),c("L2","L1"),c("L1","L3"),
#c("L3","L6"),c("L1","L6"),c("L1","L4"),
#c("L1","L7"),c('L4','L7'),
ggboxplot(meta,x="seurat_clusters",y="scores",size =  0.8, width = 0.6,outlier.shape=NA,color = "seurat_clusters", add = "jitter",add.params= list(size=0.3),legend="none")+
  theme(axis.text.x=element_text(angle=30,size=8,vjust=1,hjust=1))+
  stat_compare_means(comparisons = my_comparisons1, method= "t.test",size=2.3) +
  labs(y = "EXTEND Scores") 

sub_meta <- meta %>% filter(seurat_clusters %in% c("L3","L4","L7","C3","C5"))
head(sub_meta)
ggboxplot(sub_meta,x="seurat_clusters",y="scores",size =  0.8, width = 0.6,outlier.shape=NA,color = "seurat_clusters", add = "jitter",add.params= list(size=0.3),legend="none")+
  theme(axis.text.x=element_text(angle=30,size=8,vjust=1,hjust=1))+
  labs(y = "EXTEND Scores") 

my_comparisons2 = list(c("L3","L4"),c("L3","L7"),c("L3","C3"),
                       c("L3","C5"),c("L4","L7"),c("L4","C3"),
                       c("L4","C5"),c("L7","C3"),c("L7","C5"),
                       c("C3","C5"))
ggboxplot(sub_meta,x="seurat_clusters",y="scores",size =  0.8, width = 0.6,outlier.shape=NA,color = "seurat_clusters", add = "jitter",add.params= list(size=0.3),legend="none")+
  theme(axis.text.x=element_text(angle=30,size=8,vjust=1,hjust=1))+
  stat_compare_means(comparisons = my_comparisons2, method= "t.test",size=2.3) +
  labs(y = "EXTEND Scores") 

my_comparisons3 = list(c("Luc_48hr_PHX","Cre_48hr_PHX"))
ggboxplot(sub_meta,x="condition",y="scores",size =  0.8, width = 0.6,outlier.shape=NA,color = "condition", add = "jitter",add.params= list(size=0.3),legend="none")+
  theme(axis.text.x=element_text(angle=30,size=8,vjust=1,hjust=1))+
  stat_compare_means(comparisons = my_comparisons3, method= "t.test",size=2.3) +
  labs(y = "EXTEND Scores") 

