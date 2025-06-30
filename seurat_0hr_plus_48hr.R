setwd("/Volumes/My_Passport/SMO_scData_analysis/")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Hmisc)
library(EnhancedVolcano)
library(magrittr)
library(tidyverse)
#library(SeuratDisk)

# Load the dataset
cre_0hr.data <- Read10X(data.dir = "./Cre_0hr/filtered_feature_bc_matrix")
cre_0hr <- CreateSeuratObject(counts = cre_0hr.data, min.cells = 3, min.features = 200)
cre_0hr
cre_0hr[["percent.mt"]] <- PercentageFeatureSet(cre_0hr, pattern = "^mt-")
head(cre_0hr@meta.data, 5)
VlnPlot(cre_0hr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cre_0hr <- subset(cre_0hr, subset = percent.mt <= 60)
cre_0hr
cre_0hr <- subset(cre_0hr, subset = nFeature_RNA > 1600 & percent.mt <= 60)
cre_0hr
summary(cre_0hr@meta.data)

cre_0hr <- NormalizeData(cre_0hr, normalization.method = "LogNormalize", scale.factor = 10000)
raw_counts <- as.data.frame(cre_0hr[['RNA']]@counts)
normalized_counts <- as.data.frame(cre_0hr[['RNA']]@data)
raw_counts[1:5,1:5]

write.table(raw_counts, file = "./Cre_0hr_PHX_umi.txt", row.names = T, col.names = T,
            quote = F, sep = "\t")
write.table(normalized_counts, file = "./Cre_0hr_logNorm_umi.txt", row.names = T, col.names = T,
            quote = F, sep = "\t")



luc_0hr.data <- Read10X(data.dir = "./Luc_0hr/filtered_feature_bc_matrix")
luc_0hr <- CreateSeuratObject(counts = luc_0hr.data, min.cells = 3, min.features = 200)
luc_0hr
luc_0hr[["percent.mt"]] <- PercentageFeatureSet(luc_0hr, pattern = "^mt-")
head(luc_0hr@meta.data, 5)
VlnPlot(luc_0hr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
luc_0hr <- subset(luc_0hr, subset = nFeature_RNA > 1600 & percent.mt <= 60)
luc_0hr
summary(luc_0hr@meta.data)

luc_48hr.data <- Read10X(data.dir = "./Luc_48hr_PHX/filtered_feature_bc_matrix")
luc_48hr <- CreateSeuratObject(counts = luc_48hr.data, min.cells = 3, min.features = 200)
luc_48hr
luc_48hr[["percent.mt"]] <- PercentageFeatureSet(luc_48hr, pattern = "^mt-")
head(luc_48hr@meta.data, 5)
VlnPlot(luc_48hr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
luc_48hr <- subset(luc_48hr, subset = nFeature_RNA > 1600 & percent.mt <= 60 & nFeature_RNA < 7000)
luc_48hr
summary(luc_48hr@meta.data)

cre_48hr.data <- Read10X(data.dir = "./Cre_48hr_PHX/filtered_feature_bc_matrix")
cre_48hr <- CreateSeuratObject(counts = cre_48hr.data, min.cells = 3, min.features = 200)
cre_48hr
cre_48hr[["percent.mt"]] <- PercentageFeatureSet(cre_48hr, pattern = "^mt-")
head(cre_48hr@meta.data, 5)
VlnPlot(cre_48hr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cre_48hr <- subset(cre_48hr, subset = nFeature_RNA > 1600 & percent.mt <= 60)
cre_48hr
summary(cre_48hr@meta.data)
###################
pbmc <- merge(cre_0hr, y = c(luc_0hr,cre_48hr,luc_48hr), add.cell.ids = c("Cre_0hr", "Luc_0hr", "Cre_48hr_PHX",
                                                                         "Luc_48hr_PHX"), project = "HEP_SMO")
pbmc # 17399 features across 17862 samples
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.6)

pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap", label = TRUE)
saveRDS(pbmc, file = "./Cre&Luc_0hr_48hr_17862.rds")

#######################################################################
#######################################################################
#######################################################################
pbmc <- readRDS("./Cre&Luc_0hr_48hr_17862.rds")
pbmc # 17399 features across 17862 samples
head(pbmc@meta.data)
pbmc@meta.data$condition <- paste0(sapply(strsplit(as.character(rownames(pbmc@meta.data)), "_"), "[",1),
                                   sapply(strsplit(as.character(rownames(pbmc@meta.data)), "_"), "[",2))
head(pbmc@meta.data)
DimPlot(pbmc, reduction = "umap", group.by = "condition")
DimPlot(pbmc, reduction = "umap", split.by = "condition", label = TRUE)
features=c('Alb','Tat','Cyp7a1','Apoa1','Ttr','Serpina1c','Hnf4a','Fga','Cebpa','Fah','Rbp4','Krt8','Krt18')
features=c('Lrat','Des') # ,'Acta2','Gfap','Msn'
features=c('Cd52','Ptprc')
features=c('Alb','Tat','Cyp7a1','Apoa1')
features=c('Axin2','Tbx3','Yap1','Sox9','Afp','Lgr5','Igf2bp3','Tnfrsf12a') # 
features=c('Mki67','Top2a','Ccnb1','Foxm1','Ccna2','Cenpm','Cdk1') # ,'Ccnd1'
DotPlot(hep, features = features, group.by = "updated_clusters") + RotatedAxis()
VlnPlot(pbmc, features = features)
FeaturePlot(pbmc, features = features, label = TRUE)
FeaturePlot(pbmc, features = c("Gli2","Smo"))
FeaturePlot(pbmc, features = "Gli2",split.by = "condition")

#### Subsetting to hepatocytes population only
hep <- subset(pbmc, subset = seurat_clusters %in% c('0','1','2','3','4','6','8','10','12','13','14','15') )
hep # 17399 features across 14842 samples 
head(hep@meta.data)
dim(hep@meta.data) # [1] 6911    7
table(hep@meta.data$seurat_clusters)
DimPlot(hep, reduction = "umap", label = TRUE)
DimPlot(hep, reduction = "umap", split.by = "condition", label = TRUE)
hep@reductions$umap %>% class()

saveRDS(hep, file = "./heps_Cre&Luc_0hr_plus_48hrPHX_14842.rds")

saveRDS(hep, file = "./heps_Cre&Luc_0hr_plus_48hrPHX_14839.rds")

#########################################################
#########################################################
#########################################################
hep <- readRDS("./heps_Cre&Luc_0hr_plus_48hrPHX_14842.rds")
hep
head(hep@meta.data)
hep_0hr <- readRDS("./heps_Cre&Luc_0hr_6911.rds")
hep_0hr
head(hep_0hr@meta.data)
hep_0hr@meta.data$new_cluster_label <- paste0("orig.",hep_0hr@meta.data$seurat_clusters)

hep@meta.data$updated_clusters <- hep@meta.data$seurat_clusters
all(rownames(hep_0hr@meta.data) %in% rownames(hep@meta.data))
length(which(rownames(hep_0hr@meta.data) %in% rownames(hep@meta.data))) # [1] 6879 out of 6911; 32 cells not included
hep@meta.data[1,"updated_clusters"]
hep_0hr@meta.data[1,"new_cluster_label"]
hep@meta.data$updated_clusters <- as.character(hep@meta.data$updated_clusters)

for(i in rownames(hep_0hr@meta.data)) {
  if(i %in% rownames(hep@meta.data)) {
  hep@meta.data[i,"updated_clusters"] <- hep_0hr@meta.data[i,"new_cluster_label"]
  } else {
    print(i)
  }
}

head(hep@meta.data)
table(hep@meta.data$updated_clusters)
'%!in%' = Negate('%in%')
hep <- subset(hep, subset = updated_clusters %!in% c('0') )

idx <- which(rownames(hep_0hr@meta.data) %!in% rownames(hep@meta.data))
cell_names <- rownames(hep_0hr@meta.data)[idx]
cell_names
hep_0hr@meta.data$missing_cells <- "0"
hep_0hr@meta.data$missing_cells[idx] <- "1"
table(hep_0hr@meta.data$missing_cells)
DimPlot(hep_0hr, reduction = "umap", label = FALSE, group.by = "missing_cells",
        label.size = 5, cols = c("grey87","red1")) + NoLegend()

idx_1 <- which(rownames(pbmc@meta.data) %in% cell_names)
idx_1
pbmc@meta.data$missing_cells <- "0"
pbmc@meta.data$missing_cells[idx_1] <- "1"
test <- subset(pbmc, subset = missing_cells %in% c('1'))
test
table(test@meta.data$seurat_clusters)

DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "missing_cells",
        label.size = 5, cols = c("grey87","red1")) + NoLegend()



DimPlot(hep, reduction = "umap", label = TRUE, group.by = "updated_clusters",
        label.size = 5) + NoLegend()
DimPlot(hep, reduction = "umap", split.by = "condition", label = TRUE,
        group.by = "updated_clusters")

genes=c("Vegfa")
VlnPlot(hep, features = genes)
VlnPlot(hep, features = genes, split.by = "condition")
FeaturePlot(hep, features = c("Ptch1"), split.by = "condition")
# Find markers
# find markers for every cluster compared to all remaining cells, report only the positive ones
hep.markers <- FindAllMarkers(hep, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.58)
# min.diff.pct = 0.25, max.cells.per.ident = 200
hep.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% view()
top10 <- hep.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p <- DoHeatmap(hep, features = top10$gene) + NoLegend()
ggsave("test.jpg", p,width=20, height=10)

markers <- FindMarkers(hep, ident.1 = 3, ident.2 = 0, logfc.threshold = 0.25, test.use = "MAST", only.pos = TRUE)
head(markers, n=20)

DEenrichRPlot(
  hep,
  ident.1 = 3,
  ident.2 = 0,
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = NULL,
  max.genes,
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = ,
  num.pathway = 10,
  return.gene.list = FALSE)

###########################

cluster.markers <- FindMarkers(hep, ident.1 = 2, ident.2 = c(3),
                               min.cells.group = 1, 
                               min.cells.feature = 1,
                               min.pct = 0,
                               logfc.threshold = 0,
                               only.pos = FALSE)
dim(cluster.markers) # [1] 13864  5
class(cluster.markers)
head(cluster.markers, n = 5)

EnhancedVolcano(cluster.markers,
                lab = rownames(cluster.markers),
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Cluster 2 vs.3',
                pCutoff = 10e-16,
                FCcutoff = 0.58,
                pointSize = 0.8,
                labSize = 3.0,
                col=c('grey', 'black', 'black', 'red3'),
                colAlpha = 0.8)
dev.off()
c1vs3 <- cluster.markers
c2vs1 <- cluster.markers
c2vs3 <- cluster.markers


head(c1vs3)
head(c2vs1)
write.table(c1vs3, file = "./DE_1_vs_3.txt", sep = "\t", quote = F,
            row.names = T, col.names = T)
write.table(c2vs1, file = "./DE_2_vs_1.txt", sep = "\t", quote = F,
            row.names = T, col.names = T)
write.table(c2vs3, file = "./DE_2_vs_3.txt", sep = "\t", quote = F,
            row.names = T, col.names = T)


BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library("AnnotationDbi")
library("org.Mm.eg.db")
keytypes(org.Mm.eg.db)

# we want the log2 fold change
df <- c2vs1
head(df)
df$ensid = mapIds(org.Mm.eg.db,
                         keys=rownames(df), 
                         column="ENSEMBL",
                         keytype="SYMBOL",
                         multiVals="first")
head(df)
dim(df) # [1] 13856     6
library(tidyr)
df <- df %>% drop_na()
dim(df) # [1] 12799     6
df <- df[!duplicated(df$ensid), ]
dim(df) # [1] 12796     6

original_gene_list <- df$avg_log2FC
names(original_gene_list) <- df$ensid
gene_list<-na.omit(original_gene_list) # # omit any NA values
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")




tmp <- c3vs1plus2 %>% filter(p_val < 10e-16) %>% filter(avg_log2FC <= -0.58) %>% rownames() %>% as.data.frame()
head(tmp)
write.table(tmp, file = "./sigDEgenes_c3vs1plus2_neg.txt", sep = "\t", quote = F, row.names = F,col.names = F)

features <- c3vs1plus2 %>% filter(p_val < 10e-16) %>% filter(abs(avg_log2FC) >= 0.58) %>% rownames()

hep@meta.data$seurat_clusters <- paste0('c',hep@meta.data$seurat_clusters)
hep <- subset(hep, subset = seurat_clusters %in% c('c0','c1','c2','c3'))
hep
head(hep@meta.data)

levels(hep)
levels(hep) <- c("3", "0", "2", "1")
p <- DoHeatmap(subset(hep, downsample = 100), features = features, size = 3)
ggsave("heatmap.jpg", p,width=13, height=10)


dim(cluster.markers) # [1] 13212     5
write.table(cluster.markers, file = "./DE_0plus3_vs_1plus2.txt", sep = "\t", quote = F,
            row.names = T, col.names = T)
View(cluster.markers)
VlnPlot(hep, features = c("Fgl1","Gm42031"))

library(cowplot)
theme_set(theme_cowplot())
t.cells <- subset(hep, subset = seurat_clusters %in% c('0','3'))
t.cells
head(t.cells@meta.data)
Idents(t.cells) <- "seurat_clusters"
avg.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
head(avg.cells)
#avg.t.cells$gene <- rownames(avg.t.cells)
colnames(avg.cells) <- c("cluster0","cluster3")
head(avg.cells)
### establish linear regression model
model <- lm(cluster3~cluster0, data = avg.cells)
model
avg.cells$resid <-  data.frame(resid=resid(model)) 
head(avg.cells)
test <- avg.cells[order(-avg.cells$resid),]
head(test)

#head(markers) %>% rownames()
genes.to.label = rownames(test)[1:12]

p1 <- ggplot(avg.cells, aes(x=cluster0, y=cluster3)) + geom_point() + ggtitle("Cluster 3 vs. 0")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
ggsave("test.jpg", p1,width=13, height=9)




## percentage histogram plot
df <- data.frame(clusters=pbmc@meta.data$seurat_clusters,
                 samples=pbmc@meta.data$condition)
df <- data.frame(clusters=hep@meta.data$updated_clusters,
                 samples=hep@meta.data$condition)

head(df)
probs = data.frame(prop.table(table(df),1))
probs
library(tidyr)
probs <- probs %>% drop_na()
probs <- probs %>% dplyr::group_by(clusters) %>%
  dplyr::mutate(pos=Freq-0.08) %>%
  dplyr::mutate(pos=ifelse(Freq==0,NA,pos))

brks <- c(0, 0.25, 0.5, 0.75, 1)
ggplot(data=probs,aes(x=clusters,y=Freq,fill=samples)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))


probs$pos[13]=0.45
ggplot(data=probs,aes(x=clusters,y=Freq,fill=samples)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  geom_text(data=probs, aes(x = clusters, y = pos,
                            label = paste0(round(100*Freq),"%")), size=3)

write.table(probs, file = "./hep_stacked_histplot.txt",sep = "\t",
            quote = F, row.names = F, col.names = T)

## plot pathway activity scores
tmp <- read.table("./reactome_hedgehog_on_state_geneset.txt",header = T)
dim(tmp) # 1] 36  1
head(tmp)
colnames(tmp) <- "name"
tmp$name <- tolower(tmp$name) 
tmp$name <- capitalize(tmp$name)
head(tmp)
length(which(tmp$name %in% rownames(hep@assays$RNA)))
# [1] 31
idx <- which(tmp$name %in% rownames(hep@assays$RNA))
idx
tmp <- tmp[idx,]
head(tmp)

head(hep@meta.data)
hep <- AddModuleScore(object = hep, features = list(tmp),
                          name = 'REACTOME_hh_ON')
head(hep@meta.data)
ggplot(hep@meta.data, aes(x=updated_clusters, y=REACTOME_hh_ON1)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="yellow")+
  labs(title="REACTOME_hedgehog_ON_state",x="Heps cluster", y = "Scores") +
  theme_classic()

## Create interactive 3D plot
hep <- RunUMAP(hep,dims = 1:20,n.components = 3L)
# Extract UMAP information from Seurat Object
umap_1 <- hep[["umap"]]@cell.embeddings[,1]
umap_2 <- hep[["umap"]]@cell.embeddings[,2]
umap_3 <- hep[["umap"]]@cell.embeddings[,3]
# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = hep, reduction = "umap")
# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = hep, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))
head(plot.data)
# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))
# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~UMAP_2, y = ~UMAP_1, z=~UMAP_3, 
               color = ~seurat_clusters,
               colors = c("pink",
                          "brown",
                          "darkorchid1",
                          "magenta",
                          "lawngreen",
                          "darkgreen",
                          "lightgreen",
                          "lightblue2",
                          "royalblue1"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 3, width=2), # controls size of points
               text=~seurat_clusters, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
fig



#pdf("./Vlnplot_hep_markers.pdf",width=3, height=3)
source("./utility.R")
p <- StackedVlnPlot(obj = pbmc, features = features)
ggsave("test.jpg", p,width=6, height=10)
#dev.off()

########################################################################################################
## Convert Seurat object to AnnData
SaveH5Seurat(pbmc, filename = "Luc0_Cre0.h5Seurat")
Convert("Luc0_Cre0.h5Seurat", dest = "h5ad")

pbmc_ad <- Convert(from = "Luc0_Cre0.h5Seurat", to = "anndata", filename = "./test.h5ad")
pbmc_ad
########################################################################################################
#install.packages('plotly')
library(plotly)




