setwd("/Volumes/My_Passport/SMO_scData_analysis/SMO_48hr_PHX")

library(RColorBrewer)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(tidyverse)
library(tradeSeq)
library(slingshot)

seuratfromsce <- readRDS("./savedRDS/Luc+Cre_seuratfromsce.rds")
seuratfromsce # 13212 features across 7258 samples
head(seuratfromsce@meta.data)
table(seuratfromsce@meta.data$seurat_clusters)
"%!in%" <- Negate("%in%")
hep <- subset(seuratfromsce, seurat_clusters %!in% c("L5","C4","C7","C8","C9"))
hep
saveRDS(hep, file = "./savedRDS/Luc+Cre_hepOnly_seuratfromsce.rds")

hep <- readRDS("./savedRDS/Luc+Cre_hepOnly_seuratfromsce.rds")
hep
idx <- which(hep@meta.data$seurat_clusters == "C3")
hep@meta.data$label <- "others"
hep@meta.data$label[idx] <- "C3"
DimPlot(hep, reduction = "UMAP", group.by = "label")
DimPlot(hep, reduction = "UMAP", group.by = "seurat_clusters", label = T)
DimPlot(hep, reduction = "UMAP", group.by = "condition")


counts <- hep[["RNA"]]@counts
counts[1:5,1:5]
logcounts <- hep[['RNA']]@data
logcounts[1:5,1:5]
df <- hep[["UMAP"]]@cell.embeddings %>% as.data.frame()
colnames(df) <- c("UMAP1","UMAP2")
df <- as.matrix(df)
head(df)
dim(df) # [1] 6553    2

cl <- hep@meta.data$seurat_clusters
head(cl)
length(cl) # [1] 6553

sce <- SingleCellExperiment(assays = List(counts = counts, logcounts = logcounts))
#pca <- prcomp(t(assays(sce)$expdata), scale. = FALSE)
#reducedDims(sce) <- list(PCA=pca_data$x, UMAP=df)
reducedDims(sce) <- SimpleList(UMAP=df)
colData(sce)$seurat_clusters <- cl

sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = 'seurat_clusters',
                 start.clus = 'L0')
sce
counts <- as.matrix(counts(sce))
counts[1:5,1:5]
crv <- SlingshotDataSet(sce)
crv

#################################################################
############# RUN GeneSwitch ####################################
#################################################################
library(GeneSwitches)

head(sce@colData)
pseudotime_idx = "slingPseudotime_1"
sce_p1 <- convert_slingshot(sce_slingshot = sce, 
                           pseudotime_idx = "slingPseudotime_1", 
                           assayname = "logcounts")
sce_p1

## Download example files to current directory
get_example_inputData()
## Load input data log-normalized gene expression
load("./logexpdata.RData")
logexpdata[1:5,1:5]
class(logexpdata)
## Load Monocle2 object with pseudo-time and dimensionality reduction
load("./cardiac_monocle2.RData")
cardiac_monocle2

# monocle:::plot_cell_trajectory(cardiac_monocle2, color_by = "State")
plot_monocle_State(cardiac_monocle2)
test <- convert_monocle2(monocle2_obj = cardiac_monocle2, 
                           states = c(3,2,1), expdata = logexpdata)
test
assay(sce_p1)[1:5,1:5]
assay(p1)[1:5,1:5]
head(sce_p1@colData)
head(p1@colData)
### binarize gene expression using gene-specific thresholds
test <- binarize_exp(test, ncores = 2)
range(test$Pseudotime)
range(assay(test))
range(sce_p1$Pseudotime)
range(assay(sce_p1))
assays(sce_p1)$expdata[1:5,1:5]
hist(assays(test)$expdata, breaks = 20, plot = TRUE)
h
plot(h, freq = FALSE, main = "Histogram of gene expression",
     xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
abline(v=0.2, col="blue")


sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = 0.2, ncores = 2)
sce_p1
assays(sce_p1)$binary[1:5,1:5]
## fit logistic regression and find the switching pseudo-time point for each gene
## with downsampling. This step takes less than 1 mins
sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = TRUE, show_warning = FALSE)
sce_p1
# DataFrame with 15 rows and 13 columns
## filter top 15 best fitting switching genes among all the genes
sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, topnum = 15)
sg_allgenes # DataFrame with 15 rows and 13 columns
colnames(sg_allgenes)
# [1] "geneID"            "zerop_gene"        "passBinary"        "switch_at_time"   
# [5] "CI"                "pvalues"           "pseudoR2s"         "estimates"        
# [9] "prd_quality"       "direction"         "FDR"               "switch_at_timeidx"
# [13] "feature_type" 
## filter top 15 best fitting switching genes among surface proteins and TFs only
sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, topnum = 20,
                                genelists = gs_genelists, genetype = c("Surface proteins", "TFs"))
sg_gtypes
## combine switching genes and remove duplicated genes from sg_allgenes
sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])
sg_vis # DataFrame with 15 rows and 13 columns

plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3)


