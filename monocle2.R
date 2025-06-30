setwd("/Volumes/My_Passport/SMO_scData_analysis/")

library(Seurat)
library(monocle)

luc_48 <- readRDS("./SMO_48hr_PHX/savedRDS/Luc_48hr_PHX_2955_PC10.rds")
luc_48
luc_48@meta.data$seurat_clusters <- paste0("L",luc_48@meta.data$seurat_clusters)
head(luc_48@meta.data)


cre_48 <- readRDS("./SMO_48hr_PHX/savedRDS/Cre_48hr_PHX_4357_PC10.rds")
cre_48
head(cre_48@meta.data)

combined <- merge(luc_48, y = cre_48, add.cell.ids = c("Luc","Cre"), project = "48hr_PHX")
combined # 7312
head(combined@meta.data)
tail(combined@meta.data)
"%!in%" <- Negate("%in%")
combined <- subset(combined, seurat_clusters %!in% c('L5','C4','C7','C8','C9'))
combined # 6606

luc_48 <- subset(luc_48, seurat_clusters %!in% c('L5'))
luc_48 <- subset(luc_48, seurat_clusters %in% c('L3','L4','L7','L2'))
luc_48 # 2861
raw_count_data <- GetAssayData(luc_48, assay = "RNA", slot = "counts")
class(raw_count_data)
raw_count_data[1:5,1:5]

cells_info <- luc_48@meta.data
genes_info <- data.frame(row.names = rownames(raw_count_data),
                         gene_id = rownames(raw_count_data),
                         gene_short_name = rownames(raw_count_data))
head(genes_info)
head(cells_info)
pd <- new("AnnotatedDataFrame", data = cells_info)
fd <- new("AnnotatedDataFrame", data = genes_info)
cds <- newCellDataSet(raw_count_data, phenoData = pd, featureData = fd)
cds

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
table(disp_table$mean_expression>=0.1) # FALSE 8655  TRUE 6075
#clustering_genes <- subset(disp_table, mean_expression >= 0.1)
hv_genes <- VariableFeatures(luc_48)
head(hv_genes)
length(hv_genes)  # 2000
cds <- setOrderingFilter(cds, hv_genes)
cds <- reduceDimension(cds, num_dim = 20, reduction_method = 'tSNE')
cds <- clusterCells(cds, method = "louvain")
plot_cell_clusters(cds, cell_size = 0.5) +
  theme(legend.position = "none") +
  labs(x = "tSNE1", y = "tSNE2")
plot_cell_clusters(cds, cell_size = 0.5, color_by = "seurat_clusters") +
  #scale_color_brewer(name = "cell type", type = "qual", palette = "Set2") +
  labs(x = "tSNE1", y = "tSNE2") +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3)))

#cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- reduceDimension(cds, norm_method = 'log',verbose = T, max_components = 10)
head(pData(cds))
#table(pData(cds)$Cluster)
#table(pData(cds)$State)
#table(pData(cds)$State, pData(cds)$seurat_clusters)[,'L2']
pData(cds)$State <- pData(cds)$seurat_clusters
cds <- orderCells(cds, root_state = c("L0","L1","L2"))
plot_cell_trajectory(cds, color_by = "seurat_clusters", cell_size = 1,
                     show_branch_points = F)
plot_complex_cell_trajectory(cds, color_by = "seurat_clusters",
                             show_branch_points = T, cell_size = 0.5,cell_link_size = 0.3,
                             root_states = 'L2')

cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points



