# Follow up with seurat_48hr_Luc to perform CV cluster correlation 
# with Cre_48hr scRNA-seq data.

setwd("/Volumes/My_Passport/SMO_scData_analysis/SMO_48hr_PHX/savedRDS")
library(Seurat)

pbmc <- readRDS("./Luc_48hr_PHX_2955_PC10.rds")
pbmc
head(pbmc@meta.data)
table(pbmc@meta.data$seurat_clusters)
#  0   1   2   3   4   5   6   7
norm_counts <- pbmc[['RNA']]@data
norm_counts[1:5,1:5]
#idx <- rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "0"),])
l0 <- norm_counts[,rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "0"),])]
l0[1:5,1:5]
dim(l0)
l1 <- norm_counts[,rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "1"),])]
l2 <- norm_counts[,rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "2"),])]
l3 <- norm_counts[,rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "3"),])]
l4 <- norm_counts[,rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "4"),])]
l5 <- norm_counts[,rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "5"),])]
l6 <- norm_counts[,rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "6"),])]
l7 <- norm_counts[,rownames(pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == "7"),])]

mean_per_gene <- apply(l0, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(l0, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene,
                          sd = sd_per_gene,
                          cv = sd_per_gene/mean_per_gene)
head(cv_per_gene)

df <- data.frame(l0 = cv_per_gene$cv)

mean_per_gene <- apply(l1, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(l1, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df$l1 <- cv_per_gene$cv

mean_per_gene <- apply(l2, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(l2, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df$l2 <- cv_per_gene$cv

mean_per_gene <- apply(l3, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(l3, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df$l3 <- cv_per_gene$cv

mean_per_gene <- apply(l4, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(l4, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df$l4 <- cv_per_gene$cv

mean_per_gene <- apply(l5, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(l5, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df$l5 <- cv_per_gene$cv

mean_per_gene <- apply(l6, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(l6, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df$l6 <- cv_per_gene$cv

mean_per_gene <- apply(l7, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(l7, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df$l7 <- cv_per_gene$cv

head(df)
dim(df)
dim(norm_counts)
rownames(df) <- rownames(norm_counts)
df <- df[complete.cases(df),]
head(df)
dim(df)


##  read in Cre_48hr scRNA-seq data
## Seurat lognormalized count data and associated meta data
cre <- readRDS("./Cre_48hr_PHX_4357_PC10.rds")
cre
head(cre@meta.data)
table(cre@meta.data$seurat_clusters)
#  0   1   2   3   4   5   6   7  8 9
norm_counts <- cre[['RNA']]@data
norm_counts[1:5,1:5]
ct <- norm_counts
ct[1:5,1:5]
meta <- cre@meta.data
meta$label <- meta$seurat_clusters
head(meta)
c0 <- ct[,rownames(meta[which(meta$label == "0"),])]
c1 <- ct[,rownames(meta[which(meta$label == "1"),])]
c2 <- ct[,rownames(meta[which(meta$label == "2"),])]
c3 <- ct[,rownames(meta[which(meta$label == "3"),])]
c4 <- ct[,rownames(meta[which(meta$label == "4"),])]
c5 <- ct[,rownames(meta[which(meta$label == "5"),])]
c6 <- ct[,rownames(meta[which(meta$label == "6"),])]
c7 <- ct[,rownames(meta[which(meta$label == "7"),])]
c8 <- ct[,rownames(meta[which(meta$label == "8"),])]
c9 <- ct[,rownames(meta[which(meta$label == "9"),])]
dim(c1)
dim(c7)
table(meta$label)

mean_per_gene <- apply(c1, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c1, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2 <- data.frame(c1 = cv_per_gene$cv)
head(df_2)

mean_per_gene <- apply(c2, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c2, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c2 <- cv_per_gene$cv

mean_per_gene <- apply(c3, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c3, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c3 <- cv_per_gene$cv

mean_per_gene <- apply(c4, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c4, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c4 <- cv_per_gene$cv

mean_per_gene <- apply(c5, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c5, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c5 <- cv_per_gene$cv

mean_per_gene <- apply(c6, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c6, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c6<- cv_per_gene$cv

mean_per_gene <- apply(c7, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c7, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c7 <- cv_per_gene$cv

mean_per_gene <- apply(c8, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c8, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c8 <- cv_per_gene$cv

mean_per_gene <- apply(c9, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c9, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c9 <- cv_per_gene$cv

mean_per_gene <- apply(c0, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(c0, 1, sd, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene, sd = sd_per_gene, cv = sd_per_gene/mean_per_gene)
df_2$c0 <- cv_per_gene$cv

head(df_2)
rownames(df_2) <- rownames(ct)
sum(complete.cases(df_2))
df_2 <- df_2[complete.cases(df_2),]
head(df_2)

df$gene <- rownames(df)
idx <- intersect(rownames(df), rownames(df_2))
length(idx) # [1] 8434
df <- df[idx,]
df_2 <- df_2[idx,]
head(df)
head(df_2)

df <- cbind(df,df_2)
#df <- df[,-9]
head(df)
dim(df) # [1] 8434   18
write.table(df, file = "../smo_Luc_plus_Cre_48hr_PHX_coefficient_of_variation.txt",
            row.names = T, col.names = T, quote = F, sep = "\t")

# calculate correlation
res <- cor(df, method = "pearson", use = "complete.obs")
round(res, 2)
dim(res) # [1] 18 18
View(res)
res <- res[1:8,9:18]
res <- t(res)
#rownames(res) <- paste0("r",seq(0,8))
#colnames(res) <- paste0("L",seq(0,7))

write.table(res, file = "./correlation_table.txt",row.names = T, col.names = T,
            quote = F)

cl.topics <- fastcluster::hclust.vector(t(res), method="ward", metric="euclidean")
dd.topics <- as.dendrogram(cl.topics)

cn = colnames(res)
library(ComplexHeatmap)
ComplexHeatmap::Heatmap(as.matrix(res), cluster_columns=dd.topics, name = "Correlation(r)",
                        show_column_names=FALSE, cluster_rows=TRUE, show_row_names = TRUE,
                        bottom_annotation = HeatmapAnnotation(
                          text = anno_text(cn, rot = 45, offset = unit(1, "npc"), just = "right"),
                          annotation_height = max_text_width(cn)
                        ))

ComplexHeatmap::Heatmap(as.matrix(res), cluster_columns=FALSE, name = "Correlation(r)",
                        show_column_names=TRUE, cluster_rows=FALSE, show_row_names = TRUE)

library("ggpubr")
head(tmp)
colnames(tmp) <- c(colnames(res),rownames(res))
head(tmp)
ggscatter(tmp, x = "L7", y = "r6", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "SMO_L7", ylab = "Previous_r6_proliferative")





