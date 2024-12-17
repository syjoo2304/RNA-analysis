#!/usr/bin/env Rscript 
library(scran)
library(celldex)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
##Scoreheatmap
#conda install -c conda-forge r-viridis
#conda install -c bioconda r-pheatmap
library(viridis)
library(pheatmap)

##V1 (rough version)
seurat_obj <-  readRDS('[path_to_rds]/[seurat].rds')
test.clusters <- seurat_obj$seurat_clusters
test <- as.SingleCellExperiment(seurat_obj)
bpe <- BlueprintEncodeData()
pred.bpe2 <- SingleR(test = test, ref = bpe, genes = "de",de.method = "wilcox",
                     labels = bpe$label.main)

#select label
pred.bpe <- SingleR(test = test, ref = bpe, genes = "de",de.method = "wilcox",
                     labels = bpe$label.main, clusters = test.clusters)
bpe.clusters <- test.clusters
levels(bpe.clusters) <- pred.bpe$labels
seurat_obj[["SingleR.labels"]] <- bpe.clusters
Idents(seurat_obj) <- seurat_obj[["SingleR.labels"]]
select.labels <- unique(seurat_obj$SingleR.labels)

#drawing heatmap
png('test1.score_hm.png')
plotScoreHeatmap(pred.bpe2,labels.use=select.labels, clusters=seurat_obj$orig.ident, show.labels=F)
dev.off()
========================================================================
##subcluster
library(Seurat)
test <- readRDS('[path_to_rds]/[singler_object].singler.rds')
#install.packages('SeuratObject')
#library(SeuratObject)
========================================================================
test_small <- subset(test, idents="Epithelial cells")
test_small <- FindVariableFeatures(test_small, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(test_small)
test_small <- ScaleData(test_small, features = all.genes)
test_small <- RunPCA(test_small, features = VariableFeatures(object = test_small))
test_small <- JackStraw(test_small, num.replicate = 100)
test_small <- ScoreJackStraw(test_small, dims = 1:20)
JackStrawPlot(test_small, dims = 1:20)
ElbowPlot(test_small)
test_small <- FindNeighbors(test_small, dims = 1:6, k.param = 5)
test_small <- FindClusters(test_small)
head(Idents(test_small), 5)

test_small <- RunUMAP(test_small, dims = 1:6)
DimPlot(test_small, reduction = "umap")

test_small <- RunTSNE(test_small, dims = 1:6)
DimPlot(test_small, reduction = "tsne")

test_small.markers <- FindAllMarkers(test_small, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
test_small.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
epithelial.marker <- test_small.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(epithelial.marker, 'epithelial_marker.txt', sep='\t',quote=F)
write.table(epithelial.marker,"[path_to_output]/epithelial.txt",col.names=T, quote=F)

DimPlot(test_small)

#panglao database
