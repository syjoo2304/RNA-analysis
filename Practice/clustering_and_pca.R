#data loading
gc_exp_tb<-read.delim("2024_Clustering_and_PCA_Expression.txt",
                      sep = "\t",
                      header = T,
                      stringsAsFactors = F)

phe_tb<-read.delim("2024_Clustering_and_PCA_Sample_info.txt",
                   sep = "\t",
                   header = T,
                   stringsAsFactors = F)


library(pheatmap)
library(RColorBrewer)

rownames(gc_exp_tb) <- gc_exp_tb[,1] #dim(gc_exp_tb)
gc_exp_tb <- gc_exp_tb[,-1]
gene_sd <- apply(gc_exp_tb, 2, sd)
top <- names(sort(gene_sd,decreasing = T)[1:100])
top_exp_tb <- gc_exp_tb[,top]

our_color<-rev(brewer.pal(n = 5, name = "RdBu"))

pheatmap(top_exp_tb,
         clustering_method = "average",
         clustering_distance_cols = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize_row = 5,
         fontsize_col = 4,
         scale = "column",
         color = our_color)


pca_results<-prcomp(top_exp_tb,
                    center = F,
                    scale = F)

PC_values<-as.data.frame(pca_results$x)
plot(x = PC_values$PC1,
     y = PC_values$PC2,
     xlab = "PC1",
     ylab = "PC2",
     pch = 19,
     main = "PCA of lung cancer cell lines")


pc_variance<-pca_results$sdev
pc_variance_norm<-pc_variance/sum(pc_variance)
names(pc_variance_norm)<-paste0("PC", 1:length(pc_variance_norm))
str(pc_variance_norm)
# Plot explained variance of top 10 principle components
barplot(pc_variance_norm[1:10],
        main = "Explained Variance (Top 10)",
        ylab = "Explained vaiance",
        ylim = c(0,0.25), las = 2) # "las" argument rotates axis labels


pc_rotation<-pca_results$rotation
# order() function gives the order of each value in a vector
pc_rotation<-pc_rotation[order(pc_rotation[,1], decreasing = T),]
pc1_top5gene <- rownames(pc_rotation[1:5,])
top5_exp_tb <- gc_exp_tb[,pc1_top5gene]
pheatmap(top5_exp_tb,
         fontsize_row = 5,
         fontsize_col = 4,
         scale = "column",
         color = our_color)
