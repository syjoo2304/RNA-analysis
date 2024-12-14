#Load required libraries
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(GO.db)
library(data.table)
library(gplots)

# Load DEGs and expression data
DEGs<-read.table("DEG_ICB_data(response).txt", sep="\t", header=T)
exp_tb<-read.table("Expression_vst_ICB_data.csv",sep=",",header=T)
head(exp_tb)
# Loaexp_tb# Load gene set information
geneset_temp<-readLines("c5.go.bp.v7.4.symbols.gmt")
geneset_temp[[1]]
# Process gene set list
geneset_list <- strsplit(geneset_temp, "\t")  
names(geneset_list) <- sapply(geneset_list, function(x) x[1]) 
head(geneset_list)
geneset_list <- lapply(geneset_list, function(x) x[-c(1:2)])
#Filter DEGs for upregulated and downregulated genes using log2 fold change and p-value criteria
upgenes<-DEGs$X[which(DEGs$log2FoldChange>1 & DEGs$pvalue<0.05)] # Genes with positive log2 fold change
downgenes<-DEGs$X[which(DEGs$log2FoldChange<1& DEGs$pvalue<0.05)] # Genes with negative log2 fold change

# GSA for upregulated genes
background<-setdiff(unique(exp_tb$X),c(upgenes)) #background means whole genes - upgene

# Initialize vectors to store results
in.input <- c()  # DEGs in geneset
out.input <- c() # DEGs not in geneset
in.background <- c() # Genes in geneset but not DEGs (background)
out.background <- c() # Genes not in geneset or DEGs (background)
p.value <- c()
hit.genes <- c()
geneset.size <- c()
fdr <- c()

# Loop through each geneset and perform hypergeometric test
for(i in 1:length(geneset_list)){
  hit.genes[i] <- paste(collapse = ",", intersect(upgenes, geneset_list[[i]])) # Genes that overlap
  in.input[i] <- length(intersect(upgenes, geneset_list[[i]])) # Overlap of DEGs with geneset
  out.input[i] <- length(setdiff(upgenes, geneset_list[[i]])) # DEGs not in geneset
  in.background[i] <- length(intersect(background, geneset_list[[i]])) # Geneset genes not in DEGs
  out.background[i] <- length(setdiff(background, geneset_list[[i]])) # Genes not in DEGs or geneset
  
  # Perform hypergeometric test
  p.value[i] <- phyper(in.input[i]-1, length(upgenes), length(background), in.input[i] + in.background[i], lower.tail = F) 
  
  # Store geneset size
  geneset.size[i] <- length(geneset_list[[i]]) 
}

# Adjust p-values for multiple testing using FDR
fdr <- p.adjust(as.numeric(p.value), "fdr") 

# Store results for upregulated DEGs
up.gsa <- data.frame("Geneset" = names(geneset_list),
                     "Geneset.Size" = geneset.size,
                     "Pvalue" = p.value,
                     "FDR" = fdr,
                     "In" = in.input,
                     "Out" = out.input,
                     "In.Background" = in.background,
                     "Out.Background" = out.background,
                     "Hits" = hit.genes)
head(up.gsa)
# Order results by p-value
up.gsa <- up.gsa[order(up.gsa$Pvalue),]

# Filter Significant Gene Sets
up.gsa_filtered <- up.gsa[up.gsa$In >= 5 &  
                            up.gsa$Geneset.Size >= 10 & up.gsa$Geneset.Size < 200, ]

# Calculate gene ratio and -log10FDR for easier interpretation
up.gsa_filtered$ratio <- up.gsa_filtered$In / up.gsa_filtered$Geneset.Size
up.gsa_filtered$log10FDR <- -log10(up.gsa_filtered$FDR)

# Final table sorted by -log10FDR
up.gsa_final <- up.gsa_filtered[order(up.gsa_filtered$log10FDR, decreasing = TRUE),]


#=============================
# Repeat for downregulated genes
#=============================
background <- setdiff(unique(exp_tb$X), c(downgenes)) 

in.input <- c()
out.input <- c()
in.background <- c()
out.background <- c()
p.value <- c()
hit.genes <- c()
geneset.size <- c()
fdr <- c()

for(i in 1:length(geneset_list)){
  hit.genes[i] <- paste(collapse = ",", intersect(downgenes, geneset_list[[i]])) 
  in.input[i] <- length(intersect(downgenes, geneset_list[[i]]))
  out.input[i] <- length(setdiff(downgenes, geneset_list[[i]]))
  in.background[i] <- length(intersect(background, geneset_list[[i]])) 
  out.background[i] <- length(setdiff(background, geneset_list[[i]])) 
  
  # Perform hypergeometric test for downregulated genes
  p.value[i] <- phyper(in.input[i]-1, length(downgenes), length(background), in.input[i] + in.background[i], lower.tail = F) 
  
  geneset.size[i] <- length(geneset_list[[i]]) 
}

fdr <- p.adjust(as.numeric(p.value), "fdr") 

# Store results for downregulated DEGs
down.gsa <- data.frame("Geneset" = names(geneset_list),
                       "Geneset.Size" = geneset.size,
                       "Pvalue" = p.value,
                       "FDR" = fdr,
                       "In" = in.input,
                       "Out" = out.input,
                       "In.Background" = in.background,
                       "Out.Background" = out.background,
                       "Hits" = hit.genes)
down.gsa <- down.gsa[order(down.gsa$Pvalue),]

# Filter significant gene sets for downregulated DEGs
down.gsa_filtered <- down.gsa[down.gsa$In >= 5 & 
                                down.gsa$Geneset.Size >= 10 & down.gsa$Geneset.Size < 200, ]

# Calculate gene ratio and -log10FDR
down.gsa_filtered$ratio <- down.gsa_filtered$In / down.gsa_filtered$Geneset.Size
down.gsa_filtered$log10FDR <- -log10(down.gsa_filtered$FDR)

# Final table for downregulated DEGs
down.gsa_final <- down.gsa_filtered[order(down.gsa_filtered$log10FDR, decreasing = TRUE),]

# Bubble plot of top 10 gene sets
ggplot(up.gsa_final[1:10,], aes(x = ratio, y = reorder(Geneset, ratio), size = ratio, color = log10FDR)) +
  geom_point() +  # Use points with varying sizes and colors
  scale_color_gradient(low = "pink", high = "red") +  # Color gradient based on -log10FDR
  labs(title = "Gene Set Analysis Results (Upregulated DEGs)",
       y = "Gene Set",
       x = "Gene Ratio",
       color = "-log10FDR", size = "Gene Ratio")   # Axis and legend labels


#===========================================
# Create Bubble Plot for Downregulated DEGs
#===========================================

# Bubble plot of top 10 gene sets
ggplot(down.gsa_final[1:10,], aes(x = ratio, y = reorder(Geneset, ratio), size = ratio, color = log10FDR)) +
  geom_point() +  # Use points with varying sizes and colors
  scale_color_gradient(low = "skyblue", high = "blue") +  # Color gradient based on -log10FDR
  labs(title = "Gene Set Analysis Results (Downregulated DEGs)",
       y = "Gene Set",
       x = "Gene Ratio",
       color = "-log10FDR", size = "Gene Ratio")   # Axis and legend labels


library(pheatmap)
library(RColorBrewer)

test <- fread('Expression_vst_ICB_data2.PROJ.gct')
head(test)
tmp <- data.frame(test[,-2])
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]

our_color<-rev(brewer.pal(n = 20, name = "RdBu"))

pheatmap(tmp,
         clustering_method = "average",
         clustering_distance_cols = "euclidean",
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         fontsize_row = 5,
         fontsize_col = 4,
         color = our_color)
