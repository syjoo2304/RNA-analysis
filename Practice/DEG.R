#====================================================================================#
# Lecture
# Differentially expressed genes (DEGs) analysis based on negative binomial distribution
# We will query and download data from The Cancer Genome Project (TCGA) and perform DEG analysis
#====================================================================================#
#=====================================================================================================================#
#=========================#
# How is count data made?
#=========================#
# Simple example script that to be run in a Linux-like environment:
# Suppose you have sequenced sample A by paired-end sequencing
# You have two files: A_1.fq, A_2.fq
# You must supply reference genome sequences (FASTA) and genome annotations (GTF)
# These files can be downloaded from ENSEMBL or UCSC

# First, you must build a genome index file with your reference files
# The genome indexes are save to the disk and only need to be generated one time for each genome/annotation combination
# STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /path/to/genomeDir --sjdbGTFfile /path/to/annotations.gtf

# # Perform alignment with 'STAR' (github.com/alexdobin/STAR)
# # Below is an example with basic options
# # Refer the STAR manual (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for more information
# STAR --runThreadN 10 --outSAMtype BAM --twopassMode Basic --readFilesIn A_1.fq A_2.fq \
#   --sjdbGTFfile /path/to/GTF_file --genomeDir /path/to/genomeDir --quantMode GeneCounts

# The output file ~.tab will contain your count data!
#=====================================================================================================================#

#====================#
# Install required packages
#====================#
# Colaprico A, et al. "TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data." Nucleic acids research (2015): gkv1507.
# Mounir et al. ?€œNew functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx.?€?
# PLoS computational biology, 15(3), e1006701.
# Install TCGAbiolinks for querying TCGA data from GDC data portal
# https://portal.gdc.cancer.gov/
# https://github.com/BioinformaticsFMRP/TCGAbiolinks
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("remotes")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

# Martin Morgan, Valerie Obenchain, Jim Hester and HervÃ© PagÃ¨s (2021). SummarizedExperiment: SummarizedExperiment container. R package version
# 1.22.0. https://bioconductor.org/packages/SummarizedExperiment
# Install SummarizedExperiment. We will use this package to process downloaded TCGA data
BiocManager::install("SummarizedExperiment",force=TRUE)

# Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
# Genome Biol 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
# Install DESeq2
BiocManager::install("DESeq2")

#====================#
# Query, download and prepare TCGA-GBM data
#
# The Cancer Genome Atlas (TCGA), a landmark cancer genomics program,
# molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types.
# CGA generated over 2.5 petabytes of genomic, epigenomic, transcriptomic, and proteomic data.
# The data, which has already led to improvements in our ability to diagnose, treat, and prevent cancer,
# will remain publicly available for anyone in the research community to use.
# All data collected and processed by the program is currently available at the Genomic Data Commons (GDC).
# https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga
#====================#
# Set working directory
setwd("C:\\Users\\Sunyoung\\OneDrive\\¹ÙÅÁ È­¸é\\2024_RNA_analysis_½Ç½À") # Change according to your working directory

# Import required packages
library(TCGAbiolinks)
library(SummarizedExperiment)

# TCGA-KIRP
# Query
query <- GDCquery(
  project = "TCGA-KIRP",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

query
str(query)

# Download
GDCdownload(query = query)
# Prepare downloaded data
tcga_obj <- GDCprepare(query = query)

# # Load TCGA-ACC data from RDS file
tcga_obj <- readRDS("./tcga_kirp_data.rds")

# Get row data (mRNA data)
tcga_rowdata <- as.data.frame(rowData(tcga_obj))

# Check mRNA data
head(tcga_rowdata)
unique(tcga_rowdata$gene_type)

# Get clinical data
clinical_data <- as.data.frame(colData(tcga_obj))

# Check clinical data
head(clinical_data)
str(clinical_data)
colnames(clinical_data)

# Get expression data
exp_data <- assay(tcga_obj)

# Check expression data
exp_data[1:10, 1:10]
str(exp_data)
# Remove expression data of non-protein coding genes
protein_coding_genes <- tcga_rowdata$gene_id[tcga_rowdata$gene_type == "protein_coding"]

pc_exp_data <- exp_data[protein_coding_genes, ]
str(pc_exp_data)
# There are 3 samples types within TCGA-GBM: TP (primary solid tumor), TR (recurrent solid tumor), TAP (Additional - New Primary)
unique(clinical_data$shortLetterCode)
unique(clinical_data$definition)

# Get Primary solid tumor data and normal solid tissue data
tp_samps <- rownames(clinical_data)[clinical_data$shortLetterCode == "TP"]
tp_exp_data <- pc_exp_data[, tp_samps]

nt_samps <- rownames(clinical_data)[clinical_data$shortLetterCode == "NT"]
nt_exp_data <- pc_exp_data[, nt_samps]

###########
# QUIZ!
# Get sample names of recurrent solid tumor samples in TCGA-GBM
###########

# Combine
dim(tp_exp_data)
dim(nt_exp_data)
tp_nt_exp_data <- cbind(tp_exp_data, nt_exp_data)
dim(tp_nt_exp_data)
# Create sample data for DESeq and other operations
samp_data <- data.frame(
  row.names = colnames(tp_nt_exp_data),
  "Type" = sapply(colnames(tp_nt_exp_data), function(x) {
    x <- ifelse(is.element(x, tp_samps), "Primary tumor", "Normal")
  }),
  stringsAsFactors = T
)

#====================#
# Perform DEG analysis with DESeq2
#
# The rapid adoption of high-throughput sequencing (HTS) technologies for genomic studies has resulted in a need
# for statistical methods to assess quantitative differences between experiments. An important task here is
# the analysis of RNA sequencing (RNA-seq) data with the aim of finding genes that are differentially expressed across groups of samples.
# DESeq2, a method for differential analysis of count data, using shrinkage estimation for dispersions and fold changes
# to improve stability and interpretability of estimates.
# This enables a more quantitative analysis focused on the strength rather than the mere presence of differential expression.
#====================#
# Import DESeq2
library(DESeq2)

# Create DESeq object
samp_data$Type
dds <- DESeqDataSetFromMatrix(
  countData = tp_nt_exp_data,
  colData = samp_data,
  design = ~Type
)

# Pre-filter data
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions,
# there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads,
# we reduce the memory size of the dds data object, and we increase the speed of the transformation and testing functions within DESeq2.
# Remove data of genes with less than 667 total reads (mean count of 1/sample)
dim(dds) # 19962 322
dim(dds[rowSums(counts(dds)) > 667, ]) 

filt_dds <- dds[rowSums(counts(dds)) > 667, ]

# Calculate results, deseq?•¨?ˆ˜ë¥? ?‹œ?–‰?•˜?Š”?‹¨ê³?
filt_dds <- DESeq(filt_dds)
deseq_res <- results(
  filt_dds,
  contrast = c("Type", "Primary tumor", "Normal")
)

# View results
deseq_res
summary(deseq_res)

# Sort results by significance (FDR)
deseq_res <- deseq_res[order(deseq_res$padj), ]

# Create a data frame version of deseq_res for easier viewing
deseq_res_df <- as.data.frame(deseq_res)
View(deseq_res_df)

# MA-plot
# MA plots display a log ratio (M) vs an average (A) in order to visualize the differences between two groups.
plotMA(deseq_res)

# Use plotCounts function to compare normalized counts between tissue types
# Top 25 differentially expressed genes (by FDR)
par(mfrow = c(3, 3))
for(i in 1:9) {
  plotCounts(filt_dds, gene = rownames(deseq_res)[i], intgroup = "Type", pch = 19)
}

# Volcano plot
par(mfrow = c(1, 1))
with(deseq_res, plot(log2FoldChange, -log10(padj), pch = 20, main = "Volcano plot"))

# Add colored points: grey if padj (adjusted p-value, FDR) > 0.00001,
# red if fold change > 2^4 and padj < 0.00001, blue if fold change < -(2^4) and padj < 0.00001
with(subset(deseq_res, padj > .00001 ), points(log2FoldChange, -log10(padj), pch = 20, col = "grey"))
with(subset(deseq_res, padj < .00001 & log2FoldChange > 4), points(log2FoldChange, -log10(padj), pch = 20, col = "red"))
with(subset(deseq_res, padj < .00001 & log2FoldChange < -4), points(log2FoldChange, -log10(padj), pch = 20, col = "blue"))

###########
# QUIZ
# Color the dots on the volcano plot with your colors of choice
###########

# Get differentially expressed gene matrix
# p-adjusted<0.00001, |log2FC|>4
dds_significant <- deseq_res[!is.na(deseq_res$padj) &
                               deseq_res$padj < 0.00001 &
                               abs(deseq_res$log2FoldChange) > 4, ]
head(dds_significant)
nrow(dds_significant) # Number of DEG

# Sort in log2(fold change) order
dds_significant <- dds_significant[order(dds_significant$log2FoldChange, decreasing = T), ]

# Add gene names to results
dds_significant$gene_name <- sapply(rownames(dds_significant), function(x) {
  tcga_rowdata$gene_name[tcga_rowdata$gene_id == x]
})

# View results
View(as.data.frame(dds_significant))

# Export results as tab delimited text file
write.table(
  dds_significant,
  file = "./genome_informatics_deg_results.txt",
  sep = "\t", quote = F
)

