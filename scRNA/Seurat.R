#scRNA data of cochlea from 1,2,5,12,15 months of aged mice( Data source: Aging atlas (https://ngdc.cncb.ac.cn/aging/index),https://doi.org/10.1093/procel/pwac058 )
#Load Pre-requsite packages
library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(sctransform)

# Load seurat object (integrated and manually cell-type-annotated)
sh.integrated <- readRDS('integrated_clustered_annotated_SCTcorrected_seurat.rds')

# Define the categories for pre-defined clusters
Corti <- c('HC', 'DC_PC', 'IPhC_IBC', 'TBC', 'Nudt4+_PC', 'EC')
Modiolus <- c('SGN', 'SGC', 'SC', 'CC', 'OB')
Reissner <- c('RMC')
Stria <- c('IC', 'MC', 'BC', 'CEC', 'PVM_M')
Spiral <- c('FB', 'FC1', 'FC2', 'FC3', 'FC4', 'SMC')
Immune <- c('M', 'T', 'B', 'Neu')

# Retrieve the existing cell types
cell_type1 <- as.vector(sh.integrated@meta.data$cell_type)

# Create a new vector to store the updated categories
cell_type2 <- cell_type1

# Update the categories based on conditions
for (i in 1:length(cell_type1)) {
  if (cell_type1[i] %in% Corti) {
    cell_type2[i] <- "Organ_of_corti"
  } else if (cell_type1[i] %in% Modiolus) {
    cell_type2[i] <- "Modiolus"
  } else if (cell_type1[i] %in% Reissner) {
    cell_type2[i] <- "Reissners_membrane"
  } else if (cell_type1[i] %in% Stria) {
    cell_type2[i] <- "Stria_vascularis"
  } else if (cell_type1[i] %in% Spiral) {
    cell_type2[i] <- "Spiral_ligament"
  } else if (cell_type1[i] %in% Immune) {
    cell_type2[i] <- "Immune_cell"
  }
}

# Assign the updated categories back to the metadata
sh.integrated@meta.data$cell_type2 <- cell_type2

# Visalize Gene of Interest (In this case, six deafness genes) using DotPlot
Idents(sh.integrated) <- "cell_type2"
pdf("Dotplot_bymonth_vcelltype.pdf",width=10,height=10)
DotPlot(sh.integrated, features = GOI, group.by = "cell_type2",split.by="age",cols=c(rep("blue",6)))+ RotatedAxis()+ coord_flip()
dev.off()

f1 <- DotPlot(sh.integrated, features = GOI, group.by="cell_type",split.by="age",idents="Organ_of_corti",cols=c(rep("blue",6)))+ RotatedAxis() + coord_flip()+
  annotate("text", x = 6.5, y = length(GOI)+10, label = "Organ of Corti", size = 5, fontface = "bold")

f2 <- DotPlot(sh.integrated, features = GOI, group.by="cell_type",split.by="age",idents="Modiolus",cols=c(rep("blue",6)))+ RotatedAxis() + coord_flip()+
  annotate("text", x = 6.5, y = length(GOI)+10, label = "Modiolus", size = 5, fontface = "bold")

f3 <- DotPlot(sh.integrated, features = GOI, group.by="cell_type",split.by="age",idents="Reissners_membrane",cols=c(rep("blue",6)))+ RotatedAxis() + coord_flip()+
  annotate("text", x = 6.5, y = length(GOI)-2.5, label = "Reissner's membrane", size = 5, fontface = "bold")

f4 <- DotPlot(sh.integrated, features = GOI, group.by="cell_type",split.by="age",idents="Stria_vascularis",cols=c(rep("blue",6)))+ RotatedAxis() + coord_flip()+
  annotate("text", x = 6.5, y = length(GOI)+10, label = "Stria vascularis", size = 5, fontface = "bold")

f5 <- DotPlot(sh.integrated, features = GOI, group.by="cell_type",split.by="age",idents="Spiral_ligament",cols=c(rep("blue",6)))+ RotatedAxis() + coord_flip()+
  annotate("text", x = 6.5, y = length(GOI)+10, label = "Spiral ligament", size = 5, fontface = "bold")

f6 <- DotPlot(sh.integrated, features = GOI, group.by="cell_type",split.by="age",idents="Immune_cell",cols=c(rep("blue",6)))+ RotatedAxis() + coord_flip()+
  annotate("text", x = 6.5, y = length(GOI)+10, label = "Immune cell", size = 5, fontface = "bold")

pdf("Dotplot_bymonth_vcelltype2.pdf",width=10,height=30)
wrap_plots(plots = list(f1,f2,f3,f4,f5,f6), nrow=6, ncol=1)
dev.off()
