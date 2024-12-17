#Install Seurat v4 using conda

conda create -n seurat r-base=4.0.5
conda install conda-forge::r-v8 
conda install conda-forge::r-seurat #4.1.1
conda install bioconda::bioconductor-singlecellexperiment
conda install bioconda::bioconductor-singler
conda install bioconda::bioconductor-celldex #1.0.0
conda install conda-forge::r-arrow
conda install conda-forge::r-ggpubr
conda install conda-forge::r-tidyverse
conda install conda-forge::r-ggthemes

#Install Seurat v5 using conda

conda create -n seurat5 r-base=4.4.1
conda install conda-forge::r-seurat #
conda install conda-forge::r-arrow
conda install conda-forge::r-ggpubr
conda install conda-forge::r-tidyverse
conda install conda-forge::r-ggthemes
