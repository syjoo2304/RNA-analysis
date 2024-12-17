#### Install 'nextflow' using mamba
mamba install nextflow 
conda activate base # I installed nextfow in my base environment.

#### Download igenome reference of your interest from Illumina. Here is an example command line for downloading Human genome GRCh37.  
axel -a -n 20 https://s3.amazonaws.com/igenomes.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz

#### Running nf-core/circRNA workflow 
nextflow run [path_to_gitcloned_folder]/circrna/ --input [path_to_sample]/samples.csv --outdir [path_to_outputdir]/[output_dirname] --genome GRCm38 --igenomes_base [path_to_igenome_reference] --phenotype  [path_to_phenotype]/phenotypes.csv --module circrna_discovery,mirna_detection -profile docker --mature [path_to_smallrna_references] --hisat2 [path_to_hisat2_ht2l_references]
