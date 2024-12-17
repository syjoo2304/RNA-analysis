#!/bin/bash

sqantdir='[path_to_sqanti3]'
isoannotdir='[path_to_isoannotatlite]'
refdir_='[path_to_reference]'
workdir_='[working_dir]'
sample_gtf='clustered_2_mc.gff'
prefix='clustered_2_mc'
python $sqantdir/sqanti3_qc.py $workdir_/$sample_gtf $refdir_/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf \
$refdir_/GRCm38.p6.genome.fa -o ${prefix} -d $workdir_/SQANTI3_QC_output_${prefix}  --cpus 4 --report html

python3 $isoannotdir/IsoAnnotLite_SQ3.py ${prefix}.gtf ${prefix}_classification.txt ${prefix}_junctions.txt  \
-gff3 ${refdir_}/Mus_musculus_GRCm38_Ensembl_86.gff3 -o ${prefix}_newGFF3 -stdout ${prefix}_statisticalResults -novel -nointronic -saveTranscriptIDs


