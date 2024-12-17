#!/bin/bash

sqantdir='/home/syjoo/Project/long-read/program/SQANTI3'
refdir_='/home/syjoo/Project/long-read/REF'
workdir_='/home/syjoo/Project/long-read/working_dir'
sample_gtf='clustered_2_mc.gff'
prefix='clustered_2_mc'
python $sqantdir/sqanti3_qc.py $workdir_/$sample_gtf $refdir_/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf \
$refdir_/GRCm38.p6.genome.fa -o ${prefix} -d $workdir_/SQANTI3_QC_output_${prefix}  --cpus 4 --report html



