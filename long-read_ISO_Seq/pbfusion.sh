#!/bin/bash

##step01
pbfusion gff-cache \
    --gtf [path_to_ref]/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf \
    --gtf-out [path_to_ref]/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.bin

##step02-option
pbfusion discover --gtf [path_to_ref]/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.bin \
       --output-prefix isoseq \
       --threads 8 \
       --min-fusion-quality LOW \
       --min-coverage 1 \
       --min-fusion-read-fraction 0.01 \
       --log-level "info" \
       mapped.bam

#step00
pbmm2 align --preset ISOSEQ --sort 4_Raw264_flnc.bam [path_to_ref]/GRCm38.p6.genome.mmi 4_Raw264_mapped.bam

#step02-default
pbfusion discover --gtf [path_to_ref]/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.bin \
        --output-prefix mapped_ \
        --threads 8 \
       --log-level "info" \
       mapped.bam
      

#step03-visualization
python3 visualize_fusion2.py \
    -o fusion_browser_shot2.png \
    -a [path_to_ref]/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf \
    -f mapped_.breakpoints.groups.bed \
    -b mapped.bam
