#!/bin/bash

############################################################################################
############################################################################################

## Assuming BAM files from WGS tumor/normal(*-t/n_SortedRmdup.bam) pairs are in $(pwd) 
## run the configuration scripts to create folders for SV calling

############################################################################################
############################################################################################

PRIMARY="/scratch/tmp/a_toen03/Genomes/hg38/primary_chromosomes.bed.gz"
MANTA="/home/a/a_toen03/software_2c/manta-1.4.0.centos6_x86_64/bin/configManta.py"
HG38="/scratch/tmp/a_toen03/Genomes/hg38/hg38_noALT_withDecoy.fa"

############################################################################################
############################################################################################

ls *-t*.bam | \
  awk -F "-t" '{print $1}' | \
  parallel "$MANTA \
    --normalBam={}-n_SortedRmdup.bam \
    --tumorBam={}-t_SortedRmdup.bam \
    --callRegions=${PRIMARY} \
    --runDir=VariantDir_${BASENAME}/Manta/Somatic/ \
    --referenceFasta=${HG38}"
    
ls *-n_SortedRmdup.bam | \
  awk -F "-n_SortedRmdup.bam" '{print $1}' | \
  parallel "$MANTA \
    --normalBam={}-n_SortedRmdup.bam \
    --callRegions=${PRIMARY} \
    --runDir=VariantDir_${BASENAME}/Manta/Germline/ \
    --referenceFasta=${HG38}"
    
############################################################################################
############################################################################################
