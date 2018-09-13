#!/bin/bash

############################################################################################
############################################################################################

## Assuming BAM files from WGS tumor/normal(*-t/n_SortedRmdup.bam) pairs are in $(pwd) 
## run the configuration scripts to create folders for SV calling

############################################################################################
############################################################################################

## Take variables from this file:
source Define_Paths.sh

############################################################################################
############################################################################################

ls *-t*.bam | \
  awk -F "-t" '{print $1}' | \
  parallel "$MANTA \
    --normalBam={}-n_SortedRmdup.bam \
    --tumorBam={}-t_SortedRmdup.bam \
    --callRegions=${PRIMARY} \
    --runDir=VariantDir_{}/Manta/Somatic/ \
    --referenceFasta=${HG38}"
    
ls *-n_SortedRmdup.bam | \
  awk -F "-n_SortedRmdup.bam" '{print $1}' | \
  parallel "$MANTA \
    --normalBam={}-n_SortedRmdup.bam \
    --callRegions=${PRIMARY} \
    --runDir=VariantDir_{}/Manta/Germline/ \
    --referenceFasta=${HG38}"
    
############################################################################################
############################################################################################
