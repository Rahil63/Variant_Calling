#!/bin/bash

############################################################################################
############################################################################################

## Assuming BAM files from WGS tumor/normal(*-t/n_SortedRmdup.bam) pairs are in $(pwd) 
## run the configuration scripts to create folders for variant calling
## See: https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md

############################################################################################
############################################################################################

## Take variables from this file:
source Define_Paths.sh

############################################################################################
############################################################################################

## Assumes that Manta has already been run to get the candidate indels.

## => Somatic workflow:
ls *-t*.bam | \
  awk -F "-t" '{print $1}' | \
    parallel "$STRELKA2/configureStrelkaSomaticWorkflow.py \
      --normalBam={}-n_SortedRmdup.bam \
      --tumorBam={}-t_SortedRmdup.bam \
      --callRegions=${PRIMARY} \
      --runDir=VariantDir_{}/Strelka2/Somatic \
      --referenceFasta=${HG38} \
      --indelCandidates=VariantDir_{}/Manta/Somatic/results/variants/candidateSmallIndels.vcf.gz"
   
## => Germline workflow:
ls *-n* | \
  awk -F "-n" '{print $1}' | \
    parallel "$STRELKA2/configureStrelkaGermlineWorkflow.py \
      --bam={}-n_SortedRmdup.bam \
      --callRegions=${PRIMARY} \
      --runDir=VariantDir_{}/Strelka2/Germline/ \
      --referenceFasta=${HG38} \
      --indelCandidates=VariantDir_{}/Manta/Germline/results/variants/candidateSmallIndels.vcf.gz"
                        
############################################################################################
############################################################################################
