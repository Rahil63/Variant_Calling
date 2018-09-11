#!/bin/bash

## Run the postprocessing scripts on all Strelka directories, assuming this script and
## the Strelka2_postprocessing.sh are in the basedir of the variant calling
## in which the Strelka2Dirs are located:

ls Strelka2Dir_RG* | grep ':' | awk -F ":" '{print $1}' | \
  parallel -j 20 \
    "cp Strelka2_postprocessing.sh ./{}/Somatic/results/variants && \
     cd ./{}/Somatic/results/variants && \
     ./Strelka2_postprocessing.sh 2> stderr.log"

## Make a new folder with the final variants as symbolic links:
mkdir Final_Variants

while read p
  do
  ln -s \
    $(realpath ./$p/Somatic/results/variants/somatic.PASSED.ExLC.Ex1KG.vcf.gz) \
    ./Final_Variants/${p#Strelka2Dir_}_somatic.PASSED.ExLC.Ex1KG.vcf.gz
  done < <(ls Strelka2Dir_RG* | grep ':' | awk -F ":" '{print $1}')
