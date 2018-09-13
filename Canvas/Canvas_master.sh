#!/bin/bash

## Run Canvas on both tumor and normal WGS samples, assuming strelka germline calls have already been produced.
## Note, requires to be run via Make-Canvas_Scripts.sh

################################################################################################################
################################################################################################################

DOTNET="/home/a/a_toen03/software_2c/dotnet"
CANVAS="/home/a/a_toen03/software_2c/Canvas-1.38.0.1554+master_x64/Canvas.dll"
#BASENAME=see Make-Scripts_Canvas.sh
################################################################################################################

KMER="/scratch/tmp/a_toen03/Genomes/hg38/Canvas/kmer.fa"
FILTER13="/scratch/tmp/a_toen03/Genomes/hg38/Canvas/filter13.bed"
GENOME_FOLDER="/scratch/tmp/a_toen03/Genomes/hg38/Canvas"

################################################################################################################
################################################################################################################

#### Somatic workflow:
echo '[MAIN] Running somatic workflow -- time:' && date
$DOTNET $CANVAS Somatic-WGS \
  --bam=${BASENAME}-t_SortedRmdup.bam \
  --somatic-vcf=VariantDir_${BASENAME}/Srelka2/Somatic/result/results/variants/somatic.snvs.vcf.gz \
  --sample-b-allele-vcf=VariantDir_${BASENAME}/Strelka2/Germline/results/variants/variants.vcf.gz \
  --sample-name=${BASENAME}-t \
  --output=VariantDir_${BASENAME}/Canvas/Somatic \
  --reference=${KMER} \
  --genome-folder=${GENOME_FOLDER} \
  --filter-bed=${FILTER13}

################################################################################################################

## Germline workflow:
echo '[MAIN] Running germline workflow -- time:' && date
$DOTNET $CANVAS Germline-WGS \
  --bam=${BASENAME}-n_SortedRmdup.bam \
  --sample-b-allele-vcf=VariantDir_${BASENAME}/Strelka2/Germline/results/variants/variants.vcf.gz \
  --sample-name=${BASENAME}-n \
  --output=VariantDir_${BASENAME}/Canvas/Germline \
  --reference=${KMER} \
  --genome-folder=${GENOME_FOLDER} \
  --filter-bed=${FILTER13}
  
################################################################################################################

echo '[MAIN] Finished on:' && date 

################################################################################################################
################################################################################################################
