#!/bin/bash

## Run Canvas on both tumor and normal WGS samples, assuming strelka germline calls have already been produced:

################################################################################################################
################################################################################################################

DOTNET="/home/a/a_toen03/software_2c/dotnet"
CANVAS="/home/a/a_toen03/software_2c/Canvas-1.38.0.1554+master_x64/Canvas.dll"

################################################################################################################

KMER="/scratch/tmp/a_toen03/Genomes/hg38/Canvas/kmer.fa"
FILTER13="/scratch/tmp/a_toen03/Genomes/hg38/Canvas/filter13.bed"
GENOME_FOLDER="/scratch/tmp/a_toen03/Genomes/hg38/Canvas"
BASENAME=$1

################################################################################################################
################################################################################################################

#### Somatic workflow:
echo '[MAIN] Running somatic workflow -- time:' && date
$DOTNET $CANVAS Somatic-WGS \
  --bam=${BASENAME}-t_SortedRmdup.bam \
  --somatic-vcf=Strelka2Dir_${BASENAME}/Somatic/results/variants/somatic.snvs.vcf.gz \
  --sample-b-allele-vcf=Strelka2Dir_${BASENAME}/Germline/results/variants/somatic.snvs.vcf.gz \
  --sample-name=${BASENAME}-t \
  --output=CanvasDir_${BASENAME}/Somatic \
  --reference=${KMER} \
  --genome-folder=${GENOME_FOLDER} \
  --filter-bed=${FILTER13}

################################################################################################################

## Germline workflow:
echo '[MAIN] Running germline workflow -- time:' && date
$DOTNET $CANVAS Germline-WGS \
  --bam=${BASENAME}-n_SortedRmdup.bam \
  --sample-b-allele-vcf=Strelka2Dir_${BASENAME}/Germline/results/variants/somatic.snvs.vcf.gz \
  --sample-name=${BASENAME}-n \
  --output=CanvasDir_${BASENAME}/Germline \
  --reference=${KMER} \
  --genome-folder=${GENOME_FOLDER} \
  --filter-bed=${FILTER13}
  
################################################################################################################

echo '[MAIN] Finished on:' && date 

################################################################################################################
################################################################################################################
