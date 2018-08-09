#!/bin/bash

## Read a VCF file (hg38), and use dbSNP151 VCF to remove variants that have an AF
## of greater 1% in 1KG (=flag "COMMON") or > 1% in TOPMED. Assumes uncompressed VCF:

#################################################################################################

## Usage: ./VCF_Filter_Common.sh in.vcf NUMBER_OF_FORKS > filtered_vcf.vcf

IN_VCF=$1
NFORK=$2

if [[ ! -e All_20180418.vcf.gz ]]; then
  echo '[INFO]: Downloading dbSNP151 VCF file and its index'
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms//human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz .
  wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms//human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi .
  fi
  
## Annotate VCF:
## COMMON (=1 if MAF 1KG > 1%), TOPMED (=AF of REF and ALT in TOPMED, exact match)
## Use filter_vep to exclude these variants from the VCF (stdout)
vep --fork $NFORK--vcf --format vcf --custom All_20180418_chr.vcf.gz,dbSNP151,vcf,exact,,COMMON,TOPMED \
  -i ${IN_VCF} -o STDOUT | \
  filter_vep --filter "not dbSNP151_COMMON = 1 and not dbSNP151_TOPMED > 0.01"
  

