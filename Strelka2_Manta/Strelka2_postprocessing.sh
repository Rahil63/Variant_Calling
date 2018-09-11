#!/bin/bash

## Script assumes to be in the basedir of the variant calling, then loops
## through the strelka dirs:
## -- adding AF to the SNVs and Indels, 
## -- merges the two files,
## -- excludes variants in low-complexity,
## -- annotates remaining stuff with VEP and 
## -- removes COMMON (1KG) variants.

LC="/scratch/tmp/a_toen03/Genomes/hg38/hg38_sdust_sorted.bed.gz"
DBSNP151="/scratch/tmp/a_toen03/Databases/dbSNP/All_20180418_chr.vcf.gz"

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

echo '[INFO] Started on:' && date && echo ''

## Add somatic variant allele frequency, see:
## https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#somatic
## for the reference manual on this:

echo '[MAIN] Getting somatic variant allele frequencies'
## Add allele frequency according to Illumina recommendation:

#- Extract allelic information with bcftools for SNPs:
bcftools view -f 'PASS' somatic.snvs.vcf.gz | \
  bcftools query -f '%REF\t%ALT\t[%AU{0}\t%CU{0}\t%GU{0}\t%TU{0}\t]\n' | cut -f1-10 > snp.allelicInfo.txt
  
#- ALT AF in Normal:
paste <(mawk '{if ($1 == "A") print $3; else if ($1 == "C") print $4; else if ($1 == "G") print $5; else if ($1 == "T") print $6}' snp.allelicInfo.txt) \
      <(mawk '{if ($2 == "A") print $3; else if ($2 == "C") print $4; else if ($2 == "G") print $5; else if ($2 == "T") print $6}' snp.allelicInfo.txt) | \
        mawk '{ printf("%3.2f\n", 100*$2 / ($2+$1)) }' > snp.AFnormal.txt

#- ALT AF in Tumor:
paste <(mawk '{if ($1 == "A") print $7; else if ($1 == "C") print $8; else if ($1 == "G") print $9; else if ($1 == "T") print $10}' snp.allelicInfo.txt) \
      <(mawk '{if ($2 == "A") print $7; else if ($2 == "C") print $8; else if ($2 == "G") print $9; else if ($2 == "T") print $10}' snp.allelicInfo.txt) | \
        mawk '{ printf("%3.2f\n", 100*$2 / ($2+$1)) }' > snp.AFtumor.txt

rm snp.allelicInfo.txt
############################################################################################################################################################

## AF for Indels:
#- and for Indels Tumor:
bcftools view -f 'PASS' somatic.indels.vcf.gz | \
  bcftools query -f '[%TAR{0}\t%TIR{0}\t]\n' | \
    mawk '{printf("%3.2f\n", 100*$4 / ($3+$4))}' > indel.AFtumor.txt

#- and for Indels Normal:
bcftools view -f 'PASS' somatic.indels.vcf.gz | \
  bcftools query -f '[%TAR{0}\t%TIR{0}\t]\n' | \
    mawk '{printf("%3.2f\n", 100*$2 / ($1+$2))}' > indel.AFnormal.txt

############################################################################################################################################################


echo '[MAIN] Writing allele frequencies into the VCF files'

##- add SNP AF to VCF:
cat <(bcftools view -h somatic.snvs.vcf.gz | grep '^##') \
    <(echo '##FORMAT=<ID=ALTAF,Number=1,Type=Float,Description="ALT allele frequency [%] added by custom scripts following the Illumina recommendation at the Strelka2 manual page">') \
    <(bcftools view -h somatic.snvs.vcf.gz | grep '^#CHROM') \
    <(paste <(bcftools view -H -f 'PASS' somatic.snvs.vcf.gz | mawk 'OFS = "\t" {$9=$9":ALTAF"; print $0}' | cut -f1-9) \
            <(bcftools view -H -f 'PASS' somatic.snvs.vcf.gz | cut -f10 | paste -d ":" /dev/stdin snp.AFnormal.txt) \
            <(bcftools view -H -f 'PASS' somatic.snvs.vcf.gz | cut -f11 | paste -d ":" /dev/stdin snp.AFtumor.txt)) | bgzip > somatic.snvs.PASSED.vcf.gz

##- add Indel AF to VCF:
cat <(bcftools view -h somatic.indels.vcf.gz | grep '^##' | grep -vE 'bcftools|chrU|_random') \
    <(echo '##FORMAT=<ID=ALTAF,Number=1,Type=Float,Description="ALT allele frequency [%] added by custom scripts following the Illumina recommendation at the Strelka2 manual page">') \
    <(bcftools view -h somatic.indels.vcf.gz | grep '^#CHROM') \
    <(paste <(bcftools view -H -f 'PASS' somatic.indels.vcf.gz | mawk 'OFS = "\t" {$9=$9":ALTAF"; print $0}' | cut -f1-9) \
            <(bcftools view -H -f 'PASS' somatic.indels.vcf.gz | cut -f10 | paste -d ":" /dev/stdin indel.AFnormal.txt) \
            <(bcftools view -H -f 'PASS' somatic.indels.vcf.gz | cut -f11 | paste -d ":" /dev/stdin indel.AFtumor.txt)) | bgzip > somatic.indels.PASSED.vcf.gz

rm snp.AFnormal.txt 
rm snp.AFtumor.txt
rm indel.AFnormal.txt 
rm indel.AFtumor.txt            

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

## Combine the two intermediate files (vcfcombine from vcflib, Github)
echo '[MAIN] Combining SNVs/InDels, excluding LC variants'
vcfcombine <(bgzip -c -d somatic.snvs.PASSED.vcf.gz) <(bgzip -c -d somatic.indels.PASSED.vcf.gz) | \
  mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | tee somatic.PASSED.vcf | \
  bedtools intersect -header -sorted -v -a - -b $LC | \
  bgzip -@ 4 > somatic.PASSED.ExLC.vcf.gz

bgzip somatic.PASSED.vcf

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

## Remove common variants (do not pipe vep, as it aslways complains that STDOUT.summary or so already exists):
echo '[MAIN] Annotating with VEP & excluding common variants:'
vep --buffer_size 50000 --no_stats --cache --everything --fork 4 --vcf --format vcf --custom $DBSNP151,dbSNP151,vcf,exact,,COMMON,TOPMED \
  -i somatic.PASSED.ExLC.vcf.gz -o  STDOUT | \
    filter_vep --filter "not dbSNP151_COMMON = 1" | bgzip -@ 6 > somatic.PASSED.ExLC.Ex1KG.vcf.gz

echo '' && echo '[INFO] Ended on:' && date
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
