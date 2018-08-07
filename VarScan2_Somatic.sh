#!/bin/bash

#####################################################################################################

export LC_ALL=en_US.UTF-8

#####################################################################################################
####
#### Bash script for variant calling with samtools mpileup | VarScan2.
#### Script assumes that a pair of tumor-normal BAM files is present in the dir where the script is started from.
#### Also assumes samtools, sambamba and bam-readcount in PATH, and path to varscan.jar specified
#### in this script as VARSCAN="java -jar /path/to/varscan.jar"
#### Raw variants are called with VarScan2 somatic, then split into somatic and germs by
#### processSomatic, selected for high-confidence variants and filtered for junk calls
#### with fpfilter. For fpfilter, bam-readcount is used to extract the necessary info
#### directly from the BAM files.
#### GNU parallel is used to parallelize the job over all chromosomes.
#### ./VarScan2.sh tumor.bam normal.bam samplename
####
#####################################################################################################
####
#### Tumor and Normal are supposed to be full paths to the bam files:
TUMOR=$(realpath $1)
NORMAL=$(realpath $2)
BASENAME=$3
####
VARSCAN2="java -jar -Xmx5g $HOME/software_2c/VarScan.v2.4.3.jar"
HG38="/scratch/tmp/a_toen03/Genomes/hg38/hg38_noALT_withDecoy.fa"
####
#####################################################################################################

## Check if programs are in path:
command -v samtools >/dev/null 2>&1 || { echo >&2 "[ERROR]: samtools is not in PATH"; exit 1; }
command -v sambamba >/dev/null 2>&1 || { echo >&2 "[ERROR]: sambamba is not in PATH"; exit 1; }
command -v bam-readcount >/dev/null 2>&1 || { echo >&2 "[ERROR]: bam-readcount is not in PATH"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo >&2 "[ERROR]: bedtools is not in PATH"; exit 1; }
command -v mawk >/dev/null 2>&1 || { echo >&2 "[ERROR]: mawk is not in PATH"; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo >&2 "[ERROR]: bgzip is not in PATH"; exit 1; }

#####################################################################################################

echo ''
echo '[INFO]: Pipeline for' $BASENAME 'started' && date && echo ''

#####################################################################################################

## Mpileup from samtools piped into VarScan:

## Write output files into a directory termed VCF:

if [[ ! -d ./VCF ]]; then mkdir VCF; fi

#####################################################################################################

if [[ ! -e $TUMOR ]]; then echo '[ERROR]: Tumor BAM is missing - exiting' && exit 1; fi
if [[ ! -e ${TUMOR}.bai ]]; then
  echo '[ERROR]: Tumor BAM is not indexed - indexing now:'
  sambamba index -t 16 $1
  fi

if [[ ! -e $NORMAL ]]; then echo '[ERROR]: Normal BAM is missing - exiting' && exit 1; fi
if [[ ! -e ${NORMAL}.bai ]]; then
  echo '[ERROR]: Normal BAM is not indexed - indexing now:'
  sambamba index -t 16 $2
  fi

#####################################################################################################

echo '[MAIN]: Calling raw variants over all chromosomes:'
samtools idxstats $TUMOR | \
  grep -vE 'chrU|chrM|random|\*' | \
  mawk '$3 != "0" {print $1":2-"$2}' |
    parallel \
      "samtools mpileup -q 20 -Q 25 -B -d 1000 -f $HG38 \
        <(sambamba view -t 2 -f bam -l 0 $NORMAL {}) \
        <(sambamba view -t 2 -f bam -l 0 $TUMOR {}) | \
          $VARSCAN2 somatic /dev/stdin ./VCF/${BASENAME}_{} -mpileup --strand-filter 1 --output-vcf"

cd ./VCF

#####################################################################################################

## Combine snp and indel calls:
echo ''
echo '[MAIN]: Combining raw variants into one file'

ls ${BASENAME}*.vcf | head -n 1 | xargs grep '#' > ${BASENAME}_header.txt

cat ${BASENAME}*.vcf | \
  mawk '$1 ~ /^#/ {next} {print $0 | "sort -k1,1 -k2,2n --parallel=8"}' | \
  cat ${BASENAME}_header.txt /dev/stdin | bgzip -@ 8 > ${BASENAME}_raw.vcf.gz

#####################################################################################################

if [[ ! -d processSomatic ]]
  then
  mkdir processSomatic
  fi

## Use processSomatic to split into LOH, Germline and SNP, with and without high-confidence.
echo ''
echo '[MAIN]: VarScan processSomatic'
bgzip -c -d ${BASENAME}_raw.vcf.gz | \
  $VARSCAN2 processSomatic ./processSomatic/${BASENAME}.vcf

## Sort high-confidence variants into separate folder and compress:
cd ./processSomatic

ls ${BASENAME}*.vcf | parallel "bgzip -@ 4 {}"

if [[ ! -d VCF_High_Confidence ]];
  then
  mkdir VCF_High_Confidence
  fi

mv ${BASENAME}*Somatic.hc* ./VCF_High_Confidence
mv ${BASENAME}*LOH.hc* ./VCF_High_Confidence
mv ${BASENAME}*Germline.hc* ./VCF_High_Confidence

cd ./VCF_High_Confidence

ls ${BASENAME}*.vcf.gz | parallel "tabix -p vcf {}"
echo ''
#####################################################################################################

#### BAM-READCOUNT & FPFILTER:
echo '[MAIN]: preparing bamRC BED files:' && echo ''

zcat ${BASENAME}.Somatic.hc.vcf.gz | 
  grep -v '^#' | mawk 'OFS="\t" {print $1, $2-1, $2+1}' | \
  sort -k1,1 -k2,2n --parallel=16 | \
  bedtools merge -i - | bgzip -@ 2 > ${BASENAME}.Somatic_bamRC.bed.gz
  
zcat ${BASENAME}.Germline.hc.vcf.gz | 
  grep -v '^#' | mawk 'OFS="\t" {print $1, $2-1, $2+1}' | \
  sort -k1,1 -k2,2n --parallel=16 | \
  bedtools merge -i - | bgzip -@ 2 > ${BASENAME}.Germline_bamRC.bed.gz
  
zcat ${BASENAME}.LOH.hc.vcf.gz | 
  grep -v '^#' | mawk 'OFS="\t" {print $1, $2-1, $2+1}' | \
  sort -k1,1 -k2,2n --parallel=16 | \
  bedtools merge -i - | bgzip -@ 2 > ${BASENAME}.LOH_bamRC.bed.gz  

#####################################################################################################

echo '[MAIN]: Bam-Readcount:' && echo ''

if [[ -e ${BASENAME}.Somatic_bamRC.bed.gz ]]; then
  bam-readcount -f $HG38 -q 20 -b 25 -d 1000 -l <(bgzip -c -d ${BASENAME}.Somatic_bamRC.bed.gz) -w 1 $TUMOR | \
    bgzip -@ 6 > ${BASENAME}-Somatic.bamRC.gz
    
  $VARSCAN2 fpfilter \
    <(bgzip -c -d -@ 2 ${BASENAME}.Somatic.hc.vcf.gz) \
    <(bgzip -c -d -@ 2 ${BASENAME}-Somatic.bamRC.gz) \
    --output-file ${BASENAME}.Somatic.hc_fppassed.vcf \
    --filtered-file ${BASENAME}.Somatic.hc_fpfailed.vcf \
    --dream3-settings && echo ''
  fi

if [[ ${BASENAME}.Germline_bamRC.bed.gz ]]; then
  bam-readcount -f $HG38 -q 20 -b 25 -d 1000 -l <(bgzip -c -d ${BASENAME}.Germline_bamRC.bed.gz) -w 1 $NORMAL | \
    bgzip -@ 6 > ${BASENAME}-Germline.bamRC.gz
  
  $VARSCAN2 fpfilter \
    <(bgzip -c -d -@ 2 ${BASENAME}.Germline.hc.vcf.gz) \
    <(bgzip -c -d -@ 2 ${BASENAME}-Germline.bamRC.gz) \
    --output-file ${BASENAME}.Germline.hc_fppassed.vcf \
    --filtered-file ${BASENAME}.Germline.hc_fpfailed.vcf \
    --dream3-settings && echo ''
  fi

if [[ -e ${BASENAME}.LOH_bamRC.bed.gz ]]; then
  bam-readcount -f $HG38 -q 20 -b 25 -d 1000 -l <(bgzip -c -d ${BASENAME}.LOH_bamRC.bed.gz) -w 1 $NORMAL | \
    bgzip -@ 6 > ${BASENAME}-LOH.bamRC.gz
  
  $VARSCAN2 fpfilter \
    <(bgzip -c -d -@ 2 ${BASENAME}.LOH.hc.vcf.gz) \
    <(bgzip -c -d -@ 2 ${BASENAME}-LOH.bamRC.gz) \
    --output-file ${BASENAME}.LOH.hc_fppassed.vcf \
    --filtered-file ${BASENAME}.LOH.hc_fpfailed.vcf \
    --dream3-settings && echo ''
  fi  
  
echo ''

if [[ ! -d fpfilter_passed ]]
  then
  mkdir fpfilter_passed
  fi

ls *.vcf | parallel "bgzip {}"
ls *fppassed.vcf.gz | parallel "tabix -p vcf {}"

mv ${BASENAME}*_fppassed.vcf* ./fpfilter_passed

#####################################################################################################

echo '[FINISHED]: Pipeline for' $BASENAME 'finished on:' && date
echo '###############################################################################################'

#####################################################################################################
