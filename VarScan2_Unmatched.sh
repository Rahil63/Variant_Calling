#!/bin/bash

export LC_ALL=en_US.UTF-8

#####################################################################################################
####
#### Bash script for variant calling with samtools mpileup | VarScan2.
#### Script assumes that a single/unmatched BAM file as input.
#### Also assumes samtools, sambamba and bam-readcount in PATH, and path to varscan.jar specified
#### in this script.
#### Raw variants are called with VarScan2 mpileup2cns, then split into somatic and germis by
#### 
#### with fpfilter. For fpfilter, bam-readcount is used to extract the necessary info
#### directly from the BAM files.
#### $4 is the region argument as chr:start-end to enable parallelization with GNU parallel,
#### USAGE: cat regions.txt | parallel "./VarScan2.sh tumor.bam normal.bam samplename_{} {}"
#### Variants per chromosome must then later be combined separately.
####
#####################################################################################################
####
#### Tumor and Normal are supposed to be full paths to the bam files:
BAM=$(realpath $1)
BASENAME=$2
####
VARSCAN2="java -jar -Xmx2g $HOME/software_2c/VarScan.v2.4.3.jar"
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

if [[ ! -e $BAM ]]; then echo '[ERROR]: BAM is missing - exiting' && exit 1; fi
if [[ ! -e ${BAM}.bai ]]; then
  echo '[INFO]: BAM is not indexed - indexing now:'
  sambamba index -t 16 $1
  fi

#####################################################################################################

## Raw variants in parallel over all chromosomes that have reads and are not U, random or M:
echo '[MAIN]: VarScan2'

samtools idxstats tmp.bam | 
  mawk '$3 != "0" {print $1":2-"$2}' | \
    parallel "samtools mpileup -q 20 -Q 25 -B -d 1000 -f $HG38 \
              <(samtools view -bu -@ 2 $BAM {}) | \
                $VARSCAN2 mpileup2cns --p-value 99e-02 --strand-filter 1 --output-vcf --variants > ./VCF/${BASENAME}_{}_raw.vcf"

cd ./VCF

## Merge variants per chomosome into one file:
echo '[MAIN]: Merging variants per chromosome into one file'
ls ${BASENAME}*.vcf | head -n 1 | xargs head -n 1000 | grep '#' | \
cat /dev/stdin ${BASENAME}*.vcf | \
  mawk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ${BASENAME}_raw.vcf


mapfile -n 1 < <(head -n 1000 ${BASENAME}_raw.vcf | grep -v '#' | awk NF)
if ((${#MAPFILE[@]}==0)); then
  echo '[WARNING]: Aborting script because no variants in file' ${BASENAME}_raw.vcf
  exit 0
fi  
  
#####################################################################################################

#### BAM-READCOUNT & FPFILTER:

## Now concatenate all hc-VCF of one sample, sort and make a bed file of it for bam-readcount:
echo ''
echo '[MAIN]: Preparing region list for bam-readcount:' && echo ''
egrep -hv "^#" ${BASENAME}_raw.vcf | \
  mawk 'OFS="\t" {print $1, $2-1, $2+1}' | \
  sort -k1,1 -k2,2n --parallel=8 | \
  bedtools merge -i - > ${BASENAME}_bamRC_template.bed

## Abort if template file is empty:
mapfile -n 1 < ${BASENAME}_bamRC_template.bed
if ((${#MAPFILE[@]}==0)); then
  echo '[WARNING]: Aborting script because no variants made it to this point (template for bamRC)'
  exit 0
  fi

## Now get data from bam-readcount for both the tumor and normal file:
echo '[MAIN]: Getting data from bam-readcount:' && echo ''
bam-readcount -f $HG38 -q 20 -b 25 -d 1000 -l ${BASENAME}_bamRC_template.bed -w 1 ${BAM} | \
  bgzip -@ 6 > ${BASENAME}.bamRC.gz
echo ''

#####################################################################################################

## Apply fpfilter in case the bamRC files are not empty due to whatever reason:
mapfile -n 1 < <(bgzip -c -d ${BASENAME}.bamRC.gz | awk NF)

if ((${#MAPFILE[@]} > 0)); then
  
  echo '[MAIN]: Applying false-positive filter:' && echo ''
  $VARSCAN2 fpfilter ${BASENAME}_raw.vcf <(bgzip -c -d -@ 6 ${BASENAME}.bamRC.gz) \
    --output-file ${BASENAME}_fppassed.vcf --filtered-file ${BASENAME}_failed.vcf
  echo ''

  fi

if [[ ! -d ./fpfilter_passed ]]; then mkdir fpfilter_passed; fi
mv ${BASENAME}_fppassed.vcf ./fpfilter_passed

#####################################################################################################

echo '[FINISHED]: Pipeline finished on:' && date
echo '###############################################################################################'

#####################################################################################################
