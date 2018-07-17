#!/bin/bash

#####################################################################################################

export LC_ALL=en_US.UTF-8

#####################################################################################################
####
#### Bash script for variant calling with samtools mpileup | VarScan2.
#### Script assumes that a pair of tumor-normal BAM files is present in the dir where the script is.
#### Also assumes samtools, sambamba and bam-readcount in PATH, and path to varscan.jar specified
#### in this script as VARSCAN="java -jar /path/to/varscan.jar"
#### Raw variants are called with VarScan2 somatic, then split into somatic and germis by
#### processSomatic, selected for high-confidence variants and filtered for junk calls
#### with fpfilter. For fpfilter, bam-readcount is used to extract the necessary info
#### directly from the BAM files.
#### $4 is the region argument as chr:start-end to enable parallelization with GNU parallel,
#### USAGE: cat regions.txt | parallel "./VarScan2.sh tumor.bam normal.bam samplename_{} {}"
#### Variants per chromosome must then later be combined separately.
####
#####################################################################################################
####
#### Tumor and Normal are supposed to be full paths to the bam files:
TUMOR=$1
NORMAL=$2
BASENAME=$3
REGION=$4
####
VARSCAN2="java -jar /path/to/jar"
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

## Raw variants:
echo '[MAIN]: VarScan mpileup/somatic:'
samtools mpileup -q 20 -Q 25 -B -d 1000 -f $HG38 \
  <(samtools view -bu -@ 2 $NORMAL $REGION) \
  <(samtools view -bu -@ 2 $TUMOR $REGION) | \
    $VARSCAN2 somatic /dev/stdin ./VCF/${BASENAME} -mpileup --strand-filter 1 --output-vcf

cd ./VCF

#####################################################################################################

if [[ ! -d processSomatic ]]
  then
  mkdir processSomatic
  fi

## Use processSomatic to split into LOH, Germline and SNP, with and without high-confidence.

echo ''
echo '[MAIN]: VarScan processSomatic'
cat ${BASENAME}.snp.vcf | $VARSCAN2 processSomatic ./processSomatic/${BASENAME}.snp.vcf --max-normal-freq 0.01 
cat ${BASENAME}.indel.vcf | $VARSCAN2 processSomatic ./processSomatic/${BASENAME}.indel.vcf --max-normal-freq 0.01 

## Sort high-confidence variants into separate folder:
cd ./processSomatic

if [[ ! -d VCF_High_Confidence ]];
  then
  mkdir VCF_High_Confidence
  fi

mv ${BASENAME}*Somatic.hc* ./VCF_High_Confidence
mv ${BASENAME}*LOH.hc* ./VCF_High_Confidence
mv ${BASENAME}*Germline.hc* ./VCF_High_Confidence

cd ./VCF_High_Confidence
if [[ ! -d Germline ]];
  then
  mkdir Germline
  fi
  
mv ${BASENAME}*Germline* ./Germline

#####################################################################################################

#### BAM-READCOUNT & FPFILTER:

## Now concatenate all hc-VCF of one sample, sort and make a bed file of it for bam-readcount:
echo ''
echo '[MAIN]: Preparing region list for bam-readcount:' && echo ''
egrep -hv "^#" ${BASENAME}*.hc.vcf | \
  mawk 'OFS="\t" {print $1, $2-10, $2+10}' | \
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
bam-readcount -f $HG38 -q 20 -b 25 -d 1000 -l ${BASENAME}_bamRC_template.bed -w 1 $TUMOR | \
  bgzip -@ 6 > ${BASENAME}-t.bamRC.gz
bam-readcount -f $HG38 -q 20 -b 25 -d 1000 -l ${BASENAME}_bamRC_template.bed -w 1 $NORMAL | \
  bgzip -@ 6 > ${BASENAME}-n.bamRC.gz
echo ''

#####################################################################################################

## Apply fpfilter in case the bamRC files are not empty due to whatever reason:
mapfile -n 1 < <(bgzip -c -d ${BASENAME}-t.bamRC.gz)

if ((${#MAPFILE[@]} > 0)); then
  
  echo '[MAIN]: Applying false-positive filter:' && echo ''
  $VARSCAN2 fpfilter ${BASENAME}.snp.Somatic.hc.vcf <(bgzip -c -d -@ 6 ${BASENAME}-t.bamRC.gz) \
    --output-file ${BASENAME}.snp.Somatic.hc_passed.vcf --filtered-file ${BASENAME}.snp.Somatic.hc_failed.vcf
  echo ''

  $VARSCAN2 fpfilter ${BASENAME}.indel.Somatic.hc.vcf <(bgzip -c -d -@ 6 ${BASENAME}-t.bamRC.gz) \
    --output-file ${BASENAME}.indel.Somatic.hc_passed.vcf --filtered-file ${BASENAME}.indel.Somatic.hc_failed.vcf
  echo ''

  fi

mapfile -n 1 < <(bgzip -c -d ${BASENAME}-t.bamRC.gz)
  
if ((${#MAPFILE[@]} > 0)); then  

  $VARSCAN2 fpfilter ${BASENAME}.snp.LOH.hc.vcf <(bgzip -c -d -@ 6 ${BASENAME}-n.bamRC.gz) \
    --output-file ${BASENAME}.snp.LOH.hc_passed.vcf --filtered-file ${BASENAME}.snp.LOH.hc_failed.vcf
  echo ''

  $VARSCAN2 fpfilter ${BASENAME}.snp.LOH.hc.vcf <(bgzip -c -d -@ 6 ${BASENAME}-n.bamRC.gz) \
    --output-file ${BASENAME}.indel.LOH.hc_passed.vcf --filtered-file ${BASENAME}.indel.LOH.hc_failed.vcf
  echo ''

  fi
  
if [[ ! -d fpfilter_passed ]]
  then
  mkdir fpfilter_passed
  fi

mv ${BASENAME}*_passed.vcf ./fpfilter_passed

#####################################################################################################

echo '[FINISHED]: Pipeline finished on:' && date
echo '###############################################################################################'

#####################################################################################################
