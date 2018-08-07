#!/bin/bash

#####################################################################################################

export LC_ALL=en_US.UTF-8

#####################################################################################################
####
#### Script calls variants from a single BAM file against reference genome using VarScan2,
#### and then applies the fpfilter.
#### ./VarScan2_Single.sh.sh in.bam samplename
####
#####################################################################################################
####
#### Tumor and Normal are supposed to be full paths to the bam files:
BAM=$(realpath $1)
BASENAME=$2
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

if [[ ! -e $BAM ]]; then echo '[ERROR]: Tumor BAM is missing - exiting' && exit 1; fi
if [[ ! -e ${BAM}.bai ]]; then
  echo '[ERROR]: BAM is not indexed - indexing now:'
  sambamba index -t 16 $1
fi

#####################################################################################################

echo '[MAIN]: Calling raw variants over all chromosomes:'
samtools idxstats $BAM | \
  grep -vE 'chrU|chrM|random|\*' | \
  mawk '$3 != "0" {print $1":2-"$2}' |
    parallel \
      "samtools mpileup -q 20 -Q 25 -B -d 1000 -f $HG38 \
        <(sambamba view -t 2 -f bam -l 0 $BAM {}) | \
          $VARSCAN2 mpileup2cns /dev/stdin ./VCF/${BASENAME}_{} --p-value 99e-02 --strand-filter 1 --output-vcf --variants"

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

#### BAM-READCOUNT & FPFILTER:
echo '[MAIN]: preparing bamRC BED files:' && echo ''

zcat ${BASENAME}_raw.vcf.gz | 
  grep -v '^#' | mawk 'OFS="\t" {print $1, $2-1, $2+1}' | \
  sort -k1,1 -k2,2n --parallel=16 | \
  bedtools merge -i - | bgzip -@ 2 > ${BASENAME}_bamRC.bed.gz

#####################################################################################################

echo '[MAIN]: Bam-Readcount:' && echo ''

if [[ -e ${BASENAME}_bamRC.bed.gz ]]; then
  bam-readcount -f $HG38 -q 20 -b 25 -d 1000 -l <(bgzip -c -d ${BASENAME}_bamRC.bed.gz) -w 1 $BAM | \
    bgzip -@ 6 > ${BASENAME}.bamRC.gz
    
  $VARSCAN2 fpfilter \
    <(bgzip -c -d -@ 2 ${BASENAME}_raw.vcf.gz) \
    <(bgzip -c -d -@ 2 ${BASENAME}.bamRC.gz) \
    --output-file ${BASENAME}_fppassed.vcf \
    --filtered-file ${BASENAME}_fpfailed.vcf \
    --dream3-settings && echo ''
  fi

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
