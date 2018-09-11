#!/bin/bash

## Use samtools mpileup, parallelized over all primary chromosomes of a BAM file.
## Usage: ./mpileup_paralle.sh foo.bam
## Output would be foo.mpileup.gz

####################################################################################################
####################################################################################################

GENOME="/scratch/tmp/a_toen03/Genomes/hg38/hg38_noALT_withDecoy.fa"

####################################################################################################
####################################################################################################

if [[ -e ${1%.bam}.mpileup.gz ]]; then
  exit
  fi

if [[ $# -ne 0 ]] ; then
  echo '[USAGE]: mpileup_parallel.sh in.bam'
  fi

## Mpileup parallel over chromosomes using GNU parallel:
samtools idxstats $1 | \
  awk 'OFS="\t" {print $1, "0", $2}' | \
  grep -vE 'chrM|chrU|random|EBV|\*' | \
  awk '{print $1":"$2+1"-"$3}' | \
    parallel \
      "samtools mpileup -Q 20 -B -d 1000 -f $GENOME \
        <(sambamba view -t 2 -f bam -l 0 $1 {}) | bgzip -@ 2 > ${1%.bam}_{}.mpileup.gz"

####################################################################################################
####################################################################################################

## Merge:
find ${1%.bam}_*.mpileup.gz | xargs cat > ${1%.bam}.mpileup.gz && rm ${1%.bam}_*.mpileup.gz

####################################################################################################
####################################################################################################
