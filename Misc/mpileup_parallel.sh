#!/bin/bash

## Use samtools mpileup, parallelized over all primary chromosomes of a BAM file,
## and further chopped into chunks of 50000kb.
## Reason is that the larger chr (1,2,3) always take very long, so chunking should give quiet a gain
## in speed. Do max. of 36 threads in parallel, otherwise I/O will go crazy.
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
  grep -vE 'chrM|chrU|random|EBV|\*' | \
  cut -f1,2 | \
  bedtools makewindows -g /dev/stdin -w 50000000 | \
  awk '{print $1":"$2+1"-"$3}' | \
    parallel \
      -j 36 "samtools mpileup -Q 20 -B -d 1000 -f $GENOME \
        <(sambamba view -t 2 -f bam -l 0 $1 {}) | bgzip -@ 2 > ${1%.bam}_{}.mpileup.gz"

####################################################################################################
####################################################################################################

## Merge:
find ${1%.bam}_*.mpileup.gz | xargs cat > ${1%.bam}.mpileup.gz && rm ${1%.bam}_*.mpileup.gz

####################################################################################################
####################################################################################################
