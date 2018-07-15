#!/bin/bash

## Helper script to launch VarScan2_Pipeline.sh, and to merge all somatic/germline variants at the end into one file per sample:

########################################################################################################################
#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=<int>
#SBATCH --partition=normal
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=foo@bar.com
#######
#######
BAMTUMOR="/path/to/tumor.bam"
BAMNORMAL="/path/to/normal.bam"
BASENAME="foobar"
#######
########################################################################################################################
########################################################################################################################

command -v samtools >/dev/null 2>&1 || { echo >&2 "[ERROR]: samtools is not in PATH"; exit 1; }

if [[ ! -e $BAMTUMOR ]]; then echo '[ERROR]: Tumor BAM is missing - exiting' && exit 1; fi
if [[ ! -e ${BAMTUMOR}.bai ]]; then
  echo '[ERROR]: Tumor BAM is not indexed - indexing now:'
  sambamba index -t 16 $1
  fi

if [[ ! -e $BAMNORMAL ]]; then echo '[ERROR]: Normal BAM is missing - exiting' && exit 1; fi
if [[ ! -e ${BAMNORMAL}.bai ]]; then
  echo '[ERROR]: Normal BAM is not indexed - indexing now:'
  sambamba index -t 16 $2
  fi
  
if [[ ! -e primary_chr.bed ]]; then
  samtools idxstats $BAMTUMOR | \
  awk 'OFS="\t" {if ($1 == "*") next}; {print $1, "0", $2}' | \
  grep -v -E 'chrU|random|chrM' > primary_chr.bed
  fi

########################################################################################################################
########################################################################################################################

## Parallelize variant calling by splitting up the task in number(chromosome) using GNU parallel:
cat primary_chr.bed | \
  awk '{print $1":"$2+1"-"$3}' | \
  parallel "./VarScan2_Pipeline.sh ${BAMTUMOR} ${BAMNORMAL} ${BASENAME}_{} {}"

########################################################################################################################
########################################################################################################################

## Merge all variants to one file per input sample:
cd ./fpfilter_passed
echo ''
echo '[INFO]: Merging variants to one file'

ls RG*.vcf | head -n 1 | xargs grep '#' > ${BASENAME}_header.txt
cat RG*.vcf | \
  mawk '$1 ~ /^#/ {next} {print $0 | "sort -k1,1 -k2,2n --parallel=8"}' | \
    cat ${BASENAME}_header.txt /dev/stdin | bgzip -@ 8 > ${BASENAME}_Somatic.vcf.gz
  
  tabix -p vcf ${BASENAME}_Somatic.vcf.gz &&
  rm ${BASENAME}*.vcf

exit

########################################################################################################################
