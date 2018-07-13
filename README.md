# VCF_related
Scripts related to variant calling, filtering and annotation of VCF files.

#####################################################################################################################
#####################################################################################################################

# Using VarScan2_Pipeline.sh
This script calls variants from a tumor-normal pair, requiring only the two BAM files as input.
It must be started with a helper script that specifies:
1) the full path to the bam files
2) the basename for the output and
3) indicate the chromosomes for the parallelization.

That means, e.g.

cat helper_script.sh
```shell
#!/bin/bash
cat primary_chr.bed | \
  awk '{print $1":"$2+1"-"$3}' | \
  parallel "./varscan2.sh /path/to/tumor.bam /path/to/normal.bam Basename_{} {}"
```  
In this, the primary_chr.bed is a bed file with every chromosome to be parallelized over, with its start and end coordinate,
so start = 0 and end = length(chr). Say the Basename was FOO, then the output files would be:
FOO_chr1:1-248956422*.vcf (etc).

The list of primary chromosomes can be made with this one-liner:
```shell
samtools idxstats RG024-n_SortedRmdup.bam | \
  awk 'OFS="\t" {if ($1 == "*") next}; {print $1, "0", $2}' | \
  grep -v -E 'chrU|random|chrM' > primary_chr.bed
```  
The varscan pipeline will then create a folder ./VCF and several subfolders, where the final variants that passed the fpfilter will be present in ./VCF/processSomatic/VCF_High_Confidence/fpfilter_passed, split by chromosome.

