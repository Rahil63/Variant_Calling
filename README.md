# VCF_related
Scripts related to variant calling, filtering and annotation of VCF files.

##########################################################################################
##########################################################################################

## Using VarScan2_Somatic.sh
This script calls variants from a tumor-normal pair, requiring only the two BAM files as input plus a samplename.
In the script itself, one must specify the path to the genome.fa and the path to the varscan.jar executable.

```shell
./VarScan2_Somatic.sh tumor.bam normal.bam samplename
```

## Filter_VCF_Common.sh
This is a script that reads a VCF file and then uses the Ensemble Variant Effect Predictor together with dbSNP151 to remove variants that are flagged as common variant in the 1KG (means MAF of at least 1% in any of the ethnicity groups, indicated by COMMON=1) or an AF for ALT in TOPMED greater 1%.

## Split_VCF_By_Chr.sh
This script uses Tabix and parallel to split a VCF file into one VCF per chromosome.
