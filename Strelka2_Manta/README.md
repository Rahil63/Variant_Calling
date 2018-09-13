# Strelka2-Manta
This repo contains scripts to run Strelka2 and Manta on tumor/normal BAM file pairs,
plus scripts to do postprocessing.

## Define_Paths.sh
Contains all the paths to executables and files needed to run Strelka and Manta + postprocessing.

## Make_MantaDirectories.sh
This script creates the directories for Manta SV calling, pulling the basename of the BAM files as dir name.
It creates a subfolder <Somatic> or <Germline> in each of these. 
Can simply be run in the same folder as the BAM files without any options.
  
## Make_StrelkaDirectories.sh  
Same as the ```Make_MantaDirectories.sh```but for Strelka2. Can only be run after Manta workflow has 
finished because the ```candidateSmallIndels.vcf.gz``` must be passed to Strelka2.

## Strelka2_postprocessing.sh
This script postprocesses the somatic calls from Strelka2.
Can be executed via ```start_Strelka2_postprocessing.sh``` from the base directory.
Will remove non-PASSED variants, variants in low-complexity regions,
variants that are COMMON=1 in dbSNP151/1KG, add allele frequency to the VCF,
and create one merged VCF.gz for each sample.
A symbolic link for each of these files is created in ```./Final_Variants```.

## To Do
Write something that postprocesses the Manta calls, excluding or blacklisting breakpoints that are
in LC regions. Problem is that standard bedtools intersect always intersects the entire interval rather than only the breakpoints.

