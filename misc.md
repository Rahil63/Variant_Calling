# some snippets (...) 

## Get the chromosomes that have entries:
´´´
tabix -l ${BASENAME}.Somatic*.gz > ${BASENAME}_chr_somatic.txt
tabix -l ${BASENAME}.LOH*.gz > ${BASENAME}_chr_loh.txt
tabix -l ${BASENAME}.Germline*.gz > ${BASENAME}_chr_germ.txt

## Split the combined VCF into VCF-by-chromosome:
cat ${BASENAME}_chr_somatic.txt | parallel "tabix -h ${BASENAME}.Somatic.hc.vcf.gz {} | bgzip -@ 2 > ${BASENAME}_{}.Somatic.vcf.gz"
cat ${BASENAME}_chr_loh.txt | parallel "tabix -h ${BASENAME}.LOH.hc.vcf.gz {} | bgzip -@ 2 > ${BASENAME}_{}.LOH.vcf.gz"
cat ${BASENAME}_chr_germ.txt | parallel "tabix -h ${BASENAME}.Germline.hc.vcf.gz {} | bgzip -@ 2 > ${BASENAME}_{}.Germline.vcf.gz"

## Template for bam-readcount:
echo ''
echo '[MAIN]: Preparing region list for bam-readcount:' && echo ''
for i in ${BASENAME}_*.vcf.gz 
  do
  BNAME=${i%.vcf.gz}
  bgzip -c -d $i | grep -v '^#' | mawk 'OFS="\t" {print $1, $2-1, $2+1}' | \
  sort -k1,1 -k2,2n --parallel=16 | \
  bedtools merge -i - | bgzip -@ 2 > ${BNAME}_template_bamRC.bed.gz
  done

## Now get data from bam-readcount for both the tumor and normal file:
if [[ ! $(ls ${BASENAME}*.Somatic*bamRC*.gz | wc -l) -eq 0 ]]; then 
  
  ls ${BASENAME}*.Somatic*bamRC*.gz | \
    parallel "bgzip -c -d {} 
