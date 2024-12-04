#!/bin/bash

PATH=/gpool/bin/bcftools-1.11/bin:/gpool/bin/bedtools2/bin:$PATH

MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.chronly.AWDELs
POSTFIX=

#bedtools subtract -a B73v5.NAM.chronly.vcf.gz -b B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.headers.vcf.gz > ${MYVCF}.bed
#cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${MYVCF}.vcf
#cat ${MYVCF}.bed >> ${MYVCF}.vcf
#gzip ${MYVCF}.bed
#gzip ${MYVCF}.vcf

# extract biallelic snsps
bcftools view --types snps -m 2 -M 2 ${MYVCF}.vcf.gz > ${MYVCF}.biallelic.vcf
gzip ${MYVCF}.biallelic.vcf

# do some weird shit Jeff wants. still have no freaking clue...


# extract biallelic indels
cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${MYVCF}.biallelicindels.vcf
bcftools view -M 3 -m 3 ${MYVCF}.vcf.gz | grep "*" >> ${MYVCF}.biallelicindels.vcf

# do some weird shit Jeff wants here too. Still have no freaking clue what he wants :(
# Should I combine the two then apply python?
