#!/bin/bash

NAMLINE=$1
XSOME=$2

#MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.singletest.vcf
#MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.chr10.vcf.gz

PATH=/gpool/bin/bcftools-1.10.2/bin:$PATH
PATH=/gpool/bin/bedtools2/bin:$PATH

#THE BIG JOB. NOT NAMLINE DEPENDENT
#bcftools view -M2 -v snps ${MYVCF} > NAM.biallelic.snpsonly.vcf
#cat ../B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > NAM.biallelic.indelsonly.vcf
#bcftools view -M3 ${MYVCF} | grep "*" >> NAM.biallelic.indelsonly.vcf

#bcftools concat NAM.biallelic.snpsonly.vcf NAM.biallelic.indelsonly.vcf > NAM.goodSNPs.vcf
#bcftools filter -S . -e 'GT=="0/1" || GT=="0/2" || GT=="1/2"' NAM.goodSNPs.vcf > NAM.goodSNPs_hetstomissing.vcf
#bcftools filter -S . -e 'GT=="het"' NAM.goodSNPs.vcf > NAM.goodSNPs_hetstomissing.vcf

#BEGIN NAMLINE DEPENDENT STUFF
bcftools view -s ${NAMLINE} NAM.goodSNPs_hetstomissing.vcf > ${NAMLINE}.vcf

bedtools intersect -a ${NAMLINE}.deletions.merged.bed -b ${NAMLINE}.vcf -wao > ${NAMLINE}.deletions.snps

grep -v "0$" ${NAMLINE}.${XSOME}.deletions.snps | grep -v "\./\." > ${NAMLINE}.overlap.deletions.snps

grep 0/0 ${NAMLINE}.overlap.deletions.snps > ${NAMLINE}.bad.deletions
grep "\,\*" ${NAMLINE}.overlap.deletions.snps | grep 1/1 >> ${NAMLINE}.bad.deletions
grep "*," ${NAMLINE}.overlap.deletions.snps| grep "2/2" >> ${NAMLINE}.bad.deletions

#BY XSOME
#bedtools intersect -a ${NAMLINE}.deletions.merged.${XSOME}.bed -b ${NAMLINE}.vcf -wao > ${NAMLINE}.${XSOME}.deletions.snps

#grep -v "0$" ${NAMLINE}.${XSOME}.deletions.snps | grep -v "\./\." > ${NAMLINE}.${XSOME}.overlap.deletions.snps

#grep 0/0 ${NAMLINE}.${XSOME}.overlap.deletions.snps > ${NAMLINE}.${XSOME}.bad.deletions
#grep "\,\*" ${NAMLINE}.${XSOME}.overlap.deletions.snps | grep 1/1 >> ${NAMLINE}.${XSOME}.bad.deletions
#grep "*," ${NAMLINE}.${XSOME}.overlap.deletions.snps| grep "2/2" >> ${NAMLINE}.${XSOME}.bad.deletions
