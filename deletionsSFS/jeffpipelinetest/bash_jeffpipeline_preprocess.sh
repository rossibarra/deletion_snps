#!/bin/bash

MYVCF=B73v5.NAM.chronly.vcf.gz
NAMLINE=Tzi8

PATH=/gpool/bin/bcftools-1.10.2/bin:$PATH
#PATH=/gpool/bin/bedtools2/bin:$PATH

bcftools view -M2 -v snps ${MYVCF} > NAM.biallelic.snpsonly.vcf
cat ../B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > NAM.biallelic.indelsonly.vcf
bcftools view -M3 ${MYVCF} | grep "*" >> NAM.biallelic.indelsonly.vcf

bcftools concat NAM.biallelic.snpsonly.vcf NAM.biallelic.indelsonly.vcf > NAM.goodSNPs.vcf
# Run one or the other. Latter is preferred
#bcftools filter -S . -e 'GT=="0/1" || GT=="0/2" || GT=="1/2"' NAM.goodSNPs.vcf > NAM.goodSNPs_hetstomissing.vcf
bcftools filter -S . -e 'GT=="het"' NAM.goodSNPs.vcf > NAM.goodSNPs_hetstomissing.vcf
