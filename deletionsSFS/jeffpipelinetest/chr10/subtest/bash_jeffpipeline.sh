#!/bin/bash

XSOME=$1
#MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.singletest.vcf
MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.jeffsubset.vcf
#MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.chr10.vcf.gz

PATH=/gpool/bin/bcftools-1.10.2/bin:$PATH

bcftools view -M2 -v snps ${MYVCF} > NAM.biallelic.snpsonly.vcf
cat ../B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > NAM.biallelic.indelsonly.vcf
bcftools view -M3 ${MYVCF} | grep "*" >> NAM.biallelic.indelsonly.vcf

bcftools concat NAM.biallelic.snpsonly.vcf NAM.biallelic.indelsonly.vcf > NAM.goodSNPs.vcf
#bcftools filter -S . -e 'GT=="0/1" || GT=="0/2" || GT=="1/2"' NAM.goodSNPs.vcf > NAM.goodSNPs_hetstomissing.vcf
bcftools filter -S . -e 'GT=="het"' NAM.goodSNPs.vcf > NAM.goodSNPs_hetstomissing.vcf

bcftools view -s B97 NAM.goodSNPs_hetstomissing.vcf > b97.vcf
bedtools intersect -a B97.deletions.merged.${XSOME}.bed -b b97.vcf -wao > b97.deletions.snps

grep -v "0$" b97.deletions.snps | grep -v "\./\." > b97.overlap.deletions.snps

grep 0/0 b97.overlap.deletions.snps > bad.deletions
grep "\,\*" b97.overlap.deletions.snps | grep 1/1 >> bad.deletions
grep "*," b97.overlap.deletions.snps| grep "2/2" >> bad.deletions
