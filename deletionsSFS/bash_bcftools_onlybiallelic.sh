#!/bin/bash

PATH=/gpool/bin/bcftools-1.11/bin:$PATH

MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf
POSTFIX=snps.vcf
#MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.vcf
#POSTFIX=snps.noAWDELs.vcf

bcftools view --types snps -m 2 -M 2 ${MYVCF}.gz > $(basename ${MYVCF} ${POSTFIX})biallelic${POSTFIX}
gzip $(basename ${MYVCF} ${POSTFIX})biallelic${POSTFIX}
