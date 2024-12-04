#!/bin/bash

PATH=/gpool/bin/bcftools-1.11/bin:$PATH

#MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf
#POSTFIX=snps.vcf
MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.nonbiallelicsnps.header.chr.vcf
POSTFIX=.nonbiallelicsnps.header.chr.vcf
#MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.vcf
#POSTFIX=snps.noAWDELs.vcf

#cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > $(basename ${MYVCF} ${POSTFIX}).biallelicindels.header.chr.bed
bcftools view -M 3 ${MYVCF} | grep "*" >> $(basename ${MYVCF} ${POSTFIX}).biallelicindels.header.chr.bed
#gzip $(basename ${MYVCF} ${POSTFIX}).biallelicindels.header.chr.bed
