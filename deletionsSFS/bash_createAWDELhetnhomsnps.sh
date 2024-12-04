#!/bin/bash

PATH=/gpool/bin/bedtools2/bin:$PATH

PREFIX=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.chronly.AWDELs.biallelic

zcat ${PREFIX}.vcf.gz | grep "0/1" > ${PREFIX}snps.hets.bed
zcat ${PREFIX}.vcf.gz | grep "1/0" >> ${PREFIX}snps.hets.bed
cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${PREFIX}snps.hets.sorted.vcf
sort -k1,1V -k2,2n ${PREFIX}snps.hets.bed >> ${PREFIX}snps.hets.sorted.vcf
bedtools subtract -a ${PREFIX}.vcf.gz -b ${PREFIX}snps.hets.sorted.vcf > ${PREFIX}snps.homozygous.bed

bash bash_reducevcftoshortvcf.sh ${PREFIX}snps.homozygous.bed .bed
