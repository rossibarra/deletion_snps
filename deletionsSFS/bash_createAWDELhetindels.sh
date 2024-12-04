#!/bin/bash

PATH=/gpool/bin/bedtools2/bin:$PATH

PREFIX=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.chronly.AWDELs.biallelicindels

zcat ${PREFIX}.vcf.gz | grep "0/1" > ${PREFIX}.hets.bed
zcat ${PREFIX}.vcf.gz | grep "1/0" >> ${PREFIX}.hets.bed
zcat ${PREFIX}.vcf.gz | grep "0/2" >> ${PREFIX}.hets.bed
zcat ${PREFIX}.vcf.gz | grep "2/0" >> ${PREFIX}.hets.bed
zcat ${PREFIX}.vcf.gz | grep "1/2" >> ${PREFIX}.hets.bed
zcat ${PREFIX}.vcf.gz | grep "2/1" >> ${PREFIX}.hets.bed

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${PREFIX}.hets.vcf
sort -k1,1V -k2,2n ${PREFIX}.hets.bed > ${PREFIX}.hets.sorted.bed
cat ${PREFIX}.hets.sorted.bed >> ${PREFIX}.hets.vcf
bedtools subtract -a ${PREFIX}.vcf.gz -b ${PREFIX}.hets.vcf > ${PREFIX}.homozygous.bed
