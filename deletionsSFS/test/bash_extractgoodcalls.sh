#!/bin/bash

MYFILEPREFIX=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps

# remove headers
#grep -v "^#" ${MYFILEPREFIX}.1000.vcf > ${MYFILEPREFIX}.1000.noheader.bed

# Grab first 10k results
zcat ${MYFILEPREFIX}.vcf.gz | grep -v "^#" | head -n 10000 > ${MYFILEPREFIX}.10k.vcf

# extract only B97 and convert to bed
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$12}' ${MYFILEPREFIX}.10k.vcf | awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5}' > ${MYFILEPREFIX}.10k.clean.bed

###
# extract deletions
grep -vP "\t0/0" ${MYFILEPREFIX}.10k.clean.bed | awk '$5 ~ /^\*\,[ATCG]$/' | grep -P "\t2/2" > ${MYFILEPREFIX}.10k.dels.bed
grep -vP "\t0/0" ${MYFILEPREFIX}.10k.clean.bed | awk '$5 ~ /^[ATCG]\,\*$/' | grep -P "\t1/1" >> ${MYFILEPREFIX}.10k.dels.bed

# extract non deletion calls
awk '$5 !~ /^\*\,[ATCG]$/' ${MYFILEPREFIX}.10k.clean.bed | awk '$5 !~ /^[ATCG]\,\*$/' > ${MYFILEPREFIX}.10k.nondels.bed

####
# extract only single allelic calls
#awk '$5=="A" || $5=="T" || $5=="C" || $5 == "G" {print $0}' ${MYFILEPREFIX}.10k.clean.bed > ${MYFILEPREFIX}.10k.singleallele.bed
awk '$5 ~ /^[ATCG]$/' ${MYFILEPREFIX}.10k.nondels.bed | grep -P "\t0/0" > ${MYFILEPREFIX}.10k.homozygous.bed
awk '$5 ~ /^[ATCG]$/' ${MYFILEPREFIX}.10k.nondels.bed | grep -P "\t1/1" >> ${MYFILEPREFIX}.10k.homozygous.bed

# extract missing data ./.
grep -P "\t\./\." ${MYFILEPREFIX}.10k.clean.bed | awk '$5 ~ /^\*\,[ATCG]$/' > ${MYFILEPREFIX}.10k.md.bed
grep -P "\t\./\." ${MYFILEPREFIX}.10k.clean.bed | awk '$5 ~ /^[ATCG]\,\*$/' >> ${MYFILEPREFIX}.10k.md.bed
grep -P "\t\./\." ${MYFILEPREFIX}.10k.nondels.bed | awk '$5 ~ /^[ATCG]$/' >> ${MYFILEPREFIX}.10k.md.bed

# sort bed files
sort -k1,1V -k2,2n ${MYFILEPREFIX}.10k.homozygous.bed > ${MYFILEPREFIX}.10k.homozygous.sorted.bed
sort -k1,1V -k2,2n ${MYFILEPREFIX}.10k.nondels.bed > ${MYFILEPREFIX}.10k.nondels.sorted.bed
sort -k1,1V -k2,2n ${MYFILEPREFIX}.10k.md.bed > ${MYFILEPREFIX}.10k.md.sorted.bed

# intersection of SNP sites with AW dels
bedtools intersect -a B97.deletions.merged.bed -b ${MYFILEPREFIX}.10k.homozygous.sorted.bed -wao > ${MYFILEPREFIX}.10k.ambiguouscalls.bed
