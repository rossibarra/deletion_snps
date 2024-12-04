#!/bin/bash

PATH=/gpool/bin/bcftools-1.10.2/bin:$PATH
PATH=/gpool/bin/bedtools2/bin:$PATH

NAMLINE=$1

#zcat B73v5.NAM.chronly.vcf.gz | grep -v "^#" | wc -l
echo "Total SNPs 27,256,799" > ${NAMLINE}.txt
echo "Next values are for line specific SNPs" >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

bcftools view -s ${NAMLINE} B73v5.NAM.chronly.vcf.gz > ${NAMLINE}.only.vcf
NAMLINESNPS=$(grep -v "^#" ${NAMLINE}.only.vcf | wc -l)
echo "${NAMLINE} only SNPs" >> ${NAMLINE}.txt
echo ${NAMLINESNPS} >> ${NAMLINE}.txt
#27,256,799
echo "" >> ${NAMLINE}.txt

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.MD.vcf
grep -v "^#" ${NAMLINE}.only.vcf | grep -a "\./." >> ${NAMLINE}.only.MD.vcf
NAMLINEMDS=$(grep -v "^#" ${NAMLINE}.only.MD.vcf | wc -l)
#5,756,977
echo "MDs" >> ${NAMLINE}.txt
echo ${NAMLINEMDS} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.vcf
grep -v "\./." ${NAMLINE}.only.vcf > ${NAMLINE}.only.noMD.vcf
NAMLINENOMDS=$(grep -v "^#" ${NAMLINE}.only.noMD.vcf | wc -l)
#21,499,822
echo "no MDs" >> ${NAMLINE}.txt
echo ${NAMLINENOMDS} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.muliallelic.vcf
grep -v "^#" ${NAMLINE}.only.noMD.vcf | grep "[*ATCG]," >> ${NAMLINE}.only.noMD.muliallelic.vcf
NAMLINEMULTIALLELIC=$(grep -v "^#" ${NAMLINE}.only.noMD.muliallelic.vcf | wc -l)
#720,833
echo "No MDs & Multiallelic. Includes indels" >> ${NAMLINE}.txt
echo ${NAMLINEMULTIALLELIC} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.biallelic.vcf
grep -v "^#" ${NAMLINE}.only.noMD.vcf | grep -v "[*ATCG]," >> ${NAMLINE}.only.noMD.biallelic.vcf
#20,778,989
NAMLINEBIALLELIC=$(grep -v "^#" ${NAMLINE}.only.noMD.biallelic.vcf | wc -l)
echo "Biallelic SNPs" >> ${NAMLINE}.txt
echo ${NAMLINEBIALLELIC} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

#bcftools view -M2 -v snps ${NAMLINE}.only.noMD.vcf > ${NAMLINE}.only.noMD.biallelic.vcf
#grep -v "^#" ${NAMLINE}.only.noMD.biallelic.vcf | wc -l
#20,778,989

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.indels.vcf
grep -v "^#" ${NAMLINE}.only.noMD.muliallelic.vcf | grep "\*" >> ${NAMLINE}.only.noMD.indels.vcf
NAMLINEINDELS=$(grep -v "^#" ${NAMLINE}.only.noMD.indels.vcf | wc -l)
#414,965
echo "Multiallelic Indels" >> ${NAMLINE}.txt
echo ${NAMLINEINDELS} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

echo "Multiallelic but not indels" >> ${NAMLINE}.txt
echo $((NAMLINEMULTIALLELIC - NAMLINEINDELS)) >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

#cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.indels.vcf
#bcftools view -M3 ${NAMLINE}.only.noMD.vcf | grep "\*" >> ${NAMLINE}.only.noMD.indels.vcf
#grep -v "^#" ${NAMLINE}.only.noMD.indels.vcf | wc -l
#405,725
#echo >> ${NAMLINE}.txt
#>> ${NAMLINE}.txt

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.biallelic.AWDelovl.vcf
bedtools intersect -a ${NAMLINE}.only.noMD.biallelic.vcf -b ../${NAMLINE}.deletions.sorted.bed >> ${NAMLINE}.only.noMD.biallelic.AWDelovl.vcf
NAMLINEBIALLELICAWDOVL=$(grep -v "^#" ${NAMLINE}.only.noMD.biallelic.AWDelovl.vcf | wc -l)
#2,068,366
echo "Biallelic SNPs overlapping with AWDels" >> ${NAMLINE}.txt
echo ${NAMLINEBIALLELICAWDOVL} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

### homozygous
grep -v "^#" ${NAMLINE}.only.noMD.biallelic.vcf | grep "0/0" > ${NAMLINE}.noMD.homo.vcf
grep -v "^#" ${NAMLINE}.only.noMD.biallelic.vcf | grep "1/1" >> ${NAMLINE}.noMD.homo.vcf
#20,004,204
NAMLINEHOMO=$(grep -v "^#" ${NAMLINE}.noMD.homo.vcf | wc -l)
echo "Homozygous biallelic SNPs" >> ${NAMLINE}.txt
echo ${NAMLINEHOMO} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.biallelic.hom.vcf
sort -k1,1V -k2,2n ${NAMLINE}.noMD.homo.vcf >> ${NAMLINE}.only.noMD.biallelic.hom.vcf
grep -v "^#" ${NAMLINE}.only.noMD.biallelic.hom.vcf | wc -l
#20,004,204

# calc intergenic and genic:
bedtools intersect -a ${NAMLINE}.only.noMD.biallelic.hom.vcf -b ../../annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.genes.gtf > ${NAMLINE}.only.noMD.hom.genic.vcf
NAMLINEHOMOGENIC=$(grep -v "^#" ${NAMLINE}.only.noMD.hom.genic.vcf | wc -l)
echo "Genic Homozygous biallelic SNPs" >> ${NAMLINE}.txt
echo ${NAMLINEHOMOGENIC} >> ${NAMLINE}.txt
#genic 2,938,176
echo "" >> ${NAMLINE}.txt

NAMLINEHOMOINTERGENIC=$((NAMLINEHOMO - NAMLINEHOMOGENIC))
echo "Intergenic Homozygous biallelic SNPs" >> ${NAMLINE}.txt
echo ${NAMLINEHOMOINTERGENIC} >> ${NAMLINE}.txt
#intergenic 17,066,028
echo "" >> ${NAMLINE}.txt

### hom awdels

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.biallelic.hom.AWDelovl.vcf
bedtools intersect -a ${NAMLINE}.only.noMD.biallelic.hom.vcf -b ../${NAMLINE}.deletions.sorted.bed >> ${NAMLINE}.only.noMD.biallelic.hom.AWDelovl.vcf
NAMLINEHOMAWDEL=$(grep -v "^#" ${NAMLINE}.only.noMD.biallelic.hom.AWDelovl.vcf | wc -l)
#2,212,624
#1,910,506
echo "Homozygous biallelic SNPs with AWDelovl" >> ${NAMLINE}.txt
echo ${NAMLINEHOMAWDEL} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.biallelic.hom.AWDelovl.genic.vcf
bedtools intersect -a ${NAMLINE}.only.noMD.biallelic.hom.AWDelovl.vcf -b ../../annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.genes.gtf >> ${NAMLINE}.only.noMD.biallelic.hom.AWDelovl.genic.vcf
NAMLINEAWDELGENIC=$(grep -v "^#" ${NAMLINE}.only.noMD.biallelic.hom.AWDelovl.genic.vcf | wc -l)
#122,872
#NAMLINEAWDELGENIC
echo "Genic Homozygous biallelic SNPs with AWDelovl" >> ${NAMLINE}.txt
echo ${NAMLINEAWDELGENIC} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

echo "Intergenic Homozygous biallelic SNPs with AWDelovl" >> ${NAMLINE}.txt
echo $((NAMLINEHOMAWDEL - NAMLINEAWDELGENIC)) >> ${NAMLINE}.txt
#intergenic
#1,910,506 - 122,872
#1,787,634
echo "" >> ${NAMLINE}.txt

### Het calls
grep -v "^#" ${NAMLINE}.only.noMD.biallelic.vcf | grep "0/1" > ${NAMLINE}.noMD.hets.vcf
grep -v "^#" ${NAMLINE}.only.noMD.biallelic.vcf | grep "1/0" >> ${NAMLINE}.noMD.hets.vcf
#774,785
NAMLINEHET=$(grep -v "^#" ${NAMLINE}.noMD.hets.vcf | wc -l)
echo "Heterozygous biallelic SNPs" >> ${NAMLINE}.txt
echo ${NAMLINEHET} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

### het awdels

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.biallelic.het.vcf
sort -k1,1V -k2,2n ${NAMLINE}.noMD.hets.vcf >> ${NAMLINE}.only.noMD.biallelic.het.vcf
#grep -v "^#" ${NAMLINE}.only.noMD.biallelic.het.vcf | wc -l
#774,785

# calc intergenic and genic:
bedtools intersect -a ${NAMLINE}.only.noMD.biallelic.het.vcf -b ../../annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.genes.gtf > ${NAMLINE}.only.noMD.het.genic.vcf
NAMLINEHETGENIC=$(grep -v "^#" ${NAMLINE}.only.noMD.het.genic.vcf | wc -l)
echo "Genic Heterozygous biallelic SNPs" >> ${NAMLINE}.txt
echo ${NAMLINEHETGENIC} >> ${NAMLINE}.txt
#genic: 90,204
echo "" >> ${NAMLINE}.txt

NAMLINEHETINTERGENIC=$((NAMLINEHET - NAMLINEHETGENIC))
echo "Intergenic Heterozygous biallelic SNPs" >> ${NAMLINE}.txt
echo ${NAMLINEHETINTERGENIC} >> ${NAMLINE}.txt
#intergenic: 684,581
echo "" >> ${NAMLINE}.txt

### AWDel calls

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.biallelic.het.AWDelovl.vcf
bedtools intersect -a ${NAMLINE}.only.noMD.biallelic.het.vcf -b ../${NAMLINE}.deletions.sorted.bed >> ${NAMLINE}.only.noMD.biallelic.het.AWDelovl.vcf
NAMLINEHETAWDEL=$(grep -v "^#" ${NAMLINE}.only.noMD.biallelic.het.AWDelovl.vcf | wc -l)
#157,860
echo "Heterozygous biallelic SNPs with AWDelovl" >> ${NAMLINE}.txt
echo ${NAMLINEHETAWDEL} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > ${NAMLINE}.only.noMD.biallelic.het.AWDelovl.genic.vcf
bedtools intersect -a ${NAMLINE}.only.noMD.biallelic.het.AWDelovl.vcf -b ../../annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.genes.gtf >> ${NAMLINE}.only.noMD.biallelic.het.AWDelovl.genic.vcf
NAMLINEHETAWDELGENIC=$(grep -v "^#" ${NAMLINE}.only.noMD.biallelic.het.AWDelovl.genic.vcf | wc -l)
#12,586
echo "Genic Heterozygous biallelic SNPs with AWDelovl" >> ${NAMLINE}.txt
echo ${NAMLINEHETAWDELGENIC} >> ${NAMLINE}.txt
echo "" >> ${NAMLINE}.txt

NAMLINEHETAWDELINTERGENIC=$((NAMLINEHETAWDEL-NAMLINEHETAWDELGENIC))
echo "Intergenic Heterozygous biallelic SNPs with AWDelovl" >> ${NAMLINE}.txt
echo ${NAMLINEHETAWDELINTERGENIC} >> ${NAMLINE}.txt
#intergenic
#157,860-12,586
#145,274
echo "" >> ${NAMLINE}.txt

