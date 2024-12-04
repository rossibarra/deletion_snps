#!/bin/bash

PATH=/gpool/bin/bedtools2/bin:$PATH

MYVCF=$1
MYCVFHEADERS=$2

bedtools subtract -a ${MYVCF} -b B97.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b CML103.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b CML228.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b CML247.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b CML277.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b CML322.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b CML333.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b CML52.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b CML69.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b HP301.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Il14H.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Ki11.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Ki3.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Ky21.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b M162W.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b M37W.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Mo18W.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Ms71.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b NC350.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b NC358.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Oh43.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Oh7B.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b P39.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Tx303.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} temp.2.vcf
cat temp.1.vcf >> temp.2.vcf
bedtools subtract -a temp.2.vcf -b Tzi8.deletions.bed > temp.1.vcf
cp ${MYCVFHEADERS} myvcfnodeletions.vcf
cat temp.1.vcf >> myvcfnodeletions.vcf
rm temp.[12].vcf
