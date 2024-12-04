#!/bin/bash

PATH=/gpool/bin/bedtools2/bin:$PATH

MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-biallelicsnps.vcf.gz

#cat ../B97.deletions.merged.bed ../CML103.deletions.merged.bed ../CML228.deletions.merged.bed ../CML247.deletions.merged.bed ../CML277.deletions.merged.bed ../CML322.deletions.merged.bed ../CML333.deletions.merged.bed ../CML52.deletions.merged.bed ../CML69.deletions.merged.bed ../HP301.deletions.merged.bed ../IL14H.deletions.merged.bed ../Ki11.deletions.merged.bed ../Ki3.deletions.merged.bed ../Ky21.deletions.merged.bed ../M162W.deletions.merged.bed ../Mo18W.deletions.merged.bed ../MS37W.deletions.merged.bed ../Ms71.deletions.merged.bed ../NC350.deletions.merged.bed ../NC358.deletions.merged.bed ../OH43.deletions.merged.bed ../OH7B.deletions.merged.bed ../P39.deletions.merged.bed ../TX303.deletions.merged.bed ../TZi8.deletions.merged.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin > NAM.deletions.merged.bed

bedtools intersect -v -a ${MYVCF} -b NAM.deletions.merged.bed > $(basename ${MYVCF} .vcf.gz).noAWDELs.bed
