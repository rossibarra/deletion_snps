#!/bin/bash

#MYVCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-biallelicsnps.noAWDELs.bed
MYVCF=$1
#.bed or .vcf
MYEXT=$2

cat ${MYVCF} | grep "^chr" | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34"\t"$35"\t"$36}' > $(basename ${MYVCF} ${MYEXT}).short.bed
