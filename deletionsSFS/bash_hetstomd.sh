#!/bin/bash

#B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.chronly.AWDELs.biallelicindels.hets.sorted.short.bed.gz
BED=$1

zcat ${BED} | sed 's/0\/1/.\/./g' | sed 's/0\/2/.\/./g' | sed 's/1\/2/.\/./g' | sed 's/1\/0/.\/./g' | sed 's/2\/0/.\/./g' | sed 's/2\/1/.\/./g' > $(basename ${BED} .hets.sorted.short.bed.gz).hetstomd.sorted.short.bed
