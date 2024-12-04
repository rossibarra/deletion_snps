#!/bin/bash

cat B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.hets.header.chr.vcf
grep -P "0/1" B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.header.chr.bed >> B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.hets.header.chr.vcf
grep -P "1/0" B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.header.chr.bed >> B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.hets.header.chr.vcf
grep -P "0/2" B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.header.chr.bed >> B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.hets.header.chr.vcf
grep -P "2/0" B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.header.chr.bed >> B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.hets.header.chr.vcf
grep -P "1/2" B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.header.chr.bed >> B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.hets.header.chr.vcf
grep -P "2/1" B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.header.chr.bed >> B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.noAWDELs.biallelicindels.hets.header.chr.vcf
