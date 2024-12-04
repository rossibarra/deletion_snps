#!/bin/bash
#SBATCH -J jeffpipelinepreprocess
#SBATCH -o jeffpipelinepreprocess.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e jeffpipelinepreprocess.e%A.%a   # error file name
###SBATCH -a 25 #1-25                # 1-36
###SBATCH -n 1                # Total number of mpi tasks requested
###SBATCH --ntasks-per-node 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=2500
#SBATCH -p p1,gcluster     # queue (partition) -- normal, development, largemem, etc.

#  $SLURM_CPUS_ON_NODE
#  $SLURM_CPUS_PER_TASK
#  $SLURM_ARRAY_TASK_ID

#CORES=${SLURM_CPUS_PER_TASK}
#CORES=24

MYVCF=B73v5.NAM.chronly.vcf.gz
NAMLINE=Tzi8

PATH=/gpool/bin/bcftools-1.10.2/bin:$PATH
#PATH=/gpool/bin/bedtools2/bin:$PATH

bcftools view -M2 -v snps ${MYVCF} > NAM.biallelic.snpsonly.vcf
cat ../B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.headers > NAM.biallelic.indelsonly.vcf
bcftools view -M3 ${MYVCF} | grep "*" >> NAM.biallelic.indelsonly.vcf

bcftools concat NAM.biallelic.snpsonly.vcf NAM.biallelic.indelsonly.vcf > NAM.goodSNPs.vcf
# Run one or the other. Latter is preferred
#bcftools filter -S . -e 'GT=="0/1" || GT=="0/2" || GT=="1/2"' NAM.goodSNPs.vcf > NAM.goodSNPs_hetstomissing.vcf
bcftools filter -S . -e 'GT=="het"' NAM.goodSNPs.vcf > NAM.goodSNPs_hetstomissing.vcf
