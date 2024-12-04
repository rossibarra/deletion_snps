#!/bin/bash
#SBATCH -J findambsnps
#SBATCH -o findambsnps.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e findambsnps.e%A.%a   # error file name
#SBATCH -a 25   #1-25                # 1-36
###SBATCH -n 1                # Total number of mpi tasks requested
###SBATCH --ntasks-per-node 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=30000
#SBATCH -p p1     # queue (partition) -- normal, development, largemem, etc.

#  $SLURM_CPUS_ON_NODE
#  $SLURM_CPUS_PER_TASK
#  $SLURM_ARRAY_TASK_ID

#CORES=${SLURM_CPUS_PER_TASK}
#CORES=24

JOBFILE=joblist.txt
TEANNOT=Zm-B73-REFERENCE-NAM-5.0.TE.gff3

NAMLINE=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)

PATH=/gpool/bin/bcftools-1.10.2/bin:$PATH
PATH=/gpool/bin/bedtools2/bin:$PATH

#BEGIN NAMLINE DEPENDENT STUFF
###bcftools view -s ${NAMLINE} NAM.goodSNPs_hetstomissing.vcf > ${NAMLINE}.vcf

bedtools intersect -a ${NAMLINE}.deletions.merged.bed -b ${NAMLINE}.vcf -wao > ${NAMLINE}.deletions.snps

grep -v "0$" ${NAMLINE}.deletions.snps | grep -v "\./\." > ${NAMLINE}.overlap.deletions.snps

bedtools intersect -v -a ${NAMLINE}.overlap.deletions.snps -b ${TEANNOT} > ${NAMLINE}.overlap.deletions.noTEovl.snps

grep 0/0 ${NAMLINE}.overlap.deletions.noTEovl.snps > ${NAMLINE}.noTEovl.bad.deletions
grep "\,\*" ${NAMLINE}.overlap.deletions.noTEovl.snps | grep 1/1 >> ${NAMLINE}.noTEovl.bad.deletions
grep "*," ${NAMLINE}.overlap.deletions.noTEovl.snps| grep "2/2" >> ${NAMLINE}.noTEovl.bad.deletions

#BY XSOME
#bedtools intersect -a ${NAMLINE}.deletions.merged.${XSOME}.bed -b ${NAMLINE}.vcf -wao > ${NAMLINE}.${XSOME}.deletions.snps

#grep -v "0$" ${NAMLINE}.${XSOME}.deletions.snps | grep -v "\./\." > ${NAMLINE}.${XSOME}.overlap.deletions.snps

#grep 0/0 ${NAMLINE}.${XSOME}.overlap.deletions.snps > ${NAMLINE}.${XSOME}.bad.deletions
#grep "\,\*" ${NAMLINE}.${XSOME}.overlap.deletions.snps | grep 1/1 >> ${NAMLINE}.${XSOME}.bad.deletions
#grep "*," ${NAMLINE}.${XSOME}.overlap.deletions.snps| grep "2/2" >> ${NAMLINE}.${XSOME}.bad.deletions
