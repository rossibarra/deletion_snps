#!/bin/bash
#SBATCH -J bedtoolsintersectm
#SBATCH -o bedtoolsintersectm.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e bedtoolsintersectm.e%A.%a   # error file name
#SBATCH -a 1-25                # 1-36
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

#bedtools
PATH=/gpool/bin/bedtools2/bin:$PATH

#JOBFILE=matchbedfiles.txt
JOBFILE=nondelbedfiles.txt

VCF=B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz
BED=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)
PREFIX=$(basename ${BED} .bed)

# find intersection
bedtools intersect -a ${VCF} -b ${BED} > ${PREFIX}.intersect.vcf

# find not intersecting: vcf\bed
#bedtools intersect -v -a ${VCF} -b ${BED} > ${PREFIX}.subtract.vcf
