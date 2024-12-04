#!/bin/bash
#SBATCH -J snpcallsbylines
#SBATCH -o snpcallsbylines.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e snpcallsbylines.e%A.%a   # error file name
#SBATCH -a 1-25                # 1-36
###SBATCH -n 1                # Total number of mpi tasks requested
###SBATCH --ntasks-per-node 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=8000
#SBATCH -p p1     # queue (partition) -- normal, development, largemem, etc.

#  $SLURM_CPUS_ON_NODE
#  $SLURM_CPUS_PER_TASK
#  $SLURM_ARRAY_TASK_ID

#CORES=${SLURM_CPUS_PER_TASK}
#CORES=24

JOBFILE=joblist.txt

NAMLINE=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)

bash bash_calculatesnpcounts.sh ${NAMLINE}
