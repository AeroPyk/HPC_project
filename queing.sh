#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A edu20.dd2356

# The name of the script is ProjectMM
#SBATCH -J ProjectMM

# 5 minutes wall-clock time will be given to this job
#SBATCH -t 00:05:00

# Number of Nodes
#SBATCH --nodes=4

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1

#SBATCH -C Haswell

#SBATCH -e error_file.e
#SBATCH -o output_file.o

export OMP_NUM_THREADS=8

# Run the executable named myexe
# and write the output into my_output_file
export MDIM="2048"
srun ./main > res2048.txt

export MDIM="4096"
srun ./main > res4096.txt