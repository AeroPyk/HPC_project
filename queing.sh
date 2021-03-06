#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A edu20.dd2356

# The name of the script is ProjectMM
#SBATCH -J ProjectMM

# 5 minutes wall-clock time will be given to this job
#SBATCH -t 01:00:00

# Number of Nodes
#SBATCH --nodes=144

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1

#SBATCH -C Haswell

#SBATCH -e error_file.e
#SBATCH -o output_file.o

mkdir "res"
chmod +x BeskowTools/*

for dim in 1200 2400
do
  MDIM=$dim
  export MDIM
  echo "Mat=$MDIM*$MDIM"

  for TH in 1 4 8 16 24 32
  do
    OMP_NUM_THREADS=$TH
    export OMP_NUM_THREADS

    srun ./main >> "res/res_mat_$dim.txt" # We get the results per matrices
  done
done

# We can't naive multiply with 1 thread and size 4096, it takes ages
for dim in 3600 4200
do
  MDIM=$dim
  export MDIM
  echo "Mat=$MDIM*$MDIM"

  for TH in 8 16 24 32
  do
    OMP_NUM_THREADS=$TH
    export OMP_NUM_THREADS

    srun ./main >> "res/res_mat_$dim.txt" # We get the results per matrices
  done
done

./BeskowTools/syncB.sh