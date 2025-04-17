#!/bin/bash

# Slurm job options
#SBATCH --account=d448
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --job-name=Channel
#SBATCH --time=04:00:00
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

# Due to LAMMPS version compiled with cMake
module load PrgEnv-gnu

# Now load the "lammps" package
export PATH="$PATH:/work/d448/d448/yichongchen/apprepo/mylammps/build/"

# Set the number of threads to 1
export OMP_NUM_THREADS=1

# Run the parallel job (input file should be in cwd)
srun -n 512 lmp < in.equil > log_equil
srun -n 512 lmp < in.meas  > log_meas
#    Using 512 MPI processes 
#    128 MPI processes per node
