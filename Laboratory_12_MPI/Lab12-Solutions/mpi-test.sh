#!/bin/bash
#SBATCH --partition=edu-short
#SBATCH --account=gpu.computing26
#SBATCH --job-name=mpi
#SBATCH --output=mpi.o
#SBATCH --error=mpi.e
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0
#SBATCH --mem=1G

# Load modules
module load OpenMPI

# Print compiler versions
gcc --version
mpicc --version


mpirun -np 4 hostname 

