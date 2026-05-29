#!/bin/bash
#SBATCH --partition=edu-short
#SBATCH --account=gpu.computing26
#SBATCH --job-name=bcastcollective
#SBATCH --output=bcastcollective.o
#SBATCH --error=bcastcollective.e
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --mem=1G

# Load modules
module load OpenMPI
module load CUDA/12.5.0

# Print compiler versions
gcc --version
mpicc --version
nvcc --version

# Compile
mpicc broadcast_collective_staging.c -o bcastcollective -lcudart

mpirun -np 4 bcastcollective 

