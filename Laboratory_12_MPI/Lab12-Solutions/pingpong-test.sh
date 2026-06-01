#!/bin/bash
#SBATCH --partition=edu-short
#SBATCH --account=gpu.computing26
#SBATCH --job-name=pingpong
#SBATCH --output=pingpong.o
#SBATCH --error=pingpong.e
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0
#SBATCH --mem=1G

# Load modules
module load OpenMPI
module load CUDA/12.5.0

# Print compiler versions
gcc --version
mpicc --version

# Compile
mpicc pingpongtest.c -o pingpong -lcudart

mpirun -np 4 pingpong 

