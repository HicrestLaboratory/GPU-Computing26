#!/bin/bash
#SBATCH --partition=edu-short
#SBATCH --account=gpu.computing26
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00

#SBATCH --job-name=test
#SBATCH --output=test-%j.out
#SBATCH --error=test-%j.err

module load CUDA/12.5.0

./PrefixSum 22 0
