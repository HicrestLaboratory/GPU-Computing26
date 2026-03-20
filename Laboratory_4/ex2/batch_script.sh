#!/bin/bash
#SBATCH --partition=edu-short
#SBATCH --account=gpu.computing26
#SBATCH --gres=gpu:0
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

#SBATCH --job-name=gemm_flops
#SBATCH --output=outputs/R-%x.%j.out
#SBATCH --error=outputs/R-%x.%j.err

./bin/gemm $1 
