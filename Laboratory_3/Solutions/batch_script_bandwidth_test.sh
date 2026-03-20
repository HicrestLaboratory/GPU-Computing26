#!/bin/bash
#SBATCH --partition=edu-short
#SBATCH --account=gpu.computing26
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0
#SBATCH --nodes=1

#SBATCH --job-name=bandwidth_test
#SBATCH --output=outputs/R-%x.%j.out
#SBATCH --error=outputs/R-%x.%j.err

sizes=(1024 4096 16384 65536 262144 4194304 16777216 67108864)

for n in "${sizes[@]}"
do
    ./bin/prefixsum_test $n
done
