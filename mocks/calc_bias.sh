#!/bin/bash -l
#SBATCH -J Bias
#SBATCH -N 1
#SBATCH -t 01:30:00
#SBATCH -o ComputeBias.out
#SBATCH -e ComputeBias.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A m68
#
date
#
source activate abacus
#
export OMP_NUM_THREADS=64
srun -n 1 -c ${OMP_NUM_THREADS} python calc_bias.py
#
date
#
