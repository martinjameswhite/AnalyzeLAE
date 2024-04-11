#!/bin/bash -l
#SBATCH -J Wt
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH -o Compute_clustering.out
#SBATCH -e Compute_clustering.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A m68
#
date
#
module load python
conda  activate abacus
#
export OMP_NUM_THREADS=64
srun -n 1 -c ${OMP_NUM_THREADS} python compute_clustering.py
#
chmod og+r *
date
#
