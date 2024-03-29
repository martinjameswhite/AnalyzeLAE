#!/bin/bash -l
#SBATCH -J FakeLAE
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH -o FakeLAE.out
#SBATCH -e FakeLAE.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A m68
#
module unload craype-hugepages2M
#
source activate abacus
module load gsl
#
date
srun -n 1 python fake_lae.py
date
#
