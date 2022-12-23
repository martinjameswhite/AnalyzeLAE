#!/bin/bash -l
#SBATCH -J LAE_map
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -o LAEmap.out
#SBATCH -e LAEmap.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A m68
#
module unload craype-hugepages2M
#
source activate abacus
export PYTHONPATH=${PYTHONPATH}:${PWD}/../Cobaya/lss_likelihood
#
date
srun -n 1 python fake_map.py
date
#
