#!/bin/bash -l
#SBATCH -J LAE_Cl
#SBATCH -t 0:25:00
#SBATCH -N 1
#SBATCH -o LAE_Cl.out
#SBATCH -e LAE_Cl.err
#SBATCH -p debug
#SBATCH -C haswell
#SBATCH -A m68
#
date
source activate cmb
#
# Now run the code.
export OMP_NUM_THREADS=16
srun -n 1 -c 16 python calc_lae_cl.py
#
date
#
