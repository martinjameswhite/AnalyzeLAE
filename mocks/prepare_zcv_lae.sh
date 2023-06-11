#!/bin/bash -l
#SBATCH -J PrepZCV
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -o PrepZCV.out
#SBATCH -e PrepZCV.err
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -A m68
#
source activate abacus
#
# This code uses abacusnbody.hod.zcv
#
yaml=./prepare_zcv_lae.yaml
#
# Run once per simulation.
#python -m abacusnbody.hod.zcv.ic_fields --path2config $yaml
# Run once per redshift and cosmology.
#python -m abacusnbody.hod.zcv.zenbu_window --path2config $yaml
#python -m abacusnbody.hod.zcv.zenbu_window --path2config $yaml --want_xi
# Run once per simulation and redshift bin.
#python -m abacusnbody.hod.zcv.advect_fields --path2config $yaml --want_rsd
#python -m abacusnbody.hod.zcv.advect_fields --path2config $yaml --want_rsd --save_3D_power
#
