#!/bin/bash -l
#SBATCH -J PrepZCV
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -o PrepZCV.out
#SBATCH -e PrepZCV.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A m68
#
source activate abacus
#
# This code uses abacusnbody.hod.zcv
# Run once per simulation.
#python -m abacusnbody.hod.zcv.ic_fields --path2config ./prepare_zcv_sml.yaml 
# Run once per redshift and cosmology.
#python -m abacusnbody.hod.zcv.zenbu_window --path2config ./prepare_zcv_sml.yaml
# Run once per simulation and redshift bin.
#python -m abacusnbody.hod.zcv.advect_fields --path2config ./prepare_zcv_sml.yaml --want_rsd
# Run once per simulation and redshift bin (needed even if want_rsd is true).
#python -m abacusnbody.hod.zcv.advect_fields --path2config ./prepare_zcv_sml.yaml
#
