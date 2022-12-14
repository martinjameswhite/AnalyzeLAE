#!/bin/bash -l
#SBATCH -J Prepare
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -o Prepare.out
#SBATCH -e Prepare.err
#SBATCH -q debug
#SBATCH -C haswell
#SBATCH -A m68
#
module unload craype-hugepages2M
#
source activate abacus
#
# This code uses abacusnbody.hod.prepare_sim, but this itself has
# subsampling hard-coded that is not appropriate for LBGs or LAEs.
# You should modify the prepare_sim code so that the "MT" else
# condition in subsample_halos has e.g.:
#        downfactors[x > 10.5] = 1
#
python -m abacusnbody.hod.prepare_sim \
  --path2config ./lae_base.yaml --alt_simname AbacusSummit_high_c000_ph100
#
