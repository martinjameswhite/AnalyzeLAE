#!/bin/bash -l
#SBATCH -J Prepare
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -o Prepare.out
#SBATCH -e Prepare.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A m68
#
conda activate abacus
#
# This code uses abacusnbody.hod.prepare_sim, but this itself has
# subsampling hard-coded that is not appropriate for LBGs or LAEs.
# You should modify the prepare_sim code so that the "MT" else
# condition in subsample_halos has e.g.:
#        downfactors[x > 10.5] = 1
# You may also want/need to modify the cutoff in the
# "submask_particles" subroutine.
#
python -m abacusnbody.hod.prepare_sim \
  --path2config ./lae_base.yaml --alt_simname AbacusSummit_high_c000_ph100
#
