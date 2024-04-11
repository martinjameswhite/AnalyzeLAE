#!/bin/bash -l
#SBATCH -J LAE_MASK
#SBATCH -N 1
#SBATCH -t 00:05:00
#SBATCH -o LAE_MASK.out
#SBATCH -e LAE_MASK.err
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -A m68
#
module load python
conda activate cmb
#
db=/global/cfs/cdirs/desi/users/adamyers/ODIN
# For the N501 sample:
fn=randoms-ODIN-2band-1.fits
# For the N419 sample:
#fn=randoms-ODIN-cosmos-1.fits
db=.
fn=cosmos_HSC_randoms.fits
#
filter=N501
filter=N419
#
export OMP_NUM_THREADS=32
#
date
python make_lae_survey_mask.py ${db}/${fn} ${filter}
date
#
for suf in ran msk ; do
  mv lae_survey_${filter}_${suf}.fits lae_cosmos_${filter}_${suf}.fits
done
#
chmod -R og+r *
#
