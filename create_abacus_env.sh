#!/bin/bash
# 
conda create --name abacus --clone lazy-mpi4py
#
source activate abacus
# Install some basic stuff
conda install numpy scipy sympy matplotlib astropy pandas cython -y
conda install -c conda-forge healpy -y
conda install -c conda-forge pyfftw -y
# and CLASS
cd /global/cfs/cdirs/m68/mwhite/class/python/
python setup.py install
cd
# Set up the environment for Jupyter.
conda install ipykernel ipython jupyter # Should already be there.
python3 -m ipykernel install --user --name abacus --display-name Abacus-env
# Install findiff (used for Taylor interpolations)
python3 -m pip install findiff
#
# Install velocileptors. 
python3 -m pip install git+https://github.com/sfschen/velocileptors
#
# Install Abacus
python3 -m pip install git+https://github.com/abacusorg/abacusutils.git
## No longer used: python3 -m pip install abacusutils[zcv]
#
# Currently Abacus needs fast-cksum and it isn't included.
# To install fast-cksum:
# make a directory for fast-cksum.
# git clone https://github.com/abacusorg/fast-cksum.git
# cd fast-cksum; make
# export PYTHONPATH=$PYTHONPATH:path-to/fast-cksum
#
#
#
