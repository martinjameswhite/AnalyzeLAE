#!/bin/bash
# 
conda create --name cmb --clone lazy-mpi4py
#
source activate cmb
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
python3 -m ipykernel install --user --name cmb --display-name CMB-env
# Install findiff (used for Taylor interpolations)
python3 -m pip install --upgrade findiff
#
# Install velocileptors. 
python3 -m pip install git+https://github.com/sfschen/velocileptors
# and Anzu
conda install -c conda-forge chaospy -y
python3 -m pip install -v git+https://github.com/kokron/anzu
# and corner for visualizing the MCMC chains:
python3 -m pip install corner
#
#python3 -m pip install cobaya --upgrade
#
# Then need to install NaMaster, using "install_nmt.sh" script
# or directly from conda-forge.
conda install -c conda-forge namaster -y
#
