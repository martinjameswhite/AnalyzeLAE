#!/bin/bash
#
# Ensure that "defaults" in the first channel in the .condarc file.
#
# I needed to use PrgEnv-gnu in order to get the Planck lensing
# likelihood to compile.
# 
# Set up an empty environment for Cobaya (cloning the lazy-mpi4py
# ensures mpi4py is properly included in its NERSC form).
#
conda create --name cobaya --clone lazy-mpi4py
#
# Switch to the environment.
source activate cobaya
#
# Install some basic stuff
conda install numpy scipy sympy matplotlib astropy pandas cython -y
conda install -c conda-forge pyfftw healpy -y
#
# Set up the environment for Jupyter.
conda install ipykernel ipython jupyter -y
python3 -m ipykernel install --user --name cobaya --display-name Cobaya-env
#
# Now install Cobaya
python3 -m pip install cobaya  # --upgrade
#
# and any "cosmo" packages it wants
cobaya-install cosmo -p $SCRATCH/Cobaya/Packages
# you may need to "upgrade" the packages:
#cobaya-install cosmo --upgrade -p $SCRATCH/Cobaya/Packages
#
# Install velocileptors,
#python3 -m pip install git+https://github.com/sfschen/velocileptors
# and Anzu if you want.
#conda install -c conda-forge pyccl chaospy -y
#python3 -m pip install -v git+https://github.com/kokron/anzu
# and findiff for the Taylor series emulators
#python3 -m pip install --upgrade findiff
#
