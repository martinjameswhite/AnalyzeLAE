#!/usr/bin/env python3
#
# Compute the angular correlation function of a
# catalog of objects.
#
import numpy as np
import os
#
from astropy.table import Table
#
from Corrfunc.mocks import DDtheta_mocks as Pairs
from Corrfunc.utils import convert_3d_counts_to_cf



def calc_wt(dat,ran,bins=None):
    """Does the work of calling CorrFunc."""
    # Get number of threads, datas and randoms.
    nthreads = int( os.getenv('OMP_NUM_THREADS','1') )
    Nd,Nr    = len(dat['RA']),len(ran['RA'])
    # RA and DEC should be in degrees, and all arrays should
    # be the same type and it seems as if they need to be
    # 'float' and not e.g. 'float32'.  Ensure this now.
    pp  = 'pair_product'
    dra = dat['RA' ].astype('float')
    ddc = dat['DEC'].astype('float')
    if not 'WT' in dat.keys():
        dwt = np.ones_like(dat['RA']).astype('float')
    else:
        dwt = dat['WT'].astype('float')
    rra = ran['RA' ].astype('float')
    rdc = ran['DEC'].astype('float')
    if not 'WT' in ran.keys():
        rwt = np.ones_like(ran['RA']).astype('float')
    else:
        rwt = ran['WT'].astype('float')
    # Bin edges are specified in degrees, if nothing
    # is passed in, do log-spaced bins.
    if bins is None:
        Nbin = 5
        bins = np.logspace(-1.5,-0.5,Nbin+1)
    # do the pair counting, then convert to w(theta).
    DD = Pairs(1,nthreads,bins,RA1=dra,DEC1=ddc,weights1=dwt,weight_type=pp)
    RR = Pairs(1,nthreads,bins,RA1=rra,DEC1=rdc,weights1=rwt,weight_type=pp)
    DR = Pairs(0,nthreads,bins,RA1=dra,DEC1=ddc,RA2=rra,DEC2=rdc,\
               weights1=dwt,weights2=rwt,weight_type=pp)
    wt = convert_3d_counts_to_cf(Nd,Nd,Nr,Nr,DD,DR,DR,RR)
    # Return the binning and w(theta).
    return( (bins,wt) )
    #


if __name__=="__main__":
    rng  = np.random.default_rng(1)
    # Load the data from file.
    dat  = Table.read('mock_lae_cat.fits')
    chi0 = dat.meta['CHI0']
    Lside= dat.meta['LSIDE']
    Lx   = Lside/chi0 * 180./np.pi
    Ly   = Lside/chi0 * 180./np.pi
    # Generate a uniform random distribution.
    nran = 100000
    ran  = {}
    ran['RA' ] = rng.uniform(low=-Lx/2.,high=Lx/2.,size=nran)
    ran['DEC'] = rng.uniform(low=-Ly/2.,high=Ly/2.,size=nran)
    # Get rid of negative RAs, just for ease.
    dat['RA'] += Lx
    ran['RA'] += Lx
    # and compute w(theta), converting theta to projected
    # distance, R, assuming log-space bins.
    bins,wt = calc_wt(dat,ran)
    rval    = chi0*np.sqrt(bins[:-1]*bins[1:])*np.pi/180.
    # Print the results.
    print("# {:>8s} {:>15s}".format("R[Mpc/h]","w_theta(R)"))
    for i in range(wt.size):
        print("{:10.4f} {:15.5e}".format(rval[i],wt[i]))
    #
