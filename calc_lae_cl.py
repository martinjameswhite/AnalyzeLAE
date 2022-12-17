#!/usr/bin/env python3
#
# Code to compute the angular power spectra and covariance
# matrices given FlatSky maps.  This is essentially a thin
# wrapper around NaMaster.
#
# THIS REQUIRES SOME THINKING IF WE'RE GOING TO BE
# USING MOCK-BASED COVARIANCES.
#
import numpy    as np
import pymaster as nmt
import json
import sys
from astropy.table import Table






def read_maps(isamp=1):
    """Returns the galaxy maps and masks."""
    # Read the galaxy density and mask.
    pref = "lae_s{:02d}".format(isamp)
    tt   = Table.read(pref+'.fits')
    Lx,Ly= tt.meta['LX'],tt.meta['LY'] # Radians.
    gals = tt['del']
    gmsk = tt['msk']
    # and return the results.
    return( (Lx,Ly,gals,gmsk) )
    #




def read_thy_cl(Lmax,isamp):
    """Read best-fit theory (incl. noise) and interpolate to full length."""
    mod = np.loadtxt("lae_s{:02d}_mod.txt".format(isamp))
    # Probably need to extrapolate this to higher
    # and lower ell.  Set high ell power to zero.
    ells = np.arange(Lmax)
    cgg  = np.interp(ells,mod[:,0],mod[:,1],right=0)
    return(cgg)
    #




def make_bins(lmin,lmax,LperBin):
    """Sets up the ell bins.
    For flat-sky fields, bandpowers are simply defined as intervals in ell, and
    NaMaster doesn't currently support any weighting scheme within each interval."""
    llo = np.arange(lmin,lmax-LperBin,LperBin)
    lhi = llo + LperBin
    # And tell NaMaster to set it up.
    bins = nmt.NmtBinFlat(llo,lhi)
    return(bins)
    #



def calc_pseudo_cl(Lmin=250,Lmax=2000,LperBin=200,isamp=1):
    """Compute the pseudo-Cl and coupling matrix."""
    # Set up an empty "output" dictionary.
    outd = {}
    # Load the maps.
    Lx,Ly,gals,gmsk = read_maps(isamp)
    # Work out bins, extending the high-ell range.
    bins   = make_bins(Lmin,2*Lmax,LperBin)
    ell    = bins.get_effective_ells()
    bmax   = np.argmax(np.nonzero(ell<Lmax)[0]) + 1
    print("Cutting at index ",bmax," with ell[bmax]=",ell[bmax])
    outd['ell'] = ell[:bmax].tolist()
    #
    # Note NaMaster takes care of multiplying our maps by masks.
    #
    # We construct the workspace specifically in order to
    # have access to the bandpower windows.
    #
    galxy  = nmt.NmtFieldFlat(Lx,Ly,gmsk,[gals])
    wsp    = nmt.NmtWorkspaceFlat()
    wsp.compute_coupling_matrix(galxy,galxy,bins)
    # In the full-sky routines we could now call:
    # wla    = wsp.get_bandpower_windows()[0,:bmax,0,:] # Just want W^{00}.
    # however there doesn't seem to be a flat-sky equivalent so we're going
    # to approximate it.  Probably should file a pull request on this.
    wla = np.zeros( (len(outd['ell']),Lmax) )
    llo = np.arange(Lmin,Lmax-LperBin,LperBin)
    for i in range(llo.size):
        for j in range(llo[i],llo[i]+LperBin): wla[i,j] = 1.0/LperBin
    #
    outd['wla'] = wla.tolist()
    outd['wla_comment'] = "C_bin = sum_{ell=0}^{Nell-1} W_bin,ell C_ell"
    # Now compute the galaxy autospectrum.
    gauto  = nmt.compute_full_master_flat(galxy,galxy,bins,workspace=wsp)
    cgg    = gauto[0][:bmax]
    outd['cgg'] = cgg.tolist()
    #
    # Correct the C_l for the pixel window function?
    #
    pass
    #
    # Compute an approximation to the Gaussian covariance matrix
    # using  https://arxiv.org/abs/astro-ph/0105302 , S 3.
    # Work out the powers of the window function, w_i.
    thresh = 0.01   # An arbitrary, small positive number.
    tmsk   = gmsk[gmsk>thresh]
    w2,w4  = np.mean(tmsk**2),np.mean(tmsk**4)
    fsky_g = np.sum(tmsk)/gmsk.size*(Lx*Ly/4/np.pi) * w2**2/w4
    print("Galaxy w2={:7.4f}, w4={:7.4f}, fsky={:10.6f}".format(w2,w4,fsky_g))
    # Need theory input spectra of full length.  These
    # spectra contain the stochastic noise terms.
    cgg = read_thy_cl(Lmax,isamp)
    # Now we have the Gaussian, diagonal covariances, per ell.
    lfact  = 1.0/(2*np.arange(Lmax)+1)
    var_aa = np.diag(2*cgg**2*lfact/fsky_g)
    # bin these using the bandpower weights.
    covar_aa = np.dot(wla,np.dot(var_aa,wla.T))
    outd['cov'] = covar_aa.tolist()
    # Pack the full Cov.
    return(outd)
    #











if __name__=="__main__":
    if len(sys.argv)==1:
        isamp = 1   # Default.
    elif len(sys.argv)==2:
        isamp = int(sys.argv[1])
    else:
        raise RuntimeError("Usage: "+sys.argv[0]+" [isamp]")
    Lmin,Lmax,LperBin = 200,2000,200
    outd = calc_pseudo_cl(Lmin,Lmax,LperBin,isamp)
    #
    with open("lae_s{:02d}.json".format(isamp),"w") as fout:
        json.dump(outd,fout,indent=2)
