#!/usr/bin/env python
#
#
# Generate a survey mask (an array of "filled" Healpix pixels)
# from a large (pre-generated) random file.  Write the result
# to a FITS file.
#
# Also contains a class that uses such a mask to return True
# for objects within the survey footprint/mask.
#
# Extinction correction coefficients for the ODIN NB filters
# a419 = 4.3238  # extinction correction for E(B-V)=1
# a501 = 3.54013 # extinction correction for E(B-V)=1
# a673 = 2.43846 # extinction correction for E(B-V)=1
# Magnitudes  ; reddening corrected mags
# m419=22.5-2.5*np.log10(od['flux_n419'])-a419*od['ebv']  # + 22.1069
# m501=22.5-2.5*np.log10(od['flux_n501'])-a501*od['ebv']  # + 22.1069
# m673=22.5-2.5*np.log10(od['flux_n673'])-a673*od['ebv']  # + 22.6585
#
#
import numpy  as np
import healpy as hp
import sys
import os


from   astropy.table import Table



def make_survey_mask(ran_fname,filter_name,nside=8192,is_nest=True):
    """Generate an astropy Table of the survey (inclusion)
    mask from the coordinates read from a file of randoms.
    We assume the file is large enough that Poisson fluctuations
    can be entirely neglected."""
    # Read RA/DEC from the file, restrict to those with observations
    # in the desired filter and generate a list of Healpix pixel numbers.
    ran   = Table.read(ran_fname)
    print(ran.keys())
    ran   = ran[ ran['MASKBITS']==0 ]
    ran   = ran[ ran['IN_ARJUN_MASK']==False ]
    ran   = ran[ ran['NOBS_'+filter_name]>10 ]
    theta = np.radians(90.-ran['DEC'])
    phi   = np.radians(ran['RA'])
    pixs  = hp.ang2pix(nside,theta,phi,nest=is_nest)
    # and generate a table containing the information to be returned.
    tt = Table({"HPXPIXEL":np.unique(pixs)})
    tt.meta["HPXNSID"] = nside
    tt.meta["HPXNEST"] = is_nest
    tt.meta["COMMENT"] = r'Mask lists HPX pixels included in survey'
    # We also want a random file -- this could be downsampled.
    rr = Table({'RA':ran['RA'][::15],'DEC':ran['DEC'][::15]})
    return( (tt,rr) )
    #



class SurveyMask:
    def __init__(self,maskfn):
        """Reads a FITS file containing an inclusion mask."""
        mskd       = Table.read(maskfn)
        self.nside = mskd.meta["HPXNSID"]
        self.nest  = mskd.meta["HPXNEST"]
        self.pixs  = mskd["HPXPIXEL"]
    def __call__(self,ras,decs):
        """Returns a boolean array of whether the points pass the mask,
        with the points given by (RA,DEC) in degrees."""
        tt   = np.radians(90.-decs)
        pp   = np.radians(ras)
        pixs = hp.ang2pix(self.nside,tt,pp,nest=self.nest)
        return(np.in1d(pixs,self.pixs))
        #




if __name__=="__main__":
    if len(sys.argv)!=3:
        raise RuntimeError("Usage: "+sys.argv[0]+" <ran-fname> <filter_name>")
    ran_fname   = sys.argv[1]
    filter_name = sys.argv[2].upper()
    #
    tt,rr = make_survey_mask(ran_fname,filter_name)
    #
    tt.write("lae_survey_msk_{:s}.fits".format(filter_name),overwrite=True)
    rr.write("lae_survey_ran_{:s}.fits".format(filter_name),overwrite=True)
    #
