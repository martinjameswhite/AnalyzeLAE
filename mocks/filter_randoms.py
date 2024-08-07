#!/usr/bin/env python
#
#
# Filter the randoms generated by Adam to include the LS/HSC
# coverage information that is used in the N419 selection.
#
import numpy  as np
import sys
import os

from   make_lae_survey_mask import SurveyMask
from   astropy.table        import Table



def filter_randoms(dbase,ranfile,mskfile,outfile):
    """Reads the randoms provided by Adam and filters them
    through a mask."""
    # Read RA/DEC from the file, restrict to those that pass
    # the mask and then write a "filtered" file.
    print("Masking ",ranfile," with ",mskfile,flush=True)
    ran = Table.read(dbase+ranfile)
    ran = ran[ ran['MASKBITS']==0 ]
    ran = ran[ ran['IN_ARJUN_MASK']==False ]
    print("Read ",len(ran)," randoms from ",dbase+ranfile,flush=True)
    mask= SurveyMask(mskfile,None)
    ran = ran[ mask(ran['RA'],ran['DEC']) ]
    print("Kept ",len(ran)," randoms.",flush=True)
    ran.write(outfile,overwrite=True)
    #






if __name__=="__main__":
    dbase   = "/global/cfs/cdirs/desi/users/adamyers/ODIN/"
    ranfile = "randoms-ODIN-cosmos-1.fits"
    for imaging in ["LS","HSC"]:
        mskfile = "/global/cfs/cdirs/desi/users/raichoor/laelbg/odin/phot/"
        mskfile+= "odin-N419-cosmos-angfoot-"+imaging+".fits"
        outfile = "cosmos_"+imaging+"_randoms.fits"
        filter_randoms(dbase,ranfile,mskfile,outfile)
        #
