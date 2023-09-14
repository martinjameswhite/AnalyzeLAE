import numpy as np
import json
import sys
import os

from astropy.table import Table

sys.path.append("../mocks")
from calc_wR import calc_wt


# Put in our fiducial distances and interloper fractions.
chi_dict = {}
chi_dict['N501'] = 4448.0

fint_dict = {}
fint_dict['N501'] = {'p0':0.155,'p1':0.019,'p2':0.019,'p3':0.027}


def compute_clustering(filter_name,field_name,selection):
    """Does the work of computing the clustering."""
    fname       = field_name+"_"+filter_name
    chi0        = chi_dict[filter_name]
    fint        = fint_dict[filter_name]['p'+str(selection)]
    #
    # Read the targets from file.
    dat = Table.read('odin_p{:d}_'.format(selection)+fname+'_dat.fits')
    print("Read information for ",len(dat)," data objects from file.")
    #
    ran = Table.read("lae_"+fname+"_ran.fits")
    print("Read information for ",len(ran)," rand objects from file.")
    #
    # Compute the clustering from the data
    bins,wx = calc_wt(dat,ran)
    rval    = chi0*np.sqrt( bins[:-1]*bins[1:] )*np.pi/180.
    wx     /= (1-fint)**2 # Interloper correction.
    #
    with open('odin_p{:d}_'.format(selection)+fname+'_wx.txt','w') as outfn:
        outfn.write("# Thin-shell angular correlation function.\n")
        outfn.write("# ODIN "+field_name+" "+filter_name+"\n")
        outfn.write("# Selection: "+str(selection)+"\n")
        outfn.write("# Using chi0={:.1f}Mpc/h.\n".format(chi0))
        outfn.write("# Corrected for f_int={:.3f}\n".format(fint))
        outfn.write("# {:>10s} {:>12s}\n".format("r[Mpc/h]","wx(R)"))
        for i in range(rval.size):
            outfn.write("{:12.5f} {:12.4e}\n".format(rval[i],wx[i]))



if __name__=="__main__":
    field = "cosmos"
    for filter in ["N501"]:
        for selection in [0,1,2,3]:
            compute_clustering(filter,field,selection)
