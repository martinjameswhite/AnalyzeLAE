#!/usr/bin/env python3
#
# Use the pre-computed data and random files to generate
# alms for the data and randoms.
# Uses the "directsht" library.
#
import numpy  as np
import healpy as hp
import json
import glob
import sys
from   astropy.table import Table


sys.path.append('/pscratch/sd/m/mwhite/direct_sht/csht/')
from sht import DirectSHT








def calc_alms(isamp,ifield,ifilt):
    # Write a log file, to get around annoying buffering issues at NERSC.
    fb   = "s{:d}_{:s}_{:s}".format(isamp,ifield,ifilt.upper())
    flog = open("calc_alms_"+fb+".log","w")
    flog.write("Running "+sys.argv[0]+"\n")
    flog.write("Generating alms for sample "+\
               str(isamp)+" of "+ifield+" "+ifilt+"\n")
    flog.flush()
    #
    sht = DirectSHT(2048,4096,xmax=0.1)
    flog.write("Using Nl={:d}, Nx={:d}, xmax={:f}.\n".\
               format(sht.Nell,sht.Nx,sht.xmax))
    flog.flush()
    # Start with the data.
    db = "../odin/"
    fn = "odin_"+fb+"_dat.fits"
    tt = Table.read(db+fn)
    nd = len(tt)
    flog.write("Data file name: "+db+fn+"\n")
    flog.write("Read data for {:d} objects.\n".format(nd))
    flog.flush()
    # generate the partial map:
    if len(tt)>0:
        print("Min DEC is ",np.min(tt['DEC']))
        print("Max DEC is ",np.max(tt['DEC']))
        theta,phi = np.radians(90-tt['DEC']),np.radians(tt['RA'])
        if 'WT' in tt.keys():
            wt = tt['WT']
        else:
            wt = np.ones_like(tt['RA'])
        dalm = sht(theta,phi,wt)
    # Print some summary statistics.
    flog.write("Done with galaxy file.\n\n")
    flog.flush()
    #
    # Now the randoms.
    #
    db = '../odin/'
    fn = "lae_{:s}_{:s}_ran.fits".format(ifield,ifilt.upper())
    tt = Table.read(db+fn)
    nr = len(tt)
    flog.write("Rand file name: "+db+fn+"\n")
    flog.write("Read data for {:d} objects.\n".format(nr))
    flog.flush()
    if len(tt)>0:
        theta,phi = np.radians(90-tt['DEC']),np.radians(tt['RA'])
        if 'WT' in tt.keys():
            wt = tt['WT']
        else:
            wt = np.ones_like(tt['RA'])
        ralm = sht(theta,phi,wt)
    flog.write("Done with random file.\n")
    flog.close()
    #
    # Renormalize the alms of the data and randoms such
    # that a_{00}=1.
    dalm /= dalm[0]
    ralm /= ralm[0]
    #
    outdict = {}
    outdict['isamp' ]=isamp
    outdict['ifield']=ifield
    outdict['ifilt' ]=ifilt.upper()
    outdict['Ndata' ]=nd
    outdict['Nrand' ]=nr
    outdict['Nl'    ]=sht.Nell
    outdict['Nx'    ]=sht.Nx
    outdict['Nlm'   ]=sht.Nlm
    outdict['xmax'  ]=sht.xmax
    outdict['dlm_re']=dalm.real.tolist()
    outdict['dlm_im']=dalm.imag.tolist()
    outdict['rlm_re']=ralm.real.tolist()
    outdict['rlm_im']=ralm.imag.tolist()
    #
    outfn = "odin_"+fb+"_alm.json"
    with open(outfn,"w") as fout:
        json.dump(outdict,fout,indent=2)
    #




            
if __name__=="__main__":
    if len(sys.argv)==4:
        isamp = int(sys.argv[1])
        ifield= sys.argv[2]
        ifilt = sys.argv[3]
    else:
        raise RuntimeError("Usage: "+sys.argv[0]+" <isamp> <field> <filter>")
    calc_alms(isamp,ifield,ifilt)
