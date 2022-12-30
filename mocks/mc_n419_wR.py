#!/usr/bin/env python3
#
# Runs a Monte-Carlo loop generating mock LAE catalogs and
# computing their clustering.
#
import numpy as np
from   fake_lae import MockLAE
from   calc_wR  import calc_wt





if __name__=="__main__":
    # Set up a random number generator with a
    # fixed seed for reproducability.
    rng    = np.random.default_rng(1)
    # Define the mock catalog, shell and HOD.
    laes   = MockLAE('lae_n419.yaml',3979.,33.)
    params = {'logM_cut':11.75,'logM1':11.75+np.log10(5.),\
              'sigma':0.66,'kappa':0.33,'alpha':0.50}
    laes.set_hod(params)
    laes.generate()
    laes.assign_lum(0.5)
    # Select a field so we have access to the Lside etc.
    diam   = 3.2 * np.pi/180.
    laes.select(diam,[0.,0.,0.])
    chi0   = laes.d['chi0']
    Lside  = laes.d['Lside']
    Lx     = Lside/chi0 * 180./np.pi
    Ly     = Lside/chi0 * 180./np.pi
    ichi   = 1.0/  chi0 * 180./np.pi
    # Now we want to determine the sampling fraction to
    # get the right angular number density.
    ntarget,nbar = 387.4,[]
    for i in range(25):
        # Generate the galaxies.
        offset = rng.uniform(low=-0.5,high=0.5,size=3)
        laes.select(diam,offset)
        nobj   = laes.d['nkeep']
        nbar.append(nobj / (Lx*Ly))
    fsamp = ntarget / np.median(nbar)
    # Generate a uniform random catalog.
    # We offset the RA to eliminate negative RAs
    # just to avoid a warning.
    nran = 100000
    ran  = {}
    ran['RA' ] = rng.uniform(low=-Lx/2.,high=Lx/2.,size=nran) + Lx
    ran['DEC'] = rng.uniform(low=-Ly/2.,high=Ly/2.,size=nran)
    # apply mask.
    rad2 = ( (ran['RA']-Lx)**2 + (ran['DEC'])**2 )*(np.pi/180.)**2
    ran['RA' ] = ran['RA' ][rad2<diam**2/4]
    ran['DEC'] = ran['DEC'][rad2<diam**2/4]
    # Now do the MC loop.
    rval,wts,ngals = None,[],[]
    for i in range(256):
        # Generate the galaxies.
        offset = rng.uniform(low=-0.5,high=0.5,size=3)
        laes.select(diam,offset)
        dat        = {}
        dat['RA' ] = laes.xpos*ichi + Lx
        dat['DEC'] = laes.ypos*ichi
        # downsample
        rand = rng.uniform(low=0,high=1,size=dat['RA'].size)
        ww   = np.nonzero( rand<fsamp )[0]
        dat['RA' ] = dat['RA' ][ww]
        dat['DEC'] = dat['DEC'][ww]
        # apply mask.
        rad2 = ( (dat['RA']-Lx)**2 + (dat['DEC'])**2 )*(np.pi/180.)**2
        dat['RA' ] = dat['RA' ][rad2<diam**2/4]
        dat['DEC'] = dat['DEC'][rad2<diam**2/4]
        # compute the clustering.
        bins,wt = calc_wt(dat,ran)
        rval    = chi0*np.sqrt(bins[:-1]*bins[1:])*np.pi/180.
        wts.append(wt)
        ngals.append(dat['RA'].size)
    wts  = np.array(wts)
    wavg = np.mean(wts,axis=0)
    werr = np.std( wts,axis=0)
    navg = np.mean(np.array(ngals,dtype='float'))
    nerr = np.std( np.array(ngals,dtype='float'))
    # Now write out some results.
    diam *= 180./np.pi
    area  = np.pi * (diam/2)**2
    with open("mc_n419_wR.txt","w") as fout:
        fout.write("# Monte-Carlo calculation of wR using {:d} mocks.\n".\
                   format(wts.shape[0]))
        fout.write("# Field diameter is {:.2f}deg, area {:.2f}deg2.\n".\
                   format(diam,area))
        fout.write("# Sampling by {:.3f} to get {:.1f} LAEs/deg2\n".\
                   format(fsamp,ntarget))
        fout.write("# Have {:.1f}+/-{:.2f} LAEs/field.\n".\
                   format(navg,nerr))
        fout.write("# {:>8s} {:>15s} {:>15s}\n".\
                   format("R[Mpc/h]","wR","dwR"))
        for i in range(rval.size):
            outstr = "{:10.3f} {:15.5e} {:15.5e}".format(rval[i],wavg[i],werr[i])
            fout.write(outstr+"\n")
    #
