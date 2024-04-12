#!/usr/bin/env python3
#
# Runs a Monte-Carlo loop generating mock LAE catalogs and
# computing their clustering.
#
import numpy  as np
import healpy as hp
import time
import sys

sys.path.append('/pscratch/sd/m/mwhite/direct_sht/csht/')
from sht import DirectSHT

from   fake_lae             import MockLAE
from   rotate_to            import rotate_to
from   make_lae_survey_mask import SurveyMask
from   astropy.table        import Table
from   scipy.interpolate    import InterpolatedUnivariateSpline as Spline




if __name__=="__main__":
    # Set up a random number generator with a
    # fixed seed for reproducability.
    rng    = np.random.default_rng(1)
    # Set the name of the field we'll work with and load
    # the mask, random catalog and radial selection fn.
    fname= "cosmos_N419"
    sfn  = 'lae_'+fname+'_sfn.txt'
    mask = SurveyMask('lae_'+fname+'_msk.fits')
    tt   = Table.read('lae_'+fname+'_ran.fits')
    # Set the center of the field.
    cra = np.median(tt['RA'])
    cdc = np.median(tt['DEC'])
    print("Setting field center to ({:f},{:f})".format(cra,cdc),flush=True)
    # Mask the randoms.
    ww,ran     = mask(tt['RA'],tt['DEC']),{}
    ran['RA' ] = tt['RA' ][ww]
    ran['DEC'] = tt['DEC'][ww]
    print("Random size ",len(ran['RA']),flush=True)
    print(time.asctime(),flush=True)
    # Set up the SHT instance, compute rlm and wl.
    sht = DirectSHT(2048,2048,xmax=0.1)
    theta,phi = np.radians(90-ran['DEC']),np.radians(ran['RA'])
    rlm = sht(theta,phi,np.ones(len(theta)))
    rlm/= rlm[0]
    wl  = hp.alm2cl(rlm)
    ells= np.arange(wl.size)
    norm= np.sum(wl*(2*ells+1))/(4*np.pi)
    print("Cl norm {:e}".format(norm),flush=True)
    print(time.asctime(),flush=True)
    # Set up a binning matrix.
    bins= np.zeros( (sht.Nell,sht.Nell) )
    ii,l0,l1 = 0,100,100+200
    while l1<=sht.Nell:
        bins[ii,l0:min(l1,sht.Nell)] = 1/float(l1-l0)
        l0,l1 = l1,l1+200
        ii   += 1
    bins = bins[:ii,:]
    ells = np.dot(bins,np.arange(sht.Nell))/np.sum(bins,axis=1)
    print("ells=",ells,flush=True)
    # Define the mock catalog, shell and HOD.
    laes   = MockLAE('lae_n419.yaml',3941.,sfn)
    params = {'logM_cut':11.10,'logM1':11.80,\
              'sigma':0.66,'kappa':0.33,'alpha':0.66}
    laes.set_hod(params)
    laes.generate()
    laes.assign_lum(0.25)
    # Select a field so we have access to the Lside etc.
    diam   = 4.5 * np.pi/180.
    laes.select(diam,[0.,0.,0.])
    chi0   = laes.d['chi0']
    Lside  = laes.d['Lside']
    Lx     = Lside/chi0 * 180./np.pi
    Ly     = Lside/chi0 * 180./np.pi
    ichi   = 1.0/  chi0 * 180./np.pi
    # Match the min/max values of chi.
    chimin,chimax = np.min(laes.zpos+chi0),np.max(laes.zpos+chi0)
    print("Raw 3D number density ",laes.d['nbar'],flush=True)
    print(time.asctime(),flush=True)
    # Now we want to determine the sampling fraction to
    # get the right angular number density.
    ntarget,nbar = 267.1*(1.-0.0276),[]
    for i in range(25):
        # Generate the galaxies.
        offset = rng.uniform(low=-0.5,high=0.5,size=3)
        laes.select(diam,offset)
        nobj   = laes.d['nkeep']
        nbar.append(nobj / (Lx*Ly))
    fsamp = ntarget / np.median(nbar)
    print("Computed fsamp=",fsamp,flush=True)
    print(time.asctime(),flush=True)
    # Now do the MC loop.
    cls,ngals = [],[]
    for i in range(256):
        # Generate the galaxies.
        offset = rng.uniform(low=-0.5,high=0.5,size=3)
        laes.select(diam,offset)
        dat        = {}
        dat['RA' ] = laes.xpos*ichi
        dat['DEC'] = laes.ypos*ichi
        dat['CHI'] = laes.zpos+chi0
        # Apply radial selection function.
        rand = rng.uniform(low=0,high=1,size=dat['RA'].size)
        ww   = np.nonzero( rand<fsamp )[0]
        dat['RA' ] = dat['RA' ][ww]
        dat['DEC'] = dat['DEC'][ww]
        dat['CHI'] = dat['CHI'][ww]
        # Rotate the objects to the field center and apply mask.
        nra,ndc    = rotate_to(dat['RA'],dat['DEC'],cra,cdc)
        ww         = mask(nra,ndc)
        dat['RA' ] = nra[ww]
        dat['DEC'] = ndc[ww]
        dat['CHI'] = dat['CHI'][ww]
        # Here we would add line-of-sight downsampling if
        # we wanted/needed it to match dN/dz.
        pass
        # compute the clustering.
        theta,phi = np.radians(90-dat['DEC']),np.radians(dat['RA'])
        dlm = sht(theta,phi,np.ones(len(theta)))
        dlm/= dlm[0]
        cl  = hp.alm2cl(dlm-rlm)
        cl  = np.dot(bins,cl)/norm
        cls.append(cl[1:-1])
        ngals.append(dat['RA'].size)
        if i%10==9: print("  ... finished ",i,time.asctime(),flush=True)
    ells = ells[1:-1]
    cls  = np.array(cls)
    cavg = np.mean(cls,axis=0)
    cerr = np.std( cls,axis=0)
    ccor = np.corrcoef(cls,rowvar=False)
    navg = np.mean(np.array(ngals,dtype='float'))
    nerr = np.std( np.array(ngals,dtype='float'))
    # Now write out some results.
    with open("mc_{:s}_cl.txt".format(fname),"w") as fout:
        fout.write("# Monte-Carlo calculation of Cl using {:d} mocks.\n".\
                   format(cls.shape[0]))
        fout.write("# Field "+fname+"\n")
        fout.write("# Centered on ({:.3f},{:.3f})\n".format(cra,cdc))
        fout.write("# Number density is {:.3e}\n".format(laes.d['nbar']))
        fout.write("# Sampling by {:.4f} to get {:.1f} LAEs/deg2\n".\
                   format(fsamp,ntarget))
        fout.write("# Have {:.1f}+/-{:.2f} LAEs/field.\n".\
                   format(navg,nerr))
        fout.write("# chi0={:f}Mpc/h.\n".format(chi0))
        fout.write("# Correlation matrix is:\n")
        for i in range(ccor.shape[0]):
            outstr = "#"
            for j in range(ccor.shape[1]): outstr += " {:8.4f}".format(ccor[i,j])
            fout.write(outstr + "\n")
        fout.write("# {:>8s} {:>15s} {:>15s}\n".\
                   format("ell","Cl","dCl"))
        for i in range(ells.size):
            outstr = "{:10.3f} {:15.5e} {:15.5e}".format(ells[i],cavg[i],cerr[i])
            fout.write(outstr+"\n")
    #
