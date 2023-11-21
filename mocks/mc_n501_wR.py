#!/usr/bin/env python3
#
# Runs a Monte-Carlo loop generating mock LAE catalogs and
# computing their clustering.
#
import numpy as np

from   fake_lae             import MockLAE
from   calc_wR              import calc_wt
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
    fname= "cosmos_N501"
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
    # Define the mock catalog, shell and HOD.
    laes   = MockLAE('lae_n501.yaml',4448.,sfn)
    params = {'logM_cut':10.80,'logM1':11.50,\
              'sigma':0.50,'kappa':1.00,'alpha':0.50}
    params = {'logM_cut':11.00,'logM1':12.30,\
              'sigma':0.66,'kappa':1.00,'alpha':0.33}
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
    # Now we want to determine the sampling fraction to
    # get the right angular number density.
    ntarget,nbar = 209.3,[]
    for i in range(25):
        # Generate the galaxies.
        offset = rng.uniform(low=-0.5,high=0.5,size=3)
        laes.select(diam,offset)
        nobj   = laes.d['nkeep']
        nbar.append(nobj / (Lx*Ly))
    fsamp = ntarget / np.median(nbar)
    # Now do the MC loop.
    rval,wxs,ngals = None,[],[]
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
        bins,wx = calc_wt(dat,ran)
        rval    = chi0*np.sqrt( bins[:-1]*bins[1:] )*np.pi/180.
        wxs.append(wx)
        ngals.append(dat['RA'].size)
    wxs  = np.array(wxs)
    wavg = np.mean(wxs,axis=0)
    werr = np.std( wxs,axis=0)
    wcor = np.corrcoef(wxs,rowvar=False)
    navg = np.mean(np.array(ngals,dtype='float'))
    nerr = np.std( np.array(ngals,dtype='float'))
    # Now write out some results.
    with open("mc_{:s}_wR.txt".format(fname),"w") as fout:
        fout.write("# Monte-Carlo calculation of wR using {:d} mocks.\n".\
                   format(wxs.shape[0]))
        fout.write("# Field "+fname+"\n")
        fout.write("# Centered on ({:.3f},{:.3f})\n".format(cra,cdc))
        fout.write("# Number density is {:.3e}\n".format(laes.d['nbar']))
        fout.write("# Sampling by {:.4f} to get {:.1f} LAEs/deg2\n".\
                   format(fsamp,ntarget))
        fout.write("# Have {:.1f}+/-{:.2f} LAEs/field.\n".\
                   format(navg,nerr))
        fout.write("# chi0={:f}Mpc/h.\n".format(chi0))
        fout.write("# Correlation matrix is:\n")
        for i in range(rval.size):
            outstr = "#"
            for j in range(rval.size): outstr += " {:8.4f}".format(wcor[i,j])
            fout.write(outstr + "\n")
        fout.write("# {:>8s} {:>15s} {:>15s}\n".\
                   format("R[Mpc/h]","wR","dwR"))
        for i in range(rval.size):
            outstr = "{:10.3f} {:15.5e} {:15.5e}".format(rval[i],wavg[i],werr[i])
            fout.write(outstr+"\n")
    #
