#!/usr/bin/env python3
#
import numpy as np
import yaml
import asdf
#
from astropy.table              import Table
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta




# These are the keys we will copy from the simulation meta-data.
simkeys = ['n_s', 'omega_b', 'omega_cdm', 'omega_ncdm', 'N_ncdm', 'N_ur',\
           'H0', 'w0', 'wa', 'w', \
           'Omega_DE', 'Omega_K', 'Omega_M', 'Omega_Smooth', \
           'Redshift', 'ScaleFactor', \
           'OmegaNow_DE', 'OmegaNow_K', 'OmegaNow_m', \
           'f_growth', 'fsmooth', 'Growth', 'Growth_on_a_n', \
           'SimComment', 'SimName', 'SimSet', \
           'BoxSize', 'NP', 'BoxSizeHMpc', 'HubbleTimeHGyr', \
           'ParticleMassHMsun']




def get_metadata(simdir,simname,redshift,proxy):
    """This routine gets the metadata for the simulation.  For simulations
    in the Abacus database this duplicates the get_meta method, but that
    database is not updated frequently and for some simulations some keys
    are missing.
    Proxy is the name of a simulation with the same cosmology that has all
    of the keys, used as a work-around for missing data."""
    # First fill in the "missing" information from the proxy simulation.
    dd = {}
    pr = get_meta(proxy,redshift)
    for k in ['n_s','omega_b','omega_cdm','omega_ncdm',\
              'N_ncdm','N_ur','SimSet']:
        dd[k] = pr[k]
    # Now get the remaining information from the ASDF files directly.
    fn = simdir+simname+\
         '/halos/z{:.3f}/halo_info/halo_info_000.asdf'.format(float(redshift))
    with asdf.open(fn) as af:
        meta = af['header']
        for k in simkeys:
            if k in meta: dd[k] = meta[k]
    return(dd)
    #


class MockLAE:
    """A class to handle making mock LAE samples from halo catalogs."""
    # We may want to do a better job of holding the information in
    # this class.  Right now it's not very well thought through what is
    # kept internally vs. what is passed to each method.
    def __init__(self,yaml_file,chi0,sfn_fname):
        """Set up the class."""
        # Load the selection function; 3 columns: z, chi, p(chi).
        self.sfn = np.loadtxt(sfn_fname)
        # Load the config file and parse in relevant parameters
        config          = yaml.safe_load(open(yaml_file))
        self.sim_params = config['sim_params']
        self.HOD_params = config['HOD_params']
        self.clustering_params = config['clustering_params']
        #
        # Get the metaparameters for the simulation -- this is a
        # workaround for the incomplete database.  It is hardcoded
        # assuming the cosmology is the same as for AbacusSummit base.
        self.meta = get_metadata(self.sim_params['sim_dir'],\
                                 self.sim_params['sim_name'],\
                                 self.sim_params['z_mock'],\
                                 'AbacusSummit_base_c000_ph000')
        #self.meta = get_meta(self.sim_params['sim_name'],\
        #                     redshift=self.sim_params['z_mock'])
        #
        # additional parameter choices -- we will apply RSD
        # ourselves so override the setting in HOD_params.
        self.want_rsd      = False
        self.write_to_disk = self.HOD_params['write_to_disk']
        #
        # Save some configuration parameters for later use.
        self.d = {}
        self.d['chi0'] = chi0
        self.d['zbox'] = self.meta['Redshift']
        self.d['abox'] = self.meta['ScaleFactor']
        self.d['Lbox'] = self.meta['BoxSizeHMpc']
        self.d['OmM' ] = self.meta['Omega_M']
        #
        # The velocity<->distance conversion for RSD,
        # assuming LCDM.
        OmM  = self.meta['Omega_M']
        Eofz = lambda zz: np.sqrt( OmM*(1+zz)**3+(1-OmM) )
        self.d['velf'] = self.d['abox']*100*Eofz(self.d['zbox'])
        #
        # Later we will want a random number generator.  Use a
        # fixed seed to make this reproducable.
        self.rng  = np.random.default_rng(1)
        #
    def periodic(self,pos):
        """Periodically wrap pos into -0.5Lbox to 0.5Lbox."""
        Lbox,Lbox2 = self.d['Lbox'],0.5*self.d['Lbox']
        wrap = np.nonzero( pos>=Lbox2 )[0]
        if len(wrap)>0: pos[wrap] -= Lbox
        wrap = np.nonzero( pos<-Lbox2 )[0]
        if len(wrap)>0: pos[wrap] += Lbox
        return(pos)
        #
    def set_hod(self,params):
        """Assign the HOD parameters."""
        # For the LRG HOD sigma is defined with natural logs,
        # with the sqrt{2}.
        # Satellite numbers are ncen times ([M-kappa.Mcut]/M1)^alpha
        for k in params.keys():
            self.HOD_params['LRG_params'][k] = params[k]
            self.d[k] = params[k]
        #
    def generate(self):
        """Calls abacusutils to generate the mock sample."""
        newBall = AbacusHOD(self.sim_params,\
                            self.HOD_params,self.clustering_params)
        laes = newBall.run_hod(newBall.tracers,self.want_rsd,\
                               self.write_to_disk,Nthread=16)
        # Generate a redshift, including peculiar velocity.
        laes['LRG']['zred'] = self.periodic(laes['LRG'][ 'z']+\
                                            laes['LRG']['vz']/self.d['velf'])
        # We want some statistics on this sample.
        self.d['nobj'] = laes['LRG']['mass'].size
        self.d['nbar'] = self.d['nobj']/self.d['Lbox']**3
        self.d['ncen'] = laes['LRG']['Ncent']
        self.d['fsat'] = 1-float(self.d['ncen'])/float(self.d['nobj'])
        # and save the sample for later processing.  Could instead
        # return this and then pass it to other methods.  TBD.
        self.laes  = laes
        # Finally let's make a bitmask for the objects.  If
        # we use a single byte FITS interprets this as a boolean,
        # so we'll use a 16-bit mask.
        self.bitmask = np.zeros(self.d['nobj'],dtype='uint16')
        #
    def assign_lum(self,bright_frac):
        """Assigns Llya to the mock objects."""
        # Eventually we could do a per-object luminosity using
        # abundance matching with scatter to a Lucy-deconvolved
        # Schechter LF.  For now all we care about is bright vs.
        # faint, which we do randomly.  We could also weight by
        # a power of halo mass.
        rr = self.rng.uniform(size=self.d['nobj'])
        self.bitmask[rr<bright_frac] |= 1
        #
    def radial_select(self,chis):
        """Radial selection function, in chi."""
        chiv = self.sfn[:,1]
        pchi = self.sfn[:,2]
        prob = np.interp(chis,chiv,pchi,left=0,right=0)
        return(prob)
        #
    def select(self,diam,offset):
        """Select a small region of the box of diameter diam (radians)
        with the center of the box shifted by offset (fraction of box)."""
        # Also, at this point we should implement a z-dependent selection
        # function, but for now I'll use a hard zcut.
        Lbox  = self.d['Lbox']
        chi0  = self.d['chi0']
        Lside = diam * chi0
        gals  = self.laes['LRG']
        self.d['Lside'] = Lside
        # Shift the center of the box and select objects.
        xpos = self.periodic(gals['x'   ]+offset[0]*Lbox)
        ypos = self.periodic(gals['y'   ]+offset[1]*Lbox)
        zpos = self.periodic(gals['zred']+offset[2]*Lbox)
        rnum = self.rng.uniform(size=zpos.size)
        rsel = self.radial_select(zpos+chi0)
        in_survey = np.nonzero( (xpos>-0.5*Lside)&\
                                (xpos< 0.5*Lside)&\
                                (ypos>-0.5*Lside)&\
                                (ypos< 0.5*Lside)&\
                                (rnum<rsel) )[0]
        # Store the results for later use.
        self.d['nkeep'] = len(in_survey)
        self.xpos = xpos[in_survey]
        self.ypos = ypos[in_survey]
        self.zpos = zpos[in_survey]
        self.bitm = self.bitmask[in_survey]
        #
    def make_map(self,Nside=512):
        """Create the map."""
        # Now bin everything into a map...there are very few objects so we
        # don't need to be particularly clever.
        Lside= self.d['Lside']
        dmap = np.zeros( (Nside,Nside) )
        for i in range(self.d['nkeep']):
            ix  = int( ((self.xpos[i])/Lside+0.5)*Nside )
            iy  = int( ((self.ypos[i])/Lside+0.5)*Nside )
            dmap[ix,iy] += 1.0
        # and convert to overdensity -- there's a question here of
        # whether to use the mean density estimated from the whole
        # box or the map itself.  The latter is more "observational".
        dmap /= np.sum(dmap)/Nside**2
        dmap -= 1.0
        return(dmap)
        #
    def make_hdr(self):
        """Puts some information into a 'header' dictionary."""
        hdr = {}
        hdr['sim'    ] = self.sim_params['sim_name']
        for k in self.d.keys(): hdr[k] = self.d[k]
        return(hdr)
        #
    def write_map(self,dmap,Nside,diam,outfn):
        """Write the map to a FITS file, outfn."""
        # Put everything in a dictionary and write it to a FITS file.
        hdr = self.make_hdr()
        hdr['Lx'] = diam
        hdr['Ly'] = diam
        hdr['COMMENT'] = 'Distances in Mpc/h, angles in radians.'
        outdict = {}
        outdict['del'] = dmap.astype('f4')
        outdict['msk'] = np.ones_like(outdict['del'])
        tt = Table(outdict)
        for k in hdr.keys(): tt.meta[k] = hdr[k]
        tt.write(outfn,overwrite=True)
        #
    def write_cat(self,outfn):
        """Write a catalog of objects in the survey."""
        # Generate the angular coordinates, in degrees.
        # Assume small angles and plane projection.
        ichi = 1.0/self.d['chi0'] * 180./np.pi
        rra  = self.xpos*ichi
        dec  = self.ypos*ichi
        # and save them in a dictionary.
        hdr,outdict    = self.make_hdr(),{}
        hdr['COMMENT'] = 'Distances in Mpc/h'
        outdict['RA' ] = rra.astype('float32')
        outdict['DEC'] = dec.astype('float32')
        outdict['BITMASK']=self.bitm
        tt = Table(outdict)
        for k in hdr.keys(): tt.meta[k] = hdr[k]
        tt.write(outfn,overwrite=True)
        #






if __name__=="__main__":
    Nside  = 512
    diam   = 3.2 * np.pi/180.
    laes   = MockLAE('lae_base.yaml',4385.,'sfn.txt')
    params = {'logM_cut':11.75,'logM1':11.75+np.log10(5.),\
              'sigma':0.66,'kappa':0.33,'alpha':0.50}
    laes.set_hod(params)
    laes.generate()
    laes.assign_lum(0.5)
    laes.select(diam,[0.,0.,0.])
    dmap = laes.make_map(Nside)
    laes.write_map(dmap,Nside,diam,"mock_lae_map.fits")
    laes.write_cat("mock_lae_cat.fits")
    #
