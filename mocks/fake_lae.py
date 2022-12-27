#!/usr/bin/env python3
#
import numpy as np
import yaml
#
from astropy.table              import Table
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta



class MockLAE:
    """A class to handle making mock LAE samples from halo catalogs."""
    # We may want to do a better job of holding the information in
    # this class.  Right now it's not very well thought through what is
    # kept internally vs. what is passed to each method.
    def __init__(self,yaml_file,chi0,dchi):
        """Set up the class."""
        # Load the config file and parse in relevant parameters
        config          = yaml.safe_load(open(yaml_file))
        self.sim_params = config['sim_params']
        self.HOD_params = config['HOD_params']
        self.clustering_params = config['clustering_params']
        #
        # Get the metaparameters for the simulation.
        self.meta = get_meta(self.sim_params['sim_name'],\
                             redshift=self.sim_params['z_mock'])
        #
        # additional parameter choices
        self.want_rsd      = self.HOD_params['want_rsd']
        self.write_to_disk = self.HOD_params['write_to_disk']
        #
        # Save some configuration parameters for later use.
        self.d = {}
        self.d['chi0'] = chi0
        self.d['dchi'] = dchi
        self.d['zcen'] = self.meta['Redshift']
        self.d['acen'] = self.meta['ScaleFactor']
        self.d['Lbox'] = self.meta['BoxSizeHMpc']
        self.d['OmM' ] = self.meta['Omega_M']
        #
        # The velocity<->distance conversion for RSD,
        # assuming LCDM.
        OmM  = self.meta['Omega_M']
        Eofz = lambda zz: np.sqrt( OmM*(1+zz)**3+(1-OmM) )
        self.d['velf'] = self.d['acen']*100*Eofz(self.d['zcen'])
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
        # Finally let's make a bitmask for the objects.
        self.bitmask = np.zeros(self.d['nobj'],dtype='b')
        #
    def assign_lum(self,bright_frac):
        """Assigns Llya to the mock objects."""
        # Eventually we could do a per-object luminosity using
        # abundance matching with scatter to a Lucy-deconvolved
        # Schechter LF.  For now all we care about is bright vs.
        # faint, which we do randomly.  We could also weight by
        # a power of halo mass.
        rr = self.rng.uniform(size=self.d['nobj'])
        self.bitmask[rr<bright_frac] |= 2
        #
    def select(self,diam):
        """Select a small region of the box of diameter diam (radians)."""
        # We may eventually want to do this with arbitrary centers or
        # shifts of the coordinates to generate many maps.  Leave that
        # for now.
        # Also, at this point we should implement a z-dependent selection
        # function, but for now I'll use a hard zcut.
        Lside = diam * self.d['chi0']
        depth = self.d['dchi']
        gals  = self.laes['LRG']
        self.d['Lside'] = Lside
        in_survey = np.nonzero( (gals['x']>-0.5*Lside)&\
                                (gals['x']< 0.5*Lside)&\
                                (gals['y']>-0.5*Lside)&\
                                (gals['y']< 0.5*Lside)&\
                                (gals['zred']>-0.5*depth)&\
                                (gals['zred']< 0.5*depth) )[0]
        self.d['nkeep'] = len(in_survey)
        self.bitmask &= 254
        self.bitmask[in_survey] |= 1
        #
    def make_map(self,Nside=512):
        """Create the map."""
        # Now bin everything into a map...there are very few objects so we
        # don't need to be particularly clever.
        gals = self.laes['LRG']
        Lside= self.d['Lside']
        dmap = np.zeros( (Nside,Nside) )
        for i in np.nonzero( self.bitmask&1 )[0]:
            ix  = int( ((gals['x'][i])/Lside+0.5)*Nside )
            iy  = int( ((gals['y'][i])/Lside+0.5)*Nside )
            dmap[ix,iy] += 1.0
        # and convert to overdensity.
        dmap /= np.sum(dmap)/Nside**2
        dmap -= 1.0
        return(dmap)
        #
    def make_hdr(self):
        """Puts some information into a 'header' dictionary."""
        hdr = {}
        hdr['sim'    ] = self.sim_params['sim_name']
        hdr['COMMENT'] = 'Distances in Mpc/h, angles in radians.'
        for k in self.d.keys(): hdr[k] = self.d[k]
        return(hdr)
        #
    def write_map(self,dmap,Nside,diam,outfn):
        """Write the map to a FITS file, outfn."""
        # Put everything in a dictionary and write it to a FITS file.
        hdr = self.make_hdr()
        hdr['Lx'] = diam
        hdr['Ly'] = diam
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
        gals = self.laes['LRG']
        rra  = []
        dec  = []
        for i in np.nonzero( self.bitmask&1 )[0]:
            rra.append( gals['x'][i]*ichi )
            dec.append( gals['y'][i]*ichi )
        # and save them in a dictionary.
        hdr,outdict    = self.make_hdr(),{}
        outdict['RA' ] = np.array(rra)
        outdict['DEC'] = np.array(dec)
        tt = Table(outdict)
        for k in hdr.keys(): tt.meta[k] = hdr[k]
        tt.write(outfn,overwrite=True)
        #






if __name__=="__main__":
    Nside  = 512
    diam   = 3.2 * np.pi/180.
    laes   = MockLAE('lae_base.yaml',4385.,30.)
    params = {'logM_cut':11.75,'logM1':11.75+np.log10(5.),\
              'sigma':0.66,'kappa':0.33,'alpha':0.33}
    laes.set_hod(params)
    laes.generate()
    laes.assign_lum(0.5) # Ignore output for now.
    laes.select(diam)
    dmap = laes.make_map(Nside)
    laes.write_map(dmap,Nside,diam,"mock_lae_map.fits")
    laes.write_cat("mock_lae_cat.fits")
    #
