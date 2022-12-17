#!/usr/bin/env python3
#
import numpy as np
import yaml
#
from lcdm                       import LCDM
from astropy.table              import Table
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta
#
# Load the config file and parse in relevant parameters
path2config= './lae_base.yaml'
config     = yaml.safe_load(open(path2config))
sim_params = config['sim_params']
HOD_params = config['HOD_params']
clustering_params = config['clustering_params']
#
print("Making a mock LAE map from ",sim_params['sim_name'],flush=True)
#
# Get the metaparameters for the simulation.
meta = get_meta(sim_params['sim_name'],redshift=sim_params['z_mock'])
#
# additional parameter choices
want_rsd      = HOD_params['want_rsd']
write_to_disk = HOD_params['write_to_disk']
#
# Set up an LCDM instance and compute the distance to the
# midpoint of the box, boxsize, etc.
zcen = meta['Redshift']
acen = meta['ScaleFactor']
Lbox = meta['BoxSizeHMpc']
cc   = LCDM(meta['Omega_M'])
chi0 = cc.chi_of_z(zcen)
velf = acen*100*cc.E_of_z(zcen)
#
print("zcen=",zcen,", acen=",acen,flush=True)
print("Lbox=",Lbox,"Mpc/h.",flush=True)
print("chi0=",chi0,"Mpc/h.",flush=True)
#
# Generate the mock objects from an HOD.
# For the LRG HOD sigma is defined with natural logs,
# with the sqrt{2}.
# Satellite numbers are ncen times ([M-kappa.Mcut]/M1)^alpha
#
HOD_params['LRG_params']['logM_cut'] = 11.75
HOD_params['LRG_params']['logM1'   ] = 11.75 + np.log10(5.)
HOD_params['LRG_params']['sigma'   ] = 1.00
HOD_params['LRG_params']['kappa'   ] = 0.33
HOD_params['LRG_params']['alpha'   ] = 0.33
newBall  = AbacusHOD(sim_params,HOD_params,clustering_params)
mock_dict= newBall.run_hod(newBall.tracers,want_rsd,\
                           write_to_disk,Nthread=16)
# We want some statistics on this sample.
nobj  = mock_dict['LRG']['mass'].size
ncen  = mock_dict['LRG']['Ncent']
fsat  = 1-float(ncen)/float(nobj)
# Go into redshift space and periodically wrap.
mock_dict['LRG']['zred'] = mock_dict['LRG']['z']+mock_dict['LRG']['vz']/velf
wrap = np.nonzero( mock_dict['LRG']['zred']>Lbox )[0]
if len(wrap)>0: mock_dict['LRG']['zred'][wrap] -= Lbox
wrap = np.nonzero( mock_dict['LRG']['zred']<0 )[0]
if len(wrap)>0: mock_dict['LRG']['zred'][wrap] += Lbox
# Now select a small region of the box.
# At this point we should implement a z-dependent selection
# function, but for now I'll use a hard zcut.
diam  = 3.25 * np.pi/180. # Diameter of field, in radians.
Lside = diam * chi0
depth = 30.0 # Mpc/h
print("Lside=",Lside,"Mpc/h.  Depth=",depth,"Mpc/h.",flush=True)
#
in_survey = np.nonzero( (mock_dict['LRG']['x']>0.5*(Lbox-Lside))&\
                        (mock_dict['LRG']['x']<0.5*(Lbox+Lside))&\
                        (mock_dict['LRG']['y']>0.5*(Lbox-Lside))&\
                        (mock_dict['LRG']['y']<0.5*(Lbox+Lside))&\
                        (mock_dict['LRG']['zred']>0.5*(Lbox-depth))&\
                        (mock_dict['LRG']['zred']<0.5*(Lbox+depth)) )[0]
print("Keeping ",len(in_survey)," objects in survey.",flush=True)
# Now bin everything into a map...there are very few objects so we
# don't need to be particularly clever.
Nside= 512
dmap = np.zeros( (Nside,Nside) )
for i in in_survey:
    ix  = int( ((mock_dict['LRG']['x'][i]-0.5*Lbox)/Lside+0.5)*Nside )
    iy  = int( ((mock_dict['LRG']['y'][i]-0.5*Lbox)/Lside+0.5)*Nside )
    dmap[ix,iy] += 1.0
# and convert to overdensity.
dmap /= np.sum(dmap)/Nside**2
dmap -= 1.0
# Put everything in a dictionary and write it to a FITS file.
hdr = {}
hdr['sim'    ] = sim_params['sim_name']
hdr['zcen'   ] = zcen
hdr['Lbox'   ] = Lbox
hdr['chi0'   ] = chi0
hdr['Lside'  ] = Lside
hdr['OmegaM' ] = meta['Omega_M']
hdr['nbar'   ] = nobj/Lbox**3
hdr['fsat'   ] = fsat
hdr['Lx'     ] = diam
hdr['Ly'     ] = diam
hdr['COMMENT'] = 'Distances in Mpc/h, angles in radians.'
outdict = {}
outdict['del'   ] = dmap.astype('f4')
outdict['msk'   ] = np.ones_like(outdict['del'])
#
print(hdr)
#
outfn = "mock_lae.fits"
tt = Table(outdict)
for k in hdr.keys(): tt.meta[k] = hdr[k]
tt.write(outfn,overwrite=True)
#
