#!/usr/bin/env python3
#
import numpy as np
import yaml
import json

from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta


# These are the keys we will copy from the simulation meta-data into
# the output JSON file.
simkeys = ['n_s', 'omega_b', 'omega_cdm', 'omega_ncdm', 'N_ncdm', 'N_ur',\
           'H0', 'w0', 'wa', 'w', 'Omega_DE', 'Omega_K', 'Omega_M', 'Omega_Smooth', \
           'Redshift', 'ScaleFactor', \
           'OmegaNow_DE', 'OmegaNow_K', 'OmegaNow_m', \
           'f_growth', 'fsmooth', 'Growth', 'Growth_on_a_n', \
           'SimComment', 'SimName', 'SimSet', \
           'BoxSize', 'NP', 'BoxSizeHMpc', 'HubbleTimeHGyr', \
           'ParticleMassHMsun']
#
# Load the config file and parse in relevant parameters
path2config= './lae_base.yaml'
config     = yaml.safe_load(open(path2config))
sim_params = config['sim_params']
HOD_params = config['HOD_params']
clustering_params = config['clustering_params']
#
# Get the metaparameters for the simulation.
meta = get_meta(sim_params['sim_name'],redshift=sim_params['z_mock'])
#
# additional parameter choices
want_rsd      = HOD_params['want_rsd']
write_to_disk = HOD_params['write_to_disk']
want_zcv      = True
#
# Load the rp pi binning from the config file. Note pimax and pi_bin_size are ints.
bin_params = clustering_params['bin_params']
rpbins     = np.logspace(bin_params['logmin'],bin_params['logmax'],bin_params['nbins']+1)
pimax,dpi  = clustering_params['pimax'],clustering_params['pi_bin_size']
# work out the rp and pi bin centers (assume log binning as above)
Rcen = np.sqrt(rpbins[1:]*rpbins[:-1])
Zcen = np.arange(0.0,float(pimax),dpi) + 0.5*dpi
#
# Now loop over HODs writing out HOD, wp(R), xi0, etc. for each.
#
hodkeys = ['logM_cut','logM1','sigma','kappa','alpha']
plateau = np.log10(5.0)
dats = []
for lgMcut in [11.50,11.75,12.00]:
    for sigm in [0.5]:
        HOD_params['LRG_params']['logM_cut'] = lgMcut
        HOD_params['LRG_params']['logM1'   ] = lgMcut + plateau
        HOD_params['LRG_params']['sigma'   ] = sigm
        HOD_params['LRG_params']['kappa'   ] = 0.333
        HOD_params['LRG_params']['alpha'   ] = 0.333
        hod      = [HOD_params['LRG_params'][k] for k in hodkeys]
        newBall  = AbacusHOD(sim_params,HOD_params,clustering_params)
        mock_dict= newBall.run_hod(newBall.tracers,want_rsd,write_to_disk,Nthread=16)
        nobj  = mock_dict['LRG']['mass'].size
        ncen  = mock_dict['LRG']['Ncent']
        fsat  = 1-float(ncen)/float(nobj)
        #wpR   = newBall.compute_wp(mock_dict,rpbins,pimax,dpi)['LRG_LRG']
        xiell = newBall.compute_multipole(mock_dict,rpbins,pimax,dpi)['LRG_LRG']
        wpR   = xiell[0*len(Rcen):1*len(Rcen)]
        xi0   = xiell[1*len(Rcen):2*len(Rcen)]
        xi2   = xiell[2*len(Rcen):3*len(Rcen)]
        #
        if want_zcv:
            # Compute variance reduced spectra.
            zcv_dict = newBall.apply_zcv(mock_dict,config)
            kk = zcv_dict['k_binc']
            pkl= zcv_dict['Pk_tr_tr_ell_zcv'] # variance-reduced multipoles
            pk0= pkl[0]
            pk2= pkl[1]
        else:
            # Just something to put in the file.
            kk,pk0,pk2 = np.zeros(1),np.zeros(1),np.zeros(1)
        #
        dats.append({'hod':hod,'nobj':nobj,'fsat':fsat,\
                     'wp':wpR.tolist(),'xi0':xi0.tolist(),'xi2':xi2.tolist(),\
                     'pk0':pk0.tolist(),'pk2':pk2.tolist()})
outdict = {}
for k in simkeys: outdict[k]=meta[k]
outdict['in_red' ] = want_rsd
outdict['R'      ] = Rcen.tolist()
outdict['k'      ] = kk.tolist()
outdict['hodkeys'] = hodkeys
outdict['mocks'  ] = dats
outsuf = 's' if want_rsd else 'r'
outsim = sim_params['sim_name']
outsim = outsim[outsim.find("_c0"):] + '_'
with open("lae_clustering"+outsim+outsuf+".json","w") as fout:
    json.dump(outdict,fout,indent=2)
