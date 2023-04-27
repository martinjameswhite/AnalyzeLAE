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
want_zcv      = False
#
# Load the rp pi binning from the config file. Note pimax and pi_bin_size are ints.
bin_params = clustering_params['bin_params']
rpbins     = np.logspace(bin_params['logmin'],\
                         bin_params['logmax'],bin_params['nbins']+1)
pimax,dpi  = clustering_params['pimax'],clustering_params['pi_bin_size']
# work out the rp and pi bin centers (assume log binning as above)
Rcen = np.sqrt(rpbins[1:]*rpbins[:-1])
Zcen = np.arange(0.0,float(pimax),dpi) + 0.5*dpi
#
# Now loop over HODs writing out HOD, wp(R), xi0, etc. for each.
# For the LRG HOD sigma is defined with natural logs,
# with the sqrt{2}.
# Satellite numbers are ncen times ([M-kappa.Mcut]/M1)^alpha
#
if want_rsd:
    lgMc_list = [11.50,11.75,12.00,12.25]
    alph_list = [0.33,0.50,0.66]
    sigm_list = [0.50,0.66]
    plat_list = [5.,10.]
else:
    lgMc_list = [11.50,11.75,12.00,12.25]
    alph_list = [0.33,0.50,0.66]
    sigm_list = [0.50,0.66]
    plat_list = [5.,10.]
#
maxobj  = 8000000 # Should be an integer.
hodkeys = ['logM_cut','logM1','sigma','kappa','alpha']
dats    = []
for lgMcut in lgMc_list:
  for alph in alph_list:
    for sigm in sigm_list:
      for plateau in plat_list:
        HOD_params['LRG_params']['logM_cut'] = lgMcut
        HOD_params['LRG_params']['logM1'   ] = lgMcut + np.log10(plateau)
        HOD_params['LRG_params']['sigma'   ] = sigm
        HOD_params['LRG_params']['kappa'   ] = 0.333
        HOD_params['LRG_params']['alpha'   ] = alph
        hod      = [HOD_params['LRG_params'][k] for k in hodkeys]
        newBall  = AbacusHOD(sim_params,HOD_params,clustering_params)
        mock_dict= newBall.run_hod(newBall.tracers,want_rsd,\
                                   write_to_disk,Nthread=16)
        nobj  = mock_dict['LRG']['mass'].size
        ncen  = mock_dict['LRG']['Ncent']
        fsat  = 1-float(ncen)/float(nobj)
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
        if False: # Old code.
            pk3d  = newBall.compute_Pkmu(mock_dict,nbins_k=50,nbins_mu=11,\
                      k_hMpc_max=0.5,logk=False,num_cells=512,\
                      paste='TSC',compensated=True,interlaced=True)
            mu    = pk3d['mu_binc']
            dmu   = 1.0/len(mu)
            kk    = pk3d['k_binc']
            pkmu  = pk3d['LRG_LRG']
            pk0   = (2*0+1)*np.dot(pkmu,1.0*(0*mu**2+1))*dmu
            pk2   = (2*2+1)*np.dot(pkmu,0.5*(3*mu**2-1))*dmu
        #
        if nobj>maxobj:
            print("Have nobj=",nobj," downsampling to ",maxobj)
            rng  = np.random.default_rng()
            inds = rng.choice(nobj,size=maxobj,replace=False)
            for k in ['x','y','z','vx','vy','vz','mass','id']:
                mock_dict['LRG'][k] = mock_dict['LRG'][k][inds]
        #wpR   = newBall.compute_wp(mock_dict,rpbins,pimax,dpi)['LRG_LRG']
        xiell = newBall.compute_multipole(mock_dict,rpbins,pimax,dpi)['LRG_LRG']
        wpR   = xiell[0*len(Rcen):1*len(Rcen)]
        xi0   = xiell[1*len(Rcen):2*len(Rcen)]
        xi2   = xiell[2*len(Rcen):3*len(Rcen)]
        #
        dats.append({'hod':hod,'nobj':nobj,'fsat':fsat,\
                     'wp':wpR.tolist(),'xi0':xi0.tolist(),'xi2':xi2.tolist(),\
                     'pk0':pk0.tolist(),'pk2':pk2.tolist()})
# Now write out our answers.
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
