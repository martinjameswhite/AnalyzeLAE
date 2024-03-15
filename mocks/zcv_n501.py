#!/usr/bin/env python3
#
# For a single mock catalog, compute the Fourier- and
# configuration-space clustering using ZCV.
#
import numpy as np
import yaml
import json

from   scipy.interpolate    import InterpolatedUnivariateSpline as Spline

from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta



if __name__=="__main__":
    # Set up a random number generator with a
    # fixed seed for reproducability.
    rng    = np.random.default_rng(1)
    # Define the mock catalog and HOD.
    # For the LRG HOD sigma is defined with natural logs,
    # with the sqrt{2}.
    # Satellite numbers are ncen times ([M-kappa.Mcut]/M1)^alpha
    hodkeys= ['logM_cut','logM1','sigma','kappa','alpha']
    hod    = {'logM_cut':11.00,'logM1':12.30,\
              'sigma':0.66,'kappa':1.00,'alpha':0.33}
    print("HOD: ",hod,flush=True)
    #
    # Load the config file and parse in relevant parameters
    path2config= './lae_n501.yaml'
    config     = yaml.safe_load(open(path2config))
    sim_params = config['sim_params']
    HOD_params = config['HOD_params']
    clustering_params = config['clustering_params']
    #
    # Get the metaparameters for the simulation.
    meta = get_meta(sim_params['sim_name'],redshift=sim_params['z_mock'])
    #
    # additional parameter choices - override YAML file.
    want_rsd      = True
    write_to_disk = False
    want_zcv      = True  # Must be true.
    #
    # Copy the HOD parameters.
    for k in hodkeys:
        HOD_params['LRG_params'][k] = hod[k]
    #
    newBall  = AbacusHOD(sim_params,HOD_params,clustering_params)
    mock_dict= newBall.run_hod(newBall.tracers,want_rsd,\
                               write_to_disk,Nthread=16)
    nobj  = mock_dict['LRG']['mass'].size
    ncen  = mock_dict['LRG']['Ncent']
    fsat  = 1-float(ncen)/float(nobj)
    print("nobj=",nobj,", ncen=",ncen,", fsat=",fsat,flush=True)
    #
    if want_zcv:
        # Compute variance reduced spectra.
        zcv_dict = newBall.apply_zcv(mock_dict,config)
        kk  = zcv_dict['k_binc']
        pkl = zcv_dict['Pk_tr_tr_ell_zcv'] # variance-reduced multipoles
        pk0 = pkl[0,:]
        pk2 = pkl[1,:]
        pk4 = pkl[2,:]
        bb  = zcv_dict['bias'][0] + 1.0 # Convert to Eulerian bias.
        print("Pk: ",pk0,flush=True)
        print("bE: ",bb,flush=True)
    #
    if False: # Old code.
        maxobj  = 8000000 # Should be an integer.
        print("Have nobj=",nobj," downsampling to ",maxobj,flush=True)
        rng  = np.random.default_rng()
        inds = rng.choice(nobj,size=maxobj,replace=False)
        for k in ['x','y','z','vx','vy','vz','mass','id']:
            mock_dict['LRG'][k] = mock_dict['LRG'][k][inds]
    if want_zcv:
        zcv_dict = newBall.apply_zcv_xi(mock_dict,\
                                        config,load_presaved=False)
        xiell    = zcv_dict['Xi_tr_tr_ell_zcv']
        rbinc    = zcv_dict['r_binc']
        wpR      = np.array([0])
        xi0      = xiell[0,:]
        xi2      = xiell[1,:]
        xi4      = xiell[2,:]
        bb       = zcv_dict['bias'][0] + 1.0 # Convert to Eulerian bias.
        print("xi: ",xi0,flush=True)
        print("bE: ",bb,flush=True)
    #
    dats={'nobj':nobj,'fsat':fsat,'bias':bb,'wp':wpR.tolist(),\
          'xi0':xi0.tolist(),'xi2':xi2.tolist(),'xi4':xi4.tolist(),\
          'pk0':pk0.tolist(),'pk2':pk2.tolist(),'pk4':pk4.tolist()}
# Now write out our answers.
outdict = {}
outdict['in_red' ] = want_rsd
outdict['r'      ] = rbinc.tolist()
outdict['k'      ] = kk.tolist()
outdict['hodkeys'] = hodkeys
outdict['hod'    ] = hod
for k in dats.keys(): outdict[k]=dats[k]
outsuf = 's' if want_rsd else 'r'
outsim = sim_params['sim_name']
outsim = outsim[outsim.find("_c0"):] + '_'
with open("lae_n501"+outsim+outsuf+"_zcv.json","w") as fout:
    json.dump(outdict,fout,indent=2)
