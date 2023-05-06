#!/usr/bin/env python3
#
# Generates a mock catalog and computes all of the
# auto- and cross-spectra for this catalog so we
# can measure the bias(es).
#
import numpy as np
import yaml
import json
import gc

import random
import asdf
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta
from abacusnbody.data.bitpacked import unpack_rvint




def load_matter(sim_dir,Lbox,z,n_chunks,N_parts,f_down,type_AB='A'):
    """Loads the matter particles from the sub-sample files."""
    # estimate total number of particles and preselect indices
    N_all = 0
    N_offset = np.zeros(n_chunks, dtype=int)
    N_file = np.zeros(n_chunks, dtype=int)
    for i_chunk in range(n_chunks):
        print(i_chunk, n_chunks)
        # halo and field particles
        fn_halo = sim_dir+f'/halos/z{z:.3f}/halo_rv_{type_AB}/halo_rv_{type_AB}_{i_chunk:03d}.asdf'
        fn_field = sim_dir+f'/halos/z{z:.3f}/field_rv_{type_AB}/field_rv_{type_AB}_{i_chunk:03d}.asdf'
        N_this = asdf.open(fn_halo)['data']['rvint'].shape[0]+asdf.open(fn_field)['data']['rvint'].shape[0]
        N_offset[i_chunk] = N_all
        N_file[i_chunk] = N_this
        N_all += N_this
    gc.collect()
    print("offsets", N_offset, flush=True)
    print("per file", N_file, flush=True)
    print("all the particles in the halo and field files", N_all, flush=True)

    # global indices to keep
    inds_keep = random.sample(range(N_all), N_all//f_down)
    N_keep = len(inds_keep)
    print("N_keep", N_keep)
    pos_down = np.zeros((N_keep, 3), dtype=np.float32)
    vel_down = np.zeros((N_keep, 3), dtype=np.float32)

    # load the matter particles
    count = 0
    for i_chunk in range(n_chunks):
        print(i_chunk, n_chunks)
        # indices to keep in this chunk
        inds_keep_this = inds_keep - N_offset[i_chunk]
        inds_keep_this = inds_keep_this[(inds_keep_this >= 0) & (inds_keep_this < N_file[i_chunk])]

        # halo and field particles
        fn_halo = sim_dir+f'/halos/z{z:.3f}/halo_rv_{type_AB}/halo_rv_{type_AB}_{i_chunk:03d}.asdf'
        fn_field = sim_dir+f'/halos/z{z:.3f}/field_rv_{type_AB}/field_rv_{type_AB}_{i_chunk:03d}.asdf'

        halo_data = (asdf.open(fn_halo)['data'])['rvint']
        pos_halo, vel_halo = unpack_rvint(halo_data, Lbox, float_dtype=np.float32, velout=None)
        print("pos_halo = ", pos_halo[:5])

        field_data = (asdf.open(fn_field)['data'])['rvint']
        pos_field, vel_field = unpack_rvint(field_data, Lbox, float_dtype=np.float32, velout=None)
        print("pos_field = ", pos_field[:5])

        # stack halo and field particles
        pos_both = np.vstack((pos_halo, pos_field))
        vel_both = np.vstack((vel_halo, vel_field))

        pos_down[count:count+len(inds_keep_this)] = pos_both[inds_keep_this]
        vel_down[count:count+len(inds_keep_this)] = vel_both[inds_keep_this]
        count += len(inds_keep_this)
        del halo_data, pos_halo, vel_halo, field_data, pos_field, vel_field
        gc.collect()
    print("these two must be the same", count, pos_down.shape[0])
    pos_down = pos_down[:count]
    vel_down = vel_down[:count]
    return pos_down, vel_down






if __name__=="__main__":
    # Set up a random number generator with a
    # fixed seed for reproducability.
    rng = np.random.default_rng(1)
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
    n_chunks = 34
    #
    # additional parameter choices
    want_rsd      = HOD_params['want_rsd']
    want_rsd      = False
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
    # For the LRG HOD sigma is defined with natural logs,
    # with the sqrt{2}.
    # Satellite numbers are ncen times ([M-kappa.Mcut]/M1)^alpha
    #
    params = {'logM_cut':11.80,'logM1':11.80+np.log10(5.),\
              'sigma':0.66,'kappa':0.33,'alpha':0.50}
    #
    hodkeys = ['logM_cut','logM1','sigma','kappa','alpha']
    for k in hodkeys:
        HOD_params['LRG_params'][k] = params[k]
    hod     = [HOD_params['LRG_params'][k] for k in hodkeys]
    newBall = AbacusHOD(sim_params,HOD_params,clustering_params)
    mock_dict= newBall.run_hod(newBall.tracers,want_rsd,\
                               write_to_disk,Nthread=16)
    nobj  = mock_dict['LRG']['mass'].size
    ncen  = mock_dict['LRG']['Ncent']
    fsat  = 1-float(ncen)/float(nobj)
    #
    sampfact = 10
    mpos,vcel = load_matter(sim_params['sim_dir']+sim_params['sim_name'],\
                  meta['BoxSize'],sim_params['z_mock'],n_chunks,\
                  meta['NP'],sampfact,type_AB='A')
    mock_dict['matter'] = {}
    mock_dict['matter']['x'] = mpos[:,0]
    mock_dict['matter']['y'] = mpos[:,1]
    mock_dict['matter']['z'] = mpos[:,2]
    ndm = mock_dict['matter']['x'].size
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
    if False:
        pk3d  = newBall.compute_Pkmu(mock_dict,nbins_k=50,nbins_mu=11,\
                  k_hMpc_max=0.5,logk=False,num_cells=512,\
                  paste='TSC',compensated=True,interlaced=True)
        mu    = pk3d['mu_binc']
        dmu   = 1.0/len(mu)
        kk    = pk3d['k_binc']
        pkmu  = pk3d['LRG_LRG']
        pkgg  = (2*0+1)*np.dot(pkmu,1.0*(0*mu**2+1))*dmu
        pkmu  = pk3d['LRG_matter']
        pkgm  = (2*0+1)*np.dot(pkmu,1.0*(0*mu**2+1))*dmu
        pkmu  = pk3d['matter_matter']
        pkmm  = (2*0+1)*np.dot(pkmu,1.0*(0*mu**2+1))*dmu
        bka   = np.sqrt(pkgg/pkmm)
        bkx   = pkgm/pkmm
        #
        with open("lae_bk.txt","w") as fout:
            fout.write("# Real-space Fourier biases.\n")
            fout.write("# "+sim_name+"\n")
            fout.write("# z={:.2f}\n".format(sim_params['z_mock']))
            fout.write("# {:>10s} {:>15s} {:>15s}\n".\
                       format("k[h/Mpc]","ba","bx","Pmm"))
            for i in range(kk.size):
                fout.write("{:12.4e} {:15.5e} {:15.5e}".\
                        format(kk[i],bka[i],bkx[i],pkmm[i]))
    #
    if True:
        # This is very slow if there are many objects and
        # the maximum lag is too large.
        maxobj  = 8000000 # Should be an integer.
        if nobj>maxobj:
            print("Have nobj=",nobj," downsampling to ",maxobj,flush=True)
            inds = rng.choice(nobj,size=maxobj,replace=False)
            samp = 'LRG'
            for k in ['x','y','z','vx','vy','vz','mass','id']:
                mock_dict[samp][k] = mock_dict[samp][k][inds]
            print("Have ndm=",ndm," downsampling to ",maxobj,flush=True)
            inds = rng.choice(ndm,size=maxobj,replace=False)
            samp = 'matter'
            for k in ['x','y','z']:
                mock_dict[samp][k] = mock_dict[samp][k][inds]
        xi3d = newBall.compute_multipole(mock_dict,rpbins,pimax,dpi)
        xiell= xi3d['LRG_LRG']
        xigg = xiell[1*len(Rcen):2*len(Rcen)]
        xiell= xi3d['LRG_matter']
        xigm = xiell[1*len(Rcen):2*len(Rcen)]
        xiell= xi3d['matter_matter']
        ximm = xiell[1*len(Rcen):2*len(Rcen)]
        bra  = np.sqrt(xigg/ximm)
        brx  = xigm/ximm
        #
        with open("lae_br.txt","w") as fout:
            fout.write("# Real-space configuration biases.\n")
            fout.write("# "+sim_name+"\n")
            fout.write("# z={:.2f}\n".format(sim_params['z_mock']))
            fout.write("# {:>10s} {:>15s} {:>15s}\n".\
                       format("r[Mpc/h]","ba","bx","xi_mm"))
            for i in range(Rcen.size):
                fout.write("{:12.4e} {:15.5e} {:15.5e}".\
                        format(Rcen[i],bra[i],brx[i],ximm[i]))
    #
