#!/usr/bin/env python
#
# Makes a table of some useful quantities for the LAE paper.
#
import numpy  as     np
from   classy import Class



# Set up a default set of parameters.
params = {'output': 'tCl lCl mPk',
          'l_max_scalars': 3000,
          'lensing': 'yes',
          'P_k_max_h/Mpc': 50.,
          'non linear':'halofit',
          'z_pk': '0.0,6',
          'A_s': 2.083e-09,
          'n_s': 0.9649,
          'alpha_s': 0.,
          'h': 0.6770,
          'N_ur': 2.0328,
          'N_ncdm': 1,
          'h': 0.6736,
          'tau_reio': 0.0568,
          'omega_b': 0.02237,
          'omega_cdm': 0.12,
          'omega_ncdm': 0.0006442,
          'Omega_k': 0.}
#
cosmo = Class()
cosmo.set(params)
cosmo.compute()
#
wb = cosmo.omega_b()
#
print("OmegaM=",cosmo.Omega_m())
print("sigma8=",cosmo.sigma8())
print("hubble=",cosmo.h())
print("omegab=",wb)
#
print(cosmo.get_current_derived_parameters(['H0','Omega_Lambda',\
                                            'age','conformal_age','Neff',\
                                            'z_reio','100*theta_s','rs_rec','rs_d']))
print("")
#
# Compute the Zeldovich displacement and hence k_{nl}.
hub= cosmo.h() # To convert to "conventional" Mpc/h units.
kk = np.logspace(-3.5,0.3,250)
#
for zz,dz in zip([2.4,3.123,4.5],[0.03,0.03,0.04]):
    pk = np.array( [cosmo.pk(k*hub,zz)*hub**3 for k in kk] )
    pl = np.array( [cosmo.pk_lin(k*hub,zz)*hub**3 for k in kk] )
    #
    Ez   = cosmo.Hubble(zz)/cosmo.Hubble(0.0)
    chi0 = cosmo.comoving_distance(zz)*hub
    dcdz = 2997.925/Ez
    dchi = dcdz*(2*dz)
    #
    # Compute the Zeldovich displacement and hence k_{nl}.
    knl= 1/np.sqrt( np.trapz(pl,x=kk)/6./np.pi**2 )
    #
    outstr = "z={:.3f}, chi0={:.0f}, dchi={:.1f}, dchi/dz={:.1f} knl={:.2f}".\
             format(zz,chi0,dchi,dcdz,knl)
    print(outstr)
#
