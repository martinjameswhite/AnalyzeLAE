#!/usr/bin/env python
#
# Use the pre-computed alms for the "data" and "randoms" to
# make a power spectrum and window function (and plot them).
# This code is very slow!
# 
#
import numpy  as np
import healpy as hp
import time
import json
import sys
#
import matplotlib.pyplot as plt

# Load the fiducial interloper fraction information.
sys.path.append('/pscratch/sd/m/mwhite/AnalyzeLAE/odin/')
from fiducial import fint_dict

# Use the "C" variant of the direcsht code.
sys.path.append('/pscratch/sd/m/mwhite/direct_sht/csht/')
from mask_deconvolution import MaskDeconvolution

# Set up the list of filters and samples to process.
flist,slist = ["N501"],["s3"]

#
# Load the pre-computed alms.
tt,Nl = {},0
for ifilt in flist:
    tt[ifilt] = json.load(open("odin_s3_cosmos_"+ifilt+"_alm.json","r"))
    Nl = tt[ifilt]['Nl']
    #
    print(ifilt,", SHT has Nl=",Nl,", Nx=",tt[ifilt]['Nx'],\
          " and xmax=",tt[ifilt]['xmax'],flush=True)
    # Compute the harmonic-space description of the data and randoms.
    hdat = np.array(tt[ifilt]['dlm_re'])+1j*np.array(tt[ifilt]['dlm_im'])
    hran = np.array(tt[ifilt]['rlm_re'])+1j*np.array(tt[ifilt]['rlm_im'])
    hran*= hdat[0]/hran[0]
    tt[ifilt]['hdat'] = hdat
    tt[ifilt]['hran'] = hran
    # Calculate the angular power spectrum of the randoms
    # After subtracting shot noise, this gives us the
    # angular power spectrum of the window function 
    sn = 0.0 # 1.0/float(nrand)/(4*np.pi) * float(ndata)**2
    tt[ifilt]['wl'] = hp.alm2cl(tt[ifilt]['hran']) - sn
#
# Plot the "raw" window, unnormalized and on a log y-axis.
#
fig,ax = plt.subplots(1,1,figsize=(4.5,3.0))
icol = 0
for ifilt in flist:
    wl   = tt[ifilt]['wl']
    ell  = np.arange(wl.size)
    ax.plot(ell[1:],wl[1:],color='C'+str(icol),alpha=0.7,label=ifilt)
    icol = (icol+1)%10
#
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$W_\ell$')
#
plt.tight_layout()
plt.savefig('harmonic_rawW.pdf')
#
#
# Plot the window, normalized to unity at low ell
# and on a linear y-axis.
#
fig,ax = plt.subplots(1,1,figsize=(4.5,3.0))
icol = 0
for ifilt in flist:
    wl   = tt[ifilt]['wl']
    ell  = np.arange(wl.size)
    ax.plot(ell,wl/wl[0],color='C'+str(icol),alpha=0.7,label=ifilt)
    icol = (icol+1)%10
#
ax.legend()
ax.set_xlim(0,750)
ax.set_ylim(1e-4,1.05)
ax.set_yscale('log')
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$W_\ell / W_0$')
#
plt.tight_layout()
plt.savefig('harmonic_window.pdf')
#
#
# Get harmonic coefficients and measure raw Cls
for ifilt in flist:
    hdcl = hp.alm2cl(tt[ifilt]['hdat'])
    hrcl = hp.alm2cl(tt[ifilt]['hran'])
    hdif = hp.alm2cl(tt[ifilt]['hdat']-tt[ifilt]['hran'])
    print(ifilt,", hdif.size=",hdif.size,flush=True)
    print(ifilt,", hatC in range [{:e},{:e}]".\
          format(np.min(hdif[1:]),np.max(hdif[1:])),flush=True)
    #
    tt[ifilt]['hdcl'] = hdcl
    tt[ifilt]['hrcl'] = hrcl
    tt[ifilt]['hdif'] = hdif
#
# Plot them.
fig,ax = plt.subplots(1,1,figsize=(8,5))
icol   = 0
for ifilt in flist:
    hdcl = tt[ifilt]['hdcl']
    hrcl = tt[ifilt]['hrcl']
    hdif = tt[ifilt]['hdif']
    ells = np.arange(1,hdif.size)
    ax.loglog(ells,hdcl[1:],ls=':' ,color='C'+str(icol))
    ax.loglog(ells,hrcl[1:],ls='--',color='C'+str(icol))
    ax.loglog(ells,hdif[1:],ls='-' ,color='C'+str(icol),label=ifilt)
    icol = (icol+1)%10
#
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$C_\ell$')
#
ax.legend()
plt.tight_layout()
plt.savefig('harmonic_allC.pdf')
#


#
# Initialize an instance of the MaskDeconvolution. This will let us deconvolve
# the mask-induced mode-coupling of the pseudo-Cls, convolved the theory to
# enable apples-to-apples comparisons, and provide binning functionality
# This is kind of slow.
MaskL,Lmin,NperBin = 2048,100,200
for ifilt in flist:
    print("Working on ",ifilt,flush=True)
    print("Intializing MaskDeconvolution class with MaskL=",MaskL,flush=True)
    MD = MaskDeconvolution(MaskL,wl)
    print("Done intializing MaskDeconvolution class.",flush=True)
    print(time.asctime(),flush=True)
    #
    bins = MD.binning_matrix('linear',Lmin,NperBin)
    print("bins.shape=",bins.shape,flush=True)
    Mbl  = MD.window_matrix(bins,mode='normalization')
    # Look at the sums over ell.
    print("\nMbl entries in range {:e} to {:e}:".\
          format(np.min(Mbl),np.max(Mbl)),flush=True)
    print("Row sums of Mbl:")
    print(Mbl.sum(axis=1),flush=True)
    #
    tt[ifilt]['Mbl'] = Mbl
    #
    print("MD.Mll.shape=",MD.Mll.shape,flush=True)
    print("MD.Mll entries in range {:e} to {:e}".\
          format(np.min(MD.Mll),np.max(MD.Mll)),flush=True)
#
#
# Plot the matrix:
fig,ax = plt.subplots(1,1,figsize=(8,6))
for i in range(Mbl.shape[0]):
    ax.plot(Mbl[i,:],alpha=0.5,label='Bin '+str(i))
ax.legend()
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$W_{b\ell}$')
plt.tight_layout()
plt.savefig('harmonic_wbl.pdf')
#
# and plot the mode-coupling matrix itself.
fig,ax = plt.subplots(1,1,figsize=(8,8))
mode_coupling = np.log( 1+MD.Mll.clip(0,1e30) )
ax.imshow(mode_coupling.T,origin='lower')
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$\ell^\prime$')
plt.tight_layout()
plt.savefig('harmonic_mll.pdf')


# Do the binned Cls.
for ifilt,isamp in zip(flist,slist):
    # Now work out the binned and normalized pseudo-spectrum.
    binned_ells,hdif_normed = MD(tt[ifilt]['hdif'],bins)
    # Correct for the interlopers.
    hdif_normed /= (1.0-fint_dict[ifilt][isamp])**2
    tt[ifilt]['fint'] = fint_dict[ifilt][isamp]
    #
    print(binned_ells,flush=True)
    print(hdif_normed,flush=True)
    #
    # The first Ndiscard bins are removed because they have
    # support to very low ell and the last bins because
    # they have support past lmax
    Ndiscard = 1
    binned_ells,hdif_normed  = binned_ells[Ndiscard:-Ndiscard],hdif_normed[Ndiscard:-Ndiscard]
    tt[ifilt]['binned_ells'] = binned_ells
    tt[ifilt]['hdif_normed'] = hdif_normed
    tt[ifilt]['Mbl']         = tt[ifilt]['Mbl'][Ndiscard:-Ndiscard,:]
    #
    print(binned_ells,flush=True)
    print(hdif_normed,flush=True)
#
#
fig,ax = plt.subplots(1,1,figsize=(4.5,3.0))
icol   = 0
for ifilt in flist:
    binned_ells = tt[ifilt]['binned_ells']
    hdif_normed = tt[ifilt]['hdif_normed']
    ax.plot(binned_ells,1e6*hdif_normed,'o:',\
            color='C'+str(icol),mfc='None',label=ifilt)
    ax.axhline(1.5,color='C'+str(icol),ls=':')
    icol = (icol+1)%10
#
mc = np.loadtxt("mc_cosmos_N501_cl.txt")
ax.errorbar(mc[:,0],1e6*mc[:,1],yerr=1e6*mc[:,2],\
            fmt='s',color='C0',mfc='None',label='Mock')
#
ax.legend()
ax.set_xlim(300,1700)
ax.set_ylim(0,9)
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$10^6\ \tilde{C}_\ell$')
#
ax.legend()
plt.tight_layout()
plt.savefig('harmonic_tildeC.pdf')
#
# Now write the key information to a JSON file.
#
outd = {}
for ifilt in flist:
    outd[ifilt] = {}
    outd[ifilt]['fint'] = tt[ifilt]['fint']
    for k in ['binned_ells','hdif_normed','Mbl','wl']:
        outd[ifilt][k] = tt[ifilt][k].tolist()
#
outfn = "harmonic.json"
with open(outfn,"w") as fout:
    json.dump(outd,fout,indent=2)
#
