#!/usr/bin/env python
#
# Use the pre-computed alms for the "data" and "randoms" to
# make a power spectrum and window function (and plot them).
# This code is very slow!
# 
#
import numpy  as np
import json
import matplotlib.pyplot as plt

# Set up the list of filters and samples to process.
flist,slist = ["N419","N501"],["s3","s3"]
# Load the data.
tt = json.load(open("harmonic.json","r"))

#
# Plot the window, normalized to unity at low ell
# and on a linear y-axis.
#
fig,ax = plt.subplots(1,1,figsize=(4.5,3.0))
icol = 0
for ifilt in flist:
    wl   = np.array(tt[ifilt]['wl'])
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
#
fig,ax = plt.subplots(1,1,figsize=(4.5,3.0))
icol   = 0
for ifilt in flist:
    binned_ells = np.array(tt[ifilt]['binned_ells'])
    hdif_normed = np.array(tt[ifilt]['hdif_normed'])
    ax.plot(binned_ells,1e6*hdif_normed,'o:',\
            color='C'+str(icol),mfc='None',label=ifilt)
    #ax.axhline(1.5,color='C'+str(icol),ls=':')
    icol = (icol+1)%10
#
icol = 0
for ifilt in flist:
    mc = np.loadtxt("mc_cosmos_"+ifilt+"_cl.txt")
    ax.errorbar(mc[:,0],1e6*mc[:,1],yerr=1e6*mc[:,2],\
                fmt='s',color='C'+str(icol),mfc='None',label='Mock')
    icol = (icol+1)%10
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
#
#
fig,ax = plt.subplots(1,1,figsize=(4.5,3.0))
icol   = 0
for ifilt in flist:
    mc          = np.loadtxt("mc_cosmos_"+ifilt+"_cl.txt")
    binned_ells = np.array(tt[ifilt]['binned_ells'])
    hdif_normed = np.array(tt[ifilt]['hdif_normed'])
    ax.errorbar(binned_ells,1e6*hdif_normed,yerr=1e6*mc[:,2],\
                fmt='o:',color='C'+str(icol),mfc='None',label=ifilt)
    icol = (icol+1)%10
#
icol = 0
for ifilt in flist:
    mc = np.loadtxt("mc_cosmos_"+ifilt+"_cl.txt")
    ax.plot(mc[:,0],1e6*mc[:,1],'s-',\
            color='C'+str(icol),mfc='None',label='Mock')
    icol = (icol+1)%10
#
ax.legend()
ax.set_xlim(300,1700)
ax.set_ylim(0,9)
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$10^6\ \tilde{C}_\ell$')
#
ax.legend()
plt.tight_layout()
plt.savefig('harmonic_tildeC2.pdf')
#
