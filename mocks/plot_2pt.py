#!/usr/bin/env python3
#
# Plot the (pre-computed) real-space 2-point functions
#
import numpy as np
import matplotlib.pyplot as plt




def make_xir_plot(filt_name,zeff):
    """Does the work."""
    # Load the data.
    br   = np.loadtxt("lae_"+filt_name+"_br.txt")
    ximm = br[:,3]
    xigg = br[:,1]**2 * ximm
    xigm = br[:,2]    * ximm
    # Make the figure.
    fig,ax = plt.subplots(1,1,figsize=(4.5,3))
    #
    ax.plot(br[:,0],xigg,'d',color='C0',mfc='None',label=r'$\xi_{gg}$')
    ax.plot(br[:,0],xigm,'s',color='C1',mfc='None',label=r'$\xi_{gm}$')
    ax.plot(br[:,0],ximm,'o',color='C2',mfc='None',label=r'$\xi_{mm}$')
    tstring = filt_name.upper()+r', $z\simeq '+'{:.1f}$'.format(zeff)
    ax.legend(title=tstring,loc=1)
    ax.set_xlim(0.5,25)
    ax.set_ylim(1e-2,25)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r\quad [h^{-1}\,{\rm Mpc}]$')
    ax.set_ylabel(r'$\xi(r)$')
    #
    plt.tight_layout()
    plt.savefig('odin_'+filt_name.upper()+'_xir.png')
    plt.savefig('odin_'+filt_name.upper()+'_xir.pdf')
    #




def make_pkr_plot(filt_name,zeff):
    """Does the work."""
    # Load the data.
    bk  = np.loadtxt("lae_"+filt_name+"_bk.txt")
    Pmm = bk[:,3]
    Pgg = bk[:,1]**2 * Pmm
    Pgm = bk[:,2]    * Pmm
    # Make the figure.
    fig,ax = plt.subplots(1,1,figsize=(4.5,3))
    #
    ax.plot(bk[:,0],Pgg,'d',color='C0',mfc='None',label=r'$P_{gg}$')
    ax.plot(bk[:,0],Pgm,'s',color='C1',mfc='None',label=r'$P_{gm}$')
    ax.plot(bk[:,0],Pmm,'o',color='C2',mfc='None',label=r'$P_{mm}$')
    tstring = filt_name.upper()+r', $z\simeq '+'{:.1f}$'.format(zeff)
    ax.legend(title=tstring,loc=1)
    ax.set_xlim(0.0,0.5)
    ax.set_ylim(20,2e4)
    ax.set_yscale('log')
    ax.set_xlabel(r'$k\quad [h\,{\rm Mpc}^{-1}]$')
    ax.set_ylabel(r'$P(k)\quad [h^{-3}{\rm Mpc}^3]$')
    #
    plt.tight_layout()
    plt.savefig('odin_'+filt_name.upper()+'_pkr.png')
    plt.savefig('odin_'+filt_name.upper()+'_pkr.pdf')
    #




if __name__=="__main__":
    for filt_name,zeff in zip(["n419","n501"],[2.4,3.1]):
        make_xir_plot(filt_name,zeff)
        make_pkr_plot(filt_name,zeff)
