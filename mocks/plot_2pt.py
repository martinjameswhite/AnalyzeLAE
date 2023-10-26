#!/usr/bin/env python3
#
# Plot the (pre-computed) real-space 2-point functions
#
import numpy as np
import matplotlib.pyplot as plt




def make_xir_plot():
    """Does the work."""
    # Load the data.
    br   = np.loadtxt("lae_br.txt")
    ximm = br[:,3]
    xigg = br[:,1]**2 * ximm
    xigm = br[:,2]    * ximm
    # Make the figure.
    fig,ax = plt.subplots(1,1,figsize=(6,3.25))
    #
    ax.plot(br[:,0],xigg,'d',color='C0',mfc='None',label=r'$\xi_{gg}$')
    ax.plot(br[:,0],xigm,'s',color='C1',mfc='None',label=r'$\xi_{gm}$')
    ax.plot(br[:,0],ximm,'o',color='C2',mfc='None',label=r'$\xi_{mm}$')
    ax.legend(title=r'N501, $z\simeq 3.1$',loc=1)
    ax.set_xlim(0.5,25)
    ax.set_ylim(1e-2,25)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r\quad [h^{-1}\,{\rm Mpc}]$')
    ax.set_ylabel(r'$\xi(r)$')
    #
    plt.tight_layout()
    plt.savefig('odin_N501_xir.png')
    plt.savefig('odin_N501_xir.pdf')
    #




def make_pkr_plot():
    """Does the work."""
    # Load the data.
    bk  = np.loadtxt("lae_bk.txt")
    Pmm = bk[:,3]
    Pgg = bk[:,1]**2 * Pmm
    Pgm = bk[:,2]    * Pmm
    # Make the figure.
    fig,ax = plt.subplots(1,1,figsize=(6,3.25))
    #
    ax.plot(bk[:,0],Pgg,'d',color='C0',mfc='None',label=r'$P_{gg}$')
    ax.plot(bk[:,0],Pgm,'s',color='C1',mfc='None',label=r'$P_{gm}$')
    ax.plot(bk[:,0],Pmm,'o',color='C2',mfc='None',label=r'$P_{mm}$')
    ax.legend(title=r'N501, $z\simeq 3.1$',loc=1)
    ax.set_xlim(0.0,0.5)
    ax.set_ylim(20,2e4)
    ax.set_yscale('log')
    ax.set_xlabel(r'$k\quad [h\,{\rm Mpc}^{-1}]$')
    ax.set_ylabel(r'$P(k)\quad [h^{-3}{\rm Mpc}^3]$')
    #
    plt.tight_layout()
    plt.savefig('odin_N501_pkr.png')
    plt.savefig('odin_N501_pkr.pdf')
    #




if __name__=="__main__":
    make_xir_plot()
    make_pkr_plot()
