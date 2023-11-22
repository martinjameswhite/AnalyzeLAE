#!/usr/bin/env python3
#
# Plot the (pre-computed) biases as a function of
# scale.
#
import numpy as np
import matplotlib.pyplot as plt




def make_plot():
    """Does the work."""
    # Load the data.
    bk = np.loadtxt("lae_n501_bk.txt")
    br = np.loadtxt("lae_n501_br.txt")
    # Make the figure.
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(6,2.50))
    # First do b(k).
    ax[0].plot(bk[1:,0],bk[1:,1],'d',color='C0',mfc='None',label='$b_a$')
    ax[0].plot(bk[1:,0],bk[1:,2],'s',color='C1',mfc='None',label=r'$b_\times$')
    ax[0].axhline(2.0,ls=':',color='k')
    ax[0].legend(loc=2)
    ax[0].set_xlim(0.0,0.5)
    ax[0].set_ylim(1.5,2.5)
    ax[0].set_xlabel(r'$k\quad [h\,{\rm Mpc}^{-1}]$')
    ax[0].set_ylabel(r'Bias')
    # then do b(r)
    ax[1].plot(br[:,0],br[:,1],'d',color='C0',mfc='None',label='$b_a$')
    ax[1].plot(br[:,0],br[:,2],'s',color='C1',mfc='None',label=r'$b_\times$')
    ax[1].axhline(2.0,ls=':',color='k')
    ax[1].legend(loc=1)
    ax[1].set_xlim(0.5,25)
    ax[1].set_xscale('log')
    ax[1].set_ylim(1.5,2.5)
    ax[1].set_xlabel(r'$r\quad [h^{-1}\,{\rm Mpc}]$')
    #
    plt.tight_layout()
    plt.savefig('odin_N501_bias.png')
    plt.savefig('odin_N501_bias.pdf')
    #



if __name__=="__main__":
    make_plot()
