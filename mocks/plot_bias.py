#!/usr/bin/env python3
#
# Plot the (pre-computed) biases as a function of
# scale.
#
import numpy as np
import matplotlib.pyplot as plt




def make_plot(filter_name):
    """Does the work."""
    # Set some filter-specific numbers.
    if filter_name.lower()=="n419":
        yr,bv = [1.6,2.0],1.74
    elif filter_name.lower()=="n501":
        yr,bv = [1.9,2.4],2.035
    else:
        raise RuntimeError("Unknown filter "+filter_name)
    # Load the data.
    bk = np.loadtxt("lae_"+filter_name.lower()+"_bk.txt")
    br = np.loadtxt("lae_"+filter_name.lower()+"_br.txt")
    # Make the figure.
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(6,2.50))
    # First do b(k).
    ax[0].plot(bk[1:,0],bk[1:,1],'d',color='C0',mfc='None',label='$b_a$')
    ax[0].plot(bk[1:,0],bk[1:,2],'s',color='C1',mfc='None',label=r'$b_\times$')
    ax[0].axhline(bv,ls=':',color='k')
    ax[0].legend(loc=2)
    ax[0].set_xlim(0.0,0.5)
    ax[0].set_ylim(yr)
    ax[0].set_xlabel(r'$k\quad [h\,{\rm Mpc}^{-1}]$')
    ax[0].set_ylabel(r'Bias')
    # then do b(r)
    ax[1].plot(br[:,0],br[:,1],'d',color='C0',mfc='None',label='$b_a$')
    ax[1].plot(br[:,0],br[:,2],'s',color='C1',mfc='None',label=r'$b_\times$')
    ax[1].axhline(bv,ls=':',color='k')
    ax[1].legend(loc=1)
    ax[1].set_xlim(0.5,25)
    ax[1].set_xscale('log')
    ax[1].set_ylim(yr)
    ax[1].set_xlabel(r'$r\quad [h^{-1}\,{\rm Mpc}]$')
    #
    plt.tight_layout()
    plt.savefig('odin_'+filter_name.upper()+'_bias.png')
    plt.savefig('odin_'+filter_name.upper()+'_bias.pdf')
    #



if __name__=="__main__":
    for fn in ["N419","N501"]:
        make_plot(fn)
