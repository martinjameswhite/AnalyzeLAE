#!/usr/bin/env python
#
# Make a plot of dN/dz, along with a fit
# and label the interloper fraction(s).
#
# Also write the selection function fit.
#
import numpy as np
import matplotlib.pyplot as plt


def plot_filter(filt_name,fname="cosmos"):
    """Make a figure for filter filt_name in field fname."""
    # We will need a random number generator.
    rng   = np.random.default_rng()
    # Set the ranges based on the filter.
    if filt_name=="N419": z0,dz,chi0,dcdz = 2.4,0.03,3941.,829.
    ###if filt_name=="N501": z0,dz,chi0,dcdz = 3.1,0.03,4448.,633.
    if filt_name=="N501": z0,dz,chi0,dcdz = 3.123,0.03,4462.,627.9
    if filt_name=="N673": z0,dz,chi0,dcdz = 4.5,0.04,5160.,411.
    # Set some limits, choosing round numbers.
    zmin  = 1e-2*int( 100*(z0-2*dz)+0 )
    zmax  = 1e-2*int( 100*(z0+2*dz)+1 )
    # Redshifts with VI_QUALITY >= 3 should be considered “good”
    zvals = np.loadtxt("VI_Redshifts_{:s}.txt".format(filt_name))
    zvals = zvals[zvals[:,1]>=3,0]
    print("Median redshift for ",filt_name," is ",np.median(zvals))
    # Work out the interloper fraction.
    ww   = np.nonzero( (zvals<zmin)|(zvals>zmax) )[0]
    fint = float(len(ww))/float(len(zvals))
    print("Interloper fraction {:.3f}".format(fint))
    # Fit to dN/dz histogram.
    zarr = np.arange(zmin,zmax,0.0025)
    carr = chi0+dcdz*(zarr-z0)
    pchi = np.exp(-np.abs( (zarr-z0)/(dz) )**6 )
    # Write the selection function to file.
    with open("lae_"+fname+"_sfn_"+filt_name+".txt","w") as fout:
        fout.write("# Radial selection function.\n")
        fout.write("# "+fname+", "+filt_name+"\n")
        fout.write("# {:>8s} {:>12s} {:>12s}\n".format("z","chi[Mpc/h]","pchi"))
        for i in range(zarr.size):
            fout.write("{:10.5f} {:12.2f} {:12.6f}\n".\
                       format(zarr[i],carr[i],pchi[i]))
    # Normalize to unit integral for comparison with the normed counts.
    fsamp = pchi / np.trapz(pchi,x=zarr)
    # Now make the figure.
    fig,ax= plt.subplots(1,1,figsize=(6,3))
    # A (normalized) histogram of the counts.
    Nbins = 30
    bins  = np.linspace(zmin,zmax,Nbins+1)
    ax.hist(zvals,bins=bins,density=True,color='C0',label='VI')
    # and the analytic fit.
    ax.plot(zarr,fsamp,'-',color='C1',label='Fit')
    # legends and labels.
    ax.text(zmin+0.01,np.max(fsamp),r'$f_{int}='+'{:.3f}'.format(fint)+'$')
    ax.legend(title=filt_name)
    ax.set_xlim(zmin,zmax)
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$d\ln N/dz$')
    # Upper x-axis.
    cmin,cmax = chi0+dcdz*(zmin-z0),chi0+dcdz*(zmax-z0)
    ax2 = ax.twiny()
    ax2.set_xlim(cmin,cmax)
    ax2.set_xlabel(r'$\chi\quad [h^{-1}{\rm Mpc}]$')
    #
    plt.tight_layout()
    plt.savefig('odin_'+filt_name+"_dndz.png")
    #


if __name__=="__main__":
    for filt in ["N501"]:
        plot_filter(filt)
