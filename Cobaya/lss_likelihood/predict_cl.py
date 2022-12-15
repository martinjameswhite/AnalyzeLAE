#
# Code to compute angular power spectra for a "thin shell" survey
#
import numpy       as np
import sys

from lcdm              import LCDM
from scipy.integrate   import simps as simpson
from scipy.special     import legendre
from scipy.special     import hyp2f1	# For D(z).
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

# If we want to model the (redshift-space) power spectrum using LPT.
from velocileptors.LPT.lpt_rsd_fftw  import LPT_RSD









class LinearPowerSpectra():
    """Computes the (redshift-space) power spectrum, P(k,mu) [Mpc/h units]."""
    def __init__(self,klin,plin):
        """klin,plin: Arrays containing Plinear [Mpc/h units]."""
        # Set up a spline for Plin(k).
        self.plin = Spline(klin,plin)
    def __call__(self,kk,mu,pars):
        """Computes the three P(k,z).  Returns k,Pgg,Pgm,Pmm [Mpc/h units]."""
        bb,sigma,sn = pars
        ff  = 1.0
        plin= self.plin
        pgg = (bb+ff*mu**2)**2*plin(kk)*np.exp(-0.5*(kk*mu*sigma)**2)+sn
        pgm = (bb+ff*mu**2)   *plin(kk)
        pmm = plin(kk)
        return((kk,pgg,pgm,pmm))
        #






class LPTPowerSpectra():
    """Computes the (redshift-space) power spectrum, P(k,mu) [Mpc/h units]."""
    def __init__(self,klin,plin):
        """klin,plin: Arrays containing Plinear [Mpc/h units]."""
        # Copy the arguments.
        self.klin = klin
        self.plin = plin
        # Set up the LPT class -- this can take a little while.
        self.lpt = LPT_RSD(klin,plin,kIR=0.2)
        # biases = [0.71,0.26,0.67,0.52]
        # cterms = [-3.4,-1.7,6.5,0]
        # stoch  = [1500.,-1900.,0]
        # pars   = biases + cterms + stoch
        self.lpt.make_pltable(f=1,kmin =5e-3,kmax=1.0,nk=60,nmax=4,apar=1,aperp=1)
    def __call__(self,kk,mu,pars):
        """Computes the three P(k,mu).  Returns k,Pgg,Pgm,Pmm [Mpc/h units]."""
        # We will approximate P(k,mu) through the multipole expansion.
        kl,p0,p2,p4 = self.lpt.combine_bias_terms_pkell(pars)
        p0  = Spline(kl,p0)
        p2  = Spline(kl,p2)
        p4  = Spline(kl,p4)
        # Return only P_{gg}, set the others to zero.
        pgg = p0(kk)*1.0 + p2(kk)*legendre(2)(mu) + p4(kk)*legendre(4)(mu)
        pgm = np.zeros_like(kk)
        pmm = np.zeros_like(kk)
        return((kk,pgg,pgm,pmm))
        #








class RedshiftDistribution:
    # Holds information about the redshift distribution.
    def __init__(self,OmM,dndz):
        cc   = LCDM(OmM)
        Nz   = 256
        Nchi = 2048
        zmin = dndz[ 0,0]
        zmax = dndz[-1,0]
        zz   = np.linspace(zmin,zmax,Nz)
        dndz = Spline(dndz[:,0],dndz[:,1])(zz)
        # Normalize dN/dz.
        dndz = dndz/simpson(dndz,x=zz)
        # Set up the chi(z) array and z(chi) spline.
        chiz = np.array([cc.chi_of_z(z) for z in zz])
        zchi = Spline(chiz,zz)
        # Work out W(chi) for the objects whose dNdz is supplied.
        chimin= np.min(chiz) + 1e-5
        chimax= np.max(chiz)
        chival= np.linspace(chimin,chimax,Nchi)
        zval  = zchi(chival)
        fchi  = Spline(zz,dndz*cc.E_of_z(zz))(zval)
        fchi /= simpson(fchi,x=chival)
        # Compute the effective redshift.
        zeff = simpson(zval*fchi**2/chival**2,x=chival)
        zeff/= simpson(     fchi**2/chival**2,x=chival)
        # Now store the things we want to keep.
        self.chis = chival
        self.zval = zval
        self.fchi = fchi
        self.zeff = zeff
    def fwhm(self):
        """Returns the FWHM of fchi."""
        fmax = np.max(self.fchi)
        i=0
        while self.fchi[i]<0.5*fmax: i+=1
        j=len(self.fchi)-1
        while self.fchi[j]<0.5*fmax: j-=1
        return(self.chis[j]-self.chis[i])
    #







class FlatSkyCl:
    # Computes the flat-sky expression for C_l.
    def __init__(self,OmM,pofk,dndz):
        """Initialize the class."""
        self.cc   = LCDM(OmM)
        self.phiz = RedshiftDistribution(OmM,dndz)
        #print("# zeff=",self.phiz.zeff)
        # Use a simple approximation to chi0, assuming this
        # is all for thin-shell computations.
        self.chi0 = np.mean(self.phiz.chis)
        # Compute the kpar window (needs chi0).  This should
        # work out kmax and Nk dynamically.  For now Nk is fixed.
        kmax = 25./self.phiz.fwhm()
        self.kpar= np.linspace(0,kmax,128)
        self.wp2 = self.window_par_brute(self.kpar)
        # Make a power-spectrum instance.
        # We want the effective redshift to scale Plin.
        Dz = self.cc.D_of_z(self.phiz.zeff)
        #self.pofk = LinearPowerSpectra(pofk[:,0],pofk[:,1]*Dz**2)
        self.pofk = LPTPowerSpectra(pofk[:,0],pofk[:,1]*Dz**2)
    def pkmu(self,kk,mu,pars):
        kv,pgg,pgm,pmm = self.pofk(kk,mu,pars)
        return(pgg)
    def window_par_brute(self,kpar):
        """Compute W_parallel(k_parallel) by brute force."""
        # Get the comoving distance grid and f(chi).
        chival = self.phiz.chis.copy() - self.chi0
        fchi   = self.phiz.fchi.copy()
        # Now brute force the Fourier transform.
        wpar   = np.zeros_like(kpar)
        for i,kk in enumerate(kpar):
            cosx    = simpson( fchi*np.cos(kk*chival),x=chival )
            sinx    = simpson( fchi*np.sin(kk*chival),x=chival )
            wpar[i] = cosx**2 + sinx**2
        return(wpar)
    def window_par_matrix(self,kpar):
        """Compute W_parallel^ab(k_parallel) by brute force."""
        # This is the generalization to cross-spectra, but it
        # requires we have multiple dN/dz's and hence phi(z)'s.
        # Brute force the Fourier transform.
        Nshell = len(self.phiz)
        wpar   = np.zeros( (len(kpar),Nshell,Nshell) )
        for i,kk in enumerate(kpar):
            for a in range(Nshell):
                phiz = self.phiz[a]
                chis = phiz.chis.copy() - self.chi0
                fchi = phiz.fchi.copy()
                cosa = simpson( fchi*np.cos(kk*chis),x=chis )
                sina = simpson( fchi*np.sin(kk*chis),x=chis )
                for b in range(a,Nshell):
                    phiz = self.phiz[b]
                    chis = phiz.chis.copy() - self.chi0
                    fchi = phiz.fchi.copy()
                    cosb = simpson( fchi*np.cos(kk*chis),x=chis )
                    sinb = simpson( fchi*np.sin(kk*chis),x=chis )
                    wpar[i,a,b] = cosa*cosb + sina*sinb
                    wpar[i,b,a] = wpar[i,a,b]
        return(wpar)
    def __call__(self,ells,pars):
        """Computes the flat-sky C_ell"""
        kpar = self.kpar
        wp2  = self.wp2
        cl   = np.zeros(len(ells),dtype='float')
        for i,ell in enumerate(ells):
            kperp = (ell+0.5)/self.chi0 # Should be a non-zero scalar!
            kval  = np.sqrt(kperp**2 + kpar**2)
            uval  = kpar/kval
            pku   = self.pkmu(kval,uval,pars)
            cl[i] = simpson(pku*wp2,x=kpar)
        cl /= np.pi * self.chi0**2
        return(cl)
    #
