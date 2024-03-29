#
# Code to compute angular power spectra for a "thin shell" survey
#
import numpy as np
import json
import sys


from scipy.integrate   import simps as simpson
from scipy.special     import legendre
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

# If we want to model the (redshift-space) power spectrum using LPT.
from velocileptors.LPT.lpt_rsd_fftw  import LPT_RSD







class LPTCorrelationFunctions():
    """Computes the (redshift-space) correlation function, xi(s,mu)."""
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
        self.lpt.make_pltable(f=1,kmin=5e-3,kmax=1.0,nk=60,nmax=4,apar=1,aperp=1)
    def __call__(self,ss,mu,pars):
        """Computes the three xi(s,mu).  Returns s,xi_gg,xi_gm,xi_mm."""
        # We will approximate xi(s,mu) through the multipole expansion.
        x0,x2,x4 = self.lpt.combine_bias_terms_xiell(pars)
        x0  = Spline(x0[0],x0[1])
        x2  = Spline(x2[0],x2[1])
        x4  = Spline(x4[0],x4[1])
        # Return only xi_{gg}, set the others to zero.
        xgg = x0(ss)*1.0 + x2(ss)*legendre(2)(mu) + x4(ss)*legendre(4)(mu)
        xgm = np.zeros_like(ss)
        xmm = np.zeros_like(ss)
        return((xgg,xgm,xmm))
        #



class NbodyCorrelationFunctions():
    """Returns N-body derived xi_ell(s) from pre-computed tables."""
    def __init__(self,json_file):
        """Reads the tables from a JSON file."""
        self.mod = json.load(open(json_file,"r"))
    def __call__(self,sval,uval,pars):
        """Reads the three xi(s,mu).  Returns s,xi_gg,xi_gm,xi_mm.
        In this case 'pars' is an integer index, giving the model number."""
        # We will approximate xi(s,mu) through the multipole expansion.
        mock= self.mod['mocks'][pars]
        ss  = np.array(self.mod['R'])
        x0  = Spline(ss,np.array(mock['xi0']))
        x2  = Spline(ss,np.array(mock['xi2']))
        x4  = Spline(ss,np.array(mock['xi4']))
        # Return only xi_{gg}, set the others to zero.
        xgg = x0(sval)*1.0 + x2(sval)*legendre(2)(uval) + x4(sval)*legendre(4)(uval)
        xgm = np.zeros_like(sval)
        xmm = np.zeros_like(sval)
        return((xgg,xgm,xmm))
        #










class ThinShellWR:
    # Computes the thin-shell expression for w_theta(R).
    #
    def __init__(self,XiModel,chibar,sfn_fname):
        """Initialize the class."""
        # Store the midpoint.
        self.chi0 = chibar
        # Load the selection function (z, chi, pchi) and
        # extract the fields.
        sfn  = np.loadtxt(sfn_fname)
        chi  = sfn[:,1]
        phi  = sfn[:,2]
        cmin = np.min(chi)
        cmax = np.max(chi)
        # Normalize phi(chi).
        phi /= simpson(phi,x=chi)
        # Make a window function.
        self.yy   = np.linspace(0.0,cmax-cmin,251)
        self.wind = np.zeros_like(self.yy)
        chim      = np.linspace(cmin,cmax,401)
        for i in range(self.yy.size):
            y2           = self.yy[i]/2.0
            win1         = np.interp(chim-y2,chi,phi,left=0,right=0)
            win2         = np.interp(chim+y2,chi,phi,left=0,right=0)
            self.wind[i] = simpson(win1*win2,x=chim)
        ## Analytic: self.wind = (self.delt-self.yy)/self.delt**2
        # Copy the correlation function model,
        # e.g. LPTCorrelationFunction(klin,plin*Dz**2)
        self.xiofs = XiModel
    def xismu(self,ss,mu,pars):
        xgg,xgm,xmm = self.xiofs(ss,mu,pars)
        return(xgg)
    def __call__(self,Rs,pars):
        """Computes the thin-shell w_theta(R)"""
        wR = np.zeros(len(Rs),dtype='float')
        for i,RR in enumerate(Rs):
            sval  = np.sqrt(RR**2 + self.yy**2)
            uval  = self.yy/sval
            xisu  = self.xismu(sval,uval,pars)
            wR[i] = 2*simpson(self.wind*xisu,x=self.yy)
        return(wR)
    #
