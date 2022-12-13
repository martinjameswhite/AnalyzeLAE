import numpy as np
import json
#
from cobaya.theory        import Theory
from cobaya.likelihood    import Likelihood
from scipy.interpolate    import InterpolatedUnivariateSpline as Spline
#


# Class for angular power spectrum likelihood.
class FullShapeLikelihood(Likelihood):
    zfid:      float
    Omfid:     float
    #
    basedir:   str
    fs_datfn:  str
    fs_covfn:  str
    fs_winfn:  str
    pklin_fn:  str
    #
    fs_lmin:   float
    fs_lmax:   float
    #
    def initialize(self):
        """Sets up the class."""
        # Set up the power spectrum class.
        plin = np.loadtxt(self.pklin_fn)
        self.pofk = LinearPowerSpectra(plin[:,0],plin[:,1])
        # and load the data, covariance, window functions, etc.
        self.loadData()
        #
    def get_requirements(self):
        # We will be working at fixed cosmology.
        req = {'b1': None,\
               'b2': None,\
               'bs': None,\
               'alpha0': None,\
               'alpha2': None,\
               'SN0': None,\
               'SN2': None }
        return(req)
        #
    def logp(self,**params_values):
        """Return a log-likelihood."""
        cl_thy  = self.cl_predict()
        cl_obs  = self.cl_observe(cl_thy)
        diff    = self.cls - cl_obs
        chi2    = np.dot(diff,np.dot(self.cinv,diff))
        self.thy= cl_thy
        self.obs= cl_obs
        return(-0.5*chi2)
        #
    def loadData(self):
        """
        Loads the required data.
        """
        # First load the data
        cl_dat    = np.loadtxt(self.basedir+self.fs_datfn)
        self.ells = cl_dat[:,0]
        self.cls  = cl_dat[:,1]
        # Now load the covariance matrix.
        cov = np.loadtxt(self.basedir+self.fs_covfn)
        # Finally load the window function matrix.
        self.wla = np.loadtxt(self.basedir+self.fs_winfn)
        #
        # We're only going to want some of the entries in computing chi^2,
        # handle this by adjusting Cov.
        lcut = (self.ells > self.fs_lmax)\
             | (self.ells < self.fs_lmin)
        for i in np.nonzero(kcut)[0]:
            cov[i,:] = 0
            cov[:,i] = 0
            cov[i,i] = 1e25
        # Copy cov and save the inverse.
        self.cov  = cov
        self.cinv = np.linalg.inv(self.cov)
        #
    def cl_predict(self):
        """Use the model to compute C_ell, given biases etc."""
        pp   = self.provider
        b1   = pp.get_param('b1')
        b2   = pp.get_param('b2')
        bs   = pp.get_param('bs')
        alp0 = pp.get_param('alpha0')
        alp2 = pp.get_param('alpha2')
        sn0  = pp.get_param('SN0')
        sn2  = pp.get_param('SN2')
        #
        tt   = np.zeros_like(self.ells)
        return(tt)
        #
    def cl_observe(self,tt):
        """Apply the window function matrix to get the binned prediction."""
        # Convolve with window.
        return(np.dot(self.wla,tt))