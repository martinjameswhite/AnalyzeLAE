#
# A lightweight class to hold distance-redshift and
# growth-of-structure predictions for LCDM cosmologies.
# There is some use in having this standalone code even
# though the functionality is available in CAMB/CLASS or
# in astropy.
#
import numpy as np
from   scipy.special import hyp2f1	# For D(z).



class LCDM:
    # Just some useful distance/redshift functions for LCDM.
    def T_AK(self,x):
        """The "T" function of Adachi & Kasai (2012), used below."""
        b1,b2,b3=2.64086441,0.883044401,0.0531249537
        c1,c2,c3=1.39186078,0.512094674,0.0394382061
        x3   = x**3
        x6   = x3*x3
        x9   = x3*x6
        tmp  = 2+b1*x3+b2*x6+b3*x9
        tmp /= 1+c1*x3+c2*x6+c3*x9
        tmp *= x**0.5
        return(tmp)
        #
    def chi_of_z(self,zz):
        """The comoving distance to redshift zz, in Mpc/h.
           Uses the Pade approximate of Adachi & Kasai (2012) to compute chi
           for a LCDM model, ignoring massive neutrinos."""
        s_ak = (self.OmX/self.OmM)**0.3333333
        tmp  = self.T_AK(s_ak)-self.T_AK(s_ak/(1+zz))
        tmp *= 2997.925/(s_ak*self.OmM)**0.5
        return(tmp)
        #
    def E_of_z(self,zz):
        """The dimensionless Hubble parameter at zz."""
        Ez = (self.OmM*(1+zz)**3 + self.OmX)**0.5
        return(Ez)
        #
    def D_of_z(self,zz):
        """Scale-independent growth factor for flat LCDM."""
        aa = 1./(1.+zz)
        rr = self.OmX/self.OmM
        t1 = hyp2f1(1./3,1,11./6,-aa**3*rr) 
        t2 = hyp2f1(1./3,1,11./6,-rr)
        return( aa * t1/t2 )
        #
    def __init__(self,OmM=0.3):
        self.OmM = OmM
        self.OmX = 1-OmM
