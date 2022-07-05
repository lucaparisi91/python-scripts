import numpy as np
import scipy as sp
from scipy import interpolate
import pandas as pd
from math import *
from scipy.optimize import brentq


class scatteringLengthGaussianPotential:
    def __init__(self,filename="/home/luca/data/droplet-finite-temperature/scatteringLength/extrapolatedScatteringLengthsGaussianPotential.dat"):
        self.scatteringLengthData=pd.read_csv(filename,delim_whitespace=True)
        
        self.interp=interpolate.interp1d(self.scatteringLengthData["V0"],self.scatteringLengthData["a"])
        
        self.minV0=np.min(self.scatteringLengthData["V0"])
        
    def aContact(self,V0,R0):
        return (V0*(2*pi*R0**2)**(3/2.) )/(4*pi)
    
    def __call__(self,V0,R0):
        V0=np.array(V0)
        R0=np.array(R0)
        res=V0 * 0
        p=V0*R0**2
        tooSmallMask=p <= self.minV0
        interpMask = p >self.minV0
        res[tooSmallMask]=self.aContact(V0[tooSmallMask],R0[tooSmallMask])
        res[interpMask]=self.interp(p[interpMask])*R0[interpMask]
        
        return res




def V0(na3,nR03,contact=False,minSearch=1e-7,maxSearch=1e+2):
    '''
    Return the amplitude of the gaussian potential in units where the the density is one and \hbar=m=1
    '''
    if contact:
        return 4*pi*na3**(1./3)/((2*pi)**(3/2.)*nR03)
    else:
        def V0FromGaussianSim(na3,nR03):
            ag=scatteringLengthGaussianPotential()
            R0=nR03**(1/3.)
            aTarget=na3**(1/3.)
            def fRoot(V0):
                return ag([V0],[ R0])[0] - aTarget
        
            #return ag([maxSearch],[R0])
            return brentq(fRoot,minSearch,maxSearch)
        V0s=[V0FromGaussianSim(na3Float,nR03Float)  for na3Float, nR03Float in zip(na3,nR03) ]
        return V0s




def V0_R0( a  ,contact=False,minSearch=1e-7,maxSearch=1e+1):
    '''
    Return the amplitude of the gaussian potential in units where R0=1 the is one and \hbar=m=1
    '''
    
    if contact:
        return 4*pi*a/( (2*pi)**(3/2.) )
    else:
        ag=scatteringLengthGaussianPotential()
        R0=1

        def fRoot( V0 ):
            return ag( [V0],[ R0])[0] - a
        
        return brentq(fRoot,minSearch,maxSearch)

V0_R0 = np.vectorize(V0_R0)


def energyT0Weak(na3):
    "In units of density=1 energy per particle of the weakly interacting ground state"
    
    return 0.5 * 4*pi*na3**(1./3)