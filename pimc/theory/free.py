import numpy as np
from math import *
from scipy.special import zeta
import mpmath
from scipy.optimize import brentq

def G(p,z):
    return float(mpmath.polylog(p, z))


def fugacity(nLanda3,zmin=0,zmax=0.999999999):
    ncLanda3=G(3/2.,1)
    if nLanda3 < ncLanda3:
        
        rootFunction = lambda x :  ( G(3/2,x) - nLanda3 )/ncLanda3
        return brentq(rootFunction,zmin,zmax)
    else:
        return 1


def criticalTemperature(n):
    return 2*pi*(n/G(3/2.,1))**(2./3)


def energySingleComponent(n , T ):
    '''
    Returns the energy per particle of a single component
    '''
    beta=1/T
    landa=sqrt(2*pi*beta)
    Tc =criticalTemperature(n)
    z=fugacity(n*landa**3)
    if n==0:
        return 0
    return 3*T*G(5./2,z)/(2*n * landa**3 )

def energyMixture(n1,n2, T ):
    
    n=n1+n2
    
    return (n1/n)*energy(n1,T) + (n2/n)*energy(n2,T)


def energy( n, T ):
    if hasattr(n, '__iter__'):
        n1,n2 = n
        return energyMixture(n1,n2)
    else:
        return energySingleComponent(n,T)




def compressibility(n,T):
    landa=np.sqrt(2*pi/T)
    beta=1/T
    nLanda3=n*landa**3
    z=fugacity(nLanda3,zmin=0,zmax=0.999999999)
    return G(1/2,z)*beta/landa**3


energy=np.vectorize(energy)
compressibility=np.vectorize(compressibility)
