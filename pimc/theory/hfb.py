import numpy as np
from math import *
from scipy.special import zeta
import mpmath
from scipy.optimize import brentq
import matplotlib.pylab as pl
from scipy.integrate import quad



def G(p,z):
    return float(mpmath.polylog(p, z))


def findZ(n,a=1e-12,b=1-1e-12):
      return brentq(lambda z : n - G(3/2,z) ,a,b)


def freeEnergy(n1,n2,T,g,g12):
    beta = 1/T
    xi=G(3/2,1)
    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0_1 = n1 - xi/landa**3
    n0_2 = n2 - xi/landa**3

    y1 = np.exp( - beta * g * n0_1)
    y2 = np.exp( - beta * g * n0_2)

    if (n1*landa**3 >= xi) and (n2*landa**3 >= xi):
        F = g/2*(n1**2 + n2**2) +g12*n1*n2 + g*xi**2/landa**6 -1/(landa**3*beta) *( G(5/2,y1) + G(5/2,y2) )
    else:
        if (n1*landa**3> xi ):
            z=findZ(n2*landa**3)
            mu=1/beta * np.log(z)
            F= g/2*(n1**2 + 2*n2**2 + xi**2/landa**6) + g12* n1 * n2 + mu*n2 -1/(landa**3*beta) *( G(5/2,y1) + G(5/2,z) )
        else:
            F=None
    


    return F

freeEnergy=np.vectorize(freeEnergy)


def energy(n1,n2,T,g,g12):
    beta = 1/T
    xi=G(3/2,1)
    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0_1 = n1 - xi/landa**3
    n0_2 = n2 - xi/landa**3

    y1 = np.exp( - beta * g * n0_1)
    y2 = np.exp( - beta * g * n0_2)

    U = 0.5 * g * ( n1**2 + n2**2 ) + g12*n1*n2 -2*g*xi**2/landa**6 + (  3/(2 * beta * landa**3) * ( G(5/2,y1 ) + G(5/2,y2)    ) +  (g * xi / landa**6 * 3/2   )* (G(3/2,y1) + G(3/2,y2)  )   + g/landa**3 * (n1*G(3/2,y1) +n2*G(3/2,y2)  )  ) 

    return U


def inverseSusceptibility(n1,n2,T,g,g12):
    xi=G(3/2,1)

    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0_1 = n1 - xi/landa**3
    n0_2 = n2 - xi/landa**3
    
    y1 = np.exp( -  beta * g * n0_1 )
    y2 = np.exp(  -  beta * g * n0_2)


    if  (n0_1 > 0 and n0_2 > 0 ):
        chi = ( g - g12)/2  - beta*g**2/(4*landa**3) * (  G(1/2,y1) + G(1/2,y2)     )
    else:
        chi=None
    
    return 2*chi


inverseSusceptibility=np.vectorize(inverseSusceptibility)



def criticalRatio(n,g,T):

    return brentq(lambda r : inverseSusceptibility(n1=n/2,n2=n/2,T=T,g=g,g12=r*g) ,0.001,2)


criticalRatio = np.vectorize(criticalRatio)