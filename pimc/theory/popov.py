from zipfile import ZIP_BZIP2
import numpy as np
from math import *
from scipy.special import zeta
import mpmath
from scipy.optimize import brentq
import matplotlib.pylab as plt
from scipy.integrate import quad


def H(y,a=1e-8,b=30,plot=False):
    
    f= lambda x : x**(3/2)/(np.exp(np.sqrt(x**2 + 2*y*x)) - 1) *(x+y)/np.sqrt(x**2 + 2*y*x) 
    
    if plot:
        x=np.linspace(a,b,num=1000)
        plt.plot(x,f(x),label=y)
    return 4/(3*np.sqrt(pi)) * quad(f,a,b,epsabs=1e-9)[0]



H=np.vectorize(H)

def G(p,z):
    return float(mpmath.polylog(p, z))


def findZ(n,a=1e-12,b=1-1e-12):
      return brentq(lambda z : n - G(3/2,z) ,a,b)



def freeEnergy(n1,n2,T,g,g12,plotH=False):
    # returns the free energy per unit volume of a mixture
    xi=G(3/2,1)
    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0_1 = n1 - xi/landa**3
    n0_2 = n2 - xi/landa**3

    n0 = n0_1 + n0_2
    m0 = n0_1 -  n0_2

    delta =np.sqrt(  (g**2 - g12**2 )*m0**2 + (g12*n0)**2)
    lambdaPlus = 0.5*(g*n0 + delta)
    lambdaMinus = 0.5*(g*n0 - delta)
    
    yPlus=beta*lambdaPlus
    yMinus=beta*lambdaMinus

    if (n1*landa**3 >= xi) and (n2*landa**3 >= xi):
        F = 0.5* g * (n1**2 +  n2**2 ) + g12 * n1 * n2 + g * xi**2/landa**6 - 1/(beta*landa**3)*(H(yPlus,plot=plotH) + H(yMinus,plot=plotH) )
        F+= 1/(2*pi)**(3/2)*4/(15*sqrt(pi) ) *(    2**(5/2)* (lambdaMinus**(5/2) + lambdaPlus**(5/2))  )

    else:

        if (n1*landa**3> xi ):
            z=findZ(n2*landa**3)
            mu=1/beta * np.log(z)
            F= g/2*(n1**2 + 2*n2**2 + xi**2/landa**6) + g12* n1 * n2 + mu*n2 + 1/(2*pi)**(3/2) * 4/(15*np.sqrt(pi)) *(2*g*n0_1)**(5/2)  - 1/(beta*landa**3)*H( beta*g*n0_1,plot=plotH) -1/(landa**3*beta) * G(5/2,z)
        
        else:
            if  ( (n2*landa**3 <  xi ) ):
                z2=findZ(n2*landa**3)
                z1=findZ(n1*landa**3)

                mu2=1/beta * np.log(z2)
                mu1=1/beta * np.log(z1)

                F= g*(n1**2 + n2**2) + g12*(n1*n2)   + mu1*n1 + mu2*n2 -1/(landa**3*beta) * G(5/2,z1) -1/(landa**3*beta) * G(5/2,z2) 
                

                


    
    return F



freeEnergy=np.vectorize(freeEnergy)

