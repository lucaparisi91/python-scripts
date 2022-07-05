import numpy as np
from math import *
from scipy.special import zeta
import mpmath
from scipy.optimize import brentq
import matplotlib.pylab as pl
from scipy.integrate import quad


def G(p,z):
    return float(mpmath.polylog(p, z))

G=np.vectorize(G)

def findFugacity(nLanda3,zmin=0,zmax=0.999999999):
    ncLanda3=G(3/2.,1)
    if nLanda3 < ncLanda3:
        
        rootFunction = lambda x :  ( G(3/2,x) - nLanda3 )/ncLanda3
        return brentq(rootFunction,zmin,zmax)
    else:
        return 1
def temperatureBEC(n):
    return 2*pi*(n/zeta(3/2.))**(2./3)

def energy(n, T ):
    beta=1/T
    landa=sqrt(2*pi*beta)
    Tc = temperatureBEC(n)
    z=findFugacity(n*landa**3)
    if n==0:
        return 0
    return 3*T*G(5./2,z)/(2*n * landa**3 )
    
energy=np.vectorize(energy)
        
def energyMixture(n1,n2, T ):
    
    n=n1+n2
    
    return (n1/n)*energy(n1,T) + (n2/n)*energy(n2,T)


def freeEnergyMixtureHF(n1,n2,T,g,g12):
    # returns the free energy per unit volume of a mixture
    xi=G(3/2,1)
    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0_1 = n1 - xi/landa**3
    n0_2 = n2 - xi/landa**3

    y1 = np.exp(  - beta * g * n0_1)
    y2 = np.exp( - beta * g * n0_2)

    if (n1*landa**3 >= xi) and (n2*landa**3 >= xi):
        F = 0.5* g * (n1**2 +  n2**2 ) + g12 * n1 * n2 + g * xi**2/landa**6 - 1/(beta*landa**3)*(G(5/2,y1) + G(5/2,y2) )
    else:
        F=None
    
    
    return F

def susceptibilityMixtureHF(n1,n2,T,g,g12):
    # returns the free energy per unit volume of a mixture
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
    
    return chi



susceptibilityMixtureHF=np.vectorize(susceptibilityMixtureHF)


def chemicalPotentialMixtureHF(n1,n2,T,g,g12):
    # returns the free energy per unit volume of a mixture
    xi=G(3/2,1)
    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0_1 = n1 - xi/landa**3
    n0_2 = n2 - xi/landa**3

    y1 = np.exp( - beta * g * n0_1)
    y2 = np.exp( - beta * g * n0_2)


    if n0_1 > 0:
        mu1= g*n1 + g12*n2 + 1/( landa**3) * g * G(3/2,y1)
    else:
        mu1=None
    
    if n0_2 > 0:
        mu2= g*n2 + g12*n1 + 1/( landa**3) * g * G(3/2,y2)
    else:
        mu2=None
    
    return (mu1,mu2)

chemicalPotentialMixtureHF=np.vectorize(chemicalPotentialMixtureHF)




def freeEnergyHF(n1,T,g):
    xi=G(3/2,1)
    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0 = n1 - xi/landa**3

    y = np.exp( - beta * g * n0)
    
    if (n1*landa**3 >= xi):
        F = 0.5* g * (n1**2 ) + g * xi**2/landa**6 - 1/(beta*landa**3)*(G(5/2,y)  )
    else:
        F=None
    
    
    return F



def chemicalPotentialHF(n1,T,g):
    xi=G(3/2,1)
    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0 = n1 - xi/landa**3

    y = np.exp( - beta * g * n0)
    
    if (n1*landa**3 >= xi):
        mu = g*n1 * (1 - G(3/2,y)/(n1*landa**3)  )
    else:
        mu=None

    return mu


chemicalPotentialHF=np.vectorize(chemicalPotentialHF)


freeEnergyHF=np.vectorize(freeEnergyHF)
freeEnergyMixtureHF=np.vectorize(freeEnergyMixtureHF)


def internalEnergyMixtureHF(n1,n2,T,g,g12):
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



def internalEnergyHF(n1,T,g):
    beta = 1/T
    xi=G(3/2,1)
    landa = np.sqrt(2*pi/T) 
    beta=1/T
    n0 = n1 - xi/landa**3

    y = np.exp( - beta * g * n0)
    
    U = 0.5 * g * ( n1**2 ) -2*g*xi**2/landa**6 + (  3/(2 * beta * landa**3) * ( G(5/2,y )     ) +  (g * xi / landa**6 * 3/2   )* G(3/2,y)   + g/landa**3 * (n1*G(3/2,y)  )  ) 


    return U




def createMatrixMixedStates( ratios, densities ,na3=1e-4,model='HF'):
    # using units of nlanda3==1

    isMiscibleState=np.zeros( (len(ratios),len(densities)))
    for i,ratio in enumerate(ratios):
        for j,density in enumerate(densities):
            P=np.linspace(0,1,num=10)
            n=density
            g=4*pi*(na3/n)**(1/3.)
            g12=ratio*g
            f=freeEnergyHF(n/2*(1+P),n/2*(1-P), 2*pi,g,g12)
            
            P=P[~np.isnan(f)]
            f=f[~np.isnan(f)]
            pMin=P[np.argmin(f)]

            isMiscibleState[i,j]= (1 if pMin == 0 else 0) 
    return isMiscibleState


def createPhaseDiagramPlot( miscibilityMatrix, ratios, n,  model="HF"):
    nc=2*zeta(3/2)
    X,Y=np.meshgrid(ratios,(n/nc)**(-3/2),indexing="ij")
    plt.pcolormesh(X,Y,miscibilityMatrix,linewidth=1,cmap=plt.get_cmap("viridis"))
    plt.ylabel(r"$T/T_c$")
    plt.xlabel(r"$g_{12}/g$")



def H(y,a=1e-8,b=20,plot=False):
    
    f= lambda x : x**(3/2)/(np.exp(np.sqrt(x**2 + 2*y*x)) - 1) *(x+y)/np.sqrt(x**2 + 2*y*x) 
    
    if plot:
        x=np.linspace(a,b,num=1000)
        plt.plot(x,f(x),label=y)
    return 4/(3*np.sqrt(pi)) * quad(f,a,b)[0]


H=np.vectorize(H)


def freeEnergyPopov(n1,n2,T,g,g12):
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
        F = 0.5* g * (n1**2 +  n2**2 ) + g12 * n1 * n2 + g * xi**2/landa**6 - 1/(beta*landa**3)*(H(yPlus) + H(yMinus) )
        F+= 1/(2*pi)**(3/2)*4/(15*sqrt(pi) ) *(    2**(5/2)* (lambdaMinus**(5/2) + lambdaPlus**(5/2))  )

    else:
        F=None

    
    return F


freeEnergyPopov=np.vectorize(freeEnergyPopov)