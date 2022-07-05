from math import cos
import numpy as np
import copy
from numpy.linalg import norm


class caoBernePropagator:

    def __init__(self,a):
        self.a=a
        self.D=1

    def __call__( self, x , y , tau ):
        r=norm(x)
        rp=norm(y)
        cosTeta=np.dot(x,y)/(norm(x)*norm(y) )
        
        return 1 - (self.a*(r+rp) - self.a**2 )/ (r*rp) * np.exp( -( r * rp + self.a**2 -self.a*(r + rp))*(1+cosTeta)/(2*self.D*tau))
    
    def timeDerivative(self,x,y,tau,dt=1e-3):
        return (self.__call__(x,y,tau + dt) - self.__call__(x,y,tau - dt) )/(2*dt)


if __name__ == "__main__":
    x=np.array([1.3,1.5,0.05])
    y =np.array( [1.6, 1.7,0.1] )
    #x=[2,0,0]
    #y=[-2,0,0]

    tau=0.1

    G = caoBernePropagator(a=1)
    
    g0=G(x,y ,tau=tau)
    #print( np.log(g0) )


    #print ( G.timeDerivative(x,y,tau=0.1,dt=1e-6)/g0 )

    dx=np.array([1e-6,0,0])


    print(  (G(x,y+dx,tau) - G(x,y-dx,tau) )/(2*dx[0]*G(x,y,tau))  )
    print(G(x,y,tau))
    print(G(y,x,tau))
    
