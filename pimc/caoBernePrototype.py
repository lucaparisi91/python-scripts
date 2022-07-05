import scipy as sp
import numpy as np
from numpy.random import default_rng
import tqdm
from scipy import sparse
import scipy.linalg
import caoBerne as cao


"""
Prototype for sampling two particles without exchanges using the Cao Berne propagator
"""


def createJVector(M,tau,x0):
    J=np.zeros(M-1)
    J[0]=x0/(tau)
    J[-1]=x0/(tau)
    return J


def createKineticGaussianMatrix(M,tau,fixed=True):
    if fixed:
        A=sp.sparse.diags( 
            [ np.ones( M-1)*2, -np.ones(M-2),-np.ones(M-2) ], 
            offsets=[0,-1,1], shape=(M-1,M-1)  
        )
        A/=tau
        
        return A.todense()
    else:
        A=sp.sparse.diags( 
            [ np.ones( M)*2, -np.ones(M-1),-np.ones(M-1) ], 
            offsets=[0,-1,1], shape=(M,M)  
        )
        A=A.todense()
        A[0,-1]=-1
        A[-1,0]=-1
        
        A/=tau
    
        return A 


def sampleGaussianChain(M,timeStep, lBox,rng=default_rng()):
    A = createKineticGaussianMatrix(M,timeStep,fixed=True)
    Sigma=sp.linalg.inv(A )
    D=len(lBox)
    xs=np.zeros( (M,D) )
    for d in range(D):
        x0=(rng.random()-0.5)*lBox[d]
        
        J= createJVector(M,timeStep,x0)
        mu=np.dot(Sigma,J)

        x=np.random.multivariate_normal(mu,Sigma)
        xs[1:,d]=x
        xs[0,d]=x0
    return xs

def distance(deltax,lBox):
    deltaxPBC=np.zeros(deltax.shape)
    D=deltax.shape[1]
    for d in range(D):
        deltaxPBC[:,d] = deltax[:,d] -  lBox[d] * np.round(deltax[:,d] * (1./lBox[d]));
    return deltaxPBC

def radius(deltax):
    r2=np.sum(deltax**2,axis=1)
    return np.sqrt(r2)

def averageLengthSquare(x):
    M=x.shape[0]
    A = createKineticGaussianMatrix(M,1,fixed= False )
    return np.sum(np.tensordot(  A , x , axes=(1,0) )*x)


class caoBernePropagator:

    def __init__(self,a):
        self.a=a
        self.D=1
        self._G=cao.caoBernePropagator(a)
        
    def __call__( self, x , tau ):
        M,D=x.shape
        y=np.zeros(x.shape)
        y[0:-1,: ]=x[1:,:]
        y[-1,:]=x[0,:]

        rx=radius(x)
        ry=radius(y)
        
        cosTeta=np.sum(x*y,axis=1)/(rx*ry)
        
        return 1 - (self.a*(rx+ry) - self.a**2 )/ (rx*ry) *np.exp( -( rx * ry + self.a**2 -self.a*(rx + ry))*(1+cosTeta)/(2*self.D*tau))


        
    def timeDerivative(self,x,tau,dt=1e-3):
        return (self.__call__(x,tau + dt) - self.__call__(x,tau - dt) )/(2*dt)

    def logTimeDerivative(self,x,tau,dt=1e-3):
        return self.timeDerivative(x,tau,dt=dt)/self.__call__(x,tau)

    
    def force(self,x,tau,dx=1e-4):
        M,D = x.shape
        f=np.zeros(x.shape)

        for d in range(D):
            for t in range(M):
                dx_=[0,0,0]
                dx_[d]=dx

                f[t,d]+=(self._G(x[t,:] + dx_ , x[(t+1)%M,:],tau) - self._G(x[t,:] - dx_ , x[(t+1)%M,:] , tau ) ) /(2*dx * self._G(x[t,:] , x[(t+1)%M,:] , tau ) )
                f[t,d]+=(self._G(x[(t-1)%M,:]  , x[t,:] + dx_, tau ) - self._G(x[(t-1)%M,:]   , x[t,:] - dx_ ,tau ) )/( 2*dx *self._G( x[(t-1)%M,:] , x[t,:] , tau) )

        return -f



def energy( prop, x1, x2 , tau, lBox):
    M,_=x1.shape
    T1=averageLengthSquare(x1)/(2*tau)
    T2=averageLengthSquare(x2)/(2*tau)
    T=T1 + T2
    deltax=distance(x2-x1,lBox)
    V = - np.sum(prop.logTimeDerivative(deltax,tau,dt=1e-6) )/M
    beta=M*tau
    T/=beta

    return -T + V +  3./(2*tau)*2  

def virialEnergy(G,x1,x2,tau,lBox):
    M,D=x1.shape
    beta=M*tau
    deltax=distance(x2-x1,lBox)
    force=G.force(deltax,tau,dx=1e-6)
    V = - np.sum( G.logTimeDerivative(deltax,tau,dt=1e-6) )/M

    Rc1=np.mean(x1,axis=0)
    Rc2=np.mean(x2,axis=0)

    e3 = -   np.sum( (x1 - Rc1) * force )
    e3 += np.sum( (x2 - Rc2) * force)



    e3/=(2*beta)

    return D/(2*beta)*2 + e3  + V
    

    
def sampleInitialCondition(a,M,timeStep,lBox,maxIteration=1000):
    i=0
    while True:
        x1=sampleGaussianChain(M=M,timeStep=timeStep,lBox=lBox)
        x2=sampleGaussianChain(M=M,timeStep=timeStep,lBox=lBox)
        deltaxPBC=distance(x2-x1,lBox)
        aMin=np.min(radius(deltaxPBC))
        if aMin > a:
            break
        i+=1
        if i> maxIteration:
            raise RuntimeError("Max iteration reached")
    return x1,x2
    
    
def metropolis_log(p,rng=default_rng()):
    accept=False
    if p >= 0:
            accept=True
    else:
            r=rng.random()
            if p> np.log(r):
                accept=True
            else:
                accept=False
    return accept



if __name__ == "__main__":

    M=10
    timeStep=0.1
    a=0.1
    lBox=[1,1,1]
    nSteps=100000
    subSteps=100

    G=caoBernePropagator(a)
    x1,x2=sampleInitialCondition(a,M,timeStep,lBox)
    deltaxPBC=distance(x2-x1,lBox)
    rng=default_rng()
    accept=True
    p=np.sum(np.log(G(deltaxPBC,timeStep) ))

    for i in range(nSteps):
        acceptanceRatio=0
        l2_acc=0
        dis2_acc=0
        e_acc=0
        eV_acc=0
        for ii in range(subSteps) :
            x1New=sampleGaussianChain(M=M,timeStep=timeStep,lBox=lBox,rng=rng)
            x2New=sampleGaussianChain(M=M,timeStep=timeStep,lBox=lBox,rng=rng)
            deltaxPBC=distance(x2New-x1New,lBox)
            aMin=np.min(radius(deltaxPBC))
            p2=p
            
            if aMin > a:
                p2=np.sum(np.log(G(deltaxPBC,timeStep) ))
                accept=metropolis_log(p2 - p ,rng=rng)
            else:
                accept=False

            if accept:
                acceptanceRatio+=1
                x1=x1New
                x2=x2New
                p=p2
                
            # make measurements
            
            l2_acc+=averageLengthSquare(x1)
            dis2_acc+= np.sum( distance(x2-x1,lBox)**2 )/(M*2) 
            e_acc+=energy( G, x1, x2 , timeStep, lBox)
            eV_acc+=virialEnergy( G, x1, x2 , timeStep, lBox)


        print("------")
        print ("{}/{}".format(i,nSteps))
        acceptanceRatio/=subSteps
        l2_acc/=subSteps
        dis2_acc/=subSteps
        e_acc/=subSteps
        eV_acc/=subSteps

        with open("l2.dat","a") as f:
            f.write("{} {}\n".format(i,l2_acc)  )
        with open("dis2.dat","a") as f:
            f.write("{} {}\n".format(i,dis2_acc)  )
        with open("energy.dat","a") as f:
            f.write("{} {}\n".format(i,e_acc)  )
        with open("eV.dat","a") as f:
            f.write("{} {}\n".format(i,eV_acc)  )
        
        
        print("Acceptance ratio: {}".format(acceptanceRatio))