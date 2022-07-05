from os import system
import numpy as np
from numpy.lib.function_base import average, select
import scipy as sp
from math import *
from scipy import linalg
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags
import copy
import operator
import tqdm
from functools import reduce


def averagePosition(A,J):
    x=spsolve(A,J)

    return x

def correlationMatrix(A):
    return linalg.inv(A.todense())

def averageSquareDistance(i,j,corr=None):
    return corr[i,i] + corr[j,j] - corr[i,j] - corr[j,i]

def averageEnergyHarmonic(A,beta):
    corr=correlationMatrix(A)
    l=A.shape[0]
    timeStep=beta/l
    K=-averageSquareDistance(0,1,corr=corr)*3*l*l/(2*beta**2) + 3/(2 * beta)*l
    V=corr[1,1]*0.5*3

    return K+V

def getCovarianceMean(A,J):
    Sigma=correlationMatrix(A)
    mu=Sigma.dot(J)
    return Sigma,mu

def getAJ(sigma,l,order=2,boundaries="open",harmonic=True):
    timeStep=sigma**2
    A=freeEnergyMatrix(sigma,l,boundaries=boundaries)   
    if harmonic:
        H= harmonicOscillatorMatrix(timeStep,l,boundaries=boundaries,order=order)
        A=A+H
    
    J=freeCurrentVector(sigma,l,boundaries=boundaries)
    return A,J


def getMatrixShape(l,boundaries):
    if not isinstance(boundaries, str):
        l = l + 1 - np.sum([  1 if bound is not None else 0  for bound in boundaries ])
    else:
        l= l if boundaries == "closed" else l + 1
    return l


def freeEnergyMatrix(sigma,l,boundaries="closed"):
    l=getMatrixShape(l,boundaries)
    diagonals=[ np.zeros(l) + 2/sigma**2 , np.zeros(l-1) - 1/sigma**2 , np.zeros(l-1) - 1/sigma**2 ]
    A=diags(diagonals,[0,-1,1],format="csr")    
    
    if boundaries == "open":
        A[0,0]=1/sigma**2
        #A[0,1]=-2/sigma**2
        
        A[-1,-1]=1/sigma**2
        #A[-1,-2]=-2/sigma**2
    if boundaries == "closed":
        A[0,-1]=-1/sigma**2
        A[-1,0]=-1/sigma**2
    
    if not isinstance(boundaries, str):
        leftBound,rightBound=boundaries
        if leftBound is None:
            A[0,0]/=2
        if rightBound is None:
            A[-1,-1]/=2
    
    
    return A

def freeCurrentVector(sigma,l,boundaries="open"):
    l=getMatrixShape(l,boundaries)
    J=np.zeros(l)
    
    if not isinstance(boundaries, str):
        leftBound,rightBound=boundaries

        if (leftBound is not None) :
            J[0]=leftBound/sigma**2
        if (rightBound is not None) :
            J[-1]=rightBound/sigma**2
    return J

            
def harmonicOscillatorMatrix(timeStep,l,boundaries="open",order=2):
    sigma=np.sqrt(1/timeStep)
    l=getMatrixShape(l,boundaries)
    diagonal= np.zeros(l) + 1/sigma**2 
    
    if boundaries == "open":
        if (order == 2):
            diagonal[0]=0.5/sigma**2
            diagonal[-1]=0.5/sigma**2
        else:
            if (order == 1 ):
                diagonal[-1]=0
    if not isinstance(boundaries, str):
        leftBound,rightBound=boundaries

        if leftBound is None:
            diagonal[0]=0.5 if order==2 else 1
        if rightBound is None:
            diagonal[-1]=0.5 if order==2 else 0
    
    H=diags([diagonal],[0])
    return H

def normalization(A):
    return ( (2*pi)**A.shape[0] /linalg.det(A))**(1./2)





def build2BHarmonic(N, ranges=None):
    if ranges is None:
        ranges=(range(0,N) , range(0,N) )
    
    rangeA,rangeB=ranges
    A=np.zeros((N,N))
    
    if rangeA == rangeB:
        N1= len(rangeA)
        A[ rangeA[0]:rangeA[-1] + 1,rangeA[0] : rangeA[-1] + 1 ]=-1
        A[ rangeA, rangeB  ]=N1 - 1
    else:
        N1= len(rangeA)
        N2= len(rangeB)
        
        A[ rangeA[0]:rangeA[-1] + 1,rangeB[0] : rangeB[-1] + 1 ]=-1
        A[ rangeB[0]:rangeB[-1] + 1,rangeA[0] : rangeA[-1] + 1 ]=-1
        
        A[ rangeA, rangeA  ]=N2
        A[ rangeB, rangeB  ]=N1
        
    
    return A


class builderAMatrixGaussianHamiltion:

    def __init__(self,N,M,timeStep):
        '''
        N : number of chains
        M : numer of slices
        head=( iParticle, iTime )
        tail=( iParticle, iTime )
        
        '''
        self.N = N
        self.M = M
        self.A=np.zeros((M,N,M,N))
        self.timeStep=timeStep


        self.next=[ i for i in range(self.N) ]
        self.prev=[ i for i in range(self.N) ]
        self.head=[M-1 for i in range(self.N)]
        self.tail=[ 0 for j in range(self.N)]


    def setHead(self,i,n):
        self.head[n]=i
        oldNext=self.next[n]

        self.next[n]=None

        if oldNext is not  None:
            self.setTail(self.tail[oldNext],oldNext)
        
        

    def setTail(self,i,n):
        self.tail[n]=i

        oldPrev=self.prev[n]
        self.prev[n]=None

        if oldPrev is not  None:
            self.setHead(self.head[oldPrev],oldPrev)

    def join(self,i,j):
        if self.next[i] is not None:
            self.setTail(0, self.next[i] )
        if self.prev[j] is not None:
            self.setHead(self.M-1, self.prev[j])

        
        self.next[i]=j
        self.prev[j]=i
        self.head[i]=self.M-1
        self.tail[j]=0
        


    def addKineticMatrix(self):
        N=self.N
        M=self.M
        sigma=np.sqrt(self.timeStep)
        Ts=[freeEnergyMatrix(sigma,l=self.M,boundaries="closed").toarray() for n in range(N) ]
        # copy in the bulk of the T matrices
        T=np.zeros( (M,N,M,N) )

        for n in range(N):
            T[ 1:self.M - 1  , n, :,n] = Ts[n][1:self.M-1,:]
       
    
        kineticMatrixCoefficient=1./(self.timeStep)


        # apply periodic  boundary conditions in time
        for n in range(N):

            diagonalCoefficient=2
            if self.next[n] is not None:
                # no pbc
                T[ M - 1  , n , 0 ,self.next[n]] = - 1 * kineticMatrixCoefficient
            else:
                # pbc
                diagonalCoefficient=1

            if self.head[n] == self.tail[n]:
                diagonalCoefficient=1

            # set diagonal and link with previous time slice
            T[ self.head[n] , n , self.head[n] - 1 , n ]=-1 * kineticMatrixCoefficient
            T[ self.head[n] , n , self.head[n] , n ]=diagonalCoefficient* kineticMatrixCoefficient

            if self.prev[n] is not None:
                T[ 0  , n , M-1 ,self.prev[n] ] = - 1* kineticMatrixCoefficient
                diagonalCoefficient=2
            else:
                diagonalCoefficient=1
            
            if self.head[n] == self.tail[n]:
                diagonalCoefficient=1

            if self.tail[n] +1 >= M:
                T[ self.tail[n] , n , 0 , self.next[n] ]=-1* kineticMatrixCoefficient
            else:
                T[ self.tail[n] , n , self.tail[n] + 1 , n ]=-1* kineticMatrixCoefficient
            
            
            T[ self.tail[n] , n , self.tail[n] , n ]=diagonalCoefficient* kineticMatrixCoefficient


        # erase where head and tail are present

            
            
            

        
        
        
        self.A+=T

        self._eraseHeadTail()



    def _eraseHeadTail(self):

        for n in range(self.N):
            
            tail = self.tail[n]
            head = self.head[n]
            
            if (tail > 0):
                self.A[0:tail, n, : , n ]=0
                self.A[ tail , n, tail - 1 , n ]=0
            if (head < self.M-1):
                self.A[head+1: , n, : , n ]=0
                self.A[ head , n, head +1 :, n ]=0
        






    def add2BHarmonicInteraction(self,ranges=None):
        N=self.N
        M=self.M
        timeStep=self.timeStep
        As=[ build2BHarmonic(N,ranges=ranges) for t in range(M) ]
        
        A=np.zeros((M,N,M,N))

        for t in range(M):
            A[t,:,t,:] = As[t]

        mask = self.make2BMask()
        
        self.A+=A*mask*timeStep


    def makeBeadsMask(self):
        '''
        if extenedTail = false
        return a mask(M, N) with 0 (deleted bead or head) , 1 (active bead or tail)
 
        '''
        mask=np.zeros( (self.M,self.N))
        for n in range(self.N):
            mask[self.tail[n]:self.head[n] + 1,n]=1
            if self.next[n] is None:
                mask[self.head[n] ,n ]=0
            if self.prev[n] is None:
                mask[self.tail[n] ,n ]=1
        
        return mask
    

    def make2BMask(self):
        # extends the beads mask for imaginary boundary conditions at t=0 extremity
        beadsMaskExtended=np.zeros((self.M+1,self.N))
        beadsMaskExtended[1:,:]=self.makeBeadsMask()

        for n in range(self.N):
            beadsMaskExtended[0,n]= 0 if  self.prev[n] is None else 1 

        
        mask2B=(np.multiply.outer(beadsMaskExtended[1:,:],beadsMaskExtended[1:,:]) + np.multiply.outer(beadsMaskExtended[:self.M,:],beadsMaskExtended[:self.M,:]) )*0.5

        # adjust the diagonal elements
        for t in range(self.M):
            for n in range(self.N):
                mask2B[t,n,t,n]= (np.sum(mask2B[t,n,t,:]) - mask2B[t,n,t,n])/(self.N-1)

        return mask2B


    
    def add1BHarmonicInteraction(self):
        Vs =[  harmonicOscillatorMatrix(self.timeStep,l=self.M, boundaries="closed" ,order=2).todense() for n in range(self.N) ]


        for n in range(self.N):
            tail = self.tail[n]
            head = self.head[n]

            V=Vs[n]

            if self.prev[n] is None:
                    
                V[tail,tail]*=0.5

            if self.next[n] is None:    
                V[head,head]*=0.5


            self.A[:,n,:,n]+=V
        
        self._eraseHeadTail()

    

class gaussianObservables:

    def __init__(self,system,dimensions=1):
        self.system=system
        A=self.system.A

        self.N=A.shape[-1]
        self.M=A.shape[0]
        self.A=A
        self.timeStep=self.system.timeStep
        
        self.computeInverse()
        self.dimensions=dimensions


    def computeInverse(self):
        self.shape=(self.M,self.N,self.M,self.N)
        matrixShape=(self.shape[0]*self.shape[1],self.shape[2]*self.shape[3] )
        
        A_matrix=self.A.reshape(matrixShape)
        row_mask=np.all(A_matrix==0 , axis=0)
        
        reducedRows= matrixShape[0] - np.count_nonzero(row_mask)
        reducedShape=(reducedRows,reducedRows)

        self.mask=np.logical_or.outer(row_mask,row_mask).reshape(self.shape)

        A_reduced=self.A[~self.mask].reshape(reducedShape)

        corr_reduced=linalg.inv(A_reduced)

        self.corr=np.zeros(self.shape)
        self.corr[~self.mask]=corr_reduced.ravel()

        self.A_reduced=A_reduced

    def numberOfLinks(self,n=None):

        if n is None:
            l=0
            for n in range(self.N):
                l+=self.numberOfLinks(n)
            return l

        l=self.system.head[n] - self.system.tail[n]
        if self.system.next[n] is not None:
            l+=1
        return l

    def averageSquareDistance(self, beadA , beadB ):
        i,nA= beadA
        j,nB = beadB
        return (self.corr[i,nA,i,nA] + self.corr[j,nB,j,nB] - self.corr[i,nA,j,nB ] - self.corr[j,nB,i,nA])*self.dimensions



    
    def centerOfMassSquared(self):
        cm2=0
        for t in range(self.M):
            for i in range(self.N):
                cm2+=self.corr[t,i,t,i]
        return cm2/(self.M*self.N)



    def averageTwoBodySquareDistance(self):
        dis2=0
        for i in range(self.N):
            for j in range(0,  i):
                for t in range(self.M):
                    dis2+=self.averageSquareDistance(  (t,i) , (t,j)  )
        return dis2


    def energy(self,oneBodyHarmonic=True,twoBodyHarmonic=False):

        T = self.averageLengthSquare() / (2* self.timeStep)
        V=0

        if oneBodyHarmonic:
            V+= 0.5 * self.centerOfMassSquared() * self.N * self.M
        

        if twoBodyHarmonic:
            V+= 0.5 * self.averageTwoBodySquareDistance()
        
        
        beta=self.timeStep * self.M

        return V/self.M - T/(beta) + self.dimensions/(2*self.timeStep)*self.N

        



    def averageChainLengthSquare(self,n, t0=None , l=None):

        head = self.system.head[n] 
        tail = self.system.tail[n] 

        if l==0:
            return 0

        if t0 is None:
            t0=tail
        if l is None:
            t1=head
        else:
            t1=min(head,t0+l)
        
        t0=max(tail,t0)
        l2=0
        for t in range(t0, t1 ):
            l2+=self.averageSquareDistance( (t,n),(t+1,n) )

        if self.system.next[n] is not None and t0+l>head:        
            # boundary condition
            l2+=self.averageSquareDistance( (self.M-1,n),(0,self.system.next[n]) )
            l2+=self.averageChainLengthSquare(self.system.next[n],l=t0+ l - (t1-t0+1)     )
        return l2

    

    def averageLengthSquare(self,n=None):
        l2=0
        for n in range(self.N):
            l2+=self.averageChainLengthSquare(n=n,t0=None,l=None)
        return l2

    def normalization(self):
        l=self.numberOfLinks()
        Z=np.exp(  0.5 * ( np.log(2*pi)*(self.A_reduced.shape[0]) -np.linalg.slogdet(self.A_reduced))[1]   - np.log(2*pi*self.timeStep)*(l/2.) )

        return Z**(self.dimensions)

def createAMatrix(M, N , timeStep, trapped=True , twoBodyHarmonicInteraction=True):
    builder = builderAMatrixGaussianHamiltion(N=N,M=M,timeStep=timeStep)
    builder.addKineticMatrix()
    if trapped:
        builder.add1BHarmonicInteraction()
        
    if twoBodyHarmonicInteraction:
        builder.add2BHarmonicInteraction()
    return builder



class testConfigurations:
    def __init__(self,M,timeStep,trapped=True,twoBody=False,dimensions=1):
        self.M=M
        self.timeStep=timeStep
        self.trapped=trapped
        self.twoBody=twoBody
        self.dimensions=dimensions

    def build(self,N):
        return builderAMatrixGaussianHamiltion(N,M=self.M,timeStep=self.timeStep)
    def createObservable(self,b):
        b.addKineticMatrix()
        if self.trapped:
            b.add1BHarmonicInteraction()
        if self.twoBody and b.N>1:
            b.add2BHarmonicInteraction()

        return gaussianObservables(b,dimensions=self.dimensions)


class testAdvanceRecedeWorm(testConfigurations):
    def __init__(self,lengths,*args,**kwds) : 
        super(testAdvanceRecedeWorm,self).__init__(*args,**kwds)
        self.lengths=lengths
        
    def __iter__( self):
        for l in tqdm.tqdm(self.lengths):
            N=l//self.M + 1
            b=self.build(N)
            b.setTail(0,0)
            b.setHead(l%self.M,N-1)
            for i in range(0,N-1):
                b.join(i,i+1)
            
            yield self.createObservable(b)


class testAdvanceSwappedRecedeWorm(testConfigurations):
    def __init__(self,*args,**kwds) : 
        super(testAdvanceSwappedRecedeWorm,self).__init__(*args,**kwds)
    
    def __iter__( self):
        # unswapped
        lengths=range(1,2*self.M)

        for l in lengths:

            N = l//(self.M) + 1
            b=self.build(N)
            b.setTail(0,0)
            b.setHead(l%self.M,N-1)
            for i in range(0,N-1):
                b.join(i,i+1)
            
            yield self.createObservable(b) 
        
         
        #swapped
        for l in range(1,self.M):
            N = 2
            b=self.build(N)
            b.join(1,1)
            b.setHead(l,0)
            b.setTail(0,0)
            yield self.createObservable(b)



def nLinksInChain(b,n=None):

  

    l=0
    nCurrent=n        
    while True:
        l+=b.head[nCurrent] - b.tail[nCurrent] 
        nCurrent=b.next[nCurrent]
        
        if nCurrent is not None:
            l+=1

        if (nCurrent is None) or (nCurrent==n):
            break

    
    return l

def nLinks(b,n=None):

    if n is None:
        return reduce(operator.add, map( lambda n : nLinks(b,n) ,range(0,b.N)   ) )

    l=b.head[n] - b.tail[n] 
    if b.next[n] is not None:
        l+=1
    return l
    
def windingCoeff(beta=1,nMax=100,L=1):
    n=np.arange(1,nMax)
    return 1 + 2*np.sum( np.exp(-L**2*n**2/(2*beta)) )

def eWindingCoeff(beta,nMax=100,L=1):
    n=np.arange(1,nMax)
    x=n**2*L**2/(2*beta**2) * 2
    Zs=np.exp(-L**2*n**2/(2*beta))
    return np.sum( x*Zs)/(1 + 2*np.sum(Zs))


def _ZFree(M=10,beta=1,V=1,winding=True):
    Z=V/(2*np.pi*beta)**(3/2)
    L=V**(1/3.)
    if winding:
        Z*= (windingCoeff(beta=beta,nMax=200,L=L)**3 )
    return Z


def ZFree(beta,V,N):
    Z=np.zeros(N+1)
    Z[0]=1

    for n in range(1,N+1):
        for k in range(1,n+1):
            Z[n]+=_ZFree(beta=k*beta,V=V)*Z[n-k]
        Z[n]/=n
    return Z[N]




def singleParticleEnergy(M=10,beta=1,V=1,winding=True):
    e=1/(2*beta)
    L=V**(1/3)
    if (winding):
        e-=eWindingCoeff(beta=beta,nMax=400,L=L)
    return e*3

def energyFree(beta,V,N=1):
    Z=np.zeros(N+1)
    e=np.zeros(N+1)
    Z[0]=1

    for n in range(1,N+1):
        for k in range(1,n+1):
            Z[n]+=_ZFree(beta=k*beta,V=V)*Z[n-k]
        for k in range(1,n+1):    
            Z0=_ZFree(beta=k*beta,V=V)
            e[n]+=k*singleParticleEnergy(V=V,beta=k*beta)*Z0*Z[n-k] + e[n-k]*Z[n-k] *Z0
    e/=Z[N]

    return e[N]/N








