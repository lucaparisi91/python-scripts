import numpy as np
import matplotlib.pylab as plt
import h5py

class pbc:
    @staticmethod
    def warp(x,L):
       return x- np.round(x/L) * L
    @staticmethod
    def difference(x1,x2,L):
        dx = x1 - x2;
        if L is not None:
            dx -= L * np.round(dx * (1./L));
        return dx





def _pairDistancesAA( x , L ):
    N=len(x)
    dis=np.zeros( (N*(N-1) )//2 )
    k=0
    for i in range( N) :
        for j in range( i ):
            dis[k]=pbc.difference(x[i],x[j],L)
            k+=1
    return dis

def _pairDistancesAB( x , y, L ):
    N1=len(x)
    N2=len(y)
    dis=np.zeros( N1*N2 )
    k=0
    for i in range( N1 ) :
        for j in range( N2 ):
            dis[k]=pbc.difference(x[i],y[j],L)
            k+=1
    return dis


def pairDistances( x , y=None, L=None ):
    if y is None:
        return _pairDistancesAA(x,L)
    else:
        return _pairDistancesAB(x,y,L)



def distance(x):
    return np.sqrt(x[:,0]**2 + x[:,1]**2 + x[:,2]**2)


class particles:

    def __init__(self,filename,lBox):
        self.file=h5py.File(filename)
        self.lBox=lBox
        self.M=self.file["set0"]["particles"].shape[0]
    

    def pairDistances(self, groupA,groupB=None):
        deltas=[]
        M = self.M

        for t in range(M):
            for d in range(3):
                x=self.file[groupA]["particles"][t,d,:]
                if groupB is not None:
                    y=self.file[groupB]["particles"][t,d,:]
                else:
                    y=None
                
                deltax=pairDistances(x,y=y,L=self.lBox[0])

                deltas.append(deltax)
        
        delta=np.zeros( (M, 3,len(deltas[0]) )  )

        k=0
        for t in range(M):
            for d in range(3):
                delta[t,d,:]=deltas[k]
                k+=1       
        return delta
    







class vis:
    def __init__(self,lBox, file):
        self.file=file
        self.lBox=lBox
        self.fig=plt.figure()
        self.ax = self.fig.add_subplot(projection='3d')

    def scatterPlot( self, t , group,*args,**kwds):
        lBox=self.lBox
        ax=self.ax
        ps=self.file[group]["particles"]
        xs=pbc.warp(ps[t,0,:],lBox[0])
        ys=pbc.warp(ps[t,1,:],lBox[1])
        zs=pbc.warp(ps[t,2,:],lBox[2])
        ax.scatter(xs,ys,zs,*args,**kwds)
        ax.set_xlim3d(-self.lBox[0]/2,self.lBox[0]/2)
        ax.set_ylim3d(-self.lBox[1]/2,self.lBox[1]/2)
        ax.set_zlim3d(-self.lBox[2]/2,self.lBox[2]/2)