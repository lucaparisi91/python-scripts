import numpy as np
import json
import random


class semiCanonicalEnsamble:
    def __init__(self,N,boxSize,T,deltaMu):
        self.N=N
        self.deltaMu=deltaMu
        self.T=T
        self.boxSize=boxSize


class canonicalEnsamble:
    def __init__(self,N, boxSize, T):
        self.N=N
        self.boxSize=boxSize
        self.T=T

class run:
    def __init__(self, nBlocks,stepsPerBlock,correlationSteps=None,maxTime=None):
        self.nBlocks=nBlocks
        self.stepsPerBlock=stepsPerBlock
        self.correlationSteps=correlationSteps
        self.saveConfigurations=True
        self.maxTime=maxTime


class model:
    def __init__(self, ensamble , nBeads , actions,nCells=None ):
        self.ensamble=ensamble
        self.actions=actions
        self.nBeads=nBeads
        self.nCells=nCells
    


class simulation:
    def __init__( self, model, run , moves,observables=None,folder=".",N0=None):
        self.model=model
        self.run=run
        self.observables=observables
        self.moves=moves
        self.seed=567
        self.N0=N0

    def toJson(self ):
        j={}
        N0=self.N0
        N0_max=[]
        
        if type(self.model.ensamble) == semiCanonicalEnsamble:
            if self.N0 is None:
                N0=[self.model.ensamble.N/2,self.model.ensamble.N - self.model.ensamble.N/2 ]
            
            N0_max = [ np.sum(N0) + 2 for n in N0 ]
        else:
            if type(self.model.ensamble) == canonicalEnsamble:
                if N0 is not None:
                    if N0 != self.model.ensamble.N:
                        raise RuntimeError("Invalid initial atom number.")
                else:
                    N0=self.model.ensamble.N
            
                N0_max = [ n+2 for n in N0 ]
            else:
                raise RuntimeError( "Unkown ensamble {}".format(type(self.model.ensamble)) )


        j["inverseTemperature"]=float(1/self.model.ensamble.T)
        j["nBeads"]=int(self.model.nBeads)
        j["stepsPerBlock"]=int(self.run.stepsPerBlock)
        j["correlationSteps"]=int(self.run.correlationSteps)
        j["nBlocks"]=int(self.run.nBlocks)
        j["particles"]=[int(n) for n in N0]
        j["ensamble"]="grandCanonical"
        j["checkPointFile"]="latest.hdf5"
        j["saveConfigurations"]=self.run.saveConfigurations
        j["chemicalPotential"]=[0 for n in range(len(N0))]

        if type(self.model.ensamble) == semiCanonicalEnsamble:
            j["chemicalPotential"][0]=self.model.ensamble.deltaMu
        
        j["maxParticles"]=[ int(n) for n in N0_max ]
        j["checkPointFile"]="latest.hdf5"
        j["lBox"]=[ float(L) for L in self.model.ensamble.boxSize ]

        j["action"]= [action.toJson() for action in self.model.actions]
        j["observables"] = [ ob.toJson() for ob in self.observables ]

        if len(self.model.actions) !=0:
            minimumDistance=max( [S.minimumDistance for S in self.model.actions ]  )
        else:
            minimumDistance=0
        
        j["movesTable"]=self.moves.toJson()
        j["randomInitialCondition"]= {
            "minimumDistance" :minimumDistance
        }
        j["seed"]=int(random.randint(1, 1000000))


        if self.run.maxTime is not None:
            j["maxTime"]=int(self.run.maxTime)
        

        if self.model.nCells is not  None:
            if len(N0) > 1 :
                raise RuntimeError("Mesh action only implemented for one component(so far)")
            j["nCells"]=[self.model.nCells for d in range(3)]

            for action in self.model.actions:
                if not action.mesh:
                    raise RuntimeError("Not all actions are compatible with meshed configurations")



        return j









