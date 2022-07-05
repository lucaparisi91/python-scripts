from pytest import param
from . import simulation
from . import action
from . import observable
from . import moves
import argparse
import pandas as pd
import numpy as np
from .theory import free
from . import inputFileTools
import os

def createSim( a,ratio, boxSize, N, T, C , nBeads,deltaMu,polarization_range=None,nBlocks=10000):
    ensamble=simulation.semiCanonicalEnsamble(N=N,boxSize=boxSize,T=T,deltaMu=deltaMu)
    run=simulation.run(nBlocks=nBlocks,stepsPerBlock=1000,correlationSteps=N)
    run.saveConfigurations=False
    a12=a*ratio
    actions=[action.caoBerneAction(a=a,groups=[0,0] ),action.caoBerneAction(a=a,groups=[1,1]),action.caoBerneAction(a=a12,groups=[0,1] ) ]
        
    observables=[observable.thermodynamicEnergy(label="energy",magnetizationDistribution=True,N=N) , observable.virialEnergy(label="eV",magnetizationDistribution=True,N=N) , observable.magnetizationDistribution(groups=[0,1],label="M",magRange=[-N,N]    ) ]

    model=simulation.model( ensamble=ensamble,actions=actions,nBeads=nBeads)

    restriction=None
    if polarization_range is not None:
        pMin=polarization_range[0]
        pMax=polarization_range[1]
        nA_range=[ int( N/2*(1+pMin) ), int( N/2*(1+pMax))]
        nB_range=[ N  - nA_range[1], N - nA_range[0] ]
        restriction=moves.restriction( sets=[0,1],particleRanges=[nA_range, nB_range ])
    

    tab= moves.createTableSemiCanonical(C=C,l=int(0.6*nBeads),lShort=int(0.3*nBeads),groups=[0,1],delta=0.3*boxSize[0],restriction=restriction)
    
    sim=simulation.simulation(model=model,run=run, observables=observables,moves=tab,N0=[ nA_range[0],N-nA_range[0]  ] )
    
    return(sim)



def generateInputFiles(data):
    js=[]
    for i,row in data.iterrows():
        N=int(row["N"])
        T=float(row["T"])
        na3=float(row["na3"])
        C=[ float(row["CA"]),float(row["CB"]),float(row["CAB"]) ]
        nBeads=int(row["nBeads"])
        ratio=int(row["ratio"])
        nBlocks=int(row["nBlocks"])
        polarization_range=None

        if "pMin" in row.keys():
            polarization_range=[float(row["pMin"]),float(row["pMax"])]
        
        landaC=np.sqrt(2*np.pi)
        L=( (N/2) * landaC**3/free.G(3/2,1) )**(1/3)
        n=N/L**3
        a=(na3/n)**(1/3)

        sim=createSim(a=a,N=N,T=T,boxSize=[L,L,L] ,C=C,nBeads=nBeads,ratio=ratio,deltaMu=0,polarization_range=polarization_range,nBlocks=nBlocks)

        js.append(sim.toJson())



    return(js)



def generateLabels(data):
    labels=[]
    for i,row in data.iterrows():
        values=[ value for value in row.values ]

        formatStrings=["{}{{}}".format( key) for key in row.keys() ]   
        formatString="_".join(formatStrings)
        #return values
        labels.append(formatString.format(*values) )
    return (labels)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Generate list of json files')
    parser.add_argument('parameters_file', type=str )

    parser.add_argument('--folder', type=str , default="." )

    args = parser.parse_args()
        
    parameters=args.parameters_file

    data=pd.read_csv(parameters,delim_whitespace=True)

    print(data)

    js=generateInputFiles(data)

    labels=generateLabels(data)


    settings=[ {"folder" : os.path.join(args.folder,label) , "jSon" : [ ["input.json", j  ] ] } for label,j in zip(labels,js) ]

    inputFileTools.createSimFolders(settings)








